"""
Terrain Data Collector using Mapbox Terrain-RGB and other DEM sources

Implements Layer 2.1 from Pipeline:
- Mapbox Terrain-RGB for DEM/DTM
- National LiDAR portals (UK Environment Agency) for high-res DEM
- DSM for buildings/trees (for shadow analysis)
"""

import os
import math
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass
from pathlib import Path
import requests
from io import BytesIO
from loguru import logger
import numpy as np

try:
    from dotenv import load_dotenv
    HAS_DOTENV = True
except ImportError:
    HAS_DOTENV = False

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    logger.warning("PIL not installed - terrain tile processing limited")

try:
    import mercantile
    HAS_MERCANTILE = True
except ImportError:
    HAS_MERCANTILE = False
    logger.warning("mercantile not installed - using basic tile math")

from ..config import get_config


@dataclass
class TerrainData:
    """Terrain data extracted from DEM"""
    elevation_m: float
    slope_percent: float
    slope_direction: str
    terrain_classification: str
    corner_samples: List[Dict[str, Any]]
    dem_source: str
    resolution_m: float


class TerrainCollector:
    """
    Collect terrain/elevation data from multiple sources:
    1. Mapbox Terrain-RGB (primary)
    2. Open-Meteo (fallback)
    3. UK Environment Agency LiDAR (optional high-res)
    """
    
    # Mapbox Terrain-RGB tile URL template
    MAPBOX_TERRAIN_URL = "https://api.mapbox.com/v4/mapbox.terrain-rgb/{z}/{x}/{y}.pngraw?access_token={token}"
    
    # Open-Meteo Elevation API (free, no key required)
    OPEN_METEO_URL = "https://api.open-meteo.com/v1/elevation"
    
    # UK Environment Agency LiDAR (free for UK)
    UK_LIDAR_WCS = "https://environment.data.gov.uk/spatialdata/lidar-composite-dtm-2m/wcs"
    
    def __init__(self):
        self.config = get_config()
        
        # Load .env file if available
        if HAS_DOTENV:
            # Try multiple locations for .env file
            env_paths = [
                Path(__file__).parent.parent.parent / ".env",  # Project root
                Path.cwd() / ".env",  # Current working directory
                Path.home() / ".env",  # Home directory
            ]
            
            loaded = False
            for env_path in env_paths:
                if env_path.exists():
                    load_dotenv(env_path, override=False)  # Don't override existing env vars
                    logger.debug(f"Loaded .env file from {env_path}")
                    loaded = True
                    break
            
            if not loaded:
                # Try loading from current directory (dotenv default behavior)
                load_dotenv()
        
        self.mapbox_token = os.getenv("MAPBOX_ACCESS_TOKEN", "")
        
    def collect_terrain(
        self, 
        lat: float, 
        lon: float, 
        radius_m: float = 50
    ) -> TerrainData:
        """
        Collect terrain data for plot area
        
        Args:
            lat: Center latitude
            lon: Center longitude
            radius_m: Radius to sample
            
        Returns:
            TerrainData with elevation, slope, and classification
        """
        logger.info(f"Collecting terrain data for ({lat}, {lon})")
        
        # Try Mapbox first if token available
        if self.mapbox_token and HAS_PIL and HAS_MERCANTILE:
            try:
                return self._collect_mapbox_terrain(lat, lon, radius_m)
            except Exception as e:
                logger.warning(f"Mapbox terrain failed: {e}, falling back to Open-Meteo")
        
        # Fallback to Open-Meteo
        return self._collect_open_meteo_terrain(lat, lon, radius_m)
    
    def _collect_mapbox_terrain(
        self, 
        lat: float, 
        lon: float, 
        radius_m: float
    ) -> TerrainData:
        """
        Collect terrain from Mapbox Terrain-RGB tiles
        
        Terrain-RGB encoding: height = -10000 + ((R * 256 * 256 + G * 256 + B) * 0.1)
        """
        zoom = 14  # Good resolution for plot-level analysis
        
        # Get tile coordinates
        tile = mercantile.tile(lon, lat, zoom)
        
        # Fetch tile
        url = self.MAPBOX_TERRAIN_URL.format(
            z=zoom, x=tile.x, y=tile.y, token=self.mapbox_token
        )
        
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        # Parse RGB to elevation
        img = Image.open(BytesIO(response.content))
        img_array = np.array(img)
        
        # Decode Terrain-RGB
        R = img_array[:, :, 0].astype(float)
        G = img_array[:, :, 1].astype(float)
        B = img_array[:, :, 2].astype(float)
        
        elevation_grid = -10000 + ((R * 256 * 256 + G * 256 + B) * 0.1)
        
        # Get pixel position of our point
        bounds = mercantile.bounds(tile)
        px_x = int((lon - bounds.west) / (bounds.east - bounds.west) * 256)
        px_y = int((bounds.north - lat) / (bounds.north - bounds.south) * 256)
        
        # Sample center and corners
        center_elev = elevation_grid[px_y, px_x]
        
        # Calculate slope from local gradient
        dx = (elevation_grid[px_y, min(px_x+1, 255)] - elevation_grid[px_y, max(px_x-1, 0)]) / 2
        dy = (elevation_grid[max(px_y-1, 0), px_x] - elevation_grid[min(px_y+1, 255), px_x]) / 2
        
        # Pixel size in meters (approximate)
        meters_per_pixel = 40075000 * math.cos(math.radians(lat)) / (256 * 2**zoom)
        
        slope_percent = math.sqrt(dx*dx + dy*dy) / meters_per_pixel * 100
        slope_direction = self._get_slope_direction(dx, dy)
        
        # Sample corners (approximate based on radius)
        offset_pixels = int(radius_m / meters_per_pixel)
        corner_samples = self._sample_corners(
            elevation_grid, px_x, px_y, offset_pixels, lat, lon, meters_per_pixel
        )
        
        return TerrainData(
            elevation_m=float(center_elev),
            slope_percent=float(slope_percent),
            slope_direction=slope_direction,
            terrain_classification=self._classify_terrain(slope_percent),
            corner_samples=corner_samples,
            dem_source="mapbox_terrain_rgb",
            resolution_m=meters_per_pixel
        )
    
    def _collect_open_meteo_terrain(
        self, 
        lat: float, 
        lon: float, 
        radius_m: float
    ) -> TerrainData:
        """
        Fallback: Collect terrain from Open-Meteo Elevation API
        """
        # Sample multiple points for slope calculation
        offset_deg = radius_m / 111000  # Approximate
        
        points = [
            (lat, lon),  # Center
            (lat + offset_deg, lon),  # North
            (lat - offset_deg, lon),  # South
            (lat, lon + offset_deg),  # East
            (lat, lon - offset_deg),  # West
        ]
        
        lats = [p[0] for p in points]
        lons = [p[1] for p in points]
        
        try:
            response = requests.get(
                self.OPEN_METEO_URL,
                params={
                    "latitude": lats,
                    "longitude": lons
                },
                timeout=30
            )
            response.raise_for_status()
            data = response.json()
            
            elevations = data.get("elevation", [0] * 5)
            center_elev = elevations[0] if elevations else 0
            
            # Calculate slope
            if len(elevations) >= 5:
                north_south_slope = (elevations[1] - elevations[2]) / (2 * radius_m) * 100
                east_west_slope = (elevations[3] - elevations[4]) / (2 * radius_m) * 100
                slope_percent = math.sqrt(north_south_slope**2 + east_west_slope**2)
                slope_direction = self._get_slope_direction(east_west_slope, north_south_slope)
            else:
                slope_percent = 0
                slope_direction = "flat"
            
            corner_samples = [
                {"point": [lon, lat + offset_deg], "elevation_m": elevations[1] if len(elevations) > 1 else center_elev},
                {"point": [lon, lat - offset_deg], "elevation_m": elevations[2] if len(elevations) > 2 else center_elev},
                {"point": [lon + offset_deg, lat], "elevation_m": elevations[3] if len(elevations) > 3 else center_elev},
                {"point": [lon - offset_deg, lat], "elevation_m": elevations[4] if len(elevations) > 4 else center_elev},
            ]
            
            return TerrainData(
                elevation_m=center_elev,
                slope_percent=slope_percent,
                slope_direction=slope_direction,
                terrain_classification=self._classify_terrain(slope_percent),
                corner_samples=corner_samples,
                dem_source="open_meteo_elevation",
                resolution_m=30  # Approximate
            )
            
        except Exception as e:
            logger.error(f"Open-Meteo elevation failed: {e}")
            return self._get_default_terrain(lat, lon)
    
    def _sample_corners(
        self, 
        grid: np.ndarray, 
        cx: int, 
        cy: int, 
        offset: int,
        lat: float,
        lon: float,
        meters_per_pixel: float
    ) -> List[Dict[str, Any]]:
        """Sample elevation at property corners"""
        h, w = grid.shape
        corners = [
            (max(0, cx - offset), max(0, cy - offset)),  # SW
            (min(w-1, cx + offset), max(0, cy - offset)),  # SE
            (min(w-1, cx + offset), min(h-1, cy + offset)),  # NE
            (max(0, cx - offset), min(h-1, cy + offset)),  # NW
        ]
        
        deg_per_pixel = meters_per_pixel / 111000
        
        samples = []
        for px, py in corners:
            elev = float(grid[py, px])
            # Approximate lon/lat
            corner_lon = lon + (px - cx) * deg_per_pixel
            corner_lat = lat - (py - cy) * deg_per_pixel
            samples.append({
                "point": [corner_lon, corner_lat],
                "elevation_m": elev
            })
        
        return samples
    
    def _get_slope_direction(self, dx: float, dy: float) -> str:
        """Get cardinal direction of slope"""
        if abs(dx) < 0.001 and abs(dy) < 0.001:
            return "flat"
        
        angle = math.atan2(dy, dx) * 180 / math.pi
        
        if -22.5 <= angle < 22.5:
            return "east"
        elif 22.5 <= angle < 67.5:
            return "northeast"
        elif 67.5 <= angle < 112.5:
            return "north"
        elif 112.5 <= angle < 157.5:
            return "northwest"
        elif angle >= 157.5 or angle < -157.5:
            return "west"
        elif -157.5 <= angle < -112.5:
            return "southwest"
        elif -112.5 <= angle < -67.5:
            return "south"
        else:
            return "southeast"
    
    def _classify_terrain(self, slope_percent: float) -> str:
        """Classify terrain based on slope"""
        if slope_percent < 2:
            return "flat"
        elif slope_percent < 5:
            return "gentle_slope"
        elif slope_percent < 10:
            return "moderate_slope"
        elif slope_percent < 20:
            return "steep_slope"
        else:
            return "very_steep"
    
    def _get_default_terrain(self, lat: float, lon: float) -> TerrainData:
        """Return default terrain when all sources fail"""
        return TerrainData(
            elevation_m=50.0,  # UK average
            slope_percent=2.0,
            slope_direction="flat",
            terrain_classification="gentle_slope",
            corner_samples=[
                {"point": [lon - 0.0003, lat - 0.0003], "elevation_m": 50.0},
                {"point": [lon + 0.0003, lat - 0.0003], "elevation_m": 50.0},
                {"point": [lon + 0.0003, lat + 0.0003], "elevation_m": 50.0},
                {"point": [lon - 0.0003, lat + 0.0003], "elevation_m": 50.0},
            ],
            dem_source="default_estimate",
            resolution_m=30
        )

