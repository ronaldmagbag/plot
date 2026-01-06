"""
Elevation data collector using Open-Elevation API and Open-Meteo API
Collects elevation data for plot corners and calculates slope
"""

import time
import math
from typing import List, Dict, Any, Tuple, Optional
import requests
from loguru import logger

from ..config import get_config


class ElevationCollector:
    """Collect elevation data from free APIs"""
    
    def __init__(self):
        self.config = get_config()
        self._last_request_time = 0
        self._min_request_interval = 1.0
        
        # Primary API: Open-Meteo (more reliable, free)
        self.open_meteo_url = "https://api.open-meteo.com/v1/elevation"
        
        # Fallback: Open-Elevation
        self.open_elevation_url = self.config.api.open_elevation_url
    
    def _rate_limit(self):
        """Ensure we don't exceed rate limits"""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()
    
    def get_elevation_single(self, lat: float, lon: float) -> Optional[float]:
        """
        Get elevation for a single point
        Returns elevation in meters or None if failed
        """
        logger.debug(f"Fetching elevation for single point ({lat:.6f}, {lon:.6f})")
        # Try Open-Meteo first (more reliable)
        elevation = self._get_elevation_open_meteo(lat, lon)
        if elevation is not None:
            logger.info(f"✓ Open-Meteo elevation at ({lat:.6f}, {lon:.6f}): {elevation}m")
            return elevation
        
        # Fallback to Open-Elevation
        logger.debug(f"Open-Meteo failed, trying Open-Elevation for ({lat:.6f}, {lon:.6f})")
        elevation = self._get_elevation_open_elevation(lat, lon)
        if elevation is not None:
            logger.info(f"✓ Open-Elevation at ({lat:.6f}, {lon:.6f}): {elevation}m")
        else:
            logger.warning(f"✗ Failed to get elevation for ({lat:.6f}, {lon:.6f}) from both APIs")
        return elevation
    
    def get_elevation_multi(self, points: List[Tuple[float, float]]) -> List[Optional[float]]:
        """
        Get elevation for multiple points
        points: List of (lat, lon) tuples
        Returns list of elevations in meters
        """
        if not points:
            return []
        
        logger.info(f"Fetching elevation for {len(points)} points using Open-Meteo batch API")
        # Try batch request with Open-Meteo
        elevations = self._get_elevation_open_meteo_batch(points)
        if elevations and any(e is not None for e in elevations):
            success_count = sum(1 for e in elevations if e is not None)
            logger.info(f"✓ Open-Meteo batch: {success_count}/{len(points)} points succeeded")
            if success_count < len(points):
                failed_indices = [i for i, e in enumerate(elevations) if e is None]
                logger.warning(f"  Failed points: {failed_indices}")
            return elevations
        
        # Fallback to individual requests
        logger.warning(f"Open-Meteo batch failed, falling back to individual requests (Open-Meteo → Open-Elevation)")
        return [self.get_elevation_single(lat, lon) for lat, lon in points]
    
    def _get_elevation_open_meteo(self, lat: float, lon: float) -> Optional[float]:
        """Get elevation using Open-Meteo API (single point)"""
        self._rate_limit()
        
        try:
            logger.debug(f"Requesting Open-Meteo API: ({lat:.6f}, {lon:.6f})")
            response = requests.get(
                self.open_meteo_url,
                params={"latitude": lat, "longitude": lon},
                headers={"User-Agent": self.config.api.user_agent},
                timeout=self.config.api.request_timeout
            )
            response.raise_for_status()
            data = response.json()
            
            elevation = data.get("elevation", [None])[0]
            if elevation is not None:
                logger.debug(f"Open-Meteo API response: {elevation}m")
                return float(elevation)
            else:
                logger.debug(f"Open-Meteo API returned no elevation data")
        except requests.exceptions.RequestException as e:
            logger.debug(f"Open-Meteo API request error: {type(e).__name__}: {e}")
        except Exception as e:
            logger.debug(f"Open-Meteo API unexpected error: {type(e).__name__}: {e}")
        
        return None
    
    def _get_elevation_open_meteo_batch(self, points: List[Tuple[float, float]]) -> List[Optional[float]]:
        """Get elevation using Open-Meteo API (batch)"""
        self._rate_limit()
        
        if not points:
            return []
        
        # Open-Meteo accepts comma-separated latitudes and longitudes
        lats = ",".join([str(p[0]) for p in points])
        lons = ",".join([str(p[1]) for p in points])
        
        try:
            logger.debug(f"Requesting Open-Meteo batch API for {len(points)} points")
            response = requests.get(
                self.open_meteo_url,
                params={"latitude": lats, "longitude": lons},
                headers={"User-Agent": self.config.api.user_agent},
                timeout=self.config.api.request_timeout
            )
            response.raise_for_status()
            data = response.json()
            
            elevations = data.get("elevation", [])
            if isinstance(elevations, list) and len(elevations) == len(points):
                elevations_float = [float(e) if e is not None else None for e in elevations]
                logger.debug(f"Open-Meteo batch API response: {elevations_float}")
                return elevations_float
            else:
                logger.warning(f"Open-Meteo batch API returned {len(elevations)} elevations, expected {len(points)}")
        except requests.exceptions.RequestException as e:
            logger.debug(f"Open-Meteo batch API request error: {type(e).__name__}: {e}")
        except Exception as e:
            logger.debug(f"Open-Meteo batch API unexpected error: {type(e).__name__}: {e}")
        
        return [None] * len(points)
    
    def _get_elevation_open_elevation(self, lat: float, lon: float) -> Optional[float]:
        """Get elevation using Open-Elevation API (fallback)"""
        self._rate_limit()
        
        try:
            logger.debug(f"Requesting Open-Elevation API: ({lat:.6f}, {lon:.6f})")
            # Open-Elevation API format
            response = requests.get(
                self.open_elevation_url,
                params={"locations": f"{lat},{lon}"},
                headers={"User-Agent": self.config.api.user_agent},
                timeout=self.config.api.request_timeout
            )
            response.raise_for_status()
            data = response.json()
            
            results = data.get("results", [])
            if results and "elevation" in results[0]:
                elevation = results[0]["elevation"]
                logger.debug(f"Open-Elevation API response: {elevation}m")
                return float(elevation)
            else:
                logger.debug(f"Open-Elevation API returned no elevation data")
        except requests.exceptions.RequestException as e:
            logger.debug(f"Open-Elevation API request error: {type(e).__name__}: {e}")
        except Exception as e:
            logger.debug(f"Open-Elevation API unexpected error: {type(e).__name__}: {e}")
        
        return None
    
    def get_elevation_map(
        self,
        center_lat: float,
        center_lon: float,
        width_m: float,
        height_m: float
    ) -> Dict[str, Any]:
        """
        Get elevation map for a rectangular area
        
        Returns elevation data for corners and calculated slope info
        """
        # Convert meters to approximate degrees
        # At ~51° latitude (UK): 1 degree lat ≈ 111km, 1 degree lon ≈ 69km
        lat_per_m = 1 / 111000
        lon_per_m = 1 / (111000 * math.cos(math.radians(center_lat)))
        
        half_width = (width_m / 2) * lon_per_m
        half_height = (height_m / 2) * lat_per_m
        
        # Define corner points (SW, SE, NE, NW)
        corners = [
            (center_lat - half_height, center_lon - half_width),  # SW
            (center_lat - half_height, center_lon + half_width),  # SE
            (center_lat + half_height, center_lon + half_width),  # NE
            (center_lat + half_height, center_lon - half_width),  # NW
        ]
        
        # Also get center elevation
        all_points = corners + [(center_lat, center_lon)]
        
        logger.info(f"Fetching elevation for {len(all_points)} points around ({center_lat:.6f}, {center_lon:.6f})")
        logger.info(f"  Area: {width_m:.1f}m × {height_m:.1f}m")
        elevations = self.get_elevation_multi(all_points)
        
        # Process results
        corner_samples = []
        corner_labels = ["SW", "SE", "NE", "NW"]
        
        logger.info("Elevation results:")
        for i, (corner, elev) in enumerate(zip(corners, elevations[:4])):
            elevation_value = elev if elev is not None else 0.0
            status = "✓" if elev is not None else "✗ (using 0.0)"
            logger.info(f"  {corner_labels[i]}: {elevation_value}m {status}")
            corner_samples.append({
                "point": [corner[1], corner[0]],  # [lon, lat]
                "elevation_m": elevation_value,
                "label": corner_labels[i]
            })
        
        center_elevation = elevations[4] if len(elevations) > 4 else None
        if center_elevation is not None:
            logger.info(f"  Center: {center_elevation}m ✓")
        else:
            logger.warning(f"  Center: failed to get elevation")
        
        # Calculate slope and direction
        slope_info = self._calculate_slope(corner_samples, width_m, height_m)
        
        return {
            "corner_samples": corner_samples,
            "slope_percent": slope_info["slope_percent"],
            "slope_direction": slope_info["slope_direction"],
            "average_elevation_m": slope_info["average_elevation"],
            "terrain_classification": self._classify_terrain(slope_info["slope_percent"]),
            "dem_reference": {
                "source": "open-meteo",
                "resolution_m": 90.0,  # SRTM-based
                "file_path": None
            }
        }
    
    def _calculate_slope(
        self,
        corner_samples: List[Dict[str, Any]],
        width_m: float,
        height_m: float
    ) -> Dict[str, Any]:
        """Calculate slope percentage and direction from corner elevations"""
        
        # Get elevations (SW, SE, NE, NW)
        elevations = [s["elevation_m"] for s in corner_samples]
        
        if not all(e is not None for e in elevations):
            return {
                "slope_percent": 0.0,
                "slope_direction": "flat",
                "average_elevation": sum(e for e in elevations if e) / max(1, len([e for e in elevations if e]))
            }
        
        sw, se, ne, nw = elevations
        
        # Calculate E-W slope (average of north and south edges)
        ew_slope_south = (se - sw) / width_m * 100  # positive = rising to east
        ew_slope_north = (ne - nw) / width_m * 100
        ew_slope = (ew_slope_south + ew_slope_north) / 2
        
        # Calculate N-S slope (average of east and west edges)
        ns_slope_west = (nw - sw) / height_m * 100  # positive = rising to north
        ns_slope_east = (ne - se) / height_m * 100
        ns_slope = (ns_slope_west + ns_slope_east) / 2
        
        # Overall slope magnitude
        overall_slope = math.sqrt(ew_slope**2 + ns_slope**2)
        
        # Determine primary direction
        direction = self._get_slope_direction(ew_slope, ns_slope)
        
        average_elevation = sum(elevations) / 4
        
        return {
            "slope_percent": round(overall_slope, 2),
            "slope_direction": direction,
            "average_elevation": round(average_elevation, 2)
        }
    
    def _get_slope_direction(self, ew_slope: float, ns_slope: float) -> str:
        """Determine slope direction from E-W and N-S components"""
        if abs(ew_slope) < 0.5 and abs(ns_slope) < 0.5:
            return "flat"
        
        # Calculate angle
        angle = math.degrees(math.atan2(ns_slope, ew_slope))
        
        # Map angle to compass direction (rising towards)
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
        """Classify terrain based on slope percentage"""
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

