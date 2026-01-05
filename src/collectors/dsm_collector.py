"""
DSM (Digital Surface Model) Data Collector

Provides DSM data for calculating real building heights.
Building Height = DSM - DTM

Data Sources:
1. UK Environment Agency LiDAR (for UK locations)
2. Copernicus DEM (global coverage)
3. Mapbox (if available - note: Mapbox Terrain-RGB is DTM only)
"""

import os
import math
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import requests
from loguru import logger
import numpy as np

try:
    from dotenv import load_dotenv
    HAS_DOTENV = True
except ImportError:
    HAS_DOTENV = False

try:
    import mercantile
    HAS_MERCANTILE = True
except ImportError:
    HAS_MERCANTILE = False
    logger.warning("mercantile not installed - DSM tile processing limited")

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    logger.warning("PIL not installed - DSM tile processing limited")

try:
    from owslib.wcs import WebCoverageService
    OWSLIB_AVAILABLE = True
except ImportError:
    OWSLIB_AVAILABLE = False
    logger.warning("OWSLib not available. Install with: pip install owslib")

try:
    from pyproj import Transformer
    HAS_PYPROJ = True
except ImportError:
    HAS_PYPROJ = False
    logger.warning("pyproj not available - coordinate transformation limited")

try:
    from rasterio.io import MemoryFile
    import rasterio
    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False
    logger.warning("rasterio not available - DSM raster processing limited")

from ..config import get_config


class DSMCollector:
    """
    Collect DSM (Digital Surface Model) data from multiple sources
    
    DSM includes buildings and vegetation, unlike DTM which is ground only.
    Building Height = DSM - DTM
    """
    
    # UK Environment Agency LiDAR DSM (Web Coverage Service)
    # Note: UK EA provides both DTM and DSM via WCS
    UK_LIDAR_DSM_WCS = "https://environment.data.gov.uk/spatialdata/lidar-composite-dsm-2m/wcs"
    
    # Copernicus DEM (30m resolution, global coverage)
    # Note: Copernicus provides both DSM and DTM
    COPERNICUS_DSM_URL = "https://portal.opentopography.org/raster"
    
    # OpenTopography API (provides access to various DEM datasets)
    OPENTOPOGRAPHY_API = "https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster"
    
    def __init__(self):
        self.config = get_config()
        
        # Load .env file if available
        if HAS_DOTENV:
            env_paths = [
                Path(__file__).parent.parent.parent / ".env",
                Path.cwd() / ".env",
                Path.home() / ".env",
            ]
            
            loaded = False
            for env_path in env_paths:
                if env_path.exists():
                    load_dotenv(env_path, override=False)
                    logger.debug(f"Loaded .env file from {env_path}")
                    loaded = True
                    break
            
            if not loaded:
                load_dotenv()
    
    def get_dsm_elevation_at_point(
        self, 
        lat: float, 
        lon: float,
        country: Optional[str] = None
    ) -> Optional[float]:
        """
        Get DSM elevation at a specific point
        
        Args:
            lat: Latitude
            lon: Longitude
            country: Country code (e.g., 'GB' for UK) - helps select best data source
            
        Returns:
            DSM elevation in meters, or None if unavailable
        """
        # Try UK LiDAR first if in UK
        if country == "GB" or country is None:
            try:
                dsm_elev = self._get_uk_lidar_dsm(lat, lon)
                if dsm_elev is not None:
                    logger.debug(f"DSM elevation at ({lat:.6f}, {lon:.6f}) = {dsm_elev:.1f}m from UK LiDAR")
                    return dsm_elev
            except Exception as e:
                logger.debug(f"UK LiDAR DSM failed: {e}")
        
        # Try Copernicus DEM (global coverage)
        try:
            dsm_elev = self._get_copernicus_dsm(lat, lon)
            if dsm_elev is not None:
                logger.debug(f"DSM elevation at ({lat:.6f}, {lon:.6f}) = {dsm_elev:.1f}m from Copernicus DEM")
                return dsm_elev
        except Exception as e:
            logger.debug(f"Copernicus DSM failed: {e}")
        
        logger.warning(f"Could not get DSM elevation at ({lat:.6f}, {lon:.6f}) from any source")
        return None
    
    def _get_uk_lidar_dsm(self, lat: float, lon: float) -> Optional[float]:
        """
        Get DSM from UK Environment Agency LiDAR using WCS GetCoverage
        
        UK EA LiDAR provides 2m resolution DSM data via WCS
        """
        logger.debug(f"Attempting to get UK LiDAR DSM for ({lat:.6f}, {lon:.6f})")
        
        if not OWSLIB_AVAILABLE:
            logger.debug("OWSLib not available for UK LiDAR WCS access")
            return None
        
        if not HAS_PYPROJ:
            logger.debug("pyproj not available for coordinate transformation")
            return None
        
        try:
            # Transform WGS84 (lat/lon) to British National Grid (EPSG:27700)
            transformer = Transformer.from_crs("EPSG:4326", "EPSG:27700", always_xy=True)
            x, y = transformer.transform(lon, lat)
            
            # Create a small bounding box around the point (10m x 10m)
            buffer_m = 5.0
            bbox = (x - buffer_m, y - buffer_m, x + buffer_m, y + buffer_m)
            
            logger.debug(f"Querying UK LiDAR DSM WCS at BNG coordinates: {x:.1f}, {y:.1f}")
            
            # Connect to WCS
            wcs = WebCoverageService(self.UK_LIDAR_DSM_WCS, version='2.0.1')
            
            # Get available coverages
            coverages = list(wcs.contents.keys())
            if not coverages:
                logger.debug("No coverages found in UK LiDAR DSM WCS")
                return None
            
            # Use first available coverage
            coverage_id = coverages[0]
            logger.debug(f"Using coverage: {coverage_id}")
            
            # Get coverage for bounding box
            # Note: WCS 2.0 uses subset syntax
            response = wcs.getCoverage(
                identifier=coverage_id,
                bbox=bbox,
                crs='EPSG:27700',
                format='image/tiff',
                width=10,
                height=10
            )
            
            if HAS_RASTERIO:
                # Read the raster data
                with MemoryFile(response.read()) as memfile:
                    with memfile.open() as dataset:
                        data = dataset.read(1)
                        # Get center pixel value (DSM elevation)
                        center_y, center_x = data.shape[0] // 2, data.shape[1] // 2
                        dsm_elev = float(data[center_y, center_x])
                        
                        # Check for no-data values
                        if dataset.nodata is not None and dsm_elev == dataset.nodata:
                            logger.debug(f"No data value at point, DSM unavailable")
                            return None
                        
                        logger.info(f"UK LiDAR DSM elevation: {dsm_elev:.1f}m at ({lat:.6f}, {lon:.6f})")
                        return dsm_elev
            else:
                logger.debug("rasterio not available for DSM raster processing")
                return None
                
        except Exception as e:
            logger.debug(f"UK LiDAR DSM WCS query failed: {e}")
            return None
    
    def _get_copernicus_dsm(self, lat: float, lon: float) -> Optional[float]:
        """
        Get DSM from Copernicus DEM
        
        Note: Copernicus DEM is available via OpenTopography
        """
        logger.debug(f"Attempting to get Copernicus DSM for ({lat:.6f}, {lon:.6f})")
        
        # Copernicus DEM is available at 30m resolution globally
        # We can use OpenTopography API or download tiles
        
        # For now, this is a placeholder
        # In production, you'd:
        # 1. Determine which Copernicus tile covers the point
        # 2. Download or query the tile
        # 3. Extract elevation value
        
        # TODO: Implement Copernicus DEM access
        return None
    
    def calculate_building_height(
        self,
        dsm_elevation: float,
        dtm_elevation: float
    ) -> float:
        """
        Calculate building height from DSM and DTM
        
        Building Height = DSM - DTM
        
        Args:
            dsm_elevation: DSM elevation (includes building)
            dtm_elevation: DTM elevation (ground only)
            
        Returns:
            Building height in meters
        """
        height = dsm_elevation - dtm_elevation
        
        # Ensure height is non-negative (DSM should be >= DTM)
        if height < 0:
            logger.warning(f"Negative building height calculated: {height:.1f}m (DSM={dsm_elevation:.1f}m, DTM={dtm_elevation:.1f}m)")
            return 0.0
        
        logger.info(f"Calculated building height: {height:.1f}m (DSM={dsm_elevation:.1f}m - DTM={dtm_elevation:.1f}m)")
        return height

