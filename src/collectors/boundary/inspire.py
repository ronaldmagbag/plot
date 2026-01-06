"""
INSPIRE GML data handler

Handles loading and querying INSPIRE cadastral data
"""

from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from loguru import logger

try:
    import geopandas as gpd
    GEOPANDAS_AVAILABLE = True
except ImportError:
    GEOPANDAS_AVAILABLE = False
    logger.debug("Geopandas not available - INSPIRE GML support disabled")

try:
    from shapely.geometry import Point, Polygon
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - INSPIRE queries will be limited")


class InspireHandler:
    """Handles INSPIRE GML cadastral data"""
    
    def __init__(self, inspire_gml_path: Optional[str] = None):
        self.inspire_gml_path = inspire_gml_path
        self._inspire_gdf = None
        self._inspire_gdf_wgs84 = None
        self._loaded_gml_files = []  # Track which files were loaded
        
        if GEOPANDAS_AVAILABLE:
            if self.inspire_gml_path:
                # Load specific GML file
                self._load_inspire_gml(self.inspire_gml_path)
            else:
                # Load all GML files from default directory
                self._load_all_inspire_gml()
    
    def _find_all_inspire_gml(self) -> List[Path]:
        """Find all INSPIRE GML files in default data/inspires directory"""
        # Try to find the project root (assuming we're in src/collectors/boundary/)
        current_file = Path(__file__)
        project_root = current_file.parent.parent.parent.parent
        
        # Look for INSPIRE GML files in data/inspires/
        inspire_dir = project_root / "data" / "inspires"
        
        if inspire_dir.exists():
            # Look for all .gml files
            gml_files = sorted(inspire_dir.glob("*.gml"))
            if gml_files:
                logger.debug(f"Found {len(gml_files)} INSPIRE GML file(s): {[f.name for f in gml_files]}")
                return gml_files
        
        logger.debug("No INSPIRE GML files found in data/inspires/ - will use OSM data")
        return []
    
    def _load_all_inspire_gml(self):
        """Load all INSPIRE GML files from default directory and merge them"""
        if not GEOPANDAS_AVAILABLE:
            logger.warning("Geopandas not available - cannot load INSPIRE GML")
            return
        
        gml_files = self._find_all_inspire_gml()
        if not gml_files:
            return
        
        all_gdfs = []
        self._loaded_gml_files = []
        
        for gml_path in gml_files:
            try:
                logger.info(f"Loading INSPIRE GML file: {gml_path.name}")
                gdf = gpd.read_file(str(gml_path))
                
                # Ensure CRS is set (should be EPSG:27700 for UK INSPIRE data)
                if gdf.crs is None:
                    logger.warning(f"No CRS found in {gml_path.name}, assuming EPSG:27700")
                    gdf.set_crs("EPSG:27700", inplace=True)
                
                # Add source column to track which file each parcel came from
                gdf['_source_file'] = gml_path.name
                
                all_gdfs.append(gdf)
                self._loaded_gml_files.append(str(gml_path))
                logger.info(f"  Loaded {len(gdf)} parcels from {gml_path.name}")
                
            except Exception as e:
                logger.warning(f"Failed to load {gml_path.name}: {e}")
                continue
        
        if not all_gdfs:
            logger.warning("No GML files could be loaded")
            return
        
        # Merge all GeoDataFrames
        try:
            import pandas as pd
            self._inspire_gdf = gpd.GeoDataFrame(pd.concat(all_gdfs, ignore_index=True))
            
            # Ensure CRS is consistent
            if self._inspire_gdf.crs is None:
                self._inspire_gdf.set_crs("EPSG:27700", inplace=True)
            
            # Convert to WGS84 for spatial queries
            self._inspire_gdf_wgs84 = self._inspire_gdf.to_crs("EPSG:4326")
            
            logger.info(f"Successfully merged {len(self._loaded_gml_files)} GML file(s): {len(self._inspire_gdf)} total cadastral parcels")
            
        except Exception as e:
            logger.warning(f"Failed to merge GML files: {e}")
            self._inspire_gdf = None
            self._inspire_gdf_wgs84 = None
    
    def _load_inspire_gml(self, gml_path: str):
        """Load a specific INSPIRE GML file for cadastral boundary queries"""
        if not GEOPANDAS_AVAILABLE:
            logger.warning("Geopandas not available - cannot load INSPIRE GML")
            return
        
        try:
            gml_path_obj = Path(gml_path)
            if not gml_path_obj.exists():
                logger.warning(f"INSPIRE GML file not found: {gml_path}")
                return
            
            logger.info(f"Loading INSPIRE GML file: {gml_path_obj.name}")
            self._inspire_gdf = gpd.read_file(str(gml_path_obj))
            
            # Ensure CRS is set (should be EPSG:27700 for UK INSPIRE data)
            if self._inspire_gdf.crs is None:
                logger.warning("No CRS found in GML, assuming EPSG:27700")
                self._inspire_gdf.set_crs("EPSG:27700", inplace=True)
            
            # Add source column
            self._inspire_gdf['_source_file'] = gml_path_obj.name
            self._loaded_gml_files = [str(gml_path_obj)]
            
            # Convert to WGS84 for spatial queries
            self._inspire_gdf_wgs84 = self._inspire_gdf.to_crs("EPSG:4326")
            
            logger.info(f"Loaded {len(self._inspire_gdf)} INSPIRE cadastral parcels from {gml_path_obj.name}")
            
        except Exception as e:
            logger.warning(f"Failed to load INSPIRE GML file: {e}")
            self._inspire_gdf = None
            self._inspire_gdf_wgs84 = None
    
    def find_boundary(
        self,
        lat: float,
        lon: float
    ) -> Optional[Dict[str, Any]]:
        """
        Find boundary from INSPIRE cadastral data
        
        Returns:
            Dict with coordinates, area, perimeter, source, inspire_id, or None
        """
        if not SHAPELY_AVAILABLE:
            logger.debug("Shapely not available - cannot query INSPIRE GML data")
            return None
        
        if self._inspire_gdf_wgs84 is None:
            if not self._loaded_gml_files:
                logger.debug("No INSPIRE GML data loaded - no GML files found or failed to load")
            else:
                logger.debug(f"INSPIRE GML data loaded from {len(self._loaded_gml_files)} file(s) but GeoDataFrame is None")
            return None
        
        try:
            point = Point(lon, lat)
            
            # Find polygons containing the point
            matches = []
            for idx, row in self._inspire_gdf_wgs84.iterrows():
                geom = row.geometry
                if geom is None:
                    continue
                
                if geom.contains(point):
                    matches.append((idx, row, geom))
            
            if not matches:
                logger.debug(f"No INSPIRE parcel found containing point ({lat}, {lon}) - point is outside GML coverage")
                return None
            
            # If multiple matches, prefer the smallest polygon (most specific)
            if len(matches) > 1:
                matches.sort(key=lambda x: x[2].area)
            
            # Get the best match
            idx, row, geom = matches[0]
            
            # Convert polygon to coordinates
            if geom.geom_type == "Polygon":
                coords = [[coord[0], coord[1]] for coord in geom.exterior.coords]
                if coords[0] != coords[-1]:
                    coords.append(coords[0])
            elif geom.geom_type == "MultiPolygon":
                # Take the largest polygon
                largest = max(geom.geoms, key=lambda g: g.area)
                coords = [[coord[0], coord[1]] for coord in largest.exterior.coords]
                if coords[0] != coords[-1]:
                    coords.append(coords[0])
            else:
                return None
            
            # Get area in original CRS (EPSG:27700) for accurate area calculation
            original_geom = self._inspire_gdf.iloc[idx].geometry
            area_sqm = original_geom.area
            
            # Get INSPIRE ID if available
            inspire_id = None
            for col in ['inspire_id', 'INSPIREID', 'inspireId', 'id']:
                if col in row:
                    inspire_id = row[col]
                    break
            
            # Get source file name if available
            source_file = row.get('_source_file', 'unknown')
            
            return {
                "coordinates": coords,
                "area_sqm": area_sqm,
                "source": "inspire_cadastral",
                "accuracy_m": 2.0,  # INSPIRE cadastral data is accurate(2~3m)
                "inspire_id": inspire_id,
                "inspire_source_file": source_file
            }
            
        except Exception as e:
            logger.warning(f"Error querying INSPIRE data: {e}")
            return None
    
    def find_nearby_parcels(
        self,
        lat: float,
        lon: float,
        radius_m: float
    ) -> List[Dict[str, Any]]:
        """
        Find nearby parcels (for merging small properties)
        
        Returns:
            List of parcel dicts with coordinates, area, inspire_id
        """
        if self._inspire_gdf_wgs84 is None or not SHAPELY_AVAILABLE:
            return []
        
        try:
            point = Point(lon, lat)
            radius_deg = radius_m / 111000  # Rough conversion
            
            nearby_parcels = []
            for idx, row in self._inspire_gdf_wgs84.iterrows():
                geom = row.geometry
                if geom is None:
                    continue
                
                # Check if polygon is within radius
                if point.distance(geom) * 111000 <= radius_m:
                    if geom.geom_type == "Polygon":
                        coords = [[coord[0], coord[1]] for coord in geom.exterior.coords]
                        if coords[0] != coords[-1]:
                            coords.append(coords[0])
                    elif geom.geom_type == "MultiPolygon":
                        largest = max(geom.geoms, key=lambda g: g.area)
                        coords = [[coord[0], coord[1]] for coord in largest.exterior.coords]
                        if coords[0] != coords[-1]:
                            coords.append(coords[0])
                    else:
                        continue
                    
                    # Get area in original CRS
                    original_geom = self._inspire_gdf.iloc[idx].geometry
                    area_sqm = original_geom.area
                    
                    # Get INSPIRE ID
                    inspire_id = None
                    for col in ['inspire_id', 'INSPIREID', 'inspireId', 'id']:
                        if col in row:
                            inspire_id = row[col]
                            break
                    
                    nearby_parcels.append({
                        "coordinates": coords,
                        "area_sqm": area_sqm,
                        "inspire_id": inspire_id,
                        "distance_m": point.distance(geom) * 111000
                    })
            
            return nearby_parcels
            
        except Exception as e:
            logger.warning(f"Error finding nearby parcels: {e}")
            return []

