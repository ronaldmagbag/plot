"""
Property line detection and processing

Handles property line detection from INSPIRE and OSM data,
validation (size checks), and merging of small adjacent properties
"""

import math
from typing import List, Dict, Any, Optional
from loguru import logger

try:
    from shapely.geometry import Polygon, Point
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - property line processing will be limited")

from .models import PropertyLine
from .inspire import InspireHandler
from .utils import (
    calculate_polygon_area,
    calculate_perimeter,
    close_polygon,
    haversine_distance
)
from ...config import get_config


class PropertyLineProcessor:
    """Processes property line detection, validation, and merging"""
    
    def __init__(self, inspire_handler: Optional[InspireHandler] = None):
        self.config = get_config()
        self.inspire_handler = inspire_handler or InspireHandler()
        
        # Get thresholds from config
        self.max_area_sqm = self.config.uk_regulatory.property_line_max_area_sqm
        self.max_perimeter_m = self.config.uk_regulatory.property_line_max_perimeter_m
        self.min_side_distance_m = self.config.uk_regulatory.property_line_min_side_distance_m
    
    def get_property_line(
        self,
        lat: float,
        lon: float,
        search_radius_m: float,
        osm_buildings: Optional[List[Dict[str, Any]]] = None,
        osm_roads: Optional[List[Dict[str, Any]]] = None,
        osm_landuse: Optional[List[Dict[str, Any]]] = None
    ) -> Optional[PropertyLine]:
        """
        Get property line at given location
        
        Priority:
        1. INSPIRE GML data (if available)
        2. Validate size - if too big, use OSM logic
        3. If too small, merge with adjacent properties
        4. Fallback to OSM-based detection
        
        Returns:
            PropertyLine object or None
        """
        logger.info(f"Detecting property line at ({lat}, {lon})")
        
        # Try INSPIRE first
        logger.info(f"Attempting to find property line from INSPIRE GML data...")
        inspire_data = self.inspire_handler.find_boundary(lat, lon)
        
        if inspire_data:
            logger.info(f"Successfully found property line from INSPIRE GML data (area: {inspire_data.get('area_sqm', 0):.1f} m²)")
            # Check if property is too big
            if (inspire_data["area_sqm"] > self.max_area_sqm or
                inspire_data.get("perimeter_m", 0) > self.max_perimeter_m):
                logger.info(f"INSPIRE property too large ({inspire_data['area_sqm']:.1f} sqm), using OSM logic")
                # Fall through to OSM logic
            else:
                # Check if property is too small and should be merged
                property_line = self._create_property_line_from_data(inspire_data, lat)
                
                # Check if too small
                if self._is_too_small(property_line):
                    logger.info("Property too small, attempting to merge with adjacent properties")
                    merged = self._merge_adjacent_properties(
                        property_line,
                        lat,
                        lon,
                        search_radius_m,
                        osm_buildings
                    )
                    if merged:
                        return merged
                
                return property_line
        else:
            # INSPIRE failed - log the reason
            if self.inspire_handler._inspire_gdf_wgs84 is None:
                if not self.inspire_handler._loaded_gml_files:
                    logger.info("No INSPIRE GML files loaded - no GML files found in data/inspires/ directory")
                else:
                    logger.info(f"INSPIRE GML files loaded but no parcel found at ({lat}, {lon}) - point is outside GML coverage area")
            else:
                logger.info(f"INSPIRE GML query failed - no parcel found containing point ({lat}, {lon})")
        
        # Fallback to OSM-based detection
        logger.info("Using OSM-based property line detection")
        osm_data = self._get_osm_property_line(
            lat, lon, search_radius_m, osm_buildings, osm_roads, osm_landuse
        )
        
        if osm_data:
            return self._create_property_line_from_data(osm_data, lat)
        
        return None
    
    def _create_property_line_from_data(
        self,
        data: Dict[str, Any],
        center_lat: float
    ) -> PropertyLine:
        """Create PropertyLine object from data dict"""
        coords = close_polygon(data["coordinates"])
        area = data.get("area_sqm", calculate_polygon_area(coords, center_lat))
        perimeter = data.get("perimeter_m", calculate_perimeter(coords, center_lat))
        
        return PropertyLine(
            coordinates=coords,
            area_sqm=area,
            perimeter_m=perimeter,
            source=data.get("source", "unknown"),
            accuracy_m=data.get("accuracy_m", 5.0),
            inspire_id=data.get("inspire_id"),
            metadata=data.get("metadata", {})
        )
    
    def _is_too_small(self, property_line: PropertyLine) -> bool:
        """Check if property line is too small (needs merging)"""
        if not SHAPELY_AVAILABLE:
            return False
        
        try:
            # Check minimum side distance
            prop_poly = Polygon(property_line.coordinates)
            
            # Get bounding box dimensions
            bounds = prop_poly.bounds
            width_deg = bounds[2] - bounds[0]  # lon difference
            height_deg = bounds[3] - bounds[1]  # lat difference
            
            # Convert to meters (rough approximation)
            center_lat = (bounds[1] + bounds[3]) / 2
            width_m = width_deg * 111000 * abs(math.cos(math.radians(center_lat)))
            height_m = height_deg * 111000
            
            min_dimension = min(width_m, height_m)
            
            if min_dimension < self.min_side_distance_m:
                logger.debug(f"Property too small: min dimension = {min_dimension:.1f}m < {self.min_side_distance_m}m")
                return True
            
            return False
            
        except Exception as e:
            logger.warning(f"Error checking property size: {e}")
            return False
    
    def _merge_adjacent_properties(
        self,
        property_line: PropertyLine,
        lat: float,
        lon: float,
        search_radius_m: float,
        osm_buildings: Optional[List[Dict[str, Any]]]
    ) -> Optional[PropertyLine]:
        """
        Merge adjacent small properties
        
        Checks if buildings span multiple property zones before merging
        """
        if not SHAPELY_AVAILABLE:
            return None
        
        try:
            # Find nearby parcels from INSPIRE
            nearby_parcels = self.inspire_handler.find_nearby_parcels(
                lat, lon, search_radius_m
            )
            
            if not nearby_parcels:
                return None
            
            # Check if any buildings span multiple properties
            if osm_buildings:
                for building in osm_buildings:
                    footprint = building.get("footprint", {})
                    if not footprint:
                        continue
                    b_coords = footprint.get("coordinates", [[]])[0]
                    if not b_coords or len(b_coords) < 4:
                        continue
                    
                    building_poly = Polygon(b_coords)
                    
                    # Check if building overlaps with current property
                    current_poly = Polygon(property_line.coordinates)
                    if not building_poly.intersects(current_poly):
                        continue
                    
                    # Check if building also overlaps with nearby parcels
                    for parcel in nearby_parcels:
                        parcel_poly = Polygon(parcel["coordinates"])
                        if building_poly.intersects(parcel_poly):
                            # Building spans multiple properties - don't merge
                            logger.info("Building spans multiple properties - skipping merge")
                            return None
            
            # Merge with closest small parcel
            current_poly = Polygon(property_line.coordinates)
            best_merge = None
            min_distance = float('inf')
            
            for parcel in nearby_parcels:
                # Check if parcel is also small
                parcel_area = parcel["area_sqm"]
                if parcel_area > self.max_area_sqm:
                    continue
                
                parcel_poly = Polygon(parcel["coordinates"])
                distance = current_poly.distance(parcel_poly) * 111000
                
                # Only merge if very close (within 2m)
                if distance < 2.0 and distance < min_distance:
                    min_distance = distance
                    best_merge = parcel
            
            if best_merge:
                # Merge polygons
                merged_poly = unary_union([current_poly, Polygon(best_merge["coordinates"])])
                
                if isinstance(merged_poly, Polygon) and merged_poly.is_valid:
                    merged_coords = [[c[0], c[1]] for c in merged_poly.exterior.coords]
                    merged_coords = close_polygon(merged_coords)
                    
                    center_lat = sum(c[1] for c in merged_coords) / len(merged_coords)
                    area = calculate_polygon_area(merged_coords, center_lat)
                    perimeter = calculate_perimeter(merged_coords, center_lat)
                    
                    logger.info(f"Merged property: {property_line.area_sqm:.1f} + {best_merge['area_sqm']:.1f} = {area:.1f} sqm")
                    
                    return PropertyLine(
                        coordinates=merged_coords,
                        area_sqm=area,
                        perimeter_m=perimeter,
                        source="inspire_cadastral_merged",
                        accuracy_m=1.0,
                        inspire_id=property_line.inspire_id,
                        metadata={
                            "merged": True,
                            "original_area_sqm": property_line.area_sqm,
                            "merged_parcel_area_sqm": best_merge["area_sqm"]
                        }
                    )
            
            return None
            
        except Exception as e:
            logger.warning(f"Failed to merge adjacent properties: {e}")
            return None
    
    def _get_osm_property_line(
        self,
        lat: float,
        lon: float,
        search_radius_m: float,
        osm_buildings: Optional[List[Dict[str, Any]]],
        osm_roads: Optional[List[Dict[str, Any]]],
        osm_landuse: Optional[List[Dict[str, Any]]]
    ) -> Optional[Dict[str, Any]]:
        """
        Get property line from OSM data (fallback when INSPIRE not available or too large)
        
        This uses the existing OSM-based boundary detection logic
        Returns None to let the main collector handle OSM fallback
        """
        # Return None - the main collector will handle OSM fallback
        # This allows the collector to use its comprehensive OSM detection strategies
        return None

