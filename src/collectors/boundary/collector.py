"""
Main Boundary Collector

Orchestrates property line, setback line, and buildable envelope detection
"""

import time
import math
from typing import List, Dict, Any, Optional
from loguru import logger

try:
    from shapely.geometry import Polygon, Point, LineString
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - boundary detection will be limited")

from .models import PropertyLine, SetbackLine, BuildableEnvelope
from .inspire import InspireHandler
from .property_line import PropertyLineProcessor
from .classifier import PropertyLineClassifier
from .setback_line import SetbackLineProcessor
from .buildable_envelope import BuildableEnvelopeProcessor
from .utils import (
    calculate_polygon_area,
    calculate_perimeter,
    close_polygon,
    point_roughly_in_polygon,
    haversine_distance
)
from ...config import get_config


class BoundaryCollector:
    """
    Collect and estimate plot boundaries
    
    Uses modular components:
    - PropertyLineProcessor: Detects property boundaries
    - PropertyLineClassifier: Classifies into front/rear/sides
    - SetbackLineProcessor: Calculates setbacks (to be implemented)
    - BuildableEnvelopeProcessor: Calculates buildable area (to be implemented)
    """
    
    def __init__(self, inspire_gml_path: Optional[str] = None):
        self.config = get_config()
        self.overpass_url = self.config.api.overpass_url
        self._last_request_time = 0
        self._min_request_interval = 1.0
        
        # Typical UK residential plot dimensions (meters)
        self.default_plot_width = 12.0   # Front width
        self.default_plot_depth = 25.0   # Garden depth
        
        # Initialize processors
        self.inspire_handler = InspireHandler(inspire_gml_path)
        self.property_line_processor = PropertyLineProcessor(self.inspire_handler)
        self.classifier = PropertyLineClassifier()
        self.setback_processor = SetbackLineProcessor()
        self.buildable_processor = BuildableEnvelopeProcessor()
        
        # Expose inspire_gml_path for backward compatibility
        self.inspire_gml_path = self.inspire_handler.inspire_gml_path
    
    def _rate_limit(self):
        """Ensure we don't exceed rate limits"""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()
    
    def get_property_line(
        self,
        lat: float,
        lon: float,
        search_radius_m: float,
        osm_buildings: Optional[List[Dict[str, Any]]] = None,
        osm_roads: Optional[List[Dict[str, Any]]] = None,
        classify: bool = True
    ) -> Optional[PropertyLine]:
        """
        Get property line with optional classification
        
        Args:
            lat: Center latitude
            lon: Center longitude
            search_radius_m: Search radius in meters
            osm_buildings: Optional OSM building data
            osm_roads: Optional OSM road data
            classify: Whether to classify into front/rear/sides
            
        Returns:
            PropertyLine object with classified segments
        """
        # Get property line
        property_line = self.property_line_processor.get_property_line(
            lat, lon, search_radius_m, osm_buildings, osm_roads  # OSM landuse removed
        )
        
        if not property_line:
            return None
        
        # Classify if requested
        if classify and osm_roads and osm_buildings:
            property_line = self.classifier.classify(
                property_line, osm_roads, osm_buildings
            )
        
        return property_line
    
    def get_setback_line(
        self,
        property_line: PropertyLine
    ) -> Optional[SetbackLine]:
        """
        Get setback line from property line
        
        Args:
            property_line: PropertyLine object
            
        Returns:
            SetbackLine object or None
        """
        return self.setback_processor.calculate_setback_line(property_line)
    
    def get_buildable_envelope(
        self,
        property_line: PropertyLine,
        setback_line: Optional[SetbackLine] = None
    ) -> Optional[BuildableEnvelope]:
        """
        Get buildable envelope from property line and setback line
        
        Args:
            property_line: PropertyLine object
            setback_line: Optional SetbackLine (will be calculated if not provided)
            
        Returns:
            BuildableEnvelope object or None
        """
        if not setback_line:
            setback_line = self.get_setback_line(property_line)
            if not setback_line:
                logger.warning("Failed to calculate setback line")
                return None
        
        return self.buildable_processor.calculate_buildable_envelope(
            property_line, setback_line
        )
    
    def get_plot_boundary(
        self,
        lat: float,
        lon: float,
        search_radius_m: float,
        osm_buildings: Optional[List[Dict[str, Any]]] = None,
        osm_roads: Optional[List[Dict[str, Any]]] = None,
    ) -> Dict[str, Any]:
        """
        Get or estimate plot boundary at given location (backward compatibility)
        
        This method maintains the old API while using the new modular components
        
        Returns:
            Dict with boundary data in old format
        """
        # Try to get property line first
        property_line = self.get_property_line(
            lat, lon, search_radius_m, osm_buildings, osm_roads, None, classify=False  # OSM landuse removed
        )
        
        if property_line:
            # Convert to old format
            return {
                "type": "Polygon",
                "coordinates": [property_line.coordinates],
                "area_sqm": property_line.area_sqm,
                "perimeter_m": property_line.perimeter_m,
                "source": property_line.source,
                "accuracy_m": property_line.accuracy_m,
                "inspire_id": property_line.inspire_id
            }
        
        # Fallback to OSM-based detection (original logic)
        logger.info("Using OSM-based boundary detection (fallback)")
        
        # Strategy 2: Derive from roads
        if osm_roads:
            road_boundary = self._derive_from_roads(lat, lon, osm_roads, search_radius_m)
            if road_boundary:
                return road_boundary
        
        # Strategy 3: Derive from neighbor buildings
        if osm_buildings:
            neighbor_boundary = self._derive_from_neighbor_buildings(lat, lon, osm_buildings, search_radius_m)
            if neighbor_boundary:
                return neighbor_boundary
        
        # Strategy 4: Estimate from building
        building_boundary = self._estimate_from_building(lat, lon, search_radius_m, osm_buildings)
        if building_boundary:
            return building_boundary
        
        # Strategy 5: Default plot
        return self._create_default_plot(lat, lon)
    
    
    def _derive_from_roads(
        self,
        lat: float,
        lon: float,
        osm_roads: List[Dict[str, Any]],
        radius_m: float
    ) -> Optional[Dict[str, Any]]:
        """Derive property boundary from road network"""
        if not SHAPELY_AVAILABLE:
            return None
        
        try:
            point = Point(lon, lat)
            nearby_roads = []
            
            for road in osm_roads:
                centerline = road.get("centerline", {})
                if not centerline:
                    continue
                coords = centerline.get("coordinates", [])
                if len(coords) < 2:
                    continue
                
                road_line = LineString(coords)
                road_point = Point(coords[len(coords)//2])
                if point.distance(road_point) * 111000 < radius_m:
                    nearby_roads.append({
                        "line": road_line,
                        "width_m": road.get("width_m", 5.0),
                        "coords": coords
                    })
            
            if len(nearby_roads) < 2:
                return None
            
            buffered_roads = []
            for road in nearby_roads:
                buffer_m = (road["width_m"] / 2) + 1.5
                buffered = road["line"].buffer(
                    buffer_m / 111000,
                    cap_style=3,
                    join_style=2
                )
                buffered_roads.append(buffered)
            
            if len(buffered_roads) >= 2:
                union = unary_union(buffered_roads)
                if hasattr(union, 'geoms'):
                    for geom in union.geoms:
                        if isinstance(geom, Polygon) and geom.contains(point):
                            coords = [[c[0], c[1]] for c in geom.exterior.coords]
                            return self._create_boundary_from_coords(
                                coords, "derived_from_roads", 3.0
                            )
        
        except Exception as e:
            logger.warning(f"Failed to derive boundary from roads: {e}")
        
        return None
    
    def _derive_from_neighbor_buildings(
        self,
        lat: float,
        lon: float,
        osm_buildings: List[Dict[str, Any]],
        radius_m: float
    ) -> Optional[Dict[str, Any]]:
        """Derive property boundary from neighbor building boundaries"""
        if not SHAPELY_AVAILABLE:
            return None
        
        try:
            point = Point(lon, lat)
            nearby_buildings = []
            
            for building in osm_buildings:
                footprint = building.get("footprint", {})
                if not footprint:
                    continue
                coords = footprint.get("coordinates", [[]])[0]
                if not coords or len(coords) < 4:
                    continue
                
                building_poly = Polygon(coords)
                if point.distance(building_poly) * 111000 < radius_m:
                    nearby_buildings.append(building_poly)
            
            if len(nearby_buildings) < 1:
                return None
            
            buffered_buildings = []
            for building in nearby_buildings:
                buffered = building.buffer(
                    1.0 / 111000,
                    cap_style=3,
                    join_style=2
                )
                buffered_buildings.append(buffered)
            
            if len(buffered_buildings) >= 1:
                search_radius_deg = radius_m / 111000
                search_area = point.buffer(search_radius_deg)
                
                result = search_area
                for buffered in buffered_buildings:
                    result = result.difference(buffered)
                
                if hasattr(result, 'geoms'):
                    largest = None
                    largest_area = 0
                    for geom in result.geoms:
                        if isinstance(geom, Polygon):
                            area = geom.area
                            if area > largest_area and (geom.contains(point) or geom.distance(point) * 111000 < 10):
                                largest = geom
                                largest_area = area
                    
                    if largest:
                        coords = [[c[0], c[1]] for c in largest.exterior.coords]
                        return self._create_boundary_from_coords(
                            coords, "derived_from_neighbor_buildings", 5.0
                        )
        
        except Exception as e:
            logger.warning(f"Failed to derive boundary from neighbor buildings: {e}")
        
        return None
    
    def _estimate_from_building(
        self,
        lat: float,
        lon: float,
        radius_m: float,
        osm_buildings: Optional[List[Dict[str, Any]]] = None
    ) -> Optional[Dict[str, Any]]:
        """Estimate plot boundary from nearby building footprint"""
        if not osm_buildings:
            return None
        
        closest_building = None
        min_distance = float('inf')
        
        for building in osm_buildings:
            footprint = building.get("footprint", {})
            if not footprint:
                continue
            b_coords_raw = footprint.get("coordinates", [])
            if not b_coords_raw:
                logger.debug(f"Building {building.get('osm_id')} has no coordinates in footprint")
                continue
            
            # Extract coordinates - handle GeoJSON Polygon format [[[lon, lat], ...]]
            if isinstance(b_coords_raw[0], list) and len(b_coords_raw[0]) > 0:
                if isinstance(b_coords_raw[0][0], list):
                    # Already in format [[[lon, lat], ...]] - extract first ring
                    b_coords = b_coords_raw[0]
                else:
                    # Format [[lon, lat], ...] - use as is
                    b_coords = b_coords_raw[0]
            else:
                logger.debug(f"Building {building.get('osm_id')} has invalid coordinate format")
                continue
            
            if not b_coords or len(b_coords) < 4:
                logger.debug(f"Building {building.get('osm_id')} has insufficient coordinates: {len(b_coords) if b_coords else 0} points")
                continue
            
            centroid_lon = sum(c[0] for c in b_coords) / len(b_coords)
            centroid_lat = sum(c[1] for c in b_coords) / len(b_coords)
            
            dist = haversine_distance(lat, lon, centroid_lat, centroid_lon)
            if dist < min_distance:
                min_distance = dist
                closest_building = {
                    "coords": b_coords,
                    "tags": building.get("tags", {}),
                    "id": building.get("osm_id")
                }
        
        if closest_building and min_distance < radius_m:
            building_coords = closest_building["coords"]
            logger.debug(f"Building coords for buffering: {len(building_coords)} points")
            
            plot_coords = self._buffer_building_to_plot(building_coords, lat)
            if plot_coords:
                logger.debug(f"_buffer_building_to_plot returned: type={type(plot_coords)}, length={len(plot_coords) if isinstance(plot_coords, list) else 'N/A'}")
                
                # Validate plot_coords - must have at least 3 points for a polygon
                if not isinstance(plot_coords, list):
                    logger.error(f"Buffered plot coords is not a list: type={type(plot_coords)}")
                    logger.error(f"Building coords that failed: {len(building_coords)} points")
                    return None
                
                if len(plot_coords) < 3:
                    logger.error(f"Buffered plot coords has insufficient points: {len(plot_coords)} (expected >= 3)")
                    logger.error(f"Plot coords content: {plot_coords}")
                    logger.error(f"Building coords that failed: {len(building_coords)} points")
                    return None
                
                # Validate each coordinate is a list of 2 numbers
                for i, coord in enumerate(plot_coords):
                    if not isinstance(coord, list) or len(coord) != 2:
                        logger.error(f"Invalid coordinate at index {i}: {coord} (type={type(coord)})")
                        logger.error(f"Full plot_coords: {plot_coords}")
                        return None
                
                try:
                    area = calculate_polygon_area(plot_coords, lat)
                    perimeter = calculate_perimeter(plot_coords, lat)
                    
                    logger.info(f"Successfully estimated plot from building {closest_building['id']}: {len(plot_coords)} points, {area:.1f} m²")
                    
                    # Create result with proper GeoJSON format
                    # plot_coords is [[lon, lat], [lon, lat], ...] (list of coordinate pairs)
                    # GeoJSON Polygon format requires: [[[lon, lat], [lon, lat], ...]] (list containing one ring)
                    result = {
                        "type": "Polygon",
                        "coordinates": [plot_coords],  # GeoJSON format: [[[lon, lat], ...]]
                        "area_sqm": area,
                        "perimeter_m": perimeter,
                        "source": "estimated_from_building",
                        "accuracy_m": 5.0,
                        "building_osm_id": closest_building["id"]
                    }
                    
                    # Validate the result structure before returning
                    result_coords = result.get("coordinates", [])
                    if not result_coords or len(result_coords) == 0:
                        logger.error(f"Result coordinates is empty after creation!")
                        logger.error(f"plot_coords had {len(plot_coords)} points before wrapping")
                        return None
                    
                    if not isinstance(result_coords, list):
                        logger.error(f"Result coordinates is not a list: {type(result_coords)}")
                        return None
                    
                    if len(result_coords) == 0:
                        logger.error(f"Result coordinates list is empty")
                        return None
                    
                    first_ring = result_coords[0]
                    if not isinstance(first_ring, list):
                        logger.error(f"First ring is not a list: {type(first_ring)}")
                        return None
                    
                    if len(first_ring) < 3:
                        logger.error(f"Result first ring has insufficient points: {len(first_ring)} (expected >= 3)")
                        logger.error(f"First ring content: {first_ring[:3]}...")
                        return None
                    
                    logger.info(f"Returning property line with {len(first_ring)} points in first ring (source: estimated_from_building)")
                    return result
                    
                except Exception as e:
                    logger.error(f"Failed to calculate area/perimeter for plot coords: {e}")
                    logger.error(f"Plot coords: {plot_coords[:3]}... (showing first 3)")
                    import traceback
                    logger.debug(traceback.format_exc())
                    return None
            else:
                logger.warning(f"_buffer_building_to_plot returned None for building {closest_building.get('id', 'unknown')}")
        
        return None
    
    def _buffer_building_to_plot(
        self,
        building_coords: List[List[float]],
        center_lat: float
    ) -> Optional[List[List[float]]]:
        """Create plot boundary by buffering building footprint
        
        Uses a simple offset method instead of Shapely buffering for reliability.
        """
        if not building_coords or len(building_coords) < 4:
            logger.warning(f"Building coordinates invalid: {len(building_coords) if building_coords else 0} points")
            return None
        
        logger.debug(f"Buffering building with {len(building_coords)} points")
        
        # Use simple offset method - more reliable than Shapely buffering
        # Calculate bounding box and expand with buffer
        try:
            # Calculate meters per degree at this latitude
            lat_per_m = 1 / 111000
            lon_per_m = 1 / (111000 * math.cos(math.radians(center_lat)))
            
            # Find bounding box
            min_lon = min(c[0] for c in building_coords)
            max_lon = max(c[0] for c in building_coords)
            min_lat = min(c[1] for c in building_coords)
            max_lat = max(c[1] for c in building_coords)
            
            # Calculate buffer distance (2-5m depending on building size)
            building_width = (max_lon - min_lon) * 111000 * math.cos(math.radians(center_lat))
            building_depth = (max_lat - min_lat) * 111000
            building_perimeter = 2 * (building_width + building_depth)
            
            # Buffer: minimum 2m, or 5% of perimeter
            buffer_m = max(2.0, building_perimeter * 0.05)
            logger.debug(f"Calculated buffer: {buffer_m:.2f}m (building size: {building_width:.1f}m x {building_depth:.1f}m)")
            
            # Convert buffer to degrees
            buffer_lon = buffer_m * lon_per_m
            buffer_lat = buffer_m * lat_per_m
            
            # Create buffered rectangular plot
            plot_coords = [
                [min_lon - buffer_lon, min_lat - buffer_lat],
                [max_lon + buffer_lon, min_lat - buffer_lat],
                [max_lon + buffer_lon, max_lat + buffer_lat],
                [min_lon - buffer_lon, max_lat + buffer_lat],
                [min_lon - buffer_lon, min_lat - buffer_lat],  # Close polygon
            ]
            
            # Validate result - should always have 5 points (4 corners + closing)
            if len(plot_coords) < 3:
                logger.error(f"Generated plot coords only has {len(plot_coords)} points - this should never happen!")
                return None
            
            logger.debug(f"Successfully created buffered plot: {len(plot_coords)} points")
            return plot_coords
            
        except Exception as e:
            logger.error(f"Failed to buffer building coordinates: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return None
    
    def _create_default_plot(self, lat: float, lon: float) -> Dict[str, Any]:
        """Create default rectangular plot at location"""
        lat_per_m = 1 / 111000
        lon_per_m = 1 / (111000 * math.cos(math.radians(lat)))
        
        half_width = (self.default_plot_width / 2) * lon_per_m
        half_depth = (self.default_plot_depth / 2) * lat_per_m
        
        coords = [
            [lon - half_width, lat - half_depth],
            [lon + half_width, lat - half_depth],
            [lon + half_width, lat + half_depth],
            [lon - half_width, lat + half_depth],
            [lon - half_width, lat - half_depth],
        ]
        
        area = self.default_plot_width * self.default_plot_depth
        perimeter = 2 * (self.default_plot_width + self.default_plot_depth)
        
        return {
            "type": "Polygon",
            "coordinates": [coords],
            "area_sqm": area,
            "perimeter_m": perimeter,
            "source": "default_estimate",
            "accuracy_m": 10.0,
            "note": "Default rectangular plot - verify with cadastral data"
        }
    
    def _create_boundary_from_coords(
        self,
        coords: List[List[float]],
        source: str,
        accuracy_m: float
    ) -> Dict[str, Any]:
        """Create boundary dict from coordinates"""
        if not coords:
            return None
        
        coords = close_polygon(coords)
        
        # Validate coordinates
        if len(coords) < 3:
            logger.warning(f"Coordinates have too few points ({len(coords)}), cannot create boundary")
            return None
        
        center_lat = sum(c[1] for c in coords) / len(coords)
        area = calculate_polygon_area(coords, center_lat)
        perimeter = calculate_perimeter(coords, center_lat)
        
        return {
            "type": "Polygon",
            "coordinates": [coords],
            "area_sqm": area,
            "perimeter_m": perimeter,
            "source": source,
            "accuracy_m": accuracy_m
        }

