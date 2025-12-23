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
from ..config import get_config


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
        search_radius_m: float = 30,
        osm_buildings: Optional[List[Dict[str, Any]]] = None,
        osm_roads: Optional[List[Dict[str, Any]]] = None,
        osm_landuse: Optional[List[Dict[str, Any]]] = None,
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
            osm_landuse: Optional OSM landuse data
            classify: Whether to classify into front/rear/sides
            
        Returns:
            PropertyLine object with classified segments
        """
        # Get property line
        property_line = self.property_line_processor.get_property_line(
            lat, lon, search_radius_m, osm_buildings, osm_roads, osm_landuse
        )
        
        if not property_line:
            return None
        
        # Classify if requested
        if classify and osm_roads and osm_buildings:
            property_line = self.classifier.classify(
                property_line, osm_roads, osm_buildings
            )
        
        return property_line
    
    def get_plot_boundary(
        self,
        lat: float,
        lon: float,
        search_radius_m: float = 30,
        osm_buildings: Optional[List[Dict[str, Any]]] = None,
        osm_roads: Optional[List[Dict[str, Any]]] = None,
        osm_landuse: Optional[List[Dict[str, Any]]] = None
    ) -> Dict[str, Any]:
        """
        Get or estimate plot boundary at given location (backward compatibility)
        
        This method maintains the old API while using the new modular components
        
        Returns:
            Dict with boundary data in old format
        """
        # Try to get property line first
        property_line = self.get_property_line(
            lat, lon, search_radius_m, osm_buildings, osm_roads, osm_landuse, classify=False
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
        
        # Strategy 1: Look for landuse polygon
        if osm_landuse:
            for landuse in osm_landuse:
                geom = landuse.get("geometry", {})
                if geom.get("type") == "Polygon":
                    coords = geom.get("coordinates", [[]])[0]
                    if coords and len(coords) >= 4:
                        if point_roughly_in_polygon(lon, lat, coords):
                            coords = close_polygon(coords)
                            area = calculate_polygon_area(coords, lat)
                            perimeter = calculate_perimeter(coords, lat)
                            
                            return {
                                "type": "Polygon",
                                "coordinates": [coords],
                                "area_sqm": area,
                                "perimeter_m": perimeter,
                                "source": "openstreetmap_landuse",
                                "accuracy_m": 2.0,
                                "osm_id": landuse.get("osm_id")
                            }
        else:
            landuse_boundary = self._find_landuse_boundary(lat, lon, search_radius_m)
            if landuse_boundary:
                return landuse_boundary
        
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
    
    def _find_landuse_boundary(
        self,
        lat: float,
        lon: float,
        radius_m: float
    ) -> Optional[Dict[str, Any]]:
        """Find landuse=residential polygon containing the point"""
        import requests
        
        logger.debug(f"Searching for landuse boundary at ({lat}, {lon}) within {radius_m}m")
        self._rate_limit()
        
        query = f"""
        [out:json][timeout:30];
        (
            way["landuse"="residential"](around:{radius_m},{lat},{lon});
            way["landuse"="allotments"](around:{radius_m},{lat},{lon});
            relation["landuse"="residential"](around:{radius_m},{lat},{lon});
        );
        out body;
        >;
        out skel qt;
        """
        
        try:
            response = requests.post(
                self.overpass_url,
                data={"data": query},
                headers={"User-Agent": self.config.api.user_agent},
                timeout=30
            )
            response.raise_for_status()
            data = response.json()
            
            # Parse nodes
            nodes = {}
            for element in data.get("elements", []):
                if element["type"] == "node":
                    nodes[element["id"]] = (element["lon"], element["lat"])
            
            # Find ways
            for element in data.get("elements", []):
                if element["type"] == "way":
                    way_nodes = element.get("nodes", [])
                    if len(way_nodes) >= 4:
                        coords = []
                        for node_id in way_nodes:
                            if node_id in nodes:
                                coords.append(list(nodes[node_id]))
                        
                        if len(coords) >= 4:
                            if point_roughly_in_polygon(lon, lat, coords):
                                coords = close_polygon(coords)
                                area = calculate_polygon_area(coords, lat)
                                perimeter = calculate_perimeter(coords, lat)
                                
                                return {
                                    "type": "Polygon",
                                    "coordinates": [coords],
                                    "area_sqm": area,
                                    "perimeter_m": perimeter,
                                    "source": "openstreetmap_landuse",
                                    "accuracy_m": 2.0,
                                    "osm_id": element["id"]
                                }
        
        except Exception as e:
            logger.warning(f"Landuse boundary search failed: {e}")
        
        return None
    
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
            b_coords = footprint.get("coordinates", [[]])[0]
            if not b_coords or len(b_coords) < 4:
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
            plot_coords = self._buffer_building_to_plot(closest_building["coords"], lat)
            if plot_coords:
                area = calculate_polygon_area(plot_coords, lat)
                perimeter = calculate_perimeter(plot_coords, lat)
                
                return {
                    "type": "Polygon",
                    "coordinates": [plot_coords],
                    "area_sqm": area,
                    "perimeter_m": perimeter,
                    "source": "estimated_from_building",
                    "accuracy_m": 5.0,
                    "building_osm_id": closest_building["id"]
                }
        
        return None
    
    def _buffer_building_to_plot(
        self,
        building_coords: List[List[float]],
        center_lat: float
    ) -> Optional[List[List[float]]]:
        """Create plot boundary by buffering building footprint"""
        if len(building_coords) < 4:
            return None
        
        if SHAPELY_AVAILABLE:
            try:
                coords = close_polygon(building_coords)
                building_poly = Polygon(coords)
                
                min_buffer_m = max(2.0, building_poly.length * 111000 * 0.05)
                buffered = building_poly.buffer(
                    min_buffer_m / 111000,
                    cap_style=3,
                    join_style=2
                )
                
                if isinstance(buffered, Polygon):
                    return [[c[0], c[1]] for c in buffered.exterior.coords]
            except Exception as e:
                logger.warning(f"Shapely buffering failed: {e}")
        
        # Fallback
        lat_per_m = 1 / 111000
        lon_per_m = 1 / (111000 * math.cos(math.radians(center_lat)))
        
        min_lon = min(c[0] for c in building_coords)
        max_lon = max(c[0] for c in building_coords)
        min_lat = min(c[1] for c in building_coords)
        max_lat = max(c[1] for c in building_coords)
        
        buffer_lon = 2.0 * lon_per_m
        buffer_lat = 5.0 * lat_per_m
        
        return [
            [min_lon - buffer_lon, min_lat - buffer_lat],
            [max_lon + buffer_lon, min_lat - buffer_lat],
            [max_lon + buffer_lon, max_lat + buffer_lat],
            [min_lon - buffer_lon, max_lat + buffer_lat],
            [min_lon - buffer_lon, min_lat - buffer_lat],
        ]
    
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

