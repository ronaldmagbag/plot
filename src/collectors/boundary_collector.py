"""
Boundary collector for detecting plot boundaries
Uses OSM landuse data, building footprints, roads, and neighbor parcels to estimate plot boundaries
Supports arbitrary polygon shapes (not just rectangles)
"""

import time
import math
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
import requests
from loguru import logger

try:
    from shapely.geometry import Polygon, Point, LineString
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - polygon buffering will be limited")

try:
    import geopandas as gpd
    GEOPANDAS_AVAILABLE = True
except ImportError:
    GEOPANDAS_AVAILABLE = False
    logger.debug("Geopandas not available - INSPIRE GML support disabled")

from ..config import get_config


class BoundaryCollector:
    """
    Collect and estimate plot boundaries
    
    Uses multiple strategies (in priority order):
    0. INSPIRE cadastral data (automatically loaded from data/inspires/ if available)
    1. OSM landuse=residential polygons
    2. Road-derived boundaries
    3. Neighbor building boundaries
    4. Building footprint with buffer
    5. Estimated rectangular plot based on typical UK plot sizes
    """
    
    def __init__(self, inspire_gml_path: Optional[str] = None):
        self.config = get_config()
        self.overpass_url = self.config.api.overpass_url
        self._last_request_time = 0
        self._min_request_interval = 1.0
        
        # Typical UK residential plot dimensions (meters)
        self.default_plot_width = 12.0   # Front width
        self.default_plot_depth = 25.0   # Garden depth
        
        # INSPIRE GML data (default: try to find in data/inspires directory)
        if inspire_gml_path is None:
            # Try to find INSPIRE GML files in default location
            inspire_gml_path = self._find_default_inspire_gml()
        
        self.inspire_gml_path = inspire_gml_path
        self._inspire_gdf = None
        self._inspire_gdf_wgs84 = None
        if inspire_gml_path and GEOPANDAS_AVAILABLE:
            self._load_inspire_gml()
    
    def _find_default_inspire_gml(self) -> Optional[str]:
        """Find INSPIRE GML file in default data/inspires directory"""
        # Try to find the project root (assuming we're in src/collectors/)
        current_file = Path(__file__)
        project_root = current_file.parent.parent.parent
        
        # Look for INSPIRE GML files in data/inspires/
        inspire_dir = project_root / "data" / "inspires"
        
        if inspire_dir.exists():
            # Look for .gml files
            gml_files = list(inspire_dir.glob("*.gml"))
            if gml_files:
                # Use the first .gml file found
                default_gml = str(gml_files[0])
                logger.debug(f"Found default INSPIRE GML file: {default_gml}")
                return default_gml
        
        logger.debug("No INSPIRE GML file found in data/inspires/ - will use OSM data")
        return None
    
    def _rate_limit(self):
        """Ensure we don't exceed rate limits"""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()
    
    def _load_inspire_gml(self):
        """Load INSPIRE GML file for cadastral boundary queries"""
        if not GEOPANDAS_AVAILABLE:
            logger.warning("Geopandas not available - cannot load INSPIRE GML")
            return
        
        try:
            gml_path = Path(self.inspire_gml_path)
            if not gml_path.exists():
                logger.warning(f"INSPIRE GML file not found: {gml_path}")
                return
            
            logger.info(f"Loading INSPIRE GML file: {gml_path}")
            self._inspire_gdf = gpd.read_file(str(gml_path))
            
            # Ensure CRS is set (should be EPSG:27700 for UK INSPIRE data)
            if self._inspire_gdf.crs is None:
                logger.warning("No CRS found in GML, assuming EPSG:27700")
                self._inspire_gdf.set_crs("EPSG:27700", inplace=True)
            
            # Convert to WGS84 for spatial queries
            self._inspire_gdf_wgs84 = self._inspire_gdf.to_crs("EPSG:4326")
            
            logger.info(f"Loaded {len(self._inspire_gdf)} INSPIRE cadastral parcels")
            
        except Exception as e:
            logger.warning(f"Failed to load INSPIRE GML file: {e}")
            self._inspire_gdf = None
            self._inspire_gdf_wgs84 = None
    
    def _find_inspire_boundary(
        self,
        lat: float,
        lon: float
    ) -> Optional[Dict[str, Any]]:
        """Find boundary from INSPIRE cadastral data"""
        if self._inspire_gdf_wgs84 is None:
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
            
            # Calculate perimeter
            perimeter = self._calculate_perimeter(coords, lat)
            
            # Get INSPIRE ID if available
            inspire_id = None
            for col in ['inspire_id', 'INSPIREID', 'inspireId', 'id']:
                if col in row:
                    inspire_id = row[col]
                    break
            
            return {
                "type": "Polygon",
                "coordinates": [coords],
                "area_sqm": area_sqm,
                "perimeter_m": perimeter,
                "source": "inspire_cadastral",
                "accuracy_m": 0.5,  # INSPIRE cadastral data is very accurate
                "inspire_id": inspire_id
            }
            
        except Exception as e:
            logger.warning(f"Error querying INSPIRE data: {e}")
            return None
    
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
        Get or estimate plot boundary at given location
        
        Priority order:
        0. INSPIRE cadastral data (official UK land registry - highest accuracy)
        1. OSM landuse=residential polygons (actual parcel boundaries)
        2. Derive from roads (property line follows road centerlines)
        3. Derive from neighbor building boundaries
        4. Building footprint with proper polygon buffer (not rectangle)
        5. Default estimate (last resort)
        
        Returns property boundary polygon (any shape, not just rectangle) and metadata
        """
        logger.info(f"Detecting plot boundary at ({lat}, {lon})")
        
        # Strategy 0: Check INSPIRE cadastral data (highest priority - official UK land registry)
        if self._inspire_gdf_wgs84 is not None:
            logger.debug("Strategy 0: Checking INSPIRE cadastral data...")
            inspire_boundary = self._find_inspire_boundary(lat, lon)
            if inspire_boundary:
                logger.info("Found boundary from INSPIRE cadastral data (official UK land registry)")
                return inspire_boundary
        logger.debug("Strategy 0: No INSPIRE boundary found")
        
        # Strategy 1: Look for landuse polygon containing the point (best - actual parcel)
        logger.debug("Strategy 1: Checking for landuse polygon...")
        # Use provided landuse data if available (from batch query), otherwise query API
        if osm_landuse:
            # Check provided landuse polygons
            for landuse in osm_landuse:
                geom = landuse.get("geometry", {})
                if geom.get("type") == "Polygon":
                    coords = geom.get("coordinates", [[]])[0]
                    if coords and len(coords) >= 4:
                        # Check if point is in polygon
                        if self._point_roughly_in_polygon(lon, lat, coords):
                            # Close polygon if not closed
                            if coords[0] != coords[-1]:
                                coords.append(coords[0])
                            
                            area = self._calculate_polygon_area(coords, lat)
                            perimeter = self._calculate_perimeter(coords, lat)
                            
                            logger.info("Found landuse polygon for boundary (actual parcel)")
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
            # Fallback to API query if no landuse data provided
            landuse_boundary = self._find_landuse_boundary(lat, lon, search_radius_m)
            if landuse_boundary:
                logger.info("Found landuse polygon for boundary (actual parcel)")
                return landuse_boundary
        logger.debug("Strategy 1: No landuse polygon found")
        
        # Strategy 2: Derive from roads (property line follows road boundaries)
        logger.debug("Strategy 2: Deriving from roads...")
        if osm_roads:
            road_boundary = self._derive_from_roads(lat, lon, osm_roads, search_radius_m)
            if road_boundary:
                logger.info("Derived boundary from road network")
                return road_boundary
        logger.debug("Strategy 2: No road-derived boundary found")
        
        # Strategy 3: Derive from neighbor building boundaries
        logger.debug("Strategy 3: Deriving from neighbor buildings...")
        if osm_buildings:
            neighbor_boundary = self._derive_from_neighbor_buildings(lat, lon, osm_buildings, search_radius_m)
            if neighbor_boundary:
                logger.info("Derived boundary from neighbor building boundaries")
                return neighbor_boundary
        logger.debug("Strategy 3: No neighbor-derived boundary found")
        
        # Strategy 4: Find building and create plot around it (proper polygon buffer, not rectangle)
        logger.debug("Strategy 4: Estimating from building footprint...")
        try:
            building_boundary = self._estimate_from_building(lat, lon, search_radius_m, osm_buildings)
            if building_boundary:
                logger.info("Estimated boundary from building footprint with polygon buffer")
                return building_boundary
        except Exception as e:
            logger.warning(f"Strategy 4 failed: {e}")
        logger.debug("Strategy 4: No building-based boundary found")
        
        # Strategy 5: Create default plot (last resort - still rectangular but noted as estimate)
        logger.info("Using default plot estimate (verify with cadastral data)")
        return self._create_default_plot(lat, lon)
    
    def _find_landuse_boundary(
        self,
        lat: float,
        lon: float,
        radius_m: float
    ) -> Optional[Dict[str, Any]]:
        """Find landuse=residential polygon containing the point"""
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
            logger.debug(f"Querying Overpass API for landuse data...")
            try:
                response = requests.post(
                    self.overpass_url,
                    data={"data": query},
                    headers={"User-Agent": self.config.api.user_agent},
                    timeout=30
                )
                logger.debug(f"Overpass API response status: {response.status_code}")
                response.raise_for_status()
                data = response.json()
            except requests.exceptions.Timeout:
                logger.warning("Landuse boundary query timed out after 30s")
                return None
            except requests.exceptions.RequestException as e:
                logger.warning(f"Landuse boundary query failed: {e}")
                return None
            
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
                            # Check if point is roughly within this polygon
                            if self._point_roughly_in_polygon(lon, lat, coords):
                                # Close polygon if not closed
                                if coords[0] != coords[-1]:
                                    coords.append(coords[0])
                                
                                area = self._calculate_polygon_area(coords, lat)
                                perimeter = self._calculate_perimeter(coords, lat)
                                
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
            # This is a fallback strategy, so failure is acceptable
    
    def _estimate_from_building(
        self,
        lat: float,
        lon: float,
        radius_m: float,
        osm_buildings: Optional[List[Dict[str, Any]]] = None
    ) -> Optional[Dict[str, Any]]:
        """Estimate plot boundary from nearby building footprint"""
        
        # Use provided buildings if available, otherwise query OSM
        if osm_buildings:
            buildings_data = osm_buildings
        else:
            self._rate_limit()
            query = f"""
            [out:json][timeout:30];
            way["building"](around:{radius_m},{lat},{lon});
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
                
                buildings_data = []
                for element in data.get("elements", []):
                    if element["type"] == "way":
                        way_nodes = element.get("nodes", [])
                        if len(way_nodes) >= 4:
                            coords = []
                            for node_id in way_nodes:
                                if node_id in nodes:
                                    coords.append(list(nodes[node_id]))
                            
                            if len(coords) >= 4:
                                buildings_data.append({
                                    "coords": coords,
                                    "tags": element.get("tags", {}),
                                    "id": element["id"]
                                })
            except Exception as e:
                logger.warning(f"Failed to query buildings: {e}")
                # This is a fallback strategy, so failure is acceptable
                return None
        
        # Find closest building
        closest_building = None
        min_distance = float('inf')
        
        for building in buildings_data:
            if isinstance(building, dict) and "footprint" in building:
                # OSM building format
                footprint = building.get("footprint", {})
                b_coords = footprint.get("coordinates", [[]])[0]
            else:
                # Direct coords format
                b_coords = building.get("coords", [])
            
            if not b_coords or len(b_coords) < 4:
                continue
            
            # Calculate centroid
            centroid_lon = sum(c[0] for c in b_coords) / len(b_coords)
            centroid_lat = sum(c[1] for c in b_coords) / len(b_coords)
            
            dist = self._haversine_distance(lat, lon, centroid_lat, centroid_lon)
            if dist < min_distance:
                min_distance = dist
                closest_building = {
                    "coords": b_coords,
                    "tags": building.get("tags", {}),
                    "id": building.get("id", building.get("osm_id"))
                }
        
        if closest_building and min_distance < radius_m:
            # Create plot boundary as building + proper polygon buffer (not rectangle)
            plot_coords = self._buffer_building_to_plot(
                closest_building["coords"],
                lat,
                use_shapely=True  # Use proper polygon buffering
            )
            
            if plot_coords:
                area = self._calculate_polygon_area(plot_coords, lat)
                perimeter = self._calculate_perimeter(plot_coords, lat)
                
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
    
    def _create_default_plot(self, lat: float, lon: float) -> Dict[str, Any]:
        """Create default rectangular plot at location"""
        
        # Convert meters to degrees (approximate)
        lat_per_m = 1 / 111000
        lon_per_m = 1 / (111000 * math.cos(math.radians(lat)))
        
        half_width = (self.default_plot_width / 2) * lon_per_m
        half_depth = (self.default_plot_depth / 2) * lat_per_m
        
        # Create default rectangular plot (orientation may vary - this is a fallback estimate)
        coords = [
            [lon - half_width, lat - half_depth],  # SW
            [lon + half_width, lat - half_depth],  # SE
            [lon + half_width, lat + half_depth],  # NE
            [lon - half_width, lat + half_depth],  # NW
            [lon - half_width, lat - half_depth],  # Close
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
    
    def _derive_from_roads(
        self,
        lat: float,
        lon: float,
        osm_roads: List[Dict[str, Any]],
        radius_m: float
    ) -> Optional[Dict[str, Any]]:
        """
        Derive property boundary from road network
        Property line typically follows road centerlines with setback
        """
        if not SHAPELY_AVAILABLE:
            return None
        
        try:
            # Find roads near the point
            point = Point(lon, lat)
            nearby_roads = []
            
            for road in osm_roads:
                centerline = road.get("centerline", {})
                if not centerline:
                    continue
                coords = centerline.get("coordinates", [])
                if len(coords) < 2:
                    continue
                
                # Create LineString from road centerline
                road_line = LineString(coords)
                
                # Check if road is within radius
                road_point = Point(coords[len(coords)//2])  # Midpoint
                if point.distance(road_point) * 111000 < radius_m:  # Rough distance check
                    nearby_roads.append({
                        "line": road_line,
                        "width_m": road.get("width_m", 5.0),
                        "coords": coords
                    })
            
            if len(nearby_roads) < 2:
                return None  # Need at least 2 roads to form a boundary
            
            # Buffer roads to create potential property boundaries
            # Property line is typically road centerline + half road width + 1-2m setback from road edge
            buffered_roads = []
            for road in nearby_roads:
                # Buffer: half road width + 1-2m setback from road edge (not 2-3m)
                buffer_m = (road["width_m"] / 2) + 1.5  # 1.5m from road edge (within 1-2m range)
                # Use square cap and miter join to create straight edges (not rounded corners)
                buffered = road["line"].buffer(
                    buffer_m / 111000,
                    cap_style=3,  # SQUARE cap (flat ends)
                    join_style=2  # MITRE join (sharp corners, not rounded)
                )
                buffered_roads.append(buffered)
            
            # Find intersection area (where property might be)
            if len(buffered_roads) >= 2:
                # Union all buffered roads
                union = unary_union(buffered_roads)
                
                # Find polygon containing the point
                if hasattr(union, 'geoms'):
                    for geom in union.geoms:
                        if isinstance(geom, Polygon) and geom.contains(point):
                            coords = list(geom.exterior.coords)
                            # Convert to [lon, lat] format
                            plot_coords = [[c[0], c[1]] for c in coords]
                            return self._create_boundary_from_coords(
                                plot_coords, "derived_from_roads", 3.0
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
        """
        Derive property boundary from neighbor building boundaries
        Property line typically follows shared boundaries with neighbors
        """
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
                
                # Create polygon from building footprint
                building_poly = Polygon(coords)
                
                # Check if building is nearby
                if point.distance(building_poly) * 111000 < radius_m:
                    nearby_buildings.append(building_poly)
            
            if len(nearby_buildings) < 1:
                return None
            
            # Buffer buildings outward to create property boundaries
            # Property line follows neighbor's setback line (typically 1m side setback)
            # Use neighbor building's setback as property boundary
            buffered_buildings = []
            for building in nearby_buildings:
                # Buffer by neighbor's side setback (typically 1m for sides)
                # This creates the neighbor's setback line, which becomes our property line
                side_setback_m = 1.0  # Standard UK side setback
                buffered = building.buffer(
                    side_setback_m / 111000,
                    cap_style=3,  # SQUARE cap (straight edges)
                    join_style=2  # MITRE join (sharp corners, not rounded)
                )
                buffered_buildings.append(buffered)
            
            # Find area that's not covered by buffered buildings (potential property)
            if len(buffered_buildings) >= 1:
                # Create a search area
                search_radius_deg = radius_m / 111000
                search_area = point.buffer(search_radius_deg)
                
                # Subtract buffered buildings from search area
                result = search_area
                for buffered in buffered_buildings:
                    result = result.difference(buffered)
                
                # Find largest polygon containing or near the point
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
                        coords = list(largest.exterior.coords)
                        plot_coords = [[c[0], c[1]] for c in coords]
                        return self._create_boundary_from_coords(
                            plot_coords, "derived_from_neighbor_buildings", 5.0
                        )
        
        except Exception as e:
            logger.warning(f"Failed to derive boundary from neighbor buildings: {e}")
        
        return None
    
    def _buffer_building_to_plot(
        self,
        building_coords: List[List[float]],
        center_lat: float,
        use_shapely: bool = True
    ) -> List[List[float]]:
        """
        Create plot boundary by buffering building footprint
        
        Uses proper polygon buffering (not just rectangle) to preserve shape
        Uses typical UK setbacks: 5m front/rear, 2m sides
        Ensures property is always larger than building
        """
        if len(building_coords) < 4:
            return None
        
        # Use Shapely for proper polygon buffering if available
        if use_shapely and SHAPELY_AVAILABLE:
            try:
                # Close polygon if needed
                if building_coords[0] != building_coords[-1]:
                    building_coords = building_coords + [building_coords[0]]
                
                # Create Shapely polygon
                building_poly = Polygon(building_coords)
                
                # Calculate appropriate buffer distance
                # Use variable buffer: larger for front/rear, smaller for sides
                # For now, use uniform buffer, but ensure minimum sizes
                building_area = building_poly.area * (111000 ** 2) * math.cos(math.radians(center_lat))
                building_perimeter = building_poly.length * 111000
                
                # Minimum buffers
                min_buffer_m = max(2.0, building_perimeter * 0.05)  # At least 2m or 5% of perimeter
                
                # Buffer the polygon (preserves shape, not just bounding box)
                # Use square cap and miter join for straight edges (not rounded corners)
                buffered = building_poly.buffer(
                    min_buffer_m / 111000,
                    cap_style=3,  # SQUARE cap (flat ends)
                    join_style=2  # MITRE join (sharp corners, not rounded)
                )
                
                # Extract coordinates
                if isinstance(buffered, Polygon):
                    coords = list(buffered.exterior.coords)
                    return [[c[0], c[1]] for c in coords]
                
            except Exception as e:
                logger.warning(f"Shapely buffering failed, using fallback: {e}")
        
        # Fallback: use bounding box with buffers (less ideal but works)
        lat_per_m = 1 / 111000
        lon_per_m = 1 / (111000 * math.cos(math.radians(center_lat)))
        
        min_lon = min(c[0] for c in building_coords)
        max_lon = max(c[0] for c in building_coords)
        min_lat = min(c[1] for c in building_coords)
        max_lat = max(c[1] for c in building_coords)
        
        building_width_m = (max_lon - min_lon) * (111000 * math.cos(math.radians(center_lat)))
        building_depth_m = (max_lat - min_lat) * 111000
        
        front_rear_buffer = max(5.0, building_depth_m * 0.1)
        side_buffer = max(2.0, building_width_m * 0.1)
        
        buffer_lon = side_buffer * lon_per_m
        buffer_lat = front_rear_buffer * lat_per_m
        
        # Create buffered rectangle (fallback)
        return [
            [min_lon - buffer_lon, min_lat - buffer_lat],
            [max_lon + buffer_lon, min_lat - buffer_lat],
            [max_lon + buffer_lon, max_lat + buffer_lat],
            [min_lon - buffer_lon, max_lat + buffer_lat],
            [min_lon - buffer_lon, min_lat - buffer_lat],
        ]
    
    def _point_roughly_in_polygon(
        self,
        x: float,
        y: float,
        polygon: List[List[float]]
    ) -> bool:
        """Simple point-in-polygon check (ray casting)"""
        n = len(polygon)
        inside = False
        
        j = n - 1
        for i in range(n):
            xi, yi = polygon[i]
            xj, yj = polygon[j]
            
            if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
                inside = not inside
            j = i
        
        return inside
    
    def _calculate_polygon_area(
        self,
        coords: List[List[float]],
        center_lat: float
    ) -> float:
        """Calculate polygon area in square meters using shoelace formula"""
        if len(coords) < 3:
            return 0.0
        
        # Convert to local meters
        m_per_deg_lat = 111000
        m_per_deg_lon = 111000 * math.cos(math.radians(center_lat))
        
        # Convert coordinates to meters
        ref_lon = coords[0][0]
        ref_lat = coords[0][1]
        
        meters_coords = []
        for lon, lat in coords:
            x = (lon - ref_lon) * m_per_deg_lon
            y = (lat - ref_lat) * m_per_deg_lat
            meters_coords.append((x, y))
        
        # Shoelace formula
        n = len(meters_coords)
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += meters_coords[i][0] * meters_coords[j][1]
            area -= meters_coords[j][0] * meters_coords[i][1]
        
        return abs(area) / 2.0
    
    def _calculate_perimeter(
        self,
        coords: List[List[float]],
        center_lat: float
    ) -> float:
        """Calculate polygon perimeter in meters"""
        if len(coords) < 2:
            return 0.0
        
        perimeter = 0.0
        for i in range(len(coords) - 1):
            perimeter += self._haversine_distance(
                coords[i][1], coords[i][0],
                coords[i+1][1], coords[i+1][0]
            )
        
        return perimeter
    
    def _haversine_distance(
        self,
        lat1: float,
        lon1: float,
        lat2: float,
        lon2: float
    ) -> float:
        """Calculate distance between two points in meters"""
        R = 6371000  # Earth radius in meters
        
        phi1 = math.radians(lat1)
        phi2 = math.radians(lat2)
        delta_phi = math.radians(lat2 - lat1)
        delta_lambda = math.radians(lon2 - lon1)
        
        a = math.sin(delta_phi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda/2)**2
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        
        return R * c
    
    def _create_boundary_from_coords(
        self,
        coords: List[List[float]],
        source: str,
        accuracy_m: float
    ) -> Dict[str, Any]:
        """Create boundary dict from coordinates"""
        if not coords:
            return None
        
        # Close polygon if not closed
        if coords[0] != coords[-1]:
            coords = coords + [coords[0]]
        
        # Calculate area and perimeter
        center_lat = sum(c[1] for c in coords) / len(coords)
        area = self._calculate_polygon_area(coords, center_lat)
        perimeter = self._calculate_perimeter(coords, center_lat)
        
        return {
            "type": "Polygon",
            "coordinates": [coords],
            "area_sqm": area,
            "perimeter_m": perimeter,
            "source": source,
            "accuracy_m": accuracy_m
        }

