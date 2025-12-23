"""
OpenStreetMap data collector using Overpass API
Collects buildings, roads, water features, and vegetation data

Uses BATCH QUERY approach: single request for all features to reduce API calls
"""

import time
import os
import json
import hashlib
import math
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
import requests
from loguru import logger

try:
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - building merging will be limited")

from ..config import get_config


@dataclass
class OSMNode:
    id: int
    lat: float
    lon: float
    tags: Dict[str, str]


@dataclass
class OSMWay:
    id: int
    nodes: List[OSMNode]
    tags: Dict[str, str]
    geometry: Optional[List[List[float]]] = None  # Direct geometry from Overpass
    
    def get_coordinates(self) -> List[List[float]]:
        """Get coordinates as [lon, lat] list"""
        # Prefer direct geometry if available (from 'out geom')
        if self.geometry:
            return self.geometry
        # Fallback to node-based coordinates
        return [[n.lon, n.lat] for n in self.nodes]


class OSMCollector:
    """
    Collect data from OpenStreetMap via Overpass API
    
    Uses batch query approach: single request for all features
    to minimize API calls and avoid rate limiting.
    
    Supports caching to disk for debugging and reuse.
    """
    
    def __init__(self, cache_dir: Optional[str] = None):
        self.config = get_config()
        self.overpass_url = self.config.api.overpass_url
        self.timeout = self.config.api.overpass_timeout
        self._last_request_time = 0
        self._min_request_interval = 2.0  # Increased delay between requests
        self.cache_dir = cache_dir
        
    def _rate_limit(self):
        """Ensure we don't exceed rate limits"""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()
    
    def _query_overpass(self, query: str, retry_delay: float = 5.0) -> Dict[str, Any]:
        """
        Execute Overpass API query with retry logic
        
        Args:
            query: Overpass QL query string
            retry_delay: Initial delay between retries (increases with attempts)
        """
        self._rate_limit()
        
        headers = {
            "User-Agent": self.config.api.user_agent,
            "Content-Type": "application/x-www-form-urlencoded"
        }
        
        for attempt in range(self.config.api.max_retries):
            try:
                response = requests.post(
                    self.overpass_url,
                    data={"data": query},
                    headers=headers,
                    timeout=self.timeout
                )
                response.raise_for_status()
                return response.json()
            except requests.exceptions.Timeout:
                wait_time = retry_delay * (attempt + 1)
                logger.warning(f"Overpass timeout (attempt {attempt + 1}/{self.config.api.max_retries}). Retrying in {wait_time}s...")
                if attempt < self.config.api.max_retries - 1:
                    time.sleep(wait_time)
                else:
                    logger.error(f"OSM API failed: Overpass timeout after {self.config.api.max_retries} attempts")
                    raise RuntimeError(f"Overpass API timeout after {self.config.api.max_retries} attempts") from e
            except requests.exceptions.HTTPError as e:
                if e.response.status_code in [429, 504] and attempt < self.config.api.max_retries - 1:
                    wait_time = retry_delay * (attempt + 1)
                    logger.warning(f"Overpass {e.response.status_code} (attempt {attempt + 1}/{self.config.api.max_retries}). Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    logger.error(f"OSM API failed: HTTP {e.response.status_code} after {self.config.api.max_retries} attempts")
                    raise RuntimeError(f"Overpass API HTTP error {e.response.status_code} after {self.config.api.max_retries} attempts") from e
            except requests.exceptions.RequestException as e:
                logger.warning(f"Overpass request failed (attempt {attempt + 1}): {e}")
                if attempt < self.config.api.max_retries - 1:
                    time.sleep(retry_delay * (attempt + 1))
                else:
                    logger.error(f"OSM API failed: Request exception after {self.config.api.max_retries} attempts: {e}")
                    raise RuntimeError(f"Overpass API request failed after {self.config.api.max_retries} attempts: {e}") from e
        
        return {"elements": []}
    
    def _parse_elements(self, data: Dict[str, Any]) -> Tuple[Dict[int, OSMNode], List[OSMWay]]:
        """
        Parse Overpass response into nodes and ways
        
        Handles both 'out body' (node references) and 'out geom' (direct geometry) formats
        """
        nodes = {}
        ways = []
        
        for element in data.get("elements", []):
            if element["type"] == "node":
                nodes[element["id"]] = OSMNode(
                    id=element["id"],
                    lat=element["lat"],
                    lon=element["lon"],
                    tags=element.get("tags", {})
                )
            elif element["type"] == "way":
                # Check if geometry is directly provided (from 'out geom')
                geometry = None
                if "geometry" in element:
                    # Overpass 'out geom' provides geometry as list of {lat, lon} objects
                    # Convert to [lon, lat] format for consistency
                    geometry = []
                    for node in element["geometry"]:
                        if isinstance(node, dict):
                            # Format: {"lat": ..., "lon": ...}
                            geometry.append([node.get("lon"), node.get("lat")])
                        elif isinstance(node, list) and len(node) >= 2:
                            # Format: [lon, lat] (already correct)
                            geometry.append(node)
                
                # Build node list (for 'out body' format)
                way_nodes = []
                if "nodes" in element:
                    for node_id in element.get("nodes", []):
                        if node_id in nodes:
                            way_nodes.append(nodes[node_id])
                
                # Create way (works with or without resolved nodes)
                ways.append(OSMWay(
                    id=element["id"],
                    nodes=way_nodes,
                    tags=element.get("tags", {}),
                    geometry=geometry
                ))
        
        return nodes, ways
    
    def _get_cache_path(self, lat: float, lon: float, radius_m: float) -> Optional[str]:
        """Get cache file path for OSM query"""
        if not self.cache_dir:
            return None
        # Create hash of query parameters for cache key
        cache_key = f"{lat:.6f}_{lon:.6f}_{radius_m:.0f}"
        cache_hash = hashlib.md5(cache_key.encode()).hexdigest()[:8]
        return os.path.join(self.cache_dir, f"osm_{cache_hash}.json")
    
    def _load_from_cache(self, cache_path: str) -> Optional[Dict[str, Any]]:
        """Load OSM data from cache if exists"""
        if os.path.exists(cache_path):
            try:
                with open(cache_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    logger.info(f"Loaded OSM data from cache: {cache_path}")
                    return data
            except Exception as e:
                logger.warning(f"Failed to load cache {cache_path}: {e}")
        return None
    
    def _save_to_cache(self, cache_path: str, data: Dict[str, Any]):
        """Save OSM data to cache"""
        if not self.cache_dir:
            return
        try:
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            with open(cache_path, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            logger.info(f"Saved OSM data to cache: {cache_path}")
        except Exception as e:
            logger.warning(f"Failed to save cache {cache_path}: {e}")
    
    def _parse_cached_osm_data(self, cached_data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse cached Overpass response into the same format as fresh data"""
        nodes, ways = self._parse_elements(cached_data)
        
        # Parse features by type (same logic as fetch_all_features)
        result = {
            "buildings": [],
            "roads": [],
            "water": [],
            "vegetation_zones": [],
            "trees": [],
            "landuse": []
        }
        
        # Parse buildings
        for way in ways:
            tags = way.tags
            if "building" in tags:
                coords = way.get_coordinates()
                if len(coords) >= 4:
                    if coords[0] != coords[-1]:
                        coords.append(coords[0])
                    
                    result["buildings"].append({
                        "id": f"building_{way.id}",
                        "footprint": {
                            "type": "Polygon",
                            "coordinates": [coords]
                        },
                        "height_m": self._estimate_height(tags),
                        "stories": self._estimate_stories(tags),
                        "building_type": self._classify_building_type(tags),
                        "usage": self._get_building_usage(tags),
                        "osm_tags": tags,
                        "osm_id": way.id
                    })
        
        # Parse roads
        for way in ways:
            tags = way.tags
            if "highway" in tags:
                coords = way.get_coordinates()
                if len(coords) >= 2:
                    highway_type = tags.get("highway", "road")
                    result["roads"].append({
                        "id": f"road_{way.id}",
                        "name": tags.get("name", "Unnamed Road"),
                        "type": self._classify_road_type(highway_type),
                        "centerline": {
                            "type": "LineString",
                            "coordinates": coords
                        },
                        "width_m": self._estimate_road_width(highway_type),
                        "traffic_level": self._estimate_traffic_level(highway_type),
                        "speed_limit_mph": self._parse_speed_limit(tags.get("maxspeed")),
                        "has_sidewalk": tags.get("sidewalk", "no") != "no",
                        "osm_tags": tags,
                        "osm_id": way.id
                    })
        
        # Parse water features
        for way in ways:
            tags = way.tags
            water_type = tags.get("waterway") or tags.get("natural") or tags.get("landuse")
            if water_type in ["water", "river", "stream", "canal", "ditch", "drain", "basin", "reservoir"]:
                coords = way.get_coordinates()
                if len(coords) >= 2:
                    if coords[0] == coords[-1] and len(coords) >= 4:
                        geometry = {"type": "Polygon", "coordinates": [coords]}
                    else:
                        geometry = {"type": "LineString", "coordinates": coords}
                    
                    result["water"].append({
                        "id": f"water_{way.id}",
                        "type": water_type,
                        "geometry": geometry,
                        "name": tags.get("name"),
                        "osm_tags": tags,
                        "osm_id": way.id
                    })
        
        # Parse vegetation zones
        for way in ways:
            tags = way.tags
            veg_type = tags.get("landuse") or tags.get("natural") or tags.get("leisure")
            if veg_type in ["forest", "wood", "grass", "garden", "park", "meadow"]:
                coords = way.get_coordinates()
                if len(coords) >= 4:
                    if coords[0] != coords[-1]:
                        coords.append(coords[0])
                    
                    result["vegetation_zones"].append({
                        "id": f"vegetation_{way.id}",
                        "type": veg_type,
                        "geometry": {"type": "Polygon", "coordinates": [coords]},
                        "osm_tags": tags,
                        "osm_id": way.id
                    })
        
        # Parse landuse
        for way in ways:
            tags = way.tags
            if "landuse" in tags:
                coords = way.get_coordinates()
                if len(coords) >= 4:
                    if coords[0] != coords[-1]:
                        coords.append(coords[0])
                    
                    result["landuse"].append({
                        "id": f"landuse_{way.id}",
                        "type": tags.get("landuse"),
                        "geometry": {"type": "Polygon", "coordinates": [coords]},
                        "osm_tags": tags,
                        "osm_id": way.id
                    })
        
        # Parse individual trees
        for node_id, node in nodes.items():
            if node.tags.get("natural") == "tree":
                result["trees"].append({
                    "id": f"tree_{node_id}",
                    "location": [node.lon, node.lat],
                    "species": node.tags.get("species") or node.tags.get("genus"),
                    "height_m": self._parse_tree_height(node.tags.get("height")),
                    "protected": "protected" in str(node.tags).lower() or "heritage" in str(node.tags).lower(),
                    "osm_tags": node.tags,
                    "osm_id": node_id
                })
        
        logger.info(f"Cached data parsed: {len(result['buildings'])} buildings, "
                   f"{len(result['roads'])} roads, {len(result['water'])} water features, "
                   f"{len(result['trees'])} trees, {len(result['vegetation_zones'])} vegetation zones")
        
        return result
    
    def fetch_all_features(
        self,
        lat: float,
        lon: float,
        radius_m: float = 150
    ) -> Dict[str, Any]:
        """
        Fetch ALL OSM features in a SINGLE batch query
        
        This reduces API calls from 4+ separate requests to just 1,
        avoiding rate limiting and improving reliability.
        
        Supports caching to disk for debugging and reuse.
        
        Args:
            lat: Center latitude
            lon: Center longitude
            radius_m: Search radius in meters
            
        Returns:
            Dict with keys: 'buildings', 'roads', 'water', 'vegetation', 'landuse'
        """
        # Check cache first
        cache_path = self._get_cache_path(lat, lon, radius_m)
        if cache_path:
            cached_data = self._load_from_cache(cache_path)
            if cached_data:
                # Parse cached data into the same format
                return self._parse_cached_osm_data(cached_data)
        
        logger.info(f"Fetching ALL OSM features within {radius_m}m of ({lat}, {lon}) in single batch query")
        
        # Single comprehensive Overpass query
        query = f"""
        [out:json][timeout:{self.timeout}];
        (
            // Buildings
            way["building"](around:{radius_m},{lat},{lon});
            relation["building"](around:{radius_m},{lat},{lon});
            
            // Roads (all highway types)
            way["highway"](around:{radius_m},{lat},{lon});
            
            // Water features
            way["natural"="water"](around:{radius_m},{lat},{lon});
            way["waterway"](around:{radius_m},{lat},{lon});
            relation["natural"="water"](around:{radius_m},{lat},{lon});
            way["landuse"="basin"](around:{radius_m},{lat},{lon});
            way["landuse"="reservoir"](around:{radius_m},{lat},{lon});
            
            // Vegetation zones
            way["landuse"="forest"](around:{radius_m},{lat},{lon});
            way["natural"="wood"](around:{radius_m},{lat},{lon});
            way["landuse"="grass"](around:{radius_m},{lat},{lon});
            way["leisure"="garden"](around:{radius_m},{lat},{lon});
            
            // Individual trees
            node["natural"="tree"](around:{radius_m},{lat},{lon});
            
            // Landuse (for boundary detection)
            way["landuse"="residential"](around:{radius_m},{lat},{lon});
            relation["landuse"="residential"](around:{radius_m},{lat},{lon});
        );
        out geom;
        """
        
        try:
            data = self._query_overpass(query)
            
            # Save raw Overpass response to cache
            if cache_path:
                self._save_to_cache(cache_path, data)
            
            nodes, ways = self._parse_elements(data)
            
            # Parse features by type
            result = {
                "buildings": [],
                "roads": [],
                "water": [],
                "vegetation_zones": [],
                "trees": [],
                "landuse": []
            }
            
            # Parse buildings (filter out point buildings and merge adjacent/overlapping ones)
            building_ways = []
            for way in ways:
                tags = way.tags
                if "building" in tags:
                    coords = way.get_coordinates()
                    # Skip point buildings (less than 3 unique points)
                    if len(coords) < 3:
                        continue
                    # Skip if it's actually a point (all coordinates are the same)
                    if len(set((round(c[0], 8), round(c[1], 8)) for c in coords)) < 3:
                        continue
                    
                    # Ensure polygon is closed
                    if len(coords) >= 3:
                        if coords[0] != coords[-1]:
                            coords.append(coords[0])
                        
                        # Track nodes for this building (for merging buildings that share nodes)
                        node_set = set()
                        if way.nodes:
                            node_set = {n.id for n in way.nodes}
                        
                        building_ways.append({
                            "id": f"building_{way.id}",
                            "osm_id": way.id,
                            "coords": coords,
                            "tags": tags,
                            "nodes": node_set,  # Store node IDs for merging
                            "height_m": self._estimate_height(tags),
                            "stories": self._estimate_stories(tags),
                            "building_type": self._classify_building_type(tags),
                            "usage": self._get_building_usage(tags)
                        })
            
            # Merge adjacent/overlapping buildings (within 2m threshold)
            merged_buildings = self._merge_adjacent_buildings(building_ways, merge_threshold_m=2.0)
            
            # Add merged buildings to result
            for building in merged_buildings:
                result["buildings"].append({
                    "id": building["id"],
                    "footprint": {
                        "type": "Polygon",
                        "coordinates": [building["coords"]]
                    },
                    "height_m": building["height_m"],
                    "stories": building["stories"],
                    "building_type": building["building_type"],
                    "usage": building["usage"],
                    "osm_tags": building["tags"],
                    "osm_id": building["osm_id"]
                })
            
            # Parse roads
            for way in ways:
                tags = way.tags
                if "highway" in tags:
                    coords = way.get_coordinates()
                    if len(coords) >= 2:
                        highway_type = tags.get("highway", "road")
                        result["roads"].append({
                            "id": f"road_{way.id}",
                            "name": tags.get("name", "Unnamed Road"),
                            "type": self._classify_road_type(highway_type),
                            "centerline": {
                                "type": "LineString",
                                "coordinates": coords
                            },
                            "width_m": self._estimate_road_width(highway_type),
                            "traffic_level": self._estimate_traffic_level(highway_type),
                            "speed_limit_mph": self._parse_speed_limit(tags.get("maxspeed")),
                            "has_sidewalk": tags.get("sidewalk", "no") != "no",
                            "osm_tags": tags,
                            "osm_id": way.id
                        })
            
            # Parse water features
            for way in ways:
                tags = way.tags
                water_type = tags.get("waterway") or tags.get("natural") or tags.get("landuse")
                if water_type in ["water", "river", "stream", "canal", "ditch", "drain", "basin", "reservoir"]:
                    coords = way.get_coordinates()
                    if len(coords) >= 2:
                        # Determine geometry type
                        if coords[0] == coords[-1] and len(coords) >= 4:
                            geometry = {"type": "Polygon", "coordinates": [coords]}
                        else:
                            geometry = {"type": "LineString", "coordinates": coords}
                        
                        result["water"].append({
                            "id": f"water_{way.id}",
                            "type": water_type,
                            "geometry": geometry,
                            "name": tags.get("name"),
                            "osm_tags": tags,
                            "osm_id": way.id
                        })
            
            # Parse vegetation zones
            for way in ways:
                tags = way.tags
                veg_type = tags.get("landuse") or tags.get("natural") or tags.get("leisure")
                if veg_type in ["forest", "wood", "grass", "garden", "park", "meadow"]:
                    coords = way.get_coordinates()
                    if len(coords) >= 4:
                        if coords[0] != coords[-1]:
                            coords.append(coords[0])
                        
                        result["vegetation_zones"].append({
                            "id": f"vegetation_{way.id}",
                            "type": veg_type,
                            "geometry": {"type": "Polygon", "coordinates": [coords]},
                            "osm_tags": tags,
                            "osm_id": way.id
                        })
            
            # Parse landuse
            for way in ways:
                tags = way.tags
                if "landuse" in tags:
                    coords = way.get_coordinates()
                    if len(coords) >= 4:
                        if coords[0] != coords[-1]:
                            coords.append(coords[0])
                        
                        result["landuse"].append({
                            "id": f"landuse_{way.id}",
                            "type": tags.get("landuse"),
                            "geometry": {"type": "Polygon", "coordinates": [coords]},
                            "osm_tags": tags,
                            "osm_id": way.id
                        })
            
            # Parse individual trees (from nodes)
            for node_id, node in nodes.items():
                if node.tags.get("natural") == "tree":
                    result["trees"].append({
                        "id": f"tree_{node_id}",
                        "location": [node.lon, node.lat],
                        "species": node.tags.get("species") or node.tags.get("genus"),
                        "height_m": self._parse_tree_height(node.tags.get("height")),
                        "protected": "protected" in str(node.tags).lower() or "heritage" in str(node.tags).lower(),
                        "osm_tags": node.tags,
                        "osm_id": node_id
                    })
            
            logger.info(f"Batch query results: {len(result['buildings'])} buildings, "
                       f"{len(result['roads'])} roads, {len(result['water'])} water features, "
                       f"{len(result['trees'])} trees, {len(result['vegetation_zones'])} vegetation zones")
            
            return result
            
        except Exception as e:
            logger.error(f"OSM API query failed after all retries: {e}")
            logger.error("Cannot proceed without OSM data. Please retry the request.")
            # Raise exception instead of returning empty results
            raise RuntimeError(f"Failed to fetch OSM data from Overpass API: {e}") from e
    
    def _merge_adjacent_buildings(
        self,
        buildings: List[Dict[str, Any]],
        merge_threshold_m: float = 0.5
    ) -> List[Dict[str, Any]]:
        """
        Merge buildings that are adjacent or overlapping (within threshold distance)
        
        Args:
            buildings: List of building dicts with coords, tags, etc.
            merge_threshold_m: Maximum distance between buildings to merge (meters)
        
        Returns:
            List of merged buildings
        """
        if not buildings:
            return []
        
        if not SHAPELY_AVAILABLE:
            # Without Shapely, can't do proper merging - return as-is
            logger.warning("Shapely not available - cannot merge adjacent buildings")
            return buildings
        
        try:
            # Convert to Shapely polygons
            building_polys = []
            for b in buildings:
                try:
                    # Remove duplicate closing point if present
                    coords = b["coords"]
                    if coords[0] == coords[-1] and len(coords) > 1:
                        coords = coords[:-1]
                    
                    # Create polygon (lon, lat order for Shapely)
                    poly_coords = [(c[0], c[1]) for c in coords]
                    if len(poly_coords) >= 3:
                        poly = Polygon(poly_coords)
                        if poly.is_valid:
                            building_polys.append({
                                "poly": poly,
                                "building": b
                            })
                except Exception as e:
                    logger.warning(f"Failed to create polygon for building {b.get('id', 'unknown')}: {e}")
                    continue
            
            if not building_polys:
                return buildings
            
            # Merge overlapping/adjacent buildings
            merged = []
            used = set()
            
            for i, bp1 in enumerate(building_polys):
                if i in used:
                    continue
                
                poly1 = bp1["poly"]
                merged_group = [bp1["building"]]
                used.add(i)
                
                # Find all buildings that overlap or are adjacent to this one
                b1_nodes = bp1["building"].get("nodes", set())
                
                for j, bp2 in enumerate(building_polys):
                    if j <= i or j in used:
                        continue
                    
                    poly2 = bp2["poly"]
                    b2_nodes = bp2["building"].get("nodes", set())
                    
                    # First check: if buildings share nodes, they should definitely be merged
                    shares_nodes = bool(b1_nodes and b2_nodes and b1_nodes.intersection(b2_nodes))
                    
                    # Second check: if buildings overlap, touch, or are very close
                    distance = poly1.distance(poly2)
                    
                    # Convert distance (degrees) to meters (rough approximation)
                    # At UK latitude (~51°), 1 degree ≈ 111km
                    # Use average latitude for better accuracy
                    avg_lat = (poly1.centroid.y + poly2.centroid.y) / 2 if hasattr(poly1, 'centroid') else 51.0
                    m_per_deg_lat = 111000
                    m_per_deg_lon = 111000 * math.cos(math.radians(avg_lat))
                    
                    # For distance calculation, use the larger of lat/lon conversion
                    # Distance is in degrees, convert to meters using average
                    distance_m = distance * m_per_deg_lat  # Use lat conversion as approximation
                    
                    # Merge if: share nodes, overlap, touch, or are very close
                    if shares_nodes or distance_m <= merge_threshold_m or poly1.intersects(poly2) or poly1.touches(poly2):
                        # Merge these buildings
                        merged_group.append(bp2["building"])
                        used.add(j)
                        # Update node set for merged group
                        if b2_nodes:
                            b1_nodes = b1_nodes.union(b2_nodes)
                
                # If multiple buildings, merge them into one
                if len(merged_group) > 1:
                    # Union all polygons
                    polys_to_merge = []
                    for b in merged_group:
                        coords = b["coords"]
                        if coords[0] == coords[-1] and len(coords) > 1:
                            coords = coords[:-1]
                        poly_coords = [(c[0], c[1]) for c in coords]
                        if len(poly_coords) >= 3:
                            try:
                                polys_to_merge.append(Polygon(poly_coords))
                            except:
                                pass
                    
                    if polys_to_merge:
                        try:
                            merged_poly = unary_union(polys_to_merge)
                            if isinstance(merged_poly, Polygon) and merged_poly.is_valid:
                                # Extract coordinates from merged polygon
                                merged_coords = [[c[0], c[1]] for c in merged_poly.exterior.coords]
                                if merged_coords[0] != merged_coords[-1]:
                                    merged_coords.append(merged_coords[0])
                                
                                # Use first building's metadata, combine IDs
                                merged_ids = [b["osm_id"] for b in merged_group]
                                merged_building = {
                                    "id": f"building_merged_{'_'.join(str(oid) for oid in merged_ids)}",
                                    "osm_id": merged_ids[0],  # Primary ID
                                    "coords": merged_coords,
                                    "tags": merged_group[0]["tags"],  # Use first building's tags
                                    "height_m": max(b["height_m"] for b in merged_group),  # Use max height
                                    "stories": max(b["stories"] for b in merged_group),  # Use max stories
                                    "building_type": merged_group[0]["building_type"],
                                    "usage": merged_group[0]["usage"]
                                }
                                merged.append(merged_building)
                                continue
                        except Exception as e:
                            logger.warning(f"Failed to merge buildings: {e}")
                
                # Single building (not merged)
                merged.append(merged_group[0])
            
            return merged
            
        except Exception as e:
            logger.warning(f"Building merging failed: {e}, returning original buildings")
            return buildings
    
    # ============================================================
    # Helper Methods
    # ============================================================
    
    def _estimate_height(self, tags: Dict[str, str]) -> float:
        """Estimate building height from OSM tags"""
        # Direct height tag
        if "height" in tags:
            try:
                height_str = tags["height"].replace("m", "").replace(" ", "")
                return float(height_str)
            except ValueError:
                pass
        
        # Estimate from levels
        if "building:levels" in tags:
            try:
                levels = int(tags["building:levels"])
                return levels * 3.0  # 3m per level
            except ValueError:
                pass
        
        # Use building type defaults
        building_type = tags.get("building", "yes")
        return self.config.building_heights.get(building_type, 8.0)
    
    def _estimate_stories(self, tags: Dict[str, str]) -> int:
        """Estimate number of stories from OSM tags"""
        if "building:levels" in tags:
            try:
                return int(tags["building:levels"])
            except ValueError:
                pass
        
        # Estimate from height
        height = self._estimate_height(tags)
        return max(1, int(height / 3.0))
    
    def _classify_building_type(self, tags: Dict[str, str]) -> str:
        """Classify building type from OSM tags"""
        building = tags.get("building", "unknown")
        
        type_mapping = {
            "house": "detached_house",
            "detached": "detached_house",
            "semidetached_house": "semi_detached_house",
            "semi": "semi_detached_house",
            "terrace": "terraced_house",
            "apartments": "apartments",
            "flat": "apartments",
            "residential": "residential",
            "commercial": "commercial",
            "retail": "retail",
            "industrial": "industrial",
            "garage": "garage",
            "shed": "shed",
            "yes": "building",
        }
        
        return type_mapping.get(building, building)
    
    def _get_building_usage(self, tags: Dict[str, str]) -> str:
        """Get building usage from OSM tags"""
        if "residential" in tags.get("building", ""):
            return "residential"
        if tags.get("building") in ["house", "detached", "semi", "semidetached_house", "terrace", "apartments"]:
            return "residential"
        if tags.get("building") in ["commercial", "retail", "office"]:
            return "commercial"
        if tags.get("building") in ["industrial", "warehouse"]:
            return "industrial"
        return "unknown"
    
    def _classify_road_type(self, highway_type: str) -> str:
        """Classify road type for display"""
        type_mapping = {
            "motorway": "motorway",
            "trunk": "trunk",
            "primary": "primary",
            "secondary": "secondary",
            "tertiary": "tertiary",
            "unclassified": "unclassified",
            "residential": "local_residential",
            "living_street": "local_residential",
            "service": "service",
            "footway": "footpath",
            "path": "footpath",
            "cycleway": "cycleway",
        }
        return type_mapping.get(highway_type, "road")
    
    def _estimate_road_width(self, highway_type: str) -> float:
        """Estimate road width from type"""
        return self.config.road_widths.get(highway_type, 5.0)
    
    def _estimate_traffic_level(self, highway_type: str) -> str:
        """Estimate traffic level from road type"""
        if highway_type in ["motorway", "trunk", "primary"]:
            return "high"
        elif highway_type in ["secondary", "tertiary"]:
            return "medium"
        else:
            return "low"
    
    def _parse_speed_limit(self, maxspeed_str: Optional[str]) -> Optional[int]:
        """Parse speed limit from OSM tag"""
        if not maxspeed_str:
            return None
        try:
            # Extract number (handles "30 mph", "30", etc.)
            import re
            match = re.search(r'\d+', maxspeed_str)
            if match:
                return int(match.group())
        except:
            pass
        return None
    
    def _parse_tree_height(self, height_str: Optional[str]) -> Optional[float]:
        """Parse tree height from OSM tag"""
        if not height_str:
            return None
        try:
            height_str = height_str.replace("m", "").replace(" ", "")
            return float(height_str)
        except:
            return None
