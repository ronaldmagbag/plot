"""
Building-specific logic

Handles building parsing, classification, and merging
"""

import math
from typing import List, Dict, Any, Optional
from loguru import logger

try:
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - building merging will be limited")

from ...config import get_config
from .models import OSMWay


class BuildingProcessor:
    """Processes and classifies buildings from OSM data"""
    
    def __init__(self, terrain_collector: Optional[Any] = None, dsm_collector: Optional[Any] = None):
        self.config = get_config()
        self.terrain_collector = terrain_collector  # Optional terrain collector for DTM (ground elevation)
        self.dsm_collector = dsm_collector  # Optional DSM collector for surface elevation (includes buildings)
    
    def parse_buildings(self, ways: List[OSMWay]) -> List[Dict[str, Any]]:
        """
        Parse building ways into building dictionaries
        
        Args:
            ways: List of OSMWay objects
            
        Returns:
            List of building dictionaries with metadata
        """
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
                    
                    # Estimate height with terrain data if available
                    height_m = self._estimate_height(
                        tags, 
                        building_coords=coords, 
                        terrain_collector=self.terrain_collector,
                        building_osm_id=way.id
                    )
                    
                    building_ways.append({
                        "id": f"building_{way.id}",
                        "osm_id": way.id,
                        "coords": coords,
                        "tags": tags,
                        "nodes": node_set,  # Store node IDs for merging
                        "height_m": height_m,
                        "stories": self._estimate_stories(tags),
                        "building_type": self._classify_building_type(tags),
                        "usage": self._get_building_usage(tags)
                    })
                    
                    logger.info(f"Building OSM-{way.id}: Final height = {height_m:.1f}m, type = {self._classify_building_type(tags)}, stories = {self._estimate_stories(tags)}")
        
        return building_ways
    
    def merge_adjacent_buildings(
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
    
    def to_feature_dict(self, building: Dict[str, Any]) -> Dict[str, Any]:
        """Convert building dict to feature dictionary format"""
        return {
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
        }
    
    def _estimate_height(
        self, 
        tags: Dict[str, str],
        building_coords: Optional[List[List[float]]] = None,
        terrain_collector: Optional[Any] = None,
        building_osm_id: Optional[int] = None
    ) -> float:
        """
        Estimate building height from OSM tags
        
        Priority:
        1. Direct height tag from OSM
        2. Building levels calculation
        3. Building type defaults
        
        Note: DSM-DTM calculation is disabled for performance reasons
        (DSM fetching takes 30+ seconds per building)
        """
        building_id = f"OSM-{building_osm_id}" if building_osm_id else tags.get("id", "unknown")
        building_type = tags.get("building", "yes")
        
        logger.info(f"Building {building_id}: Estimating height (type: '{building_type}')")
        
        # 1. Direct height tag (highest priority)
        if "height" in tags:
            try:
                height_str = tags["height"].replace("m", "").replace(" ", "")
                height = float(height_str)
                logger.info(f"Building {building_id}: ✓ Using direct height tag '{tags['height']}' = {height:.1f}m")
                return height
            except ValueError:
                logger.warning(f"Building {building_id}: ✗ Failed to parse height tag '{tags.get('height')}'")
        
        # 2. Estimate from levels
        if "building:levels" in tags:
            try:
                levels = int(tags["building:levels"])
                height = levels * 3.0  # 3m per level
                logger.info(f"Building {building_id}: ✓ Calculated from {levels} levels = {height:.1f}m (3m per level)")
                return height
            except ValueError:
                logger.warning(f"Building {building_id}: ✗ Failed to parse building:levels '{tags.get('building:levels')}'")
        
        # 3. DSM-DTM calculation disabled for performance
        # Note: DSM fetching takes too long (30+ seconds per building)
        # Skipping DSM/DTM-based height calculation for now
        # if building_coords and len(building_coords) > 0:
        #     # DSM-DTM calculation code disabled
        pass
        
        # 4. Use building type defaults
        default_height = self.config.building_heights.get(building_type, 8.0)
        logger.info(f"Building {building_id}: → Using default height for type '{building_type}' = {default_height:.1f}m")
        logger.info(f"Building {building_id}: Height estimation complete: {default_height:.1f}m")
        return default_height
    
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

