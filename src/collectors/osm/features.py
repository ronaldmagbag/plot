"""
Feature parsing for water, vegetation, and trees

Handles parsing of non-building, non-road OSM features
"""

from typing import Dict, Any, List, Optional
from .models import OSMNode, OSMWay


class FeatureProcessor:
    """Processes water, vegetation, and trees features"""
    
    @staticmethod
    def parse_water_features(ways: List[OSMWay]) -> List[Dict[str, Any]]:
        """Parse water features from ways"""
        water_features = []
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
                    
                    water_features.append({
                        "id": f"water_{way.id}",
                        "type": water_type,
                        "geometry": geometry,
                        "name": tags.get("name"),
                        "osm_tags": tags,
                        "osm_id": way.id
                    })
        return water_features
    
    @staticmethod
    def parse_vegetation_zones(ways: List[OSMWay]) -> List[Dict[str, Any]]:
        """Parse vegetation zones from ways"""
        vegetation_zones = []
        for way in ways:
            tags = way.tags
            veg_type = tags.get("landuse") or tags.get("natural") or tags.get("leisure")
            if veg_type in ["forest", "wood", "grass", "garden", "park", "meadow"]:
                coords = way.get_coordinates()
                if len(coords) >= 4:
                    if coords[0] != coords[-1]:
                        coords.append(coords[0])
                    
                    vegetation_zones.append({
                        "id": f"vegetation_{way.id}",
                        "type": veg_type,
                        "geometry": {"type": "Polygon", "coordinates": [coords]},
                        "osm_tags": tags,
                        "osm_id": way.id
                    })
        return vegetation_zones
    
    @staticmethod
    def parse_trees(nodes: Dict[int, OSMNode]) -> List[Dict[str, Any]]:
        """Parse individual trees from nodes"""
        trees = []
        for node_id, node in nodes.items():
            if node.tags.get("natural") == "tree":
                trees.append({
                    "id": f"tree_{node_id}",
                    "location": [node.lon, node.lat],
                    "species": node.tags.get("species") or node.tags.get("genus"),
                    "height_m": FeatureProcessor._parse_tree_height(node.tags.get("height")),
                    "protected": "protected" in str(node.tags).lower() or "heritage" in str(node.tags).lower(),
                    "osm_tags": node.tags,
                    "osm_id": node_id
                })
        return trees
    
    @staticmethod
    def _parse_tree_height(height_str: Optional[str]) -> Optional[float]:
        """Parse tree height from OSM tag"""
        if not height_str:
            return None
        try:
            height_str = height_str.replace("m", "").replace(" ", "")
            return float(height_str)
        except:
            return None

