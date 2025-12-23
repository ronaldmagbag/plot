"""
OSM response parser

Parses Overpass API responses into OSMNode and OSMWay objects
"""

from typing import Dict, Any, Tuple, List
from .models import OSMNode, OSMWay


class OSMResponseParser:
    """Parses Overpass API responses"""
    
    @staticmethod
    def parse_elements(data: Dict[str, Any]) -> Tuple[Dict[int, OSMNode], List[OSMWay]]:
        """
        Parse Overpass response into nodes and ways
        
        Handles both 'out body' (node references) and 'out geom' (direct geometry) formats
        
        Args:
            data: JSON response from Overpass API
            
        Returns:
            Tuple of (nodes dict, ways list)
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

