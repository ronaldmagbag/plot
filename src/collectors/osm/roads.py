"""
Road-specific logic

Handles road parsing and classification
"""

from typing import Dict, Any, Optional, List
from .models import OSMWay
from ...config import get_config


class RoadProcessor:
    """Processes and classifies roads from OSM data"""
    
    def __init__(self):
        self.config = get_config()
    
    def parse_roads(self, ways: List[OSMWay]) -> List[Dict[str, Any]]:
        """
        Parse road ways into road dictionaries
        
        Args:
            ways: List of OSMWay objects
            
        Returns:
            List of road dictionaries with metadata
        """
        roads = []
        for way in ways:
            tags = way.tags
            if "highway" in tags:
                coords = way.get_coordinates()
                if len(coords) >= 2:
                    highway_type = tags.get("highway", "road")
                    roads.append({
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
        return roads
    
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

