"""
OSM data models

Data classes for representing OSM nodes and ways
"""

from typing import List, Dict, Optional
from dataclasses import dataclass


@dataclass
class OSMNode:
    """Represents an OSM node (point)"""
    id: int
    lat: float
    lon: float
    tags: Dict[str, str]


@dataclass
class OSMWay:
    """Represents an OSM way (line or polygon)"""
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

