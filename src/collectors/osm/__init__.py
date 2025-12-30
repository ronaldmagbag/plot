"""
OpenStreetMap data collection module

Modular OSM data collector with separate components for:
- API client: Overpass API communication
- Models: Data structures (OSMNode, OSMWay)
- Parser: Response parsing
- Buildings: Building-specific logic
- Roads: Road-specific logic
- Features: Water, vegetation, trees parsing
- Cache: Caching functionality
- Collector: Main orchestrator class
"""

from .models import OSMNode, OSMWay
from .collector import OSMCollector

__all__ = [
    "OSMNode",
    "OSMWay",
    "OSMCollector",
]

