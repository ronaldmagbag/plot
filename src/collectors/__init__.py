"""
Data collectors for Plot Analysis Generator

Implements Pipeline Layer 2 (Data Acquisition):
- OSMCollector: Buildings, roads, water, vegetation from OpenStreetMap
- TerrainCollector: Elevation/DEM from Mapbox Terrain-RGB or Open-Meteo
- VegetationCollector: Tree zones with RLE encoding, SAM segmentation
- SoilCollector: Soil data from SoilGrids/BGS
- BoundaryCollector: Property boundaries from cadastral sources
"""

from .osm import OSMCollector
from .elevation_collector import ElevationCollector
from .soil_collector import SoilCollector
from .boundary import BoundaryCollector
from .terrain_collector import TerrainCollector
from .vegetation_collector import VegetationCollector
from .mapbox_imagery_collector import MapboxImageryCollector

__all__ = [
    "OSMCollector",
    "ElevationCollector", 
    "SoilCollector",
    "BoundaryCollector",
    "TerrainCollector",
    "VegetationCollector",
    "MapboxImageryCollector",
]

