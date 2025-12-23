"""
Main OSM Collector

Orchestrates all OSM data collection components
"""

from typing import Dict, Any, Optional
from loguru import logger

from .api_client import OverpassAPIClient
from .cache import OSMCache
from .parser import OSMResponseParser
from .buildings import BuildingProcessor
from .roads import RoadProcessor
from .features import FeatureProcessor
from ...config import get_config


class OSMCollector:
    """
    Collect data from OpenStreetMap via Overpass API
    
    Uses batch query approach: single request for all features
    to minimize API calls and avoid rate limiting.
    
    Supports caching to disk for debugging and reuse.
    """
    
    def __init__(self, cache_dir: Optional[str] = None):
        self.config = get_config()
        self.api_client = OverpassAPIClient()
        self.cache = OSMCache(cache_dir)
        self.parser = OSMResponseParser()
        self.building_processor = BuildingProcessor()
        self.road_processor = RoadProcessor()
        self.feature_processor = FeatureProcessor()
        self.timeout = self.config.api.overpass_timeout
    
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
            Dict with keys: 'buildings', 'roads', 'water', 'vegetation_zones', 'trees', 'landuse'
        """
        # Check cache first
        cache_path = self.cache.get_cache_path(lat, lon, radius_m)
        if cache_path:
            cached_data = self.cache.load(cache_path)
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
            data = self.api_client.query(query)
            
            # Save raw Overpass response to cache
            if cache_path:
                self.cache.save(cache_path, data)
            
            nodes, ways = self.parser.parse_elements(data)
            
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
            building_ways = self.building_processor.parse_buildings(ways)
            
            # Merge adjacent/overlapping buildings (within 2m threshold)
            merged_buildings = self.building_processor.merge_adjacent_buildings(
                building_ways, 
                merge_threshold_m=2.0
            )
            
            # Add merged buildings to result
            for building in merged_buildings:
                result["buildings"].append(
                    self.building_processor.to_feature_dict(building)
                )
            
            # Parse roads
            result["roads"] = self.road_processor.parse_roads(ways)
            
            # Parse water features
            result["water"] = self.feature_processor.parse_water_features(ways)
            
            # Parse vegetation zones
            result["vegetation_zones"] = self.feature_processor.parse_vegetation_zones(ways)
            
            # Parse landuse
            result["landuse"] = self.feature_processor.parse_landuse(ways)
            
            # Parse individual trees (from nodes)
            result["trees"] = self.feature_processor.parse_trees(nodes)
            
            logger.info(f"Batch query results: {len(result['buildings'])} buildings, "
                       f"{len(result['roads'])} roads, {len(result['water'])} water features, "
                       f"{len(result['trees'])} trees, {len(result['vegetation_zones'])} vegetation zones")
            
            return result
            
        except Exception as e:
            logger.error(f"OSM API query failed after all retries: {e}")
            logger.error("Cannot proceed without OSM data. Please retry the request.")
            # Raise exception instead of returning empty results
            raise RuntimeError(f"Failed to fetch OSM data from Overpass API: {e}") from e
    
    def _parse_cached_osm_data(self, cached_data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse cached Overpass response into the same format as fresh data"""
        nodes, ways = self.parser.parse_elements(cached_data)
        
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
        building_ways = self.building_processor.parse_buildings(ways)
        merged_buildings = self.building_processor.merge_adjacent_buildings(
            building_ways, 
            merge_threshold_m=2.0
        )
        for building in merged_buildings:
            result["buildings"].append(
                self.building_processor.to_feature_dict(building)
            )
        
        # Parse roads
        result["roads"] = self.road_processor.parse_roads(ways)
        
        # Parse water features
        result["water"] = self.feature_processor.parse_water_features(ways)
        
        # Parse vegetation zones
        result["vegetation_zones"] = self.feature_processor.parse_vegetation_zones(ways)
        
        # Parse landuse
        result["landuse"] = self.feature_processor.parse_landuse(ways)
        
        # Parse individual trees
        result["trees"] = self.feature_processor.parse_trees(nodes)
        
        logger.info(f"Cached data parsed: {len(result['buildings'])} buildings, "
                   f"{len(result['roads'])} roads, {len(result['water'])} water features, "
                   f"{len(result['trees'])} trees, {len(result['vegetation_zones'])} vegetation zones")
        
        return result

