"""
OSM data caching

Handles caching of Overpass API responses to disk
"""

import os
import json
import hashlib
from typing import Dict, Any, Optional
from loguru import logger


class OSMCache:
    """Handles caching of OSM data to disk"""
    
    def __init__(self, cache_dir: Optional[str] = None):
        self.cache_dir = cache_dir
    
    def get_cache_path(self, lat: float, lon: float, radius_m: float) -> Optional[str]:
        """Get cache file path for OSM query"""
        if not self.cache_dir:
            return None
        # Create hash of query parameters for cache key
        cache_key = f"{lat:.6f}_{lon:.6f}_{radius_m:.0f}"
        cache_hash = hashlib.md5(cache_key.encode()).hexdigest()[:8]
        return os.path.join(self.cache_dir, f"osm_{cache_hash}.json")
    
    def load(self, cache_path: str) -> Optional[Dict[str, Any]]:
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
    
    def save(self, cache_path: str, data: Dict[str, Any]):
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

