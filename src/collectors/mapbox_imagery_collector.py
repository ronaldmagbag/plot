"""
Mapbox Satellite Imagery Collector

Downloads high-resolution Mapbox satellite tiles for a given area and merges them
into a single image for SAM3 segmentation.
"""

import os
import math
from typing import Dict, Any, Optional, Tuple, List
from pathlib import Path
import requests
from io import BytesIO
from loguru import logger
import numpy as np

from ..config import get_config

try:
    from dotenv import load_dotenv
    HAS_DOTENV = True
except ImportError:
    HAS_DOTENV = False

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    logger.error("PIL required for Mapbox imagery collection")

try:
    import mercantile
    HAS_MERCANTILE = True
except ImportError:
    HAS_MERCANTILE = False
    logger.error("mercantile required for Mapbox imagery collection")


class MapboxImageryCollector:
    """
    Downloads and merges Mapbox satellite imagery tiles for a given area.
    
    Uses the highest available zoom level (typically 20) and merges multiple
    tiles to cover the requested radius.
    """
    
    # Mapbox Satellite tile URL template
    MAPBOX_SATELLITE_URL = "https://api.mapbox.com/styles/v1/mapbox/satellite-v9/tiles/{z}/{x}/{y}?access_token={token}"
    
    def __init__(self, cache_dir: Optional[str] = None):
        """
        Initialize Mapbox imagery collector.
        
        Args:
            cache_dir: Directory to cache downloaded tiles and merged images
        """
        self.config = get_config()
        
        # Load .env file if available
        if HAS_DOTENV:
            # Try multiple locations for .env file
            env_paths = [
                Path(__file__).parent.parent.parent / ".env",  # Project root
                Path.cwd() / ".env",  # Current working directory
                Path.home() / ".env",  # Home directory
            ]
            
            loaded = False
            for env_path in env_paths:
                if env_path.exists():
                    load_dotenv(env_path, override=False)  # Don't override existing env vars
                    logger.info(f"Loaded .env file from {env_path}")
                    loaded = True
                    break
            
            if not loaded:
                # Try loading from current directory (dotenv default behavior)
                if load_dotenv():
                    logger.info("Loaded .env file from current directory")
                else:
                    logger.debug("No .env file found in common locations")
        else:
            logger.warning("python-dotenv not installed - .env file support unavailable")
        
        self.mapbox_token = os.getenv("MAPBOX_ACCESS_TOKEN", "")
        if not self.mapbox_token:
            logger.error("MAPBOX_ACCESS_TOKEN not set - imagery collection will fail")
            logger.info("Please set MAPBOX_ACCESS_TOKEN in .env file or environment variable")
            logger.info(f"Looking for .env in: {Path(__file__).parent.parent.parent}, {Path.cwd()}")
        
        self.cache_dir = cache_dir
        if cache_dir:
            Path(cache_dir).mkdir(parents=True, exist_ok=True)
    
    def download_imagery(
        self,
        lat: float,
        lon: float,
        radius_m: float,
        zoom: Optional[int] = None
    ) -> Tuple[Optional[Image.Image], Dict[str, Any]]:
        """
        Download and merge Mapbox satellite imagery for the specified area.
        
        Args:
            lat: Center latitude
            lon: Center longitude
            radius_m: Radius in meters to cover
            zoom: Zoom level (defaults to MAX_ZOOM for best resolution)
        
        Returns:
            Tuple of (merged_image, metadata_dict)
            merged_image: PIL Image of merged tiles, or None if failed
            metadata_dict: Contains bounds, zoom, tile_count, etc.
        """
        if not self.mapbox_token:
            logger.error("MAPBOX_ACCESS_TOKEN not set")
            return None, {}
        
        if not HAS_PIL or not HAS_MERCANTILE:
            logger.error("PIL and mercantile required for imagery collection")
            return None, {}
        
        zoom = zoom or self.config.api.mapbox_max_zoom
        
        logger.info(f"Downloading Mapbox satellite imagery for ({lat}, {lon}), radius={radius_m}m, zoom={zoom}")
        
        # Calculate bounding box in degrees
        # Approximate: 1 degree lat ≈ 111km, 1 degree lon ≈ 111km * cos(lat)
        meters_per_degree_lat = 111000
        meters_per_degree_lon = 111000 * math.cos(math.radians(lat))
        
        radius_deg_lat = radius_m / meters_per_degree_lat
        radius_deg_lon = radius_m / meters_per_degree_lon
        
        bbox = {
            "west": lon - radius_deg_lon,
            "east": lon + radius_deg_lon,
            "south": lat - radius_deg_lat,
            "north": lat + radius_deg_lat
        }
        
        # Get all tiles that cover this bounding box
        tiles = list(mercantile.tiles(bbox["west"], bbox["south"], bbox["east"], bbox["north"], zoom))
        
        if not tiles:
            logger.error("No tiles found for bounding box")
            return None, {}
        
        logger.info(f"Found {len(tiles)} tiles to download")
        
        # Download tiles
        tile_images = []
        tile_bounds = []
        
        for tile in tiles:
            tile_img, tile_bounds_data = self._download_tile(tile, zoom)
            if tile_img:
                tile_images.append((tile, tile_img))
                tile_bounds.append(tile_bounds_data)
        
        if not tile_images:
            logger.error("Failed to download any tiles")
            return None, {}
        
        logger.info(f"Downloaded {len(tile_images)} tiles successfully")
        
        # Merge tiles into single image
        merged_image = self._merge_tiles(tile_images, bbox, zoom)
        
        if merged_image is None:
            logger.error("Failed to merge tiles")
            return None, {}
        
        # Save to cache if cache_dir is set
        if self.cache_dir and merged_image:
            cache_path = self._save_to_cache(merged_image, lat, lon, radius_m, zoom)
            logger.info(f"Cached merged image to {cache_path}")
        
        metadata = {
            "bounds": bbox,
            "zoom": zoom,
            "tile_count": len(tile_images),
            "image_size": merged_image.size if merged_image else None,
            "center": [lon, lat],
            "radius_m": radius_m
        }
        
        return merged_image, metadata
    
    def _download_tile(self, tile: mercantile.Tile, zoom: int) -> Tuple[Optional[Image.Image], Dict[str, Any]]:
        """
        Download a single Mapbox tile.
        
        Args:
            tile: Mercantile tile object
            zoom: Zoom level
        
        Returns:
            Tuple of (PIL Image, bounds dict)
        """
        url = self.MAPBOX_SATELLITE_URL.format(
            z=zoom, x=tile.x, y=tile.y, token=self.mapbox_token
        )
        
        # Check cache first
        if self.cache_dir:
            cache_path = Path(self.cache_dir) / f"tile_{zoom}_{tile.x}_{tile.y}.jpg"
            if cache_path.exists():
                try:
                    img = Image.open(cache_path)
                    bounds = mercantile.bounds(tile)
                    return img, {
                        "west": bounds.west,
                        "east": bounds.east,
                        "south": bounds.south,
                        "north": bounds.north
                    }
                except Exception as e:
                    logger.warning(f"Failed to load cached tile {cache_path}: {e}")
        
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            img = Image.open(BytesIO(response.content))
            
            # Save to cache
            if self.cache_dir:
                cache_path = Path(self.cache_dir) / f"tile_{zoom}_{tile.x}_{tile.y}.jpg"
                img.save(cache_path, quality=95)
            
            bounds = mercantile.bounds(tile)
            return img, {
                "west": bounds.west,
                "east": bounds.east,
                "south": bounds.south,
                "north": bounds.north
            }
        except Exception as e:
            logger.error(f"Failed to download tile {tile.x}/{tile.y}/{zoom}: {e}")
            return None, {}
    
    def _merge_tiles(
        self,
        tile_images: List[Tuple[mercantile.Tile, Image.Image]],
        bbox: Dict[str, float],
        zoom: int
    ) -> Optional[Image.Image]:
        """
        Merge multiple tiles into a single image.
        
        Args:
            tile_images: List of (tile, image) tuples
            bbox: Target bounding box in degrees
            zoom: Zoom level
        
        Returns:
            Merged PIL Image
        """
        if not tile_images:
            return None
        
        # Each tile is 512x512 pixels (Mapbox standard)
        TILE_SIZE = 512
        
        # Sort tiles by y (north to south) then x (west to east) for correct ordering
        tile_images.sort(key=lambda t: (t[0].y, t[0].x))
        
        # Find the tile grid bounds
        min_x = min(tile.x for tile, _ in tile_images)
        max_x = max(tile.x for tile, _ in tile_images)
        min_y = min(tile.y for tile, _ in tile_images)
        max_y = max(tile.y for tile, _ in tile_images)
        
        # Calculate merged image dimensions based on tile grid
        num_tiles_x = max_x - min_x + 1
        num_tiles_y = max_y - min_y + 1
        
        merged_width = num_tiles_x * TILE_SIZE
        merged_height = num_tiles_y * TILE_SIZE
        
        logger.debug(f"Merging {len(tile_images)} tiles into {num_tiles_x}x{num_tiles_y} grid ({merged_width}x{merged_height}px)")
        
        # Create blank canvas
        merged = Image.new("RGB", (merged_width, merged_height), color=(0, 0, 0))
        
        # Place each tile at its correct grid position
        for tile, img in tile_images:
            # Calculate position in merged image based on tile grid coordinates
            # x increases from west to east, y increases from north to south
            grid_x = tile.x - min_x
            grid_y = tile.y - min_y
            
            x_offset = grid_x * TILE_SIZE
            y_offset = grid_y * TILE_SIZE
            
            # Ensure tile is exactly TILE_SIZE x TILE_SIZE
            if img.size != (TILE_SIZE, TILE_SIZE):
                logger.warning(f"Tile {tile.x}/{tile.y} has size {img.size}, resizing to {TILE_SIZE}x{TILE_SIZE}")
                img = img.resize((TILE_SIZE, TILE_SIZE), Image.Resampling.LANCZOS)
            
            # Paste tile onto merged image
            merged.paste(img, (x_offset, y_offset))
            logger.debug(f"Placed tile {tile.x}/{tile.y} at ({x_offset}, {y_offset})")
        
        # Now crop to the exact bounding box
        # Calculate pixel positions of bbox corners within the merged image
        def lon_to_pixel_in_merged(lon):
            """Convert longitude to pixel X in merged image"""
            # Find which tile contains this longitude
            tile = mercantile.tile(lon, (bbox["north"] + bbox["south"]) / 2, zoom)
            if tile.x < min_x or tile.x > max_x:
                # Outside our tile grid, clamp to edge
                return 0 if tile.x < min_x else merged_width
            
            # Get tile bounds
            tile_bounds = mercantile.bounds(tile)
            
            # Calculate pixel within the tile
            pixel_in_tile = int(((lon - tile_bounds.west) / (tile_bounds.east - tile_bounds.west)) * TILE_SIZE)
            
            # Add tile offset
            return (tile.x - min_x) * TILE_SIZE + pixel_in_tile
        
        def lat_to_pixel_in_merged(lat):
            """Convert latitude to pixel Y in merged image (Y increases downward)"""
            # Find which tile contains this latitude
            tile = mercantile.tile((bbox["west"] + bbox["east"]) / 2, lat, zoom)
            if tile.y < min_y or tile.y > max_y:
                # Outside our tile grid, clamp to edge
                return 0 if tile.y < min_y else merged_height
            
            # Get tile bounds
            tile_bounds = mercantile.bounds(tile)
            
            # Calculate pixel within the tile (Y increases downward, so north is top)
            pixel_in_tile = int(((tile_bounds.north - lat) / (tile_bounds.north - tile_bounds.south)) * TILE_SIZE)
            
            # Add tile offset
            return (tile.y - min_y) * TILE_SIZE + pixel_in_tile
        
        # Calculate crop region
        crop_left = lon_to_pixel_in_merged(bbox["west"])
        crop_right = lon_to_pixel_in_merged(bbox["east"])
        crop_top = lat_to_pixel_in_merged(bbox["north"])
        crop_bottom = lat_to_pixel_in_merged(bbox["south"])
        
        # Ensure crop coordinates are within bounds
        crop_left = max(0, min(crop_left, merged_width - 1))
        crop_right = max(crop_left + 1, min(crop_right, merged_width))
        crop_top = max(0, min(crop_top, merged_height - 1))
        crop_bottom = max(crop_top + 1, min(crop_bottom, merged_height))
        
        # Crop to exact bounding box
        if crop_right > crop_left and crop_bottom > crop_top:
            merged = merged.crop((crop_left, crop_top, crop_right, crop_bottom))
            logger.debug(f"Cropped to ({crop_left}, {crop_top}, {crop_right}, {crop_bottom}) = {merged.size}")
        
        return merged
    
    def _save_to_cache(
        self,
        image: Image.Image,
        lat: float,
        lon: float,
        radius_m: float,
        zoom: int
    ) -> Path:
        """Save merged image to cache directory."""
        cache_path = Path(self.cache_dir) / f"mapbox_merged_{lat:.6f}_{lon:.6f}_{radius_m}m_z{zoom}.jpg"
        image.save(cache_path, quality=95)
        return cache_path

