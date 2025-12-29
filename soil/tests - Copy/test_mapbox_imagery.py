#!/usr/bin/env python
"""
Test script for Mapbox imagery download and merging

Tests both MapboxImageryCollector and TerrainCollector to download
satellite tiles and create merged images for a given location.

Usage:
    python tests/test_mapbox_imagery.py --lat 51.268535 --lon -0.570979
    python tests/test_mapbox_imagery.py --lat 51.268535 --lon -0.570979 --radius 50
"""

import os
import sys
import argparse
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from loguru import logger
from src.collectors.mapbox_imagery_collector import MapboxImageryCollector
from src.collectors.terrain_collector import TerrainCollector


def setup_logging(verbose: bool = False):
    """Configure logging"""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stderr,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level=level
    )


def test_mapbox_imagery(lat: float, lon: float, radius_m: float = 50.0, verbose: bool = False):
    """
    Test Mapbox imagery download and merging.
    
    Args:
        lat: Latitude
        lon: Longitude
        radius_m: Radius in meters (default: 50m)
        verbose: Enable verbose logging
    """
    setup_logging(verbose)
    
    logger.info("=" * 60)
    logger.info("Mapbox Imagery Download Test")
    logger.info("=" * 60)
    logger.info(f"Location: ({lat}, {lon})")
    logger.info(f"Radius: {radius_m}m")
    
    # Create output directory for test results
    test_output_dir = project_root / "tests" / "output"
    test_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Test 1: Mapbox Imagery Collector
    logger.info("\n" + "-" * 60)
    logger.info("Test 1: Mapbox Satellite Imagery Download")
    logger.info("-" * 60)
    
    try:
        imagery_collector = MapboxImageryCollector(cache_dir=str(test_output_dir))
        
        if not imagery_collector.mapbox_token:
            logger.error("MAPBOX_ACCESS_TOKEN not set!")
            logger.info("Please set MAPBOX_ACCESS_TOKEN in .env file or environment variable")
            return 1
        
        logger.info(f"Mapbox token found: {imagery_collector.mapbox_token[:10]}...")
        
        # Download imagery
        merged_image, metadata = imagery_collector.download_imagery(
            lat=lat,
            lon=lon,
            radius_m=radius_m
        )
        
        if merged_image is None:
            logger.error("Failed to download/merge imagery")
            return 1
        
        logger.info(f"✓ Successfully downloaded and merged imagery")
        logger.info(f"  Image size: {merged_image.size[0]}x{merged_image.size[1]} pixels")
        logger.info(f"  Zoom level: {metadata.get('zoom', 'unknown')}")
        logger.info(f"  Tile count: {metadata.get('tile_count', 'unknown')}")
        logger.info(f"  Bounds: {metadata.get('bounds', {})}")
        
        # Note: Image is already saved by MapboxImageryCollector._save_to_cache()
        # Check for the cached merged image
        cached_merged = list(test_output_dir.glob(f"mapbox_merged_{lat:.6f}_{lon:.6f}_{radius_m}m_z*.jpg"))
        if cached_merged:
            logger.info(f"✓ Merged image cached: {cached_merged[0].name}")
        
        # List downloaded tiles
        tile_files = list(test_output_dir.glob("tile_*.jpg"))
        if tile_files:
            logger.info(f"\nDownloaded {len(tile_files)} tile files:")
            for tile_file in sorted(tile_files):
                logger.info(f"  - {tile_file.name}")
        
    except Exception as e:
        logger.error(f"Error in imagery download: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    # Test 2: Terrain Collector (uses Mapbox Terrain-RGB)
    logger.info("\n" + "-" * 60)
    logger.info("Test 2: Mapbox Terrain-RGB Download")
    logger.info("-" * 60)
    
    try:
        terrain_collector = TerrainCollector()
        
        if not terrain_collector.mapbox_token:
            logger.warning("MAPBOX_ACCESS_TOKEN not set for terrain collector")
        else:
            logger.info(f"Mapbox token found: {terrain_collector.mapbox_token[:10]}...")
            
            terrain_data = terrain_collector.collect_terrain(
                lat=lat,
                lon=lon,
                radius_m=radius_m
            )
            
            logger.info(f"✓ Terrain data collected")
            logger.info(f"  Elevation: {terrain_data.elevation_m:.2f}m")
            logger.info(f"  Slope: {terrain_data.slope_percent:.2f}%")
            logger.info(f"  Slope direction: {terrain_data.slope_direction}")
            logger.info(f"  Terrain classification: {terrain_data.terrain_classification}")
            logger.info(f"  DEM source: {terrain_data.dem_source}")
            logger.info(f"  Resolution: {terrain_data.resolution_m:.2f}m")
            
    except Exception as e:
        logger.warning(f"Error in terrain collection: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
    
    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Test Summary")
    logger.info("=" * 60)
    cached_merged = list(test_output_dir.glob(f"mapbox_merged_{lat:.6f}_{lon:.6f}_{radius_m}m_z*.jpg"))
    if cached_merged:
        logger.info(f"✓ Merged image: {cached_merged[0].name}")
    logger.info(f"✓ Tile files: {len(tile_files)} files in {test_output_dir}")
    logger.info(f"✓ Test completed successfully!")
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Test Mapbox imagery download and merging",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tests/test_mapbox_imagery.py --lat 51.268535 --lon -0.570979
  python tests/test_mapbox_imagery.py --lat 51.268535 --lon -0.570979 --radius 50 --verbose
        """
    )
    
    parser.add_argument(
        "--lat",
        type=float,
        required=True,
        help="Latitude of test location"
    )
    
    parser.add_argument(
        "--lon",
        type=float,
        required=True,
        help="Longitude of test location"
    )
    
    parser.add_argument(
        "--radius",
        type=float,
        default=50.0,
        help="Radius in meters (default: 50m)"
    )
    
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    return test_mapbox_imagery(
        lat=args.lat,
        lon=args.lon,
        radius_m=args.radius,
        verbose=args.verbose
    )


if __name__ == "__main__":
    sys.exit(main())

