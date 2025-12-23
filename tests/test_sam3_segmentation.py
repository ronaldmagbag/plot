#!/usr/bin/env python
"""
Test script for SAM3 segmentation

Tests SAM3 segmentation on Mapbox satellite imagery to segment
roads, trees, and grasses for a given location.

Usage:
    python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979
    python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979 --radius 50
    python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979 --image path/to/image.jpg
    
Note: The SAM3 segmentation script is now located at src/segmentation/sam3_segment.py
"""

import os
import sys
import json
import argparse
import subprocess
from pathlib import Path
from typing import Optional

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from loguru import logger
from src.collectors.mapbox_imagery_collector import MapboxImageryCollector


def setup_logging(verbose: bool = False):
    """Configure logging"""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stderr,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level=level
    )


def test_sam3_segmentation(
    lat: float,
    lon: float,
    radius_m: float = 50.0,
    image_path: Optional[str] = None,
    verbose: bool = False
):
    """
    Test SAM3 segmentation on satellite imagery.
    
    Args:
        lat: Latitude
        lon: Longitude
        radius_m: Radius in meters (default: 50m)
        image_path: Optional path to existing image (skips download if provided)
        verbose: Enable verbose logging
    """
    setup_logging(verbose)
    
    logger.info("=" * 60)
    logger.info("SAM3 Segmentation Test")
    logger.info("=" * 60)
    logger.info(f"Location: ({lat}, {lon})")
    logger.info(f"Radius: {radius_m}m")
    
    # Create output directory for test results
    test_output_dir = project_root / "tests" / "output"
    test_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Get or download imagery
    logger.info("\n" + "-" * 60)
    logger.info("Step 1: Preparing Satellite Imagery")
    logger.info("-" * 60)
    
    if image_path and Path(image_path).exists():
        logger.info(f"Using provided image: {image_path}")
        image_file = Path(image_path)
        # Create bbox from image metadata or use default
        bbox = {
            "west": lon - (radius_m / 111000),
            "east": lon + (radius_m / 111000),
            "south": lat - (radius_m / 111000),
            "north": lat + (radius_m / 111000)
        }
    else:
        logger.info("Downloading Mapbox satellite imagery...")
        try:
            imagery_collector = MapboxImageryCollector(cache_dir=str(test_output_dir))
            
            if not imagery_collector.mapbox_token:
                logger.error("MAPBOX_ACCESS_TOKEN not set!")
                logger.info("Please set MAPBOX_ACCESS_TOKEN in .env file or environment variable")
                return 1
            
            # Download imagery
            merged_image, metadata = imagery_collector.download_imagery(
                lat=lat,
                lon=lon,
                radius_m=radius_m
            )
            
            if merged_image is None:
                logger.error("Failed to download/merge imagery")
                return 1
            
            logger.info(f"✓ Successfully downloaded imagery")
            logger.info(f"  Image size: {merged_image.size[0]}x{merged_image.size[1]} pixels")
            logger.info(f"  Zoom level: {metadata.get('zoom', 'unknown')}")
            
            # Find cached image
            zoom = metadata.get('zoom', 20)
            cached_images = list(test_output_dir.glob(f"mapbox_merged_{lat:.6f}_{lon:.6f}_{radius_m}m_z{zoom}.jpg"))
            if cached_images:
                image_file = cached_images[0]
                bbox = metadata.get("bounds", {})
                logger.info(f"✓ Using cached image: {image_file.name}")
            else:
                logger.error("Cached image not found")
                return 1
                
        except Exception as e:
            logger.error(f"Error downloading imagery: {e}")
            if verbose:
                import traceback
                traceback.print_exc()
            return 1
    
    # Step 2: Run SAM3 segmentation
    logger.info("\n" + "-" * 60)
    logger.info("Step 2: Running SAM3 Segmentation")
    logger.info("-" * 60)
    
    try:
        # Find SAM3 script
        script_path = project_root / "src" / "segmentation" / "sam3_segment.py"
        if not script_path.exists():
            logger.error(f"SAM3 script not found at {script_path}")
            return 1
        
        # Prepare output path
        output_path = test_output_dir / f"sam3_results_{lat:.6f}_{lon:.6f}.json"
        
        # Prepare bbox JSON
        bbox_json = json.dumps(bbox)
        
        logger.info("Running SAM3 segmentation in cu126 environment...")
        logger.info(f"  Input image: {image_file}")
        logger.info(f"  Output: {output_path}")
        
        # Find micromamba
        import shutil
        micromamba_cmd = shutil.which("micromamba")
        if not micromamba_cmd:
            logger.warning("micromamba not found in PATH, trying direct command...")
            micromamba_cmd = "micromamba"
        
        # Run SAM3 segmentation
        cmd = [
            micromamba_cmd, "run", "-n", "cu126",
            "python", str(script_path),
            "--image", str(image_file),
            "--output", str(output_path),
            "--bbox", bbox_json,
            "--center-lon", str(lon),
            "--center-lat", str(lat)
        ]
        
        logger.debug(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode != 0:
            logger.error(f"SAM3 segmentation failed!")
            logger.error(f"Return code: {result.returncode}")
            logger.error(f"STDERR: {result.stderr}")
            if result.stdout:
                logger.error(f"STDOUT: {result.stdout}")
            return 1
        
        logger.info("✓ SAM3 segmentation completed")
        if result.stdout:
            logger.debug(f"Output: {result.stdout}")
        
    except subprocess.TimeoutExpired:
        logger.error("SAM3 segmentation timed out (5 minutes)")
        return 1
    except FileNotFoundError:
        logger.error("micromamba not found. Make sure micromamba is installed and in PATH")
        return 1
    except Exception as e:
        logger.error(f"Error running SAM3: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    # Step 3: Display results
    logger.info("\n" + "-" * 60)
    logger.info("Step 3: Segmentation Results")
    logger.info("-" * 60)
    
    if not output_path.exists():
        logger.error("SAM3 output file not found")
        return 1
    
    try:
        with open(output_path, "r") as f:
            results = json.load(f)
        
        roads = results.get("roads", [])
        trees = results.get("trees", [])
        grasses = results.get("grasses", [])
        
        logger.info(f"✓ Segmentation results loaded from {output_path.name}")
        logger.info(f"  Roads found: {len(roads)}")
        logger.info(f"  Trees found: {len(trees)}")
        logger.info(f"  Grasses/vegetation found: {len(grasses)}")
        
        # Show details for each category
        if roads:
            logger.info(f"\n  Road segments:")
            for i, road in enumerate(roads[:5]):  # Show first 5
                conf = road.get("confidence", 0)
                logger.info(f"    {i+1}. {road.get('id', 'unknown')} (confidence: {conf:.2f})")
            if len(roads) > 5:
                logger.info(f"    ... and {len(roads) - 5} more")
        
        if trees:
            logger.info(f"\n  Tree segments:")
            for i, tree in enumerate(trees[:5]):  # Show first 5
                conf = tree.get("confidence", 0)
                logger.info(f"    {i+1}. {tree.get('id', 'unknown')} (confidence: {conf:.2f})")
            if len(trees) > 5:
                logger.info(f"    ... and {len(trees) - 5} more")
        
        if grasses:
            logger.info(f"\n  Grass/vegetation segments:")
            for i, grass in enumerate(grasses[:5]):  # Show first 5
                conf = grass.get("confidence", 0)
                logger.info(f"    {i+1}. {grass.get('id', 'unknown')} (confidence: {conf:.2f})")
            if len(grasses) > 5:
                logger.info(f"    ... and {len(grasses) - 5} more")
        
    except Exception as e:
        logger.error(f"Error reading results: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Test Summary")
    logger.info("=" * 60)
    logger.info(f"✓ Input image: {image_file.name}")
    logger.info(f"✓ Output JSON: {output_path.name}")
    logger.info(f"✓ Roads: {len(roads)} segments")
    logger.info(f"✓ Trees: {len(trees)} segments")
    logger.info(f"✓ Grasses: {len(grasses)} segments")
    logger.info(f"✓ Test completed successfully!")
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Test SAM3 segmentation on satellite imagery",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979
  python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979 --radius 50
  python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979 --image tests/output/mapbox_merged_*.jpg
  python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979 --verbose
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
        help="Radius in meters (default: 50m). Only used if --image not provided"
    )
    
    parser.add_argument(
        "--image",
        type=str,
        default=None,
        help="Path to existing satellite image (skips download if provided)"
    )
    
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    return test_sam3_segmentation(
        lat=args.lat,
        lon=args.lon,
        radius_m=args.radius,
        image_path=args.image,
        verbose=args.verbose
    )


if __name__ == "__main__":
    sys.exit(main())

