#!/usr/bin/env python
"""
Test script for Boundary Collector

Tests the BoundaryCollector to detect plot boundaries using various strategies:
0. INSPIRE cadastral data (automatically used if available in data/inspires/)
1. OSM landuse polygons
2. Road-derived boundaries
3. Neighbor building boundaries
4. Building footprint with buffer
5. Default estimate

The collector automatically tries INSPIRE GML first (if found in data/inspires/),
then falls back to OSM strategies if no INSPIRE boundary is found.

Usage:
    python tests/test_boundary_collector.py --lat 51.268535 --lon -0.570979
    python tests/test_boundary_collector.py --lat 51.268535 --lon -0.570979 --radius 50 --verbose
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from loguru import logger
from src.collectors.boundary import BoundaryCollector
from src.collectors.osm import OSMCollector


def setup_logging(verbose: bool = False):
    """Configure logging"""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stderr,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level=level
    )


def test_boundary_collector(
    lat: float,
    lon: float,
    radius_m: float = 30.0,
    verbose: bool = False,
    with_osm_data: bool = True
):
    """
    Test Boundary Collector with various strategies.
    
    Args:
        lat: Latitude
        lon: Longitude
        radius_m: Search radius in meters (default: 30m)
        verbose: Enable verbose logging
        with_osm_data: Whether to fetch OSM data first (default: True)
    """
    setup_logging(verbose)
    
    logger.info("=" * 60)
    logger.info("Boundary Collector Test")
    logger.info("=" * 60)
    logger.info(f"Location: ({lat}, {lon})")
    logger.info(f"Search radius: {radius_m}m")
    logger.info(f"Using OSM data: {with_osm_data}")
    
    # Create output directory for test results
    test_output_dir = project_root / "tests" / "output"
    test_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize collectors (BoundaryCollector will automatically try INSPIRE GML if available)
    boundary_collector = BoundaryCollector()
    if boundary_collector.inspire_gml_path:
        logger.info(f"Using INSPIRE GML: {boundary_collector.inspire_gml_path}")
    else:
        logger.info("No INSPIRE GML found - will use OSM data")
    osm_collector = OSMCollector() if with_osm_data else None
    
    # Test 1: Get boundary with OSM data (as used in pipeline)
    logger.info("\n" + "-" * 60)
    logger.info("Test 1: Get Plot Boundary (with OSM data)")
    logger.info("-" * 60)
    
    try:
        osm_buildings = None
        osm_roads = None
        osm_landuse = None
        
        if with_osm_data and osm_collector:
            logger.info("Fetching OSM data first...")
            all_osm_features = osm_collector.fetch_all_features(lat, lon, radius_m)
            osm_buildings = all_osm_features.get("buildings", [])
            osm_roads = all_osm_features.get("roads", [])
            osm_landuse = all_osm_features.get("landuse", [])
            
            logger.info(f"  Found {len(osm_buildings)} buildings")
            logger.info(f"  Found {len(osm_roads)} roads")
            logger.info(f"  Found {len(osm_landuse)} landuse polygons")
        
        # Get boundary
        boundary_data = boundary_collector.get_plot_boundary(
            lat=lat,
            lon=lon,
            search_radius_m=radius_m,
            osm_buildings=osm_buildings,
            osm_roads=osm_roads,
            osm_landuse=osm_landuse
        )
        
        if boundary_data:
            logger.info("✓ Successfully retrieved boundary")
            logger.info(f"  Source: {boundary_data.get('source', 'unknown')}")
            logger.info(f"  Type: {boundary_data.get('type', 'unknown')}")
            logger.info(f"  Area: {boundary_data.get('area_sqm', 0):.2f} m²")
            logger.info(f"  Perimeter: {boundary_data.get('perimeter_m', 0):.2f} m")
            logger.info(f"  Accuracy: {boundary_data.get('accuracy_m', 0):.2f} m")
            
            coords = boundary_data.get("coordinates", [[]])[0]
            if coords:
                logger.info(f"  Coordinates: {len(coords)} points")
                logger.info(f"  First point: [{coords[0][0]:.6f}, {coords[0][1]:.6f}]")
                logger.info(f"  Last point: [{coords[-1][0]:.6f}, {coords[-1][1]:.6f}]")
            
            # Save boundary to JSON file
            output_file = test_output_dir / f"boundary_{lat:.6f}_{lon:.6f}.json"
            with open(output_file, 'w') as f:
                json.dump(boundary_data, f, indent=2)
            logger.info(f"  Saved to: {output_file}")
        else:
            logger.error("Failed to retrieve boundary")
            return 1
            
    except Exception as e:
        logger.error(f"Error in boundary collection: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    # Test 2: Get boundary without OSM data (direct API call)
    if with_osm_data:
        logger.info("\n" + "-" * 60)
        logger.info("Test 2: Get Plot Boundary (without pre-fetched OSM data)")
        logger.info("-" * 60)
        
        try:
            boundary_data_direct = boundary_collector.get_plot_boundary(
                lat=lat,
                lon=lon,
                search_radius_m=radius_m,
                osm_buildings=None,
                osm_roads=None,
                osm_landuse=None
            )
            
            if boundary_data_direct:
                logger.info("✓ Successfully retrieved boundary (direct API call)")
                logger.info(f"  Source: {boundary_data_direct.get('source', 'unknown')}")
                logger.info(f"  Area: {boundary_data_direct.get('area_sqm', 0):.2f} m²")
                logger.info(f"  Perimeter: {boundary_data_direct.get('perimeter_m', 0):.2f} m")
            else:
                logger.warning("Failed to retrieve boundary (direct API call)")
                
        except Exception as e:
            logger.warning(f"Error in direct boundary collection: {e}")
            if verbose:
                import traceback
                traceback.print_exc()
    
    # Test 3: Test individual strategies (if verbose)
    if verbose:
        logger.info("\n" + "-" * 60)
        logger.info("Test 3: Individual Strategy Tests")
        logger.info("-" * 60)
        
        # Test landuse boundary
        if osm_landuse:
            logger.info("\nTesting Strategy 1: Landuse polygons...")
            for i, landuse in enumerate(osm_landuse[:3]):  # Test first 3
                geom = landuse.get("geometry", {})
                if geom.get("type") == "Polygon":
                    coords = geom.get("coordinates", [[]])[0]
                    if coords and len(coords) >= 4:
                        logger.info(f"  Landuse {i+1}: {len(coords)} points, OSM ID: {landuse.get('osm_id', 'unknown')}")
        
        # Test building-based boundary
        if osm_buildings:
            logger.info("\nTesting Strategy 4: Building footprint...")
            building_boundary = boundary_collector._estimate_from_building(
                lat, lon, radius_m, osm_buildings
            )
            if building_boundary:
                logger.info(f"  Building-based boundary found: {building_boundary.get('area_sqm', 0):.2f} m²")
            else:
                logger.info("  No building-based boundary found")
        
        # Test road-derived boundary
        if osm_roads:
            logger.info("\nTesting Strategy 2: Road-derived boundary...")
            road_boundary = boundary_collector._derive_from_roads(
                lat, lon, osm_roads, radius_m
            )
            if road_boundary:
                logger.info(f"  Road-derived boundary found: {road_boundary.get('area_sqm', 0):.2f} m²")
            else:
                logger.info("  No road-derived boundary found")
        
        # Test neighbor building boundary
        if osm_buildings:
            logger.info("\nTesting Strategy 3: Neighbor building boundary...")
            neighbor_boundary = boundary_collector._derive_from_neighbor_buildings(
                lat, lon, osm_buildings, radius_m
            )
            if neighbor_boundary:
                logger.info(f"  Neighbor-derived boundary found: {neighbor_boundary.get('area_sqm', 0):.2f} m²")
            else:
                logger.info("  No neighbor-derived boundary found")
    
    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Test Summary")
    logger.info("=" * 60)
    logger.info(f"✓ Boundary collected successfully")
    logger.info(f"  Source: {boundary_data.get('source', 'unknown')}")
    logger.info(f"  Area: {boundary_data.get('area_sqm', 0):.2f} m²")
    logger.info(f"  Perimeter: {boundary_data.get('perimeter_m', 0):.2f} m")
    output_file = test_output_dir / f"boundary_{lat:.6f}_{lon:.6f}.json"
    if output_file.exists():
        logger.info(f"  Output saved: {output_file}")
    logger.info(f"✓ Test completed successfully!")
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Test Boundary Collector for plot boundary detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tests/test_boundary_collector.py --lat 51.268535 --lon -0.570979
  python tests/test_boundary_collector.py --lat 51.268535 --lon -0.570979 --radius 50 --verbose
  python tests/test_boundary_collector.py --lat 51.268535 --lon -0.570979 --no-osm-data
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
        default=30.0,
        help="Search radius in meters (default: 30m)"
    )
    
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging and individual strategy tests"
    )
    
    parser.add_argument(
        "--no-osm-data",
        action="store_true",
        help="Skip OSM data pre-fetching (test direct API calls only)"
    )
    
    args = parser.parse_args()
    
    return test_boundary_collector(
        lat=args.lat,
        lon=args.lon,
        radius_m=args.radius,
        verbose=args.verbose,
        with_osm_data=not args.no_osm_data
    )


if __name__ == "__main__":
    sys.exit(main())

