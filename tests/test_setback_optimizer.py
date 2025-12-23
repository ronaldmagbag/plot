#!/usr/bin/env python
"""
Test script for Setback Optimizer

Tests the SetbackOptimizer to calculate optimal rectangular setbacks
inside property boundaries based on rules and scoring.

Usage:
    python tests/test_setback_optimizer.py --lat 51.268535 --lon -0.570979
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
from src.analysis.setback_optimizer import SetbackOptimizer
from src.collectors.boundary_collector import BoundaryCollector
from src.collectors.osm_collector import OSMCollector


def setup_logging(verbose: bool = False):
    """Configure logging"""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stderr,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level=level
    )


def test_setback_optimizer(
    lat: float,
    lon: float,
    radius_m: float = 50.0,
    verbose: bool = False
):
    """
    Test Setback Optimizer with property boundary, roads, and buildings.
    
    Args:
        lat: Latitude
        lon: Longitude
        radius_m: Search radius in meters (default: 50m)
        verbose: Enable verbose logging
    """
    setup_logging(verbose)
    
    logger.info("=" * 60)
    logger.info("Setback Optimizer Test")
    logger.info("=" * 60)
    logger.info(f"Location: ({lat}, {lon})")
    logger.info(f"Search radius: {radius_m}m")
    
    # Create output directory for test results
    test_output_dir = project_root / "tests" / "output"
    test_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize collectors
    boundary_collector = BoundaryCollector()
    osm_collector = OSMCollector()
    setback_optimizer = SetbackOptimizer()
    
    # Get property boundary
    logger.info("\n" + "-" * 60)
    logger.info("Step 1: Getting property boundary...")
    logger.info("-" * 60)
    
    all_osm_features = osm_collector.fetch_all_features(lat, lon, radius_m)
    osm_buildings = all_osm_features.get("buildings", [])
    osm_roads = all_osm_features.get("roads", [])
    osm_landuse = all_osm_features.get("landuse", [])
    
    boundary_data = boundary_collector.get_plot_boundary(
        lat=lat,
        lon=lon,
        search_radius_m=radius_m,
        osm_buildings=osm_buildings,
        osm_roads=osm_roads,
        osm_landuse=osm_landuse
    )
    
    if not boundary_data:
        logger.error("Failed to get property boundary")
        return 1
    
    property_coords = boundary_data.get("coordinates", [[]])[0]
    logger.info(f"Property boundary source: {boundary_data.get('source', 'unknown')}")
    logger.info(f"Property area: {boundary_data.get('area_sqm', 0):.2f} m²")
    logger.info(f"Property coordinates: {len(property_coords)} points")
    
    # Calculate optimal setback
    logger.info("\n" + "-" * 60)
    logger.info("Step 2: Calculating optimal setback...")
    logger.info("-" * 60)
    
    setback_result = setback_optimizer.calculate_optimal_setback(
        property_boundary=property_coords,
        roads=osm_roads,
        neighbor_buildings=osm_buildings
    )
    
    if not setback_result:
        logger.error("Failed to calculate setback")
        return 1
    
    logger.info("✓ Setback calculated successfully")
    logger.info(f"  Setback area: {setback_result.get('area_sqm', 0):.2f} m²")
    logger.info(f"  Setback source: {setback_result.get('derived_from', 'unknown')}")
    
    # Display scoring
    scoring = setback_result.get("scoring", {})
    logger.info("\nScoring breakdown:")
    logger.info(f"  Total score: {scoring.get('total_score', 0):.2f}")
    logger.info(f"  Inside property: {scoring.get('inside_property_score', 0):.2f}")
    logger.info(f"  Rectangle shape: {scoring.get('rectangle_score', 0):.2f}")
    logger.info(f"  Alignment: {scoring.get('alignment_score', 0):.2f}")
    logger.info(f"  Distance compliance: {scoring.get('distance_compliance_score', 0):.2f}")
    logger.info(f"  Area: {scoring.get('area_score', 0):.2f}")
    
    # Display setbacks applied
    setbacks = setback_result.get("setbacks_applied", {})
    logger.info("\nSetbacks applied:")
    logger.info(f"  Front: {setbacks.get('front_m', 0):.2f} m")
    logger.info(f"  Rear: {setbacks.get('rear_m', 0):.2f} m")
    logger.info(f"  Side East: {setbacks.get('side_east_m', 0):.2f} m")
    logger.info(f"  Side West: {setbacks.get('side_west_m', 0):.2f} m")
    
    # Save results
    output_file = test_output_dir / f"setback_{lat:.6f}_{lon:.6f}.json"
    with open(output_file, 'w') as f:
        json.dump({
            "property_boundary": boundary_data,
            "setback": setback_result
        }, f, indent=2)
    logger.info(f"\n✓ Results saved to: {output_file}")
    
    # Verify setback is inside property
    logger.info("\n" + "-" * 60)
    logger.info("Step 3: Verifying setback rules...")
    logger.info("-" * 60)
    
    try:
        from shapely.geometry import Polygon, Point
        
        property_poly = Polygon(property_coords)
        setback_coords = setback_result.get("coordinates", [[]])[0]
        setback_poly = Polygon(setback_coords)
        
        is_inside = property_poly.contains(setback_poly)
        logger.info(f"  Rule 1 (Inside property): {'✓ PASS' if is_inside else '✗ FAIL'}")
        
        # Check if it's a rectangle (4 corners)
        is_rectangle = len(setback_coords) == 5  # 4 corners + closing point
        logger.info(f"  Rule 2 (Rectangle shape): {'✓ PASS' if is_rectangle else '✗ FAIL'}")
        
        if is_rectangle:
            # Check if angles are close to 90 degrees
            coords = setback_coords[:-1]  # Remove closing point
            angles_ok = True
            for i in range(4):
                p1 = coords[i]
                p2 = coords[(i + 1) % 4]
                p3 = coords[(i + 2) % 4]
                
                # Calculate angle
                v1 = (p2[0] - p1[0], p2[1] - p1[1])
                v2 = (p3[0] - p2[0], p3[1] - p2[1])
                
                import math
                dot = v1[0] * v2[0] + v1[1] * v2[1]
                len1 = math.sqrt(v1[0]**2 + v1[1]**2)
                len2 = math.sqrt(v2[0]**2 + v2[1]**2)
                
                if len1 > 0 and len2 > 0:
                    cos_angle = dot / (len1 * len2)
                    cos_angle = max(-1, min(1, cos_angle))
                    angle = math.acos(cos_angle)
                    angle_deg = math.degrees(angle)
                    
                    if abs(angle_deg - 90) > 10:  # Allow 10 degree tolerance
                        angles_ok = False
                        break
            
            logger.info(f"  Rule 2 (Right angles): {'✓ PASS' if angles_ok else '✗ FAIL'}")
        
    except Exception as e:
        logger.warning(f"Could not verify rules: {e}")
    
    logger.info("\n" + "=" * 60)
    logger.info("Test Summary")
    logger.info("=" * 60)
    logger.info(f"✓ Property boundary: {boundary_data.get('source', 'unknown')}")
    logger.info(f"✓ Setback calculated: {setback_result.get('area_sqm', 0):.2f} m²")
    logger.info(f"✓ Score: {scoring.get('total_score', 0):.2f}")
    logger.info(f"✓ Test completed successfully!")
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Test Setback Optimizer for optimal rectangular setbacks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tests/test_setback_optimizer.py --lat 51.268535 --lon -0.570979
  python tests/test_setback_optimizer.py --lat 51.268535 --lon -0.570979 --radius 50 --verbose
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
        help="Search radius in meters (default: 50m)"
    )
    
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    return test_setback_optimizer(
        lat=args.lat,
        lon=args.lon,
        radius_m=args.radius,
        verbose=args.verbose
    )


if __name__ == "__main__":
    sys.exit(main())

