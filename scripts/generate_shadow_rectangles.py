#!/usr/bin/env python3
"""
Generate shadow rectangles for surrounding buildings and trees
Input: plot.json, date/time, temp folder
Output: JSON with shadow rectangles
"""

import json
import math
import argparse
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional
import sys

# Add src to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "src"))

try:
    import pvlib
    PVLIB_AVAILABLE = True
except ImportError:
    PVLIB_AVAILABLE = False
    print("Warning: pvlib not available, using simplified sun position calculations")

# Import geometry utils
from src.analysis.geometry_utils import GeometryUtils


def calculate_sun_position(
    lat: float,
    lon: float,
    dt: datetime,
    timezone: str = "Europe/London"
) -> Dict[str, float]:
    """
    Calculate sun position (azimuth and elevation) for given date/time
    
    Returns:
        Dict with 'azimuth' (degrees, 0=North, 90=East) and 'elevation' (degrees)
    """
    if PVLIB_AVAILABLE:
        import pandas as pd
        times = pd.DatetimeIndex([dt], tz=timezone)
        solar_pos = pvlib.solarposition.get_solarposition(times, lat, lon)
        azimuth = float(solar_pos['azimuth'].iloc[0])
        elevation = float(solar_pos['apparent_elevation'].iloc[0])
    else:
        # Simplified calculation
        # Day of year
        doy = dt.timetuple().tm_yday
        
        # Solar declination
        declination = 23.45 * math.sin(math.radians(360 * (284 + doy) / 365))
        
        # Hour angle (simplified - assumes solar noon at 12:00)
        hour = dt.hour + dt.minute / 60.0
        hour_angle = 15 * (hour - 12)  # degrees
        
        # Convert to radians
        lat_rad = math.radians(lat)
        dec_rad = math.radians(declination)
        hour_rad = math.radians(hour_angle)
        
        # Elevation
        sin_elevation = (
            math.sin(lat_rad) * math.sin(dec_rad) +
            math.cos(lat_rad) * math.cos(dec_rad) * math.cos(hour_rad)
        )
        elevation = math.degrees(math.asin(max(-1, min(1, sin_elevation))))
        
        # Azimuth (simplified)
        cos_azimuth = (
            (math.sin(dec_rad) - math.sin(lat_rad) * sin_elevation) /
            (math.cos(lat_rad) * math.cos(math.radians(elevation)))
        )
        azimuth_rad = math.acos(max(-1, min(1, cos_azimuth)))
        azimuth = math.degrees(azimuth_rad)
        
        # Adjust azimuth based on time of day
        if hour_angle < 0:
            azimuth = 360 - azimuth
    
    return {
        'azimuth': azimuth,  # 0=North, 90=East, 180=South, 270=West
        'elevation': elevation
    }


def calculate_shadow_rectangle(
    footprint_coords: List[List[float]],  # List of [lon, lat] coordinates
    height_m: float,
    sun_azimuth: float,  # degrees, 0=North, 90=East
    sun_elevation: float,  # degrees
    ref_lon: float,
    ref_lat: float
) -> Dict[str, Any]:
    """
    Calculate shadow rectangle (bounding box) for a building or tree
    
    Args:
        footprint_coords: Building/tree footprint coordinates [[lon, lat], ...]
        height_m: Height of building/tree in meters
        sun_azimuth: Sun azimuth angle in degrees (0=North, 90=East)
        sun_elevation: Sun elevation angle in degrees
        ref_lon, ref_lat: Reference point for local coordinate conversion
    
    Returns:
        Dict with shadow rectangle as GeoJSON Polygon
    """
    if sun_elevation <= 0:
        # Sun below horizon, no shadow
        return None
    
    # Convert to local coordinates (meters)
    local_coords = GeometryUtils.degrees_to_local(footprint_coords, ref_lon, ref_lat)
    
    # Calculate shadow length
    # shadow_length = height * cot(elevation)
    elevation_rad = math.radians(sun_elevation)
    shadow_length = height_m / math.tan(elevation_rad) if elevation_rad > 0 else 0
    
    # Calculate shadow direction vector (in local coordinates)
    # Shadow extends opposite to sun direction
    # Azimuth: 0=North, 90=East, 180=South, 270=West
    # Shadow direction = sun_azimuth + 180°
    shadow_azimuth = (sun_azimuth + 180) % 360
    
    # In local coords: North = +Y, East = +X
    shadow_azimuth_rad = math.radians(shadow_azimuth)
    shadow_dx = shadow_length * math.sin(shadow_azimuth_rad)  # East component
    shadow_dy = shadow_length * math.cos(shadow_azimuth_rad)  # North component
    
    # Project each footprint point to shadow
    shadow_points = []
    for x, y in local_coords:
        shadow_x = x + shadow_dx
        shadow_y = y + shadow_dy
        shadow_points.append([shadow_x, shadow_y])
    
    # Combine footprint and shadow to get full shadow polygon
    # Shadow extends from footprint edge to shadow edge
    # For rectangle calculation, we'll use the bounding box of footprint + shadow
    all_points = local_coords + shadow_points
    
    # Find bounding box
    min_x = min(p[0] for p in all_points)
    max_x = max(p[0] for p in all_points)
    min_y = min(p[1] for p in all_points)
    max_y = max(p[1] for p in all_points)
    
    # Convert bounding box corners back to degrees
    bbox_corners_local = [
        (min_x, min_y),
        (max_x, min_y),
        (max_x, max_y),
        (min_x, max_y),
        (min_x, min_y)  # Close polygon
    ]
    
    bbox_corners_deg = GeometryUtils.local_to_degrees(bbox_corners_local, ref_lon, ref_lat)
    
    return {
        "type": "Polygon",
        "coordinates": [bbox_corners_deg]
    }


def process_buildings(
    buildings: List[Dict[str, Any]],
    sun_azimuth: float,
    sun_elevation: float,
    plot_centroid: Tuple[float, float]
) -> List[Dict[str, Any]]:
    """Process buildings and calculate shadow rectangles"""
    shadow_rectangles = []
    ref_lon, ref_lat = plot_centroid
    
    for building in buildings:
        building_id = building.get("id", "unknown")
        footprint = building.get("footprint", {})
        height_m = building.get("height_m", 6.0)  # Default 6m if not specified
        
        if footprint.get("type") != "Polygon":
            continue
        
        coords = footprint.get("coordinates", [])
        if not coords or len(coords) == 0:
            continue
        
        # Extract coordinates (remove nesting if needed)
        footprint_coords = coords[0] if isinstance(coords[0][0], list) else coords
        
        shadow_rect = calculate_shadow_rectangle(
            footprint_coords,
            height_m,
            sun_azimuth,
            sun_elevation,
            ref_lon,
            ref_lat
        )
        
        if shadow_rect:
            shadow_rectangles.append({
                "id": building_id,
                "type": "building",
                "height_m": height_m,
                "shadow_rectangle": shadow_rect,
                "sun_azimuth": round(sun_azimuth, 2),
                "sun_elevation": round(sun_elevation, 2)
            })
    
    return shadow_rectangles


def process_trees(
    tree_zones: Dict[str, Any],
    sun_azimuth: float,
    sun_elevation: float,
    plot_centroid: Tuple[float, float]
) -> List[Dict[str, Any]]:
    """Process trees and calculate shadow rectangles"""
    shadow_rectangles = []
    ref_lon, ref_lat = plot_centroid
    
    trees = tree_zones.get("trees", [])
    if not trees:
        return shadow_rectangles
    
    for tree in trees:
        tree_id = tree.get("id", f"tree_{len(shadow_rectangles)}")
        location = tree.get("location", {})
        height_m = tree.get("height_m", 5.0)  # Default 5m for trees
        
        if location.get("type") != "Point":
            continue
        
        coords = location.get("coordinates", [])
        if not coords or len(coords) < 2:
            continue
        
        # For point trees, create a small circular footprint
        # Use tree height as radius estimate (or fixed small radius)
        tree_radius_m = min(height_m * 0.1, 2.0)  # 10% of height, max 2m
        
        # Create circular footprint (approximate as square for simplicity)
        lon, lat = coords[0], coords[1]
        
        # Convert to local, create square footprint
        center_local = GeometryUtils.degrees_to_local([[lon, lat]], ref_lon, ref_lat)[0]
        half_size = tree_radius_m
        
        footprint_local = [
            (center_local[0] - half_size, center_local[1] - half_size),
            (center_local[0] + half_size, center_local[1] - half_size),
            (center_local[0] + half_size, center_local[1] + half_size),
            (center_local[0] - half_size, center_local[1] + half_size),
            (center_local[0] - half_size, center_local[1] - half_size)
        ]
        
        footprint_coords = GeometryUtils.local_to_degrees(footprint_local, ref_lon, ref_lat)
        
        shadow_rect = calculate_shadow_rectangle(
            footprint_coords,
            height_m,
            sun_azimuth,
            sun_elevation,
            ref_lon,
            ref_lat
        )
        
        if shadow_rect:
            shadow_rectangles.append({
                "id": tree_id,
                "type": "tree",
                "height_m": height_m,
                "shadow_rectangle": shadow_rect,
                "sun_azimuth": round(sun_azimuth, 2),
                "sun_elevation": round(sun_elevation, 2)
            })
    
    return shadow_rectangles


def main():
    parser = argparse.ArgumentParser(
        description="Generate shadow rectangles for surrounding buildings and trees"
    )
    parser.add_argument(
        "plot_json",
        type=str,
        help="Path to plot.json file"
    )
    parser.add_argument(
        "--date",
        type=str,
        required=True,
        help="Date in YYYY-MM-DD format"
    )
    parser.add_argument(
        "--time",
        type=str,
        required=True,
        help="Time in HH:MM format (24-hour)"
    )
    parser.add_argument(
        "--timezone",
        type=str,
        default="Europe/London",
        help="Timezone (default: Europe/London)"
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output JSON file path (default: temp/shadow_rectangles_<timestamp>.json)"
    )
    parser.add_argument(
        "--temp-folder",
        type=str,
        default="temp",
        help="Temp folder for output (default: temp)"
    )
    
    args = parser.parse_args()
    
    # Parse date and time
    try:
        date_obj = datetime.strptime(args.date, "%Y-%m-%d")
        time_parts = args.time.split(":")
        if len(time_parts) != 2:
            raise ValueError("Time must be in HH:MM format")
        hour = int(time_parts[0])
        minute = int(time_parts[1])
        dt = date_obj.replace(hour=hour, minute=minute)
    except ValueError as e:
        print(f"Error parsing date/time: {e}")
        print("Date format: YYYY-MM-DD (e.g., 2024-06-21)")
        print("Time format: HH:MM (e.g., 14:30)")
        return 1
    
    # Load plot.json
    plot_path = Path(args.plot_json)
    if not plot_path.exists():
        print(f"Error: plot.json not found at {plot_path}")
        return 1
    
    with open(plot_path, 'r', encoding='utf-8') as f:
        plot_data = json.load(f)
    
    # Get plot centroid
    centroid = plot_data.get("centroid", {})
    if centroid.get("type") != "Point":
        print("Error: Invalid centroid in plot.json")
        return 1
    
    plot_centroid = tuple(centroid["coordinates"])  # (lon, lat)
    lat, lon = plot_centroid[1], plot_centroid[0]
    
    # Calculate sun position
    print(f"Calculating sun position for {dt.strftime('%Y-%m-%d %H:%M')} at ({lat:.6f}, {lon:.6f})...")
    sun_pos = calculate_sun_position(lat, lon, dt, args.timezone)
    print(f"  Sun azimuth: {sun_pos['azimuth']:.2f}° (0=North, 90=East, 180=South, 270=West)")
    print(f"  Sun elevation: {sun_pos['elevation']:.2f}°")
    
    if sun_pos['elevation'] <= 0:
        print("Warning: Sun is below horizon, no shadows will be generated")
    
    # Get surrounding context
    surrounding = plot_data.get("surrounding_context", {})
    buildings = surrounding.get("buildings", [])
    tree_zones = surrounding.get("tree_zones", {})
    
    print(f"\nProcessing {len(buildings)} buildings and trees...")
    
    # Process buildings
    building_shadows = process_buildings(
        buildings,
        sun_pos['azimuth'],
        sun_pos['elevation'],
        plot_centroid
    )
    print(f"  Generated {len(building_shadows)} building shadow rectangles")
    
    # Process trees
    tree_shadows = process_trees(
        tree_zones,
        sun_pos['azimuth'],
        sun_pos['elevation'],
        plot_centroid
    )
    print(f"  Generated {len(tree_shadows)} tree shadow rectangles")
    
    # Create output structure
    output_data = {
        "timestamp": dt.isoformat(),
        "location": {
            "latitude": lat,
            "longitude": lon,
            "timezone": args.timezone
        },
        "sun_position": {
            "azimuth": round(sun_pos['azimuth'], 2),
            "elevation": round(sun_pos['elevation'], 2)
        },
        "shadow_rectangles": {
            "buildings": building_shadows,
            "trees": tree_shadows
        },
        "summary": {
            "total_buildings": len(buildings),
            "total_trees": len(tree_zones.get("trees", [])),
            "shadows_generated": len(building_shadows) + len(tree_shadows)
        }
    }
    
    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        temp_folder = Path(args.temp_folder)
        temp_folder.mkdir(exist_ok=True)
        timestamp_str = dt.strftime("%Y%m%d_%H%M%S")
        output_path = temp_folder / f"shadow_rectangles_{timestamp_str}.json"
    
    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)
    
    print(f"\n[OK] Shadow rectangles written to: {output_path}")
    print(f"  Total shadows: {len(building_shadows) + len(tree_shadows)}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

