#!/usr/bin/env python
"""
Script to draw OSM data on Mapbox merged imagery

This script:
1. Loads a Mapbox merged image from cache
2. Loads OSM data (buildings, roads, trees, water) from plot JSON
3. Projects WGS84 coordinates to image pixel coordinates
4. Draws OSM features with different colors/styles
5. Saves the annotated image
"""

import json
import argparse
import math
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any

try:
    from PIL import Image, ImageDraw, ImageFont
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    print("ERROR: PIL (Pillow) is required. Install with: pip install Pillow")
    exit(1)

try:
    import mercantile
    HAS_MERCANTILE = True
except ImportError:
    HAS_MERCANTILE = False
    print("ERROR: mercantile is required. Install with: pip install mercantile")
    exit(1)

from loguru import logger

try:
    from pyproj import Transformer
    HAS_PYPROJ = True
except ImportError:
    HAS_PYPROJ = False
    logger.warning("pyproj not available - will use approximate WGS84 interpolation (less accurate)")


def wgs84_to_web_mercator(lon: float, lat: float) -> Tuple[float, float]:
    """
    Convert WGS84 (EPSG:4326) to Web Mercator (EPSG:3857)
    
    Args:
        lon: Longitude in degrees
        lat: Latitude in degrees
    
    Returns:
        (x, y) in Web Mercator meters
    """
    if HAS_PYPROJ:
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
        x, y = transformer.transform(lon, lat)
        return x, y
    else:
        # Fallback: approximate conversion (less accurate but no dependency)
        import math
        x = lon * 111320.0 * math.cos(math.radians(lat))
        y = lat * 110540.0
        return x, y


def lon_to_pixel(lon: float, bbox: Dict[str, float], image_width: int) -> int:
    """
    Convert longitude to pixel X coordinate in image using Web Mercator projection
    
    Args:
        lon: Longitude in WGS84 (degrees)
        bbox: Bounding box dict with 'west' and 'east' keys (degrees)
        image_width: Width of image in pixels
    
    Returns:
        Pixel X coordinate (0 to image_width-1)
    """
    if bbox["east"] == bbox["west"]:
        return image_width // 2
    
    # Convert to Web Mercator for accurate projection
    x_west, _ = wgs84_to_web_mercator(bbox["west"], (bbox["north"] + bbox["south"]) / 2)
    x_east, _ = wgs84_to_web_mercator(bbox["east"], (bbox["north"] + bbox["south"]) / 2)
    x_point, _ = wgs84_to_web_mercator(lon, (bbox["north"] + bbox["south"]) / 2)
    
    # Linear interpolation in Web Mercator space
    if x_east == x_west:
        return image_width // 2
    
    ratio = (x_point - x_west) / (x_east - x_west)
    pixel_x = int(ratio * image_width)
    return max(0, min(pixel_x, image_width - 1))


def lat_to_pixel(lat: float, bbox: Dict[str, float], image_height: int) -> int:
    """
    Convert latitude to pixel Y coordinate in image using Web Mercator projection
    Note: Image Y increases downward, but latitude increases upward.
    
    Args:
        lat: Latitude in WGS84 (degrees)
        bbox: Bounding box dict with 'south' and 'north' keys (degrees)
        image_height: Height of image in pixels
    
    Returns:
        Pixel Y coordinate (0 to image_height-1)
    """
    if bbox["north"] == bbox["south"]:
        return image_height // 2
    
    # Convert to Web Mercator for accurate projection
    center_lon = (bbox["west"] + bbox["east"]) / 2
    _, y_north = wgs84_to_web_mercator(center_lon, bbox["north"])
    _, y_south = wgs84_to_web_mercator(center_lon, bbox["south"])
    _, y_point = wgs84_to_web_mercator(center_lon, lat)
    
    # Linear interpolation in Web Mercator space (inverted because Y increases downward)
    if y_north == y_south:
        return image_height // 2
    
    ratio = (y_north - y_point) / (y_north - y_south)
    pixel_y = int(ratio * image_height)
    return max(0, min(pixel_y, image_height - 1))


def coords_to_pixels(
    coords: List[List[float]],
    bbox: Dict[str, float],
    image_size: Tuple[int, int]
) -> List[Tuple[int, int]]:
    """Convert WGS84 coordinates to pixel coordinates"""
    width, height = image_size
    pixels = []
    
    for coord in coords:
        if len(coord) < 2:
            continue
        lon, lat = coord[0], coord[1]
        x = lon_to_pixel(lon, bbox, width)
        y = lat_to_pixel(lat, bbox, height)
        pixels.append((x, y))
    
    return pixels


def find_mapbox_image(
    cache_dir: Path,
    lat: float,
    lon: float,
    radius_m: Optional[float] = None
) -> Optional[Path]:
    """Find Mapbox merged image in cache directory"""
    if not cache_dir.exists():
        return None
    
    if radius_m:
        pattern = f"mapbox_merged_{lat:.6f}_{lon:.6f}_{radius_m}m_z*.jpg"
    else:
        pattern = f"mapbox_merged_{lat:.6f}_{lon:.6f}_*m_z*.jpg"
    
    matches = list(cache_dir.glob(pattern))
    if matches:
        return max(matches, key=lambda p: p.stat().st_mtime)
    
    # Find closest match
    all_images = list(cache_dir.glob("mapbox_merged_*_*_*m_z*.jpg"))
    if not all_images:
        return None
    
    best_match = None
    min_distance = float('inf')
    
    for img_path in all_images:
        parts = img_path.stem.split('_')
        if len(parts) >= 5:
            try:
                img_lat = float(parts[2])
                img_lon = float(parts[3])
                dist = math.sqrt((img_lat - lat)**2 + (img_lon - lon)**2)
                if dist < min_distance:
                    min_distance = dist
                    best_match = img_path
            except ValueError:
                continue
    
    return best_match


def extract_bbox_from_filename(filename: str) -> Optional[Dict[str, float]]:
    """Extract bounding box from image filename"""
    parts = Path(filename).stem.split('_')
    if len(parts) < 5:
        return None
    
    try:
        lat = float(parts[2])
        lon = float(parts[3])
        radius_str = parts[4].replace('m', '')
        radius_m = float(radius_str)
    except (ValueError, IndexError):
        return None
    
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
    
    return bbox


def draw_buildings(
    draw: ImageDraw.Draw,
    buildings: List[Dict[str, Any]],
    bbox: Dict[str, float],
    image_size: Tuple[int, int],
    fill_color: Tuple[int, int, int, int],
    outline_color: Tuple[int, int, int, int],
    outline_width: int = 1
):
    """Draw building footprints as polygons"""
    drawn_count = 0
    
    for building in buildings:
        footprint = building.get("footprint", {})
        if not footprint:
            continue
        
        coords = footprint.get("coordinates", [[]])[0]
        if not coords or len(coords) < 3:
            continue
        
        # Convert to pixels
        pixel_coords = coords_to_pixels(coords, bbox, image_size)
        if len(pixel_coords) < 3:
            continue
        
        # Close polygon
        if pixel_coords[0] != pixel_coords[-1]:
            pixel_coords.append(pixel_coords[0])
        
        # Draw filled polygon
        draw.polygon(pixel_coords, fill=fill_color, outline=outline_color, width=outline_width)
        drawn_count += 1
    
    logger.info(f"Drew {drawn_count} buildings")


def draw_roads(
    draw: ImageDraw.Draw,
    roads: List[Dict[str, Any]],
    bbox: Dict[str, float],
    image_size: Tuple[int, int],
    line_color: Tuple[int, int, int, int],
    line_width: int = 2
):
    """Draw road centerlines as lines"""
    drawn_count = 0
    
    for road in roads:
        centerline = road.get("centerline", {})
        if not centerline:
            continue
        
        coords = centerline.get("coordinates", [])
        if not coords or len(coords) < 2:
            continue
        
        # Convert to pixels
        pixel_coords = coords_to_pixels(coords, bbox, image_size)
        if len(pixel_coords) < 2:
            continue
        
        # Draw line
        draw.line(pixel_coords, fill=line_color, width=line_width)
        drawn_count += 1
    
    logger.info(f"Drew {drawn_count} roads")


def draw_trees(
    draw: ImageDraw.Draw,
    trees: List[Dict[str, Any]],
    bbox: Dict[str, float],
    image_size: Tuple[int, int],
    fill_color: Tuple[int, int, int, int],
    outline_color: Tuple[int, int, int, int],
    radius: int = 3
):
    """Draw trees as circles"""
    drawn_count = 0
    
    for tree in trees:
        location = tree.get("location", [])
        if not location or len(location) < 2:
            continue
        
        lon, lat = location[0], location[1]
        width, height = image_size
        
        x = lon_to_pixel(lon, bbox, width)
        y = lat_to_pixel(lat, bbox, height)
        
        # Draw circle
        bbox_circle = [x - radius, y - radius, x + radius, y + radius]
        draw.ellipse(bbox_circle, fill=fill_color, outline=outline_color, width=1)
        drawn_count += 1
    
    logger.info(f"Drew {drawn_count} trees")


def draw_water_features(
    draw: ImageDraw.Draw,
    water_features: List[Dict[str, Any]],
    bbox: Dict[str, float],
    image_size: Tuple[int, int],
    fill_color: Tuple[int, int, int, int],
    outline_color: Tuple[int, int, int, int],
    outline_width: int = 1
):
    """Draw water features as polygons"""
    drawn_count = 0
    
    for feature in water_features:
        geometry = feature.get("geometry", {})
        if not geometry:
            continue
        
        geom_type = geometry.get("type", "")
        coords = geometry.get("coordinates", [])
        
        if geom_type == "Polygon" and coords:
            # Polygon coordinates: [[[lon, lat], ...]]
            polygon_coords = coords[0] if isinstance(coords[0][0], list) else coords
            if polygon_coords and len(polygon_coords) >= 3:
                pixel_coords = coords_to_pixels(polygon_coords, bbox, image_size)
                if len(pixel_coords) >= 3:
                    if pixel_coords[0] != pixel_coords[-1]:
                        pixel_coords.append(pixel_coords[0])
                    draw.polygon(pixel_coords, fill=fill_color, outline=outline_color, width=outline_width)
                    drawn_count += 1
        elif geom_type == "MultiPolygon" and coords:
            # MultiPolygon: [[[[lon, lat], ...], ...]]
            for polygon in coords:
                if polygon and len(polygon) > 0:
                    polygon_coords = polygon[0] if isinstance(polygon[0][0], list) else polygon
                    if polygon_coords and len(polygon_coords) >= 3:
                        pixel_coords = coords_to_pixels(polygon_coords, bbox, image_size)
                        if len(pixel_coords) >= 3:
                            if pixel_coords[0] != pixel_coords[-1]:
                                pixel_coords.append(pixel_coords[0])
                            draw.polygon(pixel_coords, fill=fill_color, outline=outline_color, width=outline_width)
                            drawn_count += 1
    
    logger.info(f"Drew {drawn_count} water features")


def load_osm_data(plot_json_path: Path) -> Dict[str, Any]:
    """Load OSM data from plot JSON"""
    try:
        with open(plot_json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        surrounding_context = data.get("surrounding_context", {})
        
        return {
            "buildings": surrounding_context.get("buildings", []),
            "roads": surrounding_context.get("roads", []),
            "trees": surrounding_context.get("tree_zones", {}).get("trees", []),
            "water_features": surrounding_context.get("water_features", {}).get("features", [])
        }
    except Exception as e:
        logger.error(f"Failed to load OSM data from {plot_json_path}: {e}")
        return {
            "buildings": [],
            "roads": [],
            "trees": [],
            "water_features": []
        }


def draw_osm_on_mapbox(
    image_path: Path,
    osm_data: Dict[str, Any],
    bbox: Dict[str, float],
    output_path: Path,
    image_alpha: float = 0.7,
    draw_buildings_flag: bool = True,
    draw_roads_flag: bool = True,
    draw_trees_flag: bool = True,
    draw_water_flag: bool = True,
    building_color: Tuple[int, int, int] = (128, 128, 128),  # Gray
    road_color: Tuple[int, int, int] = (255, 255, 255),  # White
    tree_color: Tuple[int, int, int] = (34, 139, 34),  # Forest green
    water_color: Tuple[int, int, int] = (0, 100, 200),  # Blue
    building_alpha: int = 150,
    road_width: int = 3,
    tree_radius: int = 4,
    water_alpha: int = 100
) -> bool:
    """
    Draw OSM data on Mapbox image with transparency
    
    Args:
        image_path: Path to Mapbox merged image
        osm_data: Dict with 'buildings', 'roads', 'trees', 'water_features' keys
        bbox: Bounding box dict
        output_path: Path to save annotated image
        image_alpha: Alpha transparency for base image (0.0-1.0)
        draw_buildings_flag: Whether to draw buildings
        draw_roads_flag: Whether to draw roads
        draw_trees_flag: Whether to draw trees
        draw_water_flag: Whether to draw water features
        building_color: RGB color for buildings
        road_color: RGB color for roads
        tree_color: RGB color for trees
        water_color: RGB color for water
        building_alpha: Alpha for building fill (0-255)
        road_width: Width of road lines in pixels
        tree_radius: Radius of tree circles in pixels
        water_alpha: Alpha for water fill (0-255)
    """
    try:
        # Load image
        img = Image.open(image_path)
        img_width, img_height = img.size
        
        logger.info(f"Loaded image: {img_width}x{img_height} pixels")
        logger.info(f"Bounding box: {bbox}")
        
        # Convert to RGBA to support transparency
        if img.mode != 'RGBA':
            img = img.convert('RGBA')
        
        # Apply alpha to base image
        if image_alpha < 1.0:
            alpha_channel = img.split()[3]
            alpha_channel = alpha_channel.point(lambda p: int(p * image_alpha))
            img.putalpha(alpha_channel)
        
        # Create overlay for OSM features
        overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)
        
        # Draw features in order: water (background), roads, buildings, trees (foreground)
        
        # 1. Water features (background layer)
        if draw_water_flag and osm_data.get("water_features"):
            draw_water_features(
                draw, osm_data["water_features"], bbox, (img_width, img_height),
                fill_color=water_color + (water_alpha,),
                outline_color=water_color + (255,),
                outline_width=2
            )
        
        # 2. Roads
        if draw_roads_flag and osm_data.get("roads"):
            draw_roads(
                draw, osm_data["roads"], bbox, (img_width, img_height),
                line_color=road_color + (255,),
                line_width=road_width
            )
        
        # 3. Buildings
        if draw_buildings_flag and osm_data.get("buildings"):
            draw_buildings(
                draw, osm_data["buildings"], bbox, (img_width, img_height),
                fill_color=building_color + (building_alpha,),
                outline_color=building_color + (255,),
                outline_width=1
            )
        
        # 4. Trees (foreground layer)
        if draw_trees_flag and osm_data.get("trees"):
            draw_trees(
                draw, osm_data["trees"], bbox, (img_width, img_height),
                fill_color=tree_color + (200,),
                outline_color=tree_color + (255,),
                radius=tree_radius
            )
        
        # Composite overlay on base image
        result = Image.alpha_composite(img, overlay)
        
        # Convert back to RGB for saving as JPEG
        if output_path.suffix.lower() in ['.jpg', '.jpeg']:
            background = Image.new('RGB', result.size, (255, 255, 255))
            background.paste(result, mask=result.split()[3] if result.mode == 'RGBA' else None)
            result = background
        
        # Save result
        result.save(output_path, quality=95)
        logger.info(f"Saved annotated image to: {output_path}")
        
        return True
        
    except Exception as e:
        logger.error(f"Failed to draw OSM data: {e}")
        import traceback
        traceback.print_exc()
        return False


def parse_color(color_str: str) -> Tuple[int, int, int]:
    """Parse color string like '255,0,0' to RGB tuple"""
    try:
        parts = [int(x.strip()) for x in color_str.split(',')]
        if len(parts) != 3:
            raise ValueError("Color must have 3 components")
        return tuple(parts)
    except Exception as e:
        raise ValueError(f"Invalid color format: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Draw OSM data on Mapbox merged imagery",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Draw all OSM features
  python scripts/draw_osm_on_mapbox.py \\
      --plot-json output/plot_ash1.json \\
      --cache-dir output/plot_ash1 \\
      --output output/plot_ash1_osm.jpg
  
  # Draw only buildings and roads
  python scripts/draw_osm_on_mapbox.py \\
      --plot-json output/plot_ash1.json \\
      --cache-dir output/plot_ash1 \\
      --output output/plot_ash1_osm.jpg \\
      --no-trees --no-water
  
  # Custom colors
  python scripts/draw_osm_on_mapbox.py \\
      --plot-json output/plot_ash1.json \\
      --cache-dir output/plot_ash1 \\
      --output output/plot_ash1_osm.jpg \\
      --building-color "100,100,100" \\
      --road-color "255,255,0" \\
      --tree-color "0,200,0"
        """
    )
    
    parser.add_argument(
        "--plot-json",
        type=str,
        required=True,
        help="Path to plot JSON file containing OSM data"
    )
    
    parser.add_argument(
        "--image",
        type=str,
        default=None,
        help="Path to Mapbox merged image (if not provided, will search in cache-dir)"
    )
    
    parser.add_argument(
        "--cache-dir",
        type=str,
        default=None,
        help="Cache directory to search for Mapbox images"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output path for annotated image"
    )
    
    # Feature flags
    parser.add_argument(
        "--no-buildings",
        action="store_true",
        help="Don't draw buildings"
    )
    
    parser.add_argument(
        "--no-roads",
        action="store_true",
        help="Don't draw roads"
    )
    
    parser.add_argument(
        "--no-trees",
        action="store_true",
        help="Don't draw trees"
    )
    
    parser.add_argument(
        "--no-water",
        action="store_true",
        help="Don't draw water features"
    )
    
    # Colors
    parser.add_argument(
        "--building-color",
        type=str,
        default="128,128,128",
        help="Building color as RGB (default: 128,128,128)"
    )
    
    parser.add_argument(
        "--road-color",
        type=str,
        default="255,255,255",
        help="Road color as RGB (default: 255,255,255)"
    )
    
    parser.add_argument(
        "--tree-color",
        type=str,
        default="34,139,34",
        help="Tree color as RGB (default: 34,139,34)"
    )
    
    parser.add_argument(
        "--water-color",
        type=str,
        default="0,100,200",
        help="Water color as RGB (default: 0,100,200)"
    )
    
    # Alpha and sizes
    parser.add_argument(
        "--image-alpha",
        type=float,
        default=0.7,
        help="Base image transparency 0.0-1.0 (default: 0.7)"
    )
    
    parser.add_argument(
        "--building-alpha",
        type=int,
        default=150,
        help="Building fill alpha 0-255 (default: 150)"
    )
    
    parser.add_argument(
        "--road-width",
        type=int,
        default=3,
        help="Road line width in pixels (default: 3)"
    )
    
    parser.add_argument(
        "--tree-radius",
        type=int,
        default=4,
        help="Tree circle radius in pixels (default: 4)"
    )
    
    parser.add_argument(
        "--water-alpha",
        type=int,
        default=100,
        help="Water fill alpha 0-255 (default: 100)"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger.remove()
    level = "DEBUG" if args.verbose else "INFO"
    logger.add(
        lambda msg: print(msg, end=''),
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level=level
    )
    
    # Load OSM data
    plot_json_path = Path(args.plot_json)
    if not plot_json_path.exists():
        logger.error(f"Plot JSON not found: {plot_json_path}")
        return 1
    
    osm_data = load_osm_data(plot_json_path)
    logger.info(f"Loaded OSM data: {len(osm_data['buildings'])} buildings, "
                f"{len(osm_data['roads'])} roads, {len(osm_data['trees'])} trees, "
                f"{len(osm_data['water_features'])} water features")
    
    # Find or load Mapbox image
    if args.image:
        image_path = Path(args.image)
        if not image_path.exists():
            logger.error(f"Image not found: {image_path}")
            return 1
    else:
        if not args.cache_dir:
            logger.error("Either --image or --cache-dir must be provided")
            return 1
        
        cache_dir = Path(args.cache_dir)
        
        # Try to get lat/lon from plot JSON
        try:
            with open(plot_json_path, 'r') as f:
                data = json.load(f)
            centroid = data.get("centroid", {}).get("coordinates", [])
            if len(centroid) >= 2:
                lon, lat = centroid[0], centroid[1]
                image_path = find_mapbox_image(cache_dir, lat, lon)
            else:
                image_path = None
        except Exception as e:
            logger.warning(f"Failed to get centroid from JSON: {e}")
            image_path = None
        
        if not image_path:
            all_images = list(cache_dir.glob("mapbox_merged_*.jpg"))
            if all_images:
                image_path = max(all_images, key=lambda p: p.stat().st_mtime)
                logger.info(f"Using most recent image: {image_path.name}")
            else:
                logger.error(f"No Mapbox images found in {cache_dir}")
                return 1
    
    logger.info(f"Using image: {image_path}")
    
    # Extract bounding box
    bbox = extract_bbox_from_filename(image_path.name)
    if not bbox:
        logger.error(f"Failed to extract bounding box from filename: {image_path.name}")
        return 1
    
    # Parse colors
    try:
        building_color = parse_color(args.building_color)
        road_color = parse_color(args.road_color)
        tree_color = parse_color(args.tree_color)
        water_color = parse_color(args.water_color)
    except ValueError as e:
        logger.error(f"Invalid color format: {e}")
        return 1
    
    # Draw OSM data on image
    output_path = Path(args.output)
    success = draw_osm_on_mapbox(
        image_path,
        osm_data,
        bbox,
        output_path,
        image_alpha=args.image_alpha,
        draw_buildings_flag=not args.no_buildings,
        draw_roads_flag=not args.no_roads,
        draw_trees_flag=not args.no_trees,
        draw_water_flag=not args.no_water,
        building_color=building_color,
        road_color=road_color,
        tree_color=tree_color,
        water_color=water_color,
        building_alpha=args.building_alpha,
        road_width=args.road_width,
        tree_radius=args.tree_radius,
        water_alpha=args.water_alpha
    )
    
    if success:
        logger.info(f"✓ Successfully created annotated image: {output_path}")
        return 0
    else:
        logger.error("Failed to create annotated image")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())

