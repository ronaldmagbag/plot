#!/usr/bin/env python
"""
Script to draw property line on Mapbox merged imagery

This script:
1. Loads a Mapbox merged image from cache
2. Loads property line coordinates from plot JSON
3. Projects WGS84 coordinates to image pixel coordinates
4. Draws property line with transparency
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
    
    # Clamp to image bounds
    return max(0, min(pixel_x, image_width - 1))


def lat_to_pixel(lat: float, bbox: Dict[str, float], image_height: int) -> int:
    """
    Convert latitude to pixel Y coordinate in image using Web Mercator projection
    
    Note: Image Y increases downward, but latitude increases upward.
    So we need to invert: north (higher lat) = smaller Y (top of image)
    
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
    
    # Clamp to image bounds
    return max(0, min(pixel_y, image_height - 1))


def coords_to_pixels(
    coords: List[List[float]],
    bbox: Dict[str, float],
    image_size: Tuple[int, int]
) -> List[Tuple[int, int]]:
    """
    Convert WGS84 coordinates to pixel coordinates
    
    Args:
        coords: List of [lon, lat] coordinate pairs in WGS84
        bbox: Bounding box dict with 'west', 'east', 'south', 'north' keys
        image_size: (width, height) tuple in pixels
    
    Returns:
        List of (x, y) pixel coordinate tuples
    """
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
    """
    Find Mapbox merged image in cache directory
    
    Args:
        cache_dir: Cache directory path
        lat: Center latitude
        lon: Center longitude
        radius_m: Optional radius to match (if None, finds any matching lat/lon)
    
    Returns:
        Path to image file or None if not found
    """
    if not cache_dir.exists():
        return None
    
    # Try to find image with matching lat/lon
    if radius_m:
        pattern = f"mapbox_merged_{lat:.6f}_{lon:.6f}_{radius_m}m_z*.jpg"
    else:
        pattern = f"mapbox_merged_{lat:.6f}_{lon:.6f}_*m_z*.jpg"
    
    matches = list(cache_dir.glob(pattern))
    if matches:
        # Return the most recent one
        return max(matches, key=lambda p: p.stat().st_mtime)
    
    # Try without exact lat/lon match (find closest)
    all_images = list(cache_dir.glob("mapbox_merged_*_*_*m_z*.jpg"))
    if not all_images:
        return None
    
    # Find closest match by lat/lon
    best_match = None
    min_distance = float('inf')
    
    for img_path in all_images:
        # Parse filename: mapbox_merged_{lat}_{lon}_{radius}m_z{zoom}.jpg
        parts = img_path.stem.split('_')
        if len(parts) >= 5:
            try:
                img_lat = float(parts[2])
                img_lon = float(parts[3])
                # Calculate distance
                dist = math.sqrt((img_lat - lat)**2 + (img_lon - lon)**2)
                if dist < min_distance:
                    min_distance = dist
                    best_match = img_path
            except ValueError:
                continue
    
    return best_match


def extract_bbox_from_filename(filename: str) -> Optional[Dict[str, float]]:
    """
    Extract bounding box from image filename or reconstruct from lat/lon/radius
    
    Args:
        filename: Image filename like "mapbox_merged_53.039271_-1.200618_20.0m_z19.jpg"
    
    Returns:
        Bounding box dict or None
    """
    # Parse filename
    parts = Path(filename).stem.split('_')
    if len(parts) < 5:
        return None
    
    try:
        lat = float(parts[2])
        lon = float(parts[3])
        radius_str = parts[4].replace('m', '')
        radius_m = float(radius_str)
        zoom = int(parts[6].replace('z', '')) if len(parts) > 6 else 19
    except (ValueError, IndexError):
        return None
    
    # Reconstruct bounding box from lat/lon/radius
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
    
    return bbox


def load_property_line(plot_json_path: Path) -> Optional[List[List[float]]]:
    """
    Load property line coordinates from plot JSON
    
    Args:
        plot_json_path: Path to plot JSON file
    
    Returns:
        List of [lon, lat] coordinate pairs or None
    """
    try:
        with open(plot_json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        prop_line = data.get("boundaries", {}).get("property_line", {})
        # Use simplified coordinates if available, otherwise use original
        if prop_line.get("coordinates_simplified"):
            coords = prop_line.get("coordinates_simplified", [[]])[0]
        else:
        coords = prop_line.get("coordinates", [[]])[0]
        
        if not coords:
            logger.warning("No property line coordinates found in JSON")
            return None
        
        # Ensure coordinates are [lon, lat] format
        # Check if first coordinate looks like [lat, lon] (lat typically 40-70 for UK)
        if len(coords[0]) >= 2:
            first_lat = coords[0][1] if isinstance(coords[0][1], (int, float)) else coords[0][0]
            if 40 < first_lat < 70:  # UK latitude range
                # Already in [lon, lat] format
                return coords
            else:
                # Might be [lat, lon], swap
                return [[c[1], c[0]] for c in coords if len(c) >= 2]
        
        return coords
        
    except Exception as e:
        logger.error(f"Failed to load property line from {plot_json_path}: {e}")
        return None


def draw_property_on_mapbox(
    image_path: Path,
    property_coords: List[List[float]],
    bbox: Dict[str, float],
    output_path: Path,
    line_color: Tuple[int, int, int] = (255, 0, 0),  # Red
    line_width: int = 3,
    image_alpha: float = 0.7
) -> bool:
    """
    Draw property line on Mapbox image with transparency
    
    Args:
        image_path: Path to Mapbox merged image
        property_coords: List of [lon, lat] coordinate pairs
        bbox: Bounding box dict
        output_path: Path to save annotated image
        line_color: RGB color for property line
        line_width: Width of property line in pixels
        image_alpha: Alpha transparency for base image (0.0-1.0)
    
    Returns:
        True if successful, False otherwise
    """
    try:
        # Load image
        img = Image.open(image_path)
        img_width, img_height = img.size
        
        logger.info(f"Loaded image: {img_width}x{img_height} pixels")
        logger.info(f"Bounding box: {bbox}")
        
        # Convert property coordinates to pixels
        pixel_coords = coords_to_pixels(property_coords, bbox, (img_width, img_height))
        
        if len(pixel_coords) < 3:
            logger.error(f"Not enough pixel coordinates: {len(pixel_coords)}")
            return False
        
        logger.info(f"Converted {len(property_coords)} coordinates to {len(pixel_coords)} pixels")
        
        # Create a copy for drawing
        # Convert to RGBA to support transparency
        if img.mode != 'RGBA':
            img = img.convert('RGBA')
        
        # Create overlay for property line (fully opaque)
        overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)
        
        # Draw property line
        draw.line(pixel_coords, fill=line_color + (255,), width=line_width)
        
        # Also draw filled polygon with transparency for better visibility
        if len(pixel_coords) >= 3:
            # Close the polygon
            closed_coords = pixel_coords + [pixel_coords[0]]
            # Draw filled polygon with low opacity
            draw.polygon(closed_coords, fill=line_color + (50,), outline=line_color + (255,), width=line_width)
        
        # Apply alpha to base image
        if image_alpha < 1.0:
            alpha_channel = img.split()[3]  # Get alpha channel
            alpha_channel = alpha_channel.point(lambda p: int(p * image_alpha))
            img.putalpha(alpha_channel)
        
        # Composite overlay on base image
        result = Image.alpha_composite(img, overlay)
        
        # Convert back to RGB for saving as JPEG
        if output_path.suffix.lower() in ['.jpg', '.jpeg']:
            # Create white background
            background = Image.new('RGB', result.size, (255, 255, 255))
            background.paste(result, mask=result.split()[3] if result.mode == 'RGBA' else None)
            result = background
        
        # Save result
        result.save(output_path, quality=95)
        logger.info(f"Saved annotated image to: {output_path}")
        
        return True
        
    except Exception as e:
        logger.error(f"Failed to draw property line: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Draw property line on Mapbox merged imagery",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Draw property line on mapbox image using plot JSON
  python scripts/draw_property_on_mapbox.py \\
      --plot-json output/plot_ash1.json \\
      --cache-dir output/plot_ash1 \\
      --output output/plot_ash1_with_property.jpg
  
  # Specify exact image file
  python scripts/draw_property_on_mapbox.py \\
      --image output/plot_ash1/mapbox_merged_53.039271_-1.200618_20.0m_z19.jpg \\
      --plot-json output/plot_ash1.json \\
      --output output/annotated.jpg
        """
    )
    
    parser.add_argument(
        "--plot-json",
        type=str,
        required=True,
        help="Path to plot JSON file containing property line coordinates"
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
        help="Cache directory to search for Mapbox images (if --image not provided)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output path for annotated image"
    )
    
    parser.add_argument(
        "--line-color",
        type=str,
        default="255,0,0",
        help="Property line color as RGB (default: 255,0,0 for red)"
    )
    
    parser.add_argument(
        "--line-width",
        type=int,
        default=3,
        help="Property line width in pixels (default: 3)"
    )
    
    parser.add_argument(
        "--image-alpha",
        type=float,
        default=0.7,
        help="Base image transparency 0.0-1.0 (default: 0.7)"
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
    
    # Load property line from JSON
    plot_json_path = Path(args.plot_json)
    if not plot_json_path.exists():
        logger.error(f"Plot JSON not found: {plot_json_path}")
        return 1
    
    property_coords = load_property_line(plot_json_path)
    if not property_coords:
        logger.error("Failed to load property line coordinates")
        return 1
    
    logger.info(f"Loaded {len(property_coords)} property line coordinates")
    
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
            # Find any mapbox image in cache
            all_images = list(cache_dir.glob("mapbox_merged_*.jpg"))
            if all_images:
                image_path = max(all_images, key=lambda p: p.stat().st_mtime)
                logger.info(f"Using most recent image: {image_path.name}")
            else:
                logger.error(f"No Mapbox images found in {cache_dir}")
                return 1
    
    logger.info(f"Using image: {image_path}")
    
    # Extract bounding box from image filename
    bbox = extract_bbox_from_filename(image_path.name)
    if not bbox:
        logger.error(f"Failed to extract bounding box from filename: {image_path.name}")
        return 1
    
    logger.info(f"Image bounding box: {bbox}")
    
    # Parse line color
    try:
        color_parts = [int(x.strip()) for x in args.line_color.split(',')]
        if len(color_parts) != 3:
            raise ValueError("Color must have 3 components")
        line_color = tuple(color_parts)
    except Exception as e:
        logger.error(f"Invalid line color format: {e}")
        return 1
    
    # Draw property line on image
    output_path = Path(args.output)
    success = draw_property_on_mapbox(
        image_path,
        property_coords,
        bbox,
        output_path,
        line_color=line_color,
        line_width=args.line_width,
        image_alpha=args.image_alpha
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

