#!/usr/bin/env python
"""
Command-line interface for Plot Analysis Generator

Usage:
    python cli.py generate --lat 51.5074 --lon -0.1276 --output plot.json
    python cli.py batch --input locations.csv --output ./plots/
"""

import os
import sys
import json
import csv
import argparse
from datetime import datetime

# Add src to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from loguru import logger
from src.pipeline import PlotAnalysisPipeline


def setup_logging(verbose: bool = False):
    """Configure logging"""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stderr,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level=level
    )


def cmd_generate(args):
    """Generate plot analysis for a single location"""
    setup_logging(args.verbose)
    
    logger.info(f"Generating plot analysis for ({args.lat}, {args.lon})")
    
    # Determine output path and cache directory
    output_path = args.output or f"plot_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    
    # Create cache directory: output/plot_name/ (without .json extension)
    output_dir = os.path.dirname(output_path) or "."
    output_basename = os.path.basename(output_path)
    cache_dir = os.path.join(output_dir, os.path.splitext(output_basename)[0])
    os.makedirs(cache_dir, exist_ok=True)
    logger.info(f"Cache directory: {cache_dir}")
    
    pipeline = PlotAnalysisPipeline(cache_dir=cache_dir)
    
    try:
        result = pipeline.run(
            lat=args.lat,
            lon=args.lon,
            radius_m=args.radius,
            plot_id=args.plot_id
        )
        
        # Save to file
        pipeline.save(result, output_path)
        
        logger.info(f"✓ Generated: {output_path}")
        logger.info(f"  Plot ID: {result.plot_id}")
        logger.info(f"  Property Area: {result.boundaries.property_line.area_sqm} sqm")
        logger.info(f"  Buildable Area: {result.boundaries.buildable_envelope.area_sqm} sqm")
        
        # Print summary to stdout if requested
        if args.summary:
            summary = {
                "plot_id": result.plot_id,
                "centroid": result.centroid.coordinates,
                "property_area_sqm": result.boundaries.property_line.area_sqm,
                "buildable_area_sqm": result.boundaries.buildable_envelope.area_sqm,
                "neighbors": len(result.surrounding_context.buildings),
                "roads": len(result.surrounding_context.roads),
                "elevation_m": result.surrounding_context.elevation_map.average_elevation_m,
                "soil_type": result.soil.type_id,
                "best_solar_facade": result.analysis.shadow_analysis.best_solar_facade
            }
            print(json.dumps(summary, indent=2))
        
        return 0
        
    except Exception as e:
        logger.error(f"Failed to generate plot: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def cmd_batch(args):
    """Generate plot analyses for multiple locations from CSV"""
    setup_logging(args.verbose)
    
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        return 1
    
    # Read locations from CSV
    locations = []
    with open(args.input, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                locations.append({
                    "name": row.get("name", ""),
                    "lat": float(row["lat"]),
                    "lon": float(row["lon"]),
                    "plot_id": row.get("plot_id")
                })
            except (KeyError, ValueError) as e:
                logger.warning(f"Skipping invalid row: {e}")
    
    if not locations:
        logger.error("No valid locations found in CSV")
        return 1
    
    logger.info(f"Processing {len(locations)} locations...")
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    pipeline = PlotAnalysisPipeline()
    success = 0
    failed = 0
    
    for i, loc in enumerate(locations, 1):
        name = loc.get("name") or f"plot_{i:03d}"
        logger.info(f"[{i}/{len(locations)}] {name}: ({loc['lat']}, {loc['lon']})")
        
        try:
            result = pipeline.run(
                lat=loc["lat"],
                lon=loc["lon"],
                plot_id=loc.get("plot_id")
            )
            
            filename = f"{name.replace(' ', '_').lower()}.json"
            filepath = os.path.join(args.output, filename)
            pipeline.save(result, filepath)
            
            logger.info(f"  ✓ {filename}")
            success += 1
            
        except Exception as e:
            logger.error(f"  ✗ Failed: {e}")
            failed += 1
        
        # Rate limiting
        if i < len(locations):
            import time
            time.sleep(args.delay)
    
    logger.info(f"\nComplete: {success} succeeded, {failed} failed")
    return 0 if failed == 0 else 1


def cmd_visualize(args):
    """Visualize a plot.json file with full details"""
    setup_logging(args.verbose)
    
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        return 1
    
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        from matplotlib.patches import Polygon as MplPolygon, FancyBboxPatch, Rectangle
        from matplotlib.lines import Line2D
        import numpy as np
    except ImportError:
        logger.error("matplotlib is required for visualization. Install with: pip install matplotlib")
        return 1
    
    # Load plot data
    with open(args.input, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    logger.info(f"Visualizing: {data.get('plot_id', 'Unknown')}")
    
    # Create figure with better size
    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    
    # Get reference point
    centroid = data.get("centroid", {}).get("coordinates", [0, 0])
    ref_lon, ref_lat = centroid
    
    # Detect if coordinates are in degrees (lat/lon) or local meters
    # Lat/lon: latitude typically 49-61 for UK (e.g. 51.45), lon is -10 to +2
    # Local meters: coordinates would be small relative values like 0-100, not in lat/lon range
    is_degrees = True  # Default to degrees for real OSM data
    if data.get("boundaries", {}).get("property_line", {}).get("coordinates"):
        first_coord = data["boundaries"]["property_line"]["coordinates"][0][0]
        # Lat is typically 49-61 for UK. If the second value (lat) is in this range, it's degrees
        lat_val = first_coord[1]
        # If lat is in typical UK range (49-61), it's definitely degrees
        # If lat is small (e.g. 0-30), it's likely local meters (like example_plot.json)
        if lat_val > 40 and lat_val < 70:
            is_degrees = True
        elif abs(lat_val) < 100:
            is_degrees = False
    
    def to_local(coords):
        """Convert coordinates to local meters relative to centroid"""
        if is_degrees:
            # Convert from degrees to meters
            m_per_deg_lat = 111000
            m_per_deg_lon = 111000 * np.cos(np.radians(ref_lat))
            return [((c[0] - ref_lon) * m_per_deg_lon, (c[1] - ref_lat) * m_per_deg_lat) for c in coords]
        else:
            # Already in local meters, just offset from centroid
            return [(c[0] - ref_lon, c[1] - ref_lat) for c in coords]
    
    def get_polygon_center(coords):
        """Get center of polygon"""
        if not coords:
            return (0, 0)
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        return (sum(xs) / len(xs), sum(ys) / len(ys))
    
    def clip_to_bbox(coords, bbox):
        """Clip coordinates to bounding box [-bbox, +bbox] in both x and y"""
        clipped = []
        for x, y in coords:
            if -bbox <= x <= bbox and -bbox <= y <= bbox:
                clipped.append((x, y))
        return clipped
    
    def line_intersects_bbox(p1, p2, bbox):
        """Check if line segment intersects with bounding box"""
        x1, y1 = p1
        x2, y2 = p2
        
        # Check if either endpoint is inside
        if (-bbox <= x1 <= bbox and -bbox <= y1 <= bbox) or \
           (-bbox <= x2 <= bbox and -bbox <= y2 <= bbox):
            return True
        
        # Check if line segment crosses any edge of the bbox
        # This is a simplified check - for full accuracy, use line-box intersection
        min_x, max_x = min(x1, x2), max(x1, x2)
        min_y, max_y = min(y1, y2), max(y1, y2)
        
        # Check if line segment overlaps with bbox
        if max_x < -bbox or min_x > bbox or max_y < -bbox or min_y > bbox:
            return False
        
        return True
    
    def clip_line_to_bbox(coords, bbox):
        """Clip line to bounding box using Cohen-Sutherland algorithm"""
        if not coords:
            return []
        
        def compute_code(x, y, bbox):
            """Compute region code for point"""
            code = 0
            if x < -bbox:
                code |= 1  # Left
            elif x > bbox:
                code |= 2  # Right
            if y < -bbox:
                code |= 4  # Bottom
            elif y > bbox:
                code |= 8  # Top
            return code
        
        def clip_segment(p1, p2, bbox):
            """Clip line segment to bbox using Cohen-Sutherland"""
            x1, y1 = p1
            x2, y2 = p2
            
            code1 = compute_code(x1, y1, bbox)
            code2 = compute_code(x2, y2, bbox)
            
            while True:
                # Both endpoints inside
                if code1 == 0 and code2 == 0:
                    return [(x1, y1), (x2, y2)]
                
                # Both endpoints outside, same side
                if (code1 & code2) != 0:
                    return None
                
                # Pick endpoint outside
                code_out = code1 if code1 != 0 else code2
                
                # Find intersection point
                # Avoid division by zero
                if code_out & 8:  # Top
                    if abs(y2 - y1) < 1e-10:
                        return None  # Horizontal line, skip
                    x = x1 + (x2 - x1) * (bbox - y1) / (y2 - y1)
                    y = bbox
                elif code_out & 4:  # Bottom
                    if abs(y2 - y1) < 1e-10:
                        return None  # Horizontal line, skip
                    x = x1 + (x2 - x1) * (-bbox - y1) / (y2 - y1)
                    y = -bbox
                elif code_out & 2:  # Right
                    if abs(x2 - x1) < 1e-10:
                        return None  # Vertical line, skip
                    y = y1 + (y2 - y1) * (bbox - x1) / (x2 - x1)
                    x = bbox
                elif code_out & 1:  # Left
                    if abs(x2 - x1) < 1e-10:
                        return None  # Vertical line, skip
                    y = y1 + (y2 - y1) * (-bbox - x1) / (x2 - x1)
                    x = -bbox
                
                # Replace point outside with intersection point
                if code_out == code1:
                    x1, y1 = x, y
                    code1 = compute_code(x1, y1, bbox)
                else:
                    x2, y2 = x, y
                    code2 = compute_code(x2, y2, bbox)
        
        clipped = []
        for i in range(len(coords) - 1):
            p1, p2 = coords[i], coords[i + 1]
            segment = clip_segment(p1, p2, bbox)
            if segment:
                # Add points, ensuring continuity
                if not clipped:
                    # First segment - add both points
                    clipped.append(segment[0])
                    clipped.append(segment[1])
                else:
                    # Check distance to last point
                    last_point = clipped[-1]
                    dist_to_first = ((segment[0][0] - last_point[0])**2 + (segment[0][1] - last_point[1])**2)**0.5
                    
                    if dist_to_first > 1.0:  # More than 1m away, add the first point to maintain continuity
                        clipped.append(segment[0])
                    # Always add the second point
                    clipped.append(segment[1])
        
        # Remove duplicate consecutive points
        if clipped:
            cleaned = [clipped[0]]
            for pt in clipped[1:]:
                last_pt = cleaned[-1]
                dist = ((pt[0] - last_pt[0])**2 + (pt[1] - last_pt[1])**2)**0.5
                if dist > 0.1:  # Only add if more than 0.1m away
                    cleaned.append(pt)
            return cleaned
        
        return clipped
    
    def polygon_intersects_bbox(coords, bbox):
        """Check if polygon intersects or is within bounding box"""
        if not coords:
            return False
        
        # Check if any vertex is inside bbox
        for x, y in coords:
            if -bbox <= x <= bbox and -bbox <= y <= bbox:
                return True
        
        # Check if polygon might cross bbox (simplified - check bounding box of polygon)
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        poly_min_x, poly_max_x = min(xs), max(xs)
        poly_min_y, poly_max_y = min(ys), max(ys)
        
        # Check if polygon bounding box overlaps with our bbox
        if poly_max_x < -bbox or poly_min_x > bbox or poly_max_y < -bbox or poly_min_y > bbox:
            return False
        
        return True
    
    # Define 60m bounding box (centered on centroid at 0,0)
    VIEW_BBOX = 60.0  # 60 meters in each direction from centroid
    
    # Get property data
    prop_area = data.get("boundaries", {}).get("property_line", {}).get("area_sqm", 0)
    setback_area = data.get("boundaries", {}).get("setback_line", {}).get("area_sqm", 0)
    buildable_area = data.get("boundaries", {}).get("buildable_envelope", {}).get("area_sqm", 0)
    
    # ============================================================
    # LAYER 1: Water features (background) - CLIPPED to 60m box
    # ============================================================
    water_data = data.get("surrounding_context", {}).get("water_features", {})
    water_features = water_data.get("features", [])
    for wf in water_features:
        geom = wf.get("geometry", {})
        if geom.get("type") == "Polygon":
            w_coords = geom.get("coordinates", [[]])[0]
            if w_coords:
                local_w = to_local(w_coords)
                # Only show water features that intersect the 60m bounding box
                if polygon_intersects_bbox(local_w, VIEW_BBOX):
                    poly = MplPolygon(local_w, fill=True, facecolor='lightblue', alpha=0.5, 
                                     edgecolor='steelblue', linewidth=1)
                    ax.add_patch(poly)
    
    # ============================================================
    # LAYER 2: Tree zones as CIRCLES (canopy areas) - CLIPPED to 60m box
    # ============================================================
    tree_zones_data = data.get("surrounding_context", {}).get("tree_zones", {})
    trees = tree_zones_data.get("trees", [])
    for tree in trees:
        loc = tree.get("location", [])
        if loc and len(loc) == 2:
            local_t = to_local([loc])[0]
            # Only show trees within 60m bounding box
            if not (-VIEW_BBOX <= local_t[0] <= VIEW_BBOX and -VIEW_BBOX <= local_t[1] <= VIEW_BBOX):
                continue
            
            # Canopy radius based on tree height (approx radius = height * 0.4)
            height_m = tree.get("height_m") or 8  # Default to 8m if None or missing
            canopy_radius = height_m * 0.4
            
            # Draw tree canopy as filled circle
            circle = plt.Circle(local_t, canopy_radius, 
                               facecolor='forestgreen', alpha=0.4,
                               edgecolor='darkgreen', linewidth=1.5, zorder=3)
            ax.add_patch(circle)
            
            # Draw trunk marker at center
            ax.plot(local_t[0], local_t[1], 'o', color='saddlebrown', 
                   markersize=4, zorder=4)
    
    # ============================================================
    # LAYER 3: Roads as POLYGONS (not lines) - CLIPPED to 60m box
    # ============================================================
    roads = data.get("surrounding_context", {}).get("roads", [])
    for road in roads:
        r_coords = road.get("centerline", {}).get("coordinates", [])
        if r_coords and len(r_coords) >= 2:
            local_r = to_local(r_coords)
            
            # Check if road intersects the 60m box
            road_intersects = False
            for i in range(len(local_r) - 1):
                if line_intersects_bbox(local_r[i], local_r[i+1], VIEW_BBOX):
                    road_intersects = True
                    break
            
            if not road_intersects:
                continue  # Skip roads that don't intersect the view box
            
            # Simple clipping: process each segment and add intersection points
            local_r_clipped = []
            
            for i in range(len(local_r) - 1):
                p1 = local_r[i]
                p2 = local_r[i + 1]
                x1, y1 = p1
                x2, y2 = p2
                
                p1_inside = (-VIEW_BBOX <= x1 <= VIEW_BBOX and -VIEW_BBOX <= y1 <= VIEW_BBOX)
                p2_inside = (-VIEW_BBOX <= x2 <= VIEW_BBOX and -VIEW_BBOX <= y2 <= VIEW_BBOX)
                
                if p1_inside and p2_inside:
                    # Both inside - add both points
                    if not local_r_clipped or local_r_clipped[-1] != p1:
                        local_r_clipped.append(p1)
                    local_r_clipped.append(p2)
                elif p1_inside or p2_inside or line_intersects_bbox(p1, p2, VIEW_BBOX):
                    # Segment crosses boundary or has one point inside
                    # Find intersection points
                    intersections = []
                    
                    # Check each boundary
                    dx = x2 - x1
                    dy = y2 - y1
                    
                    if abs(dx) > 1e-10:
                        # Check left boundary (x = -VIEW_BBOX)
                        t = (-VIEW_BBOX - x1) / dx
                        if 0 <= t <= 1:
                            y_int = y1 + t * dy
                            if -VIEW_BBOX <= y_int <= VIEW_BBOX:
                                intersections.append(((-VIEW_BBOX, y_int), t))
                        
                        # Check right boundary (x = VIEW_BBOX)
                        t = (VIEW_BBOX - x1) / dx
                        if 0 <= t <= 1:
                            y_int = y1 + t * dy
                            if -VIEW_BBOX <= y_int <= VIEW_BBOX:
                                intersections.append(((VIEW_BBOX, y_int), t))
                    
                    if abs(dy) > 1e-10:
                        # Check bottom boundary (y = -VIEW_BBOX)
                        t = (-VIEW_BBOX - y1) / dy
                        if 0 <= t <= 1:
                            x_int = x1 + t * dx
                            if -VIEW_BBOX <= x_int <= VIEW_BBOX:
                                intersections.append(((x_int, -VIEW_BBOX), t))
                        
                        # Check top boundary (y = VIEW_BBOX)
                        t = (VIEW_BBOX - y1) / dy
                        if 0 <= t <= 1:
                            x_int = x1 + t * dx
                            if -VIEW_BBOX <= x_int <= VIEW_BBOX:
                                intersections.append(((x_int, VIEW_BBOX), t))
                    
                    # Sort intersections by t (parameter along segment)
                    intersections.sort(key=lambda x: x[1])
                    
                    # Add points in order
                    if p1_inside and (not local_r_clipped or local_r_clipped[-1] != p1):
                        local_r_clipped.append(p1)
                    
                    for int_pt, _ in intersections:
                        if not local_r_clipped or abs(local_r_clipped[-1][0] - int_pt[0]) > 0.1 or abs(local_r_clipped[-1][1] - int_pt[1]) > 0.1:
                            local_r_clipped.append(int_pt)
                    
                    if p2_inside:
                        local_r_clipped.append(p2)
            
            # Remove duplicate consecutive points
            if local_r_clipped:
                cleaned = [local_r_clipped[0]]
                for pt in local_r_clipped[1:]:
                    last_pt = cleaned[-1]
                    dist = ((pt[0] - last_pt[0])**2 + (pt[1] - last_pt[1])**2)**0.5
                    if dist > 0.1:
                        cleaned.append(pt)
                local_r_clipped = cleaned
            
            if len(local_r_clipped) < 2:
                continue  # Can't create valid road
            
            width = road.get("width_m", 6.0)
            half_width = width / 2
            
            # Create road polygon from clipped centerline
            # For curved roads, use average direction at each point for smoother offset
            road_polygon_left = []
            road_polygon_right = []
            
            for i in range(len(local_r_clipped)):
                x, y = local_r_clipped[i]
                
                # Calculate average direction at this point (for smooth curves)
                if i == 0:
                    # First point: use direction to next
                    dx = local_r_clipped[i+1][0] - x
                    dy = local_r_clipped[i+1][1] - y
                elif i == len(local_r_clipped) - 1:
                    # Last point: use direction from previous
                    dx = x - local_r_clipped[i-1][0]
                    dy = y - local_r_clipped[i-1][1]
                else:
                    # Middle point: average of incoming and outgoing directions
                    dx1 = local_r_clipped[i+1][0] - x
                    dy1 = local_r_clipped[i+1][1] - y
                    dx2 = x - local_r_clipped[i-1][0]
                    dy2 = y - local_r_clipped[i-1][1]
                    # Average the directions
                    dx = (dx1 + dx2) / 2
                    dy = (dy1 + dy2) / 2
                
                length = np.sqrt(dx*dx + dy*dy)
                if length > 0:
                    # Normalize
                    dx /= length
                    dy /= length
                    # Perpendicular direction (rotate 90 degrees counterclockwise)
                    px, py = -dy * half_width, dx * half_width
                    road_polygon_left.append((x + px, y + py))
                    road_polygon_right.append((x - px, y - py))
            
            # Combine left and right sides (reverse right side and close polygon)
            if len(road_polygon_left) >= 2 and len(road_polygon_right) >= 2:
                road_polygon = road_polygon_left + list(reversed(road_polygon_right))
                # Close the polygon if not already closed
                if road_polygon[0] != road_polygon[-1]:
                    road_polygon.append(road_polygon[0])
                
                # Draw road polygon
                try:
                    poly = MplPolygon(road_polygon, fill=True, facecolor='lightgray', 
                                     alpha=0.8, edgecolor='dimgray', linewidth=1, zorder=2)
                    ax.add_patch(poly)
                except Exception as e:
                    # If polygon creation fails, draw as thick line instead
                    logger.warning(f"Failed to create road polygon for {road.get('name', 'unknown')}: {e}")
                    xs, ys = zip(*local_r_clipped)
                    ax.plot(xs, ys, color='lightgray', linewidth=width*2, alpha=0.8, zorder=2, solid_capstyle='round')
                
                # Center line dashed
                xs, ys = zip(*local_r_clipped)
                ax.plot(xs, ys, color='white', linewidth=1.5, linestyle='--', alpha=0.9, zorder=3)
            
            # Road name and width label (positioned above road)
            mid_idx = len(local_r_clipped) // 2
            if mid_idx < len(local_r_clipped):
                mid_x, mid_y = local_r_clipped[mid_idx]
            else:
                mid_x, mid_y = local_r_clipped[0] if local_r_clipped else (0, 0)
            road_name = road.get("name", "Road")
            ax.annotate(f"{road_name}\n(width: {width}m)", 
                       (mid_x, mid_y + half_width + 2), 
                       fontsize=8, color='dimgray', ha='center', va='bottom',
                       fontweight='bold', zorder=10)
    
    # Store property polygon for distance calculations
    prop_coords_raw = data.get("boundaries", {}).get("property_line", {}).get("coordinates", [[]])[0]
    local_prop_for_dist = to_local(prop_coords_raw) if prop_coords_raw else []
    
    # ============================================================
    # LAYER 4: Neighbor buildings with height and PROPER distance labels - CLIPPED to 60m box
    # ============================================================
    buildings = data.get("surrounding_context", {}).get("buildings", [])
    for building in buildings:
        b_coords = building.get("footprint", {}).get("coordinates", [[]])[0]
        if b_coords:
            local_b = to_local(b_coords)
            
            # Only show buildings that intersect the 60m bounding box
            if not polygon_intersects_bbox(local_b, VIEW_BBOX):
                continue
            # Draw building with hatching
            poly = MplPolygon(local_b, fill=True, facecolor='slategray', alpha=0.6, 
                             edgecolor='darkslategray', linewidth=1.5, hatch='///', zorder=4)
            ax.add_patch(poly)
            
            # Get building info
            height = building.get("height_m", 8.0)
            stories = building.get("stories", 2)
            b_type = building.get("building_type", "house").replace("_", " ").title()
            dist = building.get("distance_to_property_line_m")
            wall_facing = building.get("wall_facing_plot", "")
            
            # Label with height and type at building center
            b_center = get_polygon_center(local_b)
            label_text = f"{height}m\n{b_type}"
            ax.annotate(label_text, b_center, fontsize=7, ha='center', va='center',
                       color='white', fontweight='bold', zorder=10,
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='slategray', alpha=0.8))
            
            # Draw distance indicator between building and property line
            if dist is not None and dist > 0 and wall_facing and local_prop_for_dist:
                # Find the edge of the building closest to property
                # and the corresponding property line edge
                b_xs = [p[0] for p in local_b]
                b_ys = [p[1] for p in local_b]
                
                # Determine which edge based on wall_facing
                if wall_facing == "west":
                    # Building's west edge -> property's east edge
                    b_edge_x = min(b_xs)
                    b_edge_y = b_center[1]
                    p_edge_x = max([p[0] for p in local_prop_for_dist])
                    p_edge_y = b_center[1]
                elif wall_facing == "east":
                    # Building's east edge -> property's west edge
                    b_edge_x = max(b_xs)
                    b_edge_y = b_center[1]
                    p_edge_x = min([p[0] for p in local_prop_for_dist])
                    p_edge_y = b_center[1]
                elif wall_facing == "north":
                    # Building's north edge -> property's south edge
                    b_edge_x = b_center[0]
                    b_edge_y = max(b_ys)
                    p_edge_x = b_center[0]
                    p_edge_y = min([p[1] for p in local_prop_for_dist])
                elif wall_facing == "south":
                    # Building's south edge -> property's north edge
                    b_edge_x = b_center[0]
                    b_edge_y = min(b_ys)
                    p_edge_x = b_center[0]
                    p_edge_y = max([p[1] for p in local_prop_for_dist])
                else:
                    continue
                
                # Draw distance line with arrows
                ax.annotate('', xy=(p_edge_x, p_edge_y), xytext=(b_edge_x, b_edge_y),
                           arrowprops=dict(arrowstyle='<->', color='red', lw=1.5),
                           zorder=8)
                
                # Label at midpoint of the distance line
                mid_x = (b_edge_x + p_edge_x) / 2
                mid_y = (b_edge_y + p_edge_y) / 2
                ax.annotate(f"{dist:.1f}m", (mid_x, mid_y),
                           fontsize=8, ha='center', va='center', color='red',
                           fontweight='bold', zorder=10,
                           bbox=dict(boxstyle='round,pad=0.15', facecolor='white', 
                                    edgecolor='red', alpha=0.9))
    
    # ============================================================
    # LAYER 5: Property boundary (green outline)
    # ============================================================
    prop_coords = data.get("boundaries", {}).get("property_line", {}).get("coordinates", [[]])[0]
    local_prop = None  # Initialize for later use in view limits
    if prop_coords:
        local_prop = to_local(prop_coords)
        poly = MplPolygon(local_prop, fill=False, edgecolor='darkgreen', linewidth=4, 
                         linestyle='-', zorder=5, label=f'Property Line ({prop_area:.1f}m²)')
        ax.add_patch(poly)
    
    # ============================================================
    # LAYER 6: Setback zone (orange dashed)
    # ============================================================
    setback_coords = data.get("boundaries", {}).get("setback_line", {}).get("coordinates", [[]])[0]
    setbacks = data.get("boundaries", {}).get("setback_line", {}).get("setbacks_applied", {})
    if setback_coords:
        local_setback = to_local(setback_coords)
        poly = MplPolygon(local_setback, fill=True, facecolor='moccasin', alpha=0.3, 
                          edgecolor='darkorange', linewidth=2, linestyle='--', zorder=5,
                          label=f'Setback Zone ({setback_area:.1f}m²)')
        ax.add_patch(poly)
        
        # Add setback distance labels at edges
        front = setbacks.get("front_m", 5)
        rear = setbacks.get("rear_m", 5)
        side_e = setbacks.get("side_east_m", 1)
        side_w = setbacks.get("side_west_m", 1)
    
    # ============================================================
    # LAYER 7: Buildable envelope (green fill)
    # ============================================================
    buildable_coords = data.get("boundaries", {}).get("buildable_envelope", {}).get("coordinates", [[]])[0]
    if buildable_coords:
        local_buildable = to_local(buildable_coords)
        poly = MplPolygon(local_buildable, fill=True, facecolor='palegreen', alpha=0.5,
                          edgecolor='green', linewidth=2, zorder=5,
                          label=f'Buildable Envelope ({buildable_area:.1f}m²)')
        ax.add_patch(poly)
    
    # ============================================================
    # LAYER 8: Existing structures (red hatched)
    # ============================================================
    existing = data.get("existing_structures", [])
    for structure in existing:
        s_coords = structure.get("footprint", {}).get("coordinates", [[]])[0]
        if s_coords:
            local_s = to_local(s_coords)
            poly = MplPolygon(local_s, fill=True, facecolor='lightcoral', alpha=0.7, 
                             edgecolor='darkred', linewidth=2, hatch='xx', zorder=6)
            ax.add_patch(poly)
            
            # Label existing structure
            center = get_polygon_center(local_s)
            status = structure.get("status", "existing")
            s_height = structure.get("height_m", 0)
            label = f"EXISTING\n({status})"
            if s_height > 0:
                label += f"\n{s_height}m"
            ax.annotate(label, center, fontsize=7, ha='center', va='center',
                       color='darkred', fontweight='bold', zorder=10)
    
    # ============================================================
    # LAYER 9: Access point
    # ============================================================
    access = data.get("access", {})
    primary_access = access.get("primary_access_point", {})
    access_loc = primary_access.get("location", {}).get("coordinates", [])
    if access_loc:
        local_access = to_local([access_loc])[0]
        ax.scatter(local_access[0], local_access[1], c='darkorange', s=200, 
                  marker='^', edgecolors='red', linewidths=2, zorder=8)
        ax.annotate('ACCESS', (local_access[0], local_access[1] + 5), 
                   fontsize=8, ha='center', color='red', fontweight='bold', zorder=10)
    
    # ============================================================
    # LAYER 10: Centroid marker
    # ============================================================
    ax.plot(0, 0, 'ro', markersize=12, markeredgecolor='darkred', markeredgewidth=2, 
           zorder=9, label='Centroid')
    
    # ============================================================
    # ANNOTATIONS: Scale bar, North arrow, Info box
    # ============================================================
    
    # Set view limits to exactly 60m x 60m box centered on centroid (0,0)
    # Add small buffer for labels and annotations
    buffer = 5.0  # 5m buffer for labels
    xlim = (-VIEW_BBOX - buffer, VIEW_BBOX + buffer)
    ylim = (-VIEW_BBOX - buffer, VIEW_BBOX + buffer)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Scale bar - positioned at bottom center of property area
    scale_length = 10  # 10 meters for better fit
    # Position below property, centered
    bar_y = ylim[0] + (ylim[1] - ylim[0]) * 0.08
    bar_x = -scale_length / 2  # Center it horizontally
    
    # Draw scale bar with end caps
    ax.plot([bar_x, bar_x + scale_length], [bar_y, bar_y], 'k-', linewidth=4, zorder=10)
    ax.plot([bar_x, bar_x], [bar_y - 1, bar_y + 1], 'k-', linewidth=2, zorder=10)  # Left cap
    ax.plot([bar_x + scale_length, bar_x + scale_length], [bar_y - 1, bar_y + 1], 'k-', linewidth=2, zorder=10)  # Right cap
    ax.annotate(f'{scale_length}m', (bar_x + scale_length/2, bar_y + 2), 
               fontsize=10, ha='center', va='bottom', fontweight='bold', zorder=10)
    
    # North arrow - positioned at top right area
    arrow_x = xlim[1] - (xlim[1] - xlim[0]) * 0.10
    arrow_y = ylim[1] - (ylim[1] - ylim[0]) * 0.10
    ax.annotate('N', (arrow_x, arrow_y + 4), fontsize=14, ha='center', fontweight='bold', 
               color='navy', zorder=10)
    ax.annotate('', xy=(arrow_x, arrow_y + 2), xytext=(arrow_x, arrow_y - 6),
               arrowprops=dict(arrowstyle='-|>', color='navy', lw=2, 
                              mutation_scale=15), zorder=10)
    
    # Info box (positioned at top left, outside the plot area)
    info_text = f"Plot ID: {data.get('plot_id', 'Unknown')}\n"
    info_text += f"Type: {data.get('plot_type', 'residential')}\n"
    info_text += f"Property: {prop_area:.0f} m²\n"
    info_text += f"Setback: {setback_area:.0f} m²\n"
    info_text += f"Buildable: {buildable_area:.0f} m²"
    
    # Add soil info if available
    soil = data.get("soil", {})
    if soil:
        info_text += f"\nSoil: {soil.get('type_id', 'unknown')}"
        info_text += f"\nBearing: {soil.get('bearing_capacity_kpa', 0)} kPa"
    
    # Add elevation if available
    elevation = data.get("surrounding_context", {}).get("elevation_map", {})
    if elevation:
        info_text += f"\nElevation: {elevation.get('average_elevation_m', 0):.0f}m"
        info_text += f"\nSlope: {elevation.get('slope_percent', 0):.1f}%"
    
    ax.text(xlim[0] + (xlim[1]-xlim[0])*0.02, ylim[1] - (ylim[1]-ylim[0])*0.02,
           info_text, transform=ax.transData, fontsize=8, va='top', ha='left',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'),
           zorder=10, family='monospace')
    
    # ============================================================
    # STYLING
    # ============================================================
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_xlabel('East-West (meters)', fontsize=10, fontweight='bold')
    ax.set_ylabel('North-South (meters)', fontsize=10, fontweight='bold')
    
    # Title with area info
    title = f"Plot Analysis: {data.get('plot_id', 'Unknown')}"
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Legend with areas - position outside plot
    legend_elements = []
    if prop_coords:
        legend_elements.append(Line2D([0], [0], color='darkgreen', linewidth=3, 
                              label=f'Property Line ({prop_area:.0f}m²)'))
    if setback_coords:
        legend_elements.append(patches.Patch(facecolor='moccasin', edgecolor='darkorange',
                              linewidth=2, linestyle='--', alpha=0.5,
                              label=f'Setback Zone ({setback_area:.0f}m²)'))
    if buildable_coords:
        legend_elements.append(patches.Patch(facecolor='palegreen', edgecolor='green',
                              linewidth=2, alpha=0.5,
                              label=f'Buildable Envelope ({buildable_area:.0f}m²)'))
    if buildings:
        legend_elements.append(patches.Patch(facecolor='slategray', edgecolor='darkslategray',
                              hatch='///', alpha=0.6, label='Neighbor Buildings'))
    if existing:
        legend_elements.append(patches.Patch(facecolor='lightcoral', edgecolor='darkred',
                              hatch='xx', alpha=0.7, label='Existing Structures'))
    if trees:
        legend_elements.append(plt.Circle((0,0), 1, facecolor='forestgreen', 
                              edgecolor='darkgreen', alpha=0.4, label=f'Tree Zones ({len(trees)})'))
    if roads:
        legend_elements.append(patches.Patch(facecolor='lightgray', edgecolor='dimgray',
                              alpha=0.8, label='Roads'))
    if water_features:
        legend_elements.append(patches.Patch(facecolor='lightblue', edgecolor='steelblue',
                              alpha=0.5, label='Water'))
    if access_loc:
        legend_elements.append(Line2D([0], [0], marker='^', color='w', markerfacecolor='darkorange',
                              markersize=12, markeredgecolor='red', markeredgewidth=1,
                              label='Access Point'))
    legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
                          markersize=10, markeredgecolor='darkred', label='Centroid'))
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, 
             framealpha=0.9, edgecolor='gray')
    
    # Save or show
    if args.output:
        plt.savefig(args.output, dpi=150, bbox_inches='tight', facecolor='white')
        logger.info(f"Saved visualization to: {args.output}")
    else:
        plt.show()
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Plot Analysis Generator CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Generate single plot:
    python cli.py generate --lat 51.5074 --lon -0.1276 --output my_plot.json
  
  Batch generate from CSV:
    python cli.py batch --input locations.csv --output ./plots/
  
  Visualize a plot:
    python cli.py visualize --input plot.json --output plot.png
        """
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Generate command
    gen_parser = subparsers.add_parser("generate", help="Generate plot analysis for a location")
    gen_parser.add_argument("--lat", type=float, required=True, help="Latitude")
    gen_parser.add_argument("--lon", type=float, required=True, help="Longitude")
    gen_parser.add_argument("--output", "-o", help="Output JSON file")
    gen_parser.add_argument("--radius", "-r", type=float, default=50, help="Search radius in meters")
    gen_parser.add_argument("--plot-id", help="Custom plot ID")
    gen_parser.add_argument("--summary", "-s", action="store_true", help="Print summary to stdout")
    gen_parser.set_defaults(func=cmd_generate)
    
    # Batch command
    batch_parser = subparsers.add_parser("batch", help="Batch generate from CSV file")
    batch_parser.add_argument("--input", "-i", required=True, help="Input CSV file (columns: name,lat,lon)")
    batch_parser.add_argument("--output", "-o", default="output", help="Output directory")
    batch_parser.add_argument("--delay", type=float, default=2.0, help="Delay between requests (seconds)")
    batch_parser.set_defaults(func=cmd_batch)
    
    # Visualize command
    viz_parser = subparsers.add_parser("visualize", help="Visualize a plot.json file")
    viz_parser.add_argument("--input", "-i", required=True, help="Input JSON file")
    viz_parser.add_argument("--output", "-o", help="Output image file (shows window if not specified)")
    viz_parser.set_defaults(func=cmd_visualize)
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

