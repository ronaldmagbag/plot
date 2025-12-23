"""
Debug script for Property Line Classifier

Loads existing plot data and tests only the classifier with detailed debug output
"""

import json
import sys
import math
from pathlib import Path
from loguru import logger
from typing import List, Dict, Any, Tuple

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.collectors.boundary.models import PropertyLine
from src.collectors.boundary.classifier import PropertyLineClassifier


def load_plot_data(json_path: str):
    """Load property line and roads from plot JSON file"""
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Extract property line coordinates
    boundaries = data.get("boundaries", {})
    property_line_data = boundaries.get("property_line", {})
    property_coords = property_line_data.get("coordinates", [[]])[0]
    
    if not property_coords:
        raise ValueError("No property line coordinates found in JSON")
    
    # Create PropertyLine object
    property_line = PropertyLine(
        coordinates=property_coords,
        area_sqm=property_line_data.get("area_sqm", 0),
        perimeter_m=property_line_data.get("perimeter_m", 0),
        source=property_line_data.get("source", "unknown"),
        accuracy_m=property_line_data.get("accuracy_m", 5.0)
    )
    
    # Extract roads
    surrounding_context = data.get("surrounding_context", {})
    roads_data = surrounding_context.get("roads", [])
    
    # Convert roads to the format expected by classifier
    osm_roads = []
    for road in roads_data:
        centerline = road.get("centerline", {})
        if centerline and centerline.get("coordinates"):
            osm_roads.append({
                "id": road.get("id", ""),
                "name": road.get("name", "Unknown Road"),
                "type": road.get("type", "residential"),
                "centerline": {
                    "coordinates": centerline.get("coordinates", [])
                }
            })
    
    return property_line, osm_roads


def main():
    """Main debug function"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Debug Property Line Classifier")
    parser.add_argument(
        "--input",
        type=str,
        default="output/plot_adur1/plot_adur1_1.json",
        help="Path to plot JSON file"
    )
    parser.add_argument(
        "--lat",
        type=float,
        help="Latitude for test (if not loading from JSON)"
    )
    parser.add_argument(
        "--lon",
        type=float,
        help="Longitude for test (if not loading from JSON)"
    )
    
    args = parser.parse_args()
    
    # Configure logger to show INFO level
    logger.remove()
    logger.add(
        sys.stdout,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
        level="INFO"
    )
    
    try:
        # Load data from JSON file
        logger.info(f"Loading plot data from: {args.input}")
        property_line, osm_roads = load_plot_data(args.input)
        
        logger.info(f"Property line: {len(property_line.coordinates)} points")
        logger.info(f"Roads found: {len(osm_roads)}")
        for road in osm_roads:
            logger.info(f"  - {road['name']} ({road['type']})")
        logger.info("")
        
        # Create classifier
        classifier = PropertyLineClassifier()
        
        # Run classification with debug output
        logger.info("=" * 60)
        logger.info("RUNNING PROPERTY LINE CLASSIFICATION")
        logger.info("=" * 60)
        logger.info("")
        
        classified_property_line = classifier.classify(
            property_line,
            osm_roads,
            []  # Empty buildings list
        )
        
        logger.info("")
        logger.info("=" * 60)
        logger.info("CLASSIFICATION RESULTS")
        logger.info("=" * 60)
        
        # Print segment summaries
        if classified_property_line.front:
            coords = classified_property_line.front.get_coordinates(classified_property_line.coordinates)
            logger.info(f"Front segment: {classified_property_line.front.length_m:.2f}m, "
                       f"{len(classified_property_line.front.edge_indices)} edges, "
                       f"{len(coords)} points")
        else:
            logger.info("Front segment: None")
        
        if classified_property_line.rear:
            coords = classified_property_line.rear.get_coordinates(classified_property_line.coordinates)
            logger.info(f"Rear segment: {classified_property_line.rear.length_m:.2f}m, "
                       f"{len(classified_property_line.rear.edge_indices)} edges, "
                       f"{len(coords)} points")
        else:
            logger.info("Rear segment: None")
        
        if classified_property_line.left_side:
            coords = classified_property_line.left_side.get_coordinates(classified_property_line.coordinates)
            logger.info(f"Left side segment: {classified_property_line.left_side.length_m:.2f}m, "
                       f"{len(classified_property_line.left_side.edge_indices)} edges, "
                       f"{len(coords)} points")
        else:
            logger.info("Left side segment: None")
        
        if classified_property_line.right_side:
            coords = classified_property_line.right_side.get_coordinates(classified_property_line.coordinates)
            logger.info(f"Right side segment: {classified_property_line.right_side.length_m:.2f}m, "
                       f"{len(classified_property_line.right_side.edge_indices)} edges, "
                       f"{len(coords)} points")
        else:
            logger.info("Right side segment: None")
        
        logger.info("")
        logger.info("Classification complete!")
        
        # Visualize if matplotlib is available
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            from matplotlib.patches import Polygon as MplPolygon, FancyBboxPatch
            from matplotlib.lines import Line2D
            import numpy as np
            
            visualize_classification(
                classified_property_line,
                osm_roads,
                args.input.replace('.json', '_debug.png')
            )
            logger.info(f"Visualization saved to: {args.input.replace('.json', '_debug.png')}")
        except ImportError:
            logger.warning("matplotlib not available - skipping visualization")
        except Exception as e:
            logger.warning(f"Visualization failed: {e}")
            import traceback
            traceback.print_exc()
        
    except FileNotFoundError:
        logger.error(f"File not found: {args.input}")
        logger.info("Available files in output/:")
        output_dir = Path("output")
        if output_dir.exists():
            for json_file in output_dir.rglob("*.json"):
                logger.info(f"  {json_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def visualize_classification(
    property_line: PropertyLine,
    osm_roads: List[Dict[str, Any]],
    output_path: str
):
    """Visualize property line classification with edge debug details"""
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.lines import Line2D
    import numpy as np
    plt  # Make sure plt is available
    
    # Get reference point (centroid of property)
    coords = property_line.coordinates
    if not coords:
        logger.error("No coordinates in property line")
        return
    
    # Calculate centroid
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    ref_lon = sum(xs) / len(xs)
    ref_lat = sum(ys) / len(ys)
    
    def to_local(coords_list):
        """Convert lat/lon to local meters"""
        m_per_deg_lat = 111000
        m_per_deg_lon = 111000 * np.cos(np.radians(ref_lat))
        return [((c[0] - ref_lon) * m_per_deg_lon, (c[1] - ref_lat) * m_per_deg_lat) 
                for c in coords_list]
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(16, 14))
    
    # Convert coordinates to local
    local_prop = to_local(coords)
    
    # Determine view bounds (60m box around property)
    all_xs = [p[0] for p in local_prop]
    all_ys = [p[1] for p in local_prop]
    for road in osm_roads:
        road_coords = road.get("centerline", {}).get("coordinates", [])
        if road_coords:
            local_road = to_local(road_coords)
            all_xs.extend([p[0] for p in local_road])
            all_ys.extend([p[1] for p in local_road])
    
    margin = 10
    x_min, x_max = min(all_xs) - margin, max(all_xs) + margin
    y_min, y_max = min(all_ys) - margin, max(all_ys) + margin
    
    # Draw merged road (all segments with same name)
    edge_debug_info = property_line.metadata.get("edge_debug_info", [])
    closest_road_name = property_line.metadata.get("closest_road_name", "Unknown")
    
    # Find and merge all roads with the same name
    merged_road_coords = []
    for road in osm_roads:
        if road.get("name") == closest_road_name:
            road_coords = road.get("centerline", {}).get("coordinates", [])
            if road_coords:
                merged_road_coords.extend(road_coords)
    
    # Remove duplicate consecutive points from merged road
    if merged_road_coords:
        cleaned_road_coords = [merged_road_coords[0]]
        for i in range(1, len(merged_road_coords)):
            prev = cleaned_road_coords[-1]
            curr = merged_road_coords[i]
            # If points are different (more than 1e-6 degrees apart), add it
            if abs(prev[0] - curr[0]) > 1e-6 or abs(prev[1] - curr[1]) > 1e-6:
                cleaned_road_coords.append(curr)
        
        # Draw merged road
        if len(cleaned_road_coords) >= 2:
            local_road = to_local(cleaned_road_coords)
            xs, ys = zip(*local_road)
            ax.plot(xs, ys, color='gray', linewidth=4, alpha=0.7, zorder=1)
    
    # Draw property line with colored segments (no labels)
    # Reconstruct coordinates from edge indices
    if property_line.front:
        front_coords = property_line.front.get_coordinates(property_line.coordinates)
        if len(front_coords) >= 2:
            front_coords = to_local(front_coords)
            xs, ys = zip(*front_coords)
            ax.plot(xs, ys, color=property_line.front.color, linewidth=6, 
                   linestyle='-', zorder=5, alpha=0.9)
    
    if property_line.rear:
        rear_coords = property_line.rear.get_coordinates(property_line.coordinates)
        if len(rear_coords) >= 2:
            rear_coords = to_local(rear_coords)
            xs, ys = zip(*rear_coords)
            ax.plot(xs, ys, color=property_line.rear.color, linewidth=6, 
                   linestyle='-', zorder=5, alpha=0.9)
    
    if property_line.left_side:
        left_coords = property_line.left_side.get_coordinates(property_line.coordinates)
        if len(left_coords) >= 2:
            left_coords = to_local(left_coords)
            xs, ys = zip(*left_coords)
            ax.plot(xs, ys, color=property_line.left_side.color, linewidth=6, 
                   linestyle='-', zorder=5, alpha=0.9)
    
    if property_line.right_side:
        right_coords = property_line.right_side.get_coordinates(property_line.coordinates)
        if len(right_coords) >= 2:
            right_coords = to_local(right_coords)
            xs, ys = zip(*right_coords)
            ax.plot(xs, ys, color=property_line.right_side.color, linewidth=6, 
                   linestyle='-', zorder=5, alpha=0.9)
    
    # Set axis - clean, no labels
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect('equal')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    logger.info(f"Saved visualization to {output_path}")
    plt.close()


if __name__ == "__main__":
    main()

