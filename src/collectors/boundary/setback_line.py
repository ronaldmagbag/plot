"""
Setback line calculation

Calculates setback boundaries from property line
- 4m setback from front and rear lines
- 1m setback from side lines
"""

import math
from typing import List, Dict, Any, Optional, Tuple
from loguru import logger

try:
    from shapely.geometry import Polygon, Point, LineString
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - setback line calculation will be limited")

from .models import SetbackLine, PropertyLine, SegmentType
from .utils import calculate_polygon_area, calculate_perimeter, close_polygon
from ...config import get_config


class SetbackLineProcessor:
    """Processes setback line calculation"""
    
    def __init__(self):
        self.config = get_config()
        self.front_setback_m = self.config.uk_regulatory.front_setback_m
        self.rear_setback_m = self.config.uk_regulatory.rear_setback_m
        self.side_setback_m = self.config.uk_regulatory.side_setback_m
    
    def calculate_setback_line(
        self,
        property_line: PropertyLine
    ) -> Optional[SetbackLine]:
        """
        Calculate setback line from property line
        
        Setback rules (configurable in config.py):
        - front_setback_m (default 4m) from front edges
        - rear_setback_m (default 3m) from rear edges
        - side_setback_m (default 1m) from side edges
        
        Uses classified segments to apply appropriate setback to each edge.
        
        Args:
            property_line: PropertyLine object with classified segments
            
        Returns:
            SetbackLine object or None
        """
        if not SHAPELY_AVAILABLE:
            logger.warning("Shapely not available - cannot calculate setback line")
            return None
        
        if not property_line.coordinates or len(property_line.coordinates) < 4:
            logger.warning("Invalid property line coordinates")
            return None
        
        try:
            # Create property polygon
            prop_coords = close_polygon(property_line.coordinates)
            prop_poly = Polygon(prop_coords)
            
            if not prop_poly.is_valid:
                logger.warning("Invalid property polygon")
                return None
            
            # Calculate dynamic rear setback if rear_setback_m is 0.0
            actual_rear_setback_m = self.rear_setback_m
            if self.rear_setback_m == 0.0:
                # Calculate dynamically: round(average_side_length/2) - 5
                side_lengths = []
                if property_line.left_side:
                    side_lengths.append(property_line.left_side.length_m)
                if property_line.right_side:
                    side_lengths.append(property_line.right_side.length_m)
                
                if len(side_lengths) > 0:
                    # Calculate average of side lengths
                    average_side_length = sum(side_lengths) / len(side_lengths)
                    actual_rear_setback_m = round(average_side_length / 2.0) - 5.0
                    # Ensure minimum of 0
                    actual_rear_setback_m = max(0.0, actual_rear_setback_m)
                    logger.info(f"Dynamic rear setback calculated: average_side_length={average_side_length:.2f}m, rear_setback={actual_rear_setback_m:.2f}m (round({average_side_length:.2f}/2) - 5)")
                else:
                    logger.warning("Cannot calculate dynamic rear setback: no side edges found, using 0")
                    actual_rear_setback_m = 0.0
            else:
                logger.info(f"Using static rear setback: {actual_rear_setback_m:.2f}m")
            
            # Calculate setback polygon
            setback_coords = self._calculate_setback_polygon(
                prop_poly,
                property_line,
                actual_rear_setback_m
            )
            
            if not setback_coords or len(setback_coords) < 4:
                logger.warning("Failed to calculate setback polygon")
                return None
            
            # Calculate area and perimeter
            center_lat = sum(c[1] for c in setback_coords) / len(setback_coords)
            area = calculate_polygon_area(setback_coords, center_lat)
            perimeter = calculate_perimeter(setback_coords, center_lat)
            
            return SetbackLine(
                coordinates=setback_coords,
                area_sqm=round(area, 1),
                perimeter_m=round(perimeter, 1),
                setback_type="full",
                metadata={
                    "front_setback_m": self.front_setback_m,
                    "rear_setback_m": actual_rear_setback_m,
                    "side_setback_m": self.side_setback_m
                }
            )
            
        except Exception as e:
            logger.warning(f"Failed to calculate setback line: {e}")
            return None
    
    def _calculate_setback_polygon(
        self,
        prop_poly: Polygon,
        property_line: PropertyLine,
        rear_setback_m: Optional[float] = None
    ) -> Optional[List[List[float]]]:
        """
        Calculate setback polygon by offsetting edges with variable setbacks
        
        Applies different setbacks per edge based on classification:
        - Front edges: front_setback_m (default 4m)
        - Rear edges: rear_setback_m (default 3m)
        - Side edges: side_setback_m (default 1m)
        
        Strategy:
        1. Identify which edges belong to front/rear/sides
        2. Offset each edge inward by its appropriate distance
        3. Reconstruct polygon from offset edges
        """
        try:
            # Get property coordinates
            prop_coords = close_polygon(property_line.coordinates)
            num_edges = len(prop_coords) - 1  # Exclude closing point
            
            # Build edge classification map
            edge_classification = {}  # edge_index -> setback_distance_m
            
            # Classify edges based on segments
            if property_line.front:
                for edge_idx in property_line.front.edge_indices:
                    edge_classification[edge_idx] = self.front_setback_m
            
            # Use provided rear_setback_m if given, otherwise use self.rear_setback_m
            rear_setback = rear_setback_m if rear_setback_m is not None else self.rear_setback_m
            
            if property_line.rear:
                for edge_idx in property_line.rear.edge_indices:
                    edge_classification[edge_idx] = rear_setback
            
            if property_line.left_side:
                for edge_idx in property_line.left_side.edge_indices:
                    if edge_idx not in edge_classification:  # Don't override if already classified
                        edge_classification[edge_idx] = self.side_setback_m
            
            if property_line.right_side:
                for edge_idx in property_line.right_side.edge_indices:
                    if edge_idx not in edge_classification:  # Don't override if already classified
                        edge_classification[edge_idx] = self.side_setback_m
            
            # Default to side setback for unclassified edges
            for i in range(num_edges):
                if i not in edge_classification:
                    edge_classification[i] = self.side_setback_m
            
            # Convert meters to degrees
            center_lat = prop_poly.centroid.y
            m_per_deg_lat = 111000
            m_per_deg_lon = 111000 * abs(math.cos(math.radians(center_lat)))
            m_per_deg = (m_per_deg_lat + m_per_deg_lon) / 2
            
            # Determine polygon orientation (counter-clockwise = positive area)
            # For counter-clockwise: inward is to the LEFT of edge direction (interior is left when walking boundary)
            # For clockwise: inward is to the RIGHT of edge direction (interior is right when walking boundary)
            # Use exterior.is_ccw which is more reliable than area in geographic coords
            try:
                is_ccw = prop_poly.exterior.is_ccw
            except AttributeError:
                # Fallback to area check if is_ccw not available
                polygon_area = prop_poly.area
                is_ccw = polygon_area > 0
            
            # Verify orientation by checking if centroid is inside
            # This is a sanity check to catch orientation errors
            centroid = prop_poly.centroid
            if not prop_poly.contains(centroid):
                logger.warning(f"Property polygon centroid is outside polygon - orientation may be wrong")
                logger.warning(f"  is_ccw={is_ccw}, polygon area (deg²)={prop_poly.area}")
            
            # Offset each edge and collect offset points
            offset_points = []
            
            for i in range(num_edges):
                start_point = Point(prop_coords[i])
                end_point = Point(prop_coords[i + 1])
                
                # Get setback distance for this edge
                setback_m = edge_classification.get(i, self.side_setback_m)
                setback_deg = setback_m / m_per_deg
                
                # Calculate edge vector
                dx = end_point.x - start_point.x
                dy = end_point.y - start_point.y
                edge_length = math.sqrt(dx**2 + dy**2)
                
                if edge_length < 1e-10:
                    # Degenerate edge, skip
                    continue
                
                # Normalize edge direction
                dx_norm = dx / edge_length
                dy_norm = dy / edge_length
                
                # Calculate perpendicular vector (pointing inward)
                # For counter-clockwise polygon: right of edge = (-dy, dx)
                # For clockwise polygon: left of edge = (dy, -dx)
                if is_ccw:
                    perp_x = -dy_norm
                    perp_y = dx_norm
                else:
                    perp_x = dy_norm
                    perp_y = -dx_norm
                
                # Offset both start and end points inward
                start_offset = Point(
                    start_point.x + perp_x * setback_deg,
                    start_point.y + perp_y * setback_deg
                )
                end_offset = Point(
                    end_point.x + perp_x * setback_deg,
                    end_point.y + perp_y * setback_deg
                )
                
                offset_points.append((i, start_offset, end_offset, setback_m))
            
            # Reconstruct polygon from offset points
            # Connect offset points, handling corners where edges meet
            setback_coords = []
            
            for i in range(len(offset_points)):
                curr_idx, curr_start, curr_end, curr_setback = offset_points[i]
                next_idx, next_start, next_end, next_setback = offset_points[(i + 1) % len(offset_points)]
                
                # Use the end point of current edge
                # At corners, we need to find intersection of two offset lines
                if i == 0:
                    # First point: use start of first edge
                    setback_coords.append([curr_start.x, curr_start.y])
                
                # Calculate corner intersection
                # Line 1: from curr_end to curr_end + edge_direction
                # Line 2: from next_start to next_start + edge_direction
                # Find intersection of offset lines
                
                # Get edge directions
                curr_edge_start = Point(prop_coords[curr_idx])
                curr_edge_end = Point(prop_coords[curr_idx + 1])
                next_edge_start = Point(prop_coords[next_idx])
                next_edge_end = Point(prop_coords[next_idx + 1])
                
                curr_dx = curr_edge_end.x - curr_edge_start.x
                curr_dy = curr_edge_end.y - curr_edge_start.y
                curr_len = math.sqrt(curr_dx**2 + curr_dy**2)
                
                next_dx = next_edge_end.x - next_edge_start.x
                next_dy = next_edge_end.y - next_edge_start.y
                next_len = math.sqrt(next_dx**2 + next_dy**2)
                
                if curr_len > 1e-10 and next_len > 1e-10:
                    # Normalize
                    curr_dx /= curr_len
                    curr_dy /= curr_len
                    next_dx /= next_len
                    next_dy /= next_len
                    
                    # Perpendicular vectors (inward) - use same orientation logic
                    if is_ccw:
                        curr_perp_x = -curr_dy
                        curr_perp_y = curr_dx
                        next_perp_x = -next_dy
                        next_perp_y = next_dx
                    else:
                        curr_perp_x = curr_dy
                        curr_perp_y = -curr_dx
                        next_perp_x = next_dy
                        next_perp_y = -next_dx
                    
                    # Offset distances
                    curr_offset_deg = edge_classification.get(curr_idx, self.side_setback_m) / m_per_deg
                    next_offset_deg = edge_classification.get(next_idx, self.side_setback_m) / m_per_deg
                    
                    # Offset line 1: through curr_end
                    line1_p1 = Point(
                        curr_edge_end.x + curr_perp_x * curr_offset_deg,
                        curr_edge_end.y + curr_perp_y * curr_offset_deg
                    )
                    line1_p2 = Point(
                        line1_p1.x + curr_dx,
                        line1_p1.y + curr_dy
                    )
                    
                    # Offset line 2: through next_start
                    line2_p1 = Point(
                        next_edge_start.x + next_perp_x * next_offset_deg,
                        next_edge_start.y + next_perp_y * next_offset_deg
                    )
                    line2_p2 = Point(
                        line2_p1.x + next_dx,
                        line2_p1.y + next_dy
                    )
                    
                    # Find intersection of two lines
                    intersection = self._line_intersection(
                        line1_p1, line1_p2, line2_p1, line2_p2
                    )
                    
                    if intersection:
                        # CRITICAL: Verify intersection is inside the property polygon
                        # If not, use a safer fallback (offset the original corner point)
                        intersection_point = Point(intersection.x, intersection.y)
                        if prop_poly.contains(intersection_point) or prop_poly.touches(intersection_point):
                            setback_coords.append([intersection.x, intersection.y])
                        else:
                            # Intersection is outside property - use corner point offset by minimum setback
                            corner_point = curr_edge_end  # This is the actual corner
                            min_setback = min(curr_setback, next_setback)
                            min_setback_deg = min_setback / m_per_deg
                            # Use average perpendicular direction for corner
                            avg_perp_x = (curr_perp_x + next_perp_x) / 2
                            avg_perp_y = (curr_perp_y + next_perp_y) / 2
                            # Normalize
                            avg_perp_len = math.sqrt(avg_perp_x**2 + avg_perp_y**2)
                            if avg_perp_len > 1e-10:
                                avg_perp_x /= avg_perp_len
                                avg_perp_y /= avg_perp_len
                            safe_corner = Point(
                                corner_point.x + avg_perp_x * min_setback_deg,
                                corner_point.y + avg_perp_y * min_setback_deg
                            )
                            setback_coords.append([safe_corner.x, safe_corner.y])
                            logger.debug(f"Corner intersection outside property, using safe corner offset (edge {curr_idx})")
                    else:
                        # Lines parallel - use corner point offset by minimum setback
                        corner_point = curr_edge_end
                        min_setback = min(curr_setback, next_setback)
                        min_setback_deg = min_setback / m_per_deg
                        # Use average perpendicular direction
                        avg_perp_x = (curr_perp_x + next_perp_x) / 2
                        avg_perp_y = (curr_perp_y + next_perp_y) / 2
                        avg_perp_len = math.sqrt(avg_perp_x**2 + avg_perp_y**2)
                        if avg_perp_len > 1e-10:
                            avg_perp_x /= avg_perp_len
                            avg_perp_y /= avg_perp_len
                        safe_corner = Point(
                            corner_point.x + avg_perp_x * min_setback_deg,
                            corner_point.y + avg_perp_y * min_setback_deg
                        )
                        setback_coords.append([safe_corner.x, safe_corner.y])
                        logger.debug(f"Parallel lines at corner, using safe corner offset (edge {curr_idx})")
                else:
                    # Degenerate edge - use the already-calculated offset point
                    # curr_end and next_start are the offset Points from offset_points
                    # Use the midpoint of these offset points as a safe fallback
                    mid_x = (curr_end.x + next_start.x) / 2
                    mid_y = (curr_end.y + next_start.y) / 2
                    setback_coords.append([mid_x, mid_y])
                    logger.debug(f"Degenerate edge at corner, using midpoint of offset points (edge {curr_idx})")
            
            if len(setback_coords) < 3:
                logger.warning("Not enough offset points for setback polygon")
                return None
            
            # Create polygon from offset coordinates
            setback_poly = Polygon(setback_coords)
            
            # Validate
            if not setback_poly.is_valid:
                # Try to fix - buffer(0) might return MultiPolygon
                fixed = setback_poly.buffer(0)
                
                # Handle MultiPolygon case
                if hasattr(fixed, 'geoms'):
                    # It's a MultiPolygon - take the largest polygon
                    largest = max(fixed.geoms, key=lambda g: g.area if hasattr(g, 'area') else 0)
                    setback_poly = largest if isinstance(largest, Polygon) else fixed
                else:
                    setback_poly = fixed
            
            # Check if it's still a MultiPolygon
            if hasattr(setback_poly, 'geoms'):
                # Extract largest polygon from MultiPolygon
                largest = max(setback_poly.geoms, key=lambda g: g.area if hasattr(g, 'area') else 0)
                if isinstance(largest, Polygon):
                    setback_poly = largest
                else:
                    logger.warning("Setback polygon is MultiPolygon and cannot extract single polygon")
                    return None
            
            if not isinstance(setback_poly, Polygon) or not setback_poly.is_valid or setback_poly.area <= 0:
                logger.warning("Invalid setback polygon after offsetting")
                return None
            
            # CRITICAL FIX: Clip setback polygon to property boundary to ensure containment
            # This handles cases where corner intersections or buffer(0) created geometry outside
            try:
                clipped_setback = setback_poly.intersection(prop_poly)
                
                # Handle MultiPolygon result
                if hasattr(clipped_setback, 'geoms'):
                    # Take the largest polygon from intersection
                    clipped_setback = max(clipped_setback.geoms, key=lambda g: g.area if hasattr(g, 'area') else 0)
                
                if not isinstance(clipped_setback, Polygon) or not clipped_setback.is_valid:
                    logger.warning("Clipped setback polygon is invalid")
                    return None
                
                # Check if clipped polygon is too small (less than 10% of property area)
                clipped_area_sqm = calculate_polygon_area(
                    [[c[0], c[1]] for c in clipped_setback.exterior.coords[:-1]], 
                    center_lat
                )
                prop_area_sqm = calculate_polygon_area(
                    [[c[0], c[1]] for c in prop_poly.exterior.coords[:-1]], 
                    center_lat
                )
                
                if clipped_area_sqm < prop_area_sqm * 0.1:
                    logger.warning(f"Clipped setback area ({clipped_area_sqm:.2f} m²) is too small (< 10% of property {prop_area_sqm:.2f} m²)")
                    return None
                
                # Use clipped polygon
                setback_poly = clipped_setback
                logger.debug(f"Setback polygon clipped to property boundary: {clipped_area_sqm:.2f} m²")
                
            except Exception as e:
                logger.warning(f"Failed to clip setback polygon: {e}")
                # Fallback: verify containment
                if not prop_poly.covers(setback_poly):
                    setback_coords_list = [[c[0], c[1]] for c in setback_poly.exterior.coords[:-1]]
                    prop_coords_list = [[c[0], c[1]] for c in prop_poly.exterior.coords[:-1]]
                    setback_area_sqm = calculate_polygon_area(setback_coords_list, center_lat)
                    prop_area_sqm = calculate_polygon_area(prop_coords_list, center_lat)
                    logger.warning(f"Setback polygon is not fully contained within property boundary")
                    logger.warning(f"  Setback area: {setback_area_sqm:.2f} m², Property area: {prop_area_sqm:.2f} m²")
                    logger.warning(f"  This indicates the offset calculation produced invalid geometry - likely corner intersection error")
                    return None
            
            # Additional area check with tolerance for floating point precision
            setback_coords_list = [[c[0], c[1]] for c in setback_poly.exterior.coords[:-1]]
            prop_coords_list = [[c[0], c[1]] for c in prop_poly.exterior.coords[:-1]]
            
            setback_area_sqm = calculate_polygon_area(setback_coords_list, center_lat)
            prop_area_sqm = calculate_polygon_area(prop_coords_list, center_lat)
            
            # Allow small tolerance (0.5%) for floating point precision
            tolerance = max(prop_area_sqm * 0.005, 0.1)  # At least 0.1 m² or 0.5% of property area
            
            if setback_area_sqm >= (prop_area_sqm - tolerance):
                logger.warning(f"Setback area ({setback_area_sqm:.2f} m²) >= property area ({prop_area_sqm:.2f} m², tolerance: {tolerance:.2f} m²)")
                return None
            
            # Extract coordinates
            coords = list(setback_poly.exterior.coords)
            if len(coords) > 1 and coords[0] == coords[-1]:
                coords = coords[:-1]
            
            return [[c[0], c[1]] for c in coords]
            
        except Exception as e:
            logger.warning(f"Error calculating setback polygon: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def _line_intersection(
        self,
        p1: Point, p2: Point,
        p3: Point, p4: Point
    ) -> Optional[Point]:
        """
        Find intersection point of two line segments
        
        Returns None if lines are parallel or don't intersect
        """
        try:
            # Line 1: p1 -> p2
            # Line 2: p3 -> p4
            
            # Parametric form:
            # Line 1: x = p1.x + t*(p2.x - p1.x), y = p1.y + t*(p2.y - p1.y)
            # Line 2: x = p3.x + s*(p4.x - p3.x), y = p3.y + s*(p4.y - p3.y)
            
            dx1 = p2.x - p1.x
            dy1 = p2.y - p1.y
            dx2 = p4.x - p3.x
            dy2 = p4.y - p3.y
            
            # Solve for t and s
            # p1.x + t*dx1 = p3.x + s*dx2
            # p1.y + t*dy1 = p3.y + s*dy2
            
            denominator = dx1 * dy2 - dy1 * dx2
            
            if abs(denominator) < 1e-10:
                # Lines are parallel
                return None
            
            t = ((p3.x - p1.x) * dy2 - (p3.y - p1.y) * dx2) / denominator
            
            # Calculate intersection point
            intersection_x = p1.x + t * dx1
            intersection_y = p1.y + t * dy1
            
            return Point(intersection_x, intersection_y)
            
        except Exception as e:
            logger.debug(f"Line intersection calculation failed: {e}")
            return None
    
