"""
Geometry utility functions

Common geometry calculations for boundary processing
"""

import math
from typing import List, Tuple, Optional
import numpy as np
from loguru import logger


def haversine_distance(
    lat1: float,
    lon1: float,
    lat2: float,
    lon2: float
) -> float:
    """Calculate distance between two points in meters"""
    R = 6371000  # Earth radius in meters
    
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)
    
    a = math.sin(delta_phi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    
    return R * c


def calculate_polygon_area(
    coords: List[List[float]],
    center_lat: float
) -> float:
    """Calculate polygon area in square meters using shoelace formula"""
    if len(coords) < 3:
        return 0.0
    
    # Convert to local meters
    m_per_deg_lat = 111000
    m_per_deg_lon = 111000 * math.cos(math.radians(center_lat))
    
    # Convert coordinates to meters
    ref_lon = coords[0][0]
    ref_lat = coords[0][1]
    
    meters_coords = []
    for lon, lat in coords:
        x = (lon - ref_lon) * m_per_deg_lon
        y = (lat - ref_lat) * m_per_deg_lat
        meters_coords.append((x, y))
    
    # Shoelace formula
    n = len(meters_coords)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += meters_coords[i][0] * meters_coords[j][1]
        area -= meters_coords[j][0] * meters_coords[i][1]
    
    return abs(area) / 2.0


def calculate_perimeter(
    coords: List[List[float]],
    center_lat: float
) -> float:
    """Calculate polygon perimeter in meters"""
    if len(coords) < 2:
        return 0.0
    
    perimeter = 0.0
    for i in range(len(coords) - 1):
        perimeter += haversine_distance(
            coords[i][1], coords[i][0],
            coords[i+1][1], coords[i+1][0]
        )
    
    return perimeter


def point_roughly_in_polygon(
    x: float,
    y: float,
    polygon: List[List[float]]
) -> bool:
    """Simple point-in-polygon check (ray casting)"""
    n = len(polygon)
    inside = False
    
    j = n - 1
    for i in range(n):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
            inside = not inside
        j = i
    
    return inside


def calculate_line_length(
    coords: List[List[float]]
) -> float:
    """Calculate length of a line string in meters"""
    if len(coords) < 2:
        return 0.0
    
    length = 0.0
    for i in range(len(coords) - 1):
        length += haversine_distance(
            coords[i][1], coords[i][0],
            coords[i+1][1], coords[i+1][0]
        )
    
    return length


def close_polygon(coords: List[List[float]]) -> List[List[float]]:
    """Ensure polygon is closed (first point == last point)"""
    if not coords:
        return coords
    
    if coords[0] != coords[-1]:
        return coords + [coords[0]]
    
    return coords


def get_polygon_centroid(coords: List[List[float]]) -> Tuple[float, float]:
    """Get centroid of polygon"""
    if not coords:
        return (0.0, 0.0)
    
    # Remove duplicate closing point if present
    coords_clean = coords[:-1] if coords[0] == coords[-1] and len(coords) > 1 else coords
    
    lon_sum = sum(c[0] for c in coords_clean)
    lat_sum = sum(c[1] for c in coords_clean)
    
    return (lon_sum / len(coords_clean), lat_sum / len(coords_clean))


def angle_deg(a: Tuple[float, float], b: Tuple[float, float], c: Tuple[float, float]) -> float:
    """
    Calculate angle ABC in degrees
    
    Args:
        a: First point (x, y)
        b: Vertex point (x, y)
        c: Third point (x, y)
    
    Returns:
        Angle in degrees (0-180)
    """
    ba = np.array(a) - np.array(b)
    bc = np.array(c) - np.array(b)
    
    # Handle zero-length vectors
    norm_ba = np.linalg.norm(ba)
    norm_bc = np.linalg.norm(bc)
    
    if norm_ba < 1e-10 or norm_bc < 1e-10:
        return 180.0  # Consider collinear if vector is zero-length
    
    cos_angle = np.dot(ba, bc) / (norm_ba * norm_bc)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    return math.degrees(math.acos(cos_angle))


def point_line_distance(
    p: Tuple[float, float], 
    a: Tuple[float, float], 
    b: Tuple[float, float],
    center_lat: Optional[float] = None
) -> float:
    """
    Calculate perpendicular distance from point p to line segment ab
    
    For lat/lon coordinates, this calculates distance in degrees.
    For accurate meter-based distance, use center_lat to scale the result.
    
    Args:
        p: Point (lon, lat) or (x, y)
        a: Line start point (lon, lat) or (x, y)
        b: Line end point (lon, lat) or (x, y)
        center_lat: Optional latitude for scaling (if provided, scales lon component)
    
    Returns:
        Distance in same units as coordinates (degrees for lat/lon)
    """
    p_arr = np.array(p)
    a_arr = np.array(a)
    b_arr = np.array(b)
    
    # If center_lat is provided, scale longitude component for lat/lon coordinates
    if center_lat is not None:
        # Scale longitude by cos(lat) to account for Earth's curvature
        lon_scale = math.cos(math.radians(center_lat))
        p_arr_scaled = np.array([p_arr[0] * lon_scale, p_arr[1]])
        a_arr_scaled = np.array([a_arr[0] * lon_scale, a_arr[1]])
        b_arr_scaled = np.array([b_arr[0] * lon_scale, b_arr[1]])
        
        # Vector from a to b (scaled)
        ab = b_arr_scaled - a_arr_scaled
        ab_norm = np.linalg.norm(ab)
        
        if ab_norm < 1e-10:
            return np.linalg.norm(p_arr_scaled - a_arr_scaled)
        
        # Vector from a to p (scaled)
        ap = p_arr_scaled - a_arr_scaled
        
        # Project ap onto ab
        t = np.dot(ap, ab) / (ab_norm * ab_norm)
        t = np.clip(t, 0.0, 1.0)
        
        # Closest point on line segment (scaled)
        closest_scaled = a_arr_scaled + t * ab
        
        # Convert back to original coordinates for distance calculation
        closest = np.array([closest_scaled[0] / lon_scale, closest_scaled[1]])
        
        # Calculate distance in original coordinate space
        # But account for lon scaling in the distance
        dlon = (p_arr[0] - closest[0]) * lon_scale
        dlat = p_arr[1] - closest[1]
        return math.sqrt(dlon*dlon + dlat*dlat)
    else:
        # Original Euclidean distance calculation
        ab = b_arr - a_arr
        ab_norm = np.linalg.norm(ab)
        
        if ab_norm < 1e-10:
            return np.linalg.norm(p_arr - a_arr)
        
        ap = p_arr - a_arr
        t = np.dot(ap, ab) / (ab_norm * ab_norm)
        t = np.clip(t, 0.0, 1.0)
        
        closest = a_arr + t * ab
        return np.linalg.norm(p_arr - closest)


def simplify_polygon_by_angle(
    coords: List[List[float]], 
    angle_threshold: float = 179.0
) -> List[List[float]]:
    """
    Simplify polygon by removing collinear points (angle-based)
    
    Removes points where the internal angle is close to 180° (collinear)
    
    Args:
        coords: List of [x, y] or [lon, lat] coordinate pairs
        angle_threshold: Angle threshold in degrees (default 179.0)
                        Points with angle >= threshold are removed
    
    Returns:
        Simplified coordinate list
    """
    if len(coords) < 3:
        return coords
    
    # Handle closed polygon
    closed = False
    if len(coords) > 1 and coords[0] == coords[-1]:
        closed = True
        coords = coords[:-1]
    
    if len(coords) < 3:
        return coords + ([coords[0]] if closed else [])
    
    simplified = []
    n = len(coords)
    removed_count = 0
    
    for i in range(n):
        prev = tuple(coords[i - 1])
        curr = tuple(coords[i])
        nxt = tuple(coords[(i + 1) % n])
        
        ang = angle_deg(prev, curr, nxt)
        
        # Keep point if angle is less than threshold (not collinear)
        if ang < angle_threshold:
            simplified.append(coords[i])
        else:
            removed_count += 1
            logger.debug(f"Removing point {i} (angle={ang:.2f}° >= {angle_threshold}°)")
    
    if removed_count > 0:
        logger.debug(f"Angle simplification: removed {removed_count} points")
    
    if closed and simplified:
        # Ensure polygon is closed
        if simplified[0] != simplified[-1]:
            simplified.append(simplified[0])
    
    return simplified


def simplify_polygon_by_distance(
    coords: List[List[float]], 
    distance_threshold: float = 0.02,
    center_lat: Optional[float] = None
) -> List[List[float]]:
    """
    Simplify polygon by removing points close to line segments (distance-based)
    
    Removes points that are within distance_threshold of the line between neighbors
    
    Args:
        coords: List of [x, y] or [lon, lat] coordinate pairs
        distance_threshold: Distance threshold in same units as coordinates (default 0.02)
                           Points within this distance from line are removed
        center_lat: Optional latitude for scaling lon component (for lat/lon coordinates)
    
    Returns:
        Simplified coordinate list
    """
    if len(coords) < 3:
        return coords
    
    # Handle closed polygon
    closed = False
    if len(coords) > 1 and coords[0] == coords[-1]:
        closed = True
        coords = coords[:-1]
    
    if len(coords) < 3:
        return coords + ([coords[0]] if closed else [])
    
    # Get center latitude if not provided (for lat/lon coordinates)
    if center_lat is None and coords:
        center_lat = sum(c[1] for c in coords) / len(coords)
    
    simplified = []
    n = len(coords)
    removed_count = 0
    
    for i in range(n):
        prev = tuple(coords[i - 1])
        curr = tuple(coords[i])
        nxt = tuple(coords[(i + 1) % n])
        
        # Calculate distance from current point to line between prev and next
        dist = point_line_distance(curr, prev, nxt, center_lat)
        
        # Keep point if distance is greater than threshold
        if dist >= distance_threshold:
            simplified.append(coords[i])
        else:
            removed_count += 1
            logger.debug(f"Removing point {i}: distance={dist:.8f}, threshold={distance_threshold:.8f}")
    
    if removed_count > 0:
        logger.debug(f"Distance simplification: removed {removed_count} points")
    
    if closed and simplified:
        # Ensure polygon is closed
        if simplified[0] != simplified[-1]:
            simplified.append(simplified[0])
    
    return simplified


def simplify_property_line(
    coords: List[List[float]], 
    angle_threshold: float = 179.0,
    distance_threshold: float = 0.02,
    center_lat: Optional[float] = None
) -> List[List[float]]:
    """
    Simplify property line polygon using angle-first, then distance algorithm
    
    Process:
    1. First simplify by angle (remove collinear points)
    2. Then simplify by distance (remove points close to line segments)
    
    Args:
        coords: List of [x, y] or [lon, lat] coordinate pairs
        angle_threshold: Angle threshold in degrees (default 179.0)
        distance_threshold: Distance threshold in same units as coordinates (default 0.02)
        center_lat: Optional latitude for scaling lon component (for lat/lon coordinates)
    
    Returns:
        Simplified coordinate list
    """
    if len(coords) < 3:
        return coords
    
    original_count = len(coords)
    
    # Get center latitude if not provided (for lat/lon coordinates)
    if center_lat is None and coords:
        center_lat = sum(c[1] for c in coords) / len(coords)
    
    # Step 1: Simplify by angle first
    simplified = simplify_polygon_by_angle(coords, angle_threshold)
    angle_reduced = len(coords) - len(simplified)
    
    # Step 2: Then simplify by distance
    simplified = simplify_polygon_by_distance(simplified, distance_threshold, center_lat)
    distance_reduced = (len(coords) - angle_reduced) - len(simplified)
    
    final_count = len(simplified)
    total_reduced = original_count - final_count
    
    if total_reduced > 0:
        logger.info(f"Property line simplified: {original_count} → {final_count} points "
                   f"(angle: -{angle_reduced}, distance: -{distance_reduced}, total: -{total_reduced})")
    else:
        logger.info(f"Property line simplification: no points removed (original: {original_count} points)")
    
    return simplified
