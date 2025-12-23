"""
Geometry utility functions

Common geometry calculations for boundary processing
"""

import math
from typing import List, Tuple


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

