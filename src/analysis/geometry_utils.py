"""
Geometry utilities for coordinate transformations and calculations
"""

import math
from typing import List, Tuple, Dict, Any, Optional
from dataclasses import dataclass


@dataclass
class LocalCoord:
    """Local coordinate in meters from reference point"""
    x: float
    y: float


class GeometryUtils:
    """Utility functions for geometric operations"""
    
    @staticmethod
    def degrees_to_local(
        coords: List[List[float]],
        ref_lon: float,
        ref_lat: float
    ) -> List[Tuple[float, float]]:
        """
        Convert [lon, lat] coordinates to local [x, y] meters from reference
        """
        m_per_deg_lat = 111000
        m_per_deg_lon = 111000 * math.cos(math.radians(ref_lat))
        
        local_coords = []
        for lon, lat in coords:
            x = (lon - ref_lon) * m_per_deg_lon
            y = (lat - ref_lat) * m_per_deg_lat
            local_coords.append((x, y))
        
        return local_coords
    
    @staticmethod
    def local_to_degrees(
        local_coords: List[Tuple[float, float]],
        ref_lon: float,
        ref_lat: float
    ) -> List[List[float]]:
        """
        Convert local [x, y] meters to [lon, lat] degrees
        """
        m_per_deg_lat = 111000
        m_per_deg_lon = 111000 * math.cos(math.radians(ref_lat))
        
        deg_coords = []
        for x, y in local_coords:
            lon = ref_lon + x / m_per_deg_lon
            lat = ref_lat + y / m_per_deg_lat
            deg_coords.append([lon, lat])
        
        return deg_coords
    
    @staticmethod
    def polygon_centroid(coords: List[List[float]]) -> Tuple[float, float]:
        """Calculate centroid of polygon"""
        if not coords:
            return (0, 0)
        
        n = len(coords)
        if coords[0] == coords[-1]:
            n -= 1
        
        sum_x = sum(c[0] for c in coords[:n])
        sum_y = sum(c[1] for c in coords[:n])
        
        return (sum_x / n, sum_y / n)
    
    @staticmethod
    def polygon_area_local(coords: List[Tuple[float, float]]) -> float:
        """Calculate polygon area in square meters (local coords)"""
        if len(coords) < 3:
            return 0.0
        
        n = len(coords)
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += coords[i][0] * coords[j][1]
            area -= coords[j][0] * coords[i][1]
        
        return abs(area) / 2.0
    
    @staticmethod
    def polygon_perimeter_local(coords: List[Tuple[float, float]]) -> float:
        """Calculate polygon perimeter in meters (local coords)"""
        if len(coords) < 2:
            return 0.0
        
        perimeter = 0.0
        n = len(coords)
        for i in range(n):
            j = (i + 1) % n
            dx = coords[j][0] - coords[i][0]
            dy = coords[j][1] - coords[i][1]
            perimeter += math.sqrt(dx*dx + dy*dy)
        
        return perimeter
    
    @staticmethod
    def buffer_polygon(
        coords: List[Tuple[float, float]],
        buffer_distance: float
    ) -> List[Tuple[float, float]]:
        """
        Simple polygon buffer (shrink if negative, expand if positive)
        
        Uses edge offset method - offset each edge then find intersections
        For negative buffer (shrink): offset edges INWARD
        For positive buffer (expand): offset edges OUTWARD
        """
        if len(coords) < 3:
            return coords
        
        # Ensure polygon is closed
        if coords[0] != coords[-1]:
            coords = list(coords) + [coords[0]]
        
        n = len(coords) - 1  # Don't count closing vertex
        
        # Determine polygon winding (CCW = positive area)
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += coords[i][0] * coords[j][1]
            area -= coords[j][0] * coords[i][1]
        is_ccw = area > 0
        
        # For CCW polygon: inward normal is to the RIGHT of edge direction
        # For CW polygon: inward normal is to the LEFT of edge direction
        # Negative buffer = shrink = move inward
        # Positive buffer = expand = move outward
        
        # Flip sign for CCW polygons to get correct inward direction
        effective_distance = buffer_distance if is_ccw else -buffer_distance
        
        # Calculate offset edges
        offset_edges = []
        for i in range(n):
            j = (i + 1) % n
            
            # Edge direction
            dx = coords[j][0] - coords[i][0]
            dy = coords[j][1] - coords[i][1]
            length = math.sqrt(dx*dx + dy*dy)
            
            if length < 0.0001:
                continue
            
            # Normalize
            dx /= length
            dy /= length
            
            # Perpendicular (right-hand normal for CCW = inward)
            nx = dy
            ny = -dx
            
            # Offset this edge
            offset_x = nx * effective_distance
            offset_y = ny * effective_distance
            
            p1 = (coords[i][0] + offset_x, coords[i][1] + offset_y)
            p2 = (coords[j][0] + offset_x, coords[j][1] + offset_y)
            
            offset_edges.append((p1, p2, dx, dy))
        
        if len(offset_edges) < 3:
            return coords
        
        # Find intersections of consecutive offset edges
        buffered = []
        for i in range(len(offset_edges)):
            j = (i + 1) % len(offset_edges)
            
            # Line 1: p1 + t * d1
            p1 = offset_edges[i][0]
            d1 = (offset_edges[i][2], offset_edges[i][3])
            
            # Line 2: p2 + s * d2  
            p2 = offset_edges[j][0]
            d2 = (offset_edges[j][2], offset_edges[j][3])
            
            # Find intersection
            cross = d1[0] * d2[1] - d1[1] * d2[0]
            
            if abs(cross) < 0.0001:
                # Parallel lines, use endpoint
                buffered.append(offset_edges[i][1])
            else:
                # Calculate intersection parameter
                dx = p2[0] - p1[0]
                dy = p2[1] - p1[1]
                t = (dx * d2[1] - dy * d2[0]) / cross
                
                # Limit t to prevent extreme values at sharp corners
                t = max(-100, min(100, t))
                
                ix = p1[0] + t * d1[0]
                iy = p1[1] + t * d1[1]
                buffered.append((ix, iy))
        
        # Close polygon
        if buffered:
            buffered.append(buffered[0])
        
        return buffered
    
    @staticmethod
    def get_polygon_edges(
        coords: List[Tuple[float, float]]
    ) -> List[Dict[str, Any]]:
        """
        Get edges of polygon with their properties
        
        Returns list of edges with start, end, length, and direction
        """
        if len(coords) < 2:
            return []
        
        # Ensure closed
        if coords[0] != coords[-1]:
            coords = coords + [coords[0]]
        
        edges = []
        n = len(coords) - 1
        
        for i in range(n):
            start = coords[i]
            end = coords[i + 1]
            
            dx = end[0] - start[0]
            dy = end[1] - start[1]
            length = math.sqrt(dx*dx + dy*dy)
            
            # Calculate direction (bearing from start to end)
            if length > 0:
                angle = math.degrees(math.atan2(dx, dy))  # 0 = north, 90 = east
                if angle < 0:
                    angle += 360
            else:
                angle = 0
            
            # Determine cardinal direction
            direction = GeometryUtils._angle_to_direction(angle)
            
            edges.append({
                "index": i,
                "start": start,
                "end": end,
                "length_m": length,
                "bearing": angle,
                "direction": direction
            })
        
        return edges
    
    @staticmethod
    def _angle_to_direction(angle: float) -> str:
        """Convert bearing angle to cardinal direction"""
        if angle < 22.5 or angle >= 337.5:
            return "north"
        elif angle < 67.5:
            return "northeast"
        elif angle < 112.5:
            return "east"
        elif angle < 157.5:
            return "southeast"
        elif angle < 202.5:
            return "south"
        elif angle < 247.5:
            return "southwest"
        elif angle < 292.5:
            return "west"
        else:
            return "northwest"
    
    @staticmethod
    def distance_point_to_line(
        point: Tuple[float, float],
        line_start: Tuple[float, float],
        line_end: Tuple[float, float]
    ) -> float:
        """Calculate perpendicular distance from point to line segment"""
        px, py = point
        x1, y1 = line_start
        x2, y2 = line_end
        
        dx = x2 - x1
        dy = y2 - y1
        
        if dx == 0 and dy == 0:
            # Line is a point
            return math.sqrt((px - x1)**2 + (py - y1)**2)
        
        # Parameter t for closest point on line
        t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / (dx*dx + dy*dy)))
        
        # Closest point
        closest_x = x1 + t * dx
        closest_y = y1 + t * dy
        
        return math.sqrt((px - closest_x)**2 + (py - closest_y)**2)
    
    @staticmethod
    def distance_polygon_to_polygon(
        poly1: List[Tuple[float, float]],
        poly2: List[Tuple[float, float]]
    ) -> float:
        """Calculate minimum distance between two polygons"""
        min_dist = float('inf')
        
        # Check all vertices of poly1 against all edges of poly2
        edges2 = GeometryUtils.get_polygon_edges(poly2)
        for vertex in poly1:
            for edge in edges2:
                dist = GeometryUtils.distance_point_to_line(
                    vertex, edge["start"], edge["end"]
                )
                min_dist = min(min_dist, dist)
        
        # Check all vertices of poly2 against all edges of poly1
        edges1 = GeometryUtils.get_polygon_edges(poly1)
        for vertex in poly2:
            for edge in edges1:
                dist = GeometryUtils.distance_point_to_line(
                    vertex, edge["start"], edge["end"]
                )
                min_dist = min(min_dist, dist)
        
        return min_dist

