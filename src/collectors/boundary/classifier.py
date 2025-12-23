"""
Property line classifier

Classifies property line into front, rear, left side, and right side segments
"""

import math
from typing import List, Dict, Any, Optional, Tuple
from loguru import logger

try:
    from shapely.geometry import Polygon, Point, LineString
    from shapely.ops import nearest_points
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - property line classification will be limited")

from .models import PropertyLine, PropertyLineSegment, SegmentType
from .utils import calculate_line_length, haversine_distance


class PropertyLineClassifier:
    """Classifies property line into front, rear, and side segments"""
    
    # Colors for visualization
    FRONT_COLOR = "#FF0000"  # Red
    REAR_COLOR = "#0000FF"   # Blue
    LEFT_SIDE_COLOR = "#00FF00"  # Green
    RIGHT_SIDE_COLOR = "#FFFF00"  # Yellow
    
    def classify(
        self,
        property_line: PropertyLine,
        osm_roads: List[Dict[str, Any]],
        osm_buildings: List[Dict[str, Any]]
    ) -> PropertyLine:
        """
        Classify property line into front, rear, left side, right side
        
        Args:
            property_line: PropertyLine object to classify
            osm_roads: List of OSM road features
            osm_buildings: List of OSM building features
            
        Returns:
            PropertyLine with classified segments
        """
        if not SHAPELY_AVAILABLE:
            logger.warning("Shapely not available - cannot classify property line")
            return property_line
        
        try:
            # Create Shapely polygon from property line
            prop_poly = Polygon(property_line.coordinates)
            
            # Find front line (adjacent to road)
            front_segment = self._find_front_segment(
                property_line.coordinates,
                prop_poly,
                osm_roads
            )
            
            # Find rear line (opposite of front)
            rear_segment = self._find_rear_segment(
                property_line.coordinates,
                prop_poly,
                front_segment
            )
            
            # Find side lines (connect front and rear, adjacent to neighbor buildings)
            left_side, right_side = self._find_side_segments(
                property_line.coordinates,
                prop_poly,
                front_segment,
                rear_segment,
                osm_buildings
            )
            
            # Assign segments
            property_line.front = front_segment
            property_line.rear = rear_segment
            property_line.left_side = left_side
            property_line.right_side = right_side
            
            return property_line
            
        except Exception as e:
            logger.warning(f"Failed to classify property line: {e}")
            return property_line
    
    def _find_front_segment(
        self,
        coords: List[List[float]],
        prop_poly: Polygon,
        osm_roads: List[Dict[str, Any]]
    ) -> Optional[PropertyLineSegment]:
        """Find front segment (adjacent to road)"""
        if not osm_roads:
            return None
        
        # Remove duplicate closing point
        coords_clean = coords[:-1] if coords[0] == coords[-1] and len(coords) > 1 else coords
        
        # Find the edge segment closest to roads
        min_distance = float('inf')
        best_segment = None
        best_segment_coords = None
        
        # Check each edge of the polygon
        for i in range(len(coords_clean)):
            j = (i + 1) % len(coords_clean)
            edge_coords = [coords_clean[i], coords_clean[j]]
            edge_line = LineString(edge_coords)
            
            # Find minimum distance to any road
            for road in osm_roads:
                centerline = road.get("centerline", {})
                if not centerline:
                    continue
                road_coords = centerline.get("coordinates", [])
                if len(road_coords) < 2:
                    continue
                
                road_line = LineString(road_coords)
                distance = edge_line.distance(road_line) * 111000  # Convert to meters
                
                # Consider edges within 10m of roads as potential front edges
                if distance < 10.0 and distance < min_distance:
                    min_distance = distance
                    best_segment_coords = edge_coords
                    best_segment = (i, j)
        
        if best_segment_coords:
            # Extend segment to include more points if they're also close to roads
            # For now, return the edge segment
            length = calculate_line_length(best_segment_coords)
            return PropertyLineSegment(
                segment_type=SegmentType.FRONT,
                coordinates=best_segment_coords,
                length_m=length,
                color=self.FRONT_COLOR,
                metadata={"distance_to_road_m": min_distance}
            )
        
        return None
    
    def _find_rear_segment(
        self,
        coords: List[List[float]],
        prop_poly: Polygon,
        front_segment: Optional[PropertyLineSegment]
    ) -> Optional[PropertyLineSegment]:
        """Find rear segment (opposite of front)"""
        if not front_segment:
            return None
        
        # Remove duplicate closing point
        coords_clean = coords[:-1] if coords[0] == coords[-1] and len(coords) > 1 else coords
        
        # Find the edge segment farthest from front
        if len(coords_clean) < 4:
            return None
        
        # Get centroid of property
        centroid = prop_poly.centroid
        
        # Find front edge midpoint
        front_midpoint = Point(
            (front_segment.coordinates[0][0] + front_segment.coordinates[1][0]) / 2,
            (front_segment.coordinates[0][1] + front_segment.coordinates[1][1]) / 2
        )
        
        # Find edge farthest from front midpoint
        max_distance = 0
        best_segment_coords = None
        
        for i in range(len(coords_clean)):
            j = (i + 1) % len(coords_clean)
            edge_coords = [coords_clean[i], coords_clean[j]]
            edge_midpoint = Point(
                (edge_coords[0][0] + edge_coords[1][0]) / 2,
                (edge_coords[0][1] + edge_coords[1][1]) / 2
            )
            
            distance = front_midpoint.distance(edge_midpoint) * 111000
            
            if distance > max_distance:
                max_distance = distance
                best_segment_coords = edge_coords
        
        if best_segment_coords:
            length = calculate_line_length(best_segment_coords)
            return PropertyLineSegment(
                segment_type=SegmentType.REAR,
                coordinates=best_segment_coords,
                length_m=length,
                color=self.REAR_COLOR,
                metadata={"distance_from_front_m": max_distance}
            )
        
        return None
    
    def _find_side_segments(
        self,
        coords: List[List[float]],
        prop_poly: Polygon,
        front_segment: Optional[PropertyLineSegment],
        rear_segment: Optional[PropertyLineSegment],
        osm_buildings: List[Dict[str, Any]]
    ) -> Tuple[Optional[PropertyLineSegment], Optional[PropertyLineSegment]]:
        """Find left and right side segments"""
        if not front_segment or not rear_segment:
            return None, None
        
        # Remove duplicate closing point
        coords_clean = coords[:-1] if coords[0] == coords[-1] and len(coords) > 1 else coords
        
        # Find edges that connect front to rear
        # Get front and rear edge indices
        front_start = None
        front_end = None
        rear_start = None
        rear_end = None
        
        for i in range(len(coords_clean)):
            j = (i + 1) % len(coords_clean)
            
            # Check if this edge matches front segment
            if self._edges_match(coords_clean[i], coords_clean[j], front_segment.coordinates):
                front_start = i
                front_end = j
            
            # Check if this edge matches rear segment
            if self._edges_match(coords_clean[i], coords_clean[j], rear_segment.coordinates):
                rear_start = i
                rear_end = j
        
        if front_start is None or rear_start is None:
            return None, None
        
        # Find paths from front to rear (two sides)
        # Left side: counter-clockwise from front_end to rear_start
        # Right side: clockwise from front_start to rear_end
        
        left_side_coords = self._get_path_coords(coords_clean, front_end, rear_start)
        right_side_coords = self._get_path_coords(coords_clean, rear_end, front_start)
        
        left_side = None
        right_side = None
        
        if left_side_coords and len(left_side_coords) >= 2:
            length = calculate_line_length(left_side_coords)
            left_side = PropertyLineSegment(
                segment_type=SegmentType.LEFT_SIDE,
                coordinates=left_side_coords,
                length_m=length,
                color=self.LEFT_SIDE_COLOR,
                metadata={}
            )
        
        if right_side_coords and len(right_side_coords) >= 2:
            length = calculate_line_length(right_side_coords)
            right_side = PropertyLineSegment(
                segment_type=SegmentType.RIGHT_SIDE,
                coordinates=right_side_coords,
                length_m=length,
                color=self.RIGHT_SIDE_COLOR,
                metadata={}
            )
        
        return left_side, right_side
    
    def _edges_match(
        self,
        p1: List[float],
        p2: List[float],
        edge: List[List[float]]
    ) -> bool:
        """Check if two points form the same edge"""
        if len(edge) < 2:
            return False
        
        # Check both directions
        match1 = (abs(p1[0] - edge[0][0]) < 1e-6 and abs(p1[1] - edge[0][1]) < 1e-6 and
                  abs(p2[0] - edge[1][0]) < 1e-6 and abs(p2[1] - edge[1][1]) < 1e-6)
        match2 = (abs(p1[0] - edge[1][0]) < 1e-6 and abs(p1[1] - edge[1][1]) < 1e-6 and
                  abs(p2[0] - edge[0][0]) < 1e-6 and abs(p2[1] - edge[0][1]) < 1e-6)
        
        return match1 or match2
    
    def _get_path_coords(
        self,
        coords: List[List[float]],
        start_idx: int,
        end_idx: int
    ) -> List[List[float]]:
        """Get coordinates along path from start_idx to end_idx"""
        if start_idx == end_idx:
            return [coords[start_idx]]
        
        path_coords = []
        current = start_idx
        
        while current != end_idx:
            path_coords.append(coords[current])
            current = (current + 1) % len(coords)
        
        # Include end point
        path_coords.append(coords[end_idx])
        
        return path_coords

