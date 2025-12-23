"""
Buildable envelope calculation

Calculates maximum area rectangle inscribed in setback polygon
"""

import math
from typing import List, Dict, Any, Optional, Tuple
from loguru import logger

try:
    from shapely.geometry import Polygon, Point, box
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - buildable envelope calculation will be limited")

from .models import BuildableEnvelope, PropertyLine, SetbackLine
from .utils import calculate_polygon_area, calculate_perimeter, close_polygon
from ...config import get_config


class BuildableEnvelopeProcessor:
    """Processes buildable envelope calculation"""
    
    def __init__(self):
        self.config = get_config()
    
    def calculate_buildable_envelope(
        self,
        property_line: PropertyLine,
        setback_line: Optional[SetbackLine] = None
    ) -> Optional[BuildableEnvelope]:
        """
        Calculate buildable envelope from setback line
        
        Finds the maximum area rectangle inscribed in the setback polygon
        
        Args:
            property_line: PropertyLine object
            setback_line: SetbackLine object (if None, will calculate from property_line)
            
        Returns:
            BuildableEnvelope object or None
        """
        if not SHAPELY_AVAILABLE:
            logger.warning("Shapely not available - cannot calculate buildable envelope")
            return None
        
        # Get setback polygon
        if setback_line:
            setback_coords = setback_line.coordinates
        else:
            # Calculate setback line first
            from .setback_line import SetbackLineProcessor
            setback_processor = SetbackLineProcessor()
            setback = setback_processor.calculate_setback_line(property_line)
            if not setback:
                logger.warning("Failed to calculate setback line")
                return None
            setback_coords = setback.coordinates
        
        if not setback_coords or len(setback_coords) < 4:
            logger.warning("Invalid setback coordinates")
            return None
        
        try:
            # Create setback polygon
            setback_coords_clean = close_polygon(setback_coords)
            setback_poly = Polygon(setback_coords_clean)
            
            if not setback_poly.is_valid:
                logger.warning("Invalid setback polygon")
                return None
            
            # Find largest inscribed rectangle
            rectangle = self._find_largest_inscribed_rectangle(setback_poly)
            
            if not rectangle:
                logger.warning("Failed to find inscribed rectangle")
                return None
            
            # Convert to coordinates
            rect_coords = [[c[0], c[1]] for c in rectangle.exterior.coords]
            rect_coords = close_polygon(rect_coords)
            
            # Calculate area and perimeter
            center_lat = sum(c[1] for c in rect_coords) / len(rect_coords)
            area = calculate_polygon_area(rect_coords, center_lat)
            perimeter = calculate_perimeter(rect_coords, center_lat)
            
            return BuildableEnvelope(
                coordinates=rect_coords,
                area_sqm=area,
                perimeter_m=perimeter,
                metadata={
                    "type": "maximum_inscribed_rectangle",
                    "setback_area_sqm": setback_line.area_sqm if setback_line else None
                }
            )
            
        except Exception as e:
            logger.warning(f"Failed to calculate buildable envelope: {e}")
            return None
    
    def _find_largest_inscribed_rectangle(
        self,
        polygon: Polygon
    ) -> Optional[Polygon]:
        """
        Find the largest rectangle inscribed in a polygon
        
        Uses a rotation-based approach:
        1. Try different rotation angles
        2. For each angle, find the largest axis-aligned rectangle
        3. Return the maximum area rectangle found
        """
        if not polygon.is_valid or polygon.area == 0:
            return None
        
        best_rectangle = None
        best_area = 0
        
        # Try different rotation angles (0 to 90 degrees in 10-degree steps for speed)
        for angle_deg in range(0, 91, 10):
            try:
                # Rotate polygon to align with axes
                rotated_poly = self._rotate_polygon(polygon, -angle_deg)
                
                # Find largest axis-aligned rectangle
                rect = self._find_largest_axis_aligned_rectangle(rotated_poly)
                
                if rect and rect.area > best_area:
                    best_area = rect.area
                    # Rotate back
                    best_rectangle = self._rotate_polygon(rect, angle_deg)
            except Exception as e:
                logger.debug(f"Error at angle {angle_deg}: {e}")
                continue
        
        return best_rectangle
    
    def _find_largest_axis_aligned_rectangle(
        self,
        polygon: Polygon
    ) -> Optional[Polygon]:
        """
        Find largest axis-aligned rectangle inside polygon
        
        Uses a sweep-line approach:
        1. Get bounding box
        2. Try different widths and heights
        3. Check if rectangle fits inside polygon
        """
        if not polygon.is_valid:
            return None
        
        bounds = polygon.bounds
        min_x, min_y, max_x, max_y = bounds
        
        # Get resolution (try different rectangle sizes)
        width = max_x - min_x
        height = max_y - min_y
        
        # Sample points to try (reduced for performance)
        num_samples = 15
        best_rect = None
        best_area = 0
        
        # Try different rectangle positions and sizes
        for i in range(num_samples):
            for j in range(num_samples):
                # Try different widths and heights
                w = width * (i + 1) / num_samples
                h = height * (j + 1) / num_samples
                
                if w == 0 or h == 0:
                    continue
                
                # Try different positions
                for k in range(num_samples // 2):
                    for l in range(num_samples // 2):
                        x = min_x + (max_x - min_x - w) * k / (num_samples // 2)
                        y = min_y + (max_y - min_y - h) * l / (num_samples // 2)
                        
                        # Create rectangle
                        rect = box(x, y, x + w, y + h)
                        
                        # Check if rectangle is inside polygon
                        if polygon.contains(rect) or polygon.intersects(rect):
                            # Get intersection (clipped rectangle)
                            intersection = polygon.intersection(rect)
                            if isinstance(intersection, Polygon) and intersection.area > best_area:
                                best_area = intersection.area
                                best_rect = intersection
        
        # If no rectangle found, try a simpler approach: use bounding box
        if not best_rect:
            # Try the bounding box itself
            bbox = box(min_x, min_y, max_x, max_y)
            if polygon.contains(bbox):
                return bbox
            
            # Try a smaller rectangle (half size, centered)
            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2
            half_w = width / 2
            half_h = height / 2
            small_rect = box(center_x - half_w, center_y - half_h,
                           center_x + half_w, center_y + half_h)
            
            if polygon.contains(small_rect):
                return small_rect
        
        return best_rect
    
    def _rotate_polygon(
        self,
        polygon: Polygon,
        angle_deg: float
    ) -> Polygon:
        """Rotate polygon around its centroid"""
        if angle_deg == 0:
            return polygon
        
        centroid = polygon.centroid
        angle_rad = math.radians(angle_deg)
        cos_a = math.cos(angle_rad)
        sin_a = math.sin(angle_rad)
        
        def rotate_point(x, y, cx, cy):
            # Translate to origin
            dx = x - cx
            dy = y - cy
            # Rotate
            new_x = dx * cos_a - dy * sin_a
            new_y = dx * sin_a + dy * cos_a
            # Translate back
            return (new_x + cx, new_y + cy)
        
        # Rotate exterior coordinates
        rotated_coords = [
            rotate_point(x, y, centroid.x, centroid.y)
            for x, y in polygon.exterior.coords[:-1]
        ]
        
        return Polygon(rotated_coords)
