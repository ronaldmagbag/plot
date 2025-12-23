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
            
            # Validate it's actually a rectangle (should have 4 or 5 points: 4 corners + optional closing)
            coords_list = list(rectangle.exterior.coords)
            if len(coords_list) < 4 or len(coords_list) > 5:
                logger.warning(f"Invalid rectangle: expected 4-5 points, got {len(coords_list)}")
                return None
            
            # Convert to coordinates (remove closing point if present)
            rect_coords = [[c[0], c[1]] for c in coords_list]
            if len(rect_coords) > 4 and rect_coords[0] == rect_coords[-1]:
                rect_coords = rect_coords[:-1]
            
            # Ensure we have exactly 4 corners
            if len(rect_coords) != 4:
                logger.warning(f"Rectangle has {len(rect_coords)} corners, expected 4")
                return None
            
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
        
        # Sample points to try (increased for better coverage)
        num_samples = 20
        best_rect = None
        best_area = 0
        
        # Try different rectangle sizes (from small to large)
        for i in range(1, num_samples + 1):
            for j in range(1, num_samples + 1):
                # Try different widths and heights
                w = width * i / num_samples
                h = height * j / num_samples
                
                if w == 0 or h == 0:
                    continue
                
                # Try different positions (more positions for smaller rectangles)
                num_positions = max(5, num_samples // 3)  # At least 5 positions
                for k in range(num_positions):
                    for l in range(num_positions):
                        # Calculate position ensuring rectangle stays within bounds
                        max_x_pos = max_x - w
                        max_y_pos = max_y - h
                        if max_x_pos < min_x or max_y_pos < min_y:
                            continue
                            
                        x = min_x + (max_x_pos - min_x) * k / max(num_positions - 1, 1)
                        y = min_y + (max_y_pos - min_y) * l / max(num_positions - 1, 1)
                        
                        # Create rectangle
                        rect = box(x, y, x + w, y + h)
                        
                        # Only accept rectangles that are fully contained (not just intersected)
                        # This ensures we get actual rectangles, not clipped polygons
                        if polygon.contains(rect):
                            if rect.area > best_area:
                                best_area = rect.area
                                best_rect = rect
        
        # If no rectangle found, try a simpler approach with progressively smaller rectangles
        if not best_rect:
            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2
            
            # Try progressively smaller rectangles (from 90% down to 10% of bounding box)
            for scale in [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
                w = width * scale
                h = height * scale
                
                if w == 0 or h == 0:
                    continue
                
                # Try centered rectangle
                test_rect = box(center_x - w/2, center_y - h/2,
                               center_x + w/2, center_y + h/2)
                
                if polygon.contains(test_rect):
                    logger.info(f"Found fallback rectangle at {scale*100:.0f}% scale")
                    return test_rect
            
            # Last resort: try bounding box itself
            bbox = box(min_x, min_y, max_x, max_y)
            if polygon.contains(bbox):
                logger.info("Using bounding box as buildable envelope")
                return bbox
            
            # If still nothing, log warning and return None
            logger.warning("Could not find any rectangle inside setback polygon")
            return None
        
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
