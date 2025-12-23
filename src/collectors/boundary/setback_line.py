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
        self.front_rear_setback_m = 4.0
        self.side_setback_m = 1.0
    
    def calculate_setback_line(
        self,
        property_line: PropertyLine
    ) -> Optional[SetbackLine]:
        """
        Calculate setback line from property line
        
        Setback rules:
        - 4m from front and rear lines
        - 1m from side lines
        
        Note: Current implementation uses 4m uniform buffer (conservative approach).
        This ensures front/rear meet requirements; sides will be 4m instead of 1m.
        For exact variable setbacks per edge, more complex geometric operations are needed.
        
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
            
            # Calculate setback polygon
            setback_coords = self._calculate_setback_polygon(
                prop_poly,
                property_line
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
                area_sqm=area,
                perimeter_m=perimeter,
                setback_type="full",
                metadata={
                    "front_rear_setback_m": self.front_rear_setback_m,
                    "side_setback_m": self.side_setback_m
                }
            )
            
        except Exception as e:
            logger.warning(f"Failed to calculate setback line: {e}")
            return None
    
    def _calculate_setback_polygon(
        self,
        prop_poly: Polygon,
        property_line: PropertyLine
    ) -> Optional[List[List[float]]]:
        """
        Calculate setback polygon by offsetting edges
        
        Strategy:
        1. Buffer inward by front/rear setback (4m) - this is the maximum
        2. This ensures front/rear have correct setback
        3. Sides will have more than 1m, but that's acceptable (conservative)
        
        Alternative: Use average setback (2.5m) as approximation
        """
        try:
            # Convert meters to degrees
            center_lat = prop_poly.centroid.y
            m_per_deg_lat = 111000
            m_per_deg_lon = 111000 * abs(math.cos(math.radians(center_lat)))
            # Use average for approximation
            m_per_deg = (m_per_deg_lat + m_per_deg_lon) / 2
            
            # Use front/rear setback (4m) as the buffer distance
            # This ensures front/rear are correct, sides will be conservative
            setback_deg = self.front_rear_setback_m / m_per_deg
            
            # Ensure polygon is valid and has correct orientation
            if not prop_poly.is_valid:
                # Try to fix invalid polygon
                prop_poly = prop_poly.buffer(0)
            
            # Buffer inward (negative value shrinks the polygon)
            # Use resolution parameter for better accuracy
            buffered = prop_poly.buffer(-setback_deg, cap_style=2, join_style=2, resolution=16)
            
            # Handle case where buffer creates empty or invalid geometry
            if not isinstance(buffered, Polygon):
                logger.warning(f"Buffer returned non-polygon: {type(buffered)}")
                # Try with smaller buffer
                setback_deg = self.side_setback_m / m_per_deg
                buffered = prop_poly.buffer(-setback_deg, cap_style=2, join_style=2, resolution=16)
            
            if not isinstance(buffered, Polygon) or not buffered.is_valid:
                logger.warning("Setback buffer failed - returning None")
                return None
            
            # Verify the buffered polygon is actually smaller (inside property)
            if buffered.area >= prop_poly.area:
                logger.warning(f"Setback buffer failed - buffered area ({buffered.area}) >= property area ({prop_poly.area})")
                # This shouldn't happen with negative buffer, but if it does, return None
                return None
            
            # Extract coordinates and ensure they're in correct format
            coords = list(buffered.exterior.coords)
            # Remove duplicate closing point if present
            if len(coords) > 1 and coords[0] == coords[-1]:
                coords = coords[:-1]
            
            return [[c[0], c[1]] for c in coords]
            
        except Exception as e:
            logger.warning(f"Error calculating setback polygon: {e}")
            return None
    
