"""
Buildable envelope calculation

Calculates maximum area rectangle inscribed in setback polygon
"""

import math
from typing import List, Dict, Any, Optional, Tuple
from loguru import logger

try:
    from shapely.geometry import Polygon, Point, box
    from shapely.ops import unary_union, transform
    from shapely.affinity import rotate
    from pyproj import Transformer
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely/pyproj not available - buildable envelope calculation will be limited")

from .models import BuildableEnvelope, PropertyLine, SetbackLine
from .utils import calculate_polygon_area, calculate_perimeter, close_polygon
from ...config import get_config


class BuildableEnvelopeProcessor:
    """Processes buildable envelope calculation"""
    
    def __init__(self):
        self.config = get_config()
        # Set up forward/backward transformers for metric calculations (default to UK EPSG:27700)
        self.crs_wgs84 = "EPSG:4326"
        self.crs_metric = getattr(self.config, "local_crs", "EPSG:27700")
        try:
            self.to_metric = Transformer.from_crs(self.crs_wgs84, self.crs_metric, always_xy=True).transform
            self.to_wgs84 = Transformer.from_crs(self.crs_metric, self.crs_wgs84, always_xy=True).transform
            self.transformers_ready = True
        except Exception:
            logger.warning("Failed to initialize CRS transformers; falling back to WGS84 computations (may distort rectangles)")
            self.transformers_ready = False
    
    def calculate_buildable_envelope(
        self,
        property_line: PropertyLine,
        setback_line: Optional[SetbackLine] = None
    ) -> Optional[BuildableEnvelope]:
        """
        Calculate buildable envelope from setback line
        
        Finds the maximum area rectangle inscribed in the setback polygon.
        For rectangle property lines (4 edges), uses setback line directly.
        
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
            setback_area = setback_line.area_sqm
        else:
            # Calculate setback line first
            from .setback_line import SetbackLineProcessor
            setback_processor = SetbackLineProcessor()
            setback = setback_processor.calculate_setback_line(property_line)
            if not setback:
                logger.warning("Failed to calculate setback line")
                return None
            setback_coords = setback.coordinates
            setback_area = setback.area_sqm
        
        if not setback_coords or len(setback_coords) < 4:
            logger.warning("Invalid setback coordinates")
            return None
        
        # Check if property line is a rectangle (4 edges after cleaning)
        # For rectangles, buildable area = setback area (no need to find inscribed rectangle)
        prop_coords = property_line.coordinates
        prop_coords_clean = prop_coords[:-1] if (prop_coords and len(prop_coords) > 1 and prop_coords[0] == prop_coords[-1]) else prop_coords
        
        if len(prop_coords_clean) == 4:
            # Check if it's a rectangle (2 long edges, 2 short edges)
            from .utils import calculate_line_length
            from ...config import get_config
            import math
            
            config = get_config()
            ratio_threshold = config.uk_regulatory.classifier_5point_rectangle_ratio
            
            # Calculate edge lengths
            edge_lengths = []
            for i in range(len(prop_coords_clean)):
                j = (i + 1) % len(prop_coords_clean)
                edge_start = prop_coords_clean[i]
                edge_end = prop_coords_clean[j]
                length = calculate_line_length([edge_start, edge_end])
                edge_lengths.append(length)
            
            # Sort by length
            sorted_lengths = sorted(edge_lengths, reverse=True)
            avg_long = (sorted_lengths[0] + sorted_lengths[1]) / 2
            avg_short = (sorted_lengths[2] + sorted_lengths[3]) / 2
            ratio = avg_long / avg_short if avg_short > 0 else 0
            
            if ratio >= ratio_threshold:
                # It's a rectangle - use setback line directly as buildable envelope
                logger.info(f"Property is rectangle (ratio={ratio:.2f}x >= {ratio_threshold}x) - using setback line as buildable envelope")
                from .models import BuildableEnvelope
                from .utils import calculate_perimeter
                
                # Calculate perimeter from setback coordinates
                avg_lat = sum(c[1] for c in setback_coords) / len(setback_coords) if setback_coords else 0
                perimeter_m = calculate_perimeter(setback_coords, avg_lat)
                
                return BuildableEnvelope(
                    coordinates=setback_coords,
                    area_sqm=setback_area,
                    perimeter_m=perimeter_m,
                    metadata={
                        "derived_from": "setback_line",
                        "type": "rectangle_direct",
                        "constraints_applied": []
                    }
                )
        
        logger.info("=== BUILDABLE ENVELOPE CALCULATION DEBUG ===")
        logger.info(f"Setback coordinates: {len(setback_coords)} points")
        logger.info(f"Transformers ready: {self.transformers_ready}")
        logger.info(f"Metric CRS: {self.crs_metric}")
        
        try:
            # Create setback polygon in WGS84
            setback_coords_clean = close_polygon(setback_coords)
            setback_poly_wgs = Polygon(setback_coords_clean)
            
            logger.info(f"Setback polygon (WGS84) area: {setback_poly_wgs.area:.6f} deg²")
            logger.info(f"Setback polygon (WGS84) valid: {setback_poly_wgs.is_valid}")
            
            if self.transformers_ready:
                # Reproject to metric CRS for accurate rectangle computation
                setback_poly = transform(self.to_metric, setback_poly_wgs)
                logger.info(f"Setback polygon (metric) area: {setback_poly.area:.2f} m²")
                logger.info(f"Setback polygon (metric) valid: {setback_poly.is_valid}")
                
                # Try to fix invalid polygons (common after reprojection)
                if not setback_poly.is_valid:
                    logger.warning(f"Setback polygon invalid after reprojection, attempting to fix...")
                    logger.debug(f"  Invalid reason: {setback_poly.validity_reason if hasattr(setback_poly, 'validity_reason') else 'unknown'}")
                    
                    # Try buffer(0) trick to fix self-intersections
                    fixed_poly = setback_poly.buffer(0)
                    
                    # Handle MultiPolygon case
                    if hasattr(fixed_poly, 'geoms'):
                        # It's a MultiPolygon - take the largest polygon
                        largest = max(fixed_poly.geoms, key=lambda g: g.area if hasattr(g, 'area') else 0)
                        if isinstance(largest, Polygon):
                            setback_poly = largest
                            logger.info(f"  Fixed: extracted largest polygon from MultiPolygon (area: {setback_poly.area:.2f} m²)")
                        else:
                            logger.warning(f"  Could not extract valid polygon from MultiPolygon")
                            return None
                    else:
                        setback_poly = fixed_poly
                    
                    if not setback_poly.is_valid:
                        logger.warning(f"  Could not fix invalid polygon")
                        return None
                    else:
                        logger.info(f"  Successfully fixed polygon (area: {setback_poly.area:.2f} m²)")
            else:
                setback_poly = setback_poly_wgs
                logger.warning("Using WGS84 for rectangle calculation (may cause distortion)")
            
            # Final validation check
            if not setback_poly.is_valid:
                logger.warning(f"Invalid setback polygon (final validation failed)")
                logger.debug(f"  Polygon type: {type(setback_poly)}")
                logger.debug(f"  Polygon area: {setback_poly.area}")
                # Try one more fix attempt
                try:
                    fixed = setback_poly.buffer(0)
                    if hasattr(fixed, 'geoms'):
                        fixed = max(fixed.geoms, key=lambda g: g.area if hasattr(g, 'area') else 0)
                    if isinstance(fixed, Polygon) and fixed.is_valid:
                        setback_poly = fixed
                        logger.info(f"  Final fix attempt succeeded (area: {setback_poly.area:.2f} m²)")
                    else:
                        logger.warning(f"  Final fix attempt failed - polygon still invalid")
                        return None
                except Exception as e:
                    logger.warning(f"  Final fix attempt raised exception: {e}")
                    return None
            
            logger.info(f"Setback polygon validated successfully (area: {setback_poly.area:.2f} m², valid: {setback_poly.is_valid})")
            
            # Find largest inscribed rectangle
            logger.info("Searching for largest inscribed rectangle...")
            rectangle = self._find_largest_inscribed_rectangle(setback_poly)
            
            if not rectangle:
                logger.warning("Failed to find inscribed rectangle")
                return None
            
            logger.info(f"Found rectangle area: {rectangle.area:.2f} m²")
            logger.info(f"Rectangle bounds: {rectangle.bounds}")
            
            # Validate it's actually a rectangle (orthogonal edges)
            is_rect = self._is_rectangle(rectangle)
            logger.info(f"Rectangle validation (orthogonal edges): {is_rect}")
            
            if not is_rect:
                logger.warning("Found shape is not a proper rectangle (edges not orthogonal)")
                # Log rectangle coordinates for debugging
                rect_coords = list(rectangle.exterior.coords)[:4]
                logger.debug(f"Rectangle coordinates: {rect_coords}")
                return None
            
            # Reproject rectangle back to WGS84 if needed
            if self.transformers_ready:
                rectangle_wgs = transform(self.to_wgs84, rectangle)
                logger.info(f"Reprojected rectangle back to WGS84")
            else:
                rectangle_wgs = rectangle
            
            # Validate it has exactly 4 corners
            coords_list = list(rectangle_wgs.exterior.coords)
            logger.info(f"Rectangle has {len(coords_list)} coordinate points")
            
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
            
            logger.info(f"Buildable envelope calculated:")
            logger.info(f"  Area: {area:.2f} m²")
            logger.info(f"  Perimeter: {perimeter:.2f} m")
            logger.info(f"  Setback area: {setback_line.area_sqm if setback_line else 'N/A'} m²")
            logger.info(f"  Efficiency: {(area / setback_line.area_sqm * 100) if setback_line else 0:.1f}%")
            logger.info("=== END BUILDABLE ENVELOPE DEBUG ===\n")
            
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
        best_angle = None
        
        logger.info(f"Searching for largest inscribed rectangle in polygon (area: {polygon.area:.2f} m²)")
        
        # Try different rotation angles (0 to 90 degrees in 1-degree steps)
        angles_tried = 0
        for angle_deg in range(0, 91, 5):
            try:
                # Rotate polygon to align with axes using Shapely's affine rotation
                # This preserves geometry integrity (parallel edges, right angles)
                rotated_poly = rotate(polygon, -angle_deg, origin='centroid', use_radians=False)
                
                # Find largest axis-aligned rectangle
                rect = self._find_largest_axis_aligned_rectangle(rotated_poly)
                
                if rect and rect.area > best_area:
                    best_area = rect.area
                    best_angle = angle_deg
                    # Rotate back using Shapely's affine rotation to preserve rectangle shape
                    best_rectangle = rotate(rect, angle_deg, origin='centroid', use_radians=False)
                angles_tried += 1
            except Exception as e:
                logger.debug(f"Error at angle {angle_deg}: {e}")
                continue
        
        if best_rectangle:
            logger.info(f"Best rectangle found at angle {best_angle}°: area = {best_area:.2f} m²")
            logger.info(f"  Rectangle efficiency: {(best_area / polygon.area * 100):.1f}% of setback area")
        else:
            logger.warning(f"No rectangle found after trying {angles_tried} angles")
        
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
        
        logger.debug(f"  Searching axis-aligned rectangles: bounding box {width:.2f}m × {height:.2f}m")
        
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
            logger.debug(f"  No rectangle found in main search, trying fallback strategy...")
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
    
    def _is_rectangle(self, poly: Polygon, tol: float = 1e-6) -> bool:
        """
        Validate that a polygon is a proper rectangle (orthogonal edges)
        
        Args:
            poly: Polygon to validate
            tol: Tolerance for angle checking
            
        Returns:
            True if polygon is a rectangle with orthogonal edges
        """
        coords = list(poly.exterior.coords)[:-1]
        if len(coords) != 4:
            return False
        
        def dot(a, b):
            """Dot product of two vectors"""
            return a[0]*b[0] + a[1]*b[1]
        
        # Check that adjacent edges are perpendicular (dot product should be ~0)
        for i in range(4):
            p0 = coords[i]
            p1 = coords[(i+1) % 4]
            p2 = coords[(i+2) % 4]
            v1 = (p1[0]-p0[0], p1[1]-p0[1])  # Edge from p0 to p1
            v2 = (p2[0]-p1[0], p2[1]-p1[1])  # Edge from p1 to p2
            # If edges are perpendicular, dot product should be 0
            if abs(dot(v1, v2)) > tol:
                return False
        
        return True
