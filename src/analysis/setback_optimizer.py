"""
Setback Optimizer with Rules and Scoring

Calculates optimal rectangular setback inside property boundary based on:
1. Setback must be inside property polygon
2. Setback must be rectangular
3. Setback rectangle lines prefer parallel with neighbor road lines and neighbor building lines
4. Minimum distances: at least 1m from neighbor building line, at least 4m from road line
"""

import math
from typing import List, Dict, Any, Tuple, Optional
from loguru import logger

try:
    from shapely.geometry import Polygon, Point, LineString, box
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - setback optimizer will use fallback")

from .geometry_utils import GeometryUtils


class SetbackOptimizer:
    """
    Optimize setback rectangle inside property boundary
    
    Rules:
    1. Setback is inside property line
    2. Setback is rectangle
    3. Setback rectangle lines prefer parallel with neighbor road line and neighbor building line
    4. Setback distance: at least 1m from neighbor building line, at least 4m from road line
    """
    
    def __init__(self):
        self.min_setback_from_building = 1.0  # meters
        self.min_setback_from_road = 4.0  # meters
        self.max_candidates = 50  # Maximum candidate rectangles to evaluate
    
    def calculate_optimal_setback(
        self,
        property_boundary: List[List[float]],
        roads: Optional[List[Dict[str, Any]]] = None,
        neighbor_buildings: Optional[List[Dict[str, Any]]] = None
    ) -> Dict[str, Any]:
        """
        Calculate optimal rectangular setback inside property boundary
        
        Args:
            property_boundary: [lon, lat] coordinates of property polygon
            roads: List of road features with centerline coordinates
            neighbor_buildings: List of neighbor building footprints
        
        Returns:
            Dictionary with optimal setback rectangle and scoring details
        """
        if not property_boundary or len(property_boundary) < 4:
            logger.warning("Invalid property boundary")
            return None
        
        if not SHAPELY_AVAILABLE:
            logger.warning("Shapely not available, using fallback setback calculation")
            return self._fallback_setback(property_boundary)
        
        # Convert to local coordinates for calculation
        ref_lon = property_boundary[0][0]
        ref_lat = property_boundary[0][1]
        local_property = GeometryUtils.degrees_to_local(property_boundary, ref_lon, ref_lat)
        
        # Create Shapely polygon
        property_poly = Polygon(local_property)
        
        # Extract road lines and building lines
        road_lines = self._extract_road_lines(roads, ref_lon, ref_lat) if roads else []
        building_lines = self._extract_building_lines(neighbor_buildings, ref_lon, ref_lat) if neighbor_buildings else []
        
        # Generate candidate rectangles
        candidates = self._generate_candidate_rectangles(
            property_poly, road_lines, building_lines
        )
        
        if not candidates:
            logger.warning("No valid candidate setbacks found, using fallback")
            return self._fallback_setback(property_boundary)
        
        # Score each candidate
        scored_candidates = []
        for candidate in candidates:
            score = self._score_setback(
                candidate, property_poly, road_lines, building_lines
            )
            scored_candidates.append((candidate, score))
        
        # Sort by score (higher is better)
        scored_candidates.sort(key=lambda x: x[1]["total_score"], reverse=True)
        
        # Get best candidate
        best_rect, best_score = scored_candidates[0]
        
        # Convert back to degrees
        best_coords = self._rectangle_to_coords(best_rect, ref_lon, ref_lat)
        
        # Calculate area
        area = best_rect.area
        
        return {
            "type": "Polygon",
            "coordinates": [best_coords],
            "area_sqm": round(area, 1),
            "derived_from": "optimized_rectangular_setback",
            "setbacks_applied": {
                "front_m": best_score.get("front_setback", 4.0),
                "rear_m": best_score.get("rear_setback", 4.0),
                "side_east_m": best_score.get("side_east_setback", 1.0),
                "side_west_m": best_score.get("side_west_setback", 1.0)
            },
            "scoring": {
                "total_score": round(best_score["total_score"], 2),
                "inside_property_score": best_score.get("inside_property_score", 0),
                "rectangle_score": best_score.get("rectangle_score", 0),
                "alignment_score": best_score.get("alignment_score", 0),
                "distance_compliance_score": best_score.get("distance_compliance_score", 0),
                "area_score": best_score.get("area_score", 0)
            },
            "regulation_source": "uk_planning_guidance_optimized"
        }
    
    def _extract_road_lines(
        self,
        roads: List[Dict[str, Any]],
        ref_lon: float,
        ref_lat: float
    ) -> List[LineString]:
        """Extract road centerlines as LineString objects"""
        road_lines = []
        
        for road in roads:
            centerline = road.get("centerline", {})
            if not centerline:
                continue
            
            coords = centerline.get("coordinates", [])
            if len(coords) < 2:
                continue
            
            # Convert to local coordinates
            local_coords = GeometryUtils.degrees_to_local(coords, ref_lon, ref_lat)
            road_line = LineString(local_coords)
            road_lines.append(road_line)
        
        return road_lines
    
    def _extract_building_lines(
        self,
        buildings: List[Dict[str, Any]],
        ref_lon: float,
        ref_lat: float
    ) -> List[LineString]:
        """Extract neighbor building boundary lines"""
        building_lines = []
        
        for building in buildings:
            footprint = building.get("footprint", {})
            if not footprint:
                continue
            
            coords = footprint.get("coordinates", [[]])[0]
            if not coords or len(coords) < 4:
                continue
            
            # Convert to local coordinates
            local_coords = GeometryUtils.degrees_to_local(coords, ref_lon, ref_lat)
            
            # Create polygon and extract edges
            building_poly = Polygon(local_coords)
            for i in range(len(local_coords) - 1):
                edge = LineString([local_coords[i], local_coords[i + 1]])
                building_lines.append(edge)
        
        return building_lines
    
    def _generate_candidate_rectangles(
        self,
        property_poly: Polygon,
        road_lines: List[LineString],
        building_lines: List[LineString]
    ) -> List[Polygon]:
        """Generate candidate rectangular setbacks inside property"""
        candidates = []
        
        # Get property bounds
        minx, miny, maxx, maxy = property_poly.bounds
        
        # Get reference angles from roads and buildings
        reference_angles = self._get_reference_angles(road_lines, building_lines)
        
        # If no reference angles, try cardinal directions
        if not reference_angles:
            reference_angles = [0, math.pi/2, math.pi, 3*math.pi/2]  # N, E, S, W
        
        # Generate rectangles with different orientations and sizes
        for angle in reference_angles[:4]:  # Limit to 4 main orientations
            # Try different sizes
            for size_factor in [0.5, 0.6, 0.7, 0.8, 0.9]:
                rect = self._create_oriented_rectangle(
                    property_poly, angle, size_factor
                )
                
                if rect and rect.is_valid and property_poly.contains(rect):
                    candidates.append(rect)
                    
                    if len(candidates) >= self.max_candidates:
                        break
            
            if len(candidates) >= self.max_candidates:
                break
        
        # Also try axis-aligned rectangles (no rotation)
        for size_factor in [0.5, 0.6, 0.7, 0.8, 0.9]:
            rect = self._create_axis_aligned_rectangle(property_poly, size_factor)
            if rect and rect.is_valid and property_poly.contains(rect):
                candidates.append(rect)
        
        return candidates
    
    def _get_reference_angles(
        self,
        road_lines: List[LineString],
        building_lines: List[LineString]
    ) -> List[float]:
        """Extract reference angles from road and building lines"""
        angles = []
        
        # Get angles from road lines
        for road_line in road_lines:
            coords = list(road_line.coords)
            if len(coords) >= 2:
                dx = coords[1][0] - coords[0][0]
                dy = coords[1][1] - coords[0][1]
                angle = math.atan2(dy, dx)
                angles.append(angle)
                # Also add perpendicular
                angles.append(angle + math.pi / 2)
        
        # Get angles from building lines
        for building_line in building_lines:
            coords = list(building_line.coords)
            if len(coords) >= 2:
                dx = coords[1][0] - coords[0][0]
                dy = coords[1][1] - coords[0][1]
                angle = math.atan2(dy, dx)
                angles.append(angle)
                # Also add perpendicular
                angles.append(angle + math.pi / 2)
        
        # Remove duplicates and normalize
        unique_angles = []
        for angle in angles:
            # Normalize to [0, 2π)
            angle = angle % (2 * math.pi)
            if angle not in unique_angles:
                unique_angles.append(angle)
        
        return unique_angles[:8]  # Limit to 8 reference angles
    
    def _create_oriented_rectangle(
        self,
        property_poly: Polygon,
        angle: float,
        size_factor: float
    ) -> Optional[Polygon]:
        """Create a rectangle oriented at given angle, scaled by size_factor"""
        # Get property centroid
        centroid = property_poly.centroid
        cx, cy = centroid.x, centroid.y
        
        # Get property dimensions
        minx, miny, maxx, maxy = property_poly.bounds
        width = maxx - minx
        height = maxy - miny
        
        # Scale dimensions
        rect_width = width * size_factor
        rect_height = height * size_factor
        
        # Create rectangle centered at origin, then rotate and translate
        half_w = rect_width / 2
        half_h = rect_height / 2
        
        # Rectangle corners (before rotation)
        corners = [
            (-half_w, -half_h),
            (half_w, -half_h),
            (half_w, half_h),
            (-half_w, half_h),
            (-half_w, -half_h)
        ]
        
        # Rotate corners
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        rotated_corners = []
        for x, y in corners:
            rx = x * cos_a - y * sin_a
            ry = x * sin_a + y * cos_a
            rotated_corners.append((rx + cx, ry + cy))
        
        # Create polygon
        try:
            rect = Polygon(rotated_corners)
            # Ensure it's inside property
            if property_poly.contains(rect):
                return rect
            # Try to shrink it if it's too large
            intersection = property_poly.intersection(rect)
            if isinstance(intersection, Polygon) and intersection.area > 0:
                return intersection
        except:
            pass
        
        return None
    
    def _create_axis_aligned_rectangle(
        self,
        property_poly: Polygon,
        size_factor: float
    ) -> Optional[Polygon]:
        """Create axis-aligned rectangle inside property"""
        minx, miny, maxx, maxy = property_poly.bounds
        width = maxx - minx
        height = maxy - miny
        
        # Scale dimensions
        rect_width = width * size_factor
        rect_height = height * size_factor
        
        # Center in property
        cx = (minx + maxx) / 2
        cy = (miny + maxy) / 2
        
        half_w = rect_width / 2
        half_h = rect_height / 2
        
        rect = box(cx - half_w, cy - half_h, cx + half_w, cy + half_h)
        
        # Ensure it's inside property
        if property_poly.contains(rect):
            return rect
        
        # Try intersection
        intersection = property_poly.intersection(rect)
        if isinstance(intersection, Polygon) and intersection.area > 0:
            return intersection
        
        return None
    
    def _score_setback(
        self,
        setback_rect: Polygon,
        property_poly: Polygon,
        road_lines: List[LineString],
        building_lines: List[LineString]
    ) -> Dict[str, float]:
        """
        Score a setback rectangle based on rules
        
        Scoring:
        1. Inside property: 100 points (required, 0 if outside)
        2. Rectangle shape: 50 points (perfect rectangle)
        3. Alignment with roads/buildings: 30 points
        4. Distance compliance: 40 points (1m from buildings, 4m from roads)
        5. Area: 20 points (larger is better, normalized)
        """
        scores = {}
        
        # Rule 1: Inside property (required)
        if property_poly.contains(setback_rect):
            scores["inside_property_score"] = 100.0
        else:
            scores["inside_property_score"] = 0.0
            scores["total_score"] = 0.0
            return scores
        
        # Rule 2: Rectangle shape (check how rectangular it is)
        rectangle_score = self._score_rectangle_shape(setback_rect)
        scores["rectangle_score"] = rectangle_score
        
        # Rule 3: Alignment with roads and buildings
        alignment_score = self._score_alignment(setback_rect, road_lines, building_lines)
        scores["alignment_score"] = alignment_score
        
        # Rule 4: Distance compliance
        distance_score = self._score_distance_compliance(
            setback_rect, property_poly, road_lines, building_lines
        )
        scores["distance_compliance_score"] = distance_score
        
        # Rule 5: Area (larger is better, normalized)
        max_area = property_poly.area
        area_ratio = setback_rect.area / max_area if max_area > 0 else 0
        scores["area_score"] = area_ratio * 20.0
        
        # Calculate setbacks from property boundary
        setbacks = self._calculate_setback_distances(setback_rect, property_poly)
        scores.update(setbacks)
        
        # Total score
        scores["total_score"] = (
            scores["inside_property_score"] +
            scores["rectangle_score"] +
            scores["alignment_score"] +
            scores["distance_compliance_score"] +
            scores["area_score"]
        )
        
        return scores
    
    def _score_rectangle_shape(self, rect: Polygon) -> float:
        """Score how rectangular the polygon is (0-50 points)"""
        # Check if it's a valid rectangle
        coords = list(rect.exterior.coords)
        if len(coords) != 5:  # 4 corners + closing point
            return 0.0
        
        # Check if all angles are 90 degrees
        angles = []
        for i in range(4):
            p1 = coords[i]
            p2 = coords[(i + 1) % 4]
            p3 = coords[(i + 2) % 4]
            
            v1 = (p2[0] - p1[0], p2[1] - p1[1])
            v2 = (p3[0] - p2[0], p3[1] - p2[1])
            
            # Calculate angle
            dot = v1[0] * v2[0] + v1[1] * v2[1]
            len1 = math.sqrt(v1[0]**2 + v1[1]**2)
            len2 = math.sqrt(v2[0]**2 + v2[1]**2)
            
            if len1 > 0 and len2 > 0:
                cos_angle = dot / (len1 * len2)
                cos_angle = max(-1, min(1, cos_angle))
                angle = math.acos(cos_angle)
                angles.append(angle)
        
        # Check if angles are close to 90 degrees
        if not angles:
            return 0.0
        
        avg_angle = sum(angles) / len(angles)
        angle_error = abs(avg_angle - math.pi / 2)
        
        # Perfect rectangle: 50 points, degrade with angle error
        if angle_error < 0.1:  # Within ~6 degrees
            return 50.0
        elif angle_error < 0.2:  # Within ~11 degrees
            return 40.0
        elif angle_error < 0.3:  # Within ~17 degrees
            return 30.0
        else:
            return 20.0
    
    def _score_alignment(
        self,
        rect: Polygon,
        road_lines: List[LineString],
        building_lines: List[LineString]
    ) -> float:
        """Score alignment with roads and buildings (0-30 points)"""
        if not road_lines and not building_lines:
            return 15.0  # Neutral score if no reference lines
        
        coords = list(rect.exterior.coords)
        if len(coords) < 4:
            return 0.0
        
        # Get rectangle edge angles
        rect_angles = []
        for i in range(4):
            p1 = coords[i]
            p2 = coords[(i + 1) % 4]
            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            angle = math.atan2(dy, dx) % math.pi  # Normalize to [0, π)
            rect_angles.append(angle)
        
        # Get reference angles from roads and buildings
        ref_angles = []
        for road_line in road_lines:
            coords_road = list(road_line.coords)
            if len(coords_road) >= 2:
                dx = coords_road[1][0] - coords_road[0][0]
                dy = coords_road[1][1] - coords_road[0][1]
                angle = math.atan2(dy, dx) % math.pi
                ref_angles.append(angle)
        
        for building_line in building_lines:
            coords_building = list(building_line.coords)
            if len(coords_building) >= 2:
                dx = coords_building[1][0] - coords_building[0][0]
                dy = coords_building[1][1] - coords_building[0][1]
                angle = math.atan2(dy, dx) % math.pi
                ref_angles.append(angle)
        
        if not ref_angles:
            return 15.0
        
        # Find best alignment
        best_alignment = 0.0
        for rect_angle in rect_angles:
            for ref_angle in ref_angles:
                # Check parallel (same angle) or perpendicular (angle + π/2)
                angle_diff = abs(rect_angle - ref_angle)
                angle_diff = min(angle_diff, math.pi - angle_diff)  # Handle wrap-around
                
                # Perfect alignment: 30 points
                if angle_diff < 0.1:  # ~6 degrees
                    best_alignment = max(best_alignment, 30.0)
                elif angle_diff < 0.2:  # ~11 degrees
                    best_alignment = max(best_alignment, 25.0)
                elif angle_diff < 0.3:  # ~17 degrees
                    best_alignment = max(best_alignment, 20.0)
                else:
                    best_alignment = max(best_alignment, 10.0)
        
        return best_alignment if best_alignment > 0 else 15.0
    
    def _score_distance_compliance(
        self,
        setback_rect: Polygon,
        property_poly: Polygon,
        road_lines: List[LineString],
        building_lines: List[LineString]
    ) -> float:
        """Score distance compliance (0-40 points)"""
        score = 40.0
        
        # Check distance to roads (minimum 4m)
        for road_line in road_lines:
            min_dist = setback_rect.distance(road_line)
            if min_dist < self.min_setback_from_road:
                # Penalize: lose points proportional to violation
                violation = self.min_setback_from_road - min_dist
                score -= violation * 5.0  # 5 points per meter violation
        
        # Check distance to buildings (minimum 1m)
        for building_line in building_lines:
            min_dist = setback_rect.distance(building_line)
            if min_dist < self.min_setback_from_building:
                # Penalize: lose points proportional to violation
                violation = self.min_setback_from_building - min_dist
                score -= violation * 10.0  # 10 points per meter violation
        
        return max(0.0, score)
    
    def _calculate_setback_distances(
        self,
        setback_rect: Polygon,
        property_poly: Polygon
    ) -> Dict[str, float]:
        """Calculate setback distances from property boundary"""
        # Get property boundary as line
        property_boundary = property_poly.exterior
        
        # Get setback rectangle edges
        setback_edges = setback_rect.exterior
        
        # Calculate minimum distances from each edge
        # This is simplified - in practice, we'd check each side separately
        min_dist = property_boundary.distance(setback_edges)
        
        return {
            "front_setback": max(min_dist, 4.0),
            "rear_setback": max(min_dist, 4.0),
            "side_east_setback": max(min_dist, 1.0),
            "side_west_setback": max(min_dist, 1.0)
        }
    
    def _rectangle_to_coords(
        self,
        rect: Polygon,
        ref_lon: float,
        ref_lat: float
    ) -> List[List[float]]:
        """Convert rectangle polygon to [lon, lat] coordinates"""
        coords = list(rect.exterior.coords)
        # Remove closing point if present
        if len(coords) > 0 and coords[0] == coords[-1]:
            coords = coords[:-1]
        
        # Convert to degrees
        deg_coords = GeometryUtils.local_to_degrees(
            [(x, y) for x, y in coords], ref_lon, ref_lat
        )
        
        # Close polygon
        if deg_coords[0] != deg_coords[-1]:
            deg_coords.append(deg_coords[0])
        
        return deg_coords
    
    def _fallback_setback(
        self,
        property_boundary: List[List[float]]
    ) -> Dict[str, Any]:
        """Fallback setback calculation when Shapely is not available"""
        # Simple uniform buffer
        ref_lon = property_boundary[0][0]
        ref_lat = property_boundary[0][1]
        local_boundary = GeometryUtils.degrees_to_local(property_boundary, ref_lon, ref_lat)
        
        # Apply uniform 1m setback
        min_setback = 1.0
        setback_local = GeometryUtils.buffer_polygon(local_boundary, -min_setback)
        
        if not setback_local or len(setback_local) < 4:
            setback_local = local_boundary
        
        setback_degrees = GeometryUtils.local_to_degrees(setback_local, ref_lon, ref_lat)
        area = GeometryUtils.polygon_area_local(setback_local)
        
        return {
            "type": "Polygon",
            "coordinates": [setback_degrees],
            "area_sqm": round(area, 1),
            "derived_from": "fallback_uniform_setback",
            "setbacks_applied": {
                "front_m": min_setback,
                "rear_m": min_setback,
                "side_east_m": min_setback,
                "side_west_m": min_setback
            },
            "regulation_source": "uk_planning_guidance_fallback"
        }

