"""
Setback and buildable envelope calculator
Applies UK building regulations to determine buildable area
"""

import math
from typing import List, Dict, Any, Tuple, Optional
from loguru import logger

try:
    from shapely.geometry import Polygon, LineString, Point
    from shapely.ops import unary_union
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False
    logger.warning("Shapely not available - edge-specific setbacks will use fallback")

from .geometry_utils import GeometryUtils
from ..config import get_config


class SetbackCalculator:
    """
    Calculate setback lines and buildable envelope from property boundary
    
    Applies UK residential building regulations for setbacks
    """
    
    def __init__(self):
        self.config = get_config()
        self.uk_reg = self.config.uk_regulatory
    
    def calculate_setbacks(
        self,
        property_boundary: List[List[float]],
        adjacency_info: List[Dict[str, Any]],
        custom_setbacks: Optional[Dict[str, float]] = None
    ) -> Dict[str, Any]:
        """
        Calculate setback line from property boundary
        
        Args:
            property_boundary: [lon, lat] coordinates of property line
            adjacency_info: List of edge adjacency information
            custom_setbacks: Optional override for default setbacks
        
        Returns:
            Setback line polygon and applied setbacks
        """
        if not property_boundary or len(property_boundary) < 4:
            logger.warning("Invalid property boundary for setback calculation")
            return None
        
        # Convert to local coordinates for calculation
        ref_lon = property_boundary[0][0]
        ref_lat = property_boundary[0][1]
        local_boundary = GeometryUtils.degrees_to_local(property_boundary, ref_lon, ref_lat)
        
        # Determine setbacks for each edge based on adjacency
        setbacks = self._determine_setbacks(adjacency_info, custom_setbacks)
        
        # Apply edge-specific setbacks (not average)
        # Front: 5m, Rear: 5m, Sides: 1m
        setback_local = self._apply_edge_specific_setbacks(
            local_boundary, 
            setbacks,
            adjacency_info
        )
        
        # Ensure valid polygon
        if not setback_local or len(setback_local) < 4:
            logger.warning("Setback calculation resulted in invalid polygon, using minimum setback")
            # Fallback: use minimum setback (1m) uniformly
            min_setback = min(setbacks.values())
            setback_local = GeometryUtils.buffer_polygon(local_boundary, -min_setback)
        
        # Convert back to degrees
        setback_degrees = GeometryUtils.local_to_degrees(setback_local, ref_lon, ref_lat)
        
        # Calculate area
        area = GeometryUtils.polygon_area_local(setback_local)
        
        return {
            "type": "Polygon",
            "coordinates": [setback_degrees],
            "area_sqm": round(area, 1),
            "derived_from": "property_line",
            "setbacks_applied": {
                "front_m": setbacks.get("front", self.uk_reg.front_setback_m),
                "rear_m": setbacks.get("rear", self.uk_reg.rear_setback_m),
                "side_east_m": setbacks.get("side_east", self.uk_reg.side_setback_m),
                "side_west_m": setbacks.get("side_west", self.uk_reg.side_setback_m)
            },
            "regulation_source": "uk_planning_guidance_pd"
        }
    
    def calculate_buildable_envelope(
        self,
        setback_polygon: List[List[float]],
        constraints: List[Dict[str, Any]] = None,
        property_area_sqm: float = 0,
        neighbor_buildings: List[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Calculate buildable envelope by removing constraints from setback area
        
        Args:
            setback_polygon: [lon, lat] coordinates of setback line
            constraints: List of constraint polygons to subtract
            property_area_sqm: Original property area for coverage calculation
            neighbor_buildings: List of neighboring buildings to avoid
        
        Returns:
            Buildable envelope polygon and applied constraints
        """
        if not setback_polygon or len(setback_polygon) < 4:
            logger.warning("Invalid setback polygon for envelope calculation")
            return None
        
        ref_lon = setback_polygon[0][0]
        ref_lat = setback_polygon[0][1]
        local_setback = GeometryUtils.degrees_to_local(setback_polygon, ref_lon, ref_lat)
        
        applied_constraints = []
        buildable_local = local_setback
        has_actual_constraints = False
        
        # Apply neighbor building constraints (only if neighbor is too close)
        if neighbor_buildings:
            for building in neighbor_buildings:
                footprint = building.get("footprint", {})
                if footprint.get("type") == "Polygon":
                    b_coords = footprint.get("coordinates", [[]])[0]
                    if b_coords:
                        local_building = GeometryUtils.degrees_to_local(b_coords, ref_lon, ref_lat)
                        min_dist = GeometryUtils.distance_polygon_to_polygon(buildable_local, local_building)
                        
                        # Only apply constraint if neighbor is actually too close (< 1m minimum)
                        if min_dist < 1.0:  # Less than 1m from neighbor building (violates minimum)
                            has_actual_constraints = True
                            applied_constraints.append({
                                "type": "neighbor_building_proximity",
                                "description": f"Building {building.get('id', 'unknown')} too close ({round(min_dist, 1)}m < 1m minimum)",
                                "area_removed_sqm": 0,  # Not calculating exact area
                                "min_distance_m": round(min_dist, 1)
                            })
        
        # Process explicit constraints (utility easements, TPOs, etc.)
        if constraints:
            for constraint in constraints:
                constraint_type = constraint.get("type", "unknown")
                constraint_poly = constraint.get("polygon")
                
                if constraint_poly:
                    has_actual_constraints = True
                    applied_constraints.append({
                        "type": constraint_type,
                        "description": constraint.get("description", ""),
                        "area_removed_sqm": constraint.get("area_sqm", 0)
                    })
        
        # Buildable envelope should be same as setback line unless constraints are applied
        # Only apply additional buffer if there are actual constraints that require shrinking
        if has_actual_constraints and len(buildable_local) >= 4:
            # Only apply construction tolerance if we have actual constraints
            inner_buildable = GeometryUtils.buffer_polygon(buildable_local, -0.5)
            if len(inner_buildable) >= 4:
                buildable_local = inner_buildable
                applied_constraints.append({
                    "type": "construction_tolerance",
                    "description": "0.5m buffer from setback line for construction tolerance",
                    "area_removed_sqm": 0
                })
        # If no constraints, buildable envelope = setback line (no additional shrinking)
        
        # Convert back to degrees
        buildable_degrees = GeometryUtils.local_to_degrees(buildable_local, ref_lon, ref_lat)
        
        # Calculate area
        area = GeometryUtils.polygon_area_local(buildable_local)
        
        # Determine if access corridor is needed
        access_corridor = self._determine_access_corridor(property_area_sqm, area)
        
        return {
            "type": "Polygon",
            "coordinates": [buildable_degrees],
            "area_sqm": round(area, 1),
            "derived_from": "setback_line",
            "constraints_applied": applied_constraints,
            "access_corridor": access_corridor
        }
    
    def get_building_constraints(
        self,
        property_area_sqm: float
    ) -> Dict[str, Any]:
        """
        Get building constraints based on UK regulations
        
        Returns max coverage, height restrictions, etc.
        """
        max_coverage_sqm = property_area_sqm * self.uk_reg.max_coverage_ratio
        
        return {
            "max_coverage_ratio": self.uk_reg.max_coverage_ratio,
            "max_coverage_sqm": round(max_coverage_sqm, 1),
            "max_height_m": self.uk_reg.max_height_m,
            "max_stories": self.uk_reg.max_stories,
            "max_stories_note": "Plus loft conversion allowed within height limit",
            "setbacks_m": {
                "front": self.uk_reg.front_setback_m,
                "rear": self.uk_reg.rear_setback_m,
                "side_east": self.uk_reg.side_setback_m,
                "side_west": self.uk_reg.side_setback_m,
                "notes": "Side setback increases to 2m if building height exceeds 8m"
            }
        }
    
    def get_zoning_info(self) -> Dict[str, Any]:
        """Get default UK residential zoning information"""
        return {
            "designation": self.uk_reg.default_zoning,
            "description": self.uk_reg.zoning_description,
            "permitted_uses": self.uk_reg.permitted_uses,
            "conditional_uses": self.uk_reg.conditional_uses,
            "source": "uk_planning_system",
            "last_verified": None  # Will be set by caller
        }
    
    def _determine_setbacks(
        self,
        adjacency_info: List[Dict[str, Any]],
        custom_setbacks: Optional[Dict[str, float]] = None
    ) -> Dict[str, float]:
        """
        Determine setbacks for each edge based on adjacency
        
        UK typical setbacks:
        - Front (street): 5m
        - Rear: 5m (or 1/2 of garden depth)
        - Sides: 1m (or 0 for terraced)
        """
        setbacks = {
            "front": self.uk_reg.front_setback_m,
            "rear": self.uk_reg.rear_setback_m,
            "side_east": self.uk_reg.side_setback_m,
            "side_west": self.uk_reg.side_setback_m
        }
        
        # Override with custom if provided
        if custom_setbacks:
            setbacks.update(custom_setbacks)
            return setbacks
        
        # Adjust based on adjacency
        for edge in adjacency_info:
            adjacent_to = edge.get("adjacent_to", "unknown")
            edge_id = edge.get("edge_id", "")
            
            # Street frontage usually has larger setback
            if adjacent_to == "street" and "north" in edge_id:
                setbacks["front"] = max(setbacks["front"], 5.0)
            elif adjacent_to == "street" and "south" in edge_id:
                setbacks["rear"] = max(setbacks["rear"], 5.0)
            
            # Side setbacks may be reduced if adjacent to similar building
            if adjacent_to == "neighbor_parcel":
                if edge.get("shared_wall_potential", False):
                    if "east" in edge_id:
                        setbacks["side_east"] = 0  # Party wall
                    elif "west" in edge_id:
                        setbacks["side_west"] = 0
        
        return setbacks
    
    def _apply_edge_specific_setbacks(
        self,
        local_boundary: List[Tuple[float, float]],
        setbacks: Dict[str, float],
        adjacency_info: List[Dict[str, Any]]
    ) -> List[Tuple[float, float]]:
        """
        Apply edge-specific setbacks to polygon
        Front/Rear: 5m, Sides: 1m (not average)
        
        Uses standard setbacks: front=5m, rear=5m, sides=1m
        """
        if len(local_boundary) < 4:
            return local_boundary
        
        # Close polygon if needed
        if local_boundary[0] != local_boundary[-1]:
            local_boundary = list(local_boundary) + [local_boundary[0]]
        
        # Get edges with their adjacency information
        edges = GeometryUtils.get_polygon_edges(local_boundary)
        
        # Create a mapping of edge index to setback distance
        edge_setbacks = []
        for i, edge in enumerate(edges):
            edge_id = edge.get("edge_id", "")
            direction = edge.get("direction", "").lower()
            
            # Default to side setback (1m)
            setback_dist = setbacks.get("side_east", self.uk_reg.side_setback_m)
            
            # Check adjacency info to determine front/rear/side
            for adj in adjacency_info:
                if adj.get("edge_id") == edge_id:
                    adjacent_to = adj.get("adjacent_to", "")
                    
                    if adjacent_to == "street":
                        # Street-facing: use front or rear (5m)
                        if "north" in direction or "front" in edge_id.lower():
                            setback_dist = setbacks.get("front", self.uk_reg.front_setback_m)
                        elif "south" in direction or "rear" in edge_id.lower():
                            setback_dist = setbacks.get("rear", self.uk_reg.rear_setback_m)
                        else:
                            # Default street-facing to front setback
                            setback_dist = setbacks.get("front", self.uk_reg.front_setback_m)
                    else:
                        # Neighbor-facing: use side setback (1m)
                        if "east" in direction:
                            setback_dist = setbacks.get("side_east", self.uk_reg.side_setback_m)
                        elif "west" in direction:
                            setback_dist = setbacks.get("side_west", self.uk_reg.side_setback_m)
            
            edge_setbacks.append(setback_dist)
        
        # Apply edge-specific offsets
        # Offset each edge inward by its specific setback distance
        offset_edges = []
        n = len(local_boundary) - 1  # Exclude closing vertex
        
        for i in range(n):
            j = (i + 1) % n
            p1 = local_boundary[i]
            p2 = local_boundary[j]
            
            # Edge vector
            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            length = math.sqrt(dx*dx + dy*dy)
            
            if length < 0.0001:
                continue
            
            # Normalize
            dx /= length
            dy /= length
            
            # Perpendicular (inward normal for CCW polygon)
            # Right-hand normal: (dy, -dx) for CCW
            nx = dy
            ny = -dx
            
            # Get setback for this edge
            setback = edge_setbacks[i] if i < len(edge_setbacks) else self.uk_reg.side_setback_m
            
            # Offset edge inward
            offset_x = nx * setback
            offset_y = ny * setback
            
            offset_p1 = (p1[0] + offset_x, p1[1] + offset_y)
            offset_p2 = (p2[0] + offset_x, p2[1] + offset_y)
            
            offset_edges.append((offset_p1, offset_p2, dx, dy))
        
        if len(offset_edges) < 3:
            # Fallback to uniform buffer with minimum setback
            min_setback = min(edge_setbacks) if edge_setbacks else self.uk_reg.side_setback_m
            return GeometryUtils.buffer_polygon(local_boundary, -min_setback)
        
        # Find intersections of consecutive offset edges
        setback_coords = []
        for i in range(len(offset_edges)):
            j = (i + 1) % len(offset_edges)
            
            # Line 1: offset_edges[i]
            p1 = offset_edges[i][0]
            d1 = (offset_edges[i][2], offset_edges[i][3])
            
            # Line 2: offset_edges[j]
            p2 = offset_edges[j][0]
            d2 = (offset_edges[j][2], offset_edges[j][3])
            
            # Find intersection
            cross = d1[0] * d2[1] - d1[1] * d2[0]
            
            if abs(cross) < 0.0001:
                # Parallel lines, use endpoint
                setback_coords.append(offset_edges[i][1])
            else:
                # Calculate intersection
                dx = p2[0] - p1[0]
                dy = p2[1] - p1[1]
                t = (dx * d2[1] - dy * d2[0]) / cross
                
                # Limit t to prevent extreme values
                t = max(-100, min(100, t))
                
                intersection = (p1[0] + t * d1[0], p1[1] + t * d1[1])
                setback_coords.append(intersection)
        
        # Close polygon
        if setback_coords and setback_coords[0] != setback_coords[-1]:
            setback_coords.append(setback_coords[0])
        
        return setback_coords if len(setback_coords) >= 4 else local_boundary
    
    def _determine_access_corridor(
        self,
        property_area: float,
        buildable_area: float
    ) -> Dict[str, Any]:
        """Determine if an access corridor is required"""
        
        # Larger plots typically need vehicle access
        needs_vehicle_access = property_area > 200  # sqm
        
        return {
            "required": needs_vehicle_access,
            "width_m": 3.5 if needs_vehicle_access else 0,
            "description": "Vehicle access corridor for driveway" if needs_vehicle_access else "",
            "affects_buildable_area": False  # Usually within setback zone
        }

