"""
Property line classifier

Classifies property line into front, rear, and side segments
using angle-based algorithm with road proximity
"""

import math
from typing import List, Dict, Any, Optional, Tuple, Set
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
from ...config import get_config


class PropertyLineClassifier:
    """Classifies property line into front, rear, and side segments"""
    
    # Colors for visualization
    FRONT_COLOR = "#FF0000"  # Red
    REAR_COLOR = "#0000FF"   # Blue
    SIDE_COLOR = "#00FF00"   # Green (both left and right sides use same color)
    
    def classify(
        self,
        property_line: PropertyLine,
        osm_roads: List[Dict[str, Any]],
        osm_buildings: List[Dict[str, Any]],
        debug: bool = False
    ) -> PropertyLine:
        """
        Classify property line into front, rear, and side segments
        
        Algorithm:
        1. Find closest named road from centroid point, calculate distance_from_road
        2. Merge all road segments with the same name into one road
        3. Find the closest road segment from the merged road to the centroid (once for all edges)
        4. For each edge:
           a. Use the same pre-selected road segment for all edges
           b. Calculate angle between edge and road segment (-90 to 90 degrees)
           c. Calculate minimum distance (edge line to road segment)
        5. Convert angles to absolute values (0~90 degrees):
           a. Convert all raw angles (-90~90) to absolute values (0~90)
        6. Normalize absolute angles:
           a. Find the edge with minimum absolute angle (most parallel to road, in 0~90 range)
           b. Subtract this minimum from all absolute angles (normalize so minimum becomes 0)
        7. Classification using normalized absolute angles:
           - If normalized absolute angle <= configured threshold (default 30 degrees): front or rear
           - Otherwise: side edge
        7. For front/rear edges:
           - If min_distance < distance_from_road: front
           - If min_distance >= distance_from_road: rear
        
        Args:
            property_line: PropertyLine object to classify
            osm_roads: List of OSM road features
            osm_buildings: List of OSM building features (not used but kept for compatibility)
            
        Returns:
            PropertyLine with classified segments
        """
        if not SHAPELY_AVAILABLE:
            logger.warning("Shapely not available - cannot classify property line")
            return property_line
        
        try:
            # Remove duplicate closing point for processing
            coords = property_line.coordinates
            coords_clean = coords[:-1] if (coords[0] == coords[-1] and len(coords) > 1) else coords
            
            if len(coords_clean) < 4:
                logger.warning("Property line has too few points for classification")
                return property_line
            
            # Step 1: Find closest NAMED road from centroid (skip unnamed roads)
            prop_poly = Polygon(coords_clean)
            centroid = prop_poly.centroid
            centroid_point = Point(centroid.x, centroid.y)
            
            closest_road, distance_from_road = self._find_closest_named_road(
                centroid_point, osm_roads
            )
            
            if not closest_road:
                logger.warning("No named roads found for classification")
                return property_line
            
            # Get angle threshold from config
            config = get_config()
            angle_threshold = config.uk_regulatory.classifier_parallel_angle_threshold
            
            road_name = closest_road.get('name', 'Unnamed')
            
            # Merge all road segments with the same name into one road
            merged_road = self._merge_roads_by_name(road_name, osm_roads, closest_road)
            
            # Find the closest road segment from the merged road to the centroid (once for all edges)
            closest_road_segment, closest_segment_idx = self._find_closest_road_segment_to_point(
                merged_road, centroid_point
            )
            
            if debug:
                logger.info(f"=== PROPERTY LINE CLASSIFICATION DEBUG ===")
                logger.info(f"Using angle threshold: ±{angle_threshold}° (from config)")
                logger.info(f"Closest road: {road_name}, centroid-to-road distance: {distance_from_road:.2f}m")
                logger.info(f"Merged {len([r for r in osm_roads if r.get('name') == road_name])} segments into one road")
                if closest_road_segment:
                    road_coords = merged_road.get("centerline", {}).get("coordinates", [])
                    if closest_segment_idx >= 0 and closest_segment_idx < len(road_coords) - 1:
                        seg_start = road_coords[closest_segment_idx]
                        seg_end = road_coords[closest_segment_idx + 1]
                        logger.info(f"Selected road segment {closest_segment_idx} for all edges:")
                        logger.info(f"  Segment start: ({seg_start[0]:.6f}, {seg_start[1]:.6f})")
                        logger.info(f"  Segment end:   ({seg_end[0]:.6f}, {seg_end[1]:.6f})")
                logger.info(f"Property line has {len(coords_clean)} edges")
                logger.info("")
            
            if closest_road_segment is None:
                logger.warning("No road segment found from merged road - cannot classify")
                return property_line
            
            # Step 2: First pass - calculate all angles and distances for all edges
            edge_data = []  # Store raw angle, absolute angle, distance, and edge info
            
            for i in range(len(coords_clean)):
                j = (i + 1) % len(coords_clean)
                edge_start = coords_clean[i]
                edge_end = coords_clean[j]
                edge_line = LineString([edge_start, edge_end])
                
                # Calculate angle and distance using the pre-selected road segment
                raw_angle, min_distance = self._analyze_edge_road_relationship(
                    edge_line, closest_road_segment, edge_index=i, debug=False
                )
                
                # Convert to absolute angle (0~90)
                abs_angle = abs(raw_angle)
                
                edge_data.append({
                    "index": i,
                    "start": edge_start,
                    "end": edge_end,
                    "edge_line": edge_line,
                    "raw_angle": raw_angle,
                    "abs_angle": abs_angle,
                    "min_distance": min_distance
                })
            
            # Find the minimum absolute angle (most parallel edge, now in 0~90 range)
            if not edge_data:
                logger.warning("No edge data found")
                return property_line
            
            min_abs_angle = min(edge["abs_angle"] for edge in edge_data)
            
            if debug:
                logger.info(f"Minimum absolute angle found: {min_abs_angle:.2f}°")
                logger.info(f"Normalizing all absolute angles by subtracting {min_abs_angle:.2f}°")
                logger.info("")
            
            # Step 3: Normalize all absolute angles and classify
            front_edge_indices = []
            rear_edge_indices = []
            side_edge_indices = []
            
            # Store debug info for each edge
            edge_debug_info = []
            
            for edge in edge_data:
                i = edge["index"]
                edge_start = edge["start"]
                edge_end = edge["end"]
                raw_angle = edge["raw_angle"]
                abs_angle = edge["abs_angle"]
                min_distance = edge["min_distance"]
                
                # Normalize absolute angle by subtracting the minimum absolute angle
                normalized_abs_angle = abs_angle - min_abs_angle
                
                # Calculate edge length for logging
                edge_length = math.sqrt(
                    (edge_end[0] - edge_start[0])**2 + 
                    (edge_end[1] - edge_start[1])**2
                )
                avg_lat = (edge_start[1] + edge_end[1]) / 2
                m_per_deg_lat = 111000
                m_per_deg_lon = 111000 * abs(math.cos(math.radians(avg_lat)))
                m_per_deg = (m_per_deg_lat + m_per_deg_lon) / 2
                edge_length_m = edge_length * m_per_deg
                
                if debug:
                    logger.info(f"--- Edge {i} ---")
                    logger.info(f"  Start: ({edge_start[0]:.6f}, {edge_start[1]:.6f})")
                    logger.info(f"  End:   ({edge_end[0]:.6f}, {edge_end[1]:.6f})")
                    logger.info(f"  Length: {edge_length_m:.2f}m")
                    logger.info(f"  Raw angle: {raw_angle:.2f}°")
                    logger.info(f"  Absolute angle: {abs_angle:.2f}°")
                    logger.info(f"  Normalized absolute angle: {normalized_abs_angle:.2f}° (abs - {min_abs_angle:.2f}°)")
                    logger.info(f"  Minimum distance (edge to road segment): {min_distance:.2f}m")
                    logger.info(f"  Centroid-to-road distance: {distance_from_road:.2f}m")
                
                # Classification based on normalized absolute angle (using config threshold)
                # Since we're using absolute angles, we only check if normalized_abs_angle <= threshold
                if normalized_abs_angle <= angle_threshold:
                    # Front or rear edge (parallel to road)
                    if min_distance < distance_from_road:
                        classification = "FRONT"
                        reason = f"normalized_abs_angle={normalized_abs_angle:.2f}° (parallel, <= {angle_threshold}°) AND min_distance={min_distance:.2f}m < centroid_distance={distance_from_road:.2f}m"
                        front_edge_indices.append(i)
                    else:
                        classification = "REAR"
                        reason = f"normalized_abs_angle={normalized_abs_angle:.2f}° (parallel, <= {angle_threshold}°) AND min_distance={min_distance:.2f}m >= centroid_distance={distance_from_road:.2f}m"
                        rear_edge_indices.append(i)
                else:
                    # Side edge (perpendicular to road)
                    classification = "SIDE"
                    reason = f"normalized_abs_angle={normalized_abs_angle:.2f}° (NOT parallel, > {angle_threshold}°)"
                    side_edge_indices.append(i)
                
                # Store debug info (include raw, absolute, and normalized absolute angles)
                edge_debug_info.append({
                    "edge_index": i,
                    "start": edge_start,
                    "end": edge_end,
                    "raw_angle": raw_angle,
                    "abs_angle": abs_angle,
                    "normalized_abs_angle": normalized_abs_angle,
                    "min_distance": min_distance,
                    "distance_from_road": distance_from_road,
                    "classification": classification,
                    "reason": reason
                })
                
                if debug:
                    logger.info(f"  ✅ Classification: {classification}")
                    logger.info(f"  Reason: {reason}")
                    logger.info("")
            
            if debug:
                logger.info(f"=== SUMMARY ===")
                logger.info(f"Front edges: {len(front_edge_indices)} (indices: {front_edge_indices})")
                logger.info(f"Rear edges: {len(rear_edge_indices)} (indices: {rear_edge_indices})")
                logger.info(f"Side edges: {len(side_edge_indices)} (indices: {side_edge_indices})")
                logger.info(f"Total: {len(front_edge_indices) + len(rear_edge_indices) + len(side_edge_indices)}/{len(coords_clean)} edges classified")
                logger.info("")
            
            # Step 3: Create segments
            front_segment = self._create_segment_from_indices(
                coords_clean, front_edge_indices, SegmentType.FRONT, self.FRONT_COLOR
            )
            rear_segment = self._create_segment_from_indices(
                coords_clean, rear_edge_indices, SegmentType.REAR, self.REAR_COLOR
            )
            
            # For sides, we don't distinguish left/right - just create one side segment
            # But we still need to create left_side and right_side for the model
            # Split sides into left and right based on polygon order
            left_side_indices, right_side_indices = self._split_sides(
                coords_clean, front_edge_indices, rear_edge_indices, side_edge_indices
            )
            
            left_side = self._create_segment_from_indices(
                coords_clean, left_side_indices, SegmentType.LEFT_SIDE, self.SIDE_COLOR
            )
            right_side = self._create_segment_from_indices(
                coords_clean, right_side_indices, SegmentType.RIGHT_SIDE, self.SIDE_COLOR
            )
            
            # Assign segments
            property_line.front = front_segment
            property_line.rear = rear_segment
            property_line.left_side = left_side
            property_line.right_side = right_side
            
            # Store debug info in metadata
            property_line.metadata["edge_debug_info"] = edge_debug_info
            property_line.metadata["distance_from_road"] = distance_from_road
            property_line.metadata["closest_road_name"] = merged_road.get('name', 'Unnamed')
            
            # Verify all edges are classified
            total_edges = len(coords_clean)
            classified_edges = (len(front_edge_indices) + len(rear_edge_indices) + 
                              len(side_edge_indices))
            if classified_edges != total_edges:
                logger.warning(f"Classification incomplete: {classified_edges}/{total_edges} edges classified")
            
            return property_line
            
        except Exception as e:
            logger.warning(f"Failed to classify property line: {e}")
            import traceback
            traceback.print_exc()
            return property_line
    
    def _find_closest_named_road(
        self,
        centroid_point: Point,
        osm_roads: List[Dict[str, Any]]
    ) -> Tuple[Optional[Dict[str, Any]], float]:
        """
        Find the closest NAMED road from centroid point
        Skips unnamed roads or roads without names
        
        Returns:
            Tuple of (closest_road_dict, distance_in_meters)
        """
        if not osm_roads:
            return None, float('inf')
        
        closest_road = None
        min_distance = float('inf')
        
        for road in osm_roads:
            # Skip unnamed roads or roads without names
            road_name = road.get("name", "").strip()
            road_name_lower = road_name.lower() if road_name else ""
            
            # Skip if road is unnamed, empty, or generic "road"
            if not road_name or road_name_lower in ["unnamed road", "unnamed", "road", ""]:
                continue
            
            centerline = road.get("centerline", {})
            if not centerline:
                continue
            road_coords = centerline.get("coordinates", [])
            if len(road_coords) < 2:
                continue
            
            try:
                road_line = LineString(road_coords)
                distance_deg = centroid_point.distance(road_line)
                
                # Convert to meters using centroid latitude
                avg_lat = centroid_point.y
                m_per_deg_lat = 111000
                m_per_deg_lon = 111000 * abs(math.cos(math.radians(avg_lat)))
                m_per_deg = (m_per_deg_lat + m_per_deg_lon) / 2
                distance_m = distance_deg * m_per_deg
                
                if distance_m < min_distance:
                    min_distance = distance_m
                    closest_road = road
            except Exception as e:
                logger.debug(f"Error processing road {road_name}: {e}")
                continue
        
        return closest_road, min_distance
    
    def _merge_roads_by_name(
        self,
        road_name: str,
        osm_roads: List[Dict[str, Any]],
        reference_road: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Merge all road segments with the same name into one continuous road
        
        Args:
            road_name: Name of the road to merge
            osm_roads: List of all OSM roads
            reference_road: The reference road (closest one)
        
        Returns:
            Merged road dictionary with combined centerline coordinates
        """
        # Find all roads with the same name
        same_name_roads = [
            road for road in osm_roads 
            if road.get('name') == road_name
        ]
        
        if len(same_name_roads) == 1:
            # Only one segment, return as is
            return reference_road
        
        # Collect all coordinate lists from all segments
        road_segments = []
        for road in same_name_roads:
            centerline = road.get("centerline", {})
            if not centerline:
                continue
            road_coords = centerline.get("coordinates", [])
            if road_coords and len(road_coords) >= 2:
                road_segments.append(road_coords)
        
        if not road_segments:
            return reference_road
        
        # Merge segments by connecting endpoints
        # Start with the reference road's coordinates
        merged_coords = []
        used_segments = set()
        
        # Find reference road segment
        ref_centerline = reference_road.get("centerline", {})
        ref_coords = ref_centerline.get("coordinates", []) if ref_centerline else []
        
        if ref_coords:
            merged_coords = ref_coords.copy()
            # Mark which segment was used (find by matching coordinates)
            for idx, seg in enumerate(road_segments):
                if (len(seg) >= 2 and 
                    abs(seg[0][0] - ref_coords[0][0]) < 1e-6 and 
                    abs(seg[0][1] - ref_coords[0][1]) < 1e-6):
                    used_segments.add(idx)
                    break
        
        # Try to connect other segments
        # Simple approach: just concatenate all segments, removing duplicates
        all_coords = []
        for seg in road_segments:
            all_coords.extend(seg)
        
        # Remove duplicate consecutive points
        cleaned_coords = []
        for i, coord in enumerate(all_coords):
            if i == 0:
                cleaned_coords.append(coord)
            else:
                prev = cleaned_coords[-1]
                # If points are different (more than 1e-6 degrees apart), add it
                if abs(prev[0] - coord[0]) > 1e-6 or abs(prev[1] - coord[1]) > 1e-6:
                    cleaned_coords.append(coord)
        
        # Create merged road
        merged_road = {
            "id": reference_road.get("id", ""),
            "name": road_name,
            "type": reference_road.get("type", "residential"),
            "centerline": {
                "coordinates": cleaned_coords
            }
        }
        
        # Log only if debug is enabled (but we don't have debug param here, so use logger.debug)
        logger.debug(f"Merged {len(same_name_roads)} segments of '{road_name}' into {len(cleaned_coords)} points")
        
        return merged_road
    
    def _find_closest_road_segment_to_point(
        self,
        road: Dict[str, Any],
        point: Point
    ) -> Tuple[Optional[LineString], int]:
        """
        Find the closest road segment from the merged road to a point (e.g., centroid)
        
        Args:
            road: Merged road dictionary with centerline coordinates
            point: Point to find closest segment to (e.g., property centroid)
        
        Returns:
            Tuple of (closest_road_segment_linestring, segment_index)
        """
        centerline = road.get("centerline", {})
        if not centerline:
            return None, -1
        
        road_coords = centerline.get("coordinates", [])
        if len(road_coords) < 2:
            return None, -1
        
        # Find closest road segment to the point
        closest_segment = None
        min_distance_deg = float('inf')
        closest_segment_idx = -1
        
        # Check each segment of the road
        for i in range(len(road_coords) - 1):
            road_segment = LineString([road_coords[i], road_coords[i + 1]])
            distance = point.distance(road_segment)
            
            if distance < min_distance_deg:
                min_distance_deg = distance
                closest_segment = road_segment
                closest_segment_idx = i
        
        return closest_segment, closest_segment_idx
    
    def _analyze_edge_road_relationship(
        self,
        edge_line: LineString,
        road_segment: LineString,
        edge_index: Optional[int] = None,
        debug: bool = False
    ) -> Tuple[float, float]:
        """
        Analyze relationship between edge and a pre-selected road segment
        
        Algorithm:
        1. Calculate angle between edge and road segment (-90 to 90 degrees)
        2. Calculate minimum distance (edge line to road segment)
        
        Args:
            edge_line: Property line edge as LineString
            road_segment: Pre-selected road segment as LineString (same for all edges)
            edge_index: Optional edge index for debugging
            debug: Enable debug logging
        
        Returns:
            Tuple of (angle_degrees, min_distance_meters)
        """
        if road_segment is None or len(road_segment.coords) < 2:
            return 0.0, float('inf')
        
        # Calculate minimum distance from edge to road segment
        min_distance_deg = edge_line.distance(road_segment)
        
        # Get road segment coordinates for distance conversion
        road_seg_coords = list(road_segment.coords)
        road_seg_start = road_seg_coords[0]
        road_seg_end = road_seg_coords[-1]
        avg_lat = (road_seg_start[1] + road_seg_end[1]) / 2
        m_per_deg_lat = 111000
        m_per_deg_lon = 111000 * abs(math.cos(math.radians(avg_lat)))
        m_per_deg = (m_per_deg_lat + m_per_deg_lon) / 2
        min_distance_m = min_distance_deg * m_per_deg
        
        if debug:
            logger.info(f"  Using pre-selected road segment for all edges")
            logger.info(f"    Road segment start: ({road_seg_start[0]:.6f}, {road_seg_start[1]:.6f})")
            logger.info(f"    Road segment end:   ({road_seg_end[0]:.6f}, {road_seg_end[1]:.6f})")
            logger.info(f"    Minimum distance (edge to road segment): {min_distance_m:.2f}m")
        
        # Calculate angle between edge and road segment
        angle = self._calculate_angle_between_lines(edge_line, road_segment, debug=debug)
        
        return angle, min_distance_m
    
    def _calculate_angle_between_lines(
        self,
        line1: LineString,
        line2: LineString,
        debug: bool = False
    ) -> float:
        """
        Calculate angle between two lines in degrees (-90 to 90)
        
        Returns angle in degrees, where:
        - 0 degrees = parallel
        - 90 degrees = perpendicular
        - Negative = opposite orientation
        """
        if len(line1.coords) < 2 or len(line2.coords) < 2:
            return 0.0
        
        # Get direction vectors
        dx1 = line1.coords[-1][0] - line1.coords[0][0]
        dy1 = line1.coords[-1][1] - line1.coords[0][1]
        length1 = math.sqrt(dx1**2 + dy1**2)
        
        dx2 = line2.coords[-1][0] - line2.coords[0][0]
        dy2 = line2.coords[-1][1] - line2.coords[0][1]
        length2 = math.sqrt(dx2**2 + dy2**2)
        
        if length1 == 0 or length2 == 0:
            return 0.0
        
        if debug:
            logger.info(f"    Edge direction vector: ({dx1:.6f}, {dy1:.6f}), length: {length1:.6f}")
            logger.info(f"    Road direction vector: ({dx2:.6f}, {dy2:.6f}), length: {length2:.6f}")
        
        # Normalize vectors
        dx1_norm = dx1 / length1
        dy1_norm = dy1 / length1
        dx2_norm = dx2 / length2
        dy2_norm = dy2 / length2
        
        if debug:
            logger.info(f"    Normalized edge vector: ({dx1_norm:.6f}, {dy1_norm:.6f})")
            logger.info(f"    Normalized road vector: ({dx2_norm:.6f}, {dy2_norm:.6f})")
        
        # Calculate angle using dot product
        # dot = cos(angle) when vectors are normalized
        dot_product = dx1_norm * dx2_norm + dy1_norm * dy2_norm
        
        # Clamp to [-1, 1] to avoid numerical errors
        dot_product = max(-1.0, min(1.0, dot_product))
        
        if debug:
            logger.info(f"    Dot product: {dot_product:.6f}")
        
        # Calculate angle in radians, then convert to degrees
        angle_rad = math.acos(dot_product)
        angle_deg = math.degrees(angle_rad)
        
        if debug:
            logger.info(f"    Initial angle: {angle_deg:.2f}°")
        
        # Adjust to -90 to 90 range
        # If angle > 90, subtract 180 to get the acute angle
        if angle_deg > 90:
            angle_deg = 180 - angle_deg
            if debug:
                logger.info(f"    Adjusted to acute angle: {angle_deg:.2f}°")
        
        # Determine sign based on cross product (for orientation)
        cross_product = dx1_norm * dy2_norm - dy1_norm * dx2_norm
        if cross_product < 0:
            angle_deg = -angle_deg
            if debug:
                logger.info(f"    Applied sign (cross product < 0): {angle_deg:.2f}°")
        
        return angle_deg
    
    def _split_sides(
        self,
        coords: List[List[float]],
        front_edge_indices: List[int],
        rear_edge_indices: List[int],
        side_edge_indices: List[int]
    ) -> Tuple[List[int], List[int]]:
        """
        Split side edges into left and right based on polygon order
        
        Left side: edges from front_end to rear_start
        Right side: edges from rear_end to front_start
        """
        if not side_edge_indices:
            return [], []
        
        if not front_edge_indices or not rear_edge_indices:
            # If no front/rear, split sides in half
            mid = len(side_edge_indices) // 2
            return side_edge_indices[:mid], side_edge_indices[mid:]
        
        # Get front and rear boundary points
        front_edges_list = sorted(front_edge_indices)
        rear_edges_list = sorted(rear_edge_indices)
        
        # Front end point (last point of front segment)
        front_end_idx = (front_edges_list[-1] + 1) % len(coords)
        
        # Front start point (first point of front segment)
        front_start_idx = front_edges_list[0]
        
        # Rear start point (first point of rear segment)
        rear_start_idx = rear_edges_list[0]
        
        # Rear end point (last point of rear segment)
        rear_end_idx = (rear_edges_list[-1] + 1) % len(coords)
        
        # Left side: path from front_end to rear_start
        left_side_indices = []
        current = front_end_idx
        target = rear_start_idx
        
        while current != target:
            if current in side_edge_indices:
                left_side_indices.append(current)
            current = (current + 1) % len(coords)
        
        # Right side: path from rear_end to front_start
        right_side_indices = []
        current = rear_end_idx
        target = front_start_idx
        
        while current != target:
            if current in side_edge_indices:
                right_side_indices.append(current)
            current = (current + 1) % len(coords)
        
        return left_side_indices, right_side_indices
    
    def _create_segment_from_indices(
        self,
        coords: List[List[float]],
        edge_indices: List[int],
        segment_type: SegmentType,
        color: str
    ) -> Optional[PropertyLineSegment]:
        """
        Create a PropertyLineSegment from edge indices
        
        Simply stores all edge indices without checking for consecutive groups.
        The get_coordinates() method will handle reconstruction when needed.
        """
        if not edge_indices:
            return None
        
        # Calculate total length from all edges
        total_length = 0.0
        for idx in edge_indices:
            start_point = coords[idx]
            end_idx = (idx + 1) % len(coords)
            end_point = coords[end_idx]
            
            # Calculate edge length
            dx = end_point[0] - start_point[0]
            dy = end_point[1] - start_point[1]
            edge_length = math.sqrt(dx**2 + dy**2)
            
            # Convert to meters
            avg_lat = (start_point[1] + end_point[1]) / 2
            m_per_deg_lat = 111000
            m_per_deg_lon = 111000 * abs(math.cos(math.radians(avg_lat)))
            m_per_deg = (m_per_deg_lat + m_per_deg_lon) / 2
            total_length += edge_length * m_per_deg
        
        return PropertyLineSegment(
            segment_type=segment_type,
            edge_indices=edge_indices,  # Store all edge indices as-is
            length_m=total_length,
            color=color,
            metadata={"num_edges": len(edge_indices)}
        )
    
    def _points_match(self, p1: List[float], p2: List[float], tolerance: float = 1e-6) -> bool:
        """Check if two points match"""
        return abs(p1[0] - p2[0]) < tolerance and abs(p1[1] - p2[1]) < tolerance
