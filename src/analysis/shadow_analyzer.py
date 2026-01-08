"""
Shadow analysis using sun position calculations
Calculates shadow impact on plot from neighboring buildings
"""

import math
from datetime import datetime, timedelta
from typing import List, Dict, Any, Tuple, Optional
from loguru import logger

# Lazy import pvlib to avoid scipy/numpy compatibility issues at import time
PVLIB_AVAILABLE = None
_pvlib_module = None

def _check_pvlib():
    """Lazy check for pvlib availability"""
    global PVLIB_AVAILABLE, _pvlib_module
    if PVLIB_AVAILABLE is None:
        try:
            import pvlib
            _pvlib_module = pvlib
            PVLIB_AVAILABLE = True
        except (ImportError, ValueError) as e:
            # ValueError can occur with scipy/numpy version mismatches
            PVLIB_AVAILABLE = False
            _pvlib_module = None
            logger.warning(f"pvlib not available ({type(e).__name__}: {e}) - using simplified sun position calculations")
    return PVLIB_AVAILABLE

from .geometry_utils import GeometryUtils
from ..config import get_config


class ShadowAnalyzer:
    """
    Analyze shadow patterns on a plot from neighboring buildings
    
    Uses pvlib for accurate sun position calculations when available,
    falls back to simplified calculations otherwise.
    """
    
    def __init__(self, timezone: str = "Europe/London"):
        self.timezone = timezone
    
    def analyze(
        self,
        plot_centroid: Tuple[float, float],  # (lon, lat)
        plot_boundary: List[List[float]],
        neighbor_buildings: List[Dict[str, Any]],
        property_line_obj: Optional[Any] = None
    ) -> Dict[str, Any]:
        """
        Perform comprehensive shadow analysis
        
        Returns facade scores, shadow hours, neighbor shadow angles, and sunlight score
        
        Args:
            plot_centroid: (lon, lat) of plot center
            plot_boundary: List of [lon, lat] coordinates
            neighbor_buildings: List of neighbor building dicts
            property_line_obj: Optional PropertyLine object with front/rear/side segments
        """
        lon, lat = plot_centroid
        logger.info(f"Analyzing shadows for plot at ({lat:.6f}, {lon:.6f})")
        logger.info(f"  Neighbor buildings: {len(neighbor_buildings)}")
        logger.info(f"  Plot boundary points: {len(plot_boundary)}")
        
        # Check pvlib availability
        pvlib_available = _check_pvlib()
        logger.info(f"  pvlib available: {pvlib_available}")
        
        # Key dates for analysis
        analysis_dates = {
            "winter_solstice": datetime(2024, 12, 21),
            "summer_solstice": datetime(2024, 6, 21),
            "equinox": datetime(2024, 3, 20)
        }
        
        # Calculate sun positions for each date
        logger.info("Calculating sun hours for key dates...")
        shadow_hours = {}
        for period, date in analysis_dates.items():
            hours = self._calculate_sun_hours(lat, lon, date)
            shadow_hours[period] = hours
            logger.info(f"  {period} ({date.strftime('%Y-%m-%d')}): {hours:.1f} hours")
        
        # Calculate facade scores (already accounts for neighbor shadows)
        logger.info("Calculating facade scores...")
        facade_scores = self._calculate_facade_scores(
            lat, lon, neighbor_buildings, plot_boundary
        )
        logger.info(f"  Facade scores calculated: {len(facade_scores)} facades")
        for facade, scores in facade_scores.items():
            logger.info(f"    {facade}: winter={scores['winter_avg']:.2f}, summer={scores['summer_avg']:.2f}, annual={scores['annual_avg']:.2f}")
        
        # Calculate neighbor shadow angles
        logger.info("Calculating neighbor shadow angles...")
        neighbor_angles = self._calculate_neighbor_shadow_angles(
            plot_centroid, plot_boundary, neighbor_buildings
        )
        logger.info(f"  Neighbor shadow angles: {len(neighbor_angles)} directions")
        for direction, angle in neighbor_angles.items():
            logger.info(f"    {direction}: {angle:.0f}°")
        
        # Calculate sunlight score (combines building direction + neighbor shadows)
        logger.info("Calculating sunlight score...")
        sunlight_score = self._calculate_sunlight_score(
            plot_boundary, facade_scores, property_line_obj
        )
        logger.info(f"  Sunlight score: {sunlight_score:.3f}")
        
        # Determine best solar facade
        best_facade = max(
            facade_scores.keys(),
            key=lambda f: facade_scores[f]["annual_avg"]
        )
        best_score = facade_scores[best_facade]["annual_avg"]
        logger.info(f"  Best solar facade: {best_facade} (score: {best_score:.2f})")
        
        result = {
            "facade_scores": facade_scores,
            "shadow_hours_per_day": shadow_hours,
            "neighbor_shadow_angles": neighbor_angles,
            "best_solar_facade": best_facade,
            "sunlight_score": sunlight_score,
            "computed_at": datetime.utcnow().isoformat(),
            "sun_path_params": {
                "latitude": lat,
                "longitude": lon,
                "timezone": self.timezone
            }
        }
        
        logger.info(f"Shadow analysis complete. Returning {len(result)} result keys")
        return result
    
    def _calculate_sun_hours(
        self,
        lat: float,
        lon: float,
        date: datetime
    ) -> float:
        """Calculate approximate sun hours for a given date"""
        if _check_pvlib():
            logger.debug(f"Using pvlib for sun hours calculation ({date.strftime('%Y-%m-%d')})")
            return self._calculate_sun_hours_pvlib(lat, lon, date)
        else:
            logger.debug(f"Using simplified calculation for sun hours ({date.strftime('%Y-%m-%d')})")
            return self._calculate_sun_hours_simple(lat, date)
    
    def _calculate_sun_hours_pvlib(
        self,
        lat: float,
        lon: float,
        date: datetime
    ) -> float:
        """Calculate sun hours using pvlib"""
        import pandas as pd
        
        # Create time range for the day
        times = pd.date_range(
            start=date.replace(hour=0, minute=0),
            end=date.replace(hour=23, minute=59),
            freq='15min',
            tz=self.timezone
        )
        logger.debug(f"  pvlib: Created time range with {len(times)} intervals (15min each)")
        
        # Get solar positions
        if not _check_pvlib() or _pvlib_module is None:
            logger.warning("  pvlib check failed, falling back to simplified calculation")
            return self._calculate_sun_hours_simple(lat, lon, date)
        
        logger.debug(f"  pvlib: Calling solarposition.get_solarposition(lat={lat:.6f}, lon={lon:.6f})")
        solar_pos = _pvlib_module.solarposition.get_solarposition(times, lat, lon)
        
        # Log some pvlib results
        if len(solar_pos) > 0:
            max_elevation = solar_pos['apparent_elevation'].max()
            min_elevation = solar_pos['apparent_elevation'].min()
            max_azimuth = solar_pos['azimuth'].max()
            min_azimuth = solar_pos['azimuth'].min()
            logger.debug(f"  pvlib results: elevation range [{min_elevation:.2f}°, {max_elevation:.2f}°], "
                        f"azimuth range [{min_azimuth:.2f}°, {max_azimuth:.2f}°]")
        
        # Count hours where sun is above horizon (elevation > 0)
        sun_above_horizon = solar_pos['apparent_elevation'] > 0
        sun_hours = sun_above_horizon.sum() * 0.25  # 15-minute intervals
        sun_above_count = sun_above_horizon.sum()
        logger.debug(f"  pvlib: {sun_above_count}/{len(times)} intervals above horizon = {sun_hours:.2f} hours")
        
        return round(sun_hours, 1)
    
    def _calculate_sun_hours_simple(self, lat: float, date: datetime) -> float:
        """Simplified sun hours calculation without pvlib"""
        # Day of year
        doy = date.timetuple().tm_yday
        logger.debug(f"  Simple calc: day of year = {doy}")
        
        # Solar declination angle (simplified)
        declination = 23.45 * math.sin(math.radians(360 * (284 + doy) / 365))
        logger.debug(f"  Simple calc: declination = {declination:.2f}°")
        
        # Convert to radians
        lat_rad = math.radians(lat)
        dec_rad = math.radians(declination)
        
        # Hour angle at sunrise/sunset
        cos_hour_angle = -math.tan(lat_rad) * math.tan(dec_rad)
        logger.debug(f"  Simple calc: cos(hour_angle) = {cos_hour_angle:.4f}")
        
        # Handle polar day/night
        if cos_hour_angle < -1:
            logger.debug(f"  Simple calc: polar day (24 hours)")
            return 24.0  # Midnight sun
        elif cos_hour_angle > 1:
            logger.debug(f"  Simple calc: polar night (0 hours)")
            return 0.0   # Polar night
        
        hour_angle = math.degrees(math.acos(cos_hour_angle))
        logger.debug(f"  Simple calc: hour_angle = {hour_angle:.2f}°")
        
        # Day length in hours
        day_length = 2 * hour_angle / 15
        logger.debug(f"  Simple calc: day_length = {day_length:.2f} hours")
        
        return round(day_length, 1)
    
    def _calculate_facade_scores(
        self,
        lat: float,
        lon: float,
        neighbor_buildings: List[Dict[str, Any]],
        plot_boundary: List[List[float]]
    ) -> Dict[str, Dict[str, float]]:
        """
        Calculate solar exposure scores for each facade (0-1 scale)
        
        1.0 = full sun exposure
        0.0 = fully shaded
        """
        # Convert plot boundary to local coords
        ref_lon, ref_lat = plot_boundary[0]
        local_boundary = GeometryUtils.degrees_to_local(plot_boundary, ref_lon, ref_lat)
        
        # Calculate polygon center in local coordinates
        center_deg = GeometryUtils.polygon_centroid(plot_boundary)
        center_local = GeometryUtils.degrees_to_local([list(center_deg)], ref_lon, ref_lat)[0]
        
        # Get plot edges (direction determined relative to polygon center)
        edges = GeometryUtils.get_polygon_edges(local_boundary, center=center_local)
        
        # Classify edges by cardinal direction
        facade_edges = {"north": [], "south": [], "east": [], "west": []}
        for edge in edges:
            direction = edge["direction"]
            # Facade faces opposite to edge direction
            if "north" in direction:
                facade_edges["south"].append(edge)
            elif "south" in direction:
                facade_edges["north"].append(edge)
            elif "east" in direction:
                facade_edges["west"].append(edge)
            elif "west" in direction:
                facade_edges["east"].append(edge)
        
        # Calculate base scores for each facade based on latitude
        # South facade gets best sun in northern hemisphere
        # Get base scores from config
        config = get_config()
        base_scores = config.facade_base_scores.copy()
        logger.debug(f"  Base facade scores (from config): {base_scores}")
        
        # Adjust for neighbor obstructions
        facade_scores = {}
        for facade, base in base_scores.items():
            # Check for neighbor buildings blocking this facade
            obstruction_factor = self._calculate_obstruction(
                facade, local_boundary, neighbor_buildings, ref_lon, ref_lat
            )
            logger.debug(f"  {facade} facade: obstruction_factor = {obstruction_factor:.3f}")
            
            winter_avg = base["winter"] * (1 - obstruction_factor * 0.5)
            summer_avg = base["summer"] * (1 - obstruction_factor * 0.3)
            annual_avg = (winter_avg + summer_avg) / 2
            
            logger.debug(f"  {facade} facade: base winter={base['winter']:.2f} → adjusted={winter_avg:.2f}, "
                        f"base summer={base['summer']:.2f} → adjusted={summer_avg:.2f}, annual={annual_avg:.2f}")
            
            facade_scores[facade] = {
                "winter_avg": round(winter_avg, 2),
                "summer_avg": round(summer_avg, 2),
                "annual_avg": round(annual_avg, 2)
            }
        
        logger.debug(f"  Final facade scores: {facade_scores}")
        return facade_scores
    
    def _calculate_sunlight_score(
        self,
        plot_boundary: List[List[float]],
        facade_scores: Dict[str, Dict[str, float]],
        property_line_obj: Optional[Any] = None
    ) -> float:
        """
        Calculate overall sunlight score (0-1) for building on buildable envelope
        
        Factors:
        1. Building direction: front/rear have big windows (higher score), sides have small windows (lower score)
        2. Neighbor shadows: already accounted for in facade_scores
        
        Args:
            plot_boundary: List of [lon, lat] coordinates
            facade_scores: Dict with annual_avg scores for each facade (already includes neighbor shadow effects)
            property_line_obj: Optional PropertyLine object with front/rear/side segments
        
        Returns:
            float: Sunlight score 0.0-1.0
        """
        # Get direction scores from config
        config = get_config()
        front_rear_score = config.facade_direction_scores.get("front", 0.9)
        side_score = config.facade_direction_scores.get("side", 0.4)
        
        # Default direction scores: front/rear = big windows, sides = small windows
        # These represent the window area/importance factor
        direction_scores = {
            "north": side_score,  # Default: side
            "south": side_score,  # Default: side
            "east": side_score,   # Default: side
            "west": side_score    # Default: side
        }
        
        # If property_line_obj is available, determine which facades are front/rear vs sides
        if property_line_obj and hasattr(property_line_obj, 'front') and property_line_obj.front:
            logger.debug("  Determining front/rear facades from property line segments...")
            
            # Get front and rear edge coordinates to determine their cardinal directions
            ref_lon, ref_lat = plot_boundary[0]
            local_boundary = GeometryUtils.degrees_to_local(plot_boundary, ref_lon, ref_lat)
            center_deg = GeometryUtils.polygon_centroid(plot_boundary)
            center_local = GeometryUtils.degrees_to_local([list(center_deg)], ref_lon, ref_lat)[0]
            
            # Get edges with directions
            edges = GeometryUtils.get_polygon_edges(local_boundary, center=center_local)
            
            # Map edge indices to directions
            # Edges are returned in the same order as coordinates, so index i corresponds to edge i
            edge_to_direction = {}
            for i, edge in enumerate(edges):
                direction = edge.get("direction", "")
                # Facade faces opposite to edge direction
                if "north" in direction:
                    edge_to_direction[i] = "south"
                elif "south" in direction:
                    edge_to_direction[i] = "north"
                elif "east" in direction:
                    edge_to_direction[i] = "west"
                elif "west" in direction:
                    edge_to_direction[i] = "east"
            
            # Check front edges
            if property_line_obj.front and hasattr(property_line_obj.front, 'edge_indices'):
                front_indices = property_line_obj.front.edge_indices
                for idx in front_indices:
                    facade_dir = edge_to_direction.get(idx)
                    if facade_dir:
                        direction_scores[facade_dir] = front_rear_score  # Front = big windows
                        logger.debug(f"    Front facade: {facade_dir} (score: {front_rear_score})")
            
            # Check rear edges
            if property_line_obj.rear and hasattr(property_line_obj.rear, 'edge_indices'):
                rear_indices = property_line_obj.rear.edge_indices
                for idx in rear_indices:
                    facade_dir = edge_to_direction.get(idx)
                    if facade_dir:
                        direction_scores[facade_dir] = front_rear_score  # Rear = big windows
                        logger.debug(f"    Rear facade: {facade_dir} (score: {front_rear_score})")
            
            # Sides get lower score (set explicitly from config)
            if property_line_obj.left_side and hasattr(property_line_obj.left_side, 'edge_indices'):
                side_indices = property_line_obj.left_side.edge_indices
                for idx in side_indices:
                    facade_dir = edge_to_direction.get(idx)
                    if facade_dir:
                        direction_scores[facade_dir] = side_score  # Side = small windows
                        logger.debug(f"    Side facade: {facade_dir} (score: {side_score})")
            
            if property_line_obj.right_side and hasattr(property_line_obj.right_side, 'edge_indices'):
                side_indices = property_line_obj.right_side.edge_indices
                for idx in side_indices:
                    facade_dir = edge_to_direction.get(idx)
                    if facade_dir:
                        direction_scores[facade_dir] = side_score  # Side = small windows
                        logger.debug(f"    Side facade: {facade_dir} (score: {side_score})")
        else:
            logger.debug("  No property line segments available, using default direction scores")
        
        # Calculate weighted sunlight score
        # Combine direction_score (window importance) with facade_score (solar exposure + neighbor shadows)
        # Normalize the final result to 0-1 scale based on maximum possible weighted score
        
        # Maximum possible annual_avg for each facade (base scores with no obstruction)
        max_facade_scores = {
            "north": (0.15 + 0.45) / 2,  # 0.3
            "south": (0.85 + 0.95) / 2,  # 0.9
            "east": (0.55 + 0.70) / 2,   # 0.625
            "west": (0.50 + 0.65) / 2     # 0.575
        }
        
        # Calculate actual weighted score
        total_score = 0.0
        
        # Calculate maximum possible weighted score using ACTUAL direction_scores
        # (not assuming all are 0.9, since only front/rear are 0.9, sides are 0.4)
        max_total_score = 0.0
        
        for facade in ["north", "south", "east", "west"]:
            if facade in facade_scores:
                direction_score = direction_scores.get(facade, 0.5)
                facade_score = facade_scores[facade]["annual_avg"]
                max_facade_score = max_facade_scores.get(facade, 1.0)
                
                # Actual weighted score (preserves relative importance of facades)
                weighted_score = direction_score * facade_score
                total_score += weighted_score
                
                # Maximum possible weighted score using the ACTUAL direction_score for this facade
                # This reflects the actual building orientation (front/rear vs sides)
                max_weighted_score = direction_score * max_facade_score
                max_total_score += max_weighted_score
                
                logger.debug(f"    {facade}: direction={direction_score:.2f}, facade={facade_score:.2f}, "
                           f"weighted={weighted_score:.3f}, max_weighted={max_weighted_score:.3f}")
        
        # Normalize by maximum possible score to get 0-1 range
        # This preserves the relative importance: blocking south facade hurts more than blocking north
        # And uses the actual building orientation (front/rear/sides) for maximum calculation
        if max_total_score > 0:
            sunlight_score = total_score / max_total_score
        else:
            # Fallback: average of facade scores
            sunlight_score = sum(f["annual_avg"] for f in facade_scores.values()) / len(facade_scores) if facade_scores else 0.5
        
        # Ensure score is in 0-1 range
        sunlight_score = max(0.0, min(1.0, sunlight_score))
        
        logger.debug(f"  Sunlight score calculation: total_score={total_score:.3f}, max_total_score={max_total_score:.3f}, final={sunlight_score:.3f}")
        return round(sunlight_score, 3)
    
    def _calculate_obstruction(
        self,
        facade: str,
        plot_boundary: List[Tuple[float, float]],
        neighbor_buildings: List[Dict[str, Any]],
        ref_lon: float,
        ref_lat: float
    ) -> float:
        """
        Calculate obstruction factor (0-1) from neighbor buildings for a facade
        """
        # Define direction vectors for each facade
        facade_directions = {
            "north": (0, 1),    # Check buildings to the north
            "south": (0, -1),   # Check buildings to the south
            "east": (1, 0),     # Check buildings to the east
            "west": (-1, 0)     # Check buildings to the west
        }
        
        direction = facade_directions.get(facade, (0, 0))
        if direction == (0, 0):
            return 0.0
        
        # Calculate plot centroid
        centroid = GeometryUtils.polygon_centroid(
            [(p[0], p[1]) for p in plot_boundary]
        )
        
        max_obstruction = 0.0
        building_count = 0
        
        logger.debug(f"    Calculating obstruction for {facade} facade, checking {len(neighbor_buildings)} buildings")
        
        for building in neighbor_buildings:
            # Get building footprint in local coords
            building_coords = building.get("footprint", {}).get("coordinates", [[]])
            if not building_coords or not building_coords[0]:
                continue
            
            building_local = GeometryUtils.degrees_to_local(
                building_coords[0], ref_lon, ref_lat
            )
            
            building_centroid = GeometryUtils.polygon_centroid(building_local)
            
            # Vector from plot to building
            to_building = (
                building_centroid[0] - centroid[0],
                building_centroid[1] - centroid[1]
            )
            
            # Check if building is in the direction of the facade
            dot_product = direction[0] * to_building[0] + direction[1] * to_building[1]
            
            if dot_product > 0:  # Building is in the facade direction
                distance = math.sqrt(to_building[0]**2 + to_building[1]**2)
                height = building.get("height_m", 8.0)
                building_id = building.get("id", "unknown")
                
                if distance > 0:
                    # Calculate angular size of building
                    angle = math.degrees(math.atan(height / distance))
                    
                    # Closer and taller buildings cause more obstruction
                    # Maximum obstruction at ~25m distance for 10m building
                    obstruction = min(1.0, angle / 45.0)
                    max_obstruction = max(max_obstruction, obstruction)
                    building_count += 1
                    
                    logger.debug(f"      Building {building_id}: distance={distance:.1f}m, height={height:.1f}m, "
                                f"angle={angle:.2f}°, obstruction={obstruction:.3f}")
        
        logger.debug(f"    {facade} facade: {building_count} buildings in direction, max_obstruction={max_obstruction:.3f}")
        return max_obstruction
    
    def _calculate_neighbor_shadow_angles(
        self,
        plot_centroid: Tuple[float, float],
        plot_boundary: List[List[float]],
        neighbor_buildings: List[Dict[str, Any]]
    ) -> Dict[str, float]:
        """
        Calculate the sun angle at which each neighbor starts casting shadow on plot
        
        Higher angle = less shadow impact (only affects early morning/late evening)
        Lower angle = more shadow impact (affects more of the day)
        """
        ref_lon, ref_lat = plot_centroid
        local_boundary = GeometryUtils.degrees_to_local(plot_boundary, ref_lon, ref_lat)
        plot_centroid_local = GeometryUtils.polygon_centroid(local_boundary)
        
        shadow_angles = {}
        logger.debug(f"  Calculating shadow angles for {len(neighbor_buildings)} neighbor buildings")
        
        for building in neighbor_buildings:
            building_id = building.get("id", "unknown")
            building_coords = building.get("footprint", {}).get("coordinates", [[]])
            if not building_coords or not building_coords[0]:
                continue
            
            building_local = GeometryUtils.degrees_to_local(
                building_coords[0], ref_lon, ref_lat
            )
            
            # Find closest edge of building to plot
            min_distance = float('inf')
            for i in range(len(building_local) - 1):
                dist = GeometryUtils.distance_point_to_line(
                    plot_centroid_local,
                    building_local[i],
                    building_local[i+1]
                )
                min_distance = min(min_distance, dist)
            
            if min_distance > 0:
                height = building.get("height_m", 8.0)
                # Sun angle where shadow just reaches the plot
                shadow_angle = math.degrees(math.atan(height / min_distance))
                
                # Classify by direction - use wall_facing_plot if available (more accurate)
                # wall_facing_plot indicates which wall of the building faces the plot
                # If building's east wall faces plot, building is to the west of the plot
                wall_facing = building.get("wall_facing_plot", "")
                
                if wall_facing:
                    # Convert wall_facing to building direction (inverted mapping)
                    # If building's east wall faces plot → building is to the west → "west_building"
                    # If building's west wall faces plot → building is to the east → "east_building"
                    # If building's north wall faces plot → building is to the south → "south_building"
                    # If building's south wall faces plot → building is to the north → "north_building"
                    direction_map = {
                        "north": "south_building",  # North wall faces plot → building is to the south
                        "south": "north_building",  # South wall faces plot → building is to the north
                        "east": "west_building",   # East wall faces plot → building is to the west
                        "west": "east_building"    # West wall faces plot → building is to the east
                    }
                    direction = direction_map.get(wall_facing.lower(), None)
                    
                    if direction:
                        logger.debug(f"    Building {building_id}: using wall_facing_plot='{wall_facing}' → direction={direction}")
                        shadow_angles[direction] = round(shadow_angle, 0)
                        logger.debug(f"    Building {building_id}: min_distance={min_distance:.1f}m, height={height:.1f}m, "
                                    f"shadow_angle={shadow_angle:.1f}°, direction={direction}")
                    else:
                        # Skip building if wall_facing_plot is invalid (don't use fallback)
                        logger.warning(f"    Building {building_id}: wall_facing_plot='{wall_facing}' is invalid, skipping shadow angle calculation")
                else:
                    # Skip building if wall_facing_plot not available (don't use fallback)
                    logger.warning(f"    Building {building_id}: wall_facing_plot not available, skipping shadow angle calculation")
            else:
                logger.debug(f"    Building {building_id}: min_distance=0 (overlapping or invalid), skipping")
        
        logger.debug(f"  Shadow angles result: {shadow_angles}")
        return shadow_angles

