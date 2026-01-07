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
        neighbor_buildings: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        Perform comprehensive shadow analysis
        
        Returns facade scores, shadow hours, and neighbor shadow angles
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
        
        # Calculate facade scores
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
        base_scores = {
            "north": {"winter": 0.15, "summer": 0.45},
            "south": {"winter": 0.85, "summer": 0.95},
            "east": {"winter": 0.55, "summer": 0.70},
            "west": {"winter": 0.50, "summer": 0.65}
        }
        logger.debug(f"  Base facade scores: {base_scores}")
        
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
                logger.debug(f"    Building {building_id}: min_distance={min_distance:.1f}m, height={height:.1f}m, "
                            f"shadow_angle={shadow_angle:.1f}°, direction={direction}")
            else:
                logger.debug(f"    Building {building_id}: min_distance=0 (overlapping or invalid), skipping")
        
        logger.debug(f"  Shadow angles result: {shadow_angles}")
        return shadow_angles

