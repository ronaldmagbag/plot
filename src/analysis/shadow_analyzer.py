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
        logger.info(f"Analyzing shadows for plot at ({lat}, {lon})")
        
        # Key dates for analysis
        analysis_dates = {
            "winter_solstice": datetime(2024, 12, 21),
            "summer_solstice": datetime(2024, 6, 21),
            "equinox": datetime(2024, 3, 20)
        }
        
        # Calculate sun positions for each date
        shadow_hours = {}
        for period, date in analysis_dates.items():
            shadow_hours[period] = self._calculate_sun_hours(lat, lon, date)
        
        # Calculate facade scores
        facade_scores = self._calculate_facade_scores(
            lat, lon, neighbor_buildings, plot_boundary
        )
        
        # Calculate neighbor shadow angles
        neighbor_angles = self._calculate_neighbor_shadow_angles(
            plot_centroid, plot_boundary, neighbor_buildings
        )
        
        # Determine best solar facade
        best_facade = max(
            facade_scores.keys(),
            key=lambda f: facade_scores[f]["annual_avg"]
        )
        
        return {
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
    
    def _calculate_sun_hours(
        self,
        lat: float,
        lon: float,
        date: datetime
    ) -> float:
        """Calculate approximate sun hours for a given date"""
        if _check_pvlib():
            return self._calculate_sun_hours_pvlib(lat, lon, date)
        else:
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
        
        # Get solar positions
        if not _check_pvlib() or _pvlib_module is None:
            return self._calculate_sun_hours_simple(lat, lon, date)
        solar_pos = _pvlib_module.solarposition.get_solarposition(times, lat, lon)
        
        # Count hours where sun is above horizon (elevation > 0)
        sun_above_horizon = solar_pos['apparent_elevation'] > 0
        sun_hours = sun_above_horizon.sum() * 0.25  # 15-minute intervals
        
        return round(sun_hours, 1)
    
    def _calculate_sun_hours_simple(self, lat: float, date: datetime) -> float:
        """Simplified sun hours calculation without pvlib"""
        # Day of year
        doy = date.timetuple().tm_yday
        
        # Solar declination angle (simplified)
        declination = 23.45 * math.sin(math.radians(360 * (284 + doy) / 365))
        
        # Convert to radians
        lat_rad = math.radians(lat)
        dec_rad = math.radians(declination)
        
        # Hour angle at sunrise/sunset
        cos_hour_angle = -math.tan(lat_rad) * math.tan(dec_rad)
        
        # Handle polar day/night
        if cos_hour_angle < -1:
            return 24.0  # Midnight sun
        elif cos_hour_angle > 1:
            return 0.0   # Polar night
        
        hour_angle = math.degrees(math.acos(cos_hour_angle))
        
        # Day length in hours
        day_length = 2 * hour_angle / 15
        
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
        
        # Get plot edges
        edges = GeometryUtils.get_polygon_edges(local_boundary)
        
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
        
        # Adjust for neighbor obstructions
        facade_scores = {}
        for facade, base in base_scores.items():
            # Check for neighbor buildings blocking this facade
            obstruction_factor = self._calculate_obstruction(
                facade, local_boundary, neighbor_buildings, ref_lon, ref_lat
            )
            
            winter_avg = base["winter"] * (1 - obstruction_factor * 0.5)
            summer_avg = base["summer"] * (1 - obstruction_factor * 0.3)
            annual_avg = (winter_avg + summer_avg) / 2
            
            facade_scores[facade] = {
                "winter_avg": round(winter_avg, 2),
                "summer_avg": round(summer_avg, 2),
                "annual_avg": round(annual_avg, 2)
            }
        
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
                
                if distance > 0:
                    # Calculate angular size of building
                    angle = math.degrees(math.atan(height / distance))
                    
                    # Closer and taller buildings cause more obstruction
                    # Maximum obstruction at ~25m distance for 10m building
                    obstruction = min(1.0, angle / 45.0)
                    max_obstruction = max(max_obstruction, obstruction)
        
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
                
                # Classify by direction
                building_centroid = GeometryUtils.polygon_centroid(building_local)
                dx = building_centroid[0] - plot_centroid_local[0]
                dy = building_centroid[1] - plot_centroid_local[1]
                
                if abs(dy) > abs(dx):
                    direction = "north_building" if dy > 0 else "south_building"
                else:
                    direction = "east_building" if dx > 0 else "west_building"
                
                shadow_angles[direction] = round(shadow_angle, 0)
        
        return shadow_angles

