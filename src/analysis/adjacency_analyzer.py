"""
Adjacency analysis - determines what each edge of the plot is adjacent to
"""

import math
from typing import List, Dict, Any, Tuple, Optional
from loguru import logger

from .geometry_utils import GeometryUtils


class AdjacencyAnalyzer:
    """
    Analyze what each edge of a plot boundary is adjacent to
    
    Categories: street, neighbor_parcel, water, public_land, agricultural, unknown
    """
    
    def __init__(self):
        # Threshold distances (meters)
        self.street_threshold = 15.0  # Distance to consider adjacent to street
        self.building_threshold = 30.0  # Distance to consider adjacent to building
        self.water_threshold = 20.0
    
    def analyze(
        self,
        plot_boundary: List[List[float]],
        roads: List[Dict[str, Any]],
        buildings: List[Dict[str, Any]],
        water_features: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Analyze adjacency for each edge of the plot boundary
        
        Returns list of edge analyses with adjacency information
        """
        if not plot_boundary or len(plot_boundary) < 3:
            logger.warning("Invalid plot boundary for adjacency analysis")
            return []
        
        # Convert to local coordinates
        ref_lon = plot_boundary[0][0]
        ref_lat = plot_boundary[0][1]
        local_boundary = GeometryUtils.degrees_to_local(plot_boundary, ref_lon, ref_lat)
        
        # Get edges
        edges = GeometryUtils.get_polygon_edges(local_boundary)
        
        # Prepare feature data in local coords
        local_roads = self._convert_roads_to_local(roads, ref_lon, ref_lat)
        local_buildings = self._convert_buildings_to_local(buildings, ref_lon, ref_lat)
        local_water = self._convert_water_to_local(water_features, ref_lon, ref_lat)
        
        # Analyze each edge
        adjacency_results = []
        for i, edge in enumerate(edges):
            edge_analysis = self._analyze_edge(
                edge, local_roads, local_buildings, local_water, ref_lon, ref_lat
            )
            
            # Convert geometry back to degrees
            edge_coords = GeometryUtils.local_to_degrees(
                [edge["start"], edge["end"]], ref_lon, ref_lat
            )
            
            edge_analysis["edge_id"] = self._get_edge_id(edge["direction"], i)
            edge_analysis["geometry"] = {
                "type": "LineString",
                "coordinates": edge_coords
            }
            edge_analysis["length_m"] = round(edge["length_m"], 1)
            
            adjacency_results.append(edge_analysis)
        
        logger.info(f"Analyzed {len(adjacency_results)} plot edges for adjacency")
        return adjacency_results
    
    def _analyze_edge(
        self,
        edge: Dict[str, Any],
        roads: List[Dict[str, Any]],
        buildings: List[Dict[str, Any]],
        water: List[Dict[str, Any]],
        ref_lon: float,
        ref_lat: float
    ) -> Dict[str, Any]:
        """Analyze a single edge for adjacency"""
        
        edge_midpoint = (
            (edge["start"][0] + edge["end"][0]) / 2,
            (edge["start"][1] + edge["end"][1]) / 2
        )
        
        # Check road adjacency
        road_info = self._check_road_adjacency(edge, roads)
        if road_info:
            return {
                "adjacent_to": "street",
                "street_name": road_info.get("name", "Unknown Road"),
                "street_type": road_info.get("type", "residential"),
                "primary_access": road_info.get("is_primary", True),
                "noise_level": self._estimate_noise_level(road_info.get("type")),
                "privacy_exposure": "high"  # Street-facing = public exposure
            }
        
        # Check water adjacency
        water_info = self._check_water_adjacency(edge, water)
        if water_info:
            return {
                "adjacent_to": "water_boundary",
                "water_type": water_info.get("type", "unknown"),
                "noise_level": "low",
                "privacy_exposure": "medium"
            }
        
        # Check building/neighbor adjacency
        building_info = self._check_building_adjacency(edge, buildings)
        if building_info:
            return {
                "adjacent_to": "neighbor_parcel",
                "neighbor_id": building_info.get("id"),
                "shared_wall_potential": building_info.get("distance", 100) < 2.0,
                "distance_to_neighbor_building_m": round(building_info.get("distance", 0), 1),
                "noise_level": "low",
                "privacy_exposure": "high" if building_info.get("distance", 100) < 10 else "medium"
            }
        
        # Default - adjacent to unknown parcel
        return {
            "adjacent_to": "neighbor_parcel",
            "neighbor_id": None,
            "shared_wall_potential": False,
            "distance_to_neighbor_building_m": None,
            "noise_level": "low",
            "privacy_exposure": "low"
        }
    
    def _check_road_adjacency(
        self,
        edge: Dict[str, Any],
        roads: List[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        """Check if edge is adjacent to a road"""
        edge_midpoint = (
            (edge["start"][0] + edge["end"][0]) / 2,
            (edge["start"][1] + edge["end"][1]) / 2
        )
        
        closest_road = None
        min_distance = float('inf')
        
        for road in roads:
            coords = road.get("local_coords", [])
            if len(coords) < 2:
                continue
            
            # Check distance from edge midpoint to road
            for i in range(len(coords) - 1):
                dist = GeometryUtils.distance_point_to_line(
                    edge_midpoint, coords[i], coords[i+1]
                )
                if dist < min_distance:
                    min_distance = dist
                    closest_road = road
        
        if closest_road and min_distance < self.street_threshold:
            return {
                "name": closest_road.get("name", "Unknown"),
                "type": closest_road.get("type", "residential"),
                "distance": min_distance,
                "is_primary": min_distance < 8.0  # Close enough for direct access
            }
        
        return None
    
    def _check_water_adjacency(
        self,
        edge: Dict[str, Any],
        water_features: List[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        """Check if edge is adjacent to water"""
        edge_midpoint = (
            (edge["start"][0] + edge["end"][0]) / 2,
            (edge["start"][1] + edge["end"][1]) / 2
        )
        
        for water in water_features:
            coords = water.get("local_coords", [])
            if len(coords) < 2:
                continue
            
            # Check distance
            for i in range(len(coords) - 1):
                dist = GeometryUtils.distance_point_to_line(
                    edge_midpoint, coords[i], coords[i+1]
                )
                if dist < self.water_threshold:
                    return {
                        "type": water.get("type", "water"),
                        "distance": dist
                    }
        
        return None
    
    def _check_building_adjacency(
        self,
        edge: Dict[str, Any],
        buildings: List[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        """Check if edge is adjacent to a building"""
        edge_midpoint = (
            (edge["start"][0] + edge["end"][0]) / 2,
            (edge["start"][1] + edge["end"][1]) / 2
        )
        
        closest_building = None
        min_distance = float('inf')
        
        for building in buildings:
            coords = building.get("local_coords", [])
            if len(coords) < 3:
                continue
            
            # Check distance to building polygon
            for i in range(len(coords)):
                j = (i + 1) % len(coords)
                dist = GeometryUtils.distance_point_to_line(
                    edge_midpoint, coords[i], coords[j]
                )
                if dist < min_distance:
                    min_distance = dist
                    closest_building = building
        
        if closest_building and min_distance < self.building_threshold:
            return {
                "id": closest_building.get("id"),
                "distance": min_distance,
                "type": closest_building.get("type")
            }
        
        return None
    
    def _convert_roads_to_local(
        self,
        roads: List[Dict[str, Any]],
        ref_lon: float,
        ref_lat: float
    ) -> List[Dict[str, Any]]:
        """Convert road coordinates to local meters"""
        local_roads = []
        for road in roads:
            centerline = road.get("centerline", {}).get("coordinates", [])
            if centerline:
                local_coords = GeometryUtils.degrees_to_local(centerline, ref_lon, ref_lat)
                local_roads.append({
                    **road,
                    "local_coords": local_coords
                })
        return local_roads
    
    def _convert_buildings_to_local(
        self,
        buildings: List[Dict[str, Any]],
        ref_lon: float,
        ref_lat: float
    ) -> List[Dict[str, Any]]:
        """Convert building coordinates to local meters"""
        local_buildings = []
        for building in buildings:
            footprint = building.get("footprint", {}).get("coordinates", [[]])
            if footprint and footprint[0]:
                local_coords = GeometryUtils.degrees_to_local(footprint[0], ref_lon, ref_lat)
                local_buildings.append({
                    **building,
                    "local_coords": local_coords
                })
        return local_buildings
    
    def _convert_water_to_local(
        self,
        water_features: List[Dict[str, Any]],
        ref_lon: float,
        ref_lat: float
    ) -> List[Dict[str, Any]]:
        """Convert water feature coordinates to local meters"""
        local_water = []
        for water in water_features:
            geom = water.get("geometry", {})
            coords = geom.get("coordinates", [])
            if geom.get("type") == "Polygon" and coords:
                coords = coords[0]  # Get outer ring
            if coords:
                local_coords = GeometryUtils.degrees_to_local(coords, ref_lon, ref_lat)
                local_water.append({
                    **water,
                    "local_coords": local_coords
                })
        return local_water
    
    def _get_edge_id(self, direction: str, index: int) -> str:
        """Generate edge ID from direction"""
        direction_map = {
            "north": "north_edge",
            "south": "south_edge",
            "east": "east_edge",
            "west": "west_edge",
            "northeast": "northeast_edge",
            "northwest": "northwest_edge",
            "southeast": "southeast_edge",
            "southwest": "southwest_edge"
        }
        base_id = direction_map.get(direction, f"edge_{index}")
        return f"{base_id}_{index}" if base_id in direction_map else base_id
    
    def _estimate_noise_level(self, road_type: str) -> str:
        """Estimate noise level based on road type"""
        high_noise = ["motorway", "trunk", "primary"]
        medium_noise = ["secondary", "tertiary"]
        
        if road_type in high_noise:
            return "high"
        elif road_type in medium_noise:
            return "medium"
        return "low"

