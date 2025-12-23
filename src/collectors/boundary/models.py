"""
Boundary data models

Data classes for property lines, setback lines, and buildable envelopes
"""

from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field
from enum import Enum


class SegmentType(Enum):
    """Property line segment types"""
    FRONT = "front"
    REAR = "rear"
    LEFT_SIDE = "left_side"
    RIGHT_SIDE = "right_side"
    UNKNOWN = "unknown"


@dataclass
class PropertyLineSegment:
    """A segment of a property line (front, rear, or side)"""
    segment_type: SegmentType
    edge_indices: List[int]  # Edge indices in the property line polygon
    length_m: float
    color: str  # Color for visualization
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def get_coordinates(self, property_coords: List[List[float]]) -> List[List[float]]:
        """
        Reconstruct coordinates from edge indices
        
        Returns the largest consecutive group only (for single LineString compatibility).
        Use get_coordinate_groups() to get all groups.
        
        Args:
            property_coords: Full property line coordinates
            
        Returns:
            List of coordinates for the largest consecutive group
        """
        groups = self.get_coordinate_groups(property_coords)
        if not groups:
            return []
        # Return largest group
        return max(groups, key=len) if len(groups) > 1 else groups[0]
    
    def get_coordinate_groups(self, property_coords: List[List[float]]) -> List[List[List[float]]]:
        """
        Reconstruct coordinates from edge indices, returning all consecutive groups
        
        Args:
            property_coords: Full property line coordinates
            
        Returns:
            List of coordinate groups (each group is a list of coordinates)
        """
        if not self.edge_indices:
            return []
        
        # Group consecutive edges
        sorted_indices = sorted(self.edge_indices)
        edge_groups = []
        current_group = [sorted_indices[0]]
        
        for i in range(1, len(sorted_indices)):
            prev_idx = sorted_indices[i - 1]
            curr_idx = sorted_indices[i]
            prev_end = (prev_idx + 1) % len(property_coords)
            
            if prev_end == curr_idx:
                current_group.append(curr_idx)
            else:
                edge_groups.append(current_group)
                current_group = [curr_idx]
        
        if current_group:
            edge_groups.append(current_group)
        
        # Build coordinates for each group
        coordinate_groups = []
        for group in edge_groups:
            segment_coords = []
            for i, idx in enumerate(group):
                start_point = property_coords[idx]
                end_idx = (idx + 1) % len(property_coords)
                end_point = property_coords[end_idx]
                
                if i == 0:
                    segment_coords.append(start_point)
                segment_coords.append(end_point)
            
            # Remove duplicates
            cleaned_coords = []
            for coord in segment_coords:
                if not cleaned_coords:
                    cleaned_coords.append(coord)
                else:
                    prev = cleaned_coords[-1]
                    if abs(prev[0] - coord[0]) > 1e-6 or abs(prev[1] - coord[1]) > 1e-6:
                        cleaned_coords.append(coord)
            
            # Ensure not closed
            if len(cleaned_coords) > 2:
                if (abs(cleaned_coords[0][0] - cleaned_coords[-1][0]) < 1e-6 and
                    abs(cleaned_coords[0][1] - cleaned_coords[-1][1]) < 1e-6):
                    cleaned_coords = cleaned_coords[:-1]
            
            if len(cleaned_coords) >= 2:
                coordinate_groups.append(cleaned_coords)
        
        return coordinate_groups


@dataclass
class PropertyLine:
    """
    Property boundary line
    
    Represents the actual property boundary with classification into segments
    """
    coordinates: List[List[float]]  # Polygon coordinates [lon, lat]
    area_sqm: float
    perimeter_m: float
    source: str  # e.g., "inspire_cadastral", "openstreetmap_landuse"
    accuracy_m: float
    inspire_id: Optional[str] = None
    
    # Classified segments
    front: Optional[PropertyLineSegment] = None
    rear: Optional[PropertyLineSegment] = None
    left_side: Optional[PropertyLineSegment] = None
    right_side: Optional[PropertyLineSegment] = None
    
    # Metadata
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_geojson(self) -> Dict[str, Any]:
        """Convert to GeoJSON format"""
        return {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [self.coordinates]
            },
            "properties": {
                "area_sqm": self.area_sqm,
                "perimeter_m": self.perimeter_m,
                "source": self.source,
                "accuracy_m": self.accuracy_m,
                "inspire_id": self.inspire_id,
                **self.metadata
            }
        }
    
    def get_segments_geojson(self) -> Dict[str, Any]:
        """Get all segments as GeoJSON FeatureCollection"""
        features = []
        
        for segment in [self.front, self.rear, self.left_side, self.right_side]:
            if segment:
                # Reconstruct coordinates from edge indices
                segment_coords = segment.get_coordinates(self.coordinates)
                features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": segment_coords
                    },
                    "properties": {
                        "segment_type": segment.segment_type.value,
                        "length_m": segment.length_m,
                        "color": segment.color,
                        "edge_indices": segment.edge_indices,
                        **segment.metadata
                    }
                })
        
        return {
            "type": "FeatureCollection",
            "features": features
        }


@dataclass
class SetbackLine:
    """
    Setback boundary line
    
    Represents the setback boundary (to be defined later)
    """
    coordinates: List[List[float]]  # Polygon coordinates
    area_sqm: float
    perimeter_m: float
    setback_type: str  # e.g., "front", "rear", "side"
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_geojson(self) -> Dict[str, Any]:
        """Convert to GeoJSON format"""
        return {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [self.coordinates]
            },
            "properties": {
                "area_sqm": self.area_sqm,
                "perimeter_m": self.perimeter_m,
                "setback_type": self.setback_type,
                **self.metadata
            }
        }


@dataclass
class BuildableEnvelope:
    """
    Buildable envelope
    
    Represents the buildable area within setbacks (to be defined later)
    """
    coordinates: List[List[float]]  # Polygon coordinates
    area_sqm: float
    perimeter_m: float
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_geojson(self) -> Dict[str, Any]:
        """Convert to GeoJSON format"""
        return {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [self.coordinates]
            },
            "properties": {
                "area_sqm": self.area_sqm,
                "perimeter_m": self.perimeter_m,
                **self.metadata
            }
        }

