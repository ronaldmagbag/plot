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
    coordinates: List[List[float]]  # LineString coordinates [lon, lat]
    length_m: float
    color: str  # Color for visualization
    metadata: Dict[str, Any] = field(default_factory=dict)


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
                features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": segment.coordinates
                    },
                    "properties": {
                        "segment_type": segment.segment_type.value,
                        "length_m": segment.length_m,
                        "color": segment.color,
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

