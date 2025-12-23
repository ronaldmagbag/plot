"""
Pydantic models for Plot Analysis data structure
Matches the exact JSON schema from documentation
"""

from datetime import datetime
from typing import List, Optional, Dict, Any, Literal
from pydantic import BaseModel, Field
import uuid


# ============================================================
# GeoJSON Types
# ============================================================

class GeoJSONPoint(BaseModel):
    type: Literal["Point"] = "Point"
    coordinates: List[float]  # [longitude, latitude]


class GeoJSONPolygon(BaseModel):
    type: Literal["Polygon"] = "Polygon"
    coordinates: List[List[List[float]]]  # [[[lon, lat], ...]]


class GeoJSONLineString(BaseModel):
    type: Literal["LineString"] = "LineString"
    coordinates: List[List[float]]  # [[lon, lat], ...]


# ============================================================
# Boundary Models
# ============================================================

class SetbacksApplied(BaseModel):
    front_m: float
    rear_m: float
    side_east_m: float
    side_west_m: float


class ConstraintApplied(BaseModel):
    type: str
    description: str
    area_removed_sqm: float
    polygon: Optional[GeoJSONPolygon] = None


class AccessCorridor(BaseModel):
    required: bool
    width_m: float
    description: str
    affects_buildable_area: bool


class PropertyLine(BaseModel):
    type: Literal["Polygon"] = "Polygon"
    coordinates: List[List[List[float]]]
    area_sqm: float
    perimeter_m: float
    source: str = "openstreetmap"
    source_date: str = Field(default_factory=lambda: datetime.utcnow().isoformat())
    accuracy_m: float = 2.0
    segments: Optional[Dict[str, Any]] = None  # Front, rear, left_side, right_side segments with colors


class SetbackLine(BaseModel):
    type: Literal["Polygon"] = "Polygon"
    coordinates: List[List[List[float]]]
    area_sqm: float
    derived_from: str = "property_line"
    setbacks_applied: SetbacksApplied
    regulation_source: str = "local_planning_code_R1"


class BuildableEnvelope(BaseModel):
    type: Literal["Polygon"] = "Polygon"
    coordinates: List[List[List[float]]]
    area_sqm: float
    derived_from: str = "setback_line"
    constraints_applied: List[ConstraintApplied] = Field(default_factory=list)
    access_corridor: Optional[AccessCorridor] = None


class Boundaries(BaseModel):
    property_line: PropertyLine
    setback_line: SetbackLine
    buildable_envelope: BuildableEnvelope


# ============================================================
# Surrounding Context Models
# ============================================================

class NeighborBuilding(BaseModel):
    id: str
    footprint: GeoJSONPolygon
    height_m: float
    stories: int = 2
    building_type: str
    usage: str = "residential"
    distance_to_property_line_m: float
    wall_facing_plot: str


class TreeZones(BaseModel):
    encoding: str = "rle"
    data: str = "0x00000000"
    resolution_m: float = 0.5
    bounds: GeoJSONPolygon
    coverage_percent: float = 0.0
    trees: List[Dict[str, Any]] = Field(default_factory=list)  # Individual trees


class Road(BaseModel):
    id: str
    name: str
    type: str
    centerline: GeoJSONLineString
    width_m: float
    traffic_level: str = "low"
    speed_limit_mph: Optional[int] = None
    has_sidewalk: bool = True


class ElevationSample(BaseModel):
    point: List[float]  # [x, y] in local coords
    elevation_m: float


class DEMReference(BaseModel):
    source: str
    resolution_m: float
    file_path: Optional[str] = None


class ElevationMap(BaseModel):
    corner_samples: List[ElevationSample]
    slope_percent: float
    slope_direction: str
    average_elevation_m: float
    terrain_classification: str
    dem_reference: Optional[DEMReference] = None


class NearestWaterBody(BaseModel):
    type: str
    distance_m: float
    direction: str


class WaterFeatures(BaseModel):
    encoding: str = "rle"
    data: str = "0x00000000"
    resolution_m: float = 1.0
    nearest_water_body: Optional[NearestWaterBody] = None
    features: List[Dict[str, Any]] = Field(default_factory=list)


class SurroundingContext(BaseModel):
    buildings: List[NeighborBuilding] = Field(default_factory=list)
    tree_zones: TreeZones
    roads: List[Road] = Field(default_factory=list)
    elevation_map: ElevationMap
    water_features: WaterFeatures


# ============================================================
# Analysis Models
# ============================================================

class FacadeScore(BaseModel):
    winter_avg: float
    summer_avg: float
    annual_avg: float


class ShadowHoursPerDay(BaseModel):
    winter_solstice: float
    summer_solstice: float
    equinox: float


class SunPathParams(BaseModel):
    latitude: float
    longitude: float
    timezone: str = "Europe/London"


class ShadowAnalysis(BaseModel):
    facade_scores: Dict[str, FacadeScore]
    shadow_hours_per_day: ShadowHoursPerDay
    neighbor_shadow_angles: Dict[str, float]
    best_solar_facade: str
    computed_at: str = Field(default_factory=lambda: datetime.utcnow().isoformat())
    sun_path_params: SunPathParams


class AdjacencyEdge(BaseModel):
    edge_id: str
    geometry: GeoJSONLineString
    length_m: float
    adjacent_to: str
    street_name: Optional[str] = None
    street_type: Optional[str] = None
    neighbor_id: Optional[str] = None
    primary_access: bool = False
    shared_wall_potential: bool = False
    distance_to_neighbor_building_m: Optional[float] = None
    noise_level: str = "low"
    privacy_exposure: str = "low"


class Analysis(BaseModel):
    shadow_analysis: ShadowAnalysis
    adjacency: List[AdjacencyEdge] = Field(default_factory=list)


# ============================================================
# Access Models
# ============================================================

class PrimaryAccessPoint(BaseModel):
    location: GeoJSONPoint
    side: str
    confidence: str = "high"
    determined_by: str = "street_adjacency"
    alternatives: List[Dict[str, Any]] = Field(default_factory=list)


class VehicleAccess(BaseModel):
    available: bool
    type: str = "direct_from_street"
    driveway_width_required_m: float = 3.5
    turning_radius_adequate: bool = True


class PedestrianAccess(BaseModel):
    available: bool = True
    sidewalk_present: bool = True
    grade_accessible: bool = True


class Access(BaseModel):
    primary_access_point: PrimaryAccessPoint
    vehicle_access: VehicleAccess
    pedestrian_access: PedestrianAccess


# ============================================================
# Existing Structures
# ============================================================

class ExistingStructure(BaseModel):
    id: str
    footprint: GeoJSONPolygon
    type: str
    condition: str = "unknown"
    stories: int = 2
    height_m: float
    year_built: Optional[int] = None
    status: str = "unknown"  # to_demolish, to_preserve, to_renovate, unknown
    demolition_cost_estimate_gbp: Optional[float] = None
    affects_buildable_area: bool = True


# ============================================================
# Regulatory Models
# ============================================================

class Zoning(BaseModel):
    designation: str
    description: str
    permitted_uses: List[str]
    conditional_uses: List[str] = Field(default_factory=list)
    source: str = "local_planning_authority"
    last_verified: str = Field(default_factory=lambda: datetime.utcnow().isoformat())


class SetbacksM(BaseModel):
    front: float
    rear: float
    side_east: float
    side_west: float
    notes: Optional[str] = None


class BuildingConstraints(BaseModel):
    max_coverage_ratio: float
    max_coverage_sqm: float
    max_height_m: float
    max_stories: int
    max_stories_note: Optional[str] = None
    setbacks_m: SetbacksM


class AdditionalRestriction(BaseModel):
    type: str
    status: str
    notes: Optional[str] = None


class Regulatory(BaseModel):
    zoning: Zoning
    building_constraints: BuildingConstraints
    additional_restrictions: List[AdditionalRestriction] = Field(default_factory=list)


# ============================================================
# Soil Model
# ============================================================

class Soil(BaseModel):
    type_id: str
    bearing_capacity_kpa: float
    drainage: str
    foundation_recommendation: str
    source: str = "soilgrids"


# ============================================================
# Data Quality Model
# ============================================================

class DataSource(BaseModel):
    type: str
    source: str
    date: Optional[str] = None
    resolution_m: Optional[float] = None
    accuracy_m: Optional[float] = None
    verified: Optional[str] = None


class DataQuality(BaseModel):
    data_sources: List[DataSource]
    missing_data: List[str] = Field(default_factory=list)
    last_validated: str = Field(default_factory=lambda: datetime.utcnow().isoformat())


# ============================================================
# Main Plot Analysis Model
# ============================================================

PlotType = Literal[
    "residential_building",
    "building_massive",
    "field",
    "industrial_building",
    "commercial",
    "civic_infrastructure",
    "unknown"
]


class PlotAnalysis(BaseModel):
    """Complete Plot Analysis document matching the target JSON schema"""
    
    # Core Identification
    plot_id: str = Field(default_factory=lambda: f"UKS-GB-{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())[:4]}")
    created_at: str = Field(default_factory=lambda: datetime.utcnow().isoformat())
    updated_at: str = Field(default_factory=lambda: datetime.utcnow().isoformat())
    data_version: str = "1.0"
    
    # Location & Geometry
    centroid: GeoJSONPoint
    plot_type: PlotType = "residential_building"
    
    # Boundaries
    boundaries: Boundaries
    
    # Surrounding Context
    surrounding_context: SurroundingContext
    
    # Analysis
    analysis: Analysis
    
    # Access & Existing Conditions
    access: Access
    existing_structures: List[ExistingStructure] = Field(default_factory=list)
    
    # Regulatory Constraints
    regulatory: Regulatory
    
    # Soil
    soil: Soil
    
    # Data Quality
    data_quality: DataQuality

    class Config:
        json_encoders = {
            datetime: lambda v: v.isoformat()
        }

