"""
Configuration settings for Plot Analysis Generator
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
import os


@dataclass
class UKRegulatory:
    """Default UK residential building regulations"""
    # Standard setbacks for residential (meters)
    front_setback_m: float = 5.0
    rear_setback_m: float = 5.0
    side_setback_m: float = 1.0
    
    # Building constraints
    max_height_m: float = 10.0
    max_stories: int = 2
    max_coverage_ratio: float = 0.40
    
    # Zoning
    default_zoning: str = "R1"
    zoning_description: str = "Single family residential"
    permitted_uses: List[str] = field(default_factory=lambda: [
        "single_family_dwelling", 
        "home_office"
    ])
    conditional_uses: List[str] = field(default_factory=lambda: [
        "accessory_dwelling_unit"
    ])
    
    # Property line thresholds
    property_line_max_area_sqm: float = 2000.0  # If property > this, use OSM logic
    property_line_max_perimeter_m: float = 200.0  # If perimeter > this, use OSM logic
    property_line_min_side_distance_m: float = 3.0  # If side distance < this, consider merging


@dataclass 
class APIConfig:
    """API endpoints and configuration"""
    # Overpass API (OSM) - using main endpoint with better reliability
    # Options: overpass-api.de (main), lz4.overpass-api.de, z.overpass-api.de
    overpass_url: str = "https://overpass-api.de/api/interpreter"
    overpass_timeout: int = 90  # Increased timeout for batch queries
    
    # Open-Elevation API
    open_elevation_url: str = "https://api.open-elevation.com/api/v1/lookup"
    
    # SoilGrids API
    soilgrids_url: str = "https://rest.isric.org/soilgrids/v2.0/properties/query"
    
    # UK Land Registry INSPIRE
    uk_land_registry_wfs: str = "https://use-land-property-data.service.gov.uk/api/v1/datasets/inspire/download"
    
    # Nominatim (geocoding)
    nominatim_url: str = "https://nominatim.openstreetmap.org"
    
    # Mapbox settings
    mapbox_max_zoom: int = 19  # Maximum zoom level for Mapbox satellite imagery (typically 18-22)
    
    # Request settings
    request_timeout: int = 30
    max_retries: int = 3
    retry_delay: float = 1.0
    
    # User agent for API requests
    user_agent: str = "PlotAnalysisGenerator/1.0"


@dataclass
class PipelineConfig:
    """Pipeline configuration"""
    # Search radius around center point (meters)
    search_radius_m: float = 40.0
    context_radius_m: float = 50.0  # For roads and surrounding context (neighbors, roads, etc.)
    
    # Default country
    country: str = "GB"
    
    # Output settings
    output_dir: str = "output"
    
    # API config
    api: APIConfig = field(default_factory=APIConfig)
    
    # UK regulations
    uk_regulatory: UKRegulatory = field(default_factory=UKRegulatory)
    
    # Coordinate system
    source_crs: str = "EPSG:4326"  # WGS84 lat/lon
    local_crs: str = "EPSG:27700"  # British National Grid for UK
    
    # Road width estimates by type (meters)
    road_widths: Dict[str, float] = field(default_factory=lambda: {
        "motorway": 12.0,
        "trunk": 10.0,
        "primary": 8.0,
        "secondary": 7.0,
        "tertiary": 6.0,
        "unclassified": 5.0,
        "residential": 5.0,
        "living_street": 4.0,
        "service": 3.5,
        "footway": 2.0,
        "path": 1.5,
        "cycleway": 2.0,
    })
    
    # Building height estimates by type (meters)
    building_heights: Dict[str, float] = field(default_factory=lambda: {
        "house": 7.5,
        "detached": 8.0,
        "semi": 7.5,
        "semidetached_house": 7.5,
        "terrace": 7.0,
        "apartments": 12.0,
        "residential": 8.0,
        "commercial": 10.0,
        "retail": 5.0,
        "industrial": 8.0,
        "garage": 3.0,
        "shed": 2.5,
        "yes": 8.0,  # Generic building
    })


# Global config instance
config = PipelineConfig()


def get_config() -> PipelineConfig:
    """Get global configuration"""
    return config

