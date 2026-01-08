"""
Main Pipeline Orchestrator for Plot Analysis Generation

Implements the 11-Step Automated Flow from Pipeline v2.0:

  1. Input: Lat/Lon or Parcel ID
  2. Fetch imagery + DEM (Mapbox Terrain-RGB)
  3. Fetch parcels, buildings, roads (OSM + Cadastre)
  4. Segment trees & water (SAM or landcover fallback)
  5. Compute terrain & slope
  6. Compute shadows (pvlib)
  7. Apply zoning + setbacks (rule engine)
  8. Detect access (road adjacency)
  9. Assemble plot.json

Data Sources:
  - OpenStreetMap (Overpass API): Buildings, roads, water, vegetation
  - Mapbox Terrain-RGB: Elevation/DEM
  - Open-Meteo: Elevation fallback
  - SoilGrids/BGS: Soil classification
  - pvlib: Solar/shadow calculations
"""

import json
import math
import os
import subprocess
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from loguru import logger

from .config import get_config, PipelineConfig
from .models import (
    PlotAnalysis, GeoJSONPoint, GeoJSONPolygon, GeoJSONLineString,
    Boundaries, PropertyLine, SetbackLine, BuildableEnvelope, SetbacksApplied,
    SurroundingContext, NeighborBuilding, TreeZones, Road, ElevationMap, 
    ElevationSample, DEMReference, WaterFeatures, NearestWaterBody,
    Analysis, ShadowAnalysis, FacadeScore, ShadowHoursPerDay, SunPathParams, 
    AdjacencyEdge, StreetAdjacencyEdge, BuildingAdjacencyEdge, WaterAdjacencyEdge,
    Access, PrimaryAccessPoint, VehicleAccess, PedestrianAccess,
    ExistingStructure, Regulatory, Zoning, BuildingConstraints, SetbacksM,
    Soil, DataQuality, DataSource, AccessCorridor, ConstraintApplied
)

from .collectors import (
    OSMCollector, ElevationCollector, SoilCollector, BoundaryCollector,
    TerrainCollector, DSMCollector, VegetationCollector, MapboxImageryCollector
)
from .collectors.boundary.utils import calculate_polygon_area
from .analysis import ShadowAnalyzer, AdjacencyAnalyzer, SetbackCalculator, GeometryUtils


class PlotAnalysisPipeline:
    """
    Main pipeline to generate complete plot analysis from center coordinates
    
    Usage:
        pipeline = PlotAnalysisPipeline()
        # radius_m must be provided or set in config.search_radius_m
        result = pipeline.run(lat=51.5074, lon=-0.1276, radius_m=50.0)
        pipeline.save(result, "output/plot.json")
    """
    
    def __init__(self, config: Optional[PipelineConfig] = None, cache_dir: Optional[str] = None):
        self.config = config or get_config()
        self.cache_dir = cache_dir  # Will be set based on output filename
        
        # Initialize collectors (Layer 2: Data Acquisition)
        # Pass cache directory to collectors that support it
        self.osm_collector = OSMCollector(cache_dir=cache_dir)
        self.elevation_collector = ElevationCollector()
        self.terrain_collector = TerrainCollector()
        self.dsm_collector = DSMCollector()
        self.vegetation_collector = VegetationCollector()
        self.soil_collector = SoilCollector()
        self.boundary_collector = BoundaryCollector()
        self.mapbox_imagery_collector = MapboxImageryCollector(cache_dir=cache_dir)
        
        # Connect terrain and DSM collectors to OSM collector for building height calculation
        self.osm_collector.set_terrain_collector(self.terrain_collector)
        self.osm_collector.set_dsm_collector(self.dsm_collector)
        
        # Initialize analyzers (Layers 3-8: Processing & Analysis)
        self.shadow_analyzer = ShadowAnalyzer()
        self.adjacency_analyzer = AdjacencyAnalyzer()
        self.setback_calculator = SetbackCalculator()
    
    def run(
        self,
        lat: float,
        lon: float,
        radius_m: Optional[float] = None,
        plot_id: Optional[str] = None
    ) -> PlotAnalysis:
        """
        Run the complete pipeline to generate plot analysis
        
        Args:
            lat: Latitude of plot center
            lon: Longitude of plot center
            radius_m: Search radius in meters (required - must be provided or set in config)
            plot_id: Optional custom plot ID
        
        Returns:
            Complete PlotAnalysis object
            
        Raises:
            ValueError: If required configuration values are missing
        """
        # Validate configuration first
        from src.config import validate_config
        validate_config(self.config)
        
        # Determine radius: CLI argument takes precedence, then config, otherwise error
        if radius_m is not None:
            if radius_m <= 0:
                raise ValueError(f"radius_m must be positive, got {radius_m}")
            radius = radius_m
        elif hasattr(self.config, 'search_radius_m') and self.config.search_radius_m is not None:
            if self.config.search_radius_m <= 0:
                raise ValueError(f"config.search_radius_m must be positive, got {self.config.search_radius_m}")
            radius = self.config.search_radius_m
        else:
            raise ValueError(
                "radius_m is required but not provided. "
                "Either pass radius_m parameter or set search_radius_m in config."
            )
        
        logger.info(f"Starting plot analysis pipeline for ({lat}, {lon})")
        logger.info(f"Search radius: {radius}m")
        
        # Generate plot ID
        if not plot_id:
            plot_id = f"UKS-{self.config.country}-{datetime.now().strftime('%Y%m%d%H%M%S')}-{str(uuid.uuid4())[:4]}"
        
        timestamp = datetime.utcnow().isoformat()
        
        # ============================================================
        # STAGE 1: Collect Boundary Data (use OSM data if available)
        # ============================================================
        logger.info("Stage 1: Collecting boundary data...")
        
        # First, try to get OSM data to use building footprint for boundary
        # This avoids a separate API call and uses the same data
        all_osm_features_preview = self.osm_collector.fetch_all_features(lat, lon, radius)
        osm_buildings_preview = all_osm_features_preview.get("buildings", [])
        
        # Get roads for boundary derivation
        osm_roads_preview = all_osm_features_preview.get("roads", [])
        
        # PRIORITY: Try INSPIRE/cadastral data FIRST (most accurate)
        # This ensures we use the most accurate property boundary available
        boundary_data = None
        property_line_obj_stage1 = None
        
        # Try boundary collector with OSM data (roads and buildings) for better boundary detection
        # This allows deriving boundaries from roads and neighbor buildings, not just rectangles
        # Get property line object (will be reused later for classification to avoid duplicate detection)
        # Get property line without classification first (for initial boundary)
        # This will be reused in STAGE 3 with classification to avoid duplicate INSPIRE queries
        property_line_obj_stage1 = self.boundary_collector.get_property_line(
            lat, lon, radius,
            osm_buildings=osm_buildings_preview,
            osm_roads=osm_roads_preview,
            classify=False  # Don't classify yet, will do it later with full data
        )
        
        if property_line_obj_stage1:
            # Convert to old format for compatibility
            boundary_data = {
                "type": "Polygon",
                "coordinates": [property_line_obj_stage1.coordinates],
                "area_sqm": property_line_obj_stage1.area_sqm,
                "perimeter_m": property_line_obj_stage1.perimeter_m,
                "source": property_line_obj_stage1.source,
                "accuracy_m": property_line_obj_stage1.accuracy_m,
                "inspire_id": property_line_obj_stage1.inspire_id
            }
            logger.info(f"Using property line from {property_line_obj_stage1.source} (area: {property_line_obj_stage1.area_sqm:.1f} m²)")
        else:
            # If no boundary found, use default plot estimate
            logger.info("Property line detection returned None, using default rectangular plot estimate")
            boundary_data = self.boundary_collector._create_default_plot(lat, lon)
        
        coords_raw = boundary_data.get("coordinates", [])
        source_name = boundary_data.get("source", "unknown")
        
        logger.debug(f"Extracting coordinates from source '{source_name}': type={type(coords_raw)}, length={len(coords_raw) if isinstance(coords_raw, list) else 'N/A'}")
        
        if not coords_raw:
            logger.error(f"No coordinates in boundary_data from {source_name}")
            boundary_data = self.boundary_collector._create_default_plot(lat, lon)
            coords_raw = boundary_data.get("coordinates", [])
            source_name = boundary_data.get("source", "unknown")
        
        property_coords = []
        
        if isinstance(coords_raw, list) and len(coords_raw) > 0:
            first_element = coords_raw[0]
            
            if isinstance(first_element, list) and len(first_element) > 0:
                second_element = first_element[0] if len(first_element) > 0 else None
                
                if isinstance(second_element, list) and len(second_element) > 0:
                    if len(second_element) > 0 and isinstance(second_element[0], list) and len(second_element[0]) == 2:
                        property_coords = second_element
                        logger.debug(f"Extracted coordinates from double-wrapped format: {len(property_coords)} points")
                    elif len(first_element) > 0 and isinstance(first_element[0], list) and len(first_element[0]) == 2:
                        property_coords = first_element
                        logger.debug(f"Extracted coordinates from correct nested format: {len(property_coords)} points")
                    else:
                        logger.error(f"Unable to determine coordinate structure")
                        logger.error(f"  first_element type: {type(first_element)}, length: {len(first_element)}")
                        logger.error(f"  second_element type: {type(second_element)}, length: {len(second_element) if isinstance(second_element, list) else 'N/A'}")
                        if isinstance(second_element, list) and len(second_element) > 0:
                            logger.error(f"  second_element[0] type: {type(second_element[0])}")
                        property_coords = []
                elif isinstance(second_element, (int, float)):
                    logger.error(f"Unexpected: second_element is a number: {second_element}")
                    property_coords = []
                elif len(first_element) > 0 and isinstance(first_element[0], list) and len(first_element[0]) == 2:
                    property_coords = first_element
                    logger.debug(f"Extracted coordinates from correct format: {len(property_coords)} points")
                else:
                    logger.error(f"Unexpected structure")
                    logger.error(f"  first_element type: {type(first_element)}, length: {len(first_element)}")
                    if len(first_element) > 0:
                        logger.error(f"  first_element[0] type: {type(first_element[0])}")
                    property_coords = []
            else:
                logger.error(f"First element is not a valid list: type={type(first_element)}, length={len(first_element) if isinstance(first_element, list) else 'N/A'}")
                property_coords = []
        else:
            logger.error(f"Invalid coordinates structure: type={type(coords_raw)}, length={len(coords_raw) if isinstance(coords_raw, list) else 'N/A'}")
            property_coords = []
        
        if property_coords:
            if not isinstance(property_coords, list):
                logger.error(f"property_coords is not a list: {type(property_coords)}")
                property_coords = []
            elif len(property_coords) > 0:
                first_coord = property_coords[0]
                if not isinstance(first_coord, list) or len(first_coord) != 2:
                    logger.error(f"First coordinate is invalid: {first_coord} (type: {type(first_coord)})")
                    logger.error(f"This suggests coordinate extraction failed. property_coords: {property_coords[:3]}")
                    property_coords = []
        
        # Validate property coordinates
        if not property_coords or len(property_coords) < 3:
            logger.error(f"Invalid property coordinates from {source_name}: {len(property_coords) if property_coords else 0} points")
            logger.error(f"coords_raw type: {type(coords_raw)}, length: {len(coords_raw) if isinstance(coords_raw, list) else 'N/A'}")
            if isinstance(coords_raw, list) and len(coords_raw) > 0:
                logger.error(f"coords_raw[0] type: {type(coords_raw[0])}, length: {len(coords_raw[0]) if isinstance(coords_raw[0], list) else 'N/A'}")
                if isinstance(coords_raw[0], list) and len(coords_raw[0]) > 0:
                    logger.error(f"coords_raw[0][0] type: {type(coords_raw[0][0])}, value: {coords_raw[0][0]}")
            logger.error(f"Full coords_raw structure: {coords_raw}")
            # Fallback to default plot
            logger.info("Falling back to default rectangular plot estimate")
            boundary_data = self.boundary_collector._create_default_plot(lat, lon)
            coords_raw = boundary_data.get("coordinates", [])
            if isinstance(coords_raw, list) and len(coords_raw) > 0:
                if isinstance(coords_raw[0], list):
                    property_coords = coords_raw[0]
                else:
                    property_coords = []
            else:
                property_coords = []
        
        import copy
        property_coords_original = copy.deepcopy(property_coords) if property_coords else []
        
        property_area = boundary_data.get("area_sqm", 300)
        property_perimeter = boundary_data.get("perimeter_m", 70)
        
        # ============================================================
        # STAGE 2: Collect Surrounding Context (BATCH QUERY)
        # ============================================================
        logger.info("Stage 2: Collecting surrounding context (batch OSM query)...")
        
        # Reuse OSM data if we already have it, otherwise fetch fresh
        if all_osm_features_preview and len(osm_buildings_preview) > 0:
            # Use cached/preview data
            all_osm_features = all_osm_features_preview
        else:
            # Single batch query for all OSM features (reduces API calls from 4+ to 1)
            all_osm_features = self.osm_collector.fetch_all_features(lat, lon, radius)
        
        osm_buildings = all_osm_features["buildings"]
        osm_roads = all_osm_features["roads"]
        vegetation_data = {
            "trees": all_osm_features["trees"],
            "zones": all_osm_features["vegetation_zones"],
            "coverage_percent": 0.0
        }
        water_features_raw = all_osm_features["water"]
        
        # Elevation
        elevation_data = self.elevation_collector.get_elevation_map(
            lat, lon, 
            width_m=max(30, property_perimeter / 4),
            height_m=max(30, property_perimeter / 4)
        )
        
        # Soil
        soil_data = self.soil_collector.get_soil_data(lat, lon)
        
        # Mapbox imagery download (always download if cache_dir available, regardless of SAM3 flag)
        merged_image = None
        imagery_metadata = None
        if self.cache_dir:
            # Check if merged image already exists in cache
            zoom = self.config.api.mapbox_max_zoom
            cached_image_pattern = f"mapbox_merged_{lat:.6f}_{lon:.6f}_{radius}m_z{zoom}.jpg"
            cached_images = list(Path(self.cache_dir).glob(cached_image_pattern))
            
            if cached_images:
                logger.info(f"Found cached merged image: {cached_images[0].name}, skipping download")
                try:
                    from PIL import Image
                    merged_image = Image.open(cached_images[0])
                    # Reconstruct metadata from filename and config
                    imagery_metadata = {
                        "bounds": {
                            "west": lon - radius / 111000 / abs(math.cos(math.radians(lat))),
                            "east": lon + radius / 111000 / abs(math.cos(math.radians(lat))),
                            "south": lat - radius / 111000,
                            "north": lat + radius / 111000
                        },
                        "zoom": zoom,
                        "center": [lon, lat],
                        "radius_m": radius
                    }
                except Exception as e:
                    logger.warning(f"Failed to load cached image: {e}, will re-download")
                    cached_images = []
            
            if not cached_images:
                logger.info("Downloading Mapbox satellite imagery...")
                try:
                    merged_image, imagery_metadata = self.mapbox_imagery_collector.download_imagery(
                        lat, lon, radius_m=radius
                    )
                    if merged_image:
                        logger.info(f"Successfully downloaded Mapbox imagery: {merged_image.size}")
                except Exception as e:
                    logger.warning(f"Mapbox imagery download failed: {e}")
                    merged_image = None
                    imagery_metadata = None
        
        # SAM3 segmentation (only if enabled)
        sam3_results = None
        if self.config.enable_sam3:
            if self.cache_dir:
                # Check if SAM3 results already exist
                sam3_output_path = Path(self.cache_dir) / f"sam3_results_{lat:.6f}_{lon:.6f}.json"
                if sam3_output_path.exists():
                    logger.info(f"Loading existing SAM3 results from {sam3_output_path.name}...")
                    try:
                        with open(sam3_output_path, "r") as f:
                            sam3_results = json.load(f)
                        logger.info(f"SAM3 found: {len(sam3_results.get('roads', []))} roads, "
                                  f"{len(sam3_results.get('trees', []))} trees, "
                                  f"{len(sam3_results.get('grasses', []))} grasses")
                    except Exception as e:
                        logger.warning(f"Failed to load existing SAM3 results: {e}, will regenerate")
                        sam3_results = None
                
                # If SAM3 results don't exist and we have imagery, run segmentation
                if sam3_results is None and merged_image and imagery_metadata:
                    logger.info("Running SAM3 segmentation...")
                    sam3_results = self._run_sam3_segmentation(
                        merged_image, imagery_metadata, lat, lon, radius_m=radius
                    )
                    if sam3_results:
                        logger.info(f"SAM3 found: {len(sam3_results.get('roads', []))} roads, "
                                  f"{len(sam3_results.get('trees', []))} trees, "
                                  f"{len(sam3_results.get('grasses', []))} grasses")
                    else:
                        logger.warning("SAM3 segmentation completed but returned no results")
                elif not merged_image:
                    logger.warning("No imagery available for SAM3 segmentation, skipping")
            else:
                logger.info("SAM3 segmentation requires cache_dir, skipping")
        else:
            logger.info("SAM3 segmentation is disabled in config (enable_sam3=False), skipping")
        
        # ============================================================
        # STAGE 3: Analysis
        # ============================================================
        logger.info("Stage 3: Running analysis...")
        
        # Validate property coordinates
        if not property_coords or len(property_coords) < 3:
            logger.error(f"Invalid property coordinates: {len(property_coords) if property_coords else 0} points")
            raise ValueError(f"Property coordinates must have at least 3 points, got {len(property_coords) if property_coords else 0}")
        
        # Calculate TRUE centroid of property polygon
        try:
            property_centroid = GeometryUtils.polygon_centroid(property_coords)
            if not isinstance(property_centroid, (tuple, list)) or len(property_centroid) != 2:
                logger.error(f"polygon_centroid returned invalid result: {property_centroid} (type: {type(property_centroid)})")
                raise ValueError(f"Invalid centroid from polygon_centroid: {property_centroid}")
            centroid_lon, centroid_lat = property_centroid
        except (ZeroDivisionError, ValueError, TypeError) as e:
            logger.error(f"Failed to calculate property centroid: {e}")
            logger.error(f"Property coords: {property_coords[:5] if len(property_coords) > 5 else property_coords}")
            raise ValueError(f"Failed to calculate property centroid: {e}") from e
        
        # Filter buildings: only keep those OUTSIDE the property boundary
        # Safely extract first coordinate
        if not property_coords or len(property_coords) == 0:
            logger.error("property_coords is empty, cannot extract reference point")
            raise ValueError("Property coordinates are empty")
        
        first_coord = property_coords[0]
        if not isinstance(first_coord, list) or len(first_coord) != 2:
            logger.error(f"First coordinate is invalid: {first_coord} (type: {type(first_coord)})")
            logger.error(f"property_coords structure: {property_coords[:3] if len(property_coords) > 3 else property_coords}")
            raise ValueError(f"Invalid first coordinate: {first_coord}")
        
        ref_lon, ref_lat = first_coord
        local_property = GeometryUtils.degrees_to_local(property_coords, ref_lon, ref_lat)
        
        neighbor_buildings = []
        existing_buildings = []
        
        for building in osm_buildings:
            footprint = building.get("footprint", {})
            b_coords = footprint.get("coordinates", [[]])
            if not b_coords or not b_coords[0]:
                continue
            
            b_centroid = GeometryUtils.polygon_centroid(b_coords[0])
            local_b_centroid = GeometryUtils.degrees_to_local([list(b_centroid)], ref_lon, ref_lat)[0]
            
            if self._point_in_polygon(local_b_centroid, local_property):
                # Building is INSIDE property - it's an existing structure
                existing_buildings.append(building)
            else:
                # Building is OUTSIDE property - it's a neighbor
                neighbor_buildings.append(building)
        
        logger.info(f"Separated {len(existing_buildings)} existing + {len(neighbor_buildings)} neighbor buildings")
        
        # Create PropertyLine object from original coordinates
        # Simplification will be handled by boundary collector
        from .collectors.boundary.models import PropertyLine as BoundaryPropertyLine
        # PropertyLine.coordinates expects List[List[float]] = [[lon, lat], ...] (NOT GeoJSON Polygon format)
        # property_coords is already [[lon, lat], ...], so use it directly
        logger.info(f"Creating PropertyLine object with {len(property_coords)} coordinates")
        property_line_obj_original = BoundaryPropertyLine(
            coordinates=property_coords,  # List[List[float]] = [[lon, lat], ...]
            area_sqm=property_area,
            perimeter_m=property_perimeter,
            source=boundary_data.get("source", "estimated"),
            accuracy_m=boundary_data.get("accuracy_m", 5.0)
        )
        
        logger.info("Simplifying property line using boundary collector...")
        property_line_obj = self.boundary_collector.simplify_property_line(property_line_obj_original)
        property_coords = property_line_obj.coordinates
        
        preliminary_adjacency = self._quick_adjacency(property_coords, osm_roads)
        
        logger.info(f"Classifying simplified property line with full OSM data... (coordinates: {len(property_coords)} points)")
        property_line_obj = self.boundary_collector.classifier.classify(
            property_line_obj,
            osm_roads,
            osm_buildings,
            debug=False
        )
        
        if property_line_obj:
            setback_line_obj = self.boundary_collector.get_setback_line(property_line_obj)
            buildable_envelope_obj = self.boundary_collector.get_buildable_envelope(
                property_line_obj, setback_line_obj
            )
            
            if setback_line_obj:
                # Get actual applied setbacks from metadata (includes dynamic rear setback)
                metadata = setback_line_obj.metadata
                front_setback = metadata.get("front_setback_m", self.config.uk_regulatory.front_setback_m)
                rear_setback = metadata.get("rear_setback_m", self.config.uk_regulatory.rear_setback_m)
                side_setback = metadata.get("side_setback_m", self.config.uk_regulatory.side_setback_m)
                
                logger.info(f"Setbacks applied - Front: {front_setback}m, Rear: {rear_setback}m, Side: {side_setback}m")
                
                setback_result = {
                    "coordinates": [setback_line_obj.coordinates],
                    "area_sqm": round(setback_line_obj.area_sqm, 1),
                    "setbacks_applied": {
                        "front_m": front_setback,
                        "rear_m": rear_setback,
                        "side_east_m": side_setback,
                        "side_west_m": side_setback
                    },
                    "regulation_source": "uk_planning_guidance"
                }
            else:
                setback_result = None
            
            if buildable_envelope_obj:
                buildable_result = {
                    "coordinates": [buildable_envelope_obj.coordinates],
                    "area_sqm": round(buildable_envelope_obj.area_sqm, 1),
                    "constraints_applied": []
                }
                logger.info(f"Buildable envelope calculated: {buildable_envelope_obj.area_sqm:.2f} m²")
            else:
                logger.warning("Buildable envelope calculation returned None - will fall back to setback line")
                buildable_result = None
        else:
            setback_result = self.setback_calculator.calculate_setbacks(
                property_coords, preliminary_adjacency
            )
            setback_coords = setback_result["coordinates"][0] if setback_result else property_coords
            
            buildable_result = self.setback_calculator.calculate_buildable_envelope(
                setback_coords, 
                constraints=[],
                property_area_sqm=property_area,
                neighbor_buildings=neighbor_buildings
            )
            # Note: access_corridor will be updated after adjacency_result is calculated
        
        # Use neighbor_buildings (excludes existing structures inside property)
        # Existing structures should not be considered as neighbors
        adjacency_result = self.adjacency_analyzer.analyze(
            property_coords, osm_roads, neighbor_buildings, water_features_raw
        )
        
        # Calculate access corridor after adjacency analysis is available
        if buildable_result:
            buildable_area = buildable_result.get("area_sqm", 0)
            access_corridor_data = self.setback_calculator._determine_access_corridor(
                property_area, buildable_area, adjacency_result, property_line_obj, property_coords
            )
            buildable_result["access_corridor"] = access_corridor_data
        else:
            # Fallback: calculate even if buildable_result is None
            access_corridor_data = self.setback_calculator._determine_access_corridor(
                property_area, 0, adjacency_result, property_line_obj, property_coords
            )
            buildable_result = {"access_corridor": access_corridor_data}
        
        # Calculate wall_facing_plot for neighbor buildings before shadow analysis
        # This is needed for accurate direction detection in shadow angles
        ref_lon, ref_lat = property_coords[0] if property_coords else (lon, lat)
        local_property = GeometryUtils.degrees_to_local(property_coords, ref_lon, ref_lat)
        
        for building in neighbor_buildings:
            footprint = building.get("footprint", {})
            if not footprint.get("coordinates"):
                continue
            
            b_coords = footprint["coordinates"][0]
            local_b = GeometryUtils.degrees_to_local(b_coords, ref_lon, ref_lat)
            
            # Find closest points between building and property polygons
            closest_building_pt, closest_property_pt, _ = GeometryUtils.closest_points_between_polygons(
                local_b, local_property
            )
            
            # Calculate direction of the distance line (from building to property)
            dx = closest_property_pt[0] - closest_building_pt[0]
            dy = closest_property_pt[1] - closest_building_pt[1]
            
            # Calculate angle of the distance line
            if abs(dx) < 1e-10 and abs(dy) < 1e-10:
                # Fallback to centroid-based
                b_centroid = GeometryUtils.polygon_centroid(b_coords)
                p_centroid = GeometryUtils.polygon_centroid(property_coords)
                local_b_centroid = GeometryUtils.degrees_to_local([list(b_centroid)], ref_lon, ref_lat)[0]
                local_p_centroid = GeometryUtils.degrees_to_local([list(p_centroid)], ref_lon, ref_lat)[0]
                dx = local_p_centroid[0] - local_b_centroid[0]
                dy = local_p_centroid[1] - local_b_centroid[1]
            
            angle = math.degrees(math.atan2(dx, dy))
            if angle < 0:
                angle += 360
            
            # Convert angle to cardinal direction
            if angle < 22.5 or angle >= 337.5:
                wall_facing = "north"
            elif angle < 67.5:
                wall_facing = "north"
            elif angle < 112.5:
                wall_facing = "east"
            elif angle < 157.5:
                wall_facing = "east"
            elif angle < 202.5:
                wall_facing = "south"
            elif angle < 247.5:
                wall_facing = "south"
            elif angle < 292.5:
                wall_facing = "west"
            else:
                wall_facing = "west"
            
            # Add wall_facing_plot to building dict for shadow analyzer
            building["wall_facing_plot"] = wall_facing
        
        # Use neighbor_buildings (excludes existing structures inside property)
        # Now includes wall_facing_plot for each building
        shadow_result = self.shadow_analyzer.analyze(
            (lon, lat), property_coords, neighbor_buildings, property_line_obj=property_line_obj
        )
        
        logger.info("Stage 4: Building output model...")
        
        centroid = GeoJSONPoint(coordinates=[centroid_lon, centroid_lat])
        
        boundaries = self._build_boundaries(
            property_coords, property_area, property_perimeter,
            setback_result, buildable_result, boundary_data,
            property_line_obj,  # Pass property line object for segments
            property_coords_original  # Pass original coordinates for property_line
        )
        
        # Surrounding context (use filtered neighbor_buildings, not all osm_buildings)
        surrounding_context = self._build_surrounding_context(
            property_coords, centroid_lat, centroid_lon,
            neighbor_buildings, osm_roads, vegetation_data,
            water_features_raw, elevation_data, sam3_results
        )
        
        # Analysis
        analysis = self._build_analysis(shadow_result, adjacency_result)
        
        # Access
        access = self._build_access(property_coords, osm_roads, adjacency_result, property_line_obj)
        
        # Existing structures (use pre-filtered existing_buildings)
        existing_structures = self._build_existing_structures(existing_buildings)
        
        # Regulatory
        regulatory = self._build_regulatory(property_area)
        
        # Soil
        soil = Soil(
            type_id=soil_data.get("type_id", "unknown"),
            bearing_capacity_kpa=soil_data.get("bearing_capacity_kpa", 150),
            drainage=soil_data.get("drainage", "moderate"),
            foundation_recommendation=soil_data.get("foundation_recommendation", "standard_strip"),
            source=soil_data.get("source", "soilgrids"),
            properties=soil_data.get("properties", {})  # Include raw soil properties
        )
        
        # Data quality
        data_quality = self._build_data_quality(
            boundary_data=boundary_data,
            soil_data=soil_data,
            elevation_data=elevation_data
        )
        
        plot_analysis = PlotAnalysis(
            plot_id=plot_id,
            created_at=timestamp,
            updated_at=timestamp,
            data_version="1.0",
            centroid=centroid,
            plot_type="residential_building",
            boundaries=boundaries,
            surrounding_context=surrounding_context,
            analysis=analysis,
            access=access,
            existing_structures=existing_structures,
            regulatory=regulatory,
            soil=soil,
            data_quality=data_quality
        )
        
        logger.info(f"Pipeline complete. Plot ID: {plot_id}")
        return plot_analysis
    
    def save(self, plot_analysis: PlotAnalysis, output_path: str) -> str:
        """Save plot analysis to JSON file"""
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(plot_analysis.model_dump(mode='python'), f, indent=2, ensure_ascii=False)
        
        logger.info(f"Saved plot analysis to {output_path}")
        return output_path
    
    def _run_sam3_segmentation(
        self,
        merged_image,
        imagery_metadata: Dict[str, Any],
        lat: float,
        lon: float,
        radius_m: float = None
    ) -> Optional[Dict[str, Any]]:
        """
        Run SAM3 segmentation on downloaded imagery.
        
        This method saves the image temporarily and calls the SAM3 segmentation
        script in the cu126 environment.
        
        Args:
            merged_image: PIL Image of merged Mapbox tiles
            imagery_metadata: Metadata from imagery download
            lat: Center latitude
            lon: Center longitude
        
        Returns:
            Dictionary with segmented features (roads, trees, grasses) or None if failed
        """
        if not self.cache_dir:
            return None
        
        try:
            # Reuse the already-cached merged image instead of saving again
            # The image was already saved by MapboxImageryCollector._save_to_cache()
            zoom = imagery_metadata.get("zoom", 20)
            radius = radius_m or imagery_metadata.get("radius_m", self.config.search_radius_m)
            
            # Find the cached merged image
            cached_image_pattern = f"mapbox_merged_{lat:.6f}_{lon:.6f}_{radius}m_z{zoom}.jpg"
            cached_images = list(Path(self.cache_dir).glob(cached_image_pattern))
            
            if cached_images:
                image_path = cached_images[0]
                logger.debug(f"Reusing cached merged image: {image_path.name}")
            else:
                # Fallback: save if cache file not found (shouldn't happen)
                image_path = Path(self.cache_dir) / f"mapbox_for_sam3_{lat:.6f}_{lon:.6f}.jpg"
                merged_image.save(image_path, quality=95)
                logger.debug(f"Saved image for SAM3: {image_path.name}")
            
            # Prepare output path
            output_path = Path(self.cache_dir) / f"sam3_results_{lat:.6f}_{lon:.6f}.json"
            
            # Prepare bbox JSON
            bbox = imagery_metadata.get("bounds", {})
            bbox_json = json.dumps(bbox)
            
            # Find SAM3 script path
            # __file__ is at src/pipeline.py, so go to src/segmentation/
            script_path = Path(__file__).parent / "segmentation" / "sam3_segment.py"
            if not script_path.exists():
                logger.warning(f"SAM3 script not found at {script_path}")
                logger.debug(f"Looking for script relative to: {Path(__file__).parent.parent}")
                return None
            
            # Run SAM3 segmentation in cu126 environment
            # Use micromamba to activate cu126 and run the script
            logger.info("Running SAM3 segmentation in cu126 environment...")
            
            # Use full path to python in cu126 environment to ensure correct environment
            # First try to find micromamba python
            import shutil
            micromamba_cmd = shutil.which("micromamba")
            if not micromamba_cmd:
                logger.warning("micromamba not found in PATH, trying direct command...")
                micromamba_cmd = "micromamba"
            
            cmd = [
                micromamba_cmd, "run", "-n", "cu126",
                "python", str(script_path),
                "--image", str(image_path),
                "--output", str(output_path),
                "--bbox", bbox_json,
                "--center-lon", str(lon),
                "--center-lat", str(lat)
            ]
            
            logger.debug(f"Running command: {' '.join(cmd[:4])} ... (script with args)")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                logger.warning(f"SAM3 segmentation failed: {result.stderr}")
                return None
            
            if output_path.exists():
                with open(output_path, "r") as f:
                    results = json.load(f)
                logger.info(f"SAM3 segmentation completed successfully")
                return results
            else:
                logger.warning("SAM3 output file not found")
                return None
                
        except subprocess.TimeoutExpired:
            logger.warning("SAM3 segmentation timed out")
            return None
        except Exception as e:
            logger.warning(f"SAM3 segmentation error: {e}")
            return None
    
    def _quick_adjacency(
        self,
        property_coords: List[List[float]],
        roads: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Quick adjacency check for setback determination"""
        if not property_coords or len(property_coords) < 4:
            return []
        
        ref_lon, ref_lat = property_coords[0]
        local_boundary = GeometryUtils.degrees_to_local(property_coords, ref_lon, ref_lat)
        edges = GeometryUtils.get_polygon_edges(local_boundary)
        
        result = []
        for edge in edges:
            result.append({
                "edge_id": edge["direction"],
                "adjacent_to": "neighbor_parcel",
                "shared_wall_potential": False
            })
        
        for road in roads:
            road_coords = road.get("centerline", {}).get("coordinates", [])
            if not road_coords:
                continue
            
            road_local = GeometryUtils.degrees_to_local(road_coords, ref_lon, ref_lat)
            
            for i, edge in enumerate(edges):
                edge_mid = ((edge["start"][0] + edge["end"][0]) / 2,
                           (edge["start"][1] + edge["end"][1]) / 2)
                
                for j in range(len(road_local) - 1):
                    dist = GeometryUtils.distance_point_to_line(
                        edge_mid, road_local[j], road_local[j+1]
                    )
                    if dist < 15:  # 15m threshold
                        result[i]["adjacent_to"] = "street"
                        break
        
        return result
    
    def _build_boundaries(
        self,
        property_coords: List[List[float]],
        property_area: float,
        property_perimeter: float,
        setback_result: Dict[str, Any],
        buildable_result: Dict[str, Any],
        boundary_data: Dict[str, Any],
        property_line_obj: Optional[Any] = None,
        property_coords_original: Optional[List[List[float]]] = None
    ) -> Boundaries:
        """Build Boundaries model from collected data
        
        Note: property_coords is the SIMPLIFIED version (used for all downstream operations)
        property_coords_original is the ORIGINAL version (stored as property_line)
        """
        
        # Property line - store original coordinates and simplified coordinates
        property_line_original_coords = property_coords_original if property_coords_original else property_coords
        logger.info(f"Property line: {len(property_line_original_coords)} original points, {len(property_coords)} simplified points")
        
        segments_data = None
        if property_line_obj and hasattr(property_line_obj, 'front'):
            segments_data = {}
            
            if property_line_obj.front:
                segments_data["front"] = {
                    "edges": property_line_obj.front.edge_indices,
                    "length_m": property_line_obj.front.length_m,
                    "color": property_line_obj.front.color
                }
            if property_line_obj.rear:
                segments_data["rear"] = {
                    "edges": property_line_obj.rear.edge_indices,
                    "length_m": property_line_obj.rear.length_m,
                    "color": property_line_obj.rear.color
                }
            if property_line_obj.left_side:
                segments_data["left_side"] = {
                    "edges": property_line_obj.left_side.edge_indices,
                    "length_m": property_line_obj.left_side.length_m,
                    "color": property_line_obj.left_side.color
                }
            if property_line_obj.right_side:
                segments_data["right_side"] = {
                    "edges": property_line_obj.right_side.edge_indices,
                    "length_m": property_line_obj.right_side.length_m,
                    "color": property_line_obj.right_side.color
                }
        
        property_line = PropertyLine(
            coordinates=[property_line_original_coords],
            coordinates_simplified=[property_coords],
            area_sqm=round(property_area, 1),
            perimeter_m=round(property_perimeter, 1),
            source=boundary_data.get("source", "estimated"),
            source_date=datetime.utcnow().isoformat(),
            accuracy_m=boundary_data.get("accuracy_m", 5.0),
            segments=segments_data
        )
        
        if setback_result and setback_result.get("coordinates"):
            setback_coords_raw = setback_result["coordinates"][0]
        else:
            # property_coords is already [[lon, lat], [lon, lat], ...] format
            setback_coords_raw = property_coords
        
        setbacks_applied = setback_result.get("setbacks_applied", {}) if setback_result else {}
        
        # Calculate actual area from coordinates (don't use arbitrary fallback percentages)
        if setback_result and setback_result.get("area_sqm"):
            setback_area = setback_result["area_sqm"]
        else:
            # Calculate from coordinates if available
            try:
                # setback_coords_raw should be [[lon, lat], [lon, lat], ...]
                if setback_coords_raw and len(setback_coords_raw) > 0:
                    # Check if it's already in the right format
                    if isinstance(setback_coords_raw[0], list) and len(setback_coords_raw[0]) >= 2:
                        # Calculate average latitude from coordinates for area calculation
                        avg_lat = sum(coord[1] for coord in setback_coords_raw) / len(setback_coords_raw)
                        setback_area = calculate_polygon_area(setback_coords_raw, avg_lat)
                    else:
                        # Unexpected format
                        logger.error(f"Unexpected setback coordinates format: {type(setback_coords_raw[0])}")
                        setback_area = 0.0
                else:
                    setback_area = 0.0
                if setback_result is None:
                    logger.warning(f"Setback calculation failed - calculated area from fallback coordinates: {setback_area:.2f} m²")
            except Exception as e:
                logger.error(f"Failed to calculate setback area from coordinates: {e}")
                setback_area = 0.0  # Mark as invalid/unknown
        
        # Ensure coordinates are in GeoJSON Polygon format: [[[lon, lat], ...]]
        # setback_coords_raw should be [[lon, lat], [lon, lat], ...] (list of coordinate pairs)
        try:
            if not setback_coords_raw or len(setback_coords_raw) == 0:
                logger.warning("No setback coordinates available, using property_coords")
                setback_coords_formatted = [property_coords]
            else:
                first_elem = setback_coords_raw[0]
                logger.debug(f"Setback coords format check: first element type={type(first_elem)}, value={first_elem if isinstance(first_elem, (int, float)) else 'list'}")
                
                if isinstance(first_elem, list):
                    if len(first_elem) > 0 and isinstance(first_elem[0], list):
                        setback_coords_formatted = setback_coords_raw
                        logger.debug("Setback coords already in GeoJSON format")
                    elif len(first_elem) >= 2 and isinstance(first_elem[0], (int, float)):
                        setback_coords_formatted = [setback_coords_raw]
                        logger.debug(f"Setback coords wrapped: {len(setback_coords_raw)} pairs")
                    else:
                        logger.error(f"Unexpected nested format in setback_coords_raw[0]: {type(first_elem[0]) if first_elem else 'empty'}")
                        setback_coords_formatted = [property_coords]
                elif isinstance(first_elem, (int, float)):
                    if len(setback_coords_raw) % 2 != 0:
                        logger.error(f"Flat coordinate list has odd length: {len(setback_coords_raw)}")
                        setback_coords_formatted = [property_coords]
                    else:
                        coords_pairs = [[setback_coords_raw[i], setback_coords_raw[i+1]] 
                                       for i in range(0, len(setback_coords_raw), 2)]
                        setback_coords_formatted = [coords_pairs]
                        logger.debug(f"Setback coords converted from flat list: {len(coords_pairs)} pairs")
                else:
                    logger.error(f"Unexpected setback coordinates format: first element is {type(first_elem)}")
                    setback_coords_formatted = [property_coords]
        except Exception as e:
            logger.error(f"Error formatting setback coordinates: {e}, using property_coords")
            import traceback
            logger.debug(traceback.format_exc())
            setback_coords_formatted = [property_coords]
        
        setback_line = SetbackLine(
            coordinates=setback_coords_formatted,
            area_sqm=round(setback_area, 1),
            setbacks_applied=SetbacksApplied(
                front_m=setbacks_applied.get("front_m", self.config.uk_regulatory.front_setback_m),
                rear_m=setbacks_applied.get("rear_m", self.config.uk_regulatory.rear_setback_m),
                side_east_m=setbacks_applied.get("side_east_m", self.config.uk_regulatory.side_setback_m),
                side_west_m=setbacks_applied.get("side_west_m", self.config.uk_regulatory.side_setback_m)
            ),
            regulation_source=setback_result.get("regulation_source", "uk_planning_guidance") if setback_result else "uk_planning_guidance"
        )
        
        # Buildable envelope
        if buildable_result and buildable_result.get("coordinates"):
            buildable_coords_raw = buildable_result["coordinates"][0]  # Extract first ring from GeoJSON Polygon
            logger.info(f"Using calculated buildable envelope: {buildable_result.get('area_sqm', 0):.2f} m²")
        else:
            buildable_coords_raw = setback_coords_raw  # Use same format as setback
            logger.warning("Buildable envelope not available - using setback line coordinates as fallback")
        
        constraints = buildable_result.get("constraints_applied", []) if buildable_result else []
        access_corridor_data = buildable_result.get("access_corridor", {}) if buildable_result else {}
        
        access_corridor = None
        if access_corridor_data and access_corridor_data.get("required"):
            access_corridor = AccessCorridor(
                required=True,
                width_m=access_corridor_data.get("width_m", 3.5),
                description=access_corridor_data.get("description", "Vehicle access corridor"),
                affects_buildable_area=access_corridor_data.get("affects_buildable_area", False)
            )
        
        # Calculate actual area from coordinates (don't use arbitrary fallback percentages)
        if buildable_result and buildable_result.get("area_sqm"):
            buildable_area = buildable_result["area_sqm"]
        else:
            # Calculate from coordinates if available
            try:
                # buildable_coords_raw should be [[lon, lat], [lon, lat], ...]
                if buildable_coords_raw and len(buildable_coords_raw) > 0:
                    # Check if it's already in the right format
                    if isinstance(buildable_coords_raw[0], list) and len(buildable_coords_raw[0]) >= 2:
                        # Calculate average latitude from coordinates for area calculation
                        avg_lat = sum(coord[1] for coord in buildable_coords_raw) / len(buildable_coords_raw)
                        buildable_area = calculate_polygon_area(buildable_coords_raw, avg_lat)
                    else:
                        # Unexpected format
                        logger.error(f"Unexpected buildable coordinates format: {type(buildable_coords_raw[0])}")
                        buildable_area = 0.0
                else:
                    buildable_area = 0.0
                if buildable_result is None:
                    logger.warning(f"Buildable envelope calculation failed - calculated area from fallback coordinates: {buildable_area:.2f} m²")
            except Exception as e:
                logger.error(f"Failed to calculate buildable area from coordinates: {e}")
                buildable_area = 0.0  # Mark as invalid/unknown
        
        # Ensure coordinates are in GeoJSON Polygon format: [[[lon, lat], ...]]
        # buildable_coords_raw should be [[lon, lat], [lon, lat], ...] (list of coordinate pairs)
        try:
            if not buildable_coords_raw or len(buildable_coords_raw) == 0:
                logger.warning("No buildable coordinates available, using setback_coords")
                buildable_coords_formatted = [setback_coords_raw] if setback_coords_raw else [property_coords]
            elif isinstance(buildable_coords_raw[0], list):
                # Check if it's nested lists (GeoJSON format) or list of pairs
                if len(buildable_coords_raw[0]) > 0 and isinstance(buildable_coords_raw[0][0], list):
                    # Already in GeoJSON format: [[[lon, lat], ...]]
                    buildable_coords_formatted = buildable_coords_raw
                elif len(buildable_coords_raw[0]) >= 2 and isinstance(buildable_coords_raw[0][0], (int, float)):
                    # List of pairs: [[lon, lat], [lon, lat], ...] - wrap in GeoJSON Polygon format
                    buildable_coords_formatted = [buildable_coords_raw]
                else:
                    logger.error(f"Unexpected nested format in buildable_coords_raw[0]: {type(buildable_coords_raw[0][0])}")
                    buildable_coords_formatted = [setback_coords_raw] if setback_coords_raw else [property_coords]
            elif isinstance(buildable_coords_raw[0], (int, float)):
                # Flat list: [lon, lat, lon, lat, ...] - convert to pairs then wrap
                if len(buildable_coords_raw) % 2 != 0:
                    logger.error(f"Flat coordinate list has odd length: {len(buildable_coords_raw)}")
                    buildable_coords_formatted = [setback_coords_raw] if setback_coords_raw else [property_coords]
                else:
                    coords_pairs = [[buildable_coords_raw[i], buildable_coords_raw[i+1]] 
                                   for i in range(0, len(buildable_coords_raw), 2)]
                    buildable_coords_formatted = [coords_pairs]
            else:
                logger.error(f"Unexpected buildable coordinates format: first element is {type(buildable_coords_raw[0])}")
                buildable_coords_formatted = [setback_coords_raw] if setback_coords_raw else [property_coords]
        except Exception as e:
            logger.error(f"Error formatting buildable coordinates: {e}, using setback_coords")
            buildable_coords_formatted = [setback_coords_raw] if setback_coords_raw else [property_coords]
        
        buildable_envelope = BuildableEnvelope(
            coordinates=buildable_coords_formatted,
            area_sqm=round(buildable_area, 1),
            constraints_applied=[
                ConstraintApplied(
                    type=c.get("type", "unknown"),
                    description=c.get("description", ""),
                    area_removed_sqm=c.get("area_removed_sqm", 0)
                ) for c in constraints
            ],
            access_corridor=access_corridor
        )
        
        return Boundaries(
            property_line=property_line,  # Contains original coordinates, simplified coordinates, and segments
            setback_line=setback_line,
            buildable_envelope=buildable_envelope
        )
    
    def _build_surrounding_context(
        self,
        property_coords: List[List[float]],
        lat: float,
        lon: float,
        osm_buildings: List[Dict[str, Any]],
        osm_roads: List[Dict[str, Any]],
        vegetation_data: Dict[str, Any],
        water_features_raw: List[Dict[str, Any]],
        elevation_data: Dict[str, Any],
        sam3_results: Optional[Dict[str, Any]] = None
    ) -> SurroundingContext:
        """Build SurroundingContext model"""
        
        # Neighbor buildings - calculate distances FIRST, then sort (no limit - include all that intersect 60m box)
        ref_lon, ref_lat = property_coords[0] if property_coords else (lon, lat)
        local_property = GeometryUtils.degrees_to_local(property_coords, ref_lon, ref_lat)
        
        # Calculate distance for each building to property LINE (not centroid)
        buildings_with_distances = []
        for b in osm_buildings:
            footprint = b.get("footprint", {})
            if not footprint.get("coordinates"):
                continue
            
            b_coords = footprint["coordinates"][0]
            local_b = GeometryUtils.degrees_to_local(b_coords, ref_lon, ref_lat)
            
            # Calculate MINIMUM distance from building polygon to property polygon
            distance = GeometryUtils.distance_polygon_to_polygon(local_b, local_property)
            
            # Determine wall facing plot by finding closest building edge to plot
            # Algorithm:
            # 1. Get all building edges with their directions (relative to building center)
            # 2. Find the edge with minimum distance to plot polygon
            # 3. The edge direction indicates which wall of the building faces the plot
            #
            # Edge direction angle ranges (from building center to edge midpoint):
            # - North: 337.5°-22.5° (or -22.5° to 22.5°)
            # - East: 67.5°-112.5°
            # - South: 157.5°-202.5°
            # - West: 247.5°-292.5°
            # Diagonal directions (northeast, southeast, etc.) are handled by primary component
            
            building_center = GeometryUtils.polygon_centroid(local_b)
            building_edges = GeometryUtils.get_polygon_edges(local_b, center=building_center)
            
            # Find the closest building edge to the plot polygon
            min_edge_distance = float('inf')
            closest_edge = None
            
            for edge in building_edges:
                # Calculate minimum distance from this edge to plot polygon
                edge_dist = GeometryUtils.distance_line_to_polygon(
                    edge["start"], edge["end"], local_property
                )
                if edge_dist < min_edge_distance:
                    min_edge_distance = edge_dist
                    closest_edge = edge
            
            # Determine wall_facing_plot using minimum distance line from building to property line
            # Algorithm: Find closest points between building polygon and property polygon,
            # then calculate direction of the distance line
            closest_building_pt, closest_property_pt, min_line_distance = GeometryUtils.closest_points_between_polygons(
                local_b, local_property
            )
            
            # Calculate direction of the distance line (from building to property)
            # This tells us which wall of the building faces the plot
            dx = closest_property_pt[0] - closest_building_pt[0]  # Property X - Building X
            dy = closest_property_pt[1] - closest_building_pt[1]  # Property Y - Building Y
            
            # Calculate angle of the distance line
            if abs(dx) < 1e-10 and abs(dy) < 1e-10:
                # Points are the same (shouldn't happen, but handle it)
                # Fallback to centroid-based
                b_centroid = GeometryUtils.polygon_centroid(b_coords)
                p_centroid = GeometryUtils.polygon_centroid(property_coords)
                local_b_centroid = GeometryUtils.degrees_to_local([list(b_centroid)], ref_lon, ref_lat)[0]
                local_p_centroid = GeometryUtils.degrees_to_local([list(p_centroid)], ref_lon, ref_lat)[0]
                dx = local_p_centroid[0] - local_b_centroid[0]
                dy = local_p_centroid[1] - local_b_centroid[1]
            
            angle = math.degrees(math.atan2(dx, dy))
            if angle < 0:
                angle += 360
            
            # Convert angle to cardinal direction
            # This tells us which wall of the building faces the plot
            if angle < 22.5 or angle >= 337.5:
                wall_facing = "north"  # Distance line points north → building's north wall faces plot
            elif angle < 67.5:
                wall_facing = "north"  # Northeast → use primary: north
            elif angle < 112.5:
                wall_facing = "east"  # Distance line points east → building's east wall faces plot
            elif angle < 157.5:
                wall_facing = "east"  # Southeast → use primary: east
            elif angle < 202.5:
                wall_facing = "south"  # Distance line points south → building's south wall faces plot
            elif angle < 247.5:
                wall_facing = "south"  # Southwest → use primary: south
            elif angle < 292.5:
                wall_facing = "west"  # Distance line points west → building's west wall faces plot
            else:
                wall_facing = "west"  # Northwest → use primary: west
            
            logger.debug(f"  Building {b.get('id', 'unknown')}: closest_building_pt=({closest_building_pt[0]:.1f}, {closest_building_pt[1]:.1f}), "
                        f"closest_property_pt=({closest_property_pt[0]:.1f}, {closest_property_pt[1]:.1f}), "
                        f"dx={dx:.1f}, dy={dy:.1f}, angle={angle:.1f}° → wall_facing={wall_facing}")
            
            buildings_with_distances.append({
                "building": b,
                "distance": distance,
                "wall_facing": wall_facing
            })
        
        # Sort by distance (closest first)
        # Include ALL buildings that intersect the 60m context box (no limit)
        buildings_with_distances.sort(key=lambda x: x["distance"])
        # Filter by spatial intersection with 60m box instead of limiting by count
        context_box_m = 60.0
        filtered_buildings = []
        for item in buildings_with_distances:
            b = item["building"]
            footprint = b.get("footprint", {})
            b_coords = footprint.get("coordinates", [[]])[0]
            if b_coords:
                local_b = GeometryUtils.degrees_to_local(b_coords, ref_lon, ref_lat)
                # Check if building intersects 60m box
                # A polygon intersects if any point is inside OR if any edge crosses the boundary
                intersects = False
                for pt in local_b:
                    if -context_box_m <= pt[0] <= context_box_m and -context_box_m <= pt[1] <= context_box_m:
                        intersects = True
                        break
                if not intersects:
                    # Check if any edge crosses the boundary
                    for i in range(len(local_b)):
                        p1 = local_b[i]
                        p2 = local_b[(i + 1) % len(local_b)]
                        # Simple bounding box check
                        min_x, max_x = min(p1[0], p2[0]), max(p1[0], p2[0])
                        min_y, max_y = min(p1[1], p2[1]), max(p1[1], p2[1])
                        if not (max_x < -context_box_m or min_x > context_box_m or 
                               max_y < -context_box_m or min_y > context_box_m):
                            intersects = True
                            break
                if intersects:
                    filtered_buildings.append(item)
        
        buildings_with_distances = filtered_buildings
        
        # Build NeighborBuilding objects
        buildings = []
        for item in buildings_with_distances:
            b = item["building"]
            buildings.append(NeighborBuilding(
                id=b.get("id", f"building_{len(buildings)}"),
                footprint=GeoJSONPolygon(coordinates=b.get("footprint", {}).get("coordinates", [])),
                height_m=b.get("height_m", 8.0),
                stories=b.get("stories", 2),
                building_type=b.get("building_type", "residential"),
                usage=b.get("usage", "residential"),
                distance_to_property_line_m=round(item["distance"], 1),
                wall_facing_plot=item["wall_facing"]
            ))
        
        # Roads - include ALL roads that intersect the 60m context box
        # Don't filter by count, filter by spatial intersection
        roads = []
        ref_lon, ref_lat = property_coords[0] if property_coords else (lon, lat)
        context_box_m = 60.0  # 60m box for filtering
        
        for r in osm_roads:
            centerline = r.get("centerline", {})
            if not centerline.get("coordinates"):
                continue
            
            # Check if road intersects the 60m context box
            road_coords = centerline["coordinates"]
            road_intersects = False
            
            # Convert to local meters to check intersection
            local_road = GeometryUtils.degrees_to_local(road_coords, ref_lon, ref_lat)
            
            # Check if any point is within 60m box OR if any segment crosses the box
            for pt in local_road:
                if -context_box_m <= pt[0] <= context_box_m and -context_box_m <= pt[1] <= context_box_m:
                    road_intersects = True
                    break
            
            # Also check if segments cross the boundary
            if not road_intersects:
                for i in range(len(local_road) - 1):
                    p1, p2 = local_road[i], local_road[i + 1]
                    # Simple check: if segment bounding box overlaps with context box
                    min_x, max_x = min(p1[0], p2[0]), max(p1[0], p2[0])
                    min_y, max_y = min(p1[1], p2[1]), max(p1[1], p2[1])
                    if not (max_x < -context_box_m or min_x > context_box_m or 
                           max_y < -context_box_m or min_y > context_box_m):
                        road_intersects = True
                        break
            
            # Only include roads that intersect the context box
            if road_intersects:
                roads.append(Road(
                    id=r.get("id", f"road_{len(roads)}"),
                    name=r.get("name", "Unknown Road"),
                    type=r.get("type", "residential"),
                    centerline=GeoJSONLineString(coordinates=centerline["coordinates"]),
                    width_m=r.get("width_m", 5.0),
                    traffic_level=r.get("traffic_level", "low"),
                    speed_limit_mph=r.get("speed_limit_mph"),
                    has_sidewalk=r.get("has_sidewalk", True)
                ))
        
        # Add SAM3 segmented roads (not in OSM)
        if sam3_results:
            sam3_roads = sam3_results.get("roads", [])
            osm_road_ids = {r.get("id", "") for r in osm_roads}
            
            for seg_road in sam3_roads:
                # Check if this road is already in OSM (simple check by proximity)
                road_geom = seg_road.get("geometry", {})
                if not road_geom or not road_geom.get("coordinates"):
                    continue
                
                # Convert to local meters to check if it's near any OSM road
                seg_coords = road_geom["coordinates"][0]  # Polygon coordinates
                is_duplicate = False
                
                for osm_road in osm_roads:
                    osm_centerline = osm_road.get("centerline", {})
                    if not osm_centerline.get("coordinates"):
                        continue
                    
                    # Check if segmented road is close to OSM road
                    local_seg = GeometryUtils.degrees_to_local(seg_coords, ref_lon, ref_lat)
                    local_osm = GeometryUtils.degrees_to_local(osm_centerline["coordinates"], ref_lon, ref_lat)
                    
                    # Simple distance check: if any point is within 5m, consider it duplicate
                    for seg_pt in local_seg:
                        for osm_pt in local_osm:
                            dist = math.sqrt((seg_pt[0] - osm_pt[0])**2 + (seg_pt[1] - osm_pt[1])**2)
                            if dist < 5.0:  # 5m threshold
                                is_duplicate = True
                                break
                        if is_duplicate:
                            break
                    if is_duplicate:
                        break
                
                if not is_duplicate:
                    # Extract centerline from polygon (use first ring as approximation)
                    centerline_coords = seg_coords
                    if len(centerline_coords) > 2:
                        # Simplify polygon to centerline (use first and last points for now)
                        centerline_coords = [centerline_coords[0], centerline_coords[-2]]
                    
                    roads.append(Road(
                        id=seg_road.get("id", f"road_sam3_{len(roads)}"),
                        name="Segmented Road",
                        type="residential",
                        centerline=GeoJSONLineString(coordinates=centerline_coords),
                        width_m=4.0,  # Default width for segmented roads
                        traffic_level="low",
                        has_sidewalk=False
                    ))
                    logger.debug(f"Added SAM3 segmented road: {seg_road.get('id')}")
        
        # Tree zones
        trees = vegetation_data.get("trees", [])
        zones = vegetation_data.get("zones", [])
        
        # Add SAM3 segmented trees
        if sam3_results:
            sam3_trees = sam3_results.get("trees", [])
            for seg_tree in sam3_trees:
                trees.append({
                    "id": seg_tree.get("id", f"tree_sam3_{len(trees)}"),
                    "geometry": seg_tree.get("geometry", {}),
                    "confidence": seg_tree.get("confidence", 0.5),
                    "source": "sam3_segmentation"
                })
            
            # Add SAM3 grasses as vegetation zones
            sam3_grasses = sam3_results.get("grasses", [])
            for seg_grass in sam3_grasses:
                zones.append({
                    "id": seg_grass.get("id", f"grass_sam3_{len(zones)}"),
                    "geometry": seg_grass.get("geometry", {}),
                    "type": "grass",
                    "confidence": seg_grass.get("confidence", 0.5),
                    "source": "sam3_segmentation"
                })
        
        # Create bounds polygon
        bounds_coords = [
            [lon - 0.001, lat - 0.001],
            [lon + 0.001, lat - 0.001],
            [lon + 0.001, lat + 0.001],
            [lon - 0.001, lat + 0.001],
            [lon - 0.001, lat - 0.001]
        ]
        
        tree_zones = TreeZones(
            encoding="simple",
            data=f"trees:{len(trees)},zones:{len(zones)}",
            resolution_m=0.5,
            bounds=GeoJSONPolygon(coordinates=[bounds_coords]),
            coverage_percent=min(50, len(trees) * 0.5 + len(zones) * 5),
            trees=trees
        )
        
        # Elevation map
        corner_samples = [
            ElevationSample(
                point=s.get("point", [0, 0]),
                elevation_m=s.get("elevation_m", 0)
            ) for s in elevation_data.get("corner_samples", [])
        ]
        
        dem_ref_data = elevation_data.get("dem_reference")
        dem_reference = None
        if dem_ref_data:
            dem_reference = DEMReference(
                source=dem_ref_data.get("source", "unknown"),
                resolution_m=dem_ref_data.get("resolution_m", 90.0),
                file_path=dem_ref_data.get("file_path")
            )
        
        elevation_map = ElevationMap(
            corner_samples=corner_samples,
            slope_percent=elevation_data.get("slope_percent", 0),
            slope_direction=elevation_data.get("slope_direction", "flat"),
            average_elevation_m=elevation_data.get("average_elevation_m", 0),
            terrain_classification=elevation_data.get("terrain_classification", "flat"),
            dem_reference=dem_reference
        )
        
        # Water features - include ALL water features that intersect the 60m context box
        context_box_m = 60.0
        water_features_filtered = []
        nearest_water = None
        nearest_distance = float('inf')
        
        for water in water_features_raw:
            water_geom = water.get("geometry", {})
            if not water_geom:
                continue
            
            water_intersects = False
            
            if water_geom.get("type") == "Polygon":
                water_coords = water_geom.get("coordinates", [[]])[0]
                if water_coords and len(water_coords) >= 4:
                    # Convert to local meters to check intersection
                    local_water = GeometryUtils.degrees_to_local(water_coords, ref_lon, ref_lat)
                    
                    # Check if any point is within 60m box OR if any edge crosses the box
                    for pt in local_water:
                        if -context_box_m <= pt[0] <= context_box_m and -context_box_m <= pt[1] <= context_box_m:
                            water_intersects = True
                            break
                    
                    # Also check if edges cross the boundary
                    if not water_intersects:
                        for i in range(len(local_water)):
                            p1 = local_water[i]
                            p2 = local_water[(i + 1) % len(local_water)]
                            min_x, max_x = min(p1[0], p2[0]), max(p1[0], p2[0])
                            min_y, max_y = min(p1[1], p2[1]), max(p1[1], p2[1])
                            if not (max_x < -context_box_m or min_x > context_box_m or 
                                   max_y < -context_box_m or min_y > context_box_m):
                                water_intersects = True
                                break
                    
                    if water_intersects:
                        water_features_filtered.append(water)
                        # Calculate distance for nearest water body
                        water_centroid = GeometryUtils.polygon_centroid(water_coords)
                        p_centroid = GeometryUtils.polygon_centroid(property_coords)
                        local_w = GeometryUtils.degrees_to_local([list(water_centroid)], ref_lon, ref_lat)[0]
                        local_p = GeometryUtils.degrees_to_local([list(p_centroid)], ref_lon, ref_lat)[0]
                        distance = ((local_w[0] - local_p[0])**2 + (local_w[1] - local_p[1])**2)**0.5
                        
                        if distance < nearest_distance:
                            nearest_distance = distance
                            dx = local_w[0] - local_p[0]
                            dy = local_w[1] - local_p[1]
                            if abs(dx) > abs(dy):
                                direction = "east" if dx > 0 else "west"
                            else:
                                direction = "north" if dy > 0 else "south"
                            nearest_water = NearestWaterBody(
                                type=water.get("type", "water"),
                                distance_m=round(distance, 0),
                                direction=direction
                            )
            
            elif water_geom.get("type") == "LineString":
                water_coords = water_geom.get("coordinates", [])
                if water_coords and len(water_coords) >= 2:
                    # Convert to local meters
                    local_water = GeometryUtils.degrees_to_local(water_coords, ref_lon, ref_lat)
                    
                    # Check if any point is within 60m box OR if any segment crosses the box
                    for pt in local_water:
                        if -context_box_m <= pt[0] <= context_box_m and -context_box_m <= pt[1] <= context_box_m:
                            water_intersects = True
                            break
                    
                    if not water_intersects:
                        for i in range(len(local_water) - 1):
                            p1, p2 = local_water[i], local_water[i + 1]
                            min_x, max_x = min(p1[0], p2[0]), max(p1[0], p2[0])
                            min_y, max_y = min(p1[1], p2[1]), max(p1[1], p2[1])
                            if not (max_x < -context_box_m or min_x > context_box_m or 
                                   max_y < -context_box_m or min_y > context_box_m):
                                water_intersects = True
                                break
                    
                    if water_intersects:
                        water_features_filtered.append(water)
                        # Calculate distance to nearest point for LineString
                        p_centroid = GeometryUtils.polygon_centroid(property_coords)
                        min_dist = float('inf')
                        for pt in water_coords:
                            local_pt = GeometryUtils.degrees_to_local([list(pt)], ref_lon, ref_lat)[0]
                            local_p = GeometryUtils.degrees_to_local([list(p_centroid)], ref_lon, ref_lat)[0]
                            dist = ((local_pt[0] - local_p[0])**2 + (local_pt[1] - local_p[1])**2)**0.5
                            min_dist = min(min_dist, dist)
                        
                        if min_dist < nearest_distance:
                            nearest_distance = min_dist
                            # For LineString, use direction to nearest point
                            nearest_pt = min(water_coords, key=lambda p: self._haversine_distance(
                                p_centroid[1], p_centroid[0], p[1], p[0]
                            ))
                            local_nearest = GeometryUtils.degrees_to_local([list(nearest_pt)], ref_lon, ref_lat)[0]
                            local_p = GeometryUtils.degrees_to_local([list(p_centroid)], ref_lon, ref_lat)[0]
                            dx = local_nearest[0] - local_p[0]
                            dy = local_nearest[1] - local_p[1]
                            if abs(dx) > abs(dy):
                                direction = "east" if dx > 0 else "west"
                            else:
                                direction = "north" if dy > 0 else "south"
                            nearest_water = NearestWaterBody(
                                type=water.get("type", "water"),
                                distance_m=round(min_dist, 0),
                                direction=direction
                            )
        
        water_features = WaterFeatures(
            encoding="rle",
            data="0x00000000" if not water_features_filtered else f"0x{len(water_features_filtered):08x}",
            resolution_m=1.0,
            nearest_water_body=nearest_water,
            features=[{
                "id": w.get("id"),
                "type": w.get("type"),
                "geometry": w.get("geometry")  # Include geometry for visualization
            } for w in water_features_filtered]  # Include ALL filtered water features, no limit
        )
        
        return SurroundingContext(
            buildings=buildings,
            tree_zones=tree_zones,
            roads=roads,
            elevation_map=elevation_map,
            water_features=water_features
        )
    
    def _build_analysis(
        self,
        shadow_result: Dict[str, Any],
        adjacency_result: List[Dict[str, Any]]
    ) -> Analysis:
        """Build Analysis model"""
        
        # Shadow analysis
        facade_scores = {}
        for direction, scores in shadow_result.get("facade_scores", {}).items():
            facade_scores[direction] = FacadeScore(
                winter_avg=scores.get("winter_avg", 0.5),
                summer_avg=scores.get("summer_avg", 0.7),
                annual_avg=scores.get("annual_avg", 0.6)
            )
        
        shadow_hours = shadow_result.get("shadow_hours_per_day", {})
        sun_params = shadow_result.get("sun_path_params", {})
        
        shadow_analysis = ShadowAnalysis(
            facade_scores=facade_scores,
            shadow_hours_per_day=ShadowHoursPerDay(
                winter_solstice=shadow_hours.get("winter_solstice", 4.5),
                summer_solstice=shadow_hours.get("summer_solstice", 10.0),
                equinox=shadow_hours.get("equinox", 7.5)
            ),
            neighbor_shadow_angles=shadow_result.get("neighbor_shadow_angles", {}),
            best_solar_facade=shadow_result.get("best_solar_facade", "south"),
            sunlight_score=shadow_result.get("sunlight_score", 0.5),
            computed_at=shadow_result.get("computed_at", datetime.utcnow().isoformat()),
            sun_path_params=SunPathParams(
                latitude=sun_params.get("latitude", 51.5),
                longitude=sun_params.get("longitude", -0.1),
                timezone=sun_params.get("timezone", "Europe/London")
            )
        )
        
        # Adjacency
        adjacency_edges = []
        for edge in adjacency_result:
            geom = edge.get("geometry", {})
            adjacent_to = edge.get("adjacent_to", "unknown")
            
            # Create appropriate edge type based on adjacent_to
            if adjacent_to == "street":
                adjacency_edges.append(StreetAdjacencyEdge(
                    edge_id=edge.get("edge_id", "unknown"),
                    geometry=GeoJSONLineString(coordinates=geom.get("coordinates", [])),
                    length_m=edge.get("length_m", 0),
                    street_name=edge.get("street_name", "Unknown Road"),
                    street_type=edge.get("street_type", "residential"),
                    primary_access=edge.get("primary_access", False),
                    noise_level=edge.get("noise_level", "low"),
                    privacy_exposure=edge.get("privacy_exposure", "low")
                ))
            elif adjacent_to == "water_boundary":
                adjacency_edges.append(WaterAdjacencyEdge(
                    edge_id=edge.get("edge_id", "unknown"),
                    geometry=GeoJSONLineString(coordinates=geom.get("coordinates", [])),
                    length_m=edge.get("length_m", 0),
                    water_type=edge.get("water_type", "unknown"),
                    noise_level=edge.get("noise_level", "low"),
                    privacy_exposure=edge.get("privacy_exposure", "medium")
                ))
            else:  # neighbor_parcel or unknown
                adjacency_edges.append(BuildingAdjacencyEdge(
                    edge_id=edge.get("edge_id", "unknown"),
                    geometry=GeoJSONLineString(coordinates=geom.get("coordinates", [])),
                    length_m=edge.get("length_m", 0),
                    neighbor_id=edge.get("neighbor_id"),
                    shared_wall_potential=edge.get("shared_wall_potential", False),
                    distance_to_neighbor_building_m=edge.get("distance_to_neighbor_building_m"),
                    privacy_exposure=edge.get("privacy_exposure", "low")
                ))
        
        return Analysis(
            shadow_analysis=shadow_analysis,
            adjacency=adjacency_edges
        )
    
    def _build_access(
        self,
        property_coords: List[List[float]],
        roads: List[Dict[str, Any]],
        adjacency_result: List[Dict[str, Any]],
        property_line_obj: Optional[Any] = None
    ) -> Access:
        """Build Access model"""
        
        # Try to get access point from front edges of property line
        access_coords = None
        side = "estimated"
        confidence = "low"
        determined_by = "estimation"
        
        if property_line_obj and hasattr(property_line_obj, 'front') and property_line_obj.front:
            # Calculate center of front edges (midpoint along the line by length)
            front_segment = property_line_obj.front
            front_coords = front_segment.get_coordinates(property_line_obj.coordinates)
            
            if front_coords and len(front_coords) >= 2:
                # Calculate midpoint along the line (not centroid of points)
                access_coords = self._calculate_midpoint_along_line(front_coords)
                side = "front"
                confidence = "high"
                determined_by = "front_edge_center"
        
        # Fallback to adjacency-based method if front edges not available
        if access_coords is None:
            primary_edge = None
            for edge in adjacency_result:
                if edge.get("adjacent_to") == "street":
                    primary_edge = edge
                    break
            
            if primary_edge:
                geom = primary_edge.get("geometry", {})
                coords = geom.get("coordinates", [[0, 0], [0, 0]])
                if len(coords) >= 2:
                    mid_lon = (coords[0][0] + coords[1][0]) / 2
                    mid_lat = (coords[0][1] + coords[1][1]) / 2
                    access_coords = [mid_lon, mid_lat]
                    side = "north" if "north" in primary_edge.get("edge_id", "") else "street"
                    confidence = "high"
                    determined_by = "street_adjacency"
                else:
                    access_coords = [property_coords[0][0], property_coords[0][1]] if property_coords else [0, 0]
            else:
                # Default to first point
                access_coords = [property_coords[0][0], property_coords[0][1]] if property_coords else [0, 0]
        
        primary_access = PrimaryAccessPoint(
            location=GeoJSONPoint(coordinates=access_coords),
            side=side,
            confidence=confidence,
            determined_by=determined_by,
            alternatives=[]
        )
        
        # Vehicle access
        has_road_access = any(e.get("adjacent_to") == "street" for e in adjacency_result)
        vehicle_access = VehicleAccess(
            available=has_road_access,
            type="direct_from_street" if has_road_access else "no_access",
            driveway_width_required_m=3.5,
            turning_radius_adequate=True
        )
        
        # Pedestrian access
        pedestrian_access = PedestrianAccess(
            available=True,
            sidewalk_present=has_road_access,
            grade_accessible=True
        )
        
        return Access(
            primary_access_point=primary_access,
            vehicle_access=vehicle_access,
            pedestrian_access=pedestrian_access
        )
    
    def _calculate_midpoint_along_line(self, coords: List[List[float]]) -> List[float]:
        """
        Calculate the midpoint along a polyline by length (not centroid of points)
        
        Args:
            coords: List of [lon, lat] coordinates forming a polyline
            
        Returns:
            [lon, lat] of the midpoint along the line
        """
        if not coords or len(coords) < 2:
            # Fallback: return first point or average
            if coords:
                return coords[0]
            return [0, 0]
        
        if len(coords) == 2:
            # Simple case: midpoint of single segment
            return [
                (coords[0][0] + coords[1][0]) / 2,
                (coords[0][1] + coords[1][1]) / 2
            ]
        
        # Calculate total length and segment lengths
        import math
        segment_lengths = []
        total_length = 0.0
        
        for i in range(len(coords) - 1):
            p1 = coords[i]
            p2 = coords[i + 1]
            
            # Calculate distance using Haversine (approximate for small distances)
            avg_lat = (p1[1] + p2[1]) / 2
            m_per_deg_lat = 111000
            m_per_deg_lon = 111000 * abs(math.cos(math.radians(avg_lat)))
            m_per_deg = (m_per_deg_lat + m_per_deg_lon) / 2
            
            dx = (p2[0] - p1[0]) * m_per_deg_lon
            dy = (p2[1] - p1[1]) * m_per_deg_lat
            segment_length = math.sqrt(dx**2 + dy**2)
            
            segment_lengths.append(segment_length)
            total_length += segment_length
        
        if total_length == 0:
            # All points are the same, return first point
            return coords[0]
        
        # Find midpoint: point at half the total length
        target_length = total_length / 2
        accumulated_length = 0.0
        
        for i in range(len(segment_lengths)):
            segment_len = segment_lengths[i]
            
            if accumulated_length + segment_len >= target_length:
                # Midpoint is on this segment
                remaining = target_length - accumulated_length
                ratio = remaining / segment_len if segment_len > 0 else 0
                
                p1 = coords[i]
                p2 = coords[i + 1]
                
                mid_lon = p1[0] + (p2[0] - p1[0]) * ratio
                mid_lat = p1[1] + (p2[1] - p1[1]) * ratio
                
                return [mid_lon, mid_lat]
            
            accumulated_length += segment_len
        
        # Fallback: return last point (shouldn't happen)
        return coords[-1]
    
    def _build_existing_structures(
        self,
        existing_buildings: List[Dict[str, Any]]
    ) -> List[ExistingStructure]:
        """Build ExistingStructure models from pre-filtered buildings"""
        existing = []
        
        for building in existing_buildings:
            footprint = building.get("footprint", {})
            b_coords = footprint.get("coordinates", [[]])
            if not b_coords or not b_coords[0]:
                continue
            
            existing.append(ExistingStructure(
                id=building.get("id", f"existing_{len(existing)}"),
                footprint=GeoJSONPolygon(coordinates=b_coords),
                type=building.get("building_type", "unknown"),
                condition="unknown",
                stories=building.get("stories", 2),
                height_m=building.get("height_m", 8.0),
                year_built=None,
                status="unknown",
                affects_buildable_area=True
            ))
        
        return existing
    
    def _point_in_polygon(
        self,
        point: Tuple[float, float],
        polygon: List[Tuple[float, float]]
    ) -> bool:
        """Simple point-in-polygon test"""
        x, y = point
        n = len(polygon)
        inside = False
        
        j = n - 1
        for i in range(n):
            xi, yi = polygon[i]
            xj, yj = polygon[j]
            
            if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
                inside = not inside
            j = i
        
        return inside
    
    def _point_in_polygon_degrees(
        self,
        lon: float,
        lat: float,
        polygon: List[List[float]]
    ) -> bool:
        """Point-in-polygon test for lat/lon coordinates"""
        n = len(polygon)
        inside = False
        
        j = n - 1
        for i in range(n):
            xi, yi = polygon[i]
            xj, yj = polygon[j]
            
            if ((yi > lat) != (yj > lat)) and (lon < (xj - xi) * (lat - yi) / (yj - yi) + xi):
                inside = not inside
            j = i
        
        return inside
    
    def _haversine_distance(
        self,
        lat1: float,
        lon1: float,
        lat2: float,
        lon2: float
    ) -> float:
        """Calculate distance between two points in meters"""
        import math
        R = 6371000  # Earth radius in meters
        
        phi1 = math.radians(lat1)
        phi2 = math.radians(lat2)
        delta_phi = math.radians(lat2 - lat1)
        delta_lambda = math.radians(lon2 - lon1)
        
        a = math.sin(delta_phi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda/2)**2
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        
        return R * c
    
    def _build_regulatory(self, property_area: float) -> Regulatory:
        """Build Regulatory model"""
        
        zoning = self.setback_calculator.get_zoning_info()
        constraints = self.setback_calculator.get_building_constraints(property_area)
        
        return Regulatory(
            zoning=Zoning(
                designation=zoning["designation"],
                description=zoning["description"],
                permitted_uses=zoning["permitted_uses"],
                conditional_uses=zoning["conditional_uses"],
                source="uk_planning_system",
                last_verified=datetime.utcnow().isoformat()
            ),
            building_constraints=BuildingConstraints(
                max_coverage_ratio=constraints["max_coverage_ratio"],
                max_coverage_sqm=constraints["max_coverage_sqm"],
                max_height_m=constraints["max_height_m"],
                max_stories=constraints["max_stories"],
                max_stories_note=constraints.get("max_stories_note"),
                setbacks_m=SetbacksM(
                    front=constraints["setbacks_m"]["front"],
                    rear=constraints["setbacks_m"]["rear"],
                    side_east=constraints["setbacks_m"]["side_east"],
                    side_west=constraints["setbacks_m"]["side_west"],
                    notes=constraints["setbacks_m"].get("notes")
                )
            ),
            additional_restrictions=[]
        )
    
    def _build_data_quality(
        self, 
        boundary_data: Dict[str, Any],
        soil_data: Dict[str, Any],
        elevation_data: Dict[str, Any]
    ) -> DataQuality:
        """Build DataQuality model with actual data sources used"""
        
        # Determine if sources are real data or placeholders
        def is_real_source(source: str, data_type: str) -> tuple[str, bool]:
            """Check if source is real data or placeholder. Returns (source_name, is_real)"""
            placeholder_sources = ["unknown", "default_estimate", "estimated", "placeholder"]
            
            if source in placeholder_sources:
                return f"{source}_placeholder", False
            
            # Real data sources mapping
            real_sources = {
                # Soil sources
                "local_uk_soil": ("BGS SPMM 1km (local)", True),
                "soilgrids_wcs": ("SoilGrids WCS", True),
                "soilgrids_rest": ("SoilGrids REST API", True),
                "bgs_wms": ("BGS WMS", True),
                # Boundary sources
                "inspire_cadastral": ("INSPIRE Cadastral", True),
                "inspire_cadastral_merged": ("INSPIRE Cadastral (merged)", True),
                "openstreetmap": ("OpenStreetMap", True),
                # Elevation sources
                "mapbox_terrain_rgb": ("Mapbox Terrain RGB", True),
                "open-meteo": ("Open-Meteo Elevation", True),
                "open_meteo_elevation": ("Open-Meteo Elevation", True),
            }
            
            if source in real_sources:
                return real_sources[source]
            
            # Default: assume real if not in placeholder list
            return source, True
        
        # Boundary source
        boundary_source = boundary_data.get("source", "openstreetmap")
        boundary_name, boundary_is_real = is_real_source(boundary_source, "boundary")
        
        # Soil source
        soil_source = soil_data.get("source", "unknown")
        soil_name, soil_is_real = is_real_source(soil_source, "soil")
        
        # Elevation source (check both dem_source and dem_reference.source)
        elevation_source = (
            elevation_data.get("dem_source") or 
            elevation_data.get("dem_reference", {}).get("source", "unknown")
        )
        elevation_name, elevation_is_real = is_real_source(elevation_source, "elevation")
        
        # Get elevation resolution
        elevation_resolution = (
            elevation_data.get("resolution_m") or
            elevation_data.get("dem_reference", {}).get("resolution_m")
        )
        
        sources = [
            DataSource(
                type="boundary",
                source=boundary_name if boundary_is_real else "placeholder",
                date=datetime.utcnow().strftime("%Y-%m-%d") if boundary_is_real else None,
                accuracy_m=boundary_data.get("accuracy_m") if boundary_is_real else None
            ),
            DataSource(
                type="buildings",
                source="openstreetmap",  # Always from OSM
                date=datetime.utcnow().strftime("%Y-%m-%d")
            ),
            DataSource(
                type="elevation",
                source=elevation_name if elevation_is_real else "placeholder",
                resolution_m=elevation_resolution if elevation_is_real else None
            ),
            DataSource(
                type="soil",
                source=soil_name if soil_is_real else "placeholder"
            ),
            DataSource(
                type="roads",
                source="openstreetmap",  # Always from OSM
                date=datetime.utcnow().strftime("%Y-%m-%d")
            )
        ]
        
        missing = [
            "exact_cadastral_boundary",
            "utility_locations",
            "tree_preservation_orders"
        ]
        
        return DataQuality(
            data_sources=sources,
            missing_data=missing,
            last_validated=datetime.utcnow().isoformat()
        )

