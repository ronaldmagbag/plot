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
import os
import math
import subprocess
import sys
from datetime import datetime
from typing import Dict, Any, Optional, List, Tuple
import uuid
from pathlib import Path
from loguru import logger

from .config import get_config, PipelineConfig
from .models import (
    PlotAnalysis, GeoJSONPoint, GeoJSONPolygon, GeoJSONLineString,
    Boundaries, PropertyLine, SetbackLine, BuildableEnvelope, SetbacksApplied,
    SurroundingContext, NeighborBuilding, TreeZones, Road, ElevationMap, 
    ElevationSample, WaterFeatures, NearestWaterBody,
    Analysis, ShadowAnalysis, FacadeScore, ShadowHoursPerDay, SunPathParams, AdjacencyEdge,
    Access, PrimaryAccessPoint, VehicleAccess, PedestrianAccess,
    ExistingStructure, Regulatory, Zoning, BuildingConstraints, SetbacksM,
    Soil, DataQuality, DataSource, AccessCorridor, ConstraintApplied
)

from .collectors import (
    OSMCollector, ElevationCollector, SoilCollector, BoundaryCollector,
    TerrainCollector, VegetationCollector, MapboxImageryCollector
)
from .analysis import ShadowAnalyzer, AdjacencyAnalyzer, SetbackCalculator, GeometryUtils


class PlotAnalysisPipeline:
    """
    Main pipeline to generate complete plot analysis from center coordinates
    
    Usage:
        pipeline = PlotAnalysisPipeline()
        result = pipeline.run(lat=51.5074, lon=-0.1276)
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
        self.vegetation_collector = VegetationCollector()
        self.soil_collector = SoilCollector()
        self.boundary_collector = BoundaryCollector()
        self.mapbox_imagery_collector = MapboxImageryCollector(cache_dir=cache_dir)
        
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
            radius_m: Search radius (default from config)
            plot_id: Optional custom plot ID
        
        Returns:
            Complete PlotAnalysis object
        """
        radius = radius_m or self.config.search_radius_m
        context_radius = self.config.context_radius_m
        
        logger.info(f"Starting plot analysis pipeline for ({lat}, {lon})")
        logger.info(f"Search radius: {radius}m, Context radius: {context_radius}m")
        
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
        osm_landuse_preview = all_osm_features_preview.get("landuse", [])
        
        # Try to get boundary from OSM data first
        boundary_data = None
        if osm_landuse_preview:
            # Try landuse polygons
            for landuse in osm_landuse_preview:
                geom = landuse.get("geometry", {})
                if geom.get("type") == "Polygon":
                    coords = geom.get("coordinates", [[]])[0]
                    if coords and len(coords) >= 4:
                        # Check if point is in polygon
                        if self._point_in_polygon_degrees(lon, lat, coords):
                            boundary_data = self.boundary_collector._create_boundary_from_coords(
                                coords, "openstreetmap_landuse", 2.0
                            )
                            break
        
        # Get roads for boundary derivation
        osm_roads_preview = all_osm_features_preview.get("roads", [])
        
        # Try boundary collector with OSM data (roads, buildings, and landuse) for better boundary detection
        # This allows deriving boundaries from roads and neighbor buildings, not just rectangles
        # Pass landuse data to avoid separate API call
        if not boundary_data:
            boundary_data = self.boundary_collector.get_plot_boundary(
                lat, lon, radius,
                osm_buildings=osm_buildings_preview,
                osm_roads=osm_roads_preview,
                osm_landuse=osm_landuse_preview
            )
        
        property_coords = boundary_data.get("coordinates", [[]])[0]
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
            all_osm_features = self.osm_collector.fetch_all_features(lat, lon, context_radius)
        
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
        
        # Mapbox imagery and SAM3 segmentation
        sam3_results = None
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
            
            # If SAM3 results don't exist, check for cached imagery and run segmentation
            if sam3_results is None:
                # Check if merged image already exists in cache
                zoom = self.config.api.mapbox_max_zoom
                cached_image_pattern = f"mapbox_merged_{lat:.6f}_{lon:.6f}_{context_radius}m_z{zoom}.jpg"
                cached_images = list(Path(self.cache_dir).glob(cached_image_pattern))
                
                merged_image = None
                imagery_metadata = None
                
                if cached_images:
                    logger.info(f"Found cached merged image: {cached_images[0].name}, skipping download")
                    try:
                        from PIL import Image
                        merged_image = Image.open(cached_images[0])
                        # Reconstruct metadata from filename and config
                        imagery_metadata = {
                            "bounds": {
                                "west": lon - context_radius / 111000 / abs(math.cos(math.radians(lat))),
                                "east": lon + context_radius / 111000 / abs(math.cos(math.radians(lat))),
                                "south": lat - context_radius / 111000,
                                "north": lat + context_radius / 111000
                            },
                            "zoom": zoom,
                            "center": [lon, lat],
                            "radius_m": context_radius
                        }
                    except Exception as e:
                        logger.warning(f"Failed to load cached image: {e}, will re-download")
                        cached_images = []
                
                if not cached_images:
                    logger.info("Downloading Mapbox satellite imagery for SAM3 segmentation...")
                    try:
                        merged_image, imagery_metadata = self.mapbox_imagery_collector.download_imagery(
                            lat, lon, radius_m=context_radius
                        )
                    except Exception as e:
                        logger.warning(f"Mapbox imagery download failed: {e}")
                        merged_image = None
                        imagery_metadata = None
                
                if merged_image and imagery_metadata:
                    logger.info("Running SAM3 segmentation...")
                    sam3_results = self._run_sam3_segmentation(
                        merged_image, imagery_metadata, lat, lon
                    )
                    if sam3_results:
                        logger.info(f"SAM3 found: {len(sam3_results.get('roads', []))} roads, "
                                  f"{len(sam3_results.get('trees', []))} trees, "
                                  f"{len(sam3_results.get('grasses', []))} grasses")
                elif not merged_image:
                    logger.warning("No imagery available for SAM3 segmentation, skipping")
        
        # ============================================================
        # STAGE 3: Analysis
        # ============================================================
        logger.info("Stage 3: Running analysis...")
        
        # Calculate TRUE centroid of property polygon
        property_centroid = GeometryUtils.polygon_centroid(property_coords)
        centroid_lon, centroid_lat = property_centroid
        
        # Filter buildings: only keep those OUTSIDE the property boundary
        ref_lon, ref_lat = property_coords[0] if property_coords else (lon, lat)
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
        
        # Preliminary adjacency for setback calculation
        preliminary_adjacency = self._quick_adjacency(property_coords, osm_roads)
        
        # Get property line with classification using new boundary collector
        property_line_obj = self.boundary_collector.get_property_line(
            lat, lon, radius,
            osm_buildings=osm_buildings,
            osm_roads=osm_roads,
            osm_landuse=osm_landuse_preview,
            classify=True
        )
        
        # If we got a property line object, use it; otherwise fall back to old method
        if property_line_obj:
            # Use new boundary collector methods
            setback_line_obj = self.boundary_collector.get_setback_line(property_line_obj)
            buildable_envelope_obj = self.boundary_collector.get_buildable_envelope(
                property_line_obj, setback_line_obj
            )
            
            # Convert to old format for compatibility
            if setback_line_obj:
                setback_result = {
                    "coordinates": [setback_line_obj.coordinates],
                    "area_sqm": setback_line_obj.area_sqm,
                    "setbacks_applied": {
                        "front_m": setback_line_obj.metadata.get("front_rear_setback_m", 4.0),
                        "rear_m": setback_line_obj.metadata.get("front_rear_setback_m", 4.0),
                        "side_east_m": setback_line_obj.metadata.get("side_setback_m", 1.0),
                        "side_west_m": setback_line_obj.metadata.get("side_setback_m", 1.0)
                    },
                    "regulation_source": "uk_planning_guidance"
                }
                setback_coords = setback_line_obj.coordinates
            else:
                setback_result = None
                setback_coords = property_coords
            
            if buildable_envelope_obj:
                buildable_result = {
                    "coordinates": [buildable_envelope_obj.coordinates],
                    "area_sqm": buildable_envelope_obj.area_sqm,
                    "constraints_applied": []
                }
            else:
                buildable_result = None
        else:
            # Fallback to old method
            setback_result = self.setback_calculator.calculate_setbacks(
                property_coords, preliminary_adjacency
            )
            setback_coords = setback_result["coordinates"][0] if setback_result else property_coords
            
            # Calculate buildable envelope (pass neighbor buildings for constraint checking)
            buildable_result = self.setback_calculator.calculate_buildable_envelope(
                setback_coords, 
                constraints=[],
                property_area_sqm=property_area,
                neighbor_buildings=neighbor_buildings
            )
        
        # Full adjacency analysis
        adjacency_result = self.adjacency_analyzer.analyze(
            property_coords, osm_roads, osm_buildings, water_features_raw
        )
        
        # Shadow analysis
        shadow_result = self.shadow_analyzer.analyze(
            (lon, lat), property_coords, osm_buildings
        )
        
        # ============================================================
        # STAGE 4: Build Output Model
        # ============================================================
        logger.info("Stage 4: Building output model...")
        
        # Use calculated centroid from property polygon (not input lat/lon)
        centroid = GeoJSONPoint(coordinates=[centroid_lon, centroid_lat])
        
        # Boundaries
        boundaries = self._build_boundaries(
            property_coords, property_area, property_perimeter,
            setback_result, buildable_result, boundary_data,
            property_line_obj  # Pass property line object for segments
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
        access = self._build_access(property_coords, osm_roads, adjacency_result)
        
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
            source=soil_data.get("source", "soilgrids")
        )
        
        # Data quality
        data_quality = self._build_data_quality(boundary_data)
        
        # Build complete plot analysis
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
            json.dump(plot_analysis.model_dump(), f, indent=2, ensure_ascii=False)
        
        logger.info(f"Saved plot analysis to {output_path}")
        return output_path
    
    # ============================================================
    # Helper Methods
    # ============================================================
    
    def _run_sam3_segmentation(
        self,
        merged_image,
        imagery_metadata: Dict[str, Any],
        lat: float,
        lon: float
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
            context_radius = self.config.context_radius_m
            
            # Find the cached merged image
            cached_image_pattern = f"mapbox_merged_{lat:.6f}_{lon:.6f}_{context_radius}m_z{zoom}.jpg"
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
            # __file__ is at src/pipeline.py, so go up to plot/ and then to scripts/
            script_path = Path(__file__).parent.parent / "scripts" / "sam3_segment.py"
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
            
            # Load results
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
        
        # Simple check - find the edge closest to a road
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
        
        # Check for road adjacency
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
        property_line_obj: Optional[Any] = None
    ) -> Boundaries:
        """Build Boundaries model from collected data"""
        
        # Property line - include segments if available
        property_line = PropertyLine(
            coordinates=[property_coords],
            area_sqm=round(property_area, 1),
            perimeter_m=round(property_perimeter, 1),
            source=boundary_data.get("source", "estimated"),
            source_date=datetime.utcnow().isoformat(),
            accuracy_m=boundary_data.get("accuracy_m", 5.0)
        )
        
        # Add segments if property_line_obj has them
        segments_data = None
        if property_line_obj and hasattr(property_line_obj, 'front'):
            # Store segment information for visualization
            # Store edge indices instead of coordinates to avoid duplication
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
        
        # Update property_line with segments (using model_dump and update)
        if segments_data:
            property_line_dict = property_line.model_dump()
            property_line_dict["segments"] = segments_data
            property_line = PropertyLine(**property_line_dict)
        
        # Setback line
        setback_coords = setback_result["coordinates"] if setback_result else [property_coords]
        setbacks_applied = setback_result.get("setbacks_applied", {}) if setback_result else {}
        
        setback_line = SetbackLine(
            coordinates=setback_coords,
            area_sqm=setback_result.get("area_sqm", property_area * 0.7) if setback_result else property_area * 0.7,
            setbacks_applied=SetbacksApplied(
                front_m=setbacks_applied.get("front_m", 5.0),
                rear_m=setbacks_applied.get("rear_m", 5.0),
                side_east_m=setbacks_applied.get("side_east_m", 1.0),
                side_west_m=setbacks_applied.get("side_west_m", 1.0)
            ),
            regulation_source=setback_result.get("regulation_source", "uk_planning_guidance") if setback_result else "uk_planning_guidance"
        )
        
        # Buildable envelope
        buildable_coords = buildable_result["coordinates"] if buildable_result else setback_coords
        constraints = buildable_result.get("constraints_applied", []) if buildable_result else []
        access_corridor_data = buildable_result.get("access_corridor", {}) if buildable_result else {}
        
        access_corridor = None
        if access_corridor_data.get("required"):
            access_corridor = AccessCorridor(
                required=True,
                width_m=access_corridor_data.get("width_m", 3.5),
                description=access_corridor_data.get("description", "Vehicle access corridor"),
                affects_buildable_area=access_corridor_data.get("affects_buildable_area", False)
            )
        
        buildable_envelope = BuildableEnvelope(
            coordinates=buildable_coords,
            area_sqm=buildable_result.get("area_sqm", property_area * 0.6) if buildable_result else property_area * 0.6,
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
            property_line=property_line,
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
            
            # Get building centroid for wall facing calculation
            b_centroid = GeometryUtils.polygon_centroid(b_coords)
            p_centroid = GeometryUtils.polygon_centroid(property_coords)
            local_b_centroid = GeometryUtils.degrees_to_local([list(b_centroid)], ref_lon, ref_lat)[0]
            local_p_centroid = GeometryUtils.degrees_to_local([list(p_centroid)], ref_lon, ref_lat)[0]
            
            # Determine wall facing plot
            dx = local_b_centroid[0] - local_p_centroid[0]
            dy = local_b_centroid[1] - local_p_centroid[1]
            if abs(dx) > abs(dy):
                wall_facing = "west" if dx > 0 else "east"
            else:
                wall_facing = "south" if dy > 0 else "north"
            
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
        
        elevation_map = ElevationMap(
            corner_samples=corner_samples,
            slope_percent=elevation_data.get("slope_percent", 0),
            slope_direction=elevation_data.get("slope_direction", "flat"),
            average_elevation_m=elevation_data.get("average_elevation_m", 0),
            terrain_classification=elevation_data.get("terrain_classification", "flat")
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
                "type": w.get("type")
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
            
            adjacency_edges.append(AdjacencyEdge(
                edge_id=edge.get("edge_id", "unknown"),
                geometry=GeoJSONLineString(coordinates=geom.get("coordinates", [])),
                length_m=edge.get("length_m", 0),
                adjacent_to=edge.get("adjacent_to", "unknown"),
                street_name=edge.get("street_name"),
                street_type=edge.get("street_type"),
                neighbor_id=edge.get("neighbor_id"),
                primary_access=edge.get("primary_access", False),
                shared_wall_potential=edge.get("shared_wall_potential", False),
                distance_to_neighbor_building_m=edge.get("distance_to_neighbor_building_m"),
                noise_level=edge.get("noise_level", "low"),
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
        adjacency_result: List[Dict[str, Any]]
    ) -> Access:
        """Build Access model"""
        
        # Find primary access from adjacency
        primary_edge = None
        for edge in adjacency_result:
            if edge.get("adjacent_to") == "street":
                primary_edge = edge
                break
        
        # Determine access point location
        if primary_edge:
            geom = primary_edge.get("geometry", {})
            coords = geom.get("coordinates", [[0, 0], [0, 0]])
            if len(coords) >= 2:
                mid_lon = (coords[0][0] + coords[1][0]) / 2
                mid_lat = (coords[0][1] + coords[1][1]) / 2
                access_coords = [mid_lon, mid_lat]
                side = "north" if "north" in primary_edge.get("edge_id", "") else "street"
            else:
                access_coords = [property_coords[0][0], property_coords[0][1]]
                side = "estimated"
        else:
            # Default to first point
            access_coords = [property_coords[0][0], property_coords[0][1]] if property_coords else [0, 0]
            side = "estimated"
        
        primary_access = PrimaryAccessPoint(
            location=GeoJSONPoint(coordinates=access_coords),
            side=side,
            confidence="high" if primary_edge else "low",
            determined_by="street_adjacency" if primary_edge else "estimation",
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
    
    def _build_data_quality(self, boundary_data: Dict[str, Any]) -> DataQuality:
        """Build DataQuality model"""
        
        sources = [
            DataSource(
                type="boundary",
                source=boundary_data.get("source", "openstreetmap"),
                date=datetime.utcnow().strftime("%Y-%m-%d"),
                accuracy_m=boundary_data.get("accuracy_m", 5.0)
            ),
            DataSource(
                type="buildings",
                source="openstreetmap",
                date=datetime.utcnow().strftime("%Y-%m-%d")
            ),
            DataSource(
                type="elevation",
                source="open-meteo",
                resolution_m=90.0
            ),
            DataSource(
                type="soil",
                source="soilgrids_isric"
            ),
            DataSource(
                type="roads",
                source="openstreetmap"
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

