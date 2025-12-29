"""
Soil data collector using SoilGrids WCS/API and BGS WMS
Collects soil properties for foundation recommendations
"""

import time
import json
from pathlib import Path
from typing import Dict, Any, Optional
import requests
from io import BytesIO
from loguru import logger
from rasterio.io import MemoryFile
import rasterio
from shapely.geometry import Point
import geopandas as gpd

try:
    from owslib.wcs import WebCoverageService
    from owslib.wms import WebMapService
    OWSLIB_AVAILABLE = True
except ImportError:
    OWSLIB_AVAILABLE = False
    logger.warning("OWSLib not available. Install with: pip install owslib")

from ..config import get_config


class SoilCollector:
    """Collect soil data from SoilGrids WCS/API (ISRIC) and BGS WMS (UK)"""
    
    def __init__(self):
        self.config = get_config()
        self.rest_url = self.config.api.soilgrids_url
        # SoilGrids WCS endpoint (more stable than REST API)
        self.wcs_url = self.config.api.soilgrids_wcs_url
        # BGS WMS for UK locations
        self.bgs_wms_url = self.config.api.bgs_wms_url
        self._last_request_time = 0
        self._min_request_interval = 0.5
        
        # Local UK soil data paths (if available)
        self._init_uk_soil_data()
    
    def _rate_limit(self):
        """Ensure we don't exceed rate limits"""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()
    
    def get_soil_data(self, lat: float, lon: float) -> Dict[str, Any]:
        """
        Get soil properties for a location
        
        Returns soil type, drainage, bearing capacity estimates
        
        Tries multiple data sources in order:
        1. Local UK soil data files (for UK locations, if available)
        2. SoilGrids WCS (global)
        3. SoilGrids REST API (if WCS fails)
        4. BGS WMS (for UK locations only)
        5. Unknown soil type (if all sources fail)
        """
        logger.info(f"Fetching soil data for ({lat}, {lon})")
        
        # For UK locations, try local data first (most reliable if available)
        if self._is_uk_location(lat, lon):
            properties = self._query_local_uk_soil(lat, lon)
            if properties:
                source = "local_uk_soil"
                soil_type = self._classify_soil_type(properties)
                drainage = self._estimate_drainage(properties)
                bearing_capacity = self._estimate_bearing_capacity(properties)
                foundation_rec = self._recommend_foundation(soil_type, drainage, bearing_capacity)
                
                return {
                    "type_id": soil_type,
                    "bearing_capacity_kpa": bearing_capacity,
                    "drainage": drainage,
                    "foundation_recommendation": foundation_rec,
                    "source": source,
                    "properties": properties
                }
        
        # Try SoilGrids WCS (global coverage)
        properties = self._query_soilgrids_wcs(lat, lon)
        source = "soilgrids_wcs"
        
        # Fallback to REST API if WCS fails
        if not properties:
            logger.debug("WCS failed, trying REST API")
            properties = self._query_soilgrids_rest(lat, lon)
            source = "soilgrids_rest"
        
        # For UK locations, try BGS WMS as additional fallback
        if not properties and self._is_uk_location(lat, lon):
            logger.debug("Trying BGS WMS for UK location")
            properties = self._query_bgs_wms(lat, lon)
            source = "bgs_wms"
        
        if properties:
            soil_type = self._classify_soil_type(properties)
            drainage = self._estimate_drainage(properties)
            bearing_capacity = self._estimate_bearing_capacity(properties)
            foundation_rec = self._recommend_foundation(soil_type, drainage, bearing_capacity)
            
            return {
                "type_id": soil_type,
                "bearing_capacity_kpa": bearing_capacity,
                "drainage": drainage,
                "foundation_recommendation": foundation_rec,
                "source": source,
                "properties": properties
            }
        
        # Return unknown if all data sources fail
        logger.warning("All soil data sources failed, returning unknown soil type")
        return self._get_unknown_soil()
    
    def _init_uk_soil_data(self):
        """Initialize local UK soil data paths"""
        # Look for UK soil data in multiple locations
        project_root = Path(__file__).parent.parent.parent
        soil_dir = project_root / "soil" / "data"
        data_dir = project_root / "data"
        
        self.uk_soil_data = {
            "clay": soil_dir / "ukts_clay.tif",
            "sand": soil_dir / "ukts_sand.tif",
            "silt": soil_dir / "ukts_silt.tif",
            "parent_material": soil_dir / "soil_parent_material.gpkg",
            # BGS SPMM 1km shapefile (downloaded)
            "spmm_1km": data_dir / "SPMM_1km" / "SoilParentMateriall_V1_portal1km.shp",
            # Geology GPKG (downloaded)
            "geology": data_dir / "625k_V5_Geology_UK_EPSG27700.gpkg"
        }
        
        self._uk_rasters = {}
        self._uk_soil_gdf = None
        self._uk_spmm_gdf = None
        self._uk_geology_gdf = None
        
        # Load available rasters (UKTS)
        for key, path in self.uk_soil_data.items():
            if key in ["clay", "sand", "silt"] and path.exists():
                try:
                    self._uk_rasters[key] = rasterio.open(path)
                    logger.info(f"Loaded UK soil raster: {key}")
                except Exception as e:
                    logger.debug(f"Failed to load UK soil raster {key}: {e}")
        
        # Load SPMM 1km shapefile (BGS Soil Parent Material Model)
        if self.uk_soil_data["spmm_1km"].exists():
            try:
                # Load in native CRS (EPSG:27700) for better performance
                self._uk_spmm_gdf = gpd.read_file(self.uk_soil_data["spmm_1km"])
                logger.info(f"Loaded BGS SPMM 1km shapefile: {len(self._uk_spmm_gdf)} polygons")
                logger.debug(f"SPMM CRS: {self._uk_spmm_gdf.crs}")
                logger.debug(f"SPMM columns: {list(self._uk_spmm_gdf.columns)}")
            except Exception as e:
                logger.warning(f"Failed to load BGS SPMM shapefile: {e}")
        
        # Load parent material polygons if available (alternative source)
        if self.uk_soil_data["parent_material"].exists():
            try:
                self._uk_soil_gdf = gpd.read_file(self.uk_soil_data["parent_material"]).to_crs("EPSG:4326")
                logger.info("Loaded UK soil parent material polygons")
            except Exception as e:
                logger.debug(f"Failed to load UK soil parent material: {e}")
        
        # Load geology GPKG (optional, for additional context)
        if self.uk_soil_data["geology"].exists():
            try:
                # GPKG might have multiple layers, try to load the first one
                import fiona
                layers = fiona.listlayers(str(self.uk_soil_data["geology"]))
                if layers:
                    self._uk_geology_gdf = gpd.read_file(
                        self.uk_soil_data["geology"], 
                        layer=layers[0]
                    )
                    logger.info(f"Loaded geology GPKG: {len(self._uk_geology_gdf)} features from layer '{layers[0]}'")
            except Exception as e:
                logger.debug(f"Failed to load geology GPKG: {e}")
    
    def _is_uk_location(self, lat: float, lon: float) -> bool:
        """Check if location is within UK bounds"""
        # Rough UK bounds: 49.5°N to 61°N, 8°W to 2°E
        return 49.5 <= lat <= 61.0 and -8.0 <= lon <= 2.0
    
    def _parse_soil_texture(self, texture_desc: str) -> Dict[str, float]:
        """
        Parse soil texture description and estimate clay/sand/silt percentages
        
        This is NOT a mathematical formula, but a LOOKUP TABLE based on:
        - USDA Soil Texture Triangle classifications
        - Standard soil science definitions
        
        Method:
        1. Parse keywords from texture description
        2. Match to USDA texture class
        3. Assign representative values (typically midpoints of USDA ranges)
        
        Returns values in g/kg (0-1000) to match SoilGrids format
        - 100 g/kg = 10% (percentage)
        - Example: 300 g/kg = 30% clay
        
        See docs/soil_texture_parsing.md for detailed explanation and USDA ranges.
        """
        texture_lower = texture_desc.lower()
        estimates = {}
        
        # Common UK soil texture classifications
        # Based on USDA texture triangle and UK soil classification
        
        if 'clay' in texture_lower:
            if 'sandy' in texture_lower:
                # Sandy clay
                estimates['clay'] = 400  # 40%
                estimates['sand'] = 450  # 45%
                estimates['silt'] = 150  # 15%
            elif 'silty' in texture_lower:
                # Silty clay
                estimates['clay'] = 400  # 40%
                estimates['sand'] = 100  # 10%
                estimates['silt'] = 500  # 50%
            elif 'heavy' in texture_lower or 'very' in texture_lower:
                # Heavy clay
                estimates['clay'] = 550  # 55%
                estimates['sand'] = 200  # 20%
                estimates['silt'] = 250  # 25%
            elif 'loam' in texture_lower:
                # Clay loam
                estimates['clay'] = 300  # 30%
                estimates['sand'] = 350  # 35%
                estimates['silt'] = 350  # 35%
            else:
                # Plain clay
                estimates['clay'] = 450  # 45%
                estimates['sand'] = 250  # 25%
                estimates['silt'] = 300  # 30%
        
        elif 'sand' in texture_lower:
            if 'loam' in texture_lower:
                # Sandy loam
                estimates['clay'] = 150  # 15%
                estimates['sand'] = 600  # 60%
                estimates['silt'] = 250  # 25%
            elif 'coarse' in texture_lower:
                # Coarse sand
                estimates['clay'] = 50   # 5%
                estimates['sand'] = 850  # 85%
                estimates['silt'] = 100  # 10%
            else:
                # Sand
                estimates['clay'] = 100  # 10%
                estimates['sand'] = 750  # 75%
                estimates['silt'] = 150  # 15%
        
        elif 'silt' in texture_lower:
            if 'loam' in texture_lower:
                # Silt loam
                estimates['clay'] = 150  # 15%
                estimates['sand'] = 200  # 20%
                estimates['silt'] = 650  # 65%
            else:
                # Silt
                estimates['clay'] = 100  # 10%
                estimates['sand'] = 100  # 10%
                estimates['silt'] = 800  # 80%
        
        elif 'loam' in texture_lower:
            # Loam (balanced)
            estimates['clay'] = 200  # 20%
            estimates['sand'] = 400  # 40%
            estimates['silt'] = 400  # 40%
        
        elif 'peat' in texture_lower or 'organic' in texture_lower:
            # Organic/peat soils
            estimates['clay'] = 100  # 10%
            estimates['sand'] = 200  # 20%
            estimates['silt'] = 200  # 20%
            estimates['organic'] = 500  # 50% organic matter
        
        # If we got some estimates, return them
        if estimates:
            return estimates
        
        # Default fallback - try to infer from keywords
        if 'fine' in texture_lower:
            estimates['clay'] = 300
            estimates['sand'] = 300
            estimates['silt'] = 400
        elif 'coarse' in texture_lower or 'gravel' in texture_lower:
            estimates['clay'] = 100
            estimates['sand'] = 700
            estimates['silt'] = 200
        else:
            # Unknown texture - use moderate values
            estimates['clay'] = 250
            estimates['sand'] = 400
            estimates['silt'] = 350
        
        return estimates
    
    def _query_local_uk_soil(self, lat: float, lon: float) -> Optional[Dict[str, Any]]:
        """
        Query local UK soil data files if available
        
        Uses:
        1. BGS SPMM 1km shapefile (primary - most comprehensive)
        2. UKTS rasters (if available)
        3. Other parent material datasets
        """
        results = {}
        
        # Priority 1: Query BGS SPMM 1km shapefile (most comprehensive)
        if self._uk_spmm_gdf is not None:
            try:
                from pyproj import Transformer
                # Transform from WGS84 to British National Grid (EPSG:27700)
                transformer = Transformer.from_crs("EPSG:4326", "EPSG:27700", always_xy=True)
                x, y = transformer.transform(lon, lat)
                point = Point(x, y)
                
                # Find intersecting polygon
                hits = self._uk_spmm_gdf[
                    self._uk_spmm_gdf.geometry.contains(point) | 
                    self._uk_spmm_gdf.geometry.intersects(point)
                ]
                
                if not hits.empty:
                    hit = hits.iloc[0]
                    
                    # Extract soil texture information
                    soil_tex = hit.get('SOIL_TEX', '')
                    if soil_tex:
                        results['soil_texture'] = str(soil_tex)
                        # Parse texture to estimate clay/sand/silt percentages
                        texture_estimates = self._parse_soil_texture(str(soil_tex))
                        results.update(texture_estimates)
                    
                    # Extract other useful properties
                    soil_group = hit.get('SOIL_GROUP', '')
                    if soil_group:
                        results['soil_group'] = str(soil_group)
                    
                    soil_depth = hit.get('SOIL_DEPTH', '')
                    if soil_depth:
                        results['soil_depth'] = str(soil_depth)
                    
                    esb_desc = hit.get('ESB_DESC', '')
                    if esb_desc:
                        results['esb_description'] = str(esb_desc)
                    
                    pmm_grain = hit.get('PMM_GRAIN', '')
                    if pmm_grain:
                        results['parent_material_grain'] = str(pmm_grain)
                    
                    logger.debug(f"BGS SPMM query successful: {results}")
                    
            except Exception as e:
                logger.debug(f"Failed to query BGS SPMM: {e}")
        
        # Priority 2: Sample UKTS rasters for precise texture values (if available)
        if self._uk_rasters:
            for prop in ["clay", "sand", "silt"]:
                if prop in self._uk_rasters and prop not in results:
                    try:
                        raster = self._uk_rasters[prop]
                        for value in raster.sample([(lon, lat)]):
                            val = float(value[0])
                            # UKTS values are typically in percentages
                            # But check if they need conversion (some datasets use 0-1 scale)
                            if val > 0 and val <= 1:
                                val = val * 100  # Convert to percentage
                            # Convert percentage to g/kg (multiply by 10)
                            if val <= 100:
                                results[prop] = val * 10
                            else:
                                results[prop] = val
                            break
                    except Exception as e:
                        logger.debug(f"Failed to sample UK {prop} raster: {e}")
        
        # Priority 3: Try alternative parent material source
        if self._uk_soil_gdf is not None and 'parent_material' not in results:
            try:
                pt = Point(lon, lat)
                hit = self._uk_soil_gdf[self._uk_soil_gdf.contains(pt)]
                if not hit.empty:
                    parent_data = hit.iloc[0].to_dict()
                    if "material" in parent_data or "type" in parent_data:
                        results["parent_material"] = parent_data.get("material") or parent_data.get("type")
            except Exception as e:
                logger.debug(f"Failed to query alternative parent material: {e}")
        
        if results:
            logger.info(f"Local UK soil data query successful: {list(results.keys())}")
            return results
        
        return None
    
    def _query_soilgrids_wcs(self, lat: float, lon: float) -> Optional[Dict[str, Any]]:
        """
        Query SoilGrids using WCS (Web Coverage Service) - more stable than REST API
        
        Uses OWSLib for proper WCS access, falls back to direct requests if OWSLib unavailable
        """
        self._rate_limit()
        
        # Properties and their WCS layer names
        # Using 0-5cm depth as primary (most relevant for foundations)
        wcs_layers = {
            "clay": "clay_0-5cm_mean",
            "sand": "sand_0-5cm_mean",
            "silt": "silt_0-5cm_mean",
            "bdod": "bdod_0-5cm_mean",
            "cec": "cec_0-5cm_mean",
            "phh2o": "phh2o_0-5cm_mean",
            "soc": "soc_0-5cm_mean"
        }
        
        results = {}
        
        # Try using OWSLib first (more reliable)
        if OWSLIB_AVAILABLE:
            try:
                # Try different WCS base URLs - the map parameter might need different format
                wcs_base_urls = [
                    "https://maps.isric.org/mapserv?map=/map/soilgrids.map",
                    "https://maps.isric.org/mapserv?map=soilgrids.map",
                    "https://maps.isric.org/mapserv?map=/map/soilgrids",
                ]
                
                wcs = None
                for base_url in wcs_base_urls:
                    try:
                        wcs = WebCoverageService(
                            f"{base_url}&SERVICE=WCS&VERSION=2.0.1",
                            version='2.0.1',
                            timeout=self.config.api.request_timeout
                        )
                        # If we get here, the connection worked
                        break
                    except Exception:
                        continue
                
                if wcs is None:
                    raise Exception("Could not connect to any WCS endpoint")
                
                # Small bbox around point (0.001 degrees ≈ 100m)
                bbox = (lon - 0.0005, lat - 0.0005, lon + 0.0005, lat + 0.0005)
                crs = 'EPSG:4326'
                
                # Query each property layer
                for prop_name, layer_name in wcs_layers.items():
                    try:
                        # Check if coverage exists
                        if layer_name not in wcs.contents:
                            logger.debug(f"Coverage {layer_name} not found in WCS")
                            continue
                        
                        # Get coverage data
                        response = wcs.getCoverage(
                            identifier=layer_name,
                            bbox=bbox,
                            crs=crs,
                            format='image/tiff',
                            width=1,
                            height=1
                        )
                        
                        # Read TIFF and extract value
                        with MemoryFile(BytesIO(response.read())) as memfile:
                            with memfile.open() as dataset:
                                value = dataset.read(1)[0, 0]
                                # Check for no-data values (typically -32768 or similar)
                                if value is not None and value > -1000:
                                    results[prop_name] = float(value)
                                    
                    except Exception as e:
                        logger.debug(f"WCS query failed for {prop_name}: {e}")
                        continue
                
                if results:
                    logger.debug(f"SoilGrids WCS results (OWSLib): {results}")
                    return results
                    
            except Exception as e:
                logger.debug(f"OWSLib WCS connection failed: {e}, trying direct request")
        
        # Fallback to direct HTTP request
        for prop_name, layer_name in wcs_layers.items():
            try:
                # Small bbox around point
                bbox = f"{lon-0.0005},{lat-0.0005},{lon+0.0005},{lat+0.0005}"
                
                params = {
                    "SERVICE": "WCS",
                    "VERSION": "2.0.1",
                    "REQUEST": "GetCoverage",
                    "COVERAGEID": layer_name,
                    "BBOX": bbox,
                    "FORMAT": "image/tiff",
                    "CRS": "EPSG:4326",
                    "WIDTH": "1",
                    "HEIGHT": "1"
                }
                
                response = requests.get(
                    self.wcs_url,
                    params=params,
                    headers={"User-Agent": self.config.api.user_agent},
                    timeout=self.config.api.request_timeout
                )
                response.raise_for_status()
                
                # Read TIFF and extract value
                with MemoryFile(BytesIO(response.content)) as memfile:
                    with memfile.open() as dataset:
                        value = dataset.read(1)[0, 0]
                        # Check for no-data values
                        if value is not None and value > -1000:
                            results[prop_name] = float(value)
                            
            except Exception as e:
                logger.debug(f"Direct WCS query failed for {prop_name}: {e}")
                continue
        
        if results:
            logger.debug(f"SoilGrids WCS results: {results}")
            return results
        return None
    
    def _query_soilgrids_rest(self, lat: float, lon: float) -> Optional[Dict[str, Any]]:
        """Query SoilGrids REST API for soil properties (fallback)"""
        self._rate_limit()
        
        # Properties we're interested in
        properties = ["clay", "sand", "silt", "bdod", "cec", "phh2o", "soc"]
        depths = ["0-5cm", "5-15cm", "15-30cm"]
        
        params = {
            "lon": lon,
            "lat": lat,
            "property": properties,
            "depth": depths,
            "value": "mean"
        }
        
        try:
            response = requests.get(
                self.rest_url,
                params=params,
                headers={"User-Agent": self.config.api.user_agent},
                timeout=self.config.api.request_timeout
            )
            response.raise_for_status()
            data = response.json()
            
            # Extract mean values for 0-30cm depth
            results = {}
            properties_data = data.get("properties", {})
            layers = properties_data.get("layers", [])
            
            for layer in layers:
                prop_name = layer.get("name")
                depths_data = layer.get("depths", [])
                
                # Average across depths
                values = []
                for depth in depths_data:
                    mean_val = depth.get("values", {}).get("mean")
                    if mean_val is not None:
                        values.append(mean_val)
                
                if values:
                    # Convert units where needed
                    avg = sum(values) / len(values)
                    results[prop_name] = avg
            
            logger.debug(f"SoilGrids REST results: {results}")
            return results if results else None
            
        except Exception as e:
            logger.debug(f"SoilGrids REST API request failed: {e}")
            return None
    
    def _query_bgs_wms(self, lat: float, lon: float) -> Optional[Dict[str, Any]]:
        """
        Query BGS WMS for UK soil data using GetFeatureInfo
        
        Uses the working UKSO_BGS endpoint
        """
        self._rate_limit()
        
        # Try using OWSLib first
        if OWSLIB_AVAILABLE:
            try:
                wms = WebMapService(
                    f"{self.bgs_wms_url}?service=WMS&version=1.3.0",
                    version='1.3.0',
                    timeout=self.config.api.request_timeout
                )
                
                # List available layers
                available_layers = list(wms.contents.keys())
                logger.debug(f"BGS WMS available layers: {available_layers[:10]}")  # Log first 10
                
                # Try to find soil-related layers
                # BGS layers might be named like: "0", "1", "2" or have descriptive names
                soil_layer_candidates = []
                for layer_name in available_layers:
                    layer_lower = layer_name.lower()
                    if any(term in layer_lower for term in ['soil', 'texture', 'parent', 'material', 'ukso']):
                        soil_layer_candidates.append(layer_name)
                
                # If no obvious soil layers, try first few layers
                if not soil_layer_candidates and available_layers:
                    soil_layer_candidates = available_layers[:3]
                
                for layer_name in soil_layer_candidates:
                    try:
                        # Get feature info at point
                        # Use small bbox around the point
                        bbox = (lon - 0.0001, lat - 0.0001, lon + 0.0001, lat + 0.0001)
                        
                        # Try JSON first
                        try:
                            info = wms.getfeatureinfo(
                                layers=[layer_name],
                                srs='EPSG:4326',
                                bbox=bbox,
                                size=(1, 1),
                                xy=(0, 0),
                                info_format='application/json'
                            )
                            data = json.loads(info.read())
                            if data and isinstance(data, dict):
                                logger.debug(f"BGS WMS JSON response for {layer_name}: {list(data.keys())[:5]}")
                                # Try to extract soil properties
                                return self._parse_bgs_response(data)
                        except:
                            # Try text/plain format
                            try:
                                info = wms.getfeatureinfo(
                                    layers=[layer_name],
                                    srs='EPSG:4326',
                                    bbox=bbox,
                                    size=(1, 1),
                                    xy=(0, 0),
                                    info_format='text/plain'
                                )
                                text_data = info.read().decode('utf-8')
                                logger.debug(f"BGS WMS text response for {layer_name}: {text_data[:200]}")
                                # Try to parse text response
                                return self._parse_bgs_text_response(text_data)
                            except Exception as e2:
                                logger.debug(f"BGS WMS GetFeatureInfo failed for {layer_name}: {e2}")
                                continue
                                
                    except Exception as e:
                        logger.debug(f"BGS WMS query failed for {layer_name}: {e}")
                        continue
                        
            except Exception as e:
                logger.debug(f"OWSLib BGS WMS connection failed: {e}, trying direct request")
        
        # Fallback to direct HTTP request
        try:
            # First, get capabilities to find layer names
            caps_params = {
                "SERVICE": "WMS",
                "VERSION": "1.3.0",
                "REQUEST": "GetCapabilities"
            }
            caps_response = requests.get(
                self.bgs_wms_url,
                params=caps_params,
                timeout=self.config.api.request_timeout
            )
            
            if caps_response.status_code == 200:
                # Try to find a layer name from capabilities (simplified)
                # For now, try common layer names or layer "0"
                layer_to_try = "0"
                
                # Try GetFeatureInfo
                params = {
                    "SERVICE": "WMS",
                    "VERSION": "1.3.0",
                    "REQUEST": "GetFeatureInfo",
                    "LAYERS": layer_to_try,
                    "QUERY_LAYERS": layer_to_try,
                    "CRS": "EPSG:4326",
                    "BBOX": f"{lon-0.0001},{lat-0.0001},{lon+0.0001},{lat+0.0001}",
                    "I": "0",
                    "J": "0",
                    "WIDTH": "1",
                    "HEIGHT": "1",
                    "INFO_FORMAT": "application/json"
                }
                
                response = requests.get(
                    self.bgs_wms_url,
                    params=params,
                    headers={"User-Agent": self.config.api.user_agent},
                    timeout=self.config.api.request_timeout
                )
                
                if response.status_code == 200:
                    try:
                        data = response.json()
                        return self._parse_bgs_response(data)
                    except:
                        # Try parsing as text
                        return self._parse_bgs_text_response(response.text)
                        
        except Exception as e:
            logger.debug(f"Direct BGS WMS query failed: {e}")
        
        return None
    
    def _parse_bgs_response(self, data: dict) -> Optional[Dict[str, Any]]:
        """Parse BGS WMS JSON response and extract soil properties"""
        results = {}
        
        # BGS response structure may vary - try common patterns
        # Look for soil texture values (clay, sand, silt percentages)
        if isinstance(data, dict):
            # Try to find numeric values that might be soil properties
            for key, value in data.items():
                key_lower = str(key).lower()
                if 'clay' in key_lower and isinstance(value, (int, float)):
                    results['clay'] = float(value) * 10  # Convert % to g/kg if needed
                elif 'sand' in key_lower and isinstance(value, (int, float)):
                    results['sand'] = float(value) * 10
                elif 'silt' in key_lower and isinstance(value, (int, float)):
                    results['silt'] = float(value) * 10
                elif isinstance(value, dict):
                    # Recursively search nested dicts
                    nested = self._parse_bgs_response(value)
                    if nested:
                        results.update(nested)
        
        return results if results else None
    
    def _parse_bgs_text_response(self, text: str) -> Optional[Dict[str, Any]]:
        """Parse BGS WMS text/plain response"""
        results = {}
        
        # Try to extract key-value pairs from text
        lines = text.split('\n')
        for line in lines:
            line = line.strip()
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip().lower()
                value = value.strip()
                
                # Try to parse as number
                try:
                    num_val = float(value)
                    if 'clay' in key:
                        results['clay'] = num_val * 10
                    elif 'sand' in key:
                        results['sand'] = num_val * 10
                    elif 'silt' in key:
                        results['silt'] = num_val * 10
                except:
                    pass
        
        return results if results else None
    
    def _classify_soil_type(self, properties: Dict[str, Any]) -> str:
        """Classify soil type based on texture (clay/sand/silt proportions)"""
        
        # Check if we have any texture data
        if not properties or not any(key in properties for key in ["clay", "sand", "silt"]):
            return "unknown"
        
        # SoilGrids returns values in g/kg, need to convert to percentages
        clay = properties.get("clay", 0) / 10  # g/kg to %
        sand = properties.get("sand", 0) / 10
        silt = properties.get("silt", 0) / 10
        
        # If all values are zero or missing, return unknown
        if clay == 0 and sand == 0 and silt == 0:
            return "unknown"
        
        # USDA soil texture triangle classification (simplified)
        if clay >= 40:
            if sand >= 45:
                return "sandy_clay"
            elif silt >= 40:
                return "silty_clay"
            else:
                return "clay_heavy"
        elif clay >= 27:
            if sand >= 20 and sand < 45:
                return "clay_loam"
            elif silt >= 40:
                return "silty_clay_loam"
            else:
                return "sandy_clay_loam"
        elif clay >= 12:
            if silt >= 50:
                return "silt_loam"
            elif sand >= 50:
                return "sandy_loam"
            else:
                return "loam"
        else:
            if silt >= 80:
                return "silt"
            elif sand >= 70:
                return "sand"
            else:
                return "loamy_sand"
    
    def _estimate_drainage(self, properties: Dict[str, Any]) -> str:
        """Estimate soil drainage based on texture"""
        clay = properties.get("clay", 0) / 10  # g/kg to %
        sand = properties.get("sand", 0) / 10
        
        if sand >= 60:
            return "excellent"
        elif sand >= 40:
            return "good"
        elif clay >= 40:
            return "poor"
        elif clay >= 25:
            return "moderate"
        else:
            return "good"
    
    def _estimate_bearing_capacity(self, properties: Dict[str, Any]) -> float:
        """
        Estimate soil bearing capacity in kPa
        
        This is a rough estimate - actual values require geotechnical surveys
        """
        clay = properties.get("clay", 0) / 10  # %
        sand = properties.get("sand", 0) / 10
        
        # Bulk density (g/cm³) - SoilGrids returns in cg/cm³
        bdod = properties.get("bdod", 1400) / 100
        
        # Base bearing capacity estimates (kPa) by soil type
        # These are conservative estimates for residential construction
        
        if sand >= 60:
            # Sandy soils - generally good bearing
            base_capacity = 200
        elif clay >= 40:
            # Heavy clay - can be problematic
            base_capacity = 100
        elif clay >= 25:
            # Clay loam
            base_capacity = 150
        else:
            # Loam, silt loam
            base_capacity = 150
        
        # Adjust for bulk density (higher density = more compact = better bearing)
        density_factor = bdod / 1.4  # normalized to typical value
        
        bearing_capacity = base_capacity * min(1.2, max(0.8, density_factor))
        
        return round(bearing_capacity, 0)
    
    def _recommend_foundation(
        self,
        soil_type: str,
        drainage: str,
        bearing_capacity: float
    ) -> str:
        """Recommend foundation type based on soil conditions"""
        
        # Unknown soil type requires site investigation
        if soil_type == "unknown" or drainage == "unknown":
            return "site_investigation_required"
        
        # High-risk soil types
        if soil_type in ["clay_heavy", "silty_clay"] and drainage == "poor":
            return "deep_pile_or_raft"
        
        if soil_type in ["clay_heavy", "silty_clay"]:
            return "reinforced_strip_or_raft"
        
        if bearing_capacity < 100:
            return "engineered_slab_or_pile"
        
        if drainage == "poor":
            return "raised_slab_with_drainage"
        
        if bearing_capacity >= 200 and drainage in ["excellent", "good"]:
            return "standard_strip_foundation"
        
        # Default for moderate conditions
        return "strip_foundation_with_concrete_slab"
    
    def _get_unknown_soil(self) -> Dict[str, Any]:
        """Return unknown soil type when no data is available"""
        return {
            "type_id": "unknown",
            "bearing_capacity_kpa": 100,  # Conservative default - requires site investigation
            "drainage": "unknown",
            "foundation_recommendation": "site_investigation_required",
            "source": "unknown",
            "properties": {}  # No properties available
        }

