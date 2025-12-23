"""
Soil data collector using SoilGrids API
Collects soil properties for foundation recommendations
"""

import time
from typing import Dict, Any, Optional
import requests
from loguru import logger

from ..config import get_config


class SoilCollector:
    """Collect soil data from SoilGrids API (ISRIC)"""
    
    def __init__(self):
        self.config = get_config()
        self.base_url = self.config.api.soilgrids_url
        self._last_request_time = 0
        self._min_request_interval = 0.5
    
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
        """
        logger.info(f"Fetching soil data for ({lat}, {lon})")
        
        # Query SoilGrids for key properties
        properties = self._query_soilgrids(lat, lon)
        
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
                "source": "soilgrids_isric",
                "properties": properties
            }
        
        # Return defaults for UK if API fails
        logger.warning("SoilGrids API failed, using UK defaults")
        return self._get_uk_default_soil()
    
    def _query_soilgrids(self, lat: float, lon: float) -> Optional[Dict[str, Any]]:
        """Query SoilGrids API for soil properties"""
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
                self.base_url,
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
            
            logger.debug(f"SoilGrids results: {results}")
            return results if results else None
            
        except Exception as e:
            logger.warning(f"SoilGrids API request failed: {e}")
            return None
    
    def _classify_soil_type(self, properties: Dict[str, Any]) -> str:
        """Classify soil type based on texture (clay/sand/silt proportions)"""
        
        # SoilGrids returns values in g/kg, need to convert to percentages
        clay = properties.get("clay", 0) / 10  # g/kg to %
        sand = properties.get("sand", 0) / 10
        silt = properties.get("silt", 0) / 10
        
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
    
    def _get_uk_default_soil(self) -> Dict[str, Any]:
        """Return default UK soil characteristics (London Clay typical)"""
        return {
            "type_id": "clay_loam",
            "bearing_capacity_kpa": 150,
            "drainage": "moderate",
            "foundation_recommendation": "strip_foundation_with_concrete_slab",
            "source": "uk_default",
            "properties": {
                "clay": 300,  # 30%
                "sand": 350,  # 35%
                "silt": 350,  # 35%
            }
        }

