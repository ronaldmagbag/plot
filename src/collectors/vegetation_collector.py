"""
Vegetation and Tree Zone Collector

Implements Layer 4 from Pipeline:
- SAM/SAM2 for semantic segmentation (optional)
- ESA WorldCover / CORINE as fallback
- RLE encoding for tree_zones masks
"""

import math
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass
import requests
from loguru import logger
import numpy as np

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

from ..config import get_config


@dataclass
class VegetationData:
    """Vegetation/tree zone data"""
    trees: List[Dict[str, Any]]
    zones: List[Dict[str, Any]]
    coverage_percent: float
    rle_data: Optional[str]
    resolution_m: float
    bounds: Dict[str, Any]
    source: str


class VegetationCollector:
    """
    Collect vegetation/tree data from multiple sources:
    
    1. OSM trees (node["natural"="tree"]) - point data
    2. OSM vegetation zones (landuse, forest, etc.) - polygon data
    3. ESA WorldCover 10m (fallback for coverage)
    4. SAM segmentation (optional, requires torch)
    """
    
    # ESA WorldCover API (free, global 10m land cover)
    ESA_WORLDCOVER_URL = "https://services.terrascope.be/wms/v2"
    
    def __init__(self):
        self.config = get_config()
        self._sam_model = None
    
    def collect_vegetation(
        self,
        lat: float,
        lon: float,
        radius_m: float = 100,
        trees_from_osm: Optional[List[Dict]] = None,
        zones_from_osm: Optional[List[Dict]] = None
    ) -> VegetationData:
        """
        Collect vegetation data for plot area
        
        Args:
            lat: Center latitude
            lon: Center longitude
            radius_m: Radius to analyze
            trees_from_osm: Pre-collected OSM tree data
            zones_from_osm: Pre-collected OSM vegetation zones
            
        Returns:
            VegetationData with trees, zones, coverage, and RLE mask
        """
        logger.info(f"Collecting vegetation data for ({lat}, {lon})")
        
        trees = trees_from_osm or []
        zones = zones_from_osm or []
        
        # Calculate bounds
        offset_deg = radius_m / 111000
        bounds = {
            "type": "Polygon",
            "coordinates": [[
                [lon - offset_deg, lat - offset_deg],
                [lon + offset_deg, lat - offset_deg],
                [lon + offset_deg, lat + offset_deg],
                [lon - offset_deg, lat + offset_deg],
                [lon - offset_deg, lat - offset_deg]
            ]]
        }
        
        # Calculate coverage from trees and zones
        coverage_percent = self._calculate_coverage(
            trees, zones, lat, lon, radius_m
        )
        
        # Generate RLE encoding if we have data
        rle_data = self._generate_rle_mask(
            trees, zones, lat, lon, radius_m
        )
        
        return VegetationData(
            trees=trees,
            zones=zones,
            coverage_percent=coverage_percent,
            rle_data=rle_data,
            resolution_m=1.0,  # 1m grid for RLE
            bounds=bounds,
            source="osm_vegetation"
        )
    
    def _calculate_coverage(
        self,
        trees: List[Dict],
        zones: List[Dict],
        lat: float,
        lon: float,
        radius_m: float
    ) -> float:
        """
        Calculate vegetation coverage percentage
        
        Uses tree canopy estimation + zone areas
        """
        total_area = math.pi * radius_m * radius_m  # Circular area
        
        tree_area = 0.0
        for tree in trees:
            # Estimate canopy area from height (radius ≈ height * 0.4)
            height = tree.get("height_m", 8)
            canopy_radius = height * 0.4
            tree_area += math.pi * canopy_radius * canopy_radius
        
        zone_area = 0.0
        for zone in zones:
            # For now, estimate from typical zone sizes
            zone_area += 100  # 100 sqm per zone as estimate
        
        # Avoid double-counting (trees might be in zones)
        vegetation_area = min(tree_area + zone_area * 0.5, total_area)
        
        return (vegetation_area / total_area) * 100
    
    def _generate_rle_mask(
        self,
        trees: List[Dict],
        zones: List[Dict],
        lat: float,
        lon: float,
        radius_m: float,
        resolution_m: float = 1.0
    ) -> Optional[str]:
        """
        Generate RLE-encoded vegetation mask
        
        RLE format: "count1:value1,count2:value2,..."
        value 0 = no vegetation, 1 = vegetation
        """
        if not trees and not zones:
            # No vegetation - return empty mask indicator
            return "0x00000000"
        
        # Create binary grid
        grid_size = int(2 * radius_m / resolution_m)
        grid = np.zeros((grid_size, grid_size), dtype=np.uint8)
        
        meters_per_deg = 111000 * math.cos(math.radians(lat))
        
        # Mark tree canopies
        for tree in trees:
            loc = tree.get("location", [])
            if len(loc) == 2:
                tree_lon, tree_lat = loc
                
                # Convert to grid coordinates
                dx = (tree_lon - lon) * meters_per_deg + radius_m
                dy = (tree_lat - lat) * 111000 + radius_m
                
                gx = int(dx / resolution_m)
                gy = int(dy / resolution_m)
                
                # Draw canopy circle
                height = tree.get("height_m", 8)
                canopy_radius = int(height * 0.4 / resolution_m)
                
                for i in range(-canopy_radius, canopy_radius + 1):
                    for j in range(-canopy_radius, canopy_radius + 1):
                        if i*i + j*j <= canopy_radius*canopy_radius:
                            px, py = gx + i, gy + j
                            if 0 <= px < grid_size and 0 <= py < grid_size:
                                grid[py, px] = 1
        
        # Mark vegetation zones (simplified - just mark approximate area)
        for zone in zones:
            geom = zone.get("geometry", {})
            coords = geom.get("coordinates", [[]])
            if coords and coords[0]:
                # Get zone center
                zone_coords = coords[0]
                if zone_coords:
                    center_lon = sum(c[0] for c in zone_coords) / len(zone_coords)
                    center_lat = sum(c[1] for c in zone_coords) / len(zone_coords)
                    
                    dx = (center_lon - lon) * meters_per_deg + radius_m
                    dy = (center_lat - lat) * 111000 + radius_m
                    
                    gx = int(dx / resolution_m)
                    gy = int(dy / resolution_m)
                    
                    # Mark small area around center
                    for i in range(-5, 6):
                        for j in range(-5, 6):
                            px, py = gx + i, gy + j
                            if 0 <= px < grid_size and 0 <= py < grid_size:
                                grid[py, px] = 1
        
        # RLE encode the flattened grid
        flat = grid.flatten()
        rle_parts = []
        
        if len(flat) == 0:
            return "0x00000000"
        
        current_val = flat[0]
        count = 1
        
        for val in flat[1:]:
            if val == current_val:
                count += 1
            else:
                rle_parts.append(f"{count}:{current_val}")
                current_val = val
                count = 1
        
        rle_parts.append(f"{count}:{current_val}")
        
        # Encode as hex-like string
        rle_str = ",".join(rle_parts)
        
        # If too long, compress further
        if len(rle_str) > 1000:
            total_veg = np.sum(grid)
            total_cells = grid.size
            return f"0x{total_veg:08X}_{total_cells:08X}"
        
        return f"rle:{rle_str}"
    
    def try_sam_segmentation(
        self,
        image_path: str,
        lat: float,
        lon: float
    ) -> Optional[np.ndarray]:
        """
        Try to segment vegetation using SAM model
        
        Requires: segment-anything and torch
        
        Returns:
            Binary mask of vegetation or None if SAM unavailable
        """
        try:
            from segment_anything import SamPredictor, sam_model_registry
            import torch
        except ImportError:
            logger.warning("SAM not available - install segment-anything and torch")
            return None
        
        if not HAS_PIL:
            logger.warning("PIL required for SAM segmentation")
            return None
        
        try:
            # Load image
            image = np.array(Image.open(image_path))
            
            # Load SAM model (first time only)
            if self._sam_model is None:
                # Use smallest model for speed
                sam = sam_model_registry["vit_b"](checkpoint="sam_vit_b.pth")
                if torch.cuda.is_available():
                    sam = sam.to("cuda")
                self._sam_model = SamPredictor(sam)
            
            self._sam_model.set_image(image)
            
            # Prompt for vegetation (green areas)
            # This is simplified - real implementation would use better prompts
            h, w = image.shape[:2]
            
            # Generate mask for central area
            masks, _, _ = self._sam_model.predict(
                point_coords=np.array([[w//2, h//2]]),
                point_labels=np.array([1]),
                multimask_output=True
            )
            
            return masks[0] if len(masks) > 0 else None
            
        except Exception as e:
            logger.error(f"SAM segmentation failed: {e}")
            return None
    
    def get_esa_worldcover(
        self,
        lat: float,
        lon: float,
        radius_m: float = 100
    ) -> Optional[Dict[str, float]]:
        """
        Get land cover percentages from ESA WorldCover
        
        WorldCover classes:
        - 10: Tree cover
        - 20: Shrubland  
        - 30: Grassland
        - 40: Cropland
        - 50: Built-up
        - 60: Bare/sparse vegetation
        - 70: Snow/Ice
        - 80: Water
        - 90: Herbaceous wetland
        - 95: Mangroves
        - 100: Moss/lichen
        
        Returns:
            Dict with class percentages or None if unavailable
        """
        try:
            # Calculate bounding box
            offset = radius_m / 111000
            bbox = f"{lon-offset},{lat-offset},{lon+offset},{lat+offset}"
            
            # Request WMS GetFeatureInfo (simplified)
            # Full implementation would use WCS or download tiles
            logger.info("ESA WorldCover lookup would require WCS implementation")
            
            return None  # Placeholder
            
        except Exception as e:
            logger.error(f"ESA WorldCover lookup failed: {e}")
            return None

