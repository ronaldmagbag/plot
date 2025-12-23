# Plot Analysis Data Generation Pipeline

## Overview

This document describes the complete pipeline to generate plot analysis JSON from a **center coordinate (lat/lon)**.

---

## 🎯 Required Inputs

| Input | Format | Example | Required |
|-------|--------|---------|----------|
| **Center Coordinate** | `[longitude, latitude]` | `[-0.1276, 51.5074]` | ✅ Yes |
| **Search Radius** | meters | `50` | ✅ Yes |
| **Country/Region** | ISO code | `GB`, `US`, `DE` | ✅ Yes |
| **Plot ID** (optional) | string | Cadastral reference | Optional |

---

## 🔄 Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           INPUT: Center Coordinate                          │
│                              [-0.1276, 51.5074]                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        STAGE 1: BOUNDARY DETECTION                          │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────────────────┐  │
│  │ Cadastral APIs  │  │ OpenStreetMap   │  │ Satellite Segmentation      │  │
│  │ (Land Registry) │  │ (Overpass API)  │  │ (SAM2 / DeepLab)            │  │
│  └────────┬────────┘  └────────┬────────┘  └─────────────┬───────────────┘  │
│           └────────────────────┴─────────────────────────┘                  │
│                                     │                                       │
│                          Property Line Polygon                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        STAGE 2: GEOSPATIAL DATA                             │
│  ┌───────────────┐ ┌───────────────┐ ┌───────────────┐ ┌───────────────┐    │
│  │   Buildings   │ │    Roads      │ │   Elevation   │ │  Vegetation   │    │
│  │ OSM/MS/Google │ │ OSM Overpass  │ │ SRTM/ALOS/DEM │ │ Sentinel-2    │    │
│  └───────┬───────┘ └───────┬───────┘ └───────┬───────┘ └───────┬───────┘    │
│          │                 │                 │                 │            │
│  ┌───────────────┐ ┌───────────────┐ ┌───────────────┐ ┌───────────────┐    │
│  │    Water      │ │     Soil      │ │   Zoning      │ │  Utilities    │    │
│  │ OSM/Sentinel  │ │ SSURGO/BGS    │ │ Local APIs    │ │ Local Data    │    │
│  └───────┬───────┘ └───────┬───────┘ └───────┬───────┘ └───────┬───────┘    │
│          └─────────────────┴─────────────────┴─────────────────┘            │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        STAGE 3: ANALYSIS & COMPUTATION                      │
│  ┌───────────────────┐  ┌───────────────────┐  ┌───────────────────────┐    │
│  │ Setback Calc      │  │ Shadow Analysis   │  │ Adjacency Detection   │    │
│  │ (Regulation-based)│  │ (Sun Path + DEM)  │  │ (Spatial Joins)       │    │
│  └───────────────────┘  └───────────────────┘  └───────────────────────┘    │
│  ┌───────────────────┐  ┌───────────────────┐  ┌───────────────────────┐    │
│  │ Access Points     │  │ Buildable Area    │  │ Constraint Overlay    │    │
│  │ (Road proximity)  │  │ (Boolean ops)     │  │ (Easements, TPOs)     │    │
│  └───────────────────┘  └───────────────────┘  └───────────────────────┘    │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        STAGE 4: OUTPUT GENERATION                           │
│                                                                             │
│              ┌─────────────────────────────────────────────┐                │
│              │           Plot Analysis JSON                │                │
│              │     (MongoDB-compatible document)           │                │
│              └─────────────────────────────────────────────┘                │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 📦 Stage 1: Boundary Detection

### Option A: Official Cadastral APIs (Best Accuracy)

| Country | API/Source | Access |
|---------|------------|--------|
| 🇬🇧 UK | [HM Land Registry](https://use-land-property-data.service.gov.uk/) | API (paid for bulk) |
| 🇺🇸 USA | County Assessor GIS / [Regrid](https://regrid.com/) | Varies by county |
| 🇩🇪 Germany | [ALKIS](https://www.adv-online.de/) | State-level portals |
| 🇫🇷 France | [cadastre.gouv.fr](https://cadastre.gouv.fr/) | Free API |
| 🇦🇺 Australia | State Land Registries | Varies by state |

**Example UK Request:**
```python
import requests

def get_uk_parcel(lat, lon):
    """Query UK Land Registry INSPIRE Index Polygons"""
    # WFS endpoint for INSPIRE polygons
    url = "https://data.landregistry.gov.uk/data/inspire/index-polygons"
    params = {
        "service": "WFS",
        "request": "GetFeature",
        "typeName": "inspire:index-polygons",
        "outputFormat": "json",
        "bbox": f"{lon-0.001},{lat-0.001},{lon+0.001},{lat+0.001},EPSG:4326"
    }
    response = requests.get(url, params=params)
    return response.json()
```

### Option B: OpenStreetMap (Free, Global)

```python
import overpy

def get_osm_boundaries(lat, lon, radius_m=100):
    """Query OSM for plot boundaries"""
    api = overpy.Overpass()
    
    query = f"""
    [out:json][timeout:30];
    (
      way["landuse"](around:{radius_m},{lat},{lon});
      way["boundary"="cadastral"](around:{radius_m},{lat},{lon});
      relation["boundary"](around:{radius_m},{lat},{lon});
    );
    out body;
    >;
    out skel qt;
    """
    
    result = api.query(query)
    return result
```

### Option C: AI Segmentation (When No Data Available)

Use **SAM2** (Segment Anything Model 2) or similar to detect plot boundaries from satellite imagery:

```python
from segment_anything import SamPredictor, sam_model_registry
import rasterio

def segment_plot_from_satellite(image_path, center_point):
    """Use SAM2 to segment plot boundary from satellite image"""
    sam = sam_model_registry["vit_h"](checkpoint="sam_vit_h.pth")
    predictor = SamPredictor(sam)
    
    # Load satellite image
    with rasterio.open(image_path) as src:
        image = src.read()
    
    predictor.set_image(image)
    masks, scores, _ = predictor.predict(
        point_coords=center_point,
        point_labels=[1],  # foreground
        multimask_output=True
    )
    
    return masks[scores.argmax()]  # Best mask
```

---

## 📦 Stage 2: Geospatial Data Collection

### 2.1 Building Footprints

| Source | Coverage | Quality | Access |
|--------|----------|---------|--------|
| **OpenStreetMap** | Global | Variable | Free (Overpass API) |
| **Microsoft Building Footprints** | Global | High | Free download |
| **Google Open Buildings** | Africa, Asia, Latin America | High | Free download |
| **OS MasterMap** (UK) | UK | Excellent | Licensed |

**OpenStreetMap Buildings:**
```python
def get_buildings(lat, lon, radius_m=50):
    query = f"""
    [out:json][timeout:30];
    way["building"](around:{radius_m},{lat},{lon});
    out body;
    >;
    out skel qt;
    """
    return overpy.Overpass().query(query)
```

**Microsoft Building Footprints:**
```bash
# Download from: https://github.com/microsoft/GlobalMLBuildingFootprints
# Files are organized by region (GeoJSON format)
```

### 2.2 Road Network

```python
def get_roads(lat, lon, radius_m=100):
    query = f"""
    [out:json][timeout:30];
    way["highway"](around:{radius_m},{lat},{lon});
    out body;
    >;
    out skel qt;
    """
    return overpy.Overpass().query(query)
```

**Road width estimation:**
| Road Type | Typical Width (m) |
|-----------|-------------------|
| motorway | 12-15 |
| primary | 8-10 |
| secondary | 6-8 |
| residential | 5-6 |
| service | 3-4 |

### 2.3 Elevation / DEM

| Source | Resolution | Coverage | Access |
|--------|------------|----------|--------|
| **SRTM** | 30m | Global (60°N-56°S) | Free |
| **ALOS World 3D** | 30m | Global | Free (registration) |
| **Copernicus DEM** | 30m / 90m | Global | Free |
| **LIDAR** | <1m | Country-specific | Varies |
| **Mapbox Terrain** | ~5m | Global | API (paid) |

**Using elevation-api (Python):**
```python
import elevation
import rasterio

def get_elevation(bounds, output_path):
    """Download DEM for bounding box"""
    # bounds = (west, south, east, north)
    elevation.clip(bounds=bounds, output=output_path)
    
    with rasterio.open(output_path) as src:
        return src.read(1)
```

**UK Environment Agency LIDAR:**
```python
# Free 1m resolution LIDAR for England
# https://environment.data.gov.uk/DefraDataDownload/?Mode=survey
```

### 2.4 Vegetation / Tree Zones

**Satellite-based detection:**
```python
# Using Sentinel-2 NDVI
import sentinelhub

def get_vegetation_mask(bbox, time_range):
    """Calculate NDVI from Sentinel-2"""
    config = sentinelhub.SHConfig()
    
    evalscript = """
    //VERSION=3
    function setup() {
        return {
            input: ["B04", "B08"],
            output: { bands: 1 }
        };
    }
    function evaluatePixel(sample) {
        let ndvi = (sample.B08 - sample.B04) / (sample.B08 + sample.B04);
        return [ndvi > 0.4 ? 1 : 0];  // Tree threshold
    }
    """
    
    request = sentinelhub.SentinelHubRequest(
        evalscript=evalscript,
        input_data=[sentinelhub.SentinelHubRequest.input_data(
            data_collection=sentinelhub.DataCollection.SENTINEL2_L2A,
            time_interval=time_range
        )],
        responses=[sentinelhub.SentinelHubRequest.output_response('default', 'image/tiff')],
        bbox=bbox,
        size=[512, 512],
        config=config
    )
    
    return request.get_data()[0]
```

### 2.5 Water Features

```python
def get_water_features(lat, lon, radius_m=100):
    query = f"""
    [out:json][timeout:30];
    (
      way["natural"="water"](around:{radius_m},{lat},{lon});
      way["waterway"](around:{radius_m},{lat},{lon});
      relation["natural"="water"](around:{radius_m},{lat},{lon});
    );
    out body;
    >;
    out skel qt;
    """
    return overpy.Overpass().query(query)
```

### 2.6 Soil Data

| Region | Source | Access |
|--------|--------|--------|
| 🇺🇸 USA | [USDA SSURGO](https://websoilsurvey.nrcs.usda.gov/) | Free API |
| 🇬🇧 UK | [British Geological Survey](https://www.bgs.ac.uk/) | API (some free) |
| 🇪🇺 EU | [ESDAC](https://esdac.jrc.ec.europa.eu/) | Free download |
| 🌍 Global | [SoilGrids](https://soilgrids.org/) | Free API |

```python
def get_soil_data(lat, lon):
    """Query SoilGrids for soil properties"""
    import requests
    
    url = f"https://rest.isric.org/soilgrids/v2.0/properties/query"
    params = {
        "lon": lon,
        "lat": lat,
        "property": ["clay", "sand", "silt", "phh2o"],
        "depth": "0-30cm"
    }
    
    response = requests.get(url, params=params)
    return response.json()
```

### 2.7 Zoning & Regulations

| Region | Source | Notes |
|--------|--------|-------|
| 🇺🇸 USA | [Regrid](https://regrid.com/), County GIS | Varies by jurisdiction |
| 🇬🇧 UK | [Planning Portal](https://www.planningportal.co.uk/) | Local authority APIs |
| 🇩🇪 Germany | [XPlanung](https://www.xplanung.de/) | Standardized format |

**UK Planning Data:**
```python
def get_uk_planning_constraints(lat, lon):
    """Query UK planning constraints"""
    # Each local authority has different APIs
    # Example: London DataStore
    url = "https://data.london.gov.uk/api/..."
    # Implementation varies by council
```

---

## 📦 Stage 3: Analysis & Computation

### 3.1 Setback Calculation

```python
from shapely.geometry import Polygon
from shapely.ops import unary_union

def calculate_setbacks(property_polygon, setbacks):
    """
    Calculate setback line from property boundary
    
    setbacks = {
        'front': 5.0,  # meters
        'rear': 5.0,
        'side_left': 1.0,
        'side_right': 1.0
    }
    """
    # Simple uniform buffer (negative = inward)
    avg_setback = sum(setbacks.values()) / len(setbacks)
    setback_polygon = property_polygon.buffer(-avg_setback)
    
    return setback_polygon
```

### 3.2 Shadow Analysis

```python
import pvlib
from shapely.geometry import Polygon
import numpy as np

def calculate_shadows(buildings, plot_centroid, date_times):
    """
    Calculate shadow projections from buildings
    
    Uses pvlib for sun position calculations
    """
    lat, lon = plot_centroid
    
    results = {}
    for dt in date_times:
        # Get sun position
        solar_position = pvlib.solarposition.get_solarposition(
            dt, lat, lon
        )
        
        sun_altitude = solar_position['apparent_elevation'].values[0]
        sun_azimuth = solar_position['azimuth'].values[0]
        
        shadows = []
        for building in buildings:
            height = building['height_m']
            footprint = Polygon(building['footprint'])
            
            # Calculate shadow length
            if sun_altitude > 0:
                shadow_length = height / np.tan(np.radians(sun_altitude))
                # Project shadow in opposite direction of sun
                shadow_direction = (sun_azimuth + 180) % 360
                
                # Create shadow polygon (simplified)
                shadow = project_shadow(footprint, shadow_length, shadow_direction)
                shadows.append(shadow)
        
        results[dt] = unary_union(shadows)
    
    return results
```

### 3.3 Adjacency Analysis

```python
from shapely.geometry import LineString, Polygon
from shapely.ops import nearest_points

def analyze_adjacency(property_polygon, surrounding_features):
    """
    Analyze what each edge of the property is adjacent to
    """
    edges = []
    coords = list(property_polygon.exterior.coords)
    
    for i in range(len(coords) - 1):
        edge = LineString([coords[i], coords[i + 1]])
        
        # Find nearest feature
        nearest_feature = None
        min_distance = float('inf')
        
        for feature in surrounding_features:
            distance = edge.distance(feature['geometry'])
            if distance < min_distance:
                min_distance = distance
                nearest_feature = feature
        
        edges.append({
            'geometry': edge,
            'length_m': edge.length,
            'adjacent_to': nearest_feature['type'] if nearest_feature else 'unknown',
            'distance_m': min_distance
        })
    
    return edges
```

### 3.4 Buildable Envelope

```python
from shapely.ops import unary_union

def calculate_buildable_envelope(setback_polygon, constraints):
    """
    Calculate final buildable area by subtracting all constraints
    
    constraints = [
        {'type': 'easement', 'geometry': Polygon(...)},
        {'type': 'tree_protection', 'geometry': Polygon(...)},
        ...
    ]
    """
    constraint_union = unary_union([c['geometry'] for c in constraints])
    buildable = setback_polygon.difference(constraint_union)
    
    return buildable
```

---

## 🛠️ Open Source Tools & Libraries

### Core Geospatial

| Library | Purpose | Install |
|---------|---------|---------|
| **Shapely** | Geometry operations | `pip install shapely` |
| **GeoPandas** | Spatial dataframes | `pip install geopandas` |
| **Rasterio** | Raster data I/O | `pip install rasterio` |
| **Fiona** | Vector data I/O | `pip install fiona` |
| **PyProj** | Coordinate transforms | `pip install pyproj` |
| **GDAL** | Geospatial Swiss Army knife | `conda install gdal` |

### Data Sources

| Library | Purpose | Install |
|---------|---------|---------|
| **OSMnx** | OpenStreetMap networks | `pip install osmnx` |
| **Overpy** | OSM Overpass API | `pip install overpy` |
| **elevation** | DEM downloading | `pip install elevation` |
| **sentinelhub** | Sentinel satellite data | `pip install sentinelhub` |
| **planetary-computer** | Microsoft data catalog | `pip install planetary-computer` |

### Analysis

| Library | Purpose | Install |
|---------|---------|---------|
| **pvlib** | Solar position/shadows | `pip install pvlib` |
| **suncalc** | Sun/moon calculations | `pip install suncalc` |
| **scikit-image** | Image processing | `pip install scikit-image` |
| **opencv-python** | Computer vision | `pip install opencv-python` |

### AI/ML Segmentation

| Library | Purpose | Install |
|---------|---------|---------|
| **segment-anything** | SAM2 | `pip install segment-anything` |
| **detectron2** | Object detection | [See FB docs] |
| **rasterio** + **torch** | Deep learning on rasters | Standard installs |

---

## 🌐 APIs Reference

### Free APIs

| API | Data | Rate Limits |
|-----|------|-------------|
| **Overpass API** | OSM data | ~10,000 req/day |
| **Nominatim** | Geocoding | 1 req/sec |
| **Open-Elevation** | Elevation | 1000 req/day |
| **SoilGrids** | Soil data | Unlimited |
| **Copernicus** | Satellite imagery | Registration required |

### Paid/Commercial APIs

| API | Data | Pricing |
|-----|------|---------|
| **Google Maps Platform** | Various | Pay-per-use |
| **Mapbox** | Terrain, satellite | Free tier + paid |
| **HERE** | Roads, places | Free tier + paid |
| **Regrid** | US parcels | Subscription |
| **Ordnance Survey** | UK data | Subscription |

---

## 📋 Complete Python Pipeline

```python
"""
Complete Plot Analysis Pipeline
"""

import json
from datetime import datetime
from shapely.geometry import Point, Polygon, shape
from shapely.ops import unary_union
import geopandas as gpd
import overpy
import requests

class PlotAnalysisPipeline:
    def __init__(self, center_lon, center_lat, radius_m=50, country="GB"):
        self.center = Point(center_lon, center_lat)
        self.lon = center_lon
        self.lat = center_lat
        self.radius = radius_m
        self.country = country
        self.api = overpy.Overpass()
        
    def run(self):
        """Execute full pipeline"""
        result = {
            "plot_id": self._generate_plot_id(),
            "created_at": datetime.utcnow().isoformat(),
            "centroid": {
                "type": "Point",
                "coordinates": [self.lon, self.lat]
            }
        }
        
        # Stage 1: Get boundaries
        result["boundaries"] = self._get_boundaries()
        
        # Stage 2: Get surrounding context
        result["surrounding_context"] = {
            "buildings": self._get_buildings(),
            "roads": self._get_roads(),
            "tree_zones": self._get_vegetation(),
            "water_features": self._get_water(),
            "elevation_map": self._get_elevation()
        }
        
        # Stage 3: Analysis
        result["analysis"] = {
            "shadow_analysis": self._analyze_shadows(result),
            "adjacency": self._analyze_adjacency(result)
        }
        
        # Stage 4: Access analysis
        result["access"] = self._analyze_access(result)
        
        return result
    
    def _generate_plot_id(self):
        return f"{self.country}-{datetime.now().strftime('%Y%m%d%H%M%S')}"
    
    def _get_boundaries(self):
        """Get or estimate property boundary"""
        # Try cadastral first, fall back to estimation
        # Implementation depends on country
        pass
    
    def _get_buildings(self):
        """Query OSM for buildings"""
        query = f"""
        [out:json][timeout:30];
        way["building"](around:{self.radius},{self.lat},{self.lon});
        out body;
        >;
        out skel qt;
        """
        result = self.api.query(query)
        
        buildings = []
        for way in result.ways:
            coords = [[n.lon, n.lat] for n in way.nodes]
            buildings.append({
                "id": f"building_{way.id}",
                "footprint": {"type": "Polygon", "coordinates": [coords]},
                "building_type": way.tags.get("building", "unknown"),
                "height_m": self._estimate_height(way.tags)
            })
        
        return buildings
    
    def _estimate_height(self, tags):
        """Estimate building height from OSM tags"""
        if "height" in tags:
            return float(tags["height"].replace("m", ""))
        if "building:levels" in tags:
            return int(tags["building:levels"]) * 3.0
        return 8.0  # Default 2-story
    
    def _get_roads(self):
        """Query OSM for roads"""
        query = f"""
        [out:json][timeout:30];
        way["highway"](around:{self.radius * 2},{self.lat},{self.lon});
        out body;
        >;
        out skel qt;
        """
        result = self.api.query(query)
        
        roads = []
        for way in result.ways:
            coords = [[n.lon, n.lat] for n in way.nodes]
            roads.append({
                "id": f"road_{way.id}",
                "name": way.tags.get("name", "Unnamed Road"),
                "type": way.tags.get("highway", "road"),
                "centerline": {"type": "LineString", "coordinates": coords},
                "width_m": self._estimate_road_width(way.tags.get("highway", ""))
            })
        
        return roads
    
    def _estimate_road_width(self, highway_type):
        """Estimate road width from type"""
        widths = {
            "motorway": 12,
            "trunk": 10,
            "primary": 8,
            "secondary": 7,
            "tertiary": 6,
            "residential": 5,
            "service": 4,
            "path": 2
        }
        return widths.get(highway_type, 5)
    
    def _get_vegetation(self):
        """Get tree/vegetation zones"""
        # Could use NDVI from Sentinel-2 or OSM landuse
        query = f"""
        [out:json][timeout:30];
        (
          way["landuse"="forest"](around:{self.radius},{self.lat},{self.lon});
          way["natural"="wood"](around:{self.radius},{self.lat},{self.lon});
          node["natural"="tree"](around:{self.radius},{self.lat},{self.lon});
        );
        out body;
        >;
        out skel qt;
        """
        return self.api.query(query)
    
    def _get_water(self):
        """Get water features"""
        query = f"""
        [out:json][timeout:30];
        (
          way["natural"="water"](around:{self.radius},{self.lat},{self.lon});
          way["waterway"](around:{self.radius},{self.lat},{self.lon});
        );
        out body;
        >;
        out skel qt;
        """
        return self.api.query(query)
    
    def _get_elevation(self):
        """Get elevation data"""
        # Use open-elevation API
        try:
            url = f"https://api.open-elevation.com/api/v1/lookup"
            params = {"locations": f"{self.lat},{self.lon}"}
            response = requests.get(url, params=params, timeout=10)
            data = response.json()
            return {"elevation_m": data["results"][0]["elevation"]}
        except:
            return {"elevation_m": None}
    
    def _analyze_shadows(self, data):
        """Perform shadow analysis"""
        # Simplified - would use pvlib in production
        return {
            "facade_scores": {
                "north": {"annual_avg": 0.30},
                "south": {"annual_avg": 0.90},
                "east": {"annual_avg": 0.62},
                "west": {"annual_avg": 0.58}
            }
        }
    
    def _analyze_adjacency(self, data):
        """Analyze edge adjacencies"""
        # Would analyze property boundary edges
        return []
    
    def _analyze_access(self, data):
        """Determine access points"""
        # Find nearest road
        return {
            "primary_access_point": {
                "coordinates": [self.lon, self.lat],
                "confidence": "estimated"
            }
        }


# Usage
if __name__ == "__main__":
    pipeline = PlotAnalysisPipeline(
        center_lon=-0.1276,
        center_lat=51.5074,
        radius_m=50,
        country="GB"
    )
    
    result = pipeline.run()
    print(json.dumps(result, indent=2))
```

---

## 📊 Data Quality Matrix

| Data Layer | Best Source | Accuracy | Cost |
|------------|-------------|----------|------|
| Property Boundary | Official Cadastre | ±0.1m | $$ |
| Property Boundary | OSM | ±5m | Free |
| Buildings | OS MasterMap / Local | ±0.5m | $$ |
| Buildings | OSM | ±2m | Free |
| Elevation | LIDAR | ±0.15m | Free-$$ |
| Elevation | SRTM | ±5m | Free |
| Roads | OSM | ±2m | Free |
| Vegetation | Sentinel-2 + ML | ±10m | Free |
| Zoning | Local Authority | Exact | Free-$$ |

---

## 🚀 Quick Start

### Minimal Setup (Free Data Only)

```bash
pip install shapely geopandas overpy requests elevation pvlib matplotlib
```

```python
from pipeline import PlotAnalysisPipeline

result = PlotAnalysisPipeline(
    center_lon=-0.1276,
    center_lat=51.5074
).run()
```

### Production Setup

```bash
pip install shapely geopandas overpy requests rasterio sentinelhub pvlib \
    osmnx segment-anything torch torchvision mongodb pymongo
```

---

## 📝 JSON Schema Issues in Original Document

The example JSON in the data structure document has these issues:

1. **Comments** - Uses `*// comment *` style which is invalid JSON
2. **MongoDB Types** - Uses `ObjectId()` and `ISODate()` which are BSON, not JSON
3. **Trailing commas** - Some sections have trailing commas

**Valid JSON version:**
- Replace `ObjectId("...")` with `"..."`  
- Replace `ISODate("...")` with `"..."`
- Remove all comments or use `"_comment": "..."` pattern
- Remove trailing commas

---

## 📚 Further Resources

- [OpenStreetMap Wiki](https://wiki.openstreetmap.org/)
- [Shapely Documentation](https://shapely.readthedocs.io/)
- [pvlib Documentation](https://pvlib-python.readthedocs.io/)
- [SentinelHub Documentation](https://sentinelhub-py.readthedocs.io/)
- [Microsoft Building Footprints](https://github.com/microsoft/GlobalMLBuildingFootprints)
- [Google Open Buildings](https://sites.research.google/open-buildings/)

