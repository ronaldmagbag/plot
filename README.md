# Plot Analysis Generator v2.0

Generate comprehensive plot analysis JSON files from center coordinates for UK residential properties.

## Architecture Overview

This pipeline follows an **11-step automated flow**:

```
Input: Lat/Lon or Parcel ID
  ↓
1. Fetch imagery + DEM (Mapbox Terrain-RGB)
  ↓
2. Fetch parcels, buildings, roads (OSM + Cadastre)
  ↓
3. Segment trees & water (SAM or landcover fallback)
  ↓
4. Compute terrain & slope
  ↓
5. Compute building heights (OSM tags / LiDAR)
  ↓
6. Compute shadows (pvlib)
  ↓
7. Apply zoning + setbacks (rule engine)
  ↓
8. Detect access (road adjacency)
  ↓
9. Assemble plot.json
```

## Features

### Data Acquisition (Layer 2)
- **Terrain/DEM**: Mapbox Terrain-RGB tiles, Open-Meteo fallback
- **Buildings**: OSM with height extraction from `height` and `building:levels` tags
- **Roads**: OSM highway network with width estimation
- **Vegetation**: OSM trees + SAM segmentation (optional) with RLE encoding
- **Water**: OSM water features
- **Soil**: SoilGrids ISRIC API, British Geological Survey

### Geometric Processing (Layer 3)
- Shapely/GeoPandas for polygon operations
- Centroid, area, perimeter calculations
- 50m buffer for context
- Edge iteration for adjacency

### Semantic Understanding (Layer 4)
- Tree zone detection with canopy estimation
- RLE-encoded vegetation masks
- ESA WorldCover fallback

### Shadow Analysis (Layer 6)
- pvlib for sun position calculations
- Ray-casting for shadow computation
- Facade solar scores (N/S/E/W)

### Setbacks & Buildable Envelope (Layer 7)
- Declarative zoning rules
- Polygon offset operations
- Constraint application (easements, utilities)

### Access Detection (Layer 8)
- Street adjacency detection
- Vehicle/pedestrian access analysis
- Driveway width requirements

## Installation

### Using micromamba (recommended)

```bash
# Create and activate environment
micromamba create -n plot-generator python=3.11
micromamba activate plot-generator

# Install dependencies
pip install -r requirements.txt

# Optional: Install GDAL for advanced raster processing
micromamba install gdal -c conda-forge

# Optional: Install SAM for vegetation segmentation
pip install segment-anything torch
```

### Environment Variables

```bash
# Optional: For Mapbox terrain tiles
export MAPBOX_ACCESS_TOKEN="your_mapbox_token"

# Optional: For Sentinel Hub imagery
export SENTINELHUB_CLIENT_ID="your_client_id"
export SENTINELHUB_CLIENT_SECRET="your_client_secret"
```

## Quick Start

### Generate a Single Plot

```python
from src.pipeline import PlotAnalysisPipeline

pipeline = PlotAnalysisPipeline()

# Generate plot analysis for a location in South UK
result = pipeline.run(
    lat=51.3762,  # Croydon
    lon=-0.0982,
    plot_id="MY_PLOT_001"
)

# Save to file
pipeline.save(result, "output/my_plot.json")
```

### Command Line Interface

```bash
# Generate single plot
python cli.py generate --lat 51.268535 --lon -0.570979 --output plot_gu47nn.json

# Generate 20 South UK samples
python tests/generate_samples.py --count 20 --output ./output

# Visualize a plot
python cli.py visualize --input output/plot.json --output plot.png
```

## Output JSON Structure

The generated `plot.json` follows the Plot Analysis Data Structure:

```json
{
  "plot_id": "UKS-GB-20241218-abc1",
  "centroid": {"type": "Point", "coordinates": [-0.0982, 51.3762]},
  "plot_type": "residential_building",
  
  "boundaries": {
    "property_line": {...},
    "setback_line": {...},
    "buildable_envelope": {...}
  },
  
  "surrounding_context": {
    "buildings": [...],
    "roads": [...],
    "tree_zones": {
      "encoding": "rle",
      "data": "rle:...",
      "trees": [...],
      "coverage_percent": 12.5
    },
    "elevation_map": {...},
    "water_features": {...}
  },
  
  "analysis": {
    "shadow_analysis": {...},
    "adjacency": [...]
  },
  
  "access": {...},
  "existing_structures": [...],
  "regulatory": {...},
  "soil": {...},
  "data_quality": {...}
}
```

## Data Sources

| Data Type | Primary Source | Fallback | API |
|-----------|---------------|----------|-----|
| **Terrain/DEM** | Mapbox Terrain-RGB | Open-Meteo | Free (Mapbox requires token) |
| **Buildings** | OpenStreetMap | - | Overpass API |
| **Roads** | OpenStreetMap | - | Overpass API |
| **Vegetation** | SAM Segmentation | OSM + WorldCover | Local + Overpass |
| **Water** | OpenStreetMap | - | Overpass API |
| **Soil** | SoilGrids (ISRIC) | UK defaults | Free API |
| **Parcels** | UK Land Registry | OSM estimate | WFS (planned) |

## Project Structure

```
plot/
├── src/
│   ├── __init__.py
│   ├── config.py              # Configuration settings
│   ├── models.py              # Pydantic data models
│   ├── pipeline.py            # Main 11-step orchestrator
│   ├── collectors/
│   │   ├── osm_collector.py        # OSM buildings, roads, water, vegetation
│   │   ├── terrain_collector.py    # Mapbox Terrain-RGB DEM
│   │   ├── vegetation_collector.py # Tree zones + RLE encoding
│   │   ├── elevation_collector.py  # Open-Meteo elevation fallback
│   │   ├── soil_collector.py       # SoilGrids API
│   │   └── boundary_collector.py   # Boundary detection
│   └── analysis/
│       ├── shadow_analyzer.py      # pvlib sun/shadow analysis
│       ├── adjacency_analyzer.py   # Edge adjacency detection
│       ├── setback_calculator.py   # Zoning rule engine
│       └── geometry_utils.py       # Shapely geometry helpers
├── cli.py                     # Command-line interface
├── requirements.txt
├── tests/
│   ├── generate_samples.py        # Batch sample generator
│   ├── generate_residential_samples.py  # Residential sample generator
│   ├── debug_classifier.py      # Property line classifier debug tool
│   └── test_*.py                # Unit tests
├── docs/
│   ├── Plot Analysis Data Structure.md
│   └── plot.json              # Example output
└── README.md
```

## UK Building Regulations Applied

Default setbacks (configurable):
- **Front setback**: 5m from street
- **Rear setback**: 5m from rear boundary
- **Side setbacks**: 1m from side boundaries

Default constraints:
- **Max coverage**: 40% of plot area
- **Max height**: 10m
- **Max stories**: 2 (plus loft conversion)

## Sample Locations (South UK)

The sample generator includes 20 residential locations across South UK:

- Greater London suburbs (Croydon, Bromley, Richmond, Wimbledon, Kingston)
- Surrey (Guildford, Woking, Epsom)
- Kent (Sevenoaks, Tunbridge Wells, Maidstone)
- Sussex (Brighton, Worthing, Horsham)
- Hampshire (Southampton, Winchester, Basingstoke)
- Berkshire (Reading, Windsor)
- Oxfordshire (Oxford suburb)

## API Rate Limits

The generator respects API rate limits:
- Overpass API: 1 request per second
- Open-Meteo: No strict limit
- SoilGrids: 0.5 seconds between requests
- Mapbox: 100,000 requests/month (free tier)

Default delay between samples: 3 seconds

## Optional ML Enhancements (Future)

| Task | Model |
|------|-------|
| Missing building heights | Vision-based depth models |
| Privacy scoring | Graph neural networks |
| Optimal building placement | Reinforcement learning |
| Design massing | Diffusion + constraints |

## Requirements

- Python 3.9+
- Internet connection for API access
- Optional: MAPBOX_ACCESS_TOKEN for terrain tiles
- Optional: CUDA GPU for SAM segmentation

## License

MIT License
