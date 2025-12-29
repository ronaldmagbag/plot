# soil.py
import os
import json
import requests
from pathlib import Path

import geopandas as gpd
import rasterio
from shapely.geometry import Point

# =========================
# Paths
# =========================

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
DATA_DIR.mkdir(exist_ok=True)

UKTS = {
    "clay": DATA_DIR / "ukts_clay.tif",
    "sand": DATA_DIR / "ukts_sand.tif",
    "silt": DATA_DIR / "ukts_silt.tif",
}

SOIL_POLYGONS = DATA_DIR / "soil_parent_material.gpkg"

# =========================
# Download helpers
# =========================

def download_file(url: str, dst: Path):
    if dst.exists():
        return
    print(f"Downloading {dst.name}")
    r = requests.get(url, stream=True, timeout=60)
    r.raise_for_status()
    with open(dst, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)


def download_datasets():
    """
    URLs intentionally isolated here so you can
    swap mirrors or host internally later.
    """

    # UK Compiled Topsoil (UKTS) – example mirrors
    download_file(
        "https://example-mirror.uk/ukts_clay_percent.tif",
        UKTS["clay"]
    )
    download_file(
        "https://example-mirror.uk/ukts_sand_percent.tif",
        UKTS["sand"]
    )
    download_file(
        "https://example-mirror.uk/ukts_silt_percent.tif",
        UKTS["silt"]
    )

    # Soil Parent Material (GPKG)
    download_file(
        "https://example-mirror.uk/soil_parent_material.gpkg",
        SOIL_POLYGONS
    )

# =========================
# Load datasets (lazy)
# =========================

_SOIL_GDF = None
_RASTERS = {}

def load_datasets():
    global _SOIL_GDF, _RASTERS

    if _SOIL_GDF is None and SOIL_POLYGONS.exists():
        _SOIL_GDF = gpd.read_file(SOIL_POLYGONS).to_crs("EPSG:4326")

    for k, path in UKTS.items():
        if k not in _RASTERS and path.exists():
            _RASTERS[k] = rasterio.open(path)

# =========================
# Sampling
# =========================

def sample_raster(raster, lon, lat):
    try:
        for v in raster.sample([(lon, lat)]):
            return float(v[0])
    except Exception:
        return None


def lookup_polygon(lat, lon):
    if _SOIL_GDF is None:
        return None
    pt = Point(lon, lat)
    hit = _SOIL_GDF[_SOIL_GDF.contains(pt)]
    if hit.empty:
        return None
    return hit.iloc[0].to_dict()

# =========================
# Soil classification
# =========================

def classify_texture(clay, sand, silt):
    if clay is None:
        return "unknown"
    if clay > 35:
        return "clay"
    if clay > 20 and sand > 45:
        return "clay_loam"
    if sand > 70:
        return "sandy"
    return "loam"


ENGINEERING_TABLE = {
    "clay": {
        "bearing": 100,
        "drainage": "poor",
        "foundation": "raft_or_piled"
    },
    "clay_loam": {
        "bearing": 150,
        "drainage": "moderate",
        "foundation": "strip_foundation_with_concrete_slab"
    },
    "loam": {
        "bearing": 200,
        "drainage": "good",
        "foundation": "strip_foundation"
    },
    "sandy": {
        "bearing": 250,
        "drainage": "excellent",
        "foundation": "strip_foundation"
    },
    "unknown": {
        "bearing": 100,
        "drainage": "unknown",
        "foundation": "site_investigation_required"
    }
}

# =========================
# Public API
# =========================

def resolve_soil(lat: float, lon: float) -> dict:
    load_datasets()

    clay = sample_raster(_RASTERS.get("clay"), lon, lat)
    sand = sample_raster(_RASTERS.get("sand"), lon, lat)
    silt = sample_raster(_RASTERS.get("silt"), lon, lat)

    soil_type = classify_texture(clay, sand, silt)
    eng = ENGINEERING_TABLE[soil_type]

    return {
        "soil": {
            "type_id": soil_type,
            "bearing_capacity_kpa": eng["bearing"],
            "drainage": eng["drainage"],
            "foundation_recommendation": eng["foundation"],
            "source": "bgs_ukso"
        }
    }
