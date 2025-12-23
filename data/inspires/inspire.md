
This is **UK-specific**, handles **EPSG:27700 → EPSG:4326**, and is safe for large GML files.



# 1️⃣ Dependencies

```bash
pip install geopandas shapely pyproj fiona numpy
```


# 4️⃣ `convert.py` – INSPIRE GML → plot.json

```python
import geopandas as gpd
import json
import uuid
from shapely.geometry import mapping
from pyproj import Transformer
from schema import create_plot_schema

INPUT_GML = "INSPIRE_INDEX_POLYGONS.gml"


SRC_EPSG = "EPSG:27700"   # British National Grid
DST_EPSG = "EPSG:4326"    # WGS84


def reproject_geometry(geom):
    transformer = Transformer.from_crs(SRC_EPSG, DST_EPSG, always_xy=True)
    return gpd.GeoSeries([geom], crs=SRC_EPSG).to_crs(DST_EPSG).iloc[0]


def polygon_to_coordinates(geom):
    if geom.geom_type == "Polygon":
        return [list(geom.exterior.coords)]
    elif geom.geom_type == "MultiPolygon":
        # take largest polygon
        largest = max(geom.geoms, key=lambda g: g.area)
        return [list(largest.exterior.coords)]
    else:
        raise ValueError("Unsupported geometry type")


def main():
    gdf = gpd.read_file(INPUT_GML)

    plots = []

    for _, row in gdf.iterrows():
        geom_27700 = row.geometry
        geom_4326 = reproject_geometry(geom_27700)

        plot = create_plot_schema()

        plot["plot_id"] = str(uuid.uuid4())
        plot["source_id"] = row.get("inspire_id") or row.get("INSPIREID")

        plot["geometry"]["coordinates"] = polygon_to_coordinates(geom_4326)

        plot["area_m2"] = float(geom_27700.area)

        centroid = geom_4326.centroid
        plot["centroid"] = [centroid.x, centroid.y]

        # elevation left null (to be enriched later)
        plots.append(plot)

    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(plots, f, indent=2)

    print(f"Generated {len(plots)} plots → {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
