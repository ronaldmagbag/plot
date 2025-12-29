"""
Examine the structure of downloaded UK soil data files
"""
import geopandas as gpd
from pathlib import Path
from shapely.geometry import Point

print("Examining UK Soil Data Files...\n")

# 1. Examine SPMM_1km shapefile
print("=" * 60)
print("1. Soil Parent Material Model (SPMM_1km)")
print("=" * 60)
spmm_path = Path("data/SPMM_1km/SoilParentMateriall_V1_portal1km.shp")

if spmm_path.exists():
    try:
        gdf = gpd.read_file(spmm_path)
        print(f"✓ Loaded shapefile")
        print(f"  CRS: {gdf.crs}")
        print(f"  Shape: {gdf.shape} (rows, cols)")
        print(f"  Columns: {list(gdf.columns)}")
        print(f"\n  Sample data (first row):")
        sample = gdf.iloc[0]
        for col in gdf.columns:
            if col != 'geometry':
                print(f"    {col}: {sample[col]}")
        
        # Test point query (UK location)
        print(f"\n  Testing point query at (53.025, -1.202)...")
        from pyproj import Transformer
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:27700", always_xy=True)
        lon, lat = -1.202, 53.025
        x, y = transformer.transform(lon, lat)
        print(f"    Converted to BNG: ({x:.2f}, {y:.2f})")
        
        point = Point(x, y)
        hits = gdf[gdf.geometry.contains(point) | gdf.geometry.intersects(point)]
        if not hits.empty:
            print(f"    ✓ Found {len(hits)} matching polygon(s)")
            hit = hits.iloc[0]
            print(f"    Soil Texture: {hit.get('SOIL_TEX', 'N/A')}")
            print(f"    Soil Group: {hit.get('SOIL_GROUP', 'N/A')}")
            print(f"    Soil Depth: {hit.get('SOIL_DEPTH', 'N/A')}")
            print(f"    ESB Description: {hit.get('ESB_DESC', 'N/A')}")
        else:
            print(f"    ✗ No matching polygon found")
            
    except Exception as e:
        print(f"✗ Error loading shapefile: {e}")
        import traceback
        traceback.print_exc()
else:
    print(f"✗ File not found: {spmm_path}")

print()

# 2. Examine Geology GPKG
print("=" * 60)
print("2. Geology GPKG (625k_V5_Geology_UK_EPSG27700.gpkg)")
print("=" * 60)
geology_path = Path("data/625k_V5_Geology_UK_EPSG27700.gpkg")

if geology_path.exists():
    try:
        # GPKG might have multiple layers
        import fiona
        layers = fiona.listlayers(str(geology_path))
        print(f"✓ Found {len(layers)} layer(s): {layers}")
        
        for layer_name in layers[:3]:  # Check first 3 layers
            print(f"\n  Layer: {layer_name}")
            gdf = gpd.read_file(geology_path, layer=layer_name)
            print(f"    CRS: {gdf.crs}")
            print(f"    Shape: {gdf.shape}")
            print(f"    Columns: {list(gdf.columns)[:10]}")  # First 10 columns
            if len(gdf) > 0:
                print(f"    Sample row (first column values):")
                sample = gdf.iloc[0]
                for col in list(gdf.columns)[:5]:
                    if col != 'geometry':
                        print(f"      {col}: {sample[col]}")
    except Exception as e:
        print(f"✗ Error loading GPKG: {e}")
        import traceback
        traceback.print_exc()
else:
    print(f"✗ File not found: {geology_path}")

print("\n" + "=" * 60)
print("Summary: Use this information to integrate into soil_collector.py")
print("=" * 60)

