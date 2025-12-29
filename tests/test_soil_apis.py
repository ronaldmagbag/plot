"""
Test script to verify soil data API endpoints
Run this to check which APIs are working
"""
import requests
from pathlib import Path

print("Testing Soil Data APIs...\n")

# Test 1: SoilGrids WCS GetCapabilities
print("1. Testing SoilGrids WCS GetCapabilities...")
wcs_urls = [
    "https://maps.isric.org/mapserv?map=/map/soilgrids.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCapabilities",
    "https://maps.isric.org/mapserv?SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCapabilities",
]
for url in wcs_urls:
    try:
        r = requests.get(url, timeout=10)
        print(f"   URL: {url[:80]}...")
        print(f"   Status: {r.status_code}")
        if r.status_code == 200:
            print(f"   Content-Type: {r.headers.get('Content-Type', 'unknown')}")
            if 'xml' in r.headers.get('Content-Type', '').lower() or 'text/xml' in r.text[:100].lower():
                print("   ✓ WCS GetCapabilities working!")
                break
            else:
                print(f"   ✗ Unexpected response: {r.text[:200]}")
        else:
            print(f"   ✗ Failed: {r.text[:200]}")
    except Exception as e:
        print(f"   ✗ Error: {e}")

print()

# Test 2: SoilGrids REST API
print("2. Testing SoilGrids REST API...")
rest_url = "https://rest.isric.org/soilgrids/v2.0/properties/query"
params = {
    "lon": -1.202,
    "lat": 53.025,
    "property": ["clay"],
    "depth": ["0-5cm"],
    "value": "mean"
}
try:
    r = requests.get(rest_url, params=params, timeout=10)
    print(f"   Status: {r.status_code}")
    if r.status_code == 200:
        print("   ✓ REST API working!")
        print(f"   Response: {r.json()}")
    else:
        print(f"   ✗ Failed: {r.text[:200]}")
except Exception as e:
    print(f"   ✗ Error: {e}")

print()

# Test 3: BGS WMS GetCapabilities
print("3. Testing BGS WMS GetCapabilities...")
bgs_urls = [
    "https://map.bgs.ac.uk/arcgis/services/Soil/MapServer/WMSServer?request=GetCapabilities&service=WMS",
    "https://map.bgs.ac.uk/arcgis/services/UKSO/UKSO_BGS/MapServer/WMSServer?request=GetCapabilities&service=WMS",
]
for url in bgs_urls:
    try:
        r = requests.get(url, timeout=10)
        print(f"   URL: {url[:80]}...")
        print(f"   Status: {r.status_code}")
        if r.status_code == 200:
            print(f"   Content-Type: {r.headers.get('Content-Type', 'unknown')}")
            if 'xml' in r.headers.get('Content-Type', '').lower() or 'xml' in r.text[:200].lower():
                print("   ✓ BGS WMS GetCapabilities working!")
                # Try to find layer names
                if 'Layer' in r.text or 'layer' in r.text:
                    print("   Contains layer information")
                break
            else:
                print(f"   ✗ Unexpected response: {r.text[:200]}")
        else:
            print(f"   ✗ Failed: {r.text[:200]}")
    except Exception as e:
        print(f"   ✗ Error: {e}")

print()

# Test 4: Check for local UK soil data
print("4. Checking for local UK soil data files...")
project_root = Path(__file__).parent.parent
soil_dir = project_root / "soil" / "data"
print(f"   Looking in: {soil_dir}")
if soil_dir.exists():
    files = list(soil_dir.glob("*.tif")) + list(soil_dir.glob("*.gpkg"))
    if files:
        print(f"   ✓ Found {len(files)} data file(s):")
        for f in files:
            print(f"     - {f.name}")
    else:
        print("   ✗ No data files found")
        print("   To download UK soil data, see: soil/SOIL Datasets.md")
else:
    print("   ✗ Data directory does not exist")
    print("   To download UK soil data, see: soil/SOIL Datasets.md")

print("\n" + "="*60)
print("Summary: Check which APIs are working above")
print("="*60)

