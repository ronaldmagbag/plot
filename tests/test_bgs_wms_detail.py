"""
Test BGS WMS GetFeatureInfo to see actual response format
"""
import requests
import json

# Test coordinates (UK location)
lat = 53.025325345325484
lon = -1.2020914785758279

bgs_url = "https://map.bgs.ac.uk/arcgis/services/UKSO/UKSO_BGS/MapServer/WMSServer"

print("Testing BGS WMS GetFeatureInfo...\n")

# First, get capabilities to see available layers
print("1. Getting Capabilities...")
caps_params = {
    "SERVICE": "WMS",
    "VERSION": "1.3.0",
    "REQUEST": "GetCapabilities"
}
try:
    r = requests.get(bgs_url, params=caps_params, timeout=10)
    if r.status_code == 200:
        print(f"   ✓ Got capabilities (length: {len(r.content)} bytes)")
        # Try to find layer names
        content = r.text
        if '<Layer' in content or '<layer' in content:
            print("   Contains layer information")
            # Extract layer names (simple extraction)
            import re
            layers = re.findall(r'<Name>([^<]+)</Name>', content)
            if layers:
                print(f"   Found {len(layers)} layers:")
                for i, layer in enumerate(layers[:10]):  # Show first 10
                    print(f"     {i}: {layer}")
    else:
        print(f"   ✗ Failed: {r.status_code}")
except Exception as e:
    print(f"   ✗ Error: {e}")

print()

# Try GetFeatureInfo with different formats
print("2. Testing GetFeatureInfo...")

# Try with layer "0" (first layer)
for layer_name in ["0", "1", "2"]:
    print(f"\n   Trying layer: {layer_name}")
    
    # Try JSON format
    params_json = {
        "SERVICE": "WMS",
        "VERSION": "1.3.0",
        "REQUEST": "GetFeatureInfo",
        "LAYERS": layer_name,
        "QUERY_LAYERS": layer_name,
        "CRS": "EPSG:4326",
        "BBOX": f"{lon-0.0001},{lat-0.0001},{lon+0.0001},{lat+0.0001}",
        "I": "0",
        "J": "0",
        "WIDTH": "1",
        "HEIGHT": "1",
        "INFO_FORMAT": "application/json"
    }
    
    try:
        r = requests.get(bgs_url, params=params_json, timeout=10)
        print(f"     JSON Status: {r.status_code}")
        if r.status_code == 200:
            try:
                data = json.loads(r.text)
                print(f"     ✓ JSON Response: {json.dumps(data, indent=6)[:500]}")
            except:
                print(f"     Response (not JSON): {r.text[:300]}")
        else:
            print(f"     Response: {r.text[:200]}")
    except Exception as e:
        print(f"     ✗ Error: {e}")
    
    # Try text/plain format
    params_text = params_json.copy()
    params_text["INFO_FORMAT"] = "text/plain"
    
    try:
        r = requests.get(bgs_url, params=params_text, timeout=10)
        print(f"     Text Status: {r.status_code}")
        if r.status_code == 200:
            print(f"     ✓ Text Response: {r.text[:300]}")
    except Exception as e:
        print(f"     ✗ Error: {e}")

print("\n" + "="*60)
print("Use this information to improve BGS WMS parsing")
print("="*60)

