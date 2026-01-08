#!/usr/bin/env python3
"""Merge multiple PVGIS JSON files into one"""

import json
from pathlib import Path

# Define the files with their azimuth angles and directions
files = [
    ('PVdata_51.013_-0.105_SA3_crystSi_1kWp_14_90deg_-90deg.json', -90, 'East'),
    ('PVdata_51.013_-0.105_SA3_crystSi_1kWp_14_90deg_0deg.json', 0, 'South'),
    ('PVdata_51.013_-0.105_SA3_crystSi_1kWp_14_90deg_90deg.json', 90, 'West'),
    ('PVdata_51.013_-0.105_SA3_crystSi_1kWp_14_90deg_-179deg.json', -179, 'North')
]

# Initialize merged data structure
merged_data = {
    'inputs': None,
    'outputs': {
        'by_azimuth': {}
    },
    'meta': None
}

# Process each file
for filename, azimuth, direction in files:
    filepath = Path(__file__).parent / filename
    
    if not filepath.exists():
        print(f"Warning: {filename} not found, skipping...")
        continue
    
    with open(filepath, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Use inputs and meta from first file (they should be the same)
    if merged_data['inputs'] is None:
        merged_data['inputs'] = data['inputs']
        merged_data['meta'] = data['meta']
    
    # Store outputs by azimuth
    az_key = str(azimuth)
    merged_data['outputs']['by_azimuth'][az_key] = {
        'azimuth': azimuth,
        'direction': direction,
        'monthly': data['outputs']['monthly'],
        'totals': data['outputs']['totals']
    }
    
    print(f"[OK] Processed {filename} (azimuth: {azimuth} deg, direction: {direction})")

# Write merged file
output_file = Path(__file__).parent / 'PVdata_51.013_-0.105_SA3_crystSi_1kWp_14_90deg_merged.json'
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(merged_data, f, indent=2, ensure_ascii=False)

print(f"\n[OK] Merged data written to: {output_file.name}")
print(f"  Total orientations: {len(merged_data['outputs']['by_azimuth'])}")

