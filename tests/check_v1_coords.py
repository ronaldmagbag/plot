import json

# Load the JSON file
with open('output/v1.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

# Get centroid and property line
centroid = data['centroid']['coordinates']
prop_coords = data['boundaries']['property_line']['coordinates'][0]

print(f"Centroid: {centroid}")
print(f"Property line first point: {prop_coords[0]}")
print(f"Property line last point: {prop_coords[-1]}")
print(f"Property line has {len(prop_coords)} points")

# Calculate distance from centroid
lon_diff = abs(prop_coords[0][0] - centroid[0])
lat_diff = abs(prop_coords[0][1] - centroid[1])
print(f"\nDistance from centroid (first point):")
print(f"  lon_diff: {lon_diff:.6f} degrees")
print(f"  lat_diff: {lat_diff:.6f} degrees")

# Convert to meters
m_per_deg_lat = 111000
m_per_deg_lon = 111000 * 0.6  # Approximate for UK latitude
lon_m = lon_diff * m_per_deg_lon
lat_m = lat_diff * m_per_deg_lat
print(f"  In meters: lon={lon_m:.1f}m, lat={lat_m:.1f}m")
print(f"  Total distance: {(lon_m**2 + lat_m**2)**0.5:.1f}m")

# Check all property points
prop_lons = [c[0] for c in prop_coords]
prop_lats = [c[1] for c in prop_coords]
min_lon, max_lon = min(prop_lons), max(prop_lons)
min_lat, max_lat = min(prop_lats), max(prop_lats)

print(f"\nProperty line bounds:")
print(f"  lon: {min_lon:.6f} to {max_lon:.6f} (span: {max_lon - min_lon:.6f})")
print(f"  lat: {min_lat:.6f} to {max_lat:.6f} (span: {max_lat - min_lat:.6f})")

# Calculate extent in meters
lon_span_m = (max_lon - min_lon) * m_per_deg_lon
lat_span_m = (max_lat - min_lat) * m_per_deg_lat
print(f"  Span in meters: {lon_span_m:.1f}m x {lat_span_m:.1f}m")
print(f"  Max extent: {max(lon_span_m, lat_span_m):.1f}m")

# Check if centroid is within property bounds
centroid_in_bounds = (min_lon <= centroid[0] <= max_lon and min_lat <= centroid[1] <= max_lat)
print(f"\nCentroid within property bounds: {centroid_in_bounds}")

# Check surrounding_context
surrounding = data.get('surrounding_context', {})
print(f"\nSurrounding context:")
print(f"  roads: {len(surrounding.get('roads', []))} items")
print(f"  buildings: {len(surrounding.get('buildings', []))} items")

