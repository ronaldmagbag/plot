

Got it. Given **only** the external plot inputs you listed, here’s what you can still generate **deterministically** (compute it), what you can generate **LLM-assisted** (summaries, “site archetype”, recommendations) and a concrete **Example Plot JSON** that fits your pipeline needs (plan selection \+ placement/orientation \+ scoring).

I’ll also flag a few **additional external data points that are usually easy to acquire** in the UK South and are very high leverage.

---

## **A few external data points you didn’t list that are often easy (UK) and high value**

If you truly won’t have them in 5 weeks, ignore — but these are commonly “easy wins”:

1. **Flood risk zone** (very high leverage for feasibility \+ cost \+ “architect brain” credibility).

2. **Conservation area / listed building / AONB / Green Belt flags** (planning feasibility).

3. **Road classification \+ speed** (local road vs A-road, impacts access, noise, visibility splays sometimes).

4. **Existing building footprints (if any) on the plot** (even if you’re not reusing them).

5. **Rights of way / public footpaths** crossing or bordering the parcel (can be deal-killers).

If you can’t get these, the system can still be valuable, but these 5 make your “site analysis” feel much more like an architect.

---

## **What you can generate for the plot with your listed external inputs**

### **External inputs you said you will have (raw)**

* Location (lat/lon)

* Boundary polygon

* Topography (DEM / topo map)

* Soil conditions (\~250m resolution)

* Adjacent neighbors: position, distance-to-boundary, building height

* Street frontage

* Primary access point

* Adjacent street width

* Setback requirements

* Climate zone

* Latitude (sun path)

* Prevailing winds (maybe)

* Avg precipitation/temperature

* Orientation for sun (you can compute from lat/lon \+ north; but you say you’ll have it)

Below is what you can compute on top.

---

# **1\) Deterministic derived plot metadata (compute, don’t guess)**

## **1.1 Coordinate frames & normalization (critical for everything else)**

Even if your boundary starts in WGS84 (lat/lon), you want a **local metric coordinate system**.

**Generate deterministically**

* `crs.source`: `"EPSG:4326"` (or whatever you receive)

* `crs.local`: `"EPSG:27700"` (British National Grid) or a local tangent-plane ENU system

* `transform.wgs84_to_local`: stored parameters (or just re-project on the fly)

**Why it matters:** all scoring (fit, distances, setbacks, slopes) becomes stable in meters.

---

## **1.2 Plot geometry metrics (from boundary polygon)**

**Generate deterministically**

* `plot.area_m2`

* `plot.perimeter_m`

* `plot.bounding_box` (min/max x/y)

* `plot.dimensions_approx` (length/width estimates)

* `plot.aspect_ratio` (long/short)

* `plot.compactness` (e.g. `4πA / P²`)

* `plot.irregularity_score` (1 − rectangularity)

* `plot.orthogonality_score` (% edges near 0/90/180/270°)

* `plot.primary_axis_bearing_deg` (dominant direction of long axis)

* `plot.frontage_edges` (which boundary segments touch street frontage)

* `plot.front_direction_bearing_deg` (bearing from plot interior toward street)

**Useful for:** selecting plans with compatible footprint proportions; placing building without “clipping”; choosing which façade should become “front”.

---

## **1.3 Setbacks → buildable envelope**

Given numeric setbacks, you can produce a *real buildable polygon*.

**Generate deterministically**

* `setbacks.by_edge`: attach each boundary edge with its setback (front/side/rear or explicit)

* `buildable_envelope.polygon` (offset-in polygon after setbacks)

* `buildable_envelope.area_m2`

* `buildable_envelope.coverage_ratio = envelope_area / plot_area`

* `buildable_envelope.max_inscribed_rectangle` (approx; extremely useful for plan fit)

* `buildable_envelope.clearance_map` (optional: distance-to-boundary field sampled on a grid)

**Useful for:** “Does this footprint fit at all?” and for placement optimization.

---

## **1.4 Topography summaries (from DEM/topo)**

Even if you choose **not** to model terrain in 3D, you can still compute highly valuable site metrics.

**Generate deterministically**

* `topography.dem_source` \+ resolution \+ timestamp (if known)

* `topography.elevation_stats`: min/mean/max (m)

* `topography.relief_m = max - min` within plot

* `topography.slope_stats_percent`: min/mean/max (or degrees)

* `topography.slope_direction_aspect_deg` (dominant downhill direction)

* `topography.roughness_index` (e.g. std of slope or curvature proxy)

* `topography.flatness_score` (e.g. % of plot under 5% slope)

* `topography.slope_flag`: `flat / moderate / steep` based on thresholds you set

* `drainage.expected_flow_direction` (usually downhill aspect)

* `drainage.surface_runoff_risk` (very coarse proxy: high if steep \+ high rainfall)

* `earthworks.risk_level` (LOW/MED/HIGH based on slope \+ relief)

* `earthworks.estimated_cost_allowance_range` (if you’re comfortable; otherwise output “needs QS”)

**Important:** with your inputs, you *can* do a **good high-level earthworks risk \+ allowance** without explicitly modifying the DEM.

---

## **1.5 Soil summaries (coarse)**

At \~250m resolution, treat this as advisory.

**Generate deterministically**

* `soil.source` \+ resolution

* `soil.class` (whatever classification you get)

* `soil.confidence`: `low/medium` (almost always medium→low at 250m)

* `soil.foundation_risk_proxy`: `LOW/MED/HIGH` (based on soil class)

* `soil.drainage_infiltration_proxy`: `LOW/MED/HIGH`

* `soil.disclaimer`: “needs site investigation”

**Useful for:** “architect brain” credibility and early budget allowances.

---

## **1.6 Street frontage \+ access geometry**

Given street frontage \+ access point \+ street width, you can compute a lot.

**Generate deterministically**

* `access.frontage_edges` (boundary edges that touch street)

* `access.primary_access_point_local_xy`

* `access.street_width_m`

* `access.street_axis_bearing_deg` (direction of the road segment adjacent)

* `access.front_setback_edge_id` (choose which edge is “front”)

* `access.feasible_driveway_corridor` (a polygon “cone” from access point into plot)

* `access.entry_preference_zones` (areas in buildable envelope that keep entry near access)

**Optional but valuable:**

* `access.driveway_length_est_range` (distance from access point to likely parking/garage zone)

* `access.parking_feasibility_score` (very rough; depends on buildable envelope width/depth and street width)

Even without perfect driveway design, you can produce a compelling “we understand access constraints” report.

---

## **1.7 Neighbor-derived privacy and daylight obstruction proxies**

You said you’ll have neighbor **position**, **distance-to-boundary**, and **height** (not full footprints). That’s enough to generate strong *risk maps* and summary scores even if you can’t do exact window-to-window lines.

**Generate deterministically**  
 For each neighbor “blocker”:

* `neighbors[i].bearing_deg` (direction from plot centroid or from relevant boundary)

* `neighbors[i].distance_to_boundary_m`

* `neighbors[i].height_m`

* `neighbors[i].vertical_obstruction_angle_deg = atan(height / distance)`  
   (huge—works for both shadow and privacy pressure)

* `neighbors[i].privacy_pressure_score` (a monotonic function of height/distance)

* `neighbors[i].daylight_block_risk` (primarily relevant if neighbor is south-ish)

Aggregate summaries:

* `privacy.boundary_risk_by_sector` (e.g. N/NE/E/… each with a score)

* `privacy.overall_risk_score`

* `privacy.sensitive_edges` (edges where privacy pressure is highest)

* `daylight.obstruction_by_sector` (how blocked each direction is)

* `daylight.winter_sun_risk_score` (especially if tall neighbors are to the south)

This lets you later score plan placements by:

* “Which façade should be ‘controlled’ (privacy)?”

* “Which façade should get more glazing (sun/view)?”

Even with limited neighbor data, you can still do something meaningful and defensible.

---

## **1.8 Sun path and climate-derived design parameters**

With latitude \+ climate zone (+ optional winds), you can produce architect-relevant “rules of thumb” that later guide placement and façade decisions.

**Generate deterministically**

* `solar.key_dates` (winter solstice / equinox / summer solstice)

* `solar.noon_altitude_deg_by_date`

* `solar.azimuth_range_by_date` (sunrise/sunset azimuth)

* `solar.south_bearing_deg` (local “true south” in plot coordinate frame)

* `climate.heating_dominated_proxy` (UK South is typically heating-dominated; but compute from climate zone)

* `climate.rainfall_band` (low/med/high, from precipitation)

* `wind.prevailing_dir_deg` (if you have it)

* `wind.exposure_proxy` (combine wind \+ openness, very rough)

---

# **2\) LLM-assisted plot metadata (grounded, evidence-backed)**

Given you’ll have the deterministic plot facts above, the LLM layer should NOT invent new facts. It should only:

* **Summarize**

* **Tag**

* **Recommend** *with explicit evidence references*

* **List unknowns**

### **2.1 Good LLM outputs for a plot**

**Generate via LLM (with citations back to computed fields)**

* `site_archetype_label`  
   Examples:

  * “Rectangular suburban lot with close side neighbors”

  * “Corner lot with strong street exposure”

* `constraints_summary` (bullet list)

* `opportunities_summary` (bullet list)

* `design_priorities_recommendation`  
   e.g. “Prioritize south-facing garden glazing; keep east façade controlled due to privacy pressure.”

* `data_confidence_and_unknowns`  
   e.g. “Soil is coarse (250m). Neighbor windows unknown.”

* `questions_to_ask_user` (high leverage for your chatbot)  
   e.g. “Do you care more about afternoon sun or morning sun?”  
   “Is privacy from the east neighbor critical?”

### **2.2 Plot “tags” for retrieval**

**Generate via LLM (but deterministic tags are better when possible)**

* `tags`: `["flat_site", "side_neighbors_close", "south_opportunity", "narrow_frontage", ...]`

---

# **3\) Example Plot JSON structure (aligned to your placement engine)**

Below is a concrete structure that matches your needs: keep **raw external**, **derived deterministic**, and **LLM semantics** separated.

`{`  
  `"plot_id": "UKS-EX-0001",`  
  `"region": "UK-South",`  
  `"external_raw": {`  
    `"location": {`  
      `"lat": 51.123456,`  
      `"lon": -0.123456,`  
      `"address_string": "Example Rd, Sussex, UK"`  
    `},`  
    `"boundary": {`  
      `"crs": "EPSG:4326",`  
      `"polygon_wgs84": [`  
        `[-0.12350, 51.12340],`  
        `[-0.12310, 51.12340],`  
        `[-0.12310, 51.12370],`  
        `[-0.12350, 51.12370]`  
      `]`  
    `},`  
    `"topography": {`  
      `"dem_source": "UK_LiDAR_or_DEM",`  
      `"dem_resolution_m": 2.0,`  
      `"elevation_grid_ref": "S3://.../dem.tif"`  
    `},`  
    `"soil": {`  
      `"soil_source": "SoilGrid_or_UK_dataset",`  
      `"resolution_m": 250,`  
      `"soil_class_raw": "Loamy soil, moderate drainage"`  
    `},`  
    `"neighbors": [`  
      `{`  
        `"neighbor_id": "N-E",`  
        `"bearing_deg": 90,`  
        `"distance_to_boundary_m": 4.5,`  
        `"height_m": 8.5`  
      `},`  
      `{`  
        `"neighbor_id": "N-W",`  
        `"bearing_deg": 270,`  
        `"distance_to_boundary_m": 5.0,`  
        `"height_m": 9.0`  
      `},`  
      `{`  
        `"neighbor_id": "N-S",`  
        `"bearing_deg": 180,`  
        `"distance_to_boundary_m": 18.0,`  
        `"height_m": 7.5`  
      `}`  
    `],`  
    `"street": {`  
      `"frontage_edges_hint": ["edge_0"],`  
      `"primary_access_point_wgs84": [-0.12330, 51.12340],`  
      `"adjacent_street_width_m": 6.0`  
    `},`  
    `"regulations": {`  
      `"setbacks_m": {`  
        `"front": 5.0,`  
        `"side": 1.0,`  
        `"rear": 5.0`  
      `}`  
    `},`  
    `"climate": {`  
      `"climate_zone": "UK_Temperate_Oceanic",`  
      `"latitude_deg": 51.123456,`  
      `"avg_precip_mm_per_year": 850,`  
      `"avg_temp_c": 10.5,`  
      `"prevailing_wind_bearing_deg": 225`  
    `},`  
    `"sun": {`  
      `"north_bearing_deg": 0`  
    `}`  
  `},`

  `"derived": {`  
    `"crs": {`  
      `"local_crs": "EPSG:27700",`  
      `"notes": "Reproject boundary and access point for all metric computations."`  
    `},`

    `"geometry": {`  
      `"boundary_polygon_local_xy": [[0,0],[20,0],[20,30],[0,30]],`  
      `"area_m2": 600.0,`  
      `"perimeter_m": 100.0,`  
      `"bounding_box": { "min_x": 0, "min_y": 0, "max_x": 20, "max_y": 30 },`  
      `"approx_dimensions_m": { "width": 20, "depth": 30 },`  
      `"aspect_ratio": 1.5,`  
      `"compactness": 0.75,`  
      `"orthogonality_score": 0.98,`  
      `"primary_axis_bearing_deg": 0`  
    `},`

    `"setbacks_and_envelope": {`  
      `"setbacks_by_edge_m": [`  
        `{ "edge_id": "edge_0", "type": "front", "value_m": 5.0 },`  
        `{ "edge_id": "edge_1", "type": "side", "value_m": 1.0 },`  
        `{ "edge_id": "edge_2", "type": "rear", "value_m": 5.0 },`  
        `{ "edge_id": "edge_3", "type": "side", "value_m": 1.0 }`  
      `],`  
      `"buildable_envelope_polygon_local_xy": [[1,5],[19,5],[19,25],[1,25]],`  
      `"buildable_envelope_area_m2": 360.0,`  
      `"coverage_ratio_envelope_to_plot": 0.6,`  
      `"max_inscribed_rect_approx": { "width_m": 18.0, "depth_m": 20.0 }`  
    `},`

    `"topography": {`  
      `"elevation_stats_m": { "min": 49.2, "mean": 50.1, "max": 51.3 },`  
      `"relief_m": 2.1,`  
      `"slope_stats_percent": { "min": 0.5, "mean": 3.0, "max": 6.5 },`  
      `"dominant_downhill_aspect_deg": 180,`  
      `"flatness_score": 0.85,`  
      `"earthworks_risk_level": "LOW",`  
      `"drainage_flow_direction_deg": 180`  
    `},`

    `"soil": {`  
      `"soil_class_normalized": "loam_moderate_drainage",`  
      `"soil_confidence": "LOW",`  
      `"foundation_risk_proxy": "MEDIUM",`  
      `"infiltration_proxy": "MEDIUM",`  
      `"disclaimer": "Coarse (250m) soil data; requires site investigation for foundations."`  
    `},`

    `"access": {`  
      `"street_frontage_edges": ["edge_0"],`  
      `"primary_access_point_local_xy": [10.0, 0.0],`  
      `"street_width_m": 6.0,`  
      `"front_direction_bearing_deg": 0,`  
      `"driveway_corridor_polygon_local_xy": [[8,0],[12,0],[15,10],[5,10]]`  
    `},`

    `"neighbors": {`  
      `"by_neighbor": [`  
        `{`  
          `"neighbor_id": "N-E",`  
          `"bearing_deg": 90,`  
          `"distance_to_boundary_m": 4.5,`  
          `"height_m": 8.5,`  
          `"vertical_obstruction_angle_deg": 62.0,`  
          `"privacy_pressure_score": 0.82,`  
          `"daylight_block_risk": "LOW"`  
        `},`  
        `{`  
          `"neighbor_id": "N-S",`  
          `"bearing_deg": 180,`  
          `"distance_to_boundary_m": 18.0,`  
          `"height_m": 7.5,`  
          `"vertical_obstruction_angle_deg": 22.6,`  
          `"privacy_pressure_score": 0.35,`  
          `"daylight_block_risk": "MEDIUM"`  
        `}`  
      `],`  
      `"privacy": {`  
        `"risk_by_sector": {`  
          `"N": 0.4,`  
          `"E": 0.82,`  
          `"S": 0.35,`  
          `"W": 0.78`  
        `},`  
        `"sensitive_sectors": ["E", "W"],`  
        `"overall_privacy_risk_score": 0.68`  
      `},`  
      `"daylight_obstruction": {`  
        `"risk_by_sector": {`  
          `"S": 0.55,`  
          `"E": 0.25,`  
          `"W": 0.25,`  
          `"N": 0.05`  
        `},`  
        `"winter_sun_risk_score": 0.55`  
      `}`  
    `},`

    `"solar": {`  
      `"latitude_deg": 51.123456,`  
      `"key_dates": {`  
        `"winter_solstice": { "noon_altitude_deg": 15.0 },`  
        `"equinox": { "noon_altitude_deg": 39.0 },`  
        `"summer_solstice": { "noon_altitude_deg": 62.0 }`  
      `},`  
      `"true_south_bearing_deg": 180,`  
      `"prevailing_wind_bearing_deg": 225,`  
      `"climate_summary": {`  
        `"avg_temp_c": 10.5,`  
        `"avg_precip_mm_per_year": 850,`  
        `"rainfall_band": "MEDIUM"`  
      `}`  
    `},`

    `"flags": {`  
      `"small_envelope_flag": false,`  
      `"irregular_shape_flag": false,`  
      `"steep_slope_flag": false,`  
      `"privacy_risk_flag": true,`  
      `"winter_sun_obstruction_flag": "MEDIUM"`  
    `},`

    `"matching_aids": {`  
      `"numeric_feature_vector": {`  
        `"plot_area_m2": 600,`  
        `"envelope_area_m2": 360,`  
        `"aspect_ratio": 1.5,`  
        `"mean_slope_percent": 3.0,`  
        `"privacy_risk_score": 0.68,`  
        `"winter_sun_risk_score": 0.55,`  
        `"street_width_m": 6.0`  
      `},`  
      `"tags_deterministic": [`  
        `"rectangular_plot",`  
        `"flat_to_moderate_slope",`  
        `"close_side_neighbors",`  
        `"north_frontage"`  
      `]`  
    `}`  
  `},`

  `"llm_semantics": {`  
    `"site_archetype_label": "Rectangular suburban lot with close side neighbors",`  
    `"constraints": [`  
      `{`  
        `"claim": "Side privacy pressure is high; large side windows will likely feel exposed.",`  
        `"evidence_refs": ["derived.neighbors.privacy.risk_by_sector.E", "derived.neighbors.privacy.risk_by_sector.W"],`  
        `"confidence": 0.8`  
      `}`  
    `],`  
    `"opportunities": [`  
      `{`  
        `"claim": "Rear (south) side is the best candidate for primary garden glazing and outdoor living.",`  
        `"evidence_refs": ["derived.solar.true_south_bearing_deg", "derived.neighbors.daylight_obstruction.risk_by_sector.S"],`  
        `"confidence": 0.7`  
      `}`  
    `],`  
    `"unknowns": [`  
      `"Exact neighbor window locations and overlooking angles are unknown.",`  
      `"Soil data is coarse (250m) — site investigation required."`  
    `],`  
    `"questions_to_ask_user": [`  
      `"Do you prioritize privacy from side neighbors over maximizing side windows?",`  
      `"Do you want morning sun (east) or afternoon/evening sun (west) in living spaces?"`  
    `]`  
  `}`  
`}`

---

## **What’s the minimum “plot JSON” you need to enable best placement scoring?**

If you want to keep things very lean for 5 weeks, the **must-have derived fields** (beyond raw external) are:

1. `boundary_polygon_local_xy`

2. `buildable_envelope_polygon_local_xy` (+ envelope area)

3. `topography.slope_stats + downhill aspect`

4. `neighbors[i].bearing, distance, height → vertical_obstruction_angle`

5. `solar.key_dates + true_south_bearing`

6. `access.frontage_edges + access point local xy + street width`

7. `privacy_risk_by_sector` and `daylight_obstruction_by_sector` (aggregate summaries)

That set alone will let you:

* place footprints within the envelope,

* choose orientation based on sun,

* penalize orientations that put “big glazing” toward privacy risk sectors,

* flag slope/earthworks risk,

* produce an “architect-style” explanation.

