

# **Plot Data Structure for Site Analysis**

To complement the floor plan metadata, a **plot (site) JSON** should capture key information about the land parcel and its context. Below is a structured plan (A–L) detailing what to include for a given plot, along with how each data point is obtained (external data, deterministic computation, LLM-generated, etc.). This information will enable evaluating how a floor plan fits on the site and assessing placement quality (sunlight, privacy, etc.). We then provide an example JSON snippet for a UK South plot.

## **A. Identification & Basic Geometry**

* **Plot ID:** Unique identifier for the site (e.g. `plot_001`). *Source: External/dataset input.*

* **Plot Type / Category:** Classification of the site context (e.g. *urban infill lot*, *suburban corner lot*, *rural acreage*). *Source: External or user-provided, possibly inferred from location density.*

* **Location:** Geographical data such as address or coordinates (latitude/longitude). *Source: External GIS data.*

* **Boundary Polygon:** The plot’s boundary as a polygon (list of coordinates). *Source: External survey or GIS (GeoJSON format).*

* **Area & Perimeter:** Total land area (m²) and perimeter length (m). *Source: Deterministically computed from the polygon.*

* **Dimensions & Aspect:** Characteristic length/width of the plot and aspect ratio (e.g. 30m × 20m, aspect 1.5). *Source: Computed from boundary extents.*

* **Orientation:** The bearing of the plot’s predominant axis or frontage relative to north (e.g. front side faces north). *Source: Computed from geometry or external map data.*

* **Shape Compactness:** A measure of shape regularity (e.g. plot compactness or irregularity). *Source: Computed (area ÷ area of bounding box).*

* **Orthogonality:** Percentage of plot edges aligned to cardinal directions (indicates if plot is rectangular vs. skewed). *Source: Computed from polygon angles.*

## **B. Plot Usage & Current Features**

* **Current Use:** Description of what currently exists on the plot (e.g. *vacant land*, *existing old house*, *parking lot*). *Source: External (site survey) or user input.*

* **Existing Structures:** Any structures on-site (footprint polygon, height, usage if any). *Source: External survey data.*

* **Natural Features:** Notable on-site features like large **trees**, water bodies, or terrain elements (with location on plot). *Source: External data or survey.*

* **Topography:** Ground elevation data (flat or slope). Include slope magnitude (%) and direction (e.g. slopes down toward south). *Source: External (topographic survey) and computed.*

* **Soil/Drainage:** (If relevant) Soil type or drainage conditions (for buildability). *Source: External data.*

## **C. Surrounding Context & Neighbors**

* **Adjacent Neighbors:** List of neighboring buildings/plots with their properties:

  * *Neighbor ID/Name:* Identifier or description (e.g. *house to east*).

  * *Position:* Direction relative to plot (e.g. east, west) or distance and bearing.

  * *Distance to Boundary:* How far the neighbor structure is from the plot’s edge (in meters).

  * *Building Height:* Height of neighbor building (m) or number of storeys[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=Here%E2%80%99s%20where%20you%20look%20beyond,noise%20levels%20in%20the%20neighborhood). *(External data, e.g. GIS or manual input.)*

  * *Building Type/Use:* What the neighbor is (e.g. 2-story house, apartment, school)[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=Here%E2%80%99s%20where%20you%20look%20beyond,noise%20levels%20in%20the%20neighborhood). *(External or inferred.)*

* **Sight Lines & Overlook:** Notation of sight lines from neighboring properties into the plot (which sides are overlooked). *Source: Deterministic inference from neighbor positions/heights (and assumed window orientations).*

* **Open Sides:** Which sides of the plot have no adjacent buildings (e.g. street or open land) and thus more openness. *Source: External (GIS) and computed.*

* **Street Frontage:** Which side(s) of the plot border a road or public street (e.g. *north side on Main St.*). *Source: External mapping.*

* **Neighborhood Character:** Qualitative context (e.g. “residential neighborhood with similar two-story houses”). *Source: LLM-generated summary from surrounding data.*

* **Noise Sources:** Major noise or pollution sources around (e.g. busy road on north, railway on west). *Source: External (traffic data or zoning maps) or user input.*

## **D. Access & Circulation**

* **Primary Access Point:** The main entry to the site (from which side/street). *Source: External (road mapping) or user input.*

* **Vehicular Access:** Possible driveway location or parking access on the plot. *Source: External inference from road layout.*

* **Pedestrian Access:** Sidewalks or footpaths connecting to the site. *Source: External data.*

* **Internal Circulation Space:** If applicable, area reserved for driveways or paths within the plot. *Source: Deterministic (e.g. if existing driveway given or to be allocated).*

* **Adjacent Street Width & Traffic:** Width of the fronting road and its traffic level (for entry and noise considerations). *Source: External (GIS or city data).*

* **Connectivity:** Proximity to intersections or multiple access points (e.g. corner lot with two street access). *Source: External.*

## **E. Zoning & Regulatory Constraints**

* **Zoning Designation:** Land use zoning (e.g. residential R1). *Source: External (municipal data).*

* **Allowed Building Footprint/Area:** Regulations like maximum buildable area or floor-area-ratio. *Source: External zoning code.*

* **Setback Requirements:** Required setbacks on each side (e.g. 5m front, 1m sides, 5m rear)[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=3). *(External/regulatory data.)*

* **Height Restrictions:** Maximum building height or storeys allowed on this plot[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=3). *(External/regulatory data.)*

* **Easements & Rights-of-Way:** Any portions of the plot that cannot be built on (utility lines, access easements). *Source: External (title or survey).*

* **Historical/Legal Constraints:** If in a conservation area, covenants, or future planned developments affecting the site[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=6). *Source: External.*

## **F. Utilities & Services**

* **Utility Connections:** Availability and location of utilities – water, sewer, electricity, gas, telecom. Note distance or entry points for each[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=8). *Source: External (site plans or utility maps).*

* **Drainage:** Existing drainage patterns or required on-site stormwater management (e.g. low spot on west corner). *Source: External survey.*

* **Infrastructure Constraints:** e.g. overhead power lines, manholes on site, which might affect building placement. *Source: External.*

## **G. Climate & Sun Path Data**

* **Climate Zone:** General climate information (e.g. temperate oceanic for UK South). *Source: External (based on location).*

* **Latitude (for Sun Path):** Latitude of the site (needed for solar calculations). *Source: External (from coordinates).*

* **Sun Path Characteristics:** Typical sun angles by season (e.g. low winter sun, high summer sun) relevant to design. *Source: Deterministic (compute from latitude) or external climate data.*

* **Prevailing Winds:** Main wind directions and intensity (e.g. frequent SW winds). *Source: External (weather data).*

* **Average Precipitation/Temperature:** Basic climate data for context (rainy climate, moderate temperatures). *Source: External (climate databases).*

* **Orientation for Sun:** Notation of which side is south (for northern hemisphere) and thus gets most sun. *Source: Computed from site orientation.*

* **Seasonal Sun Exposure:** Qualitative description of how sunlight reaches the site across seasons (e.g. “full southern sun year-round, eastern morning sun, western evening sun”). *Source: LLM-generated summary from latitude and neighbor info.*

## **H. Daylight & Solar Exposure Analysis**

* **Unobstructed Sky Areas:** Portions of the site open to direct sunlight (e.g. front yard, back yard) with minimal obstruction. *Source: Deterministic (analyze neighbor/building distances and heights).*

* **Neighbor Shadow Angles:** For each significant neighbor, the vertical angle up to the top of that building from the plot’s ground (indicates how much sun it blocks). \*Source: Computed (geometry). \*

* **Sun Exposure by Orientation:** Estimate of sunlight availability on each side of the plot (High/Med/Low) considering orientation and obstructions. *Source: Computed or LLM (using neighbor data and sun path)*.

* **Potential Solar Gain Zones:** Areas on the plot best suited for solar collection (e.g. the south-west corner has greatest afternoon sun). *Source: Computed (if doing solar simulation) or LLM inference.*

* **Daylight Hours Impact:** If available, data on how many hours of direct sun certain areas get (particularly in winter vs summer). *Source: External simulation or computed via sun path and obstacles.*

* **Outdoor Shadow Pattern:** Description of how shadows traverse the plot (e.g. “east neighbor casts morning shadow on site, but afternoons are clear”). *Source: LLM-generated from obstacle analysis.*

## **I. Privacy & Views Considerations**

* **Neighbor Overlook Risk:** Assessment of privacy from each side: e.g. high/medium/low risk of being overlooked by neighboring windows. *Source: Deterministic/LLM based on neighbor distance, height, likely windows.*

* **Sensitive Boundaries:** Identification of sides where privacy is most concerning (e.g. *east side very close to neighbor – sensitive for bedrooms*). *Source: Deterministic analysis.*

* **View Opportunities:** Notable outward views *from* the site (e.g. open view to park or landscape on a given side). *Source: External (GIS or site visit data) and LLM interpretation.*

* **Scenic Orientation:** Which direction offers the best vista (if any) to favor for main rooms/gardens. *Source: External/LLM (based on context features).*

* **Public Exposure:** Degree to which the site/house will be visible from public areas (street side exposure). *Source: Deterministic (frontage length, street width).*

* **Noise Exposure:** Which areas of the plot are exposed to noise or foot traffic (e.g. front facing busy road vs. quiet rear garden)[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=Here%E2%80%99s%20where%20you%20look%20beyond,noise%20levels%20in%20the%20neighborhood). *Source: External (traffic data) or LLM summary.*

## **J. Deterministic Flags & Warnings**

*(Boolean flags or alerts auto-generated from the above data to highlight design-relevant conditions)*

* **Overshadowing Risk Flag:** `true` if a neighboring structure is likely to significantly shadow the plot (e.g. tall building directly south).

* **Privacy Risk Flag:** `true` if the site is heavily overlooked (multiple close neighbors with direct sight lines).

* **Noise Risk Flag:** `true` if adjacent to high noise source (busy road, etc.).

* **Slope Flag:** `true` if the terrain slope is steep enough to impact building design.

* **Access Constraint Flag:** `true` if site access is limited (e.g. landlocked or very narrow entrance).

* **Irregular Shape Flag:** `true` if plot shape is highly irregular, which may constrain building footprint placement.

* **Small Plot/Setback Flag:** `true` if buildable area after setbacks is very limited.

* **Code Compliance Flags:** e.g. flag if *potential issue* like building footprint would exceed coverage or if a known restriction (tree preservation etc.) exists.  
   *(These flags are deterministically set based on threshold rules applied to the raw data.)*

## **K. LLM-Assisted Semantic Metadata**

*(Derived insights in natural language, generated by AI based on all the above factual data.)*

* **Site Archetype Label:** A concise label characterizing the site’s overall scenario – *“Narrow urban infill with close neighbors”*, *“Suburban lot with front street and rear garden”*, etc. *Source: LLM.*

* **Site Summary Description:** A paragraph describing the plot’s context, features, and important constraints/opportunities in plain language. *Source: LLM (using factual data).*

* **Primary Strengths:** Bullet points summarizing the site’s advantages (each tied to evidence/facts). For example: *“Excellent southern exposure for sunlight”*, *“Rear side opens to green space providing privacy and views”*. *Source: LLM (fact-based).*

* **Primary Constraints:** Bullet points for site challenges, e.g. *“Eastern neighbor is very close, posing privacy concerns”*, *“Limited frontage width restricts design footprint”*. *Source: LLM.*

* **Sunlight Potential Rating:** Qualitative rating (High/Medium/Low) of the site’s overall daylight potential for a new home, with justification (e.g. High because of unobstructed south side). *Source: LLM (based on sunlight analysis).*

* **Privacy Suitability Rating:** Rating of how well the site can provide privacy (e.g. *Moderate* – good privacy at rear, but sides are exposed). *Source: LLM.*

* **Optimal Orientation Advice:** Recommendation for how to orient/place the house on the plot for best results (e.g. *“Align living areas toward the south garden; use minimal side windows on east side”*). *Source: LLM (combining floor plan site-fit preferences with site data).*

* **Best Use Scenario:** Suggested occupant or design scenario the site is best suited for (e.g. *“Ideal for a family home that prioritizes a sunny backyard garden”*). *Source: LLM.*

* **Confidence Scores:** Confidence level for each LLM-generated claim or rating, based on data completeness. *Source: LLM (self-assessed).*

* **Explicit Unknowns:** Listing of any important unknown factors (e.g. *“Exact positions of neighbor windows not known”*, *“No data on soil stability”*). *Source: LLM.*

## **L. Retrieval & Matching Aids**

* **Feature Vector:** A numerical vector embedding summarizing key site characteristics for similarity search (combining metrics like area, openness, neighbor density, etc.). *Source: Deterministic or ML model.*

* **Tags:** A set of keywords describing the site for easy lookup (e.g. `["corner_lot", "south_facing_garden", "dense_neighbors", "flat_site"]`). *Source: Deterministic/LLM based on data.*

* **Comparable Site References:** (If applicable) Reference IDs of similar sites in the dataset (for matching floor plans or case studies). *Source: Deterministic (vector similarity).*

* **Potential Floor Plan Match Hints:** Attributes to match with floor plans, e.g. *“Fits plans with front entry on north side”*, *“Needs designs with limited side windows”*. *Source: LLM or rule-based, derived from site constraints and floor plan metadata.*

---

## **Example Plot JSON Structure (UK South)**

Below is an **example** JSON snippet for a hypothetical plot in Southern UK, illustrating how the above information could be structured:

`{`  
  `"plot_id": "UKS-001",`  
  `"type": "suburban_residential",`  
  `"location": {`  
    `"latitude": 51.0,`  
    `"longitude": -0.1,`  
    `"region": "UK-South",`  
    `"address": "1 Example Road, Townshire, UK"`  
  `},`  
  `"geometry": {`  
    `"boundary_polygon": [`  
      `[0,0], [20,0], [20,30], [0,30]  /* coordinates in meters */`  
    `],`  
    `"area_sqm": 600.0,`  
    `"perimeter_m": 100.0,`  
    `"approx_dimensions": { "length_m": 30, "width_m": 20 },`  
    `"orientation_front": "north",`   
    `"aspect_ratio": 1.5,`  
    `"compactness": 0.98,`  
    `"orthogonality": 1.0`  
  `},`  
  `"current_use": "vacant_land",`  
  `"topography": {`  
    `"slope_percent": 2.0,`  
    `"slope_direction": "south",`  
    `"elevation_avg_m": 50`  
  `},`  
  `"context": {`  
    `"neighbors": [`  
      `{`  
        `"id": "N-EAST",`  
        `"direction": "east",`  
        `"distance_to_plot_m": 5,`  
        `"height_m": 9,`  
        `"stories": 2,`  
        `"type": "detached_house",`  
        `"usage": "residential"`  
      `},`  
      `{`  
        `"id": "N-WEST",`  
        `"direction": "west",`  
        `"distance_to_plot_m": 5,`  
        `"height_m": 9,`  
        `"stories": 2,`  
        `"type": "detached_house",`  
        `"usage": "residential"`  
      `},`  
      `{`  
        `"id": "N-SOUTH",`  
        `"direction": "south",`  
        `"distance_to_plot_m": 25,`  
        `"height_m": 8,`  
        `"stories": 2,`  
        `"type": "detached_house",`  
        `"usage": "residential"`  
      `}`  
    `],`  
    `"open_sides": ["south"],`   
    `"adjacent_road": { "side": "north", "name": "Example Road", "type": "local_residential" },`  
    `"noise_levels": { "north": "medium (road traffic)", "south": "low", "east": "low", "west": "low" }`  
  `},`  
  `"access": {`  
    `"primary_access_side": "north",`  
    `"vehicle_access": true,`  
    `"pedestrian_access": true,`  
    `"corner_plot": false`  
  `},`  
  `"regulations": {`  
    `"zoning": "R1-single_family",`  
    `"setbacks_m": { "front": 5, "side": 1, "rear": 5 },`  
    `"max_height_m": 10,`  
    `"max_coverage_ratio": 0.5,`  
    `"easements": []`  
  `},`  
  `"utilities": {`  
    `"water": "available_at_street",`  
    `"sewer": "available_at_street",`  
    `"electricity": "on_site",`  
    `"gas": "available_at_street",`  
    `"telcom": "on_site"`  
  `},`  
  `"climate": {`  
    `"zone": "temperate_oceanic",`  
    `"avg_annual_temp_C": 10.5,`  
    `"avg_annual_rain_mm": 800,`  
    `"prevailing_wind": "SW",`  
    `"sun_path": "low winter sun (15° alt), high summer sun (~62° alt)"`  
  `},`  
  `"sunlight_analysis": {`  
    `"sun_exposure": { "north": "Low", "south": "High", "east": "Medium", "west": "Medium" },`  
    `"neighbor_shadow_angles": {`  
      `"east_deg": 60,   /* neighbor to east blocks sun below ~60° elevation */`  
      `"west_deg": 60,`  
      `"south_deg": 18   /* neighbor to south is farther, lower angle block */`  
    `},`  
    `"winter_sun_hours_est": 4,    /* approx direct sun hours in midwinter */`  
    `"summer_sun_hours_est": 10   /* approx direct sun hours in midsummer */`  
  `},`  
  `"privacy_analysis": {`  
    `"overlook_risk": { "north": "Medium", "south": "Medium", "east": "High", "west": "High" },`  
    `"remarks": "Close side neighbors can view into side windows and garden; rear has some distance for privacy."`  
  `},`  
  `"flags": {`  
    `"overshadow_risk": false,`  
    `"privacy_risk": true,`  
    `"noise_risk": false,`  
    `"slope_risk": false,`  
    `"access_constraint": false,`  
    `"irregular_shape": false`  
  `},`  
  `"LLM_semantics": {`  
    `"site_archetype": "Suburban lot with front street and flanking neighbors",`  
    `"summary": "Rectangular 600 m² suburban plot, fronting a quiet road to the north with adjacent houses on the east and west. The site is flat with an open south side leading to back gardens, offering excellent sunlight from the south. Neighbors on both sides are two-story houses very close to the boundaries, which raises privacy concerns for side-facing windows. The rear (south) side is more private and sunny, ideal for a garden or main living space orientation.",`  
    `"primary_strengths": [`  
      `"Unobstructed southern exposure for ample daylight on the rear side",`  
      `"Flat, rectangular lot shape is easy to build on and orient",`  
      `"Access from north side street keeps the sunny rear garden private"`  
    `],`  
    `"primary_constraints": [`  
      `"Eastern and western neighbors are in close proximity (about 5m), causing potential overlooking into side rooms or yard",`  
      `"Limited plot width (20m) and side setbacks might restrict building width and placement of windows on sides"`  
    `],`  
    `"sunlight_potential_rating": "High",`  
    `"privacy_suitability": "Moderate",`  
    `"optimal_orientation": "Place the house towards the north side, with main living areas and windows facing south into the garden. Minimize large side-facing windows or use frosted glass to mitigate direct views from east/west neighbors.",`  
    `"best_use_scenario": "Ideal for a family home design that prioritizes a south-facing backyard for outdoor living, given the excellent sunlight and the need to maximize privacy away from side neighbors.",`  
    `"confidence": {`  
      `"sunlight_potential_rating": 0.9,`  
      `"privacy_suitability": 0.7,`  
      `"optimal_orientation": 0.85`  
    `},`  
    `"unknowns_noted": [`  
      `"Exact positions of neighbor windows not known, assumptions made on typical layouts",`  
      `"No data on any large trees that might cast additional shadows"`  
    `]`  
  `},`  
  `"matching_aids": {`  
    `"feature_vector": [600, 20, 30, 2, 1, 0.5, 1, ...],`  
    `"tags": ["suburban", "flat", "south_open", "close_neighbors", "north_frontage"]`  
  `}`  
`}`

This structured data covers: basic site identity and geometry, context (neighbors, roads) and environmental factors, computed metrics for sun and privacy, regulatory limits, and high-level insights. Such a **plot JSON**, used in conjunction with a floor plan JSON, allows algorithms or an LLM to evaluate floor plan placement options on the site and to score or explain the fit across dimensions like sunlight access and privacy. The external factual fields (location, neighbors, climate, etc.) feed into deterministic computations (areas, shadow angles, flags) and LLM-generated summaries, ensuring a comprehensive understanding of the plot[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=,Views%2C%20surrounding%20architecture%20and%20landscaping)[projectmanager.com](https://www.projectmanager.com/blog/site-analysis-in-architecture#:~:text=Here%E2%80%99s%20where%20you%20look%20beyond,noise%20levels%20in%20the%20neighborhood).

