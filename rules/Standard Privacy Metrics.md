## 🪟 **Standard Privacy Metrics** in Planning

Planning authorities in England often use **separation distances** as a proxy for privacy protection in residential contexts.

### 📏 A — **Privacy Separation Distances**

There’s a well-used design standard that:

* **≥ ~20-23 m** between *habitable room windows* facing each other to avoid overlooking. ([Designing Buildings][1])
* Some design codes use **≈21 m** minimum between principal windows across gardens. ([Planning Data][2])
* Where buildings face at oblique angles or in higher-density contexts, local authorities may accept reduced distances (e.g., ~15 m). ([Trafford Design Code][3])

These distances aren’t statutory — they’re **guidance used by LPAs’ planners and design codes** to judge acceptable privacy.

### 📏 B — **Special Considerations**

* **Upper-floor living rooms** may require larger spacing (e.g., ~35 m) due to overlooking potential. ([GOV.UK][4])
* Privacy can be enhanced with **screening, design orientation, internal layout** rather than just distance. ([GOV.UK][5])

**Practical takeaway:** The most applicable metric for privacy is a **distance between habitable room windows/facades** — typically 20-25 m in suburban contexts — but this is guided by local policy rather than fixed national law.


[1]: https://www.designingbuildings.co.uk/wiki/Privacy%20distance?utm_source=chatgpt.com "Privacy distance - Designing Buildings Wiki"
[2]: https://www.planning.data.gov.uk/entity/8700292?utm_source=chatgpt.com "Separation distances | Design code rule - Planning.data.gov.uk"
[3]: https://trafforddesigncode.uk/houses/plan-and-layout/?utm_source=chatgpt.com "Plan and Layout - Trafford Design Code"
[4]: https://assets.publishing.service.gov.uk/media/652e71cd6b6fbf000db757f9/EDG_Rear_Privacy.pdf?utm_source=chatgpt.com "[PDF] Rear Privacy - GOV.UK"
[5]: https://www.gov.uk/government/collections/planning-practice-guidance?utm_source=chatgpt.com "Planning practice guidance - GOV.UK"


## 📏 Simple **Metric: Distance from Nearest Building**

This is the one you *can code directly*.

### 🎯 What It Measures

A useful planning **privacy/setback metric** is:

> **d = min(distance(footprint A to any other building footprint))**

This gives a scalar distance from a plot/building to its nearest neighbour structure. Larger distance generally implies:

* more privacy
* less overshadowing
* less acoustic impact
* more potential for permitted development without consent

### 🧠 Implementation Idea (Python / Geo)

Here’s a high-level algorithm using something like Shapely/GeoPandas:

#### 1. **Load building footprints**

```python
import geopandas as gpd
from shapely.ops import nearest_points

buildings = gpd.read_file("buildings.geojson")  # footprints
```

#### 2. **Find nearest building for each**

```python
# spatial index
sindex = buildings.sindex

def nearest_distance(bldg):
    possible = buildings.iloc[list(sindex.intersection(bldg.bounds))]
    possible = possible[possible.geometry != bldg.geometry]
    if possible.empty:
        return float('inf')
    nearest_geom = nearest_points(bldg.geometry, possible.unary_union)[1]
    return bldg.geometry.distance(nearest_geom)

buildings["nearest_dist"] = buildings.apply(nearest_distance, axis=1)
```

This produces a **distance to nearest building** for each footprint.

### 📏 Extending Metric

You could also compute:

* **privacy_score = f(nearest_dist, orientation, window locations)**
* **sky exposure index** (how open the surrounding space is)
* **context factor** (whether in AONB / Conservation area)

But the **basic `nearest_dist`** gives a useful continuous metric.