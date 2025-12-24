

**UK council systems *do* appear perfectly aligned**, even though the raw datasets are not. That’s because councils **don’t visualise data the way Mapbox + OSM does**.

Below is exactly **how councils handle (and hide) the mismatch**, step by step.

---

## 1️⃣ What councils treat as “truth”

UK local authorities do **not** treat satellite imagery or OSM as authoritative.

Their visual hierarchy is fixed:

### Legal truth (top priority)

* Ordnance Survey MasterMap
* INSPIRE cadastral polygons
* Surveyed boundaries (EPSG:27700)

These layers are:

* Survey-grade
* Internally consistent
* Never warped to fit imagery

---

## 2️⃣ Why council maps look “perfectly aligned”

### Key reason: **everything is rendered in EPSG:27700**

Councils **do not use Web Mercator**.

Typical council stack:

* CRS: **British National Grid (EPSG:27700)**
* Renderer: GeoServer / ArcGIS Server / QGIS Server
* Base map: OS MasterMap raster or vector

➡ No WGS84
➡ No EPSG:3857
➡ No satellite reprojection drift

This alone removes ~80% of the mismatch you’re seeing.

---

## 3️⃣ How imagery is handled (this is the trick)

When councils show aerial imagery:

### They **visually subordinate it**

* Lower opacity
* Slight blur
* Not used for snapping
* Often clipped under parcel boundaries

Some councils even add **disclaimers directly on the map**:

> “Imagery is indicative only and not suitable for measurement.”

This is intentional.

---

## 4️⃣ The “snap illusion” councils use

Councils **never snap parcels to imagery**.

Instead, they:

* Snap imagery *visually* under legal vectors
* Use thicker boundary strokes
* Use halo outlines

### Typical styling

* Parcel boundary: **2–3 px solid red**
* Building footprint: **bold fill**
* Imagery: muted / grayscale

Your eye aligns to the boundary — not the photo.

---

## 5️⃣ Why OSM + Mapbox can never match this

Your stack:

* OSM vectors → EPSG:4326
* Mapbox imagery → EPSG:3857
* INSPIRE → EPSG:27700 → 4326 → 3857

Council stack:

* Everything → **EPSG:27700**
* Imagery warped *once*, offline, to OS grid
* Vectors never transformed for display

➡ You are comparing **engineering-grade data** to **consumer-grade visual tiles**.

---

## 6️⃣ What councils explicitly forbid (important)

UK planning guidance explicitly states:

* ❌ Do not measure from imagery
* ❌ Do not digitise boundaries from photos
* ❌ Do not adjust legal boundaries visually

This is why:

* Planning portals never allow “edit geometry”
* Boundaries are read-only
* Downloads always come from OS / INSPIRE

---

## 7️⃣ How you should replicate council behavior (recommended)

### ✔ Internally (computation)

* Use **EPSG:27700 only**
* INSPIRE as ground truth
* All distances, setbacks, rectangles in meters

### ✔ Visually (frontend)

* Draw property lines **bold**
* Draw imagery **subtle**
* Add disclaimer:

  > “Imagery shown for context only”

### ✔ Never promise pixel-perfect alignment

Councils don’t — and legally can’t.

---

## 8️⃣ Why official apps “look identical”

Because they:

1. Use the **same OS datasets**
2. Use the **same CRS**
3. Avoid Web Mercator entirely
4. Design the UI to **hide residual mismatch**

It’s visual discipline, not better math.

---

## Final takeaway

* Your results are **technically correct**
* Councils don’t “fix” the mismatch — they **design around it**
* Pixel-perfect overlays are **not a UK planning requirement**
* Legal correctness > visual coincidence

If you want, next I can:

* Propose **Mapbox styling rules** that mimic council portals
* Show how to **lock your pipeline to EPSG:27700 end-to-end**
* Help you add a **planner-grade disclaimer layer** like councils do
