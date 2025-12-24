You’re not imagining it — **~2 m southward offset when drawing OSM vectors over Mapbox imagery is real and well-known**, and **it is not a CRS bug in your code**.

Below is the precise explanation, specific to **OSM + Mapbox imagery** and your latitude (~53° N).

---

## 1️⃣ OSM vs Mapbox imagery are *not* the same ground truth

Even though both *appear* to be “WGS84”, they come from **different pipelines**.

### OSM vectors

* Collected from:

  * GPS traces
  * Aerial imagery (various vintages)
  * Manual alignment
* Stored in **EPSG:4326**
* Often **visually nudged** by mappers to fit imagery

### Mapbox imagery

Mapbox satellite imagery:

* Comes from **multiple providers** (Maxar, Airbus, etc.)
* Each tile set has:

  * its own orthorectification model
  * its own vertical datum assumptions
* Rendered in **EPSG:3857**

➡ These two datasets are **not guaranteed to align metrically**.

---

## 2️⃣ Why the offset is consistently southward

At ~53° N (UK):

### A. Vertical datum & orthorectification bias

Satellite imagery is:

* Orthorectified using a **global DEM**
* Slight vertical errors (1–3 m) cause **horizontal shift**

Rule of thumb:

```
1 m height error ≈ 1–2 m horizontal error
```

Bias is often:

* Southward
* Slightly eastward

This matches exactly what you see.

---

### B. Web Mercator scale distortion

Mapbox imagery is rendered in **EPSG:3857**:

At latitude φ:

```
scale = 1 / cos(φ)
```

At 53°:

```
cos(53°) ≈ 0.60 → scale ≈ 1.66
```

So:

* Small reprojection rounding
* Tile snapping
* Floating-point truncation

➡ 1–2 m apparent drift is normal.

---

## 3️⃣ This is documented industry behavior (not your bug)

In UK GIS workflows:

* Ordnance Survey vs aerial imagery:

  * **±2–5 m mismatch accepted**
* Planning authorities:

  * Never approve boundaries based on imagery overlays
* Cadastral truth:

  * Comes from survey / INSPIRE / OS MasterMap

Your result is exactly in the expected range.

---

## 4️⃣ How to *prove* it to yourself (quick tests)

### Test A — draw everything in EPSG:3857

Convert both INSPIRE and OSM vectors to 3857 **before drawing**:

```python
gdf_3857 = gdf.to_crs("EPSG:3857")
```

You will still see the offset → proves it’s **imagery alignment**, not CRS math.

---

### Test B — compare OSM vs OSM

Draw:

* OSM roads
* OSM buildings

They will align perfectly.

---

### Test C — compare INSPIRE vs OS MasterMap (if available)

They will align **near perfectly** (<0.5 m).

---

## 5️⃣ What you should do in your pipeline (important)

### ✔ Treat imagery as **visual context only**

* DO NOT snap geometry to imagery
* DO NOT correct offsets manually

### ✔ Treat INSPIRE / OS-derived data as legal truth

* Property lines
* Setbacks
* Max buildable envelope

### ✔ Add tolerance buffers

For UK planning logic:

```text
legal_tolerance = 3.0 m
```

Use this when:

* testing intersection
* computing distances
* validating setbacks

---

## 6️⃣ When offsets become a *real* problem

You should worry only if:

* Offset > 5–6 m
* Direction flips between tiles
* INSPIRE vs INSPIRE mismatches

Your ~2 m southward drift is **normal**.

---

## Final verdict

* **OSM CRS:** correct
* **Mapbox CRS:** correct
* **Your reprojection:** correct
* **Offset source:** imagery orthorectification + Mercator distortion
* **Observed 2 m south shift:** ✅ expected

This is exactly what professional UK GIS pipelines see every day.
