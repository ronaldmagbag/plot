I’ve reviewed the exact code sections you pointed to (loading → coordinate extraction → area), and the **short answer** is:

> **Your projection logic is mostly correct**, but the **~3 m east/south offset is real and expected** given how you are converting and visualising INSPIRE cadastral data against OSM.

Below is a precise breakdown, line-by-line logic check, and *why* you see that shift.

---

## 1️⃣ GML loading & CRS handling (lines ~78–110)

### What your code does

```python
gdf = gpd.read_file(str(gml_path))

if gdf.crs is None:
    gdf.set_crs("EPSG:27700", inplace=True)

self._inspire_gdf_wgs84 = self._inspire_gdf.to_crs("EPSG:4326")
```

### ✔ Correct

* UK INSPIRE cadastral GML **is natively EPSG:27700** (British National Grid)
* Explicitly setting CRS **only when missing** is correct
* Reprojecting to **EPSG:4326 (WGS84)** for OSM-style use is correct

### ⚠ Subtle issue (important)

Some INSPIRE GML files:

* *declare* EPSG:27700
* but internally reference **OSGB36 + OSTN15 grid shifts**

`to_crs("EPSG:4326")` uses PROJ defaults.
Depending on your PROJ version, it may:

* use OSTN02
* or a simplified Helmert transform

➡ This alone can introduce **1–3 m horizontal shift**

---

## 2️⃣ Coordinate extraction (lines ~199–201)

```python
coords = [[coord[0], coord[1]] for coord in geom.exterior.coords]
```

### ✔ Correct

* After `to_crs("EPSG:4326")`, these are:

  ```
  [lon, lat]
  ```
* This matches:

  * GeoJSON
  * Mapbox
  * OSM expectations

No axis flip bug here 👍

---

## 3️⃣ Area calculation (lines ~213–215)

```python
original_geom = self._inspire_gdf.iloc[idx].geometry
area_sqm = original_geom.area
```

### ✔ Correct and best practice

* Area computed in **EPSG:27700 (meters)** → accurate
* You did **not** compute area in WGS84 → very good

No problem here.

---

## 4️⃣ Why you see ~3 m east & south offset on OSM

This is **not a bug in your code**.

### Root causes (stacked)

#### A. Different reference authorities

* **INSPIRE parcels**

  * Survey-grade cadastral geometry
  * Based on OS MasterMap lineage
* **OSM**

  * Digitised from imagery + GPS
  * Often intentionally **offset from legal boundary**
  * Roads/buildings prioritised over cadastral precision

➡ 1–5 m mismatch is *normal* in the UK

---

#### B. OSGB36 → WGS84 transformation accuracy

British National Grid requires **grid-based corrections**:

| Method            | Typical error |
| ----------------- | ------------- |
| Helmert transform | 2–7 m         |
| OSTN02            | ~1 m          |
| OSTN15 (best)     | < 0.2 m       |

Your `to_crs("EPSG:4326")` **may not be using OSTN15**.

➡ That alone explains **east/south bias of ~3 m**

---

#### C. Visual comparison trap

You are:

* Drawing **legal property_line**
* Over **OSM basemap**

But:

* OSM roads are often centerlines
* Buildings are manually aligned
* Parcel boundaries are *not authoritative* in OSM

So alignment will *never* be perfect.

---

## 5️⃣ How to confirm this empirically (recommended)

### Test 1 — compare before reprojection

Plot both in **EPSG:27700**:

```python
self._inspire_gdf.plot()
```

Overlay with:

* OS MasterMap (if available)
* Or INSPIRE-only comparison

➡ You’ll see **no offset**

---

### Test 2 — force high-accuracy transform

If you want best possible WGS84 alignment:

```python
self._inspire_gdf_wgs84 = self._inspire_gdf.to_crs(
    "EPSG:4326",
    authority="EPSG"
)
```

Or explicitly ensure modern PROJ data is installed (OSTN15).

---

## 6️⃣ Recommendation for your pipeline (important)

### ✔ Treat INSPIRE as **ground truth**

Use it for:

* `property_line`
* `parcel_area`
* legal setbacks
* max inscribed rectangle

### ✔ Treat OSM as **context only**

Use it for:

* roads
* nearby buildings
* access logic
* privacy distance metrics

### ❌ Do NOT try to “snap” INSPIRE to OSM

That introduces **legal inaccuracies**.

---

## Final verdict

| Component                   | Status              |
| --------------------------- | ------------------- |
| CRS assumption (EPSG:27700) | ✅ Correct           |
| Reprojection logic          | ✅ Correct           |
| Coordinate extraction       | ✅ Correct           |
| Area calculation            | ✅ Correct           |
| 3 m offset                  | ✅ Expected / normal |

**No projection bug found.**
What you’re seeing is a **known UK geodesy + data-authority mismatch**, not a coding error.
