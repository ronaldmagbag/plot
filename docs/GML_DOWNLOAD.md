You can download **INSPIRE Index Polygons** for a **specific local authority, region, or custom bounding box (BBox)** using the official INSPIRE services. Here’s the practical, no-nonsense way to do it 👇

---

## What are INSPIRE Index Polygons?

They’re cadastral **index polygons** published by **HM Land Registry** under the EU/UK **INSPIRE Directive**.
They represent **registered land parcels (titles)**, not building footprints.

---

## Official Data Access (UK)

### Service

* **INSPIRE Index Polygons (WFS / Atom)**
* Provider: **HM Land Registry**

### Main entry

* INSPIRE Index Polygons dataset (searchable via HMLR INSPIRE portal)

---

## ✅ Option 1 — Download by **Local Authority (Best & Fastest)**

HM Land Registry publishes **pre-split GML files per local authority**.

**Steps**

1. Go to INSPIRE Index Polygons download page
2. Choose **Local Authority**
3. Download:

   ```
   Index_Polygons_<LocalAuthority>.gml
   ```

**Examples**

* `Index_Polygons_West_Sussex.gml`
* `Index_Polygons_Guildford.gml`

👉 This is the **recommended approach** if your area aligns with council boundaries.

---

## ✅ Option 2 — Download by **Bounding Box (BBox)** via WFS

Use the **WFS (Web Feature Service)** to request only polygons inside a BBox.

### WFS endpoint

```
https://services.landregistry.gov.uk/INSPIRE/wfs
```

### Required layer

```
INSPIRE_CadastralParcels
```

---

### Example: BBox request (EPSG:27700 – British National Grid)

```text
SERVICE=WFS
&VERSION=2.0.0
&REQUEST=GetFeature
&TYPENAMES=INSPIRE_CadastralParcels
&BBOX=500000,120000,510000,130000,EPSG:27700
&OUTPUTFORMAT=application/gml+xml; version=3.2
```

✔ Use this if:

* Your region crosses local authorities
* You want **small, precise areas**
* You are building an automated pipeline

---

## ✅ Option 3 — Download by **Region then Clip**

1. Download **Local Authority GML**
2. Clip using:

   * **QGIS**
   * or **GDAL**

```bash
ogr2ogr clipped.gpkg source.gml \
  -clipsrc xmin ymin xmax ymax
```

✔ Best for:

* Repeatable workflows
* plot.json generation
* Offline processing

---

## Coordinate Systems (Important ⚠️)

| CRS            | Use                       |
| -------------- | ------------------------- |
| **EPSG:27700** | UK official (recommended) |
| EPSG:4326      | Lat/Lon (WGS84)           |

INSPIRE UK data is **natively EPSG:27700**.

---

## What You Should Use (Decision Guide)

| Goal             | Method                  |
| ---------------- | ----------------------- |
| One council area | **Local Authority GML** |
| Custom region    | **WFS BBox**            |
| Automation       | **WFS + script**        |
| Visualization    | **QGIS**                |

---

## Common Pitfalls

* ❌ These polygons ≠ building footprints
* ❌ Some land is **unregistered** → no polygon
* ❌ Large BBox requests may time out

---

