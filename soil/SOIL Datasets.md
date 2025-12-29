Here are the **real, authoritative UK soil datasets you can download and use in your Python workflow**, along with *where* and *how* to access them:

---

## ✅ 1) **UK Compiled Topsoil (UKTS) Raster Dataset** — GeoTIFFs

This is the best free **national topsoil data** (geochemical predictions for many elements) mapped on a ~500 m grid across much of the UK. ([Data.gov.uk][1])

### 📥 What you can get

* ~82 GeoTIFF files covering predicted soil chemistry (41 elements)
* Interpolated topsoil concentration values

### 📍 Where to download

The UKTS raster collection is listed on **data.gov.uk** with metadata and links pointing to the **BGS OpenGeoscience** portal. ([Data.gov.uk][1])

❗ **Important:** Direct “one-click” data URLs aren’t in the data.gov.uk listings — you’ll get them from the BGS *OpenGeoscience downloads* page.

👉 Go here to access the full downloadable dataset:

🔗 **[UKTS raster dataset on BGS OpenGeoscience downloads](https://www2.bgs.ac.uk/nationalgeosciencedatacentre/citedData/catalogue/76b69adf-699b-4032-a751-2db0991d55f6.html?utm_source=chatgpt.com)** ([British Geological Survey][2])

**Citation info:**
British Geological Survey. UK Compiled Topsoil raster dataset (UKTSraster). DOI: 10.5285/76b69adf-699b-4032-a751-2db0991d55f6. ([British Geological Survey][2])

You’ll find multiple GeoTIFFs for different elements — you can choose those you need (e.g., Si, Mg, P, etc.) and use them as input rasters in your Python pipeline.

---

## ✅ 2) **Soil Parent Material Model (1 km resolution)**

This dataset describes the **parent material underlying soils** (e.g., chalk, clay, sand/gravel, glacial deposits). A reduced-resolution (1 km) version is *free* to use. ([bgs.ac.uk][3])

### 📍 How to get it

🔗 **[Soil Parent Material Model (1 km free version) download page at BGS](https://www.bgs.ac.uk/download/esri-soil-parent-material-model-1km-resolution/?utm_source=chatgpt.com)** ([bgs.ac.uk][3])

This provides a **GIS dataset (e.g., ESRI shapefile/GPKG)** representing soil parent materials at regular grid resolution — good for texture baseline.

---

## 🗺️ 3) **Visual & WMS Access (optional)**

If you prefer to sample data on the fly or integrate into GIS:

📍 **BGS Soil Property WMS service**

```
https://map.bgs.ac.uk/arcgis/services/UKSO/UKSO_BGS/MapServer/WMSServer
```

Includes layers for:

* Soil parent material
* Soil texture
* Profile soil properties
  …and more. ([bgs.ac.uk][4])

This is excellent for programmatic sampling via OGC WMS GetFeatureInfo.

---

## 📌 Notes & Tips

### ⚠️ Licensing & Use

Both the **UKTS raster data** and **soil parent material model** are free for public & commercial use under UK Open Government Licence *but* require acknowledgement of BGS materials where reproduced. ([Data.gov.uk][1])

### 🧠 Combining Datasets

* **UKTS rasters:** best for element concentrations and indirectly inferring texture (if combined with machine learning or rules). ([Data.gov.uk][1])
* **Parent Material model:** provides texture-related geological context. ([bgs.ac.uk][3])

Together, they allow robust offline soil lookup in Python.

---

## 🚀 Example Download Strategy in Python

Below is a starter script outline you can use once you know the real download URLs from the **BGS OpenGeoscience** pages:

```python
import requests
from pathlib import Path

DATA_DIR = Path("soil_data")
DATA_DIR.mkdir(exist_ok=True)

urls = {
    "ukts_clay": "https://your-real-BGS-url/UKTS_clay_percent.tif",
    "ukts_sand": "https://your-real-BGS-url/UKTS_sand_percent.tif",
    "parent_material": "https://your-real-BGS-url/soil_parent_material_1km.gpkg"
}

for name, url in urls.items():
    dst = DATA_DIR / f"{name}.{url.split('.')[-1]}"
    if not dst.exists():
        print(f"Downloading {name} …")
        resp = requests.get(url, stream=True)
        resp.raise_for_status()
        with open(dst, "wb") as f:
            for chunk in resp.iter_content(8192):
                f.write(chunk)
        print(f"Saved to {dst}")
```

---

## 📌 Summary of Real Data Sources

| Dataset                     | Format          | Coverage      | Cost | Notes                                      |
| --------------------------- | --------------- | ------------- | ---- | ------------------------------------------ |
| UK Compiled Topsoil (UKTS)  | GeoTIFF rasters | UK (majority) | Free | 500 m grid geochemistry ([Data.gov.uk][1]) |
| Soil Parent Material (1 km) | GIS vector/GPKG | Great Britain | Free | Texture / parent geology ([bgs.ac.uk][3])  |
| Soil Property WMS           | WMS             | UK            | Free | Useful for live sampling ([bgs.ac.uk][4])  |

---

If you want, I can **find the explicit download links** for each UKTS GeoTIFF (e.g., clay, sand, silt) and the exact soil parent material file format so you can script the download automatically. Just let me know!

[1]: https://www.data.gov.uk/dataset/4147e16c-b6ba-48b1-834a-0cedd7c7ca0e/uk-compiled-topsoil-ukts-raster-dataset?utm_source=chatgpt.com "UK Compiled Topsoil (UKTS) raster dataset - data.gov.uk"
[2]: https://www2.bgs.ac.uk/nationalgeosciencedatacentre/citedData/catalogue/76b69adf-699b-4032-a751-2db0991d55f6.html?utm_source=chatgpt.com "UK Compiled Topsoil raster dataset (UKTSraster) | NGDC Cited Data | National Geoscience Data Centre (NGDC) | Our data | British Geological Survey (BGS)"
[3]: https://www.bgs.ac.uk/download/esri-soil-parent-material-model-1km-resolution/?utm_source=chatgpt.com "ESRI Soil Parent Material Model (1km resolution) - British Geological Survey"
[4]: https://www.bgs.ac.uk/technologies/web-map-services-wms/soil-property-data-wms/?utm_source=chatgpt.com "Soil property data WMS - British Geological Survey"
