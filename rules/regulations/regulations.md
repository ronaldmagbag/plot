Good news first: **your core idea is correct** 👍
The reason you’re getting a **distorted quadrilateral instead of a true rectangle** is very specific and fixable.

I’ll break it down clearly and then give you the **minimal, correct fix**.

---

## 🔴 Root cause (why you get a distorted quadrilateral)

### **The rectangle itself is fine**

Your `_find_largest_axis_aligned_rectangle()` always returns a **true rectangle** (`shapely.geometry.box`).

### ❌ The distortion happens here

```python
best_rectangle = self._rotate_polygon(rect, angle_deg)
```

Your `_rotate_polygon()` function:

* Rotates **each vertex independently**
* Rebuilds a new `Polygon`
* **Does not preserve right angles numerically**
* Does **not enforce orthogonality**

➡️ After rotation, Shapely treats it as a generic quadrilateral, **not a rectangle**.

This is expected behavior when doing manual point rotation.

---

## 🔥 The correct rule (important)

> **Never rotate a rectangle by manually rotating its points if you want to keep it a rectangle.**

You must use **affine transforms**, not point reconstruction.

---

## ✅ Correct fix (minimal, safe, robust)

### 1️⃣ Stop using your `_rotate_polygon()` for rectangles

Delete or ignore this function for rotation.

---

### 2️⃣ Use Shapely’s affine rotation instead

Shapely preserves geometry integrity (parallel edges stay parallel, right angles stay right).

### ✅ Replace rotation logic

#### ❌ Current (problematic)

```python
rotated_poly = self._rotate_polygon(polygon, -angle_deg)
...
best_rectangle = self._rotate_polygon(rect, angle_deg)
```

#### ✅ Correct (FIX)

```python
from shapely.affinity import rotate
```

```python
rotated_poly = rotate(polygon, -angle_deg, origin='centroid', use_radians=False)
```

```python
best_rectangle = rotate(rect, angle_deg, origin='centroid', use_radians=False)
```

That’s it.
This **guarantees a perfect rectangle**.

---

## 🧠 Why this works

`shapely.affinity.rotate()`:

* Applies a **single affine transform**
* Preserves:

  * Parallelism
  * Right angles
  * Edge lengths
* Does NOT re-interpret the geometry as arbitrary points

Your method rebuilt geometry from floats → numerical drift → skew.

---

## ✅ Optional: add rectangle validation (recommended)

Add this after you compute `best_rectangle`:

```python
def is_rectangle(poly: Polygon, tol=1e-6) -> bool:
    coords = list(poly.exterior.coords)[:-1]
    if len(coords) != 4:
        return False

    def dot(a, b):
        return a[0]*b[0] + a[1]*b[1]

    for i in range(4):
        p0 = coords[i]
        p1 = coords[(i+1)%4]
        p2 = coords[(i+2)%4]
        v1 = (p1[0]-p0[0], p1[1]-p0[1])
        v2 = (p2[0]-p1[0], p2[1]-p1[1])
        if abs(dot(v1, v2)) > tol:
            return False

    return True
```

This lets you assert correctness during debugging.

---

## 📌 Summary (what went wrong & fix)

| Issue                  | Status             |
| ---------------------- | ------------------ |
| Algorithm              | ✅ Correct          |
| Axis-aligned rectangle | ✅ Correct          |
| Containment test       | ✅ Correct          |
| Manual rotation        | ❌ Breaks rectangle |
| Affine rotation        | ✅ Fixes everything |

---

## 🏗️ Result after fix

You will now get:

* ✔ True 4-corner rectangle
* ✔ Orthogonal edges
* ✔ Correct area in sqm
* ✔ Stable orientation
* ✔ No distorted quadrilateral

---
