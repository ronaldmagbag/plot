What you want is **polygon simplification with an angle (collinearity) constraint**.
INSPIRE parcel/property-line polygons are very dense because they preserve survey accuracy, but many vertices lie on almost-straight lines.

Below are **practical, robust ways** to remove points where the internal angle is ~180° while keeping the legal boundary shape.

---

## 1️⃣ Angle-based vertex removal (exactly what you described)

### Idea

For each vertex **B** with neighbors **A → B → C**:

* Compute the angle at **B**
* If the angle is close to **180°** (collinear), remove **B**

### Geometry

Angle between vectors **BA** and **BC**:

```
cosθ = (u · v) / (|u||v|)
θ ≈ 180° → cosθ ≈ -1
```

### Python implementation (recommended)

```python
import math
import numpy as np

def angle_deg(a, b, c):
    """Angle ABC in degrees"""
    ba = np.array(a) - np.array(b)
    bc = np.array(c) - np.array(b)

    cos_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)

    return math.degrees(math.acos(cos_angle))


def simplify_polygon_by_angle(coords, angle_threshold=179.0):
    """
    Remove points where angle ~ 180°
    coords: list of (x, y), closed or open
    """
    if coords[0] == coords[-1]:
        closed = True
        coords = coords[:-1]
    else:
        closed = False

    simplified = []
    n = len(coords)

    for i in range(n):
        prev = coords[i - 1]
        curr = coords[i]
        nxt = coords[(i + 1) % n]

        ang = angle_deg(prev, curr, nxt)

        if ang < angle_threshold:
            simplified.append(curr)

    if closed:
        simplified.append(simplified[0])

    return simplified
```

### Recommended threshold

| Threshold | Effect                     |
| --------- | -------------------------- |
| `179.9°`  | Very strict (survey-grade) |
| `179°`    | Good default               |
| `178°`    | Aggressive simplification  |

---

## 2️⃣ Distance + angle (safer for legal boundaries)

Sometimes points are **not perfectly collinear** but still unnecessary.

Add a **minimum offset distance** from the line:

```python
def point_line_distance(p, a, b):
    """Distance from point p to line ab"""
    p, a, b = map(np.array, (p, a, b))
    return np.linalg.norm(np.cross(b - a, a - p)) / np.linalg.norm(b - a)


def simplify_polygon_safe(coords, angle_threshold=179, dist_threshold=0.02):
    """
    dist_threshold in meters (UK INSPIRE usually EPSG:27700)
    """
    if coords[0] == coords[-1]:
        closed = True
        coords = coords[:-1]
    else:
        closed = False

    simplified = []
    n = len(coords)

    for i in range(n):
        prev = coords[i - 1]
        curr = coords[i]
        nxt = coords[(i + 1) % n]

        ang = angle_deg(prev, curr, nxt)
        dist = point_line_distance(curr, prev, nxt)

        if not (ang > angle_threshold and dist < dist_threshold):
            simplified.append(curr)

    if closed:
        simplified.append(simplified[0])

    return simplified
```
