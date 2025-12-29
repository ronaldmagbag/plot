# Soil Texture Parsing Logic

## Overview

The `_parse_soil_texture()` function uses a **lookup table** based on the **USDA Soil Texture Triangle** classifications, not a mathematical formula. It maps texture descriptions to estimated clay/sand/silt percentages.

## Value Format

- **Units**: g/kg (grams per kilogram)
- **Range**: 0-1000 g/kg
- **Conversion**: 100 g/kg = 10% (percentage)
- **Example**: 300 g/kg = 30% clay

This format matches SoilGrids API output for consistency.

## USDA Soil Texture Triangle Definitions

The USDA defines soil texture classes with specific percentage ranges:

### Clay-based textures:
- **Clay**: ≥40% clay, <45% sand, <40% silt
- **Sandy Clay**: ≥35% clay, ≥45% sand
- **Silty Clay**: ≥40% clay, ≥40% silt
- **Clay Loam**: 27-40% clay, 20-45% sand

### Sand-based textures:
- **Sand**: ≥85% sand, <10% clay
- **Loamy Sand**: 70-90% sand, <15% clay
- **Sandy Loam**: <20% clay, ≥52% sand, <50% silt

### Silt-based textures:
- **Silt**: ≥80% silt, <12% clay
- **Silt Loam**: ≥50% silt, 12-27% clay

### Balanced textures:
- **Loam**: 7-27% clay, 28-50% silt, <52% sand

## Current Implementation

The function uses **representative values** (typically midpoints or common values) within each USDA range:

| Texture Class | Clay (g/kg) | Sand (g/kg) | Silt (g/kg) | Notes |
|--------------|-------------|-------------|-------------|-------|
| **Clay** | 450 (45%) | 250 (25%) | 300 (30%) | Midpoint of clay range |
| **Sandy Clay** | 400 (40%) | 450 (45%) | 150 (15%) | High sand, high clay |
| **Silty Clay** | 400 (40%) | 100 (10%) | 500 (50%) | High silt, high clay |
| **Heavy Clay** | 550 (55%) | 200 (20%) | 250 (25%) | Very high clay |
| **Clay Loam** | 300 (30%) | 350 (35%) | 350 (35%) | Balanced with clay dominance |
| **Sand** | 100 (10%) | 750 (75%) | 150 (15%) | High sand |
| **Coarse Sand** | 50 (5%) | 850 (85%) | 100 (10%) | Very high sand |
| **Sandy Loam** | 150 (15%) | 600 (60%) | 250 (25%) | Sand dominant |
| **Silt** | 100 (10%) | 100 (10%) | 800 (80%) | High silt |
| **Silt Loam** | 150 (15%) | 200 (20%) | 650 (65%) | Silt dominant |
| **Loam** | 200 (20%) | 400 (40%) | 400 (40%) | Balanced |

## Calculation Method

Since we only have texture **descriptions** (not lab measurements), we:

1. **Parse keywords** from the description (e.g., "sandy clay loam")
2. **Match to USDA class** using keyword hierarchy
3. **Assign representative values** based on USDA ranges

### Example:
```
Input: "Sandy Clay Loam"
→ Contains "clay" AND "sandy" AND "loam"
→ Matches: Sandy Clay Loam class
→ Assigns: clay=300, sand=500, silt=200 (representative values)
```

## Improving Accuracy

To improve accuracy, you could:

1. **Use actual USDA range midpoints**:
   ```python
   # Sandy Clay Loam: 20-35% clay, 45-65% sand
   clay = (20 + 35) / 2 * 10 = 275 g/kg
   sand = (45 + 65) / 2 * 10 = 550 g/kg
   silt = 1000 - 275 - 550 = 175 g/kg
   ```

2. **Use actual BGS data** if available in the shapefile
3. **Cross-reference with UKTS rasters** for precise values
4. **Add more texture class mappings** based on UK-specific classifications

## Current Limitations

- Uses **fixed representative values** rather than ranges
- Doesn't account for **sub-classes** (e.g., "fine sandy loam" vs "coarse sandy loam")
- May not match **UK-specific** soil classifications exactly
- **Best used as fallback** when precise raster data unavailable

## Recommended Approach

For best results:
1. **Primary**: Use UKTS raster data (if available) - provides actual measured values
2. **Secondary**: Use SPMM shapefile texture descriptions with this parser
3. **Tertiary**: Fall back to API queries or defaults

## References

- USDA Soil Texture Triangle: https://www.nrcs.usda.gov/resources/education-and-teaching-materials/soil-texture-calculator
- USDA Soil Survey Manual Chapter 3: https://www.nrcs.usda.gov/sites/default/files/2022-09/SSM-ch3.pdf

