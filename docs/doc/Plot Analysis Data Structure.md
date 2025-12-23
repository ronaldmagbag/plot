# **Plot Analysis Data Structure**

## **Core Plot Information**

**Plot ID**

* Unique identifier for the parcel

**Centroid Coordinate**

* Geographic center point (lat/lon) of the plot

**Shape Points**

* Boundary polygon vertices defining the plot perimeter

**Plot Type**

* Classification: residential building, multi-unit building, field, industrial building, commercial, civic infrastructure, unknown

---

## **Surrounding Area Context**

*(Plot boundary + 50m radius buffer zone)*

### **Geospatial Data Layers:**

**Building Array**

* List of neighboring structures with polygonal footprints, heights, and building types

**Tree Zones**

* RLE-encoded binary mask indicating vegetation coverage areas

**Road Network**

* Road type classification, centerline geometry (polyline points), and width measurements

**Elevation Map**

* Digital terrain model showing ground height variations

**Soil Type**

* Soil classification identifier (string) for the target plot only

**Water Features**

* RLE-encoded mask showing bodies of water, drainage, or flood zones

---

## **Computed Analysis Layers**

**Shadow Map**

* Solar shadow projection calculated from building heights, sun path, and time of day/season

**Adjacency Analysis**

* Plot boundary segmented into edges, each tagged with adjacent feature type:  
  * Street frontage  
  * Neighboring building  
  * Adjacent parcel  
  * Water boundary  
  * Public/governmental land  
  * Agricultural land

**Setback Lines**

* Legal building envelope polygon showing buildable area after applying required setbacks from property lines

---

## **Access & Existing Conditions**

**Primary Access Point**

* Main entry location determined from road network and plot orientation (if determinable from data)

**Vehicle Access**

* Availability and location of driveway/vehicular entry points

**Existing Structures**

* Buildings currently on the plot: footprint geometry and structural type

---

## **Regulatory Constraints**

**Zoning**

* Land use designation and permitted building types

**Allowed Building Footprint**

* Maximum ground coverage percentage or area

**Height Restrictions**

* Maximum building height (meters/stories) permitted by local code

javascript  
{  
  "_id": ObjectId("507f1f77bcf86cd799439011"),  
    
  *// ============================================================*  
  *// CORE IDENTIFICATION*  
  *// ============================================================*  
  "plot_id": "UKS-GB-2024-001",  
  "created_at": ISODate("2024-12-16T10:30:00Z"),  
  "updated_at": ISODate("2024-12-16T10:30:00Z"),  
  "data_version": "1.0",  
    
  *// ============================================================*  
  *// LOCATION & GEOMETRY*  
  *// ============================================================*  
  "centroid": {  
    "type": "Point",  
    "coordinates": [-0.1276, 51.5074]  *// [longitude, latitude] - GeoJSON format*  
  },  
    
  "plot_type": "residential_building",    
  *// Options: residential_building, building_massive, field,*   
  *//          industrial_building, commercial, civic_infrastructure, unknown*  
    
  *// ============================================================*  
  *// BOUNDARIES - THREE SEPARATE POLYGONS*  
  *// ============================================================*  
  "boundaries": {  
      
    *// 1. PROPERTY LINE - Legal cadastral boundary (outer)*  
    "property_line": {  
      "type": "Polygon",  
      "coordinates": [[  
        [0, 0],      *// southwest corner*  
        [20, 0],     *// southeast corner*  
        [20, 30],    *// northeast corner*  
        [0, 30],     *// northwest corner*  
        [0, 0]       *// close polygon*  
      ]],  
      "area_sqm": 600.0,  
      "perimeter_m": 100.0,  
      "source": "uk_land_registry",  
      "source_date": ISODate("2024-01-15T00:00:00Z"),  
      "accuracy_m": 0.5  
    },  
      
    *// 2. SETBACK LINE - Property line minus regulatory setbacks*  
    "setback_line": {  
      "type": "Polygon",  
      "coordinates": [[  
        [1, 5],      *// southwest: 5m front, 1m side setback*  
        [19, 5],     *// southeast: 5m front, 1m side setback*    
        [19, 25],    *// northeast: 5m rear, 1m side setback*  
        [1, 25],     *// northwest: 5m rear, 1m side setback*  
        [1, 5]       *// close polygon*  
      ]],  
      "area_sqm": 360.0,  *// (20-2) × (30-10) = 18 × 20 = 360*  
      "derived_from": "property_line",  
      "setbacks_applied": {  
        "front_m": 5,  
        "rear_m": 5,  
        "side_east_m": 1,  
        "side_west_m": 1  
      },  
      "regulation_source": "local_planning_code_R1"  
    },  
      
    *// 3. BUILDABLE ENVELOPE - Actual buildable area after ALL constraints*  
    "buildable_envelope": {  
      "type": "Polygon",  
      "coordinates": [[  
        [1, 5],  
        [19, 5],  
        [19, 18],    *// reduced from 25 due to utility easement*  
        [12, 18],    *// notch for utility easement*  
        [12, 25],  
        [1, 25],  
        [1, 5]  
      ]],  
      "area_sqm": 311.0,   *// reduced from setback_line due to constraints*  
      "derived_from": "setback_line",  
      "constraints_applied": [  
        {  
          "type": "utility_easement",  
          "description": "7m × 7m easement on southeast corner for underground utilities",  
          "area_removed_sqm": 49.0,  
          "polygon": {  
            "type": "Polygon",   
            "coordinates": [[[12, 18], [19, 18], [19, 25], [12, 25], [12, 18]]]  
          }  
        }  
      ],  
      "access_corridor": {  
        "required": true,  
        "width_m": 3.5,  
        "description": "Vehicle access from north side to potential garage location",  
        "affects_buildable_area": false  *// corridor is within setback zone*  
      }  
    }  
  },  
    
  *// ============================================================*  
  *// SURROUNDING AREA CONTEXT (50m radius)*  
  *// ============================================================*  
  "surrounding_context": {  
      
    "buildings": [  
      {  
        "id": "neighbor_east_001",  
        "footprint": {  
          "type": "Polygon",  
          "coordinates": [[[25, 5], [35, 5], [35, 25], [25, 25], [25, 5]]]  
        },  
        "height_m": 9.0,  
        "stories": 2,  
        "building_type": "detached_house",  
        "usage": "residential",  
        "distance_to_property_line_m": 5.0,  
        "wall_facing_plot": "west"  
      },  
      {  
        "id": "neighbor_west_001",  
        "footprint": {  
          "type": "Polygon",  
          "coordinates": [[[-15, 5], [-5, 5], [-5, 25], [-15, 25], [-15, 5]]]  
        },  
        "height_m": 8.5,  
        "stories": 2,  
        "building_type": "semi_detached_house",  
        "usage": "residential",  
        "distance_to_property_line_m": 5.0,  
        "wall_facing_plot": "east"  
      },  
      {  
        "id": "neighbor_south_001",  
        "footprint": {  
          "type": "Polygon",  
          "coordinates": [[[2, -20], [18, -20], [18, -5], [2, -5], [2, -20]]]  
        },  
        "height_m": 7.0,  
        "stories": 2,  
        "building_type": "detached_house",  
        "usage": "residential",  
        "distance_to_property_line_m": 5.0,  
        "wall_facing_plot": "north"  
      }  
    ],  
      
    "tree_zones": {  
      "encoding": "rle",  
      "data": "0x15A3B8C2...",  *// RLE compressed binary mask*  
      "resolution_m": 0.5,  
      "bounds": {  
        "type": "Polygon",  
        "coordinates": [[[-50, -50], [70, -50], [70, 80], [-50, 80], [-50, -50]]]  
      },  
      "coverage_percent": 12.5  
    },  
      
    "roads": [  
      {  
        "id": "road_north_001",  
        "name": "Example Road",  
        "type": "local_residential",  
        "centerline": {  
          "type": "LineString",  
          "coordinates": [[-100, 30], [0, 30], [20, 30], [120, 30]]  
        },  
        "width_m": 6.0,  
        "traffic_level": "low",  
        "speed_limit_mph": 20,  
        "has_sidewalk": true  
      }  
    ],  
      
    "elevation_map": {  
      "corner_samples": [  
        {"point": [0, 0], "elevation_m": 50.2},  
        {"point": [20, 0], "elevation_m": 50.5},  
        {"point": [20, 30], "elevation_m": 51.1},  
        {"point": [0, 30], "elevation_m": 50.8}  
      ],  
      "slope_percent": 2.5,  
      "slope_direction": "northeast",  
      "average_elevation_m": 50.65,  
      "terrain_classification": "gentle_slope",  
      *// Optional: full DEM reference if needed for drainage analysis*  
      "dem_reference": {  
        "source": "lidar_composite_dtm",  
        "resolution_m": 1.0,  
        "file_path": "s3://terrain-data/UKS-GB-2024-001-dem.tif"  
      }  
    },  
      
    "water_features": {  
      "encoding": "rle",  
      "data": "0x00000000",  *// no water features in 50m radius*  
      "resolution_m": 1.0,  
      "nearest_water_body": {  
        "type": "stream",  
        "distance_m": 230,  
        "direction": "west"  
      }  
    }  
  },  
    
  *// ============================================================*  
  *// COMPUTED ANALYSIS LAYERS*  
  *// ============================================================*  
  "analysis": {  
      
    "shadow_analysis": {  
      "facade_scores": {  
        *// Scale 0-1: 0 = fully shaded, 1 = full sun*  
        "north": {  
          "winter_avg": 0.15,  
          "summer_avg": 0.45,  
          "annual_avg": 0.30  
        },  
        "south": {  
          "winter_avg": 0.85,  
          "summer_avg": 0.95,  
          "annual_avg": 0.90  
        },  
        "east": {  
          "winter_avg": 0.55,  
          "summer_avg": 0.70,  
          "annual_avg": 0.62  
        },  
        "west": {  
          "winter_avg": 0.50,  
          "summer_avg": 0.65,  
          "annual_avg": 0.58  
        }  
      },  
      "shadow_hours_per_day": {  
        "winter_solstice": 4.5,  
        "summer_solstice": 10.2,  
        "equinox": 7.8  
      },  
      "neighbor_shadow_angles": {  
        "east_building": 60,   *// degrees above horizon where shadow starts*  
        "west_building": 58,  
        "south_building": 18  
      },  
      "best_solar_facade": "south",  
      "computed_at": ISODate("2024-12-16T10:30:00Z"),  
      "sun_path_params": {  
        "latitude": 51.5074,  
        "longitude": -0.1276,  
        "timezone": "Europe/London"  
      }  
    },  
      
    "adjacency": [  
      {  
        "edge_id": "north_edge",  
        "geometry": {  
          "type": "LineString",  
          "coordinates": [[0, 30], [20, 30]]  
        },  
        "length_m": 20.0,  
        "adjacent_to": "street",  
        "street_name": "Example Road",  
        "street_type": "local_residential",  
        "primary_access": true,  
        "noise_level": "medium",  
        "privacy_exposure": "high"  *// street-facing = public exposure*  
      },  
      {  
        "edge_id": "east_edge",  
        "geometry": {  
          "type": "LineString",  
          "coordinates": [[20, 0], [20, 30]]  
        },  
        "length_m": 30.0,  
        "adjacent_to": "neighbor_parcel",  
        "neighbor_id": "neighbor_east_001",  
        "shared_wall_potential": false,  
        "distance_to_neighbor_building_m": 5.0,  
        "privacy_exposure": "high"  *// close neighbor = privacy concern*  
      },  
      {  
        "edge_id": "south_edge",  
        "geometry": {  
          "type": "LineString",  
          "coordinates": [[20, 0], [0, 0]]  
        },  
        "length_m": 20.0,  
        "adjacent_to": "neighbor_parcel",  
        "neighbor_id": "neighbor_south_001",  
        "shared_wall_potential": false,  
        "distance_to_neighbor_building_m": 5.0,  
        "privacy_exposure": "medium"  
      },  
      {  
        "edge_id": "west_edge",  
        "geometry": {  
          "type": "LineString",  
          "coordinates": [[0, 0], [0, 30]]  
        },  
        "length_m": 30.0,  
        "adjacent_to": "neighbor_parcel",  
        "neighbor_id": "neighbor_west_001",  
        "shared_wall_potential": false,  
        "distance_to_neighbor_building_m": 5.0,  
        "privacy_exposure": "high"  
      }  
    ]  
  },  
    
  *// ============================================================*  
  *// ACCESS & EXISTING CONDITIONS*  
  *// ============================================================*  
  "access": {  
    "primary_access_point": {  
      "location": {  
        "type": "Point",  
        "coordinates": [10, 30]  *// midpoint of north edge*  
      },  
      "side": "north",  
      "confidence": "high",  
      "determined_by": "street_adjacency",  
      "alternatives": []  *// no alternative access points*  
    },  
      
    "vehicle_access": {  
      "available": true,  
      "type": "direct_from_street",  
      "driveway_width_required_m": 3.5,  
      "turning_radius_adequate": true  
    },  
      
    "pedestrian_access": {  
      "available": true,  
      "sidewalk_present": true,  
      "grade_accessible": true  *// slope < 5%*  
    }  
  },  
    
  "existing_structures": [  
    {  
      "id": "existing_house_001",  
      "footprint": {  
        "type": "Polygon",  
        "coordinates": [[[6, 10], [14, 10], [14, 20], [6, 20], [6, 10]]]  
      },  
      "type": "detached_house",  
      "condition": "dilapidated",  
      "stories": 2,  
      "height_m": 7.5,  
      "year_built": 1965,  
      "status": "to_demolish",  *// Options: to_demolish, to_preserve, to_renovate, unknown*  
      "demolition_cost_estimate_gbp": 8500,  
      "affects_buildable_area": false  *// will be removed*  
    }  
  ],  
    
  *// ============================================================*  
  *// REGULATORY CONSTRAINTS*  
  *// ============================================================*  
  "regulatory": {  
    "zoning": {  
      "designation": "R1",  
      "description": "Single family residential",  
      "permitted_uses": ["single_family_dwelling", "home_office"],  
      "conditional_uses": ["accessory_dwelling_unit"],  
      "source": "local_planning_authority",  
      "last_verified": ISODate("2024-12-01T00:00:00Z")  
    },  
      
    "building_constraints": {  
      "max_coverage_ratio": 0.40,  *// 40% of property_line area*  
      "max_coverage_sqm": 240.0,   *// 600 × 0.40*  
        
      "max_height_m": 10.0,  
      "max_stories": 2,  
      "max_stories_note": "Plus loft conversion allowed within height limit",  
        
      "setbacks_m": {  
        "front": 5.0,  
        "rear": 5.0,  
        "side_east": 1.0,  
        "side_west": 1.0,  
        "notes": "Side setback increases to 2m if building height exceeds 8m"  
      },  
    },  
    // For now - any structure since we don’t have good examples  
    "additional_restrictions": [  
      {  
        "type": "tree_preservation_order",  
        "status": "none_on_plot",  
        "notes": "Neighbor's protected oak tree at (-8, 15) - root protection zone extends 5m"  
      }  
    ]  
  },  
    
  *// ============================================================*  
  *// SOIL CLASSIFICATION (Need further research)*  
  *// ============================================================*  
  "soil": {  
    "type_id": "clay_heavy",  
    "bearing_capacity_kpa": 150,  
    "drainage": "poor",  
    "foundation_recommendation": "deep_pile_or_engineered_slab",  
    "source": "british_geological_survey",  
  },  
    
  *// ============================================================*  
  *// METADATA & QUALITY FLAGS*  
  *// ============================================================*  
  "data_quality": {  
    "data_sources": [  
      {"type": "cadastral_boundary", "source": "uk_land_registry"},  
      {"type": "satellite_imagery", "source": "sentinel2", "date": "2024-11-15", "resolution_m": 10},  
      {"type": "building_footprints", "source": "os_mastermap"},  
      {"type": "elevation", "source": "lidar_composite_dtm", "accuracy_m": 0.15},  
      {"type": "zoning", "source": "local_planning_authority", "verified": "2024-12-01"}  
    ],  
    "missing_data": [  
      "exact_neighbor_window_positions",  
      "underground_utility_precise_locations"  
    ],  
    "last_validated": ISODate("2024-12-16T10:30:00Z")  
  },  
    
    
  *// ============================================================*  
  *// INDEXES (for MongoDB query performance)*  
  *// ============================================================*  
  *// Create these indexes in MongoDB:*  
  *// db.plots.createIndex({"centroid": "2dsphere"})*  
  *// db.plots.createIndex({"plot_id": 1})*  
  *// db.plots.createIndex({"plot_type": 1})*  
  *// db.plots.createIndex({"boundaries.buildable_envelope.area_sqm": 1})*  
  *// db.plots.createIndex({"regulatory.zoning.designation": 1})*  
  *// db.plots.createIndex({"flags.privacy_risk": 1, "flags.overshadow_risk": 1})*

}

---

