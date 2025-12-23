#!/usr/bin/env python
"""
INSPIRE GML Boundary Finder

Finds the polygon from INSPIRE GML data that contains a given centroid point (lat/lon).
Converts from EPSG:27700 (British National Grid) to EPSG:4326 (WGS84).

Usage:
    python -m src.collectors.boundary.inspire_boundary_finder --lat 51.268535 --lon -0.570979
    python -m src.collectors.boundary.inspire_boundary_finder --lat 51.268535 --lon -0.570979 --gml data/inspires/Adur_District_Council.gml
"""

import os
import sys
import json
import math
import argparse
from pathlib import Path
from typing import Optional, Dict, Any, List

# Add project root to path
# __file__ is at src/collectors/boundary/inspire_boundary_finder.py
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

try:
    import geopandas as gpd
    from shapely.geometry import Point, Polygon
    from pyproj import Transformer
    GEOPANDAS_AVAILABLE = True
except ImportError:
    GEOPANDAS_AVAILABLE = False
    print("Error: geopandas, shapely, and pyproj are required.")
    print("Install with: pip install geopandas shapely pyproj fiona")
    sys.exit(1)

from loguru import logger


class InspireBoundaryFinder:
    """
    Find plot boundaries from INSPIRE GML data
    
    INSPIRE data uses EPSG:27700 (British National Grid) and needs to be
    converted to EPSG:4326 (WGS84) for use with lat/lon coordinates.
    """
    
    def __init__(self, gml_path: str):
        """
        Initialize with path to INSPIRE GML file
        
        Args:
            gml_path: Path to the INSPIRE GML file
        """
        self.gml_path = Path(gml_path)
        if not self.gml_path.exists():
            raise FileNotFoundError(f"GML file not found: {gml_path}")
        
        self.gdf = None
        self.transformer_27700_to_4326 = None
        self.transformer_4326_to_27700 = None
        
        # Load GML file
        self._load_gml()
    
    def _load_gml(self):
        """Load and prepare GML data"""
        logger.info(f"Loading GML file: {self.gml_path}")
        
        try:
            # Read GML file
            self.gdf = gpd.read_file(self.gml_path)
            
            logger.info(f"Loaded {len(self.gdf)} features")
            logger.info(f"CRS: {self.gdf.crs}")
            logger.info(f"Columns: {list(self.gdf.columns)}")
            
            # Ensure CRS is set (should be EPSG:27700 for UK INSPIRE data)
            if self.gdf.crs is None:
                logger.warning("No CRS found in GML, assuming EPSG:27700")
                self.gdf.set_crs("EPSG:27700", inplace=True)
            
            # Create transformers
            self.transformer_27700_to_4326 = Transformer.from_crs(
                "EPSG:27700", "EPSG:4326", always_xy=True
            )
            self.transformer_4326_to_27700 = Transformer.from_crs(
                "EPSG:4326", "EPSG:27700", always_xy=True
            )
            
            # Convert geometries to WGS84 for spatial queries
            logger.info("Converting geometries to WGS84...")
            self.gdf_wgs84 = self.gdf.to_crs("EPSG:4326")
            
            logger.info("GML data loaded successfully")
            
        except Exception as e:
            logger.error(f"Failed to load GML file: {e}")
            raise
    
    def find_boundary(
        self,
        lat: float,
        lon: float,
        buffer_m: float = 0.0
    ) -> Optional[Dict[str, Any]]:
        """
        Find the polygon that contains the given point
        
        Args:
            lat: Latitude in WGS84 (EPSG:4326)
            lon: Longitude in WGS84 (EPSG:4326)
            buffer_m: Buffer distance in meters for point-in-polygon check (default: 0)
        
        Returns:
            Dictionary with boundary data or None if not found
        """
        logger.info(f"Searching for boundary containing point ({lat}, {lon})")
        
        # Create point in WGS84
        point_wgs84 = Point(lon, lat)
        
        # If buffer is specified, create a buffer around the point
        if buffer_m > 0:
            # Convert buffer from meters to degrees (approximate)
            buffer_deg = buffer_m / 111000  # Rough conversion
            search_area = point_wgs84.buffer(buffer_deg)
        else:
            search_area = None
        
        # Find polygons containing the point
        matches = []
        
        for idx, row in self.gdf_wgs84.iterrows():
            geom = row.geometry
            
            if geom is None:
                continue
            
            # Check if point is in polygon
            if search_area:
                # Use buffer intersection
                if geom.intersects(search_area):
                    matches.append((idx, row, geom))
            else:
                # Direct point-in-polygon check
                if geom.contains(point_wgs84):
                    matches.append((idx, row, geom))
        
        if not matches:
            logger.warning(f"No polygon found containing point ({lat}, {lon})")
            return None
        
        # If multiple matches, prefer the smallest polygon (most specific)
        if len(matches) > 1:
            logger.info(f"Found {len(matches)} polygons containing the point, selecting smallest")
            matches.sort(key=lambda x: x[2].area)
        
        # Get the best match
        idx, row, geom = matches[0]
        
        # Convert polygon to coordinates
        coords = self._geometry_to_coordinates(geom)
        
        # Get area in original CRS (EPSG:27700) for accurate area calculation
        original_geom = self.gdf.iloc[idx].geometry
        area_sqm = original_geom.area
        
        # Calculate perimeter
        perimeter_m = self._calculate_perimeter_wgs84(coords, lat)
        
        # Get INSPIRE ID if available
        inspire_id = None
        for col in ['inspire_id', 'INSPIREID', 'inspireId', 'id']:
            if col in row:
                inspire_id = row[col]
                break
        
        logger.info(f"Found boundary: {area_sqm:.2f} m², perimeter: {perimeter_m:.2f} m")
        
        return {
            "type": "Polygon",
            "coordinates": [coords],
            "area_sqm": area_sqm,
            "perimeter_m": perimeter_m,
            "source": "inspire_cadastral",
            "accuracy_m": 0.5,  # INSPIRE cadastral data is very accurate
            "inspire_id": inspire_id,
            "centroid": [lon, lat]
        }
    
    def _geometry_to_coordinates(self, geom) -> List[List[float]]:
        """
        Convert Shapely geometry to coordinate list
        
        Args:
            geom: Shapely geometry (Polygon or MultiPolygon)
        
        Returns:
            List of [lon, lat] coordinate pairs
        """
        if geom.geom_type == "Polygon":
            # Return exterior coordinates as [lon, lat] pairs
            coords = [[coord[0], coord[1]] for coord in geom.exterior.coords]
            # Ensure polygon is closed
            if coords[0] != coords[-1]:
                coords.append(coords[0])
            return coords
        elif geom.geom_type == "MultiPolygon":
            # Take the largest polygon
            largest = max(geom.geoms, key=lambda g: g.area)
            coords = [[coord[0], coord[1]] for coord in largest.exterior.coords]
            if coords[0] != coords[-1]:
                coords.append(coords[0])
            return coords
        else:
            raise ValueError(f"Unsupported geometry type: {geom.geom_type}")
    
    def _calculate_perimeter_wgs84(
        self,
        coords: List[List[float]],
        center_lat: float
    ) -> float:
        """Calculate polygon perimeter in meters using haversine formula"""
        if len(coords) < 2:
            return 0.0
        
        def haversine_distance(lat1, lon1, lat2, lon2):
            """Calculate distance between two points in meters"""
            R = 6371000  # Earth radius in meters
            phi1 = math.radians(lat1)
            phi2 = math.radians(lat2)
            delta_phi = math.radians(lat2 - lat1)
            delta_lambda = math.radians(lon2 - lon1)
            
            a = (math.sin(delta_phi/2)**2 + 
                 math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda/2)**2)
            c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
            
            return R * c
        
        perimeter = 0.0
        for i in range(len(coords) - 1):
            lon1, lat1 = coords[i]
            lon2, lat2 = coords[i + 1]
            perimeter += haversine_distance(lat1, lon1, lat2, lon2)
        
        return perimeter
    
    def get_all_polygons_info(self) -> Dict[str, Any]:
        """Get summary information about all polygons in the GML file"""
        if self.gdf is None:
            return {}
        
        return {
            "total_features": len(self.gdf),
            "crs": str(self.gdf.crs),
            "columns": list(self.gdf.columns),
            "bounds": self.gdf_wgs84.total_bounds.tolist() if hasattr(self, 'gdf_wgs84') else None
        }


def find_boundary_from_gml(
    lat: float,
    lon: float,
    gml_path: str,
    buffer_m: float = 0.0
) -> Optional[Dict[str, Any]]:
    """
    Convenience function to find boundary from GML file
    
    Args:
        lat: Latitude in WGS84
        lon: Longitude in WGS84
        gml_path: Path to INSPIRE GML file
        buffer_m: Buffer distance in meters (default: 0)
    
    Returns:
        Boundary dictionary or None
    """
    finder = InspireBoundaryFinder(gml_path)
    return finder.find_boundary(lat, lon, buffer_m)


def main():
    parser = argparse.ArgumentParser(
        description="Find plot boundary from INSPIRE GML data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m src.collectors.boundary.inspire_boundary_finder --lat 51.268535 --lon -0.570979
  python -m src.collectors.boundary.inspire_boundary_finder --lat 51.268535 --lon -0.570979 --gml data/inspires/Adur_District_Council.gml
  python -m src.collectors.boundary.inspire_boundary_finder --lat 51.268535 --lon -0.570979 --buffer 10
        """
    )
    
    parser.add_argument(
        "--lat",
        type=float,
        required=True,
        help="Latitude of point (WGS84)"
    )
    
    parser.add_argument(
        "--lon",
        type=float,
        required=True,
        help="Longitude of point (WGS84)"
    )
    
    parser.add_argument(
        "--gml",
        type=str,
        default="data/inspires/Adur_District_Council.gml",
        help="Path to INSPIRE GML file (default: data/inspires/Adur_District_Council.gml)"
    )
    
    parser.add_argument(
        "--buffer",
        type=float,
        default=0.0,
        help="Buffer distance in meters for point search (default: 0)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        help="Output JSON file path (optional)"
    )
    
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger.remove()
    level = "DEBUG" if args.verbose else "INFO"
    logger.add(
        sys.stderr,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level=level
    )
    
    try:
        # Find boundary
        boundary = find_boundary_from_gml(
            lat=args.lat,
            lon=args.lon,
            gml_path=args.gml,
            buffer_m=args.buffer
        )
        
        if boundary:
            logger.info("\n" + "=" * 60)
            logger.info("Boundary Found")
            logger.info("=" * 60)
            logger.info(f"Source: {boundary['source']}")
            logger.info(f"Area: {boundary['area_sqm']:.2f} m²")
            logger.info(f"Perimeter: {boundary['perimeter_m']:.2f} m")
            logger.info(f"Accuracy: {boundary['accuracy_m']:.2f} m")
            if boundary.get('inspire_id'):
                logger.info(f"INSPIRE ID: {boundary['inspire_id']}")
            
            coords = boundary.get("coordinates", [[]])[0]
            if coords:
                logger.info(f"Coordinates: {len(coords)} points")
                logger.info(f"First point: [{coords[0][0]:.6f}, {coords[0][1]:.6f}]")
                logger.info(f"Last point: [{coords[-1][0]:.6f}, {coords[-1][1]:.6f}]")
            
            # Save to file if requested
            if args.output:
                output_path = Path(args.output)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                with open(output_path, 'w') as f:
                    json.dump(boundary, f, indent=2)
                logger.info(f"\nSaved to: {output_path}")
            else:
                # Print JSON to stdout
                print("\n" + json.dumps(boundary, indent=2))
            
            return 0
        else:
            logger.error("No boundary found")
            return 1
            
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

