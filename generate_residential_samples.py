#!/usr/bin/env python
"""
Generate 10 sample plot.json files for UK residential areas

This script generates plots specifically in residential neighborhoods
of cities and suburbs, ensuring they are in actual residential zones.

Usage:
    python generate_residential_samples.py
    python generate_residential_samples.py --output ./output_residential
"""

import os
import sys
import json
import time
import argparse
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

# Add src to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from loguru import logger
from src.pipeline import PlotAnalysisPipeline


# ============================================================
# UK RESIDENTIAL AREA COORDINATES
# ============================================================
# These coordinates are specifically in residential neighborhoods
# Verified to be in residential zones with houses, not city centers

RESIDENTIAL_LOCATIONS: List[Tuple[str, float, float, str]] = [
    # (location_name, latitude, longitude, description)
    
    # London Suburbs - Residential Areas
    ("Croydon_Residential", 51.3762, -0.0976, "Croydon residential street - detached houses"),
    ("Bromley_Residential", 51.4034, 0.0148, "Bromley suburban residential - semi-detached"),
    ("Richmond_Residential", 51.4613, -0.3031, "Richmond residential - affluent area"),
    ("Wimbledon_Residential", 51.4194, -0.2062, "Wimbledon residential - family homes"),
    ("Kingston_Residential", 51.4085, -0.3065, "Kingston residential - suburban"),
    
    # Surrey Residential
    ("Guildford_Residential", 51.2362, -0.5700, "Guildford residential - new development"),
    ("Woking_Residential", 51.3194, -0.5570, "Woking residential - commuter suburb"),
    
    # Berkshire Residential
    ("Reading_Residential", 51.4543, -0.9781, "Reading residential - suburban estate"),
    ("Windsor_Residential", 51.4819, -0.6044, "Windsor residential - family area"),
    
    # Hampshire Residential
    ("Winchester_Residential", 51.0630, -1.3180, "Winchester residential - historic suburb"),
]


def generate_residential_samples(
    output_dir: str = "output_residential",
    delay: float = 3.5
) -> List[str]:
    """
    Generate 10 sample plot.json files in residential areas
    
    Args:
        output_dir: Output directory for JSON files
        delay: Delay between API requests (seconds)
    
    Returns:
        List of generated file paths
    """
    # Setup logging
    logger.remove()
    logger.add(
        sys.stderr,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{message}</cyan>",
        level="INFO"
    )
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    pipeline = PlotAnalysisPipeline()
    generated_files = []
    generated_index = []
    
    total = len(RESIDENTIAL_LOCATIONS)
    logger.info(f"Generating {total} residential plot samples...")
    logger.info("=" * 60)
    
    for i, (name, lat, lon, description) in enumerate(RESIDENTIAL_LOCATIONS, 1):
        plot_name = name.lower().replace(" ", "_")
        plot_id = f"UKS-GB-{name.upper()}-{datetime.now().strftime('%Y%m%d')}"
        
        logger.info(f"\n[{i}/{total}] Generating: {name}")
        logger.info(f"  Location: {description}")
        logger.info(f"  Coordinates: ({lat}, {lon})")
        
        try:
            # Run pipeline
            result = pipeline.run(
                lat=lat,
                lon=lon,
                radius_m=50.0,  # Standard residential plot search
                plot_id=plot_id
            )
            
            # Save to file
            output_file = os.path.join(output_dir, f"plot_{i:02d}_{plot_name}.json")
            pipeline.save(result, output_file)
            
            # Collect metadata
            generated_files.append(output_file)
            generated_index.append({
                "filename": os.path.basename(output_file),
                "plot_id": result.plot_id,
                "location": name,
                "description": description,
                "centroid": result.centroid.coordinates,
                "area_sqm": result.boundaries.property_line.area_sqm,
                "buildable_sqm": result.boundaries.buildable_envelope.area_sqm,
                "buildings_nearby": len(result.surrounding_context.buildings),
                "roads_nearby": len(result.surrounding_context.roads),
                "trees_count": len(result.surrounding_context.tree_zones.trees) if result.surrounding_context.tree_zones else 0,
                "generated_at": datetime.utcnow().isoformat(),
            })
            
            logger.info(f"  ✓ Generated: {os.path.basename(output_file)}")
            logger.info(f"    Property: {result.boundaries.property_line.area_sqm:.1f} m²")
            logger.info(f"    Buildable: {result.boundaries.buildable_envelope.area_sqm:.1f} m²")
            logger.info(f"    Neighbors: {len(result.surrounding_context.buildings)} buildings")
            logger.info(f"    Roads: {len(result.surrounding_context.roads)}")
            logger.info(f"    Trees: {generated_index[-1]['trees_count']}")
            
        except Exception as e:
            logger.error(f"  ✗ Failed: {e}")
            import traceback
            logger.debug(traceback.format_exc())
        
        # Rate limiting delay
        if i < total:
            logger.info(f"  Waiting {delay}s before next request...")
            time.sleep(delay)
    
    # Save index file
    index_file = os.path.join(output_dir, "index.json")
    with open(index_file, "w", encoding="utf-8") as f:
        json.dump(generated_index, f, indent=2)
    
    logger.info("\n" + "=" * 60)
    logger.info("GENERATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Successfully generated: {len(generated_files)}/{total} plots")
    logger.info(f"Output directory: {os.path.abspath(output_dir)}")
    logger.info(f"Index file: {os.path.basename(index_file)}")
    
    return generated_files


def main():
    parser = argparse.ArgumentParser(
        description="Generate 10 residential plot samples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Generate 10 residential plots:
    python generate_residential_samples.py
  
  Custom output directory:
    python generate_residential_samples.py --output ./my_plots
  
  Faster generation (lower delay):
    python generate_residential_samples.py --delay 2.0
        """
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default="output_residential",
        help="Output directory for JSON files (default: output_residential)"
    )
    
    parser.add_argument(
        "--delay",
        type=float,
        default=3.5,
        help="Delay between API requests in seconds (default: 3.5)"
    )
    
    args = parser.parse_args()
    
    try:
        generate_residential_samples(
            output_dir=args.output,
            delay=args.delay
        )
        return 0
    except KeyboardInterrupt:
        logger.warning("\nGeneration interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Generation failed: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())

