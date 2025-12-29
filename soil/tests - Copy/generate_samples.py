#!/usr/bin/env python
"""
Generate sample plot.json files for South UK residential areas

This script generates 20 sample plots from real coordinates in South UK
suitable for residential building development.

Usage:
    python tests/generate_samples.py
    python tests/generate_samples.py --count 10 --output ./output
"""

import os
import sys
import json
import time
import argparse
from datetime import datetime
from typing import List, Tuple

# Add src to path
# __file__ is at tests/generate_samples.py, so go up one level to project root
from pathlib import Path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from loguru import logger
from src.pipeline import PlotAnalysisPipeline


# ============================================================
# SOUTH UK RESIDENTIAL SAMPLE LOCATIONS
# ============================================================
# These are real coordinates in residential areas across South UK
# Selected for typical residential building plots

SOUTH_UK_SAMPLES: List[Tuple[str, float, float, str]] = [
    # (location_name, latitude, longitude, description)
    
    # Greater London suburbs
    ("Croydon", 51.3762, -0.0982, "South London suburban residential"),
    ("Bromley", 51.4039, 0.0198, "Southeast London residential area"),
    ("Richmond", 51.4613, -0.3037, "Southwest London affluent residential"),
    ("Wimbledon", 51.4214, -0.2064, "South London residential"),
    ("Kingston", 51.4123, -0.3007, "Surrey border residential"),
    
    # Surrey
    ("Guildford", 51.2362, -0.5704, "Surrey market town residential"),
    ("Woking", 51.3162, -0.5600, "Surrey commuter town"),
    ("Epsom", 51.3360, -0.2680, "Surrey residential suburb"),
    
    # Kent
    ("Sevenoaks", 51.2731, 0.1907, "Kent residential town"),
    ("Tunbridge_Wells", 51.1320, 0.2630, "Kent spa town residential"),
    ("Maidstone", 51.2720, 0.5290, "Kent county town"),
    
    # Sussex
    ("Brighton_Hove", 50.8225, -0.1372, "Coastal residential"),
    ("Worthing", 50.8147, -0.3714, "West Sussex coastal"),
    ("Horsham", 51.0640, -0.3270, "West Sussex market town"),
    
    # Hampshire
    ("Southampton", 50.9097, -1.4044, "Hampshire port city residential"),
    ("Winchester", 51.0632, -1.3080, "Hampshire historic city"),
    ("Basingstoke", 51.2667, -1.0876, "Hampshire commuter town"),
    
    # Berkshire
    ("Reading", 51.4543, -0.9781, "Berkshire large town"),
    ("Windsor", 51.4839, -0.6044, "Royal Borough residential"),
    
    # Oxfordshire
    ("Oxford_suburb", 51.7520, -1.2577, "University city suburb"),
]


def generate_samples(
    output_dir: str = "output",
    count: int = 20,
    delay: float = 3.0  # Increased delay to avoid API rate limiting
) -> List[str]:
    """
    Generate sample plot.json files
    
    Args:
        output_dir: Output directory for JSON files
        count: Number of samples to generate (max 20)
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
    
    # Limit count to available samples
    count = min(count, len(SOUTH_UK_SAMPLES))
    samples = SOUTH_UK_SAMPLES[:count]
    
    logger.info(f"=" * 60)
    logger.info(f"PLOT ANALYSIS GENERATOR - South UK Samples")
    logger.info(f"=" * 60)
    logger.info(f"Generating {count} sample plots...")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"API delay: {delay}s between requests")
    logger.info(f"=" * 60)
    
    # Initialize pipeline
    pipeline = PlotAnalysisPipeline()
    
    generated_files = []
    failed = []
    
    for i, (name, lat, lon, desc) in enumerate(samples, 1):
        logger.info(f"\n[{i}/{count}] Processing: {name}")
        logger.info(f"  Location: ({lat}, {lon})")
        logger.info(f"  Description: {desc}")
        
        try:
            # Generate plot analysis
            plot_id = f"UKS-GB-{name.upper().replace(' ', '_')}-{datetime.now().strftime('%Y%m%d')}"
            
            result = pipeline.run(
                lat=lat,
                lon=lon,
                plot_id=plot_id
            )
            
            # Save to file
            filename = f"plot_{i:02d}_{name.lower().replace(' ', '_')}.json"
            filepath = os.path.join(output_dir, filename)
            
            pipeline.save(result, filepath)
            generated_files.append(filepath)
            
            logger.info(f"  ✓ Generated: {filename}")
            logger.info(f"    Plot ID: {result.plot_id}")
            logger.info(f"    Area: {result.boundaries.property_line.area_sqm} sqm")
            logger.info(f"    Buildings nearby: {len(result.surrounding_context.buildings)}")
            logger.info(f"    Roads nearby: {len(result.surrounding_context.roads)}")
            
        except Exception as e:
            logger.error(f"  ✗ Failed: {e}")
            failed.append((name, str(e)))
        
        # Delay to avoid rate limiting
        if i < count:
            logger.info(f"  Waiting {delay}s before next request...")
            time.sleep(delay)
    
    # Summary
    logger.info(f"\n" + "=" * 60)
    logger.info(f"GENERATION COMPLETE")
    logger.info(f"=" * 60)
    logger.info(f"Successfully generated: {len(generated_files)}/{count}")
    
    if failed:
        logger.warning(f"Failed: {len(failed)}")
        for name, error in failed:
            logger.warning(f"  - {name}: {error}")
    
    logger.info(f"\nOutput files in: {os.path.abspath(output_dir)}")
    
    # Create index file
    index_path = os.path.join(output_dir, "index.json")
    index_data = {
        "generated_at": datetime.utcnow().isoformat(),
        "count": len(generated_files),
        "region": "South UK",
        "files": [os.path.basename(f) for f in generated_files],
        "samples": [
            {
                "name": samples[i][0],
                "lat": samples[i][1],
                "lon": samples[i][2],
                "description": samples[i][3],
                "file": os.path.basename(generated_files[i]) if i < len(generated_files) else None
            }
            for i in range(len(samples))
        ]
    }
    
    with open(index_path, "w") as f:
        json.dump(index_data, f, indent=2)
    
    logger.info(f"Index file: {index_path}")
    
    return generated_files


def main():
    parser = argparse.ArgumentParser(
        description="Generate sample plot.json files for South UK residential areas"
    )
    parser.add_argument(
        "--output", "-o",
        default="output",
        help="Output directory (default: output)"
    )
    parser.add_argument(
        "--count", "-c",
        type=int,
        default=20,
        help="Number of samples to generate (default: 20, max: 20)"
    )
    parser.add_argument(
        "--delay", "-d",
        type=float,
        default=2.0,
        help="Delay between API requests in seconds (default: 2.0)"
    )
    
    args = parser.parse_args()
    
    try:
        generated = generate_samples(
            output_dir=args.output,
            count=args.count,
            delay=args.delay
        )
        
        if generated:
            print(f"\n[OK] Successfully generated {len(generated)} plot files")
            sys.exit(0)
        else:
            print("\n[FAILED] No files were generated")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n\nGeneration interrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"\n[ERROR] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

