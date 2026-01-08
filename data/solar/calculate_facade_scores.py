#!/usr/bin/env python3
"""Calculate facade_base_scores from merged PVGIS data"""

import json
from pathlib import Path
from typing import Dict, List

def calculate_facade_scores(merged_json_path: str) -> Dict[str, Dict[str, float]]:
    """
    Calculate facade base scores from PVGIS data.
    
    For each direction:
    - Winter: Average of 3 months with minimum E_m values
    - Summer: Average of 3 months with maximum E_m values
    - Annual: Average of all 12 months E_m values
    
    Then normalize all values by the maximum value so the best facade gets 1.0.
    
    Returns:
        Dict with structure: {direction: {"winter": float, "summer": float}}
    """
    with open(merged_json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Map azimuth to direction name
    azimuth_to_direction = {
        "-90": "east",
        "0": "south",
        "90": "west",
        "-179": "north"
    }
    
    # Collect all E_m values for each direction
    direction_data = {}
    
    for az_key, direction in azimuth_to_direction.items():
        if az_key not in data['outputs']['by_azimuth']:
            continue
        
        monthly_data = data['outputs']['by_azimuth'][az_key]['monthly']['fixed']
        
        # Extract E_m values for all months
        e_m_values = [month['E_m'] for month in monthly_data]
        
        # Winter: average of 3 months with minimum E_m
        sorted_indices = sorted(range(len(e_m_values)), key=lambda i: e_m_values[i])
        winter_months = sorted_indices[:3]
        winter_avg = sum(e_m_values[i] for i in winter_months) / 3
        
        # Summer: average of 3 months with maximum E_m
        summer_months = sorted_indices[-3:]
        summer_avg = sum(e_m_values[i] for i in summer_months) / 3
        
        # Annual: average of all 12 months
        annual_avg = sum(e_m_values) / 12
        
        direction_data[direction] = {
            'winter': winter_avg,
            'summer': summer_avg,
            'annual': annual_avg,
            'all_values': {
                'winter': winter_avg,
                'summer': summer_avg,
                'annual': annual_avg
            }
        }
    
    # Find overall maximum value across all directions and seasons for normalization
    all_values = []
    for d in direction_data.values():
        all_values.extend([d['winter'], d['summer'], d['annual']])
    max_overall = max(all_values) if all_values else 1.0
    
    # Normalize all values by the overall maximum
    facade_scores = {}
    for direction, values in direction_data.items():
        facade_scores[direction] = {
            'winter': values['winter'] / max_overall if max_overall > 0 else 0.0,
            'summer': values['summer'] / max_overall if max_overall > 0 else 0.0,
            'annual': values['annual'] / max_overall if max_overall > 0 else 0.0
        }
    
    return facade_scores, direction_data, {
        'max_overall': max_overall
    }

def main():
    """Main function to calculate and display facade scores"""
    merged_file = Path(__file__).parent / 'PVdata_51.013_-0.105_SA3_crystSi_1kWp_14_90deg_merged.json'
    
    if not merged_file.exists():
        print(f"Error: {merged_file} not found!")
        return
    
    facade_scores, raw_data, max_values = calculate_facade_scores(str(merged_file))
    
    print("=" * 70)
    print("Facade Base Scores Calculation from PVGIS Data")
    print("=" * 70)
    print("\nRaw Values (kWh/month):")
    print("-" * 70)
    for direction in ['north', 'south', 'east', 'west']:
        if direction in raw_data:
            print(f"{direction:8s}: Winter={raw_data[direction]['winter']:6.2f}, "
                  f"Summer={raw_data[direction]['summer']:6.2f}, "
                  f"Annual={raw_data[direction]['annual']:6.2f}")
    
    print(f"\nMaximum Value (for normalization):")
    print(f"  Overall Maximum: {max_values['max_overall']:.2f} kWh/month")
    
    print("\n" + "=" * 70)
    print("Normalized Facade Base Scores (0-1 scale):")
    print("=" * 70)
    print("\nFor config.py:")
    print("-" * 70)
    print("facade_base_scores: Dict[str, Dict[str, float]] = field(default_factory=lambda: {")
    
    for direction in ['north', 'south', 'east', 'west']:
        if direction in facade_scores:
            winter = facade_scores[direction]['winter']
            summer = facade_scores[direction]['summer']
            print(f'    "{direction}": {{"winter": {winter:.3f}, "summer": {summer:.3f}}},')
    
    print("})")
    
    print("\n" + "=" * 70)
    print("Detailed Scores:")
    print("=" * 70)
    for direction in ['north', 'south', 'east', 'west']:
        if direction in facade_scores:
            print(f"{direction:8s}: Winter={facade_scores[direction]['winter']:.3f}, "
                  f"Summer={facade_scores[direction]['summer']:.3f}, "
                  f"Annual={facade_scores[direction]['annual']:.3f}")

if __name__ == "__main__":
    main()

