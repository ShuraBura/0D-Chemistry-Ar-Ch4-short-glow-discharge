#!/usr/bin/env python3
"""
Find results where CH is in the 0.6-1.4× range (6e8 to 1.4e9)
"""

import os
import json
import glob

CH_TARGET = 1.0e9
CH_MIN = 0.6 * CH_TARGET  # 6e8
CH_MAX = 1.4 * CH_TARGET  # 1.4e9

H_TARGET = 5.18e13
C2_TARGET = 1.3e11
C2H2_TARGET = 4.0e12

result_dirs = glob.glob("optimization_results_*")
results_ch_in_range = []

for result_dir in result_dirs:
    json_files = glob.glob(os.path.join(result_dir, "best_*.json"))
    
    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            ch_density = data.get('all_densities', {}).get('CH', data.get('target_densities', {}).get('CH', 0))
            
            if CH_MIN <= ch_density <= CH_MAX:
                h_density = data.get('all_densities', {}).get('H', data.get('target_densities', {}).get('H', 0))
                c2_density = data.get('all_densities', {}).get('C2', data.get('target_densities', {}).get('C2', 0))
                c2h2_density = data.get('all_densities', {}).get('C2H2', data.get('target_densities', {}).get('C2H2', 0))
                
                results_ch_in_range.append({
                    'file': json_file,
                    'H': h_density,
                    'CH': ch_density,
                    'C2': c2_density,
                    'C2H2': c2h2_density,
                    'H_frac': h_density / H_TARGET,
                    'CH_frac': ch_density / CH_TARGET,
                    'C2_frac': c2_density / C2_TARGET,
                    'C2H2_frac': c2h2_density / C2H2_TARGET if c2h2_density else 0,
                })
        except Exception as e:
            pass

print(f"Found {len(results_ch_in_range)} results with CH in range [0.6-1.4×]")
print()

# Sort by C2H2 density
results_ch_in_range.sort(key=lambda x: x['C2H2'], reverse=True)

print("="*80)
print("RESULTS WITH CH IN RANGE [0.6-1.4×], SORTED BY C2H2:")
print("="*80)

for i, result in enumerate(results_ch_in_range[:15]):
    print(f"\n{i+1}. {result['file']}")
    print(f"   C2H2: {result['C2H2']:.2e} ({result['C2H2_frac']:.2f}×)")
    print(f"   H:    {result['H']:.2e} ({result['H_frac']:.2f}×)")
    print(f"   CH:   {result['CH']:.2e} ({result['CH_frac']:.2f}×) ✓")
    print(f"   C2:   {result['C2']:.2e} ({result['C2_frac']:.2f}×)")
    
    # Check if H and C2 are also in range
    h_ok = 0.6 <= result['H_frac'] <= 1.4
    c2_ok = 0.6 <= result['C2_frac'] <= 1.4
    
    if h_ok and c2_ok:
        print(f"   ★★★ ALL THREE (H, CH, C2) IN RANGE! ★★★")

print()
print("="*80)
print("SUMMARY:")
print("="*80)
print(f"When CH is in range (0.6-1.4×), the highest C2H2 achieved is:")
if results_ch_in_range:
    best = results_ch_in_range[0]
    print(f"  C2H2: {best['C2H2']:.2e} ({best['C2H2_frac']:.2f}× of 4e12 target)")
    print(f"  File: {best['file']}")
else:
    print("  No results found!")

