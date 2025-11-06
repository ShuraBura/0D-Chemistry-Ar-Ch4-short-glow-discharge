#!/usr/bin/env python3
"""Find the H~0.49×, C2~0.49×, CH~8× result"""

import os
import json
import glob

H_TARGET = 5.18e13
CH_TARGET = 1.0e9
C2_TARGET = 1.3e11

result_dirs = glob.glob("optimization_results_*")
candidates = []

for result_dir in result_dirs:
    json_files = glob.glob(os.path.join(result_dir, "best_*.json"))
    
    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            h_density = data.get('all_densities', {}).get('H', data.get('target_densities', {}).get('H', 0))
            ch_density = data.get('all_densities', {}).get('CH', data.get('target_densities', {}).get('CH', 0))
            c2_density = data.get('all_densities', {}).get('C2', data.get('target_densities', {}).get('C2', 0))
            
            h_frac = h_density / H_TARGET
            ch_frac = ch_density / CH_TARGET
            c2_frac = c2_density / C2_TARGET
            
            # Look for H and C2 around 0.4-0.5×, CH around 5-12×
            if 0.4 <= h_frac <= 0.55 and 0.4 <= c2_frac <= 0.55 and 5 <= ch_frac <= 12:
                candidates.append({
                    'file': json_file,
                    'H': h_density,
                    'CH': ch_density,
                    'C2': c2_density,
                    'H_frac': h_frac,
                    'CH_frac': ch_frac,
                    'C2_frac': c2_frac,
                })
        except Exception as e:
            pass

print(f"Found {len(candidates)} candidates with H~0.4-0.5×, C2~0.4-0.5×, CH~5-12×")
print()

if candidates:
    # Sort by how close to the described values
    candidates.sort(key=lambda x: abs(x['H_frac'] - 0.49) + abs(x['C2_frac'] - 0.49) + abs(x['CH_frac'] - 8))
    
    print("="*80)
    print("CANDIDATES (sorted by closeness to H~0.49×, C2~0.49×, CH~8×):")
    print("="*80)
    
    for i, result in enumerate(candidates[:10]):
        print(f"\n{i+1}. {result['file']}")
        print(f"   H:  {result['H']:.2e} ({result['H_frac']:.2f}×)")
        print(f"   CH: {result['CH']:.2e} ({result['CH_frac']:.2f}×)")
        print(f"   C2: {result['C2']:.2e} ({result['C2_frac']:.2f}×)")
    
    print()
    print("="*80)
    print("BEST MATCH:")
    print("="*80)
    best = candidates[0]
    print(f"\nFile: {best['file']}")
    print(f"\nDensities:")
    print(f"  H:  {best['H']:.2e} ({best['H_frac']:.2f}× target)")
    print(f"  CH: {best['CH']:.2e} ({best['CH_frac']:.2f}× target)")
    print(f"  C2: {best['C2']:.2e} ({best['C2_frac']:.2f}× target)")
    print()
    print("Distance from target range [0.6-1.4×]:")
    print(f"  H:  {best['H_frac']:.2f}× → need {0.6/best['H_frac']:.2f}× boost to reach 0.6×")
    print(f"  C2: {best['C2_frac']:.2f}× → need {0.6/best['C2_frac']:.2f}× boost to reach 0.6×")
    print(f"  CH: {best['CH_frac']:.2f}× → need {best['CH_frac']/1.4:.2f}× reduction to reach 1.4×")
    
    # Load full details
    with open(best['file'], 'r') as f:
        data = json.load(f)
    
    print(f"\nPlasma parameters:")
    print(f"  Te: {data.get('Te', 'N/A')} eV")
    print(f"  Ne: {data.get('Ne', 'N/A'):.2e} cm^-3")
    print(f"  E_field: {data.get('E_field', 'N/A')} V/cm")
    print(f"  H_drift_gain: {data.get('H_drift_gain', 'N/A'):.2e}")
    
    print(f"\nKey wall sticking rates:")
    rate_values = data.get('rate_values', {})
    for key in ['stick_C2H2_9_11', 'stick_C2_9_9', 'stick_CH_9_3', 'stick_H_9_1']:
        if key in rate_values:
            print(f"  {key}: {rate_values[key]:.1f} /s")
    
    print()
    print("="*80)
    print("STRATEGY TO IMPROVE:")
    print("="*80)
    print()
    print("This result is VERY close! We need:")
    print(f"  1. Increase H by ~{0.6/best['H_frac']:.2f}× (boost H production or reduce H loss)")
    print(f"  2. Increase C2 by ~{0.6/best['C2_frac']:.2f}× (boost C2 production or reduce C2 loss)")
    print(f"  3. Decrease CH by ~{best['CH_frac']/1.4:.2f}× (suppress CH production or increase CH loss)")
    print()
    print("Key actions:")
    print("  - Increase H_drift_gain to boost H")
    print("  - Reduce stick_C2 to keep more C2 in gas phase")
    print("  - Increase stick_CH to remove more CH")
    print("  - Fine-tune chemical reaction rates")

else:
    print("No candidates found. Searching with broader ranges...")

