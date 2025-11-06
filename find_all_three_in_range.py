#!/usr/bin/env python3
"""
Find results where H, CH, AND C2 are ALL in the 0.6-1.4× range
Since C2H2 was NOT measured experimentally, we ignore it!
"""

import os
import json
import glob

H_TARGET = 5.18e13
CH_TARGET = 1.0e9
C2_TARGET = 1.3e11

H_MIN, H_MAX = 0.6 * H_TARGET, 1.4 * H_TARGET
CH_MIN, CH_MAX = 0.6 * CH_TARGET, 1.4 * CH_TARGET
C2_MIN, C2_MAX = 0.6 * C2_TARGET, 1.4 * C2_TARGET

result_dirs = glob.glob("optimization_results_*")
results_all_in_range = []

for result_dir in result_dirs:
    json_files = glob.glob(os.path.join(result_dir, "best_*.json"))
    
    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            h_density = data.get('all_densities', {}).get('H', data.get('target_densities', {}).get('H', 0))
            ch_density = data.get('all_densities', {}).get('CH', data.get('target_densities', {}).get('CH', 0))
            c2_density = data.get('all_densities', {}).get('C2', data.get('target_densities', {}).get('C2', 0))
            
            h_ok = H_MIN <= h_density <= H_MAX
            ch_ok = CH_MIN <= ch_density <= CH_MAX
            c2_ok = C2_MIN <= c2_density <= C2_MAX
            
            if h_ok and ch_ok and c2_ok:
                c2h2_density = data.get('all_densities', {}).get('C2H2', data.get('target_densities', {}).get('C2H2', 0))
                charge = data.get('charge_balance', {}).get('imbalance_percent', 999)
                
                results_all_in_range.append({
                    'file': json_file,
                    'H': h_density,
                    'CH': ch_density,
                    'C2': c2_density,
                    'C2H2': c2h2_density,
                    'H_frac': h_density / H_TARGET,
                    'CH_frac': ch_density / CH_TARGET,
                    'C2_frac': c2_density / C2_TARGET,
                    'charge_imbalance': charge,
                })
        except Exception as e:
            pass

print("="*80)
print("SEARCHING FOR SOLUTIONS WHERE H, CH, AND C2 ARE ALL IN RANGE [0.6-1.4×]")
print("="*80)
print()
print(f"Found {len(results_all_in_range)} results with ALL THREE in range!")
print()

if results_all_in_range:
    # Sort by charge balance (best first)
    results_all_in_range.sort(key=lambda x: x['charge_imbalance'])
    
    print("="*80)
    print("★★★ SUCCESS CASES ★★★")
    print("="*80)
    
    for i, result in enumerate(results_all_in_range):
        print(f"\n{i+1}. {result['file']}")
        print(f"   H:    {result['H']:.2e} ({result['H_frac']:.2f}×) ✓")
        print(f"   CH:   {result['CH']:.2e} ({result['CH_frac']:.2f}×) ✓")
        print(f"   C2:   {result['C2']:.2e} ({result['C2_frac']:.2f}×) ✓")
        print(f"   C2H2: {result['C2H2']:.2e} (unmeasured)")
        print(f"   Charge imbalance: {result['charge_imbalance']:.1f}%")
        
        if result['charge_imbalance'] < 100:
            print(f"   ★★★ EXCELLENT CHARGE BALANCE! ★★★")
    
    print()
    print("="*80)
    print("BEST RESULT:")
    print("="*80)
    best = results_all_in_range[0]
    print(f"\nFile: {best['file']}")
    print(f"\nALL TARGETS MET:")
    print(f"  H:  {best['H']:.2e} ({best['H_frac']:.2f}× target) ✓")
    print(f"  CH: {best['CH']:.2e} ({best['CH_frac']:.2f}× target) ✓")
    print(f"  C2: {best['C2']:.2e} ({best['C2_frac']:.2f}× target) ✓")
    print(f"  Charge imbalance: {best['charge_imbalance']:.1f}%")
    
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
    
else:
    print("NO RESULTS FOUND where all three (H, CH, C2) are in range!")
    print()
    print("This confirms the fundamental coupling constraint:")
    print("  - When CH is in range, H and C2 are too low")
    print("  - When H and C2 are in range, CH is too high")
    print()
    print("Next step: Optimize with wall sticking as free parameters")

