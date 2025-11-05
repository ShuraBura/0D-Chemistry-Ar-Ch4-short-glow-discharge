#!/usr/bin/env python3
"""
Find the last working version of C2H2 optimization that actually produced results.

The user said: "It was done before and it worked" - referring to achieving C2H2 of 3-5e12.

Let's search through optimization results to find when C2H2 was last successfully high.
"""

import os
import json
import glob

print("="*80)
print("SEARCHING FOR PREVIOUS C2H2 SUCCESS")
print("="*80)
print()

# Find all result directories
result_dirs = glob.glob("optimization_results_*")
print(f"Found {len(result_dirs)} result directories:")
for d in result_dirs:
    print(f"  - {d}")
print()

# Search for best C2H2 results
best_c2h2_results = []

for result_dir in result_dirs:
    json_files = glob.glob(os.path.join(result_dir, "best_*.json"))
    
    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Check if we have C2H2 density
            c2h2_density = None
            if 'all_densities' in data and 'C2H2' in data['all_densities']:
                c2h2_density = data['all_densities']['C2H2']
            elif 'target_densities' in data and 'C2H2' in data['target_densities']:
                c2h2_density = data['target_densities']['C2H2']
            
            if c2h2_density is not None:
                h_density = data.get('all_densities', {}).get('H', data.get('target_densities', {}).get('H', 0))
                ch_density = data.get('all_densities', {}).get('CH', data.get('target_densities', {}).get('CH', 0))
                c2_density = data.get('all_densities', {}).get('C2', data.get('target_densities', {}).get('C2', 0))
                
                best_c2h2_results.append({
                    'file': json_file,
                    'C2H2': c2h2_density,
                    'H': h_density,
                    'CH': ch_density,
                    'C2': c2_density,
                    'C2H2_fraction': c2h2_density / 4.0e12 if c2h2_density else 0,
                })
        except Exception as e:
            pass

# Sort by C2H2 density
best_c2h2_results.sort(key=lambda x: x['C2H2'], reverse=True)

print(f"Found {len(best_c2h2_results)} results with C2H2 data")
print()

print("="*80)
print("TOP 10 C2H2 RESULTS:")
print("="*80)

for i, result in enumerate(best_c2h2_results[:10]):
    print(f"\n{i+1}. {result['file']}")
    print(f"   C2H2: {result['C2H2']:.2e} ({result['C2H2_fraction']:.2f}× of 4e12 target)")
    print(f"   H:    {result['H']:.2e} ({result['H']/5.18e13:.2f}×)")
    print(f"   CH:   {result['CH']:.2e} ({result['CH']/1.0e9:.2f}×)")
    print(f"   C2:   {result['C2']:.2e} ({result['C2']/1.3e11:.2f}×)")

# Find the highest C2H2 that was actually achieved
if best_c2h2_results:
    best = best_c2h2_results[0]
    print()
    print("="*80)
    print(f"HIGHEST C2H2 ACHIEVED: {best['C2H2']:.2e}")
    print("="*80)
    
    # Load full details
    with open(best['file'], 'r') as f:
        data = json.load(f)
    
    print(f"\nFile: {best['file']}")
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
print("CONCLUSION")
print("="*80)
print()
print("The user said they achieved C2H2 of 3-5e12 experimentally.")
print(f"Our best model result: {best_c2h2_results[0]['C2H2']:.2e}")
print(f"That's a gap of {(3e12 / best_c2h2_results[0]['C2H2']):.1f}× to get to experimental levels")
print()
print("Possible explanations:")
print("1. Experimental wall sticking is 10-100× lower than modeled")
print("2. Gas flow allows C2H2 to accumulate (but flow rate 5.77/s << wall sticking)")
print("3. Different chemistry pathways or rate constants in real experiment")
print("4. Temperature gradients or spatial inhomogeneity not captured in 0D model")

