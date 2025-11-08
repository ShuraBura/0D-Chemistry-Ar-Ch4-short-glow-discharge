"""
analyze_species_balance.py - Complete production/consumption analysis for key species

Usage: python3 analyze_species_balance.py <result.json>

Analyzes H, CH, and C2 production and consumption pathways to understand
chemistry balance and identify tuning opportunities.
"""

import json
import sys
import numpy as np

def load_result(filename):
    """Load optimization result JSON"""
    with open(filename, 'r') as f:
        return json.load(f)

def analyze_H_balance(result):
    """Analyze all H production and consumption pathways"""
    n = result['all_densities']
    k_opt = result['rate_values']
    Te = result['Te']
    Ne = n['e']

    # Estimated rate constants (from define_rates.py hardcoded values)
    # Electron-impact rates (Te-dependent, using ~1 eV estimates)
    k_e_CH4_CH3_H = 6e-11        # e + CH4 → CH3 + H
    k_e_CH4_CH_H2_H = 3e-11      # e + CH4 → CH + H2 + H
    k_e_CH4_CH_H_H2 = 2e-11      # e + CH4 → CH + H + H2
    k_e_H2_2H = 1e-11            # e + H2 → 2H
    k_e_CH3_CH2_H = 3e-11        # e + CH3 → CH2 + H
    k_e_C2H2_C2H_H = 1e-11       # e + C2H2 → C2H + H
    k_e_C2H4_C2H3_H = 1e-11      # e + C2H4 → C2H3 + H
    k_e_C2H6_C2H5_H = 1e-11      # e + C2H6 → C2H5 + H

    # Neutral reactions
    k_CH2_H_CH_H2 = 1e-11        # CH2 + H → CH + H2
    k_CH_H_C_H2 = 1.2e-10        # CH + H → C + H2
    k_CH_H_CH2 = 1e-10           # CH + H → CH2
    k_CH_H2_CH2_H = 1e-11        # CH + H2 → CH2 + H
    k_H_CH4_CH3_H2 = 1e-30       # H + CH4 → CH3 + H2 (with Ea=0.5eV barrier, ~0 at 400K)
    k_H_C2H4_C2H3_H2 = 1e-30     # H + C2H4 → C2H3 + H2 (with Ea=0.6eV barrier, ~0)

    # CH reactions producing H
    k_CH_CH4_C2H4_H = 1.5e-10    # CH + CH4 → C2H4 + H
    k_CH_CH3_C2H3_H = 1.5e-10    # CH + CH3 → C2H3 + H
    k_CH_CH2_C2H2_H = 1e-11      # CH + CH2 → C2H2 + H
    k_CH_C_C2_H = 1e-10          # CH + C → C2 + H
    k_CH_C2H2_C3H2_H = 1e-10     # CH + C2H2 → C3H2 + H

    # ArStar reactions
    k_ArStar_CH4_CH3_H = 5e-10   # ArStar + CH4 → CH3 + H + Ar
    k_ArStar_H2_2H = 1e-10       # ArStar + H2 → 2H + Ar
    k_ArStar_CH2_CH_H = 8e-11    # ArStar + CH2 → CH + H + Ar

    H_production = {}
    H_consumption = {}

    # ========== H PRODUCTION ==========
    # Electron-impact dissociation
    H_production['e+CH4→CH3+H'] = Ne * n['CH4'] * k_e_CH4_CH3_H
    H_production['e+CH4→CH+H2+H'] = Ne * n['CH4'] * k_e_CH4_CH_H2_H
    H_production['e+CH4→CH+H+H2'] = Ne * n['CH4'] * k_e_CH4_CH_H_H2
    H_production['e+H2→2H'] = 2 * Ne * n['H2'] * k_e_H2_2H  # Produces 2 H
    H_production['e+CH3→CH2+H'] = Ne * n['CH3'] * k_e_CH3_CH2_H
    H_production['e+C2H2→C2H+H'] = Ne * n['C2H2'] * k_e_C2H2_C2H_H
    H_production['e+C2H4→C2H3+H'] = Ne * n['C2H4'] * k_e_C2H4_C2H3_H
    H_production['e+C2H6→C2H5+H'] = Ne * n['C2H6'] * k_e_C2H6_C2H5_H

    # CH reactions (major H source!)
    H_production['CH+CH4→C2H4+H'] = n['CH'] * n['CH4'] * k_CH_CH4_C2H4_H
    H_production['CH+CH3→C2H3+H'] = n['CH'] * n['CH3'] * k_CH_CH3_C2H3_H
    H_production['CH+CH2→C2H2+H'] = n['CH'] * n['CH2'] * k_CH_CH2_C2H2_H
    H_production['CH+C→C2+H'] = n['CH'] * n['C'] * k_CH_C_C2_H
    H_production['CH+C2H2→C3H2+H'] = n['CH'] * n['C2H2'] * k_CH_C2H2_C3H2_H
    H_production['CH+H2→CH2+H'] = n['CH'] * n['H2'] * k_CH_H2_CH2_H

    # ArStar reactions
    H_production['ArStar+CH4→CH3+H+Ar'] = n['ArStar'] * n['CH4'] * k_ArStar_CH4_CH3_H
    H_production['ArStar+H2→2H+Ar'] = 2 * n['ArStar'] * n['H2'] * k_ArStar_H2_2H
    H_production['ArStar+CH2→CH+H+Ar'] = n['ArStar'] * n['CH2'] * k_ArStar_CH2_CH_H

    # ========== H CONSUMPTION ==========
    # Wall sticking
    H_consumption['H_wall_sticking'] = n['H'] * k_opt.get('stick_H_9_1', 0)

    # H abstraction (suppressed by activation barriers)
    H_consumption['H+CH4→CH3+H2'] = n['H'] * n['CH4'] * k_H_CH4_CH3_H2
    H_consumption['H+C2H4→C2H3+H2'] = n['H'] * n['C2H4'] * k_H_C2H4_C2H3_H2

    # H + CH reactions
    H_consumption['CH+H→C+H2'] = n['CH'] * n['H'] * k_CH_H_C_H2
    H_consumption['CH+H→CH2'] = n['CH'] * n['H'] * k_CH_H_CH2
    H_consumption['CH2+H→CH+H2'] = n['CH2'] * n['H'] * k_CH2_H_CH_H2

    return H_production, H_consumption

def analyze_CH_balance(result):
    """Analyze all CH production and consumption pathways"""
    n = result['all_densities']
    k_opt = result['rate_values']
    Te = result['Te']
    Ne = n['e']

    # Rate constants
    k_e_CH4_CH_H2_H = 3e-11      # e + CH4 → CH + H2 + H
    k_e_CH4_CH_H_H2 = 2e-11      # e + CH4 → CH + H + H2
    k_e_CH3_CH2_H = 3e-11        # e + CH3 → CH2 + H (then CH2 → CH)
    k_e_CH_C_H = 6e-11           # e + CH → C + H

    k_CH_CH_C2_H2 = 2.16e-10     # CH + CH → C2 + H2
    k_CH_H_C_H2 = 1.2e-10        # CH + H → C + H2
    k_CH_H_CH2 = 1e-10           # CH + H → CH2
    k_CH_H2_CH2_H = 1e-11        # CH + H2 → CH2 + H
    k_CH_CH4_C2H4_H = 1.5e-10    # CH + CH4 → C2H4 + H
    k_CH_CH4_CH2_CH3 = 1e-11     # CH + CH4 → CH2 + CH3
    k_CH_CH3_C2H3_H = 1.5e-10    # CH + CH3 → C2H3 + H
    k_CH_CH3_C2H2_H2 = 1e-10     # CH + CH3 → C2H2 + H2
    k_CH_CH2_C2H2_H = 1e-11      # CH + CH2 → C2H2 + H
    k_CH_C_C2_H = 1e-10          # CH + C → C2 + H
    k_CH_C2H2_C3H2_H = 1e-10     # CH + C2H2 → C3H2 + H
    k_CH_C2H4_C3H4_H = 1e-10     # CH + C2H4 → C3H4 + H

    k_CH2_H_CH_H2 = 1e-11        # CH2 + H → CH + H2
    k_ArStar_CH2_CH_H = 8e-11    # ArStar + CH2 → CH + H + Ar
    k_CHPlus_e_CH = 3.5e-7       # CHPlus + e → CH

    CH_production = {}
    CH_consumption = {}

    # ========== CH PRODUCTION ==========
    # Electron-impact
    CH_production['e+CH4→CH+H2+H'] = Ne * n['CH4'] * k_e_CH4_CH_H2_H
    CH_production['e+CH4→CH+H+H2'] = Ne * n['CH4'] * k_e_CH4_CH_H_H2

    # From CH2
    CH_production['CH2+H→CH+H2'] = n['CH2'] * n['H'] * k_CH2_H_CH_H2

    # ArStar
    CH_production['ArStar+CH2→CH+H+Ar'] = n['ArStar'] * n['CH2'] * k_ArStar_CH2_CH_H

    # Recombination
    CH_production['CHPlus+e→CH'] = n.get('CHPlus', 0) * Ne * k_CHPlus_e_CH

    # ========== CH CONSUMPTION ==========
    # Wall/drift
    CH_consumption['CH_wall_sticking'] = n['CH'] * k_opt.get('stick_CH_9_3', 0)
    CH_consumption['CH_ambipolar_loss'] = n['CH'] * k_opt.get('loss_CH_11_9', 0)

    # Ionization
    CH_consumption['e+CH→CHPlus+2e'] = Ne * n['CH'] * k_opt.get('e_CH_CHPlus_2e_cm3_2_10', 0)
    CH_consumption['e+CH→C+H'] = Ne * n['CH'] * k_e_CH_C_H

    # CH + CH
    CH_consumption['CH+CH→C2+H2'] = 2 * k_CH_CH_C2_H2 * n['CH']**2  # Consumes 2 CH

    # CH + H
    CH_consumption['CH+H→C+H2'] = n['CH'] * n['H'] * k_CH_H_C_H2
    CH_consumption['CH+H→CH2'] = n['CH'] * n['H'] * k_CH_H_CH2

    # CH + H2
    CH_consumption['CH+H2→CH2+H'] = n['CH'] * n['H2'] * k_CH_H2_CH2_H

    # CH + CH4 (HUGE!)
    CH_consumption['CH+CH4→C2H4+H'] = n['CH'] * n['CH4'] * k_CH_CH4_C2H4_H
    CH_consumption['CH+CH4→CH2+CH3'] = n['CH'] * n['CH4'] * k_CH_CH4_CH2_CH3

    # CH + CH3
    CH_consumption['CH+CH3→C2H3+H'] = n['CH'] * n['CH3'] * k_CH_CH3_C2H3_H
    CH_consumption['CH+CH3→C2H2+H2'] = n['CH'] * n['CH3'] * k_CH_CH3_C2H2_H2

    # CH + CH2
    CH_consumption['CH+CH2→C2H2+H'] = n['CH'] * n['CH2'] * k_CH_CH2_C2H2_H

    # CH + C
    CH_consumption['CH+C→C2+H'] = n['CH'] * n['C'] * k_CH_C_C2_H

    # CH + C2H2
    CH_consumption['CH+C2H2→C3H2+H'] = n['CH'] * n['C2H2'] * k_CH_C2H2_C3H2_H

    # CH + C2H4
    CH_consumption['CH+C2H4→C3H4+H'] = n['CH'] * n['C2H4'] * k_CH_C2H4_C3H4_H

    return CH_production, CH_consumption

def analyze_C2_balance(result):
    """Analyze all C2 production and consumption pathways"""
    n = result['all_densities']
    k_opt = result['rate_values']
    Te = result['Te']
    Ne = n['e']

    # Rate constants
    k_e_C2H2_C2_H2 = 1e-11       # e + C2H2 → C2 + H2
    k_CH_CH_C2_H2 = 2.16e-10     # CH + CH → C2 + H2
    k_CH_C_C2_H = 1e-10          # CH + C → C2 + H
    k_CH_C_C2_H2 = 1e-10         # CH + C → C2 + H2
    k_C2_H_CH_C = 1e-10          # C2 + H → CH + C (reverse)
    k_CH2_CH_C2_H2_H = 1e-11     # CH2 + CH → C2 + H2 + H
    k_CH2_CH2_C2_H2_H2 = 1e-11   # CH2 + CH2 → C2 + 2H2

    C2_production = {}
    C2_consumption = {}

    # ========== C2 PRODUCTION ==========
    # Electron-impact
    C2_production['e+C2H2→C2+H2'] = Ne * n['C2H2'] * k_e_C2H2_C2_H2

    # CH + CH (KEY!)
    C2_production['CH+CH→C2+H2'] = k_CH_CH_C2_H2 * n['CH']**2

    # CH + C
    C2_production['CH+C→C2+H'] = n['CH'] * n['C'] * k_CH_C_C2_H
    C2_production['CH+C→C2+H2'] = n['CH'] * n['C'] * k_CH_C_C2_H2

    # CH2 reactions
    C2_production['CH2+CH→C2+H2+H'] = n['CH2'] * n['CH'] * k_CH2_CH_C2_H2_H
    C2_production['CH2+CH2→C2+2H2'] = k_CH2_CH2_C2_H2_H2 * n['CH2']**2

    # ========== C2 CONSUMPTION ==========
    # Wall/drift
    C2_consumption['C2_wall_sticking'] = n['C2'] * k_opt.get('stick_C2_9_9', 0)
    C2_consumption['C2_loss'] = n['C2'] * k_opt.get('loss_C2_11_3', 0)

    # C2 + H → CH + C (reverse of production)
    C2_consumption['C2+H→CH+C'] = n['C2'] * n['H'] * k_C2_H_CH_C

    # CH + C2 → C3 + H
    C2_consumption['CH+C2→C3+H'] = n['CH'] * n['C2'] * 1e-10

    return C2_production, C2_consumption

def print_balance_table(species, production, consumption, target=None):
    """Print formatted production/consumption table"""
    print(f"\n{'='*80}")
    print(f"{species} BALANCE ANALYSIS")
    print(f"{'='*80}")

    total_prod = sum(production.values())
    total_cons = sum(consumption.values())

    # Current density
    if target:
        print(f"\nTarget: {target:.2e} cm⁻³")

    print(f"\n{'='*80}")
    print(f"{species} PRODUCTION")
    print(f"{'='*80}")
    print(f"{'Pathway':<45s} {'Rate (cm⁻³/s)':>15s} {'%':>8s}")
    print("-" * 80)

    for pathway, rate in sorted(production.items(), key=lambda x: x[1], reverse=True):
        if rate > total_prod * 0.001:  # Show if > 0.1%
            pct = rate / total_prod * 100 if total_prod > 0 else 0
            print(f"{pathway:<45s} {rate:15.3e} {pct:7.1f}%")

    print("-" * 80)
    print(f"{'TOTAL PRODUCTION':<45s} {total_prod:15.3e}")

    print(f"\n{'='*80}")
    print(f"{species} CONSUMPTION")
    print(f"{'='*80}")
    print(f"{'Pathway':<45s} {'Rate (cm⁻³/s)':>15s} {'%':>8s}")
    print("-" * 80)

    for pathway, rate in sorted(consumption.items(), key=lambda x: x[1], reverse=True):
        if rate > total_cons * 0.001:  # Show if > 0.1%
            pct = rate / total_cons * 100 if total_cons > 0 else 0
            print(f"{pathway:<45s} {rate:15.3e} {pct:7.1f}%")

    print("-" * 80)
    print(f"{'TOTAL CONSUMPTION':<45s} {total_cons:15.3e}")

    print(f"\n{'='*80}")
    print(f"STEADY STATE CHECK")
    print(f"{'='*80}")
    print(f"Production:  {total_prod:.3e} cm⁻³/s")
    print(f"Consumption: {total_cons:.3e} cm⁻³/s")
    print(f"Net rate:    {total_prod - total_cons:.3e} cm⁻³/s")
    print(f"P/C ratio:   {total_prod/total_cons:.4f}" if total_cons > 0 else "P/C ratio:   undefined")

    if abs(total_prod/total_cons - 1.0) < 0.05:
        print("✓ At steady state (within 5%)")
    else:
        print(f"⚠ NOT at steady state ({abs(total_prod/total_cons - 1.0)*100:.1f}% imbalance)")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 analyze_species_balance.py <result.json>")
        sys.exit(1)

    filename = sys.argv[1]
    result = load_result(filename)

    # Print header
    print("=" * 80)
    print("SPECIES BALANCE ANALYSIS")
    print("=" * 80)
    print(f"\nFile: {filename}")
    print(f"Te: {result['Te']:.3f} eV")
    print(f"Ne: {result['Ne']:.3e} cm⁻³")
    print(f"E: {result['E_field']:.1f} V/cm")

    n = result['all_densities']
    print(f"\nCurrent densities:")
    print(f"  H:  {n['H']:.3e} cm⁻³ (target: 2.40e14)")
    print(f"  CH: {n['CH']:.3e} cm⁻³ (target: 1.00e09)")
    print(f"  C2: {n['C2']:.3e} cm⁻³ (target: 5.60e11)")

    # Analyze each species
    H_prod, H_cons = analyze_H_balance(result)
    print_balance_table("H", H_prod, H_cons, target=2.4e14)

    CH_prod, CH_cons = analyze_CH_balance(result)
    print_balance_table("CH", CH_prod, CH_cons, target=1.0e9)

    C2_prod, C2_cons = analyze_C2_balance(result)
    print_balance_table("C2", C2_prod, C2_cons, target=5.6e11)

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")

    total_H_prod = sum(H_prod.values())
    total_CH_prod = sum(CH_prod.values())
    total_C2_prod = sum(C2_prod.values())

    print(f"\nMajor H production:   {max(H_prod.items(), key=lambda x: x[1])[0]}")
    print(f"Major H consumption:  {max(H_cons.items(), key=lambda x: x[1])[0]}")
    print(f"\nMajor CH production:  {max(CH_prod.items(), key=lambda x: x[1])[0]}")
    print(f"Major CH consumption: {max(CH_cons.items(), key=lambda x: x[1])[0]}")
    print(f"\nMajor C2 production:  {max(C2_prod.items(), key=lambda x: x[1])[0]}")
    print(f"Major C2 consumption: {max(C2_cons.items(), key=lambda x: x[1])[0]}")

if __name__ == '__main__':
    main()
