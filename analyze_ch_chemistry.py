#!/usr/bin/env python3
"""
Analyze CH chemistry to understand why CH is 53x too high

Target: CH = 1.0e9 cm^-3
Current: CH = 5.3e10 cm^-3 (53x too high!)

Strategy:
1. Identify all CH production pathways
2. Identify all CH loss pathways
3. Check which rates are at min/max of literature range
4. Quantify production vs loss rates from simulation
"""

import numpy as np
from rate_database_complete import get_complete_rate_database, get_tunable_rates_for_target
from define_rates import define_rates
from build_reactions import build_reactions


def analyze_ch_pathways():
    """Analyze CH production and loss pathways."""

    # Load complete database
    db = get_complete_rate_database()

    print("=" * 80)
    print("CH CHEMISTRY ANALYSIS")
    print("=" * 80)
    print(f"\nTarget: CH = 1.0e9 cm^-3")
    print(f"Current: CH = 5.3e10 cm^-3 (53x TOO HIGH!)")
    print("\nGoal: Reduce CH by factor of 53 while keeping H and C2 near target")

    # Get CH-relevant rates from database
    print("\n" + "=" * 80)
    print("KEY RATES FOR CH CONTROL (from literature)")
    print("=" * 80)

    ch_tunable = get_tunable_rates_for_target('CH')

    for rate_name, description in ch_tunable.items():
        if rate_name in db:
            rate = db[rate_name]

            # Check if at boundary
            is_production = 'source' in description or 'production' in description
            is_loss = 'loss' in description or 'consumption' in description

            current_val = rate.value
            at_min = abs(current_val - rate.min) / rate.min < 0.01
            at_max = abs(current_val - rate.max) / rate.max < 0.01

            boundary_str = ""
            if at_min:
                boundary_str = " [AT MINIMUM]"
            elif at_max:
                boundary_str = " [AT MAXIMUM]"

            range_factor = rate.max / rate.min if rate.min > 0 else 1.0

            print(f"\n{rate_name}:")
            print(f"  Description: {description}")
            print(f"  Current: {current_val:.2e}")
            print(f"  Range: [{rate.min:.2e}, {rate.max:.2e}] ({range_factor:.1f}x range){boundary_str}")
            print(f"  Source: {rate.source}")
            if rate.flag:
                print(f"  FLAG: {rate.flag}")

    # Detailed analysis of key reactions
    print("\n" + "=" * 80)
    print("DETAILED ANALYSIS OF KEY CH REACTIONS")
    print("=" * 80)

    key_reactions = {
        'CH_CH_C2_H2_cm3_5_4': {
            'reaction': 'CH + CH → C2 + H2',
            'impact': 'CRITICAL: Removes 2 CH, produces 1 C2',
            'notes': 'This is THE key reaction - quadratic in [CH]!'
        },
        'CH_CH3_C2H2_H2_cm3_7_23': {
            'reaction': 'CH + CH3 → C2H2 + H2',
            'impact': 'CH loss pathway, produces C2H2 (precursor to C2)',
            'notes': 'Linear in [CH], depends on [CH3]'
        },
        'CH_CH3_C2H3_H_cm3_7_10': {
            'reaction': 'CH + CH3 → C2H3 + H',
            'impact': 'Alternative CH + CH3 pathway',
            'notes': 'Competes with CH + CH3 → C2H2'
        },
        'loss_CH_11_9': {
            'reaction': 'CH → wall',
            'impact': 'CRITICAL wall loss - currently at MAX of 10x range!',
            'notes': 'Range: 1e3 to 1e4, currently 1e4'
        },
        'e_CH4_CH_H2_H_vib_cm3_1_3': {
            'reaction': 'e + CH4 → CH + H2 + H',
            'impact': 'PRIMARY CH source from electron impact',
            'notes': 'Depends on electron density and CH4'
        }
    }

    for rate_name, info in key_reactions.items():
        if rate_name in db:
            rate = db[rate_name]
            print(f"\n{info['reaction']}")
            print(f"  Rate constant: {rate_name}")
            print(f"  Current value: {rate.value:.2e}")
            print(f"  Literature range: [{rate.min:.2e}, {rate.max:.2e}]")
            print(f"  Range factor: {rate.max/rate.min if rate.min > 0 else 1:.1f}x")
            print(f"  Source: {rate.source}")
            print(f"  Impact: {info['impact']}")
            print(f"  Notes: {info['notes']}")

    # Calculate approximate rates from last simulation
    print("\n" + "=" * 80)
    print("QUANTITATIVE ANALYSIS (from last simulation)")
    print("=" * 80)

    # Load final densities from run_optimized.py output
    # H: 3.48e13, CH: 5.30e10, C2: 7.59e11, CH3: ~1e13, e: ~1e10

    n_H = 3.48e13
    n_CH = 5.30e10
    n_C2 = 7.59e11
    n_CH3 = 1e13  # Approximate
    n_e = 1e10    # Approximate
    n_CH4 = 1e14  # Approximate

    # Get current rate constants
    params = {
        'pressure': 400e-3,  # Torr
        'T': 400,            # K
        'L': 0.05,           # m
        'R': 0.05,           # m
    }
    k = define_rates(params)

    # Calculate reaction rates (cm^-3 s^-1)
    print("\nCH PRODUCTION RATES:")

    r_e_CH4_CH = k['e_CH4_CH_H2_H_vib_cm3_1_3'] * n_e * n_CH4
    print(f"  e + CH4 → CH + H2 + H:  {r_e_CH4_CH:.2e} cm^-3 s^-1")

    r_CH2_H_CH = k['CH2_H_CH_H2_cm3_7_1'] * 1e12 * n_H  # Assume [CH2] ~ 1e12
    print(f"  CH2 + H → CH + H2:      {r_CH2_H_CH:.2e} cm^-3 s^-1 (est.)")

    total_production = r_e_CH4_CH + r_CH2_H_CH
    print(f"  TOTAL PRODUCTION:       {total_production:.2e} cm^-3 s^-1")

    print("\nCH LOSS RATES:")

    r_CH_CH_C2 = k['CH_CH_C2_H2_cm3_5_4'] * n_CH * n_CH
    print(f"  CH + CH → C2 + H2:      {r_CH_CH_C2:.2e} cm^-3 s^-1 ***QUADRATIC***")

    r_CH_CH3_C2H2 = k['CH_CH3_C2H2_H2_cm3_7_23'] * n_CH * n_CH3
    print(f"  CH + CH3 → C2H2 + H2:   {r_CH_CH3_C2H2:.2e} cm^-3 s^-1")

    r_CH_wall = k['loss_CH_11_9'] * n_CH
    print(f"  CH → wall:              {r_CH_wall:.2e} cm^-3 s^-1")

    total_loss = r_CH_CH_C2 + r_CH_CH3_C2H2 + r_CH_wall
    print(f"  TOTAL LOSS:             {total_loss:.2e} cm^-3 s^-1")

    print(f"\nProduction/Loss ratio: {total_production/total_loss:.2f}")
    print(f"CH lifetime (from loss): {n_CH/total_loss:.2e} seconds")

    # Strategic recommendations
    print("\n" + "=" * 80)
    print("STRATEGIC RECOMMENDATIONS")
    print("=" * 80)

    print("""
To reduce CH by 53x, we have these options:

1. REDUCE CH PRODUCTION (electron-impact from CH4):
   - Decrease e_CH4_CH_H2_H_vib_cm3_1_3
   - Range: [3e-11, 9e-11], currently 6e-11 (middle)
   - Could go down to 3e-11 (2x reduction in production)

2. INCREASE CH + CH → C2 + H2 (reaction 5.4):
   - Increase CH_CH_C2_H2_cm3_5_4
   - Range: [1.2e-10, 1.8e-10], currently 2.16e-10 (ABOVE MAX!)
   - Wait, this is ALREADY above literature! Should check MATLAB original
   - QUADRATIC in [CH] - very powerful for high CH densities

3. DECREASE CH WALL LOSS (counterintuitive!):
   - loss_CH_11_9: Range [1e3, 1e4], currently 1e4 (at MAX)
   - Could DECREASE to 1e3
   - Why? Because CH is so high, wall loss saturates
   - Better to let CH + CH → C2 convert it to useful C2

4. INCREASE CH + CH3 → C2H2:
   - CH_CH3_C2H2_H2_cm3_7_23
   - This also helps because C2H2 can dissociate to C2

5. COMPLEX STRATEGY:
   - Need to balance C2 production (target: 1.3e11, current: 7.6e11)
   - CH + CH → C2 is BOTH a CH sink AND C2 source
   - If we reduce CH too much, C2 might drop below target
   - Need multi-parameter optimization

CRITICAL FINDING:
- CH_CH_C2_H2_cm3_5_4 = 2.16e-10 is ABOVE literature max of 1.8e-10!
- Need to check if this is correct in original MATLAB
- If incorrect, fixing this alone might help significantly
""")


if __name__ == '__main__':
    analyze_ch_pathways()
