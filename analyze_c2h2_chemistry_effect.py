#!/usr/bin/env python3
"""
Analyze how boosting C2H2 affects H, CH, and C2

User's hypothesis: Increasing C2H2 might:
1. Increase C2 via C2H2 + H → C2 + H2 + H
2. Decrease CH (or at least not increase it as much)
3. Increase H (the reaction produces H!)
"""

import json

# Load the best H~0.49×, C2~0.46×, CH~8× result
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    data = json.load(f)

print("="*80)
print("ANALYZING C2H2 CHEMISTRY IN BEST H~0.49×, C2~0.46×, CH~8× RESULT")
print("="*80)
print()

# Current state
H = data['all_densities']['H']
CH = data['all_densities']['CH']
C2 = data['all_densities']['C2']
C2H2 = data['all_densities']['C2H2']

print("Current state:")
print(f"  H:    {H:.2e} (0.49× target)")
print(f"  CH:   {CH:.2e} (8.03× target)")
print(f"  C2:   {C2:.2e} (0.46× target)")
print(f"  C2H2: {C2H2:.2e}")
print()

# Analyze C2H2 chemistry
c2h2_chem = data['chemistry_analysis']['C2H2']

print("C2H2 chemistry:")
print(f"  Density: {c2h2_chem['density']:.2e}")
print(f"  Production: {c2h2_chem['production']:.2e} /cm³/s")
print(f"  Loss: {c2h2_chem['loss']:.2e} /cm³/s")
print()

print("Top C2H2 production pathways:")
for rxn, rate in c2h2_chem['top_production'][:5]:
    pct = rate / c2h2_chem['production'] * 100
    print(f"  {rxn:50s} {rate:.2e} ({pct:.1f}%)")
print()

print("Top C2H2 loss pathways:")
for rxn, rate in c2h2_chem['top_loss'][:10]:
    pct = rate / c2h2_chem['loss'] * 100
    print(f"  {rxn:50s} {rate:.2e} ({pct:.1f}%)")
print()

# Key reaction: C2H2 + H → C2 + H2 + H
c2h2_to_c2_rate = None
for rxn, rate in c2h2_chem['top_loss']:
    if 'C2H2_H_C2' in rxn:
        c2h2_to_c2_rate = rate
        print(f"★ KEY REACTION: {rxn}")
        print(f"  Rate: {rate:.2e} /cm³/s")
        print(f"  Fraction of C2H2 loss: {rate/c2h2_chem['loss']*100:.1f}%")
        break
print()

# Wall sticking loss
wall_loss = None
for rxn, rate in c2h2_chem['top_loss']:
    if 'stick_C2H2' in rxn:
        wall_loss = rate
        print(f"WALL STICKING: {rxn}")
        print(f"  Rate: {rate:.2e} /cm³/s")
        print(f"  Fraction of C2H2 loss: {rate/c2h2_chem['loss']*100:.1f}%")
        break
print()

# Gas-phase losses
gas_loss = None
for rxn, rate in c2h2_chem['top_loss']:
    if 'loss_C2H2' in rxn and 'stick' not in rxn:
        gas_loss = rate
        print(f"GAS-PHASE LOSS: {rxn}")
        print(f"  Rate: {rate:.2e} /cm³/s")
        print(f"  Fraction of C2H2 loss: {rate/c2h2_chem['loss']*100:.1f}%")
        break
print()

print("="*80)
print("ANALYSIS: WHAT HAPPENS IF WE BOOST C2H2?")
print("="*80)
print()

if c2h2_to_c2_rate and wall_loss and gas_loss:
    print("Current C2H2 loss breakdown:")
    print(f"  Gas-phase loss: {gas_loss:.2e} ({gas_loss/c2h2_chem['loss']*100:.1f}%)")
    print(f"  Wall sticking:  {wall_loss:.2e} ({wall_loss/c2h2_chem['loss']*100:.1f}%)")
    print(f"  C2H2 + H → C2:  {c2h2_to_c2_rate:.2e} ({c2h2_to_c2_rate/c2h2_chem['loss']*100:.1f}%)")
    print()
    
    print("To boost C2H2, we need to:")
    print("  1. REDUCE wall sticking (currently 32.5% of loss)")
    print("  2. REDUCE gas-phase loss (currently 59.7% of loss)")
    print("  3. This will INCREASE C2H2 density")
    print()
    
    print("Effect of higher C2H2:")
    print("  → C2H2 + H → C2 + H2 + H rate increases")
    print("  → MORE C2 produced ✓")
    print("  → MORE H produced (reaction makes H!) ✓")
    print("  → But C2 + H → CH + C might increase CH ✗")
    print()
    
    # Calculate what happens if C2H2 goes up 3×
    print("SCENARIO: If C2H2 increases 3× (from {:.2e} to {:.2e}):".format(
        c2h2_chem['density'], c2h2_chem['density']*3))
    print(f"  C2H2 + H → C2 rate increases ~3×")
    print(f"  New C2 from C2H2: {c2h2_to_c2_rate*3:.2e} /cm³/s")
    print(f"  Current C2 production: {data['chemistry_analysis']['C2']['production']:.2e}")
    print(f"  Boost factor: {c2h2_to_c2_rate*3/data['chemistry_analysis']['C2']['production']:.2f}×")
    print()
    
    # Check C2 + H → CH
    c2_chem = data['chemistry_analysis']['C2']
    c2_to_ch_rate = None
    for rxn, rate in c2_chem['top_loss']:
        if 'C2_H_CH' in rxn:
            c2_to_ch_rate = rate
            print(f"COUPLING: {rxn}")
            print(f"  Current rate: {rate:.2e} /cm³/s")
            print(f"  Fraction of C2 loss: {rate/c2_chem['loss']*100:.1f}%")
            break
    
    if c2_to_ch_rate:
        print()
        print("CRITICAL COUPLING ANALYSIS:")
        print(f"  If C2 increases, C2 + H → CH rate increases")
        print(f"  Current C2 + H → CH: {c2_to_ch_rate:.2e} ({c2_to_ch_rate/c2_chem['loss']*100:.1f}% of C2 loss)")
        
        ch_chem = data['chemistry_analysis']['CH']
        for rxn, rate in ch_chem['top_production']:
            if 'C2_H_CH' in rxn:
                print(f"  This is {rate/ch_chem['production']*100:.1f}% of CH production")
                print()
                print("  ⚠ WARNING: Boosting C2 will boost CH via this pathway!")
                print("  Need to suppress CH consumption to compensate")
                break

print()
print("="*80)
print("RECOMMENDATION:")
print("="*80)
print()
print("Try moderate C2H2 boost (2-4×) by:")
print("  1. Reduce stick_C2H2 by 2-5×")
print("  2. Reduce loss_C2H2 by 2-3×")
print("  3. Increase stick_CH to suppress CH production")
print("  4. Monitor C2 + H → CH coupling carefully")
print()
print("Expected outcome:")
print("  H:  0.49× → 0.6-0.7× (boosted by C2H2 + H → C2 + H2 + H)")
print("  C2: 0.46× → 0.6-0.8× (boosted by C2H2 + H → C2)")
print("  CH: 8.0× → 2-4× (need to suppress via increased wall sticking)")

