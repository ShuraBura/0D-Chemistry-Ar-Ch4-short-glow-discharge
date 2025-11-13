"""
Investigate heterogeneous CH loss mechanisms

Questions:
1. Are we missing CH loss to nanoparticles/dust?
2. Are we missing CH loss to polymerization/film formation?
3. What about carbon cluster formation?
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions

def pressure_to_density(pressure_mTorr, T_K=300):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

print("="*80)
print("INVESTIGATE HETEROGENEOUS CH LOSS MECHANISMS")
print("="*80)

print("\n" + "="*80)
print("CURRENT LOSS MECHANISMS IN MODEL:")
print("="*80 + "\n")

P = 500.0
n_total = pressure_to_density(P)

params = {
    'P': P,
    'Te': baseline['Te'],
    'ne': 2.3e9,
    'E_field': baseline['E_field'],
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4',
                'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3',
                'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3',
                'C3H5', 'C4H', 'C3H6', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
}

k = define_rates(params)

# Check for loss mechanisms
print("1. WALL STICKING:")
wall_loss_keys = [key for key in k.keys() if 'stick_CH' in key]
for key in wall_loss_keys:
    print(f"   {key}: {k[key]:.2e} s⁻¹")

print("\n2. GENERAL LOSS TERMS:")
general_loss_keys = [key for key in k.keys() if 'loss_CH' in key]
for key in general_loss_keys:
    print(f"   {key}: {k[key]:.2e} s⁻¹")

print("\n" + "="*80)
print("MISSING HETEROGENEOUS LOSS MECHANISMS:")
print("="*80 + "\n")

print("In real Ar/CH4 plasmas at 500 mTorr with high carbon radical densities,")
print("we expect several heterogeneous loss processes:\n")

print("1. NANOPARTICLE/DUST FORMATION:")
print("   - CH radicals can stick to growing carbon nanoparticles")
print("   - Typical particle densities: 1e6 - 1e10 cm⁻³")
print("   - Typical particle radii: 1-100 nm")
print("   - Sticking coefficient: 0.1 - 1.0")
print("   - Loss rate formula: k_dust = n_dust × π × r² × v_th × α")
print("     where v_th = sqrt(8kT/πm) ≈ 1e5 cm/s for CH")
print()

# Calculate potential dust loss
T = 300  # K
m_CH = 13 * 1.66e-27  # kg (13 amu)
k_B = 1.38e-23  # J/K
v_th = np.sqrt(8 * k_B * T / (np.pi * m_CH)) * 100  # cm/s

print(f"   CH thermal velocity: {v_th:.2e} cm/s")

# Different dust scenarios
dust_scenarios = [
    ("Low dust", 1e6, 10e-7, 0.3),
    ("Moderate dust", 1e8, 50e-7, 0.5),
    ("High dust", 1e10, 100e-7, 0.8),
]

print("\n   Dust loss rate scenarios:")
for name, n_dust, r_dust, alpha in dust_scenarios:
    k_dust = n_dust * np.pi * r_dust**2 * v_th * alpha
    print(f"   {name:20} n={n_dust:.0e} cm⁻³, r={r_dust*1e7:.0f} nm, α={alpha:.1f}")
    print(f"                       → k_dust = {k_dust:.2e} s⁻¹")

print("\n2. POLYMERIZATION/FILM DEPOSITION:")
print("   - CH can incorporate into growing polymer film on walls")
print("   - This is BEYOND simple wall sticking (physisorption)")
print("   - Film growth rate: 1-100 nm/min typical in CVD")
print("   - Effective loss coefficient: 0.01 - 0.1 (reactive sticking)")
print()

# Current wall sticking calculation
v_th_Bohm = np.sqrt(k_B * baseline['Te'] * 1.6e-19 / m_CH) * 100  # cm/s (Bohm velocity)
A_wall = 2 * np.pi * 0.1 * 0.45 + 2 * np.pi * 0.1**2  # cylinder surface area (m²) -> cm²
A_wall *= 1e4
V_discharge = np.pi * 0.1**2 * 0.45 * 1e6  # cm³

print(f"   Discharge geometry:")
print(f"     Radius: 10 cm, Length: 45 cm")
print(f"     Volume: {V_discharge:.2e} cm³")
print(f"     Wall area: {A_wall:.2e} cm²")
print(f"     A/V ratio: {A_wall/V_discharge:.2e} cm⁻¹")

# Current wall loss
stick_CH_current = k.get('stick_CH_9_3', 0.0)
print(f"\n   Current wall sticking: {stick_CH_current:.2e} s⁻¹")

# Enhanced polymerization loss
gamma_poly = 0.01  # reactive sticking coefficient for polymerization
k_poly = gamma_poly * v_th * (A_wall / V_discharge) / 4
print(f"   Polymerization loss (γ={gamma_poly}): {k_poly:.2e} s⁻¹")

gamma_poly_high = 0.1
k_poly_high = gamma_poly_high * v_th * (A_wall / V_discharge) / 4
print(f"   Polymerization loss (γ={gamma_poly_high}): {k_poly_high:.2e} s⁻¹")

print("\n3. CARBON CLUSTER FORMATION:")
print("   - CH + CₙHₘ → Cₙ₊₁Hₘ₊₁ (growth of larger clusters)")
print("   - Model includes C3, C4 clusters but may miss larger ones")
print("   - Once clusters reach C₆-C₁₀, they can nucleate into nanoparticles")
print("   - Loss to large clusters: effectively permanent")
print()

print("4. ION-INDUCED LOSS:")
print("   - Ions hitting walls can catalyze CH incorporation into films")
print("   - Ion-enhanced sticking coefficient: 0.1 - 1.0")
print("   - Ion flux to walls: Γ_i ≈ 0.5 × n_i × u_Bohm")
print()

# Calculate ion flux
n_i_total = 2.3e9 * 3.84  # Ni/Ne ratio from previous results
u_Bohm = np.sqrt(k_B * baseline['Te'] * 1.6e-19 / (1.67e-27 * 40)) * 100  # Ar⁺ Bohm velocity (cm/s)
Gamma_i = 0.5 * n_i_total * u_Bohm  # ions/cm²/s
print(f"   Total ion density: {n_i_total:.2e} cm⁻³")
print(f"   Bohm velocity (Ar⁺): {u_Bohm:.2e} cm/s")
print(f"   Ion flux to walls: {Gamma_i:.2e} cm⁻²s⁻¹")
print(f"   Ion-induced loss (if 1 CH per ion): {Gamma_i * A_wall / V_discharge:.2e} s⁻¹")

print("\n" + "="*80)
print("IMPACT ANALYSIS:")
print("="*80 + "\n")

CH = 6.36e10  # cm⁻³ from previous results
CH_deficit = 3.11e15  # cm⁻³/s

print(f"CH density: {CH:.2e} cm⁻³")
print(f"CH accumulation rate: {CH_deficit:.2e} cm⁻³/s")
print(f"Required additional loss coefficient: {CH_deficit/CH:.2e} s⁻¹")
print()

# Evaluate each mechanism
print("Can these mechanisms close the gap?\n")

print("1. Dust loss:")
for name, n_dust, r_dust, alpha in dust_scenarios:
    k_dust = n_dust * np.pi * r_dust**2 * v_th * alpha
    rate_dust = k_dust * CH
    impact = rate_dust / CH_deficit * 100
    print(f"   {name:20} k={k_dust:.2e} s⁻¹ → rate={rate_dust:.2e} cm⁻³/s ({impact:.2f}% of deficit)")

print("\n2. Enhanced polymerization:")
rate_poly = k_poly * CH
impact_poly = rate_poly / CH_deficit * 100
print(f"   γ=0.01:              k={k_poly:.2e} s⁻¹ → rate={rate_poly:.2e} cm⁻³/s ({impact_poly:.2f}% of deficit)")

rate_poly_high = k_poly_high * CH
impact_poly_high = rate_poly_high / CH_deficit * 100
print(f"   γ=0.1:               k={k_poly_high:.2e} s⁻¹ → rate={rate_poly_high:.2e} cm⁻³/s ({impact_poly_high:.2f}% of deficit)")

print("\n3. Ion-induced loss:")
k_ion = Gamma_i * A_wall / V_discharge
rate_ion = k_ion * CH
impact_ion = rate_ion / CH_deficit * 100
print(f"   1 CH/ion:            k={k_ion:.2e} s⁻¹ → rate={rate_ion:.2e} cm⁻³/s ({impact_ion:.2f}% of deficit)")

print("\n" + "="*80)
print("CONCLUSION:")
print("="*80 + "\n")

print("The current model has:")
print(f"  - Wall sticking: {stick_CH_current:.2e} s⁻¹")
print(f"  - General loss: Limited to simple wall processes")
print()

total_potential = 0
for name, n_dust, r_dust, alpha in dust_scenarios:
    k_dust = n_dust * np.pi * r_dust**2 * v_th * alpha
    total_potential = max(total_potential, k_dust)

total_potential += k_poly_high + k_ion

print(f"Maximum plausible additional loss: {total_potential:.2e} s⁻¹")
print(f"Required to balance: {CH_deficit/CH:.2e} s⁻¹")
print(f"Gap: {(CH_deficit/CH) / total_potential:.1f}× still too low!")
print()

if total_potential * CH / CH_deficit > 0.1:
    print("✓ These mechanisms COULD significantly help (>10% of deficit)")
    print("  → Worth implementing dust/nanoparticle loss terms")
else:
    print("✗ Even with maximum plausible values, only accounts for")
    print(f"  {total_potential * CH / CH_deficit * 100:.1f}% of the deficit")

print("\nRECOMMENDATIONS:")
print("1. Add CH + Ar⁺ → products (ion-neutral reaction)")
print("2. Add dust/nanoparticle loss term with adjustable n_dust, r_dust, α")
print("3. Add enhanced polymerization loss (beyond simple physisorption)")
print("4. Consider that 63× excess CH may be physical equilibrium")
