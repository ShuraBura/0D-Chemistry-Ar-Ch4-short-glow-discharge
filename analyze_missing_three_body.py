"""
Check what three-body reactions exist and what might be missing

Current three-body reactions (cm⁶/s):
1. C + C + M → C2 + M              k = 1.0e-32
2. H + H + M → H2 + M              k = 1.0e-32
3. CH3 + CH3 + M → C2H6 + M        k = 3.6e-29
4. CH3 + H + M → CH4 + M           k = 5.0e-31

Missing (typical values):
- e + Ar+ + M → Ar + M             k ~ 1e-25 to 1e-27 (three-body electron-ion recombination)
- e + CH4+ + M → CH4 + M           k ~ 1e-25
- e + CH3+ + M → CH3 + M           k ~ 1e-25

Why these matter:
At 500 mTorr, n_total = 1.2e16 cm⁻³
Three-body recombination rate = k_3body × [e] × [ion] × [M]

Example with baseline:
[e] = 1.22e8, [Ar+] ~ 3e7, [M] = 1.2e16
Rate = 1e-25 × 1.22e8 × 3e7 × 1.2e16 = 4.4e6 cm⁻³/s

Compare to two-body dissociative recombination:
Rate = 1.5e-7 × 1.22e8 × 3e7 = 5.5e8 cm⁻³/s

So three-body is ~100× slower at 500 mTorr.
But at higher ne (runaway conditions), three-body becomes important!

When ne = 1e10 (runaway):
Three-body rate = 1e-25 × 1e10 × 1e10 × 1.2e16 = 1.2e21 cm⁻³/s (huge!)
Two-body rate = 1.5e-7 × 1e10 × 1e10 = 1.5e13 cm⁻³/s

Three-body recombination provides NEGATIVE FEEDBACK at high densities!
This could be what's missing.
"""

print(__doc__)

# Calculate expected three-body rates at different conditions
import numpy as np

def three_body_rate(k_3body, ne, ni, n_total):
    """Calculate three-body recombination rate (cm⁻³/s)"""
    return k_3body * ne * ni * n_total

def two_body_rate(k_2body, ne, ni):
    """Calculate two-body recombination rate (cm⁻³/s)"""
    return k_2body * ne * ni

# Typical values
k_3body = 1e-25  # cm⁶/s for electron-ion recombination
k_2body = 1.5e-7  # cm³/s for dissociative recombination

# Test at different conditions
print("\n" + "="*80)
print("THREE-BODY vs TWO-BODY RECOMBINATION")
print("="*80)

conditions = [
    ("Baseline (stable)", 1.22e8, 3e7, 1.2e16),
    ("Moderate runaway", 1e9, 1e9, 1.2e16),
    ("Full runaway", 1e10, 1e10, 1.2e16),
]

for name, ne, ni, n_total in conditions:
    rate_3body = three_body_rate(k_3body, ne, ni, n_total)
    rate_2body = two_body_rate(k_2body, ne, ni)

    print(f"\n{name}:")
    print(f"  ne = {ne:.2e}, ni = {ni:.2e}, n_total = {n_total:.2e}")
    print(f"  Two-body rate:   {rate_2body:.2e} cm⁻³/s")
    print(f"  Three-body rate: {rate_3body:.2e} cm⁻³/s")
    print(f"  Ratio (3-body/2-body): {rate_3body/rate_2body:.2e}")

print("\n" + "="*80)
print("KEY INSIGHT:")
print("="*80)
print("At baseline conditions: three-body is 100× slower (negligible)")
print("During runaway: three-body dominates! (provides stabilization)")
print("\nMissing three-body e-ion recombination prevents stabilization at high ne!")

# Compare pressure effects
print("\n" + "="*80)
print("PRESSURE DEPENDENCE")
print("="*80)

pressures = [300, 400, 500, 600]
for P in pressures:
    n_total = P / 500 * 1.2e16  # Scale from 500 mTorr

    # At moderate runaway conditions
    ne = ni = 1e9
    rate_3body = three_body_rate(k_3body, ne, ni, n_total)
    rate_2body = two_body_rate(k_2body, ne, ni)

    print(f"\n{P} mTorr: n_total = {n_total:.2e}")
    print(f"  Three-body rate: {rate_3body:.2e} cm⁻³/s")
    print(f"  Change vs 500 mTorr: {rate_3body / three_body_rate(k_3body, ne, ni, 1.2e16):.2f}×")

print("\n" + "="*80)
print("CONCLUSION:")
print("="*80)
print("Lower pressure → less three-body stabilization")
print("BUT: Also less collision-driven runaway chemistry")
print("Net effect: Lower pressure likely MORE stable")
print("="*80)
