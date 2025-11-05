#!/usr/bin/env python3
"""
MANUAL PARAMETER SWEEP: Very conservative C2H2 boost

Starting from best_f259.0_Te1.09.json (H=0.49×, C2=0.46×, CH=8×)
Make SMALL adjustments to avoid ODE solver instabilities
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
import os

from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions

TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

os.makedirs('manual_sweep_results', exist_ok=True)

# Load starting point
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)

print("="*80)
print("MANUAL C2H2 BOOST PARAMETER SWEEP")
print("="*80)
print()
print("Starting from: optimization_results_C2H2_boost/best_f259.0_Te1.09.json")
print("  H:    2.53e+13 (0.49×)")
print("  C2:   6.04e+10 (0.46×)")
print("  CH:   8.03e+09 (8.03×)")
print("  C2H2: 2.16e+10")
print()

baseline_rates = start_data['rate_values'].copy()
print("Baseline rates:")
val = baseline_rates.get('stick_C2H2_9_11')
print(f"  stick_C2H2_9_11:  {val:.1f}" if val is not None else "  stick_C2H2_9_11:  N/A")
val = baseline_rates.get('loss_C2H2_11_19')
print(f"  loss_C2H2_11_19:  {val:.1f}" if val is not None else "  loss_C2H2_11_19:  N/A")
val = baseline_rates.get('stick_CH_9_3')
print(f"  stick_CH_9_3:     {val:.1f}" if val is not None else "  stick_CH_9_3:     N/A")
print()


class PlasmaODE:
    def __init__(self, params):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.k = params['k']
        self.tags = params['tags']
        self.H_drift_gain = params.get('H_drift_gain', 3.2e17)
        self.e_idx = self.species.index('e')
        self.H_idx = self.species.index('H')
        self._build_stoich()

    def _build_stoich(self):
        from scipy import sparse
        rows_prod, cols_prod, vals_prod = [], [], []
        rows_react, cols_react, vals_react = [], [], []
        for rxn_idx, reaction in enumerate(self.R):
            prod_species = np.where(reaction.products > 0)[0]
            for sp_idx in prod_species:
                rows_prod.append(sp_idx)
                cols_prod.append(rxn_idx)
                vals_prod.append(reaction.products[sp_idx])
            react_species = np.where(reaction.reactants > 0)[0]
            for sp_idx in react_species:
                rows_react.append(sp_idx)
                cols_react.append(rxn_idx)
                vals_react.append(reaction.reactants[sp_idx])
        self.nr = len(self.R)
        self.prod_matrix = sparse.csr_matrix((vals_prod, (rows_prod, cols_prod)), shape=(self.ns, self.nr))
        self.react_matrix = sparse.csr_matrix((vals_react, (rows_react, cols_react)), shape=(self.ns, self.nr))

    def __call__(self, t, y):
        y = np.maximum(y, 1e-10)
        rates = np.zeros(len(self.R))
        for rxn_idx, reaction in enumerate(self.R):
            rate_constant = self.k.get(reaction.rate, 0.0)
            if 'drift' in self.tags[rxn_idx]:
                rates[rxn_idx] = rate_constant * y[self.e_idx]
            elif 'loss' in self.tags[rxn_idx] or 'stick' in self.tags[rxn_idx]:
                react_species = np.where(reaction.reactants > 0)[0]
                if len(react_species) > 0:
                    rates[rxn_idx] = rate_constant * y[react_species[0]]
            else:
                rate = rate_constant
                react_species = np.where(reaction.reactants > 0)[0]
                for sp_idx in react_species:
                    rate *= y[sp_idx] ** reaction.reactants[sp_idx]
                rates[rxn_idx] = rate
        dydt = self.prod_matrix.dot(rates) - self.react_matrix.dot(rates)
        dydt[self.H_idx] += self.H_drift_gain
        return dydt


def run_simulation(params_dict):
    try:
        params = {'species': [], 'P': 0.4, 'Tgas': 300}
        params.update(params_dict)

        k = define_rates_tunable(params)
        db = get_complete_rate_database()

        for name, val in params_dict.get('rate_values', {}).items():
            if name in k and name in db:
                k[name] = np.clip(val, db[name].min, db[name].max)

        params['k'] = k
        params['R'], params['tags'] = build_reactions(params)

        species = params['species']
        y0 = np.ones(len(species)) * 1e3

        def set_density(name, value):
            try:
                y0[species.index(name)] = value
            except ValueError:
                pass

        set_density('e', params['ne'])
        set_density('Ar', 0.85 * 9.66e15)
        set_density('CH4', 0.15 * 9.66e15)
        set_density('ArPlus', 1e7)
        set_density('CH4Plus', 1e5)
        set_density('H2', 1e12)
        set_density('H', 1e11)
        set_density('C2', 5e7)
        set_density('CH', 5e4)
        set_density('CH3', 5e7)

        sol = solve_ivp(PlasmaODE(params), (0, 100), y0, method='BDF', rtol=1e-5, atol=1e-6, max_step=10.0)

        if not sol.success:
            return None, f"ODE solver failed: {sol.message}"

        y_final = sol.y[:, -1]

        def get_density(name):
            try:
                return y_final[species.index(name)]
            except ValueError:
                return 0.0

        ion_species = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus']
        n_i = sum(get_density(sp) for sp in ion_species)
        n_e = get_density('e')
        charge_imbalance = abs(n_i - n_e) / n_e * 100 if n_e > 0 else 999

        return {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2'),
            'C2H2': get_density('C2H2'),
            'charge_imbalance': charge_imbalance,
        }, None
    except Exception as e:
        return None, str(e)


# First, verify baseline works
print("="*80)
print("STEP 1: VERIFY BASELINE")
print("="*80)
print()

baseline_params = {
    'Te': start_data['Te'],
    'ne': start_data['Ne'],  # Note: capital N in JSON
    'E_field': start_data['E_field'],
    'H_drift_gain': start_data.get('H_drift_gain', 3.2e17),
    'rate_values': baseline_rates,
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
    'L_discharge': 0.45,  # cm (from documentation)
    'P': 0.4,
    'Tgas': 300,
    'mobilities': {
        'ArPlus': 3057.28,
        'CH4Plus': 6432,
        'CH3Plus': 4949.6,
        'CH5Plus': 4761.6,
        'ArHPlus': 2969.6,
        'CH2Plus': 4949.6,
        'C2H5Plus': 4949.6,
        'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6,
        'C2HPlus': 5000,
        'H3Plus': 5000,
        'CHPlus': 5000,
        'H2Plus': 5000,
        'CH3Minus': 3000,
        'HMinus': 3000
    }
}

results, error = run_simulation(baseline_params)
if results:
    print("✓ Baseline simulation SUCCESSFUL")
    print(f"  H:    {results['H']:.2e} ({results['H']/TARGETS['H']:.2f}×)")
    print(f"  CH:   {results['CH']:.2e} ({results['CH']/TARGETS['CH']:.2f}×)")
    print(f"  C2:   {results['C2']:.2e} ({results['C2']/TARGETS['C2']:.2f}×)")
    print(f"  C2H2: {results['C2H2']:.2e}")
    print(f"  Charge imbalance: {results['charge_imbalance']:.1f}%")
else:
    print(f"✗ Baseline simulation FAILED: {error}")
    print("Cannot proceed with sweep if baseline fails!")
    exit(1)

print()
print("="*80)
print("STEP 2: SMALL C2H2 BOOST SWEEP")
print("="*80)
print()
print("Strategy: Reduce stick_C2H2 and loss_C2H2 in small steps")
print("          Increase stick_CH in small steps")
print()

# Very conservative sweep: small steps
stick_C2H2_factors = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]  # 10-50% reduction
loss_C2H2_factors = [1.0, 0.9, 0.8, 0.7]              # 10-30% reduction
stick_CH_factors = [1.0, 1.1, 1.2, 1.3]               # 10-30% increase

test_count = 0
success_count = 0
best_result = None
best_score = 1e10

for sc2h2 in stick_C2H2_factors:
    for lc2h2 in loss_C2H2_factors:
        for sch in stick_CH_factors:
            test_count += 1

            # Build adjusted rates
            test_rates = baseline_rates.copy()
            if 'stick_C2H2_9_11' in test_rates:
                test_rates['stick_C2H2_9_11'] *= sc2h2
            if 'loss_C2H2_11_19' in test_rates:
                test_rates['loss_C2H2_11_19'] *= lc2h2
            if 'stick_CH_9_3' in test_rates:
                test_rates['stick_CH_9_3'] *= sch

            test_params = {
                'Te': start_data['Te'],
                'ne': start_data['Ne'],  # Note: capital N in JSON
                'E_field': start_data['E_field'],
                'H_drift_gain': start_data.get('H_drift_gain', 3.2e17),
                'rate_values': test_rates,
                'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                            'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                            'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                            'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                            'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
                'L_discharge': 0.45,  # cm (from documentation)
                'P': 0.4,
                'Tgas': 300,
                'mobilities': {
                    'ArPlus': 3057.28,
                    'CH4Plus': 6432,
                    'CH3Plus': 4949.6,
                    'CH5Plus': 4761.6,
                    'ArHPlus': 2969.6,
                    'CH2Plus': 4949.6,
                    'C2H5Plus': 4949.6,
                    'C2H4Plus': 4949.6,
                    'C2H3Plus': 4949.6,
                    'C2HPlus': 5000,
                    'H3Plus': 5000,
                    'CHPlus': 5000,
                    'H2Plus': 5000,
                    'CH3Minus': 3000,
                    'HMinus': 3000
                }
            }

            results, error = run_simulation(test_params)

            if results is None:
                continue

            success_count += 1

            # Calculate score
            score = 0.0
            h_ratio = results['H'] / TARGETS['H']
            ch_ratio = results['CH'] / TARGETS['CH']
            c2_ratio = results['C2'] / TARGETS['C2']

            score += abs(h_ratio - 1.0) ** 2
            score += abs(ch_ratio - 1.0) ** 2
            score += abs(c2_ratio - 1.0) ** 2
            score += 0.1 * (results['charge_imbalance'] / 100.0) ** 2

            if score < best_score:
                best_score = score
                best_result = {
                    'params': test_params,
                    'results': results,
                    'factors': (sc2h2, lc2h2, sch),
                    'ratios': (h_ratio, ch_ratio, c2_ratio)
                }

                print(f"Test {test_count}: stick_C2H2={sc2h2:.1f}×, loss_C2H2={lc2h2:.1f}×, stick_CH={sch:.1f}×")
                print(f"  ★ NEW BEST (score={score:.3f})")
                print(f"    H:    {results['H']:.2e} ({h_ratio:.2f}×)")
                print(f"    CH:   {results['CH']:.2e} ({ch_ratio:.2f}×)")
                print(f"    C2:   {results['C2']:.2e} ({c2_ratio:.2f}×)")
                print(f"    C2H2: {results['C2H2']:.2e}")
                print(f"    Charge: {results['charge_imbalance']:.1f}%")

                # Check if in range
                h_ok = 0.6 <= h_ratio <= 1.4
                ch_ok = 0.6 <= ch_ratio <= 1.4
                c2_ok = 0.6 <= c2_ratio <= 1.4

                if h_ok and ch_ok and c2_ok:
                    print()
                    print("  ★★★ ALL THREE IN RANGE! ★★★")
                    print()

print()
print("="*80)
print("SWEEP COMPLETE")
print("="*80)
print()
print(f"Tests run: {test_count}")
print(f"Successful: {success_count} ({success_count/test_count*100:.1f}%)")
print()

if best_result:
    r = best_result['results']
    h_ratio, ch_ratio, c2_ratio = best_result['ratios']
    sc2h2, lc2h2, sch = best_result['factors']

    print("BEST RESULT:")
    print(f"  stick_C2H2 factor: {sc2h2:.1f}×")
    print(f"  loss_C2H2 factor:  {lc2h2:.1f}×")
    print(f"  stick_CH factor:   {sch:.1f}×")
    print()
    print(f"  H:    {r['H']:.2e} ({h_ratio:.2f}×) {'✓' if 0.6 <= h_ratio <= 1.4 else '✗'}")
    print(f"  CH:   {r['CH']:.2e} ({ch_ratio:.2f}×) {'✓' if 0.6 <= ch_ratio <= 1.4 else '✗'}")
    print(f"  C2:   {r['C2']:.2e} ({c2_ratio:.2f}×) {'✓' if 0.6 <= c2_ratio <= 1.4 else '✗'}")
    print(f"  C2H2: {r['C2H2']:.2e} ({r['C2H2']/2.16e10:.2f}× baseline)")
    print(f"  Charge imbalance: {r['charge_imbalance']:.1f}%")

    # Save best result
    save_data = {
        'factors': {
            'stick_C2H2': sc2h2,
            'loss_C2H2': lc2h2,
            'stick_CH': sch
        },
        'params': best_result['params'],
        'densities': r,
        'ratios': {
            'H': h_ratio,
            'CH': ch_ratio,
            'C2': c2_ratio
        }
    }

    with open('manual_sweep_results/best_result.json', 'w') as f:
        json.dump(save_data, f, indent=2)

    print()
    print("Saved to: manual_sweep_results/best_result.json")
else:
    print("NO SUCCESSFUL RESULTS!")
