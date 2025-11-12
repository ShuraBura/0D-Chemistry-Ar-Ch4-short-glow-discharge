"""
Run optimizer WITHOUT H_drift_gain to match all three targets

Now testing pure chemistry - no artificial H source
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=300):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

targets = {'H': 2.52e14, 'CH': 1.0e9, 'C2': 5.6e11}

# Simplified: optimize ne, C2+H suppression, and dust
PARAMS = [
    ('ne', 1e9, 1e7, 1e10),
    ('C2_H_mult', 1.0, 0.001, 1.0),
    ('dust_dens', 0.0, 0.0, 1e9),
    ('CH_H_mult', 1.0, 0.1, 10.0),
]

def run_sim(params_array):
    ne_val, c2h_mult, dust_dens, ch_h_mult = params_array

    P = 500.0
    n_total = pressure_to_density(P)

    params = {
        'P': P, 'Te': baseline['Te'], 'ne': ne_val,
        'E_field': baseline['E_field'], 'L_discharge': 0.45,
        'mobilities': {'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6, 'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6, 'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000, 'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000},
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
    }

    if dust_dens > 1e6:
        params['enable_dust_loss'] = True
        params['dust_density'] = dust_dens
        params['dust_radius'] = 50e-7
        params['dust_sticking'] = 0.5

    k = define_rates(params)
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    rate_mults = {
        'e_CH4_CH3_HMinus_cm3_8_1': 10.0, 'ArStar_CH4_CH3_H_cm3_3_1': 1.4,
        'e_CH4Plus_CH3_H_cm3_6_4': 1.5, 'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01, 'loss_C2H2_11_19': 0.01,
        'C2_H_CH_C_cm3_7_6': c2h_mult,
        'CH_H_C_H2_cm3_7_3': ch_h_mult,
    }

    for rate_name, mult in rate_mults.items():
        if rate_name in k:
            k[rate_name] *= mult

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.85
    y0[species.index('CH4')] = n_total * 0.15
    y0[species.index('e')] = ne_val

    try:
        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0)
        if not sol.success:
            return None

        y_final = sol.y[:, -1]
        H = y_final[species.index('H')]
        CH = y_final[species.index('CH')]
        C2 = y_final[species.index('C2')]

        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        Ni = sum(y_final[species.index(ion)] for ion in ions if ion in species)
        Ne = y_final[species.index('e')]

        if Ni / Ne > 10:
            return None

        return {'H': H, 'CH': CH, 'C2': C2, 'Ni_Ne': Ni/Ne}
    except:
        return None

def objective(params_array):
    result = run_sim(params_array)
    if result is None:
        return 1e10

    h_err = abs(result['H'] / targets['H'] - 1.0) * 100
    ch_err = abs(result['CH'] / targets['CH'] - 1.0) * 100
    c2_err = abs(result['C2'] / targets['C2'] - 1.0) * 100

    return np.sqrt((h_err**2 + ch_err**2 + c2_err**2) / 3)

if __name__ == '__main__':
    print('='*80)
    print('OPTIMIZE WITHOUT H_drift_gain (pure chemistry)')
    print('='*80)
    print(f'\nTargets: H={targets["H"]:.2e}, CH={targets["CH"]:.2e}, C2={targets["C2"]:.2e}\n')

    bounds = [(p[2], p[3]) for p in PARAMS]
    x0 = [p[1] for p in PARAMS]

    print('Baseline test:')
    baseline_err = objective(np.array(x0))
    result = run_sim(np.array(x0))
    if result:
        print(f'  H={result["H"]/targets["H"]*100:.1f}%, CH={result["CH"]/targets["CH"]*100:.0f}%, C2={result["C2"]/targets["C2"]*100:.1f}%, RMS={baseline_err:.1f}%\n')

    print('Running optimizer (this will take a few minutes)...\n')

    res = differential_evolution(
        objective, bounds, strategy='best1bin', maxiter=30, popsize=10,
        tol=0.01, seed=42, polish=True, workers=1,
        callback=lambda xk, convergence: print(f'  Iteration, best RMS: {convergence:.1f}%'),
    )

    print(f'\n{"="*80}')
    print('OPTIMIZATION COMPLETE')
    print(f'{"="*80}\n')

    if res.fun < baseline_err:
        print(f'✓ Improvement: {baseline_err:.1f}% → {res.fun:.1f}%\n')

        print('Optimal parameters:')
        for i, (name, _, _, _) in enumerate(PARAMS):
            print(f'  {name}: {res.x[i]:.3e}')

        final = run_sim(res.x)
        if final:
            print(f'\nFinal results:')
            print(f'  H:  {final["H"]:.2e} cm⁻³ ({final["H"]/targets["H"]*100:6.2f}%)')
            print(f'  CH: {final["CH"]:.2e} cm⁻³ ({final["CH"]/targets["CH"]*100:6.2f}%)')
            print(f'  C2: {final["C2"]:.2e} cm⁻³ ({final["C2"]/targets["C2"]*100:6.2f}%)')
            print(f'  Ni/Ne: {final["Ni_Ne"]:.2f}')
            print(f'  RMS error: {res.fun:.1f}%')

            # Save
            opt_data = {
                'ne': float(res.x[0]),
                'C2_H_CH_C_multiplier': float(res.x[1]),
                'dust_density': float(res.x[2]),
                'CH_H_C_H2_multiplier': float(res.x[3]),
                'results': {
                    'H': float(final['H']), 'CH': float(final['CH']), 'C2': float(final['C2']),
                    'H_pct': float(final['H']/targets['H']*100),
                    'CH_pct': float(final['CH']/targets['CH']*100),
                    'C2_pct': float(final['C2']/targets['C2']*100),
                    'Ni_Ne': float(final['Ni_Ne']),
                    'rms_error': float(res.fun),
                }
            }
            with open('optimized_rates_no_H_drift.json', 'w') as f:
                json.dump(opt_data, f, indent=2)
            print('\nSaved to optimized_rates_no_H_drift.json')
    else:
        print(f'✗ No improvement: RMS = {res.fun:.1f}%')
