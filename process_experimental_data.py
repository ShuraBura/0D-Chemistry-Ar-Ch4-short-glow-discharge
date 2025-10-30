#!/usr/bin/env python3
"""
Process experimental density data from combined_density_matrix.csv
Extracts CG region (0-1 mm), calculates spatial averages, and estimates diffusion parameters.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate, optimize

def load_experimental_data(filename='combined_density_matrix.csv'):
    """Load and parse experimental density data."""
    # Read CSV, skipping first row (header), using second row for column names
    df = pd.read_csv(filename, skiprows=1, header=0)

    # Clean column names
    df.columns = ['distance_mm', 'nH_m3', 'nH_err', 'nCH_m3', 'nCH_err', 'nC2_m3', 'nC2_err']

    # Convert units from m^-3 to cm^-3
    df['nH_cm3'] = df['nH_m3'] / 1e6
    df['nCH_cm3'] = df['nCH_m3'] / 1e6
    df['nC2_cm3'] = df['nC2_m3'] / 1e6
    df['nH_err_cm3'] = df['nH_err'] / 1e6
    df['nCH_err_cm3'] = df['nCH_err'] / 1e6
    df['nC2_err_cm3'] = df['nC2_err'] / 1e6

    print(f"✓ Loaded {len(df)} data points from {filename}")
    print(f"  Distance range: {df['distance_mm'].min():.3f} - {df['distance_mm'].max():.3f} mm")

    return df

def extract_cg_region(df, z_min=0.0, z_max=1.0):
    """Extract data in CG (Cathode Glow) region."""
    mask = (df['distance_mm'] >= z_min) & (df['distance_mm'] <= z_max)
    df_cg = df[mask].copy()

    print(f"\n✓ Extracted CG region ({z_min}-{z_max} mm): {len(df_cg)} points")

    return df_cg

def calculate_spatial_average(df_cg, species='H'):
    """
    Calculate spatial average density using trapezoidal integration.

    <n> = (1/L) ∫ n(z) dz  from z=0 to z=L
    """
    z = df_cg['distance_mm'].values / 10.0  # Convert mm to cm
    n = df_cg[f'n{species}_cm3'].values
    n_err = df_cg[f'n{species}_err_cm3'].values

    # Remove any NaN or zero values for valid integration
    valid = ~np.isnan(n) & (n > 0)
    z_valid = z[valid]
    n_valid = n[valid]
    n_err_valid = n_err[valid]

    if len(z_valid) < 2:
        print(f"  ⚠️  {species}: Insufficient valid data points")
        return np.nan, np.nan, np.nan

    # Integrate using trapezoidal rule
    L = z_valid[-1] - z_valid[0]  # cm
    integral = integrate.trapezoid(n_valid, z_valid)
    n_avg = integral / L

    # Estimate uncertainty (simplified - assumes uncorrelated errors)
    integral_err = integrate.trapezoid(n_err_valid, z_valid)
    n_avg_err = integral_err / L

    # Also report peak value
    n_peak = np.max(n_valid)
    z_peak = z_valid[np.argmax(n_valid)]

    print(f"  {species}: <n> = {n_avg:.2e} ± {n_avg_err:.2e} cm⁻³")
    print(f"      Peak: {n_peak:.2e} cm⁻³ at z={z_peak*10:.3f} mm")
    print(f"      Peak/Avg ratio: {n_peak/n_avg:.2f}")

    return n_avg, n_avg_err, n_peak

def fit_exponential_decay(df_cg, species='H', plot=False):
    """
    Fit exponential decay model: n(z) = n0 * exp(-z/L_diff)
    to estimate diffusion length L_diff
    """
    z = df_cg['distance_mm'].values / 10.0  # mm to cm
    n = df_cg[f'n{species}_cm3'].values

    # Remove invalid points
    valid = ~np.isnan(n) & (n > 0)
    z_valid = z[valid]
    n_valid = n[valid]

    if len(z_valid) < 5:
        print(f"  ⚠️  {species}: Insufficient data for exponential fit")
        return None, None

    # Find peak position and fit from peak onwards (decay region)
    idx_peak = np.argmax(n_valid)
    z_decay = z_valid[idx_peak:]
    n_decay = n_valid[idx_peak:]

    if len(z_decay) < 3:
        print(f"  ⚠️  {species}: Insufficient decay data")
        return None, None

    # Fit: ln(n) = ln(n0) - z/L_diff
    try:
        # Shift z to start from peak
        z_shifted = z_decay - z_decay[0]
        popt, pcov = optimize.curve_fit(
            lambda z, n0, L: n0 * np.exp(-z/L),
            z_shifted, n_decay,
            p0=[n_decay[0], 0.1],  # Initial guess: n0=peak, L_diff=0.1 cm
            maxfev=5000
        )
        n0_fit, L_diff = popt
        L_diff_err = np.sqrt(pcov[1, 1])

        print(f"  {species}: L_diff = {L_diff:.4f} ± {L_diff_err:.4f} cm")

        if plot:
            plt.figure(figsize=(10, 6))
            plt.plot(z_valid*10, n_valid, 'o', label=f'{species} data', markersize=4)
            z_fit = np.linspace(z_decay[0], z_decay[-1], 100)
            n_fit = n0_fit * np.exp(-(z_fit - z_decay[0])/L_diff)
            plt.plot(z_fit*10, n_fit, '-', label=f'Exp fit: L={L_diff:.3f} cm')
            plt.axvline(1.0, color='r', linestyle='--', label='CG boundary (1 mm)')
            plt.xlabel('Distance from cathode (mm)')
            plt.ylabel(f'n_{species} (cm⁻³)')
            plt.yscale('log')
            plt.legend()
            plt.title(f'{species} Density Profile with Exponential Fit')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(f'n{species}_profile_fit.png', dpi=150)
            print(f"  Saved plot: n{species}_profile_fit.png")
            plt.close()

        return L_diff, L_diff_err

    except Exception as e:
        print(f"  ⚠️  {species}: Fit failed: {e}")
        return None, None

def plot_profiles(df_cg, save=True):
    """Plot all three species profiles."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    species_list = ['H', 'CH', 'C2']

    for ax, species in zip(axes, species_list):
        z = df_cg['distance_mm'].values
        n = df_cg[f'n{species}_cm3'].values
        n_err = df_cg[f'n{species}_err_cm3'].values

        # Remove invalid points
        valid = ~np.isnan(n) & (n > 0)
        z_valid = z[valid]
        n_valid = n[valid]
        n_err_valid = n_err[valid]

        ax.errorbar(z_valid, n_valid, yerr=n_err_valid, fmt='o-',
                   markersize=4, capsize=3, label=f'n_{species}')
        ax.axvline(1.0, color='r', linestyle='--', alpha=0.5, label='CG edge')
        ax.set_xlabel('Distance from cathode (mm)')
        ax.set_ylabel(f'n_{species} (cm⁻³)')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{species} Density Profile')

    plt.tight_layout()
    if save:
        plt.savefig('all_species_profiles.png', dpi=150)
        print("\n✓ Saved plot: all_species_profiles.png")
    plt.close()

def main():
    """Main processing workflow."""
    print("="*70)
    print("EXPERIMENTAL DATA PROCESSING - CG REGION")
    print("="*70)

    # Load data
    df = load_experimental_data('combined_density_matrix.csv')

    # Extract CG region (0-1 mm)
    df_cg = extract_cg_region(df, z_min=0.0, z_max=1.0)

    # Calculate spatial averages
    print("\n" + "="*70)
    print("SPATIAL AVERAGES (CG Region: 0-1 mm)")
    print("="*70)

    results = {}
    for species in ['H', 'CH', 'C2']:
        n_avg, n_err, n_peak = calculate_spatial_average(df_cg, species)
        results[species] = {'avg': n_avg, 'err': n_err, 'peak': n_peak}

    # Estimate diffusion lengths
    print("\n" + "="*70)
    print("DIFFUSION LENGTH ESTIMATES")
    print("="*70)

    L_diff = {}
    for species in ['H', 'CH', 'C2']:
        L, L_err = fit_exponential_decay(df_cg, species, plot=True)
        L_diff[species] = L

    # Plot all profiles
    plot_profiles(df_cg)

    # Summary and recommendations
    print("\n" + "="*70)
    print("SUMMARY & RECOMMENDATIONS")
    print("="*70)

    print("\n1. SPATIAL AVERAGES (for 0-D model validation):")
    print(f"   H:  {results['H']['avg']:.2e} cm⁻³")
    print(f"   CH: {results['CH']['avg']:.2e} cm⁻³")
    print(f"   C₂: {results['C2']['avg']:.2e} cm⁻³")

    print("\n2. COMPARISON WITH TARGETS:")
    targets = {'H': 8.58e15, 'CH': 4.6e8, 'C2': 1.44e11}
    for species in ['H', 'CH', 'C2']:
        ratio = results[species]['avg'] / targets[species]
        print(f"   {species}: Measured/Target = {ratio:.2f}")

    print("\n3. DIFFUSION LENGTHS (for model input):")
    for species in ['H', 'CH', 'C2']:
        if L_diff[species] is not None:
            print(f"   L_diff_{species} = {L_diff[species]:.4f} cm")
        else:
            print(f"   L_diff_{species} = Not determined (use default 0.1 cm)")

    print("\n4. NEXT STEPS:")
    print("   a) Use spatial averages above for 0-D model comparison")
    print("   b) Implement Option A (dynamic H) OR Option B (fixed H)")
    print("   c) Run parameter sweep with updated targets")
    print("   d) Validate model predictions against measured densities")

    print("\n" + "="*70)
    print("✓ Processing complete!")
    print("="*70)

    # Save results to file
    with open('experimental_averages.txt', 'w') as f:
        f.write("# Experimental Spatial Averages (CG Region: 0-1 mm)\n")
        f.write("# Generated by process_experimental_data.py\n\n")
        f.write(f"nH_avg  = {results['H']['avg']:.6e}  # cm^-3\n")
        f.write(f"nCH_avg = {results['CH']['avg']:.6e}  # cm^-3\n")
        f.write(f"nC2_avg = {results['C2']['avg']:.6e}  # cm^-3\n")
        f.write(f"\n# Diffusion lengths\n")
        for species in ['H', 'CH', 'C2']:
            if L_diff[species] is not None:
                f.write(f"L_diff_{species} = {L_diff[species]:.6f}  # cm\n")

    print("\n✓ Saved: experimental_averages.txt")

    return results, L_diff

if __name__ == '__main__':
    results, L_diff = main()
