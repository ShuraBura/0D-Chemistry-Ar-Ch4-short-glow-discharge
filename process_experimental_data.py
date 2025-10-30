#!/usr/bin/env python3
"""
Process actual experimental nH data from combined_density_matrix.csv

File format (as described):
- Column 1: Distance from cathode (assuming meters or mm)
- Column 2: Measured nH (m⁻³)

This script will:
1. Load the experimental data
2. Convert units to cm and cm⁻³
3. Calculate spatial average for 0-1 mm region (CG)
4. Estimate diffusion length from profile decay
5. Provide values to use in 0-D model
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def load_experimental_nH(filepath='combined_density_matrix.csv'):
    """
    Load experimental H density profile.

    Args:
        filepath: Path to CSV file

    Returns:
        x_cm: Position from cathode (cm)
        nH_cm3: H density (cm⁻³)
    """
    print(f"Loading data from {filepath}...")

    # Load CSV (assuming no header, space or comma separated)
    try:
        data = np.loadtxt(filepath, delimiter=',')
    except:
        # Try space-separated
        data = np.loadtxt(filepath)

    # Extract columns
    x_raw = data[:, 0]  # Position (need to determine units)
    nH_raw = data[:, 1]  # Density in m⁻³

    # Determine units of x by examining range
    if x_raw.max() < 0.1:
        # Already in meters or cm (< 10 cm makes sense)
        if x_raw.max() < 0.01:
            # Likely meters (< 1 cm)
            x_cm = x_raw * 100  # m → cm
            print(f"Position interpreted as meters, converted to cm")
        else:
            # Likely already cm
            x_cm = x_raw
            print(f"Position interpreted as cm")
    else:
        # Likely mm (0-10 mm range is typical)
        x_cm = x_raw / 10  # mm → cm
        print(f"Position interpreted as mm, converted to cm")

    # Convert density m⁻³ → cm⁻³
    nH_cm3 = nH_raw / 1e6

    print(f"Loaded {len(x_cm)} data points")
    print(f"Position range: {x_cm.min():.4f} - {x_cm.max():.4f} cm")
    print(f"Density range: {nH_cm3.min():.2e} - {nH_cm3.max():.2e} cm⁻³")

    return x_cm, nH_cm3


def extract_CG_region(x_cm, nH_cm3, x_min=0.0, x_max=0.1):
    """
    Extract data for Cathode Glow region (0-1 mm = 0-0.1 cm).

    Args:
        x_cm: Position array (cm)
        nH_cm3: Density array (cm⁻³)
        x_min: Start of CG region (cm)
        x_max: End of CG region (cm)

    Returns:
        x_cg: Position in CG region (cm)
        nH_cg: Density in CG region (cm⁻³)
    """
    mask = (x_cm >= x_min) & (x_cm <= x_max)
    x_cg = x_cm[mask]
    nH_cg = nH_cm3[mask]

    print(f"\nCG region (x = {x_min:.3f} - {x_max:.3f} cm):")
    print(f"  {len(x_cg)} data points")
    print(f"  Density range: {nH_cg.min():.2e} - {nH_cg.max():.2e} cm⁻³")

    return x_cg, nH_cg


def calculate_spatial_average(x, nH):
    """
    Calculate spatial average: <nH> = (1/L) ∫ nH(x) dx

    This is what your TALIF measurement represents.
    """
    if len(x) < 2:
        return nH[0] if len(nH) > 0 else 0

    L = x[-1] - x[0]
    nH_avg = np.trapezoid(nH, x) / L if L > 0 else np.mean(nH)

    return nH_avg


def fit_exponential_decay(x, nH):
    """
    Fit exponential decay to extract diffusion length.

    Model: nH(x) = A * exp(-x / L_diff) + B

    Returns:
        L_diff: Diffusion length (cm)
        fit_params: All fit parameters
    """
    # Define exponential decay model
    def exp_decay(x, A, L_diff, B):
        return A * np.exp(-x / L_diff) + B

    try:
        # Initial guess
        A0 = nH[0] - nH[-1]
        L_diff0 = (x[-1] - x[0]) / 2
        B0 = nH[-1]

        # Fit
        popt, pcov = curve_fit(
            exp_decay,
            x,
            nH,
            p0=[A0, L_diff0, B0],
            bounds=([0, 0.001, 0], [np.inf, 1.0, np.inf])
        )

        A_fit, L_diff_fit, B_fit = popt

        # Calculate R²
        nH_fit = exp_decay(x, *popt)
        ss_res = np.sum((nH - nH_fit)**2)
        ss_tot = np.sum((nH - np.mean(nH))**2)
        r_squared = 1 - (ss_res / ss_tot)

        print(f"\nExponential decay fit:")
        print(f"  nH(x) = {A_fit:.2e} * exp(-x/{L_diff_fit:.4f}) + {B_fit:.2e}")
        print(f"  L_diff = {L_diff_fit:.4f} cm")
        print(f"  R² = {r_squared:.4f}")

        return L_diff_fit, popt, r_squared

    except Exception as e:
        print(f"\nExponential fit failed: {e}")
        print("Using simple estimate instead...")

        # Simple estimate: L_diff ~ x where n drops to 1/e
        n_target = nH[0] / np.e
        idx = np.argmin(np.abs(nH - n_target))
        L_diff_simple = x[idx] - x[0]

        return L_diff_simple, None, 0


def plot_profile(x_cm, nH_cm3, x_cg, nH_cg, nH_avg, L_diff,
                 savefig='experimental_nH_profile.png'):
    """
    Plot the experimental profile with analysis.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left plot: Full profile
    ax1.plot(x_cm * 10, nH_cm3, 'o-', label='Full profile', markersize=4)
    ax1.axvspan(0, 1.0, alpha=0.2, color='orange', label='CG region (0-1 mm)')
    ax1.axhline(nH_avg, color='red', linestyle='--',
                label=f'Spatial avg: {nH_avg:.2e} cm⁻³')
    ax1.set_xlabel('Distance from cathode (mm)')
    ax1.set_ylabel('H density (cm⁻³)')
    ax1.set_title('Experimental H Density Profile')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Right plot: CG region with fit
    ax2.plot(x_cg * 10, nH_cg, 'o', label='CG data', markersize=6)
    ax2.axhline(nH_avg, color='red', linestyle='--',
                label=f'<nH> = {nH_avg:.2e}')

    # Plot exponential fit if available
    if L_diff > 0:
        x_fit = np.linspace(x_cg[0], x_cg[-1], 100)
        nH_fit = nH_cg[0] * np.exp(-(x_fit - x_cg[0]) / L_diff)
        ax2.plot(x_fit * 10, nH_fit, 'g--',
                 label=f'Exp fit (L_diff={L_diff:.3f} cm)')

    ax2.set_xlabel('Distance from cathode (mm)')
    ax2.set_ylabel('H density (cm⁻³)')
    ax2.set_title('CG Region (0-1 mm)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(savefig, dpi=150)
    print(f"\nPlot saved to {savefig}")

    return fig


def main():
    """
    Main processing workflow.
    """
    print("="*70)
    print("EXPERIMENTAL H PROFILE ANALYSIS")
    print("="*70)

    # Check if file exists
    import os
    if not os.path.exists('combined_density_matrix.csv'):
        print("\nERROR: combined_density_matrix.csv not found!")
        print("\nPlease add the file to the repository with format:")
        print("  Column 1: Distance from cathode (m or mm)")
        print("  Column 2: Measured nH (m⁻³)")
        print("\nExample:")
        print("  0.0000, 8.58e21")
        print("  0.0002, 7.45e21")
        print("  0.0004, 6.82e21")
        print("  ...")
        return

    # Load data
    x_cm, nH_cm3 = load_experimental_nH('combined_density_matrix.csv')

    # Extract CG region (0-1 mm)
    x_cg, nH_cg = extract_CG_region(x_cm, nH_cm3, x_min=0.0, x_max=0.1)

    # Calculate spatial average
    nH_avg = calculate_spatial_average(x_cg, nH_cg)
    print(f"\n{'='*70}")
    print(f"SPATIAL AVERAGE (CG region): {nH_avg:.3e} cm⁻³")
    print(f"{'='*70}")

    # Find peak
    peak_idx = np.argmax(nH_cg)
    nH_peak = nH_cg[peak_idx]
    x_peak = x_cg[peak_idx]
    print(f"\nPeak density: {nH_peak:.3e} cm⁻³ at x = {x_peak*10:.2f} mm")
    print(f"Peak/Average ratio: {nH_peak/nH_avg:.2f}")

    # Estimate diffusion length
    L_diff, fit_params, r_squared = fit_exponential_decay(x_cg, nH_cg)

    # Plot
    try:
        fig = plot_profile(x_cm, nH_cm3, x_cg, nH_cg, nH_avg, L_diff)
    except Exception as e:
        print(f"\nPlotting failed: {e}")

    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS FOR 0-D MODEL")
    print("="*70)
    print(f"\n1. FIXED H DENSITY (Option B):")
    print(f"   Use in your simulation:")
    print(f"   set_density('H', {nH_avg:.3e})  # Spatial average")
    print(f"   # In odefun.py: dydt[H_idx] = 0")

    print(f"\n2. DIFFUSION LENGTH:")
    print(f"   Update in define_rates_tunable.py:")
    print(f"   params['L_diff'] = {L_diff:.4f}  # cm (from profile fit)")

    print(f"\n3. COMPARISON:")
    print(f"   Target H (CG):   {nH_avg:.3e} cm⁻³ (from experiment)")
    print(f"   Target CH (CG):  4.6e8 cm⁻³ (from EXPERIMENTAL_TARGETS.md)")
    print(f"   Target C₂ (CG):  1.44e11 cm⁻³")

    print(f"\n4. VALIDATION:")
    print(f"   Run: python3 run_with_fixed_H.py")
    print(f"   Check if CH and C₂ match targets within factor 2-5")

    # Save processed data
    output_file = 'nH_profile_processed.csv'
    header = "x_cm,nH_cm3,region"
    data_out = np.column_stack((
        x_cm,
        nH_cm3,
        ['CG' if 0 <= x <= 0.1 else 'bulk' for x in x_cm]
    ))
    np.savetxt(output_file, data_out, delimiter=',', header=header,
               fmt=['%.6f', '%.6e', '%s'], comments='')
    print(f"\nProcessed data saved to {output_file}")

    # Save summary
    summary_file = 'nH_profile_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("Experimental H Profile Summary\n")
        f.write("="*50 + "\n\n")
        f.write(f"CG Region (0-1 mm):\n")
        f.write(f"  Spatial average: {nH_avg:.3e} cm⁻³\n")
        f.write(f"  Peak density:    {nH_peak:.3e} cm⁻³\n")
        f.write(f"  Peak location:   {x_peak*10:.2f} mm\n")
        f.write(f"  Peak/Avg ratio:  {nH_peak/nH_avg:.2f}\n\n")
        f.write(f"Diffusion length: {L_diff:.4f} cm\n")
        f.write(f"Fit quality (R²): {r_squared:.4f}\n\n")
        f.write("Use in 0-D model:\n")
        f.write(f"  nH_fixed = {nH_avg:.3e}  # cm⁻³\n")
        f.write(f"  L_diff = {L_diff:.4f}     # cm\n")

    print(f"Summary saved to {summary_file}")
    print("="*70)


if __name__ == '__main__':
    main()
