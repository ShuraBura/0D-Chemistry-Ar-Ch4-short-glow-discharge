"""
Utility to load and process experimental nH(x) profiles for 0-D model.

Since 0-D models predict spatially uniform densities, we need to reduce
the spatial profile to a single representative value.
"""

import numpy as np
from typing import Tuple, Optional


def load_nH_profile(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load experimental H density profile from file.

    Expected format: CSV with columns [x_cm, nH_cm3]

    Args:
        filepath: Path to CSV file

    Returns:
        x: Position array (cm)
        nH: H density array (cm⁻³)
    """
    data = np.loadtxt(filepath, delimiter=',', skiprows=1)
    x = data[:, 0]
    nH = data[:, 1]
    return x, nH


def calculate_spatial_average(x: np.ndarray, nH: np.ndarray) -> float:
    """
    Calculate spatial average of H density profile.

    Uses trapezoidal integration:
        <nH> = (1/L) ∫ nH(x) dx

    Args:
        x: Position array (cm)
        nH: H density array (cm⁻³)

    Returns:
        nH_avg: Spatially-averaged H density (cm⁻³)
    """
    L = x[-1] - x[0]
    nH_avg = np.trapz(nH, x) / L
    return nH_avg


def calculate_peak_density(x: np.ndarray, nH: np.ndarray) -> Tuple[float, float]:
    """
    Find peak H density and its location.

    Args:
        x: Position array (cm)
        nH: H density array (cm⁻³)

    Returns:
        x_peak: Position of peak (cm)
        nH_peak: Peak H density (cm⁻³)
    """
    peak_idx = np.argmax(nH)
    x_peak = x[peak_idx]
    nH_peak = nH[peak_idx]
    return x_peak, nH_peak


def estimate_diffusion_length(x: np.ndarray, nH: np.ndarray,
                              x0: Optional[float] = None) -> float:
    """
    Estimate diffusion length from exponential decay of profile.

    Fits: nH(x) = nH0 * exp(-x/L_diff)

    Args:
        x: Position array (cm)
        nH: H density array (cm⁻³)
        x0: Starting position for fit (default: first point)

    Returns:
        L_diff: Estimated diffusion length (cm)
    """
    if x0 is None:
        x0 = x[0]

    # Select data for fitting
    mask = x >= x0
    x_fit = x[mask]
    nH_fit = nH[mask]

    # Log-linear fit
    log_nH = np.log(nH_fit + 1e-10)  # Avoid log(0)

    # Linear fit: log(nH) = log(nH0) - x/L_diff
    coeffs = np.polyfit(x_fit, log_nH, 1)
    slope = coeffs[0]

    # L_diff = -1/slope
    L_diff = -1.0 / slope if slope != 0 else np.inf

    return abs(L_diff)


def create_example_profile(savepath: str = 'nH_profile_example.csv'):
    """
    Create example nH(x) profile for testing.

    Simulates a Gaussian-like profile with exponential decay.
    """
    x = np.linspace(0, 0.1, 50)  # 0-1 mm in cm

    # Gaussian peak with decay
    x_peak = 0.03  # Peak at 0.3 mm
    width = 0.02   # Width ~ 0.2 mm
    nH_peak = 1e16  # Peak density

    nH = nH_peak * np.exp(-(x - x_peak)**2 / (2 * width**2))

    # Add some noise
    nH += np.random.normal(0, nH_peak * 0.05, len(x))
    nH = np.maximum(nH, 1e14)  # Floor

    # Save
    header = "x_cm,nH_cm3"
    data = np.column_stack((x, nH))
    np.savetxt(savepath, data, delimiter=',', header=header, comments='')

    print(f"Example profile saved to {savepath}")
    print(f"Spatial average: {calculate_spatial_average(x, nH):.2e} cm⁻³")
    print(f"Peak: {np.max(nH):.2e} cm⁻³ at x = {x[np.argmax(nH)]:.4f} cm")


def use_profile_in_simulation(profile_file: str, method: str = 'average'):
    """
    Process experimental profile for use in 0-D simulation.

    Args:
        profile_file: Path to nH(x) CSV file
        method: How to reduce profile to single value
                'average': Spatial average (recommended)
                'peak': Peak value
                'rms': RMS average

    Returns:
        nH_fixed: H density value to use in 0-D model (cm⁻³)
    """
    x, nH = load_nH_profile(profile_file)

    if method == 'average':
        nH_fixed = calculate_spatial_average(x, nH)
        print(f"Using spatial average: {nH_fixed:.2e} cm⁻³")

    elif method == 'peak':
        x_peak, nH_fixed = calculate_peak_density(x, nH)
        print(f"Using peak value: {nH_fixed:.2e} cm⁻³ at x = {x_peak:.4f} cm")

    elif method == 'rms':
        nH_fixed = np.sqrt(np.trapz(nH**2, x) / (x[-1] - x[0]))
        print(f"Using RMS average: {nH_fixed:.2e} cm⁻³")

    else:
        raise ValueError(f"Unknown method: {method}")

    # Also calculate diffusion length for reference
    L_diff = estimate_diffusion_length(x, nH)
    print(f"Estimated L_diff from profile: {L_diff:.4f} cm")

    return nH_fixed, L_diff


if __name__ == '__main__':
    # Example usage
    print("Creating example nH(x) profile...")
    create_example_profile()

    print("\nProcessing profile for 0-D model...")
    nH_fixed, L_diff = use_profile_in_simulation('nH_profile_example.csv',
                                                  method='average')

    print("\nRecommendation:")
    print(f"- Set H density in model: set_density('H', {nH_fixed:.2e})")
    print(f"- Use L_diff = {L_diff:.4f} cm in define_rates_tunable.py")
    print(f"- If using fixed H, set dydt[H_idx] = 0 in odefun.py")
