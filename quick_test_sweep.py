#!/usr/bin/env python3
"""
Quick test of parameter sweep with reduced grid
"""

import sys
sys.path.insert(0, '/home/user/0D-Chemistry-Ar-Ch4-short-glow-discharge')

# Monkey-patch the sweep parameters for quick test
import parameter_sweep_cg as psweep

# Override main to use smaller grid
original_main = psweep.main

def quick_main():
    # Smaller grid for testing
    print("Running QUICK TEST with reduced parameter grid...")
    print("(For full sweep, run: python3 parameter_sweep_cg.py)\n")

    # Temporarily replace parameter ranges
    old_code = psweep.__dict__.copy()

    # Run with reduced grid
    psweep.Te_values = [1.0, 2.0, 3.0]  # Just 3 values
    psweep.ne_values = [1e8, 1e9]       # Just 2 values

    # Modify main to use our ranges
    exec(open('parameter_sweep_cg.py').read().replace(
        'Te_values = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0]',
        'Te_values = [1.0, 2.0, 3.0]'
    ).replace(
        'ne_values = [1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10]',
        'ne_values = [1e8, 1e9]'
    ), psweep.__dict__)

if __name__ == '__main__':
    quick_main()
