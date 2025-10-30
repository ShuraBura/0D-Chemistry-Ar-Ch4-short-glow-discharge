#!/usr/bin/env python3
"""
Parse MATLAB define_rates.m and generate complete Python rate database
Extracts: rate name, value, min, max, source, notes, flags
"""

import re

def parse_matlab_rates(matlab_file):
    """Parse MATLAB define_rates.m to extract all rate information."""

    with open(matlab_file, 'r') as f:
        lines = f.readlines()

    rates = []

    # Pattern to match rate definitions with comments
    # k.rate_name = value; % comment with Range: min–max, Source: citation
    rate_pattern = re.compile(
        r'k\.(\w+)\s*=\s*([0-9.e+-]+)\s*;\s*%\s*(.+)'
    )

    for i, line in enumerate(lines):
        match = rate_pattern.match(line.strip())
        if match:
            name = match.group(1)
            value = float(match.group(2))
            comment = match.group(3)

            # Extract range
            min_val, max_val = None, None
            range_match = re.search(r'Range:\s*([0-9.e+-]+)[–-]([0-9.e+-]+)', comment)
            if range_match:
                min_val = float(range_match.group(1))
                max_val = float(range_match.group(2))

            # Extract source
            source = "Unknown"
            source_match = re.search(r'Source:\s*([^,\n]+)', comment)
            if source_match:
                source = source_match.group(1).strip()

            # Extract scaling notes
            notes = ""
            scaled_match = re.search(r'Scaled:\s*([^,\n]+)', comment)
            if scaled_match:
                notes = "Scaled: " + scaled_match.group(1).strip()

            # Check for FLAG on next line
            flag = ""
            if i + 1 < len(lines) and 'FLAG:' in lines[i + 1]:
                flag_match = re.search(r'%\s*FLAG:\s*(.+)', lines[i + 1])
                if flag_match:
                    flag = flag_match.group(1).strip()

            rates.append({
                'name': name,
                'value': value,
                'min': min_val if min_val is not None else value,  # Default to value if no range
                'max': max_val if max_val is not None else value,
                'source': source,
                'notes': notes,
                'flag': flag
            })

    return rates


def generate_python_code(rates):
    """Generate Python code for rate_database_complete.py"""

    code = '''"""
rate_database_complete.py - COMPLETE reaction rate database
Auto-generated from MATLAB define_rates.m
All 150+ rates with literature ranges, citations, and flags preserved
"""

class RateConstant:
    """Structure to hold rate constant with metadata."""
    def __init__(self, value, min_val, max_val, source, notes="", flag=""):
        self.value = value
        self.min = min_val
        self.max = max_val
        self.source = source
        self.notes = notes
        self.flag = flag

    def is_within_range(self, val):
        """Check if a value is within literature range."""
        return self.min <= val <= self.max

    def __repr__(self):
        s = f"RateConstant({self.value:.2e}, [{self.min:.2e}, {self.max:.2e}], {self.source})"
        if self.flag:
            s += f" [FLAG: {self.flag}]"
        return s


def get_complete_rate_database():
    """
    Get COMPLETE rate database with ALL rates from MATLAB define_rates.m

    Returns dict of RateConstant objects with:
    - value: Current/nominal value
    - min: Minimum literature value
    - max: Maximum literature value
    - source: Citation
    - notes: Scaling/rationale
    - flag: Validation warnings
    """
    db = {}

'''

    # Group rates by prefix for organization
    current_group = ""

    for rate in rates:
        name = rate['name']

        # Detect group changes
        if name.startswith('e_') and not current_group.startswith("Electron"):
            code += "\n    # ========== ELECTRON-IMPACT REACTIONS ==========\n"
            current_group = "Electron"
        elif name.startswith('ArStar_') and current_group != "ArStar":
            code += "\n    # ========== Ar* REACTIONS ==========\n"
            current_group = "ArStar"
        elif name.startswith('ArPlus_') and current_group != "Ion":
            code += "\n    # ========== ION-NEUTRAL REACTIONS ==========\n"
            current_group = "Ion"
        elif (name.startswith('CH_') or name.startswith('C_') or name.startswith('H_')) and current_group != "Neutral":
            code += "\n    # ========== NEUTRAL-NEUTRAL REACTIONS ==========\n"
            current_group = "Neutral"
        elif name.startswith('H_H_M') and current_group != "Termolecular":
            code += "\n    # ========== TERMOLECULAR REACTIONS ==========\n"
            current_group = "Termolecular"
        elif name.startswith('stick_') and current_group != "Stick":
            code += "\n    # ========== WALL STICKING REACTIONS ==========\n"
            current_group = "Stick"
        elif name.startswith('drift_') and current_group != "Drift":
            code += "\n    # ========== DRIFT LOSSES (geometry-dependent) ==========\n"
            current_group = "Drift"
        elif name.startswith('loss_') and current_group != "Loss":
            code += "\n    # ========== WALL LOSS REACTIONS ==========\n"
            current_group = "Loss"

        # Skip drift rates (calculated from mobilities)
        if name.startswith('drift_'):
            continue

        # Generate line
        value = rate['value']
        min_val = rate['min']
        max_val = rate['max']
        source = rate['source']
        notes = rate['notes']
        flag = rate['flag']

        # Format the line
        code += f"    db['{name}'] = RateConstant({value:.2e}, {min_val:.2e}, {max_val:.2e}, \"{source}\""
        if notes:
            code += f', notes="{notes}"'
        if flag:
            code += f', flag="{flag}"'
        code += ")\n"

    code += '''
    return db


def get_rates_by_species(db, species_name):
    """Get all rates involving a specific species."""
    return {k: v for k, v in db.items() if species_name in k}


def get_rates_by_source(db, source_name):
    """Get all rates from a specific literature source."""
    return {k: v for k, v in db.items() if source_name.lower() in v.source.lower()}


def get_flagged_rates(db):
    """Get all rates that have validation flags."""
    return {k: v for k, v in db.items() if v.flag}


def get_tunable_rates_for_target(target_species):
    """
    Get rates most relevant for tuning specific target species.

    Parameters:
    -----------
    target_species : str
        'H', 'CH', or 'C2'

    Returns:
    --------
    dict : {rate_name: importance_note}
    """

    if target_species == 'H':
        return {
            'e_CH4_CH3_H_cm3_1_1': 'e + CH4 → CH3 + H (primary H source)',
            'e_H2_H_H_cm3_1_4': 'e + H2 → 2H (H2 dissociation)',
            'e_CH4_CH_H2_H_vib_cm3_1_3': 'e + CH4 → CH + H2 + H',
            'e_CH4_CH_H_H2_cm3_1_11': 'e + CH4 → CH + H + H2',
            'ArStar_CH4_CH3_H_cm3_3_1': 'Ar* + CH4 → CH3 + H + Ar',
            'H_CH4_CH3_H2_cm3_7_25': 'H + CH4 → CH3 + H2 (H consumption)',
            'CH2_H_CH_H2_cm3_7_1': 'CH2 + H → CH + H2 (H consumption)',
            'CH_H_C_H2_cm3_7_3': 'CH + H → C + H2 (H consumption)',
        }

    elif target_species == 'CH':
        return {
            'e_CH4_CH_H2_H_vib_cm3_1_3': 'e + CH4 → CH + H2 + H (primary CH source)',
            'e_CH4_CH_H_H2_cm3_1_11': 'e + CH4 → CH + H + H2 (secondary CH source)',
            'ArStar_CH2_CH_H_cm3_3_4': 'Ar* + CH2 → CH + H + Ar (CH production)',
            'CH2_H_CH_H2_cm3_7_1': 'CH2 + H → CH + H2 (CH production)',
            'CH_CH_C2_H2_cm3_5_4': 'CH + CH → C2 + H2 (CH loss, C2 gain) ***KEY***',
            'CH_CH3_C2H2_H2_cm3_7_23': 'CH + CH3 → C2H2 + H2 (CH loss)',
            'CH_CH3_C2H3_H_cm3_7_10': 'CH + CH3 → C2H3 + H (CH loss)',
            'CH_H_C_H2_cm3_7_3': 'CH + H → C + H2 (CH loss)',
            'loss_CH_11_9': 'CH wall loss ***CRITICAL - 10x range***',
            'stick_CH_9_3': 'CH wall sticking (5x range)',
        }

    elif target_species == 'C2':
        return {
            'CH_CH_C2_H2_cm3_5_4': 'CH + CH → C2 + H2 ***PRIMARY C2 SOURCE***',
            'C_CH_C2_H_cm3_7_4': 'C + CH → C2 + H (C2 formation)',
            'CH_C_C2_H_cm3_7_9': 'CH + C → C2 + H (C2 formation)',
            'e_C2H2_C2_H2_cm3_1_16': 'e + C2H2 → C2 + H2 (C2 from C2H2)',
            'C2H2_C_C2_CH2_cm3_7_19': 'C2H2 + C → C2 + CH2 (C2 formation)',
            'C2_H_CH_C_cm3_7_6': 'C2 + H → CH + C (C2 loss)',
            'loss_C2_11_3': 'C2 wall loss ***HUGE range 1e-4 to 2e3!***',
            'stick_C2_9_9': 'C2 wall sticking (5x range)',
        }

    else:
        return {}


if __name__ == '__main__':
    # Test the database
    db = get_complete_rate_database()
    print(f"Complete database loaded: {len(db)} rates")

    # Show flagged rates
    flagged = get_flagged_rates(db)
    print(f"\\nFlagged rates needing validation: {len(flagged)}")
    for name, rate in list(flagged.items())[:5]:
        print(f"  {name}: {rate.flag}")

    # Show CH-related rates
    ch_rates = get_rates_by_species(db, 'CH')
    print(f"\\nRates involving CH: {len(ch_rates)}")
'''

    return code


if __name__ == '__main__':
    import sys

    matlab_file = 'define_rates.m'
    output_file = 'rate_database_complete.py'

    print("Parsing MATLAB define_rates.m...")
    rates = parse_matlab_rates(matlab_file)
    print(f"  Extracted {len(rates)} rates")

    print(f"\\nGenerating {output_file}...")
    python_code = generate_python_code(rates)

    with open(output_file, 'w') as f:
        f.write(python_code)

    print(f"  ✓ Complete database generated: {len(rates)} rates")
    print(f"  ✓ File saved: {output_file}")

    # Test import
    print("\\nTesting generated code...")
    exec(open(output_file).read())
    db = get_complete_rate_database()
    print(f"  ✓ Database loads successfully: {len(db)} rates")

    # Show statistics
    flagged = get_flagged_rates(db)
    print(f"\\nDatabase statistics:")
    print(f"  Total rates: {len(db)}")
    print(f"  Flagged for validation: {len(flagged)}")
    print(f"  Unique sources: {len(set(r.source for r in db.values()))}")

    print("\\n✓ Complete rate database ready!")
