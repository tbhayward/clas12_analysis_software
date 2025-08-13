#!/usr/bin/env python3
"""
Combine asymmetry text files by inverse-variance weighting, preserving the
section order from the first input file.

Usage:
  python combine_asymmetry_texts.py input1.txt input2.txt ... output.txt

All args except the last are inputs; the last arg is the output filename.

Each input file may contain multiple lines of the form:
  NAME = {{mean, value, error}, {mean, value, error}, ...};

For each NAME present in any input, the script aligns points by index across
the input files and computes a weighted average:

  value_comb = sum_i( value_i / sigma_i^2 ) / sum_i( 1 / sigma_i^2 )
  sigma_comb = 1 / sqrt( sum_i( 1 / sigma_i^2 ) )
  mean_comb  = sum_i( mean_i  / sigma_i^2 ) / sum_i( 1 / sigma_i^2 )

Points with NaN or nonpositive sigma are ignored for that index.
If an index has no valid contributions across inputs, it is skipped.
"""

import sys
import re
import math

def parse_file(path, return_order=False):
    """
    Returns:
      - if return_order=False: dict name -> list of (mean, value, sigma)
      - if return_order=True:  (dict, order_list)
    """
    with open(path, 'r') as f:
        text = f.read()

    # Match blocks like: name = { ... };   (across lines)
    assign_re = re.compile(r'([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};', re.DOTALL)
    # Match inner { a, b, c } where there are no further braces inside
    triple_re = re.compile(r'\{([^{}]+)\}')

    out = {}
    order = []

    for m in assign_re.finditer(text):
        name = m.group(1)
        content = m.group(2)
        order.append(name)

        triples = []
        for t in triple_re.findall(content):
            parts = [p.strip() for p in t.split(',')]
            if len(parts) != 3:
                continue
            try:
                mean  = float(parts[0])
                value = float(parts[1])
                sigma = float(parts[2])
            except ValueError:
                continue
            triples.append((mean, value, sigma))
        if triples:
            out[name] = triples

    return (out, order) if return_order else out

def is_valid(meas):
    """Check value tuple (mean, val, sig) has finite numbers and positive sigma."""
    mean, val, sig = meas
    if not (math.isfinite(val) and math.isfinite(sig) and math.isfinite(mean)):
        return False
    if sig <= 0.0:
        return False
    return True

def combine_lists(lists_by_file):
    """
    Combine multiple lists of (mean, value, sigma) aligned by index.

    lists_by_file: list of lists (one per input file)
    Returns a combined list of (mean, value, sigma), skipping indices with no valid contributions.
    """
    max_len = max((len(L) for L in lists_by_file if L is not None), default=0)
    combined = []

    for idx in range(max_len):
        values = []
        for L in lists_by_file:
            if L is None or idx >= len(L):
                continue
            if is_valid(L[idx]):
                values.append(L[idx])

        if not values:
            # skip this index if no valid contributions
            continue

        # inverse-variance weighting
        wsum = 0.0
        vwsum = 0.0
        mwsum = 0.0
        for mean, val, sig in values:
            w = 1.0 / (sig * sig)
            wsum  += w
            vwsum += val  * w
            mwsum += mean * w

        val_comb = vwsum / wsum
        sig_comb = 1.0 / math.sqrt(wsum)
        mean_comb = mwsum / wsum

        combined.append((mean_comb, val_comb, sig_comb))

    return combined

def main():
    if len(sys.argv) < 3:
        print("Usage: python combine_asymmetry_texts.py input1.txt [input2.txt ...] output.txt")
        sys.exit(1)

    *inputs, output_path = sys.argv[1:]
    if not inputs:
        print("Error: Need at least one input file.")
        sys.exit(1)

    # Parse first file (get order) and the rest
    first_dict, first_order = parse_file(inputs[0], return_order=True)
    parsed_rest = [parse_file(p) for p in inputs[1:]]

    # Collect union of all names
    all_names = set(first_dict.keys())
    for d in parsed_rest:
        all_names.update(d.keys())

    # Build final ordered list:
    # 1) in the order they appear in the first file
    # 2) then any names not in the first file, appended in alphabetical order
    ordered_names = []
    seen = set()
    for name in first_order:
        if name in all_names and name not in seen:
            ordered_names.append(name)
            seen.add(name)
    remaining = sorted(n for n in all_names if n not in seen)
    ordered_names.extend(remaining)

    # Combine per name
    combined_results = {}
    for name in ordered_names:
        lists_by_file = []
        # First file contribution
        lists_by_file.append(first_dict.get(name))
        # Remaining files
        for d in parsed_rest:
            lists_by_file.append(d.get(name))
        combo = combine_lists(lists_by_file)
        if combo:
            combined_results[name] = combo

    # Write output in the same format
    with open(output_path, 'w') as out:
        for name in ordered_names:
            triples = combined_results.get(name)
            if not triples:
                continue
            out.write(f"{name} = {{")
            rows = []
            for mean, val, sig in triples:
                rows.append(f"{{{mean:.9f}, {val:.9f}, {sig:.9f}}}")
            out.write(", ".join(rows))
            out.write("};\n")

if __name__ == "__main__":
    main()