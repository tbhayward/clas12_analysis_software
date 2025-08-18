#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine asymmetry text files by inverse-variance weighting, aligning points by their
'mean' coordinate within a numeric tolerance (rather than assuming identical indexing).

Usage:
  python combine_asymmetry_texts.py [--tol 5e-4] input1.txt input2.txt ... output.txt

Notes:
  • Section order in the output follows the order in the FIRST input file; any sections
    that exist only in later inputs are appended alphabetically (same as your original).
  • Within each section (NAME), target 'means' come from:
      1) all 'means' present in the FIRST file, in their original order;
      2) plus any extra 'means' found only in later files (appended sorted by mean).
  • For each target mean T, we gather the nearest point (within |Δ| ≤ tol) from each file
    and perform inverse-variance weighting across those contributions:
       value_comb = Σ (v_i / σ_i²) / Σ (1 / σ_i²)
       sigma_comb = 1 / sqrt(Σ (1 / σ_i²))
       mean_comb  = Σ (m_i / σ_i²) / Σ (1 / σ_i²)
    Points with NaN or nonpositive σ are ignored. If a target has no valid contributions,
    it is skipped.

This makes the combiner robust to tiny bin-center drifts or rounding differences between periods.
"""

import sys
import re
import math
import argparse

def parse_file(path, return_order=False):
    """
    Returns:
      - if return_order=False: dict name -> list of (mean, value, sigma)
      - if return_order=True:  (dict, order_list)
    """
    with open(path, 'r') as f:
        text = f.read()  #endif

    # Match blocks like: name = { ... };
    assign_re = re.compile(r'([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};', re.DOTALL)
    # Match inner { a, b, c } where there are no further braces inside
    triple_re = re.compile(r'\{([^{}]+)\}')  #endif

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
                continue  # endif
            try:
                mean  = float(parts[0])
                value = float(parts[1])
                sigma = float(parts[2])
            except ValueError:
                continue  # endif
            triples.append((mean, value, sigma))
        #endfor
        if triples:
            out[name] = triples
    #endfor

    return (out, order) if return_order else out  #endif
#endfor

def is_valid(meas):
    """Check value tuple (mean, val, sig) has finite numbers and positive sigma."""
    mean, val, sig = meas
    if not (math.isfinite(val) and math.isfinite(sig) and math.isfinite(mean)):
        return False  # endif
    if sig <= 0.0:
        return False  # endif
    return True  # endif
#endfor

def nearest_index_by_mean(target_mean, triples, tol):
    """
    Find index of the element in 'triples' whose mean is nearest to target_mean,
    provided |Δ| ≤ tol. Returns None if no element is within tol.
    """
    best_idx = None
    best_delta = None
    for i, (m, v, s) in enumerate(triples):
        delta = abs(m - target_mean)
        if delta <= tol and (best_delta is None or delta < best_delta):
            best_idx = i
            best_delta = delta
        #endif
    #endfor
    return best_idx  # may be None
#endfor

def build_target_means(anchor_triples, other_lists, tol):
    """
    Create the set of target means for alignment:
      - Start with the means (in order) from the first file (anchor).
      - Add any means from other files that are not within tol of an existing target.
    Returns (targets_in_order, is_anchor_mask) where the mask is a list[bool]
    indicating which targets came from the anchor.
    """
    targets = []
    is_anchor = []

    # 1) Anchor means, keep order
    if anchor_triples is not None:
        for m, _, _ in anchor_triples:
            targets.append(m)
            is_anchor.append(True)
        #endfor
    #endif

    # 2) Extras from other files
    def exists_close(x):
        for t in targets:
            if abs(x - t) <= tol:
                return True  # endif
        #endfor
        return False  # endif
    #endfor

    extras = []
    for L in other_lists:
        if not L:
            continue  # endif
        for m, _, _ in L:
            if not exists_close(m):
                extras.append(m)
            #endif
        #endfor
    #endfor

    if extras:
        extras_sorted = sorted(set(extras))
        for m in extras_sorted:
            targets.append(m)
            is_anchor.append(False)
        #endfor
    #endif

    return targets, is_anchor  #endif
#endfor

def combine_by_mean(lists_by_file, tol, want_summary=False, name_for_log=""):
    """
    Combine multiple lists of (mean, value, sigma) by aligning on the mean with tolerance.

    lists_by_file: [list_or_None, list_or_None, ...] parallel to the input files.
                   The FIRST element is the anchor list that sets the base order.

    Returns a combined list of (mean, value, sigma).
    """
    anchor = lists_by_file[0]
    others = lists_by_file[1:]

    # Targets: anchor means in order + unmatched means from others (sorted)
    targets, is_anchor = build_target_means(anchor, others, tol)

    combined = []
    matched_stats = []  # for optional logging: how many files contributed per target

    for ti, T in enumerate(targets):
        contrib = []
        # Examine each file's list, find nearest point to T within tol
        for L in lists_by_file:
            if not L:
                continue  # endif
            k = nearest_index_by_mean(T, L, tol)
            if k is None:
                continue  # endif
            if is_valid(L[k]):
                contrib.append(L[k])
            #endif
        #endfor

        if not contrib:
            continue  # skip targets with no valid contributions  # endif

        # Inverse-variance weighting on contrib
        wsum = 0.0
        vwsum = 0.0
        mwsum = 0.0
        for m, v, s in contrib:
            w = 1.0 / (s * s)
            wsum  += w
            vwsum += v * w
            mwsum += m * w
        #endfor

        val_comb  = vwsum / wsum
        sig_comb  = 1.0 / math.sqrt(wsum)
        mean_comb = mwsum / wsum

        combined.append((mean_comb, val_comb, sig_comb))
        matched_stats.append(len(contrib))
    #endfor

    if want_summary and targets:
        # Simple one-line summary for quick sanity checking
        # e.g., "NAME: 12 targets, contributions per target (min/median/max) = 2/3/3"
        if matched_stats:
            srt = sorted(matched_stats)
            mn, md, mx = srt[0], srt[len(srt)//2], srt[-1]
            print(f"[COMBINE] {name_for_log}: {len(matched_stats)} targets, per-target contributions (min/med/max) = {mn}/{md}/{mx}")
        else:
            print(f"[COMBINE] {name_for_log}: no matched targets")
        #endif
    #endif

    return combined  #endif
#endfor

def main():
    ap = argparse.ArgumentParser(description="Combine asymmetry text files with mean-based alignment.")
    ap.add_argument("--tol", type=float, default=5e-4, help="tolerance for matching means across files (default: 5e-4)")
    ap.add_argument("files", nargs="+", help="input1.txt input2.txt ... output.txt (last arg is output)")
    args = ap.parse_args()

    if len(args.files) < 2:
        print("Usage: python combine_asymmetry_texts.py [--tol 5e-4] input1.txt [input2.txt ...] output.txt")
        sys.exit(1)  # endif

    *inputs, output_path = args.files
    if not inputs:
        print("Error: Need at least one input file.")
        sys.exit(1)  # endif

    # Parse first file (get order) and the rest
    first_dict, first_order = parse_file(inputs[0], return_order=True)
    parsed_rest = [parse_file(p) for p in inputs[1:]]  #endif

    # Collect union of all names
    all_names = set(first_dict.keys())
    for d in parsed_rest:
        all_names.update(d.keys())
    #endfor

    # Build final ordered list:
    # 1) in the order they appear in the first file
    # 2) then any names not in the first file, appended in alphabetical order
    ordered_names = []
    seen = set()
    for name in first_order:
        if name in all_names and name not in seen:
            ordered_names.append(name)
            seen.add(name)
        #endif
    #endfor
    remaining = sorted(n for n in all_names if n not in seen)
    ordered_names.extend(remaining)

    # Combine per name
    combined_results = {}
    for name in ordered_names:
        lists_by_file = []
        # First file contribution (may be None if name absent)
        lists_by_file.append(first_dict.get(name))
        # Remaining files
        for d in parsed_rest:
            lists_by_file.append(d.get(name))
        #endfor

        combo = combine_by_mean(lists_by_file, tol=args.tol, want_summary=True, name_for_log=name)
        if combo:
            combined_results[name] = combo
        #endif
    #endfor

    # Write output in the same format
    with open(output_path, 'w') as out:
        for name in ordered_names:
            triples = combined_results.get(name)
            if not triples:
                continue  # endif
            out.write(f"{name} = {{")
            rows = []
            for mean, val, sig in triples:
                rows.append(f"{{{mean:.9f}, {val:.9f}, {sig:.9f}}}")
            #endfor
            out.write(", ".join(rows))
            out.write("};\n")
        #endfor
    #endif
#endfor

if __name__ == "__main__":
    main()
#endif