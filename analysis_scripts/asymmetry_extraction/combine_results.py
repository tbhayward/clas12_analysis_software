#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine asymmetry text files by inverse-variance weighting.

Two alignment modes:
  1) Default (mean-based): align points by their 'mean' coordinate within a numeric
     tolerance (--tol). Robust to tiny bin-center drifts across periods.
  2) Index-based (--by-index): align strictly by list index (0..N-1). Assumes the
     same bin order/structure across inputs.

Usage:
  python combine_asymmetry_texts.py [--tol 5e-4] [--by-index] input1.txt input2.txt ... output.txt

Notes:
  • Section order in the output follows the order in the FIRST input file; any sections
    that only appear in later inputs are appended alphabetically.
  • In mean-based mode, target 'means' come from:
      1) all 'means' in the FIRST file (in order), plus
      2) any extra 'means' found only in later files (appended sorted), provided they
         are not within tol of an already chosen target.
  • Points with NaN or nonpositive σ are ignored. If a target/index has no valid
    contributions, it is skipped.
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
        text = f.read()

    assign_re = re.compile(r'([A-Za-z0-9_]+)\s*=\s*\{(.*?)\};', re.DOTALL)
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
    return (math.isfinite(mean) and math.isfinite(val) and math.isfinite(sig) and sig > 0.0)

# ── Mean-based alignment helpers ─────────────────────────────────────

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
    return best_idx  # may be None

def build_target_means(anchor_triples, other_lists, tol):
    """
    Create the set of target means for alignment:
      - Start with the means (in order) from the first file (anchor).
      - Add any means from other files that are not within tol of an existing target.
    Returns (targets_in_order, is_anchor_mask) where the mask marks anchor-origin.
    """
    targets = []
    is_anchor = []

    if anchor_triples is not None:
        for m, _, _ in anchor_triples:
            targets.append(m)
            is_anchor.append(True)

    def exists_close(x):
        for t in targets:
            if abs(x - t) <= tol:
                return True
        return False

    extras = []
    for L in other_lists:
        if not L:
            continue
        for m, _, _ in L:
            if not exists_close(m):
                extras.append(m)

    if extras:
        for m in sorted(set(extras)):
            targets.append(m)
            is_anchor.append(False)

    return targets, is_anchor

def combine_by_mean(lists_by_file, tol, want_summary=False, name_for_log=""):
    """
    Combine multiple lists of (mean, value, sigma) by aligning on mean with tolerance.
    lists_by_file: [list_or_None, list_or_None, ...], first element is anchor.
    Returns a combined list of (mean, value, sigma).
    """
    anchor = lists_by_file[0]
    others = lists_by_file[1:]
    targets, _ = build_target_means(anchor, others, tol)

    combined = []
    matched_stats = []

    for T in targets:
        contrib = []
        for L in lists_by_file:
            if not L:
                continue
            k = nearest_index_by_mean(T, L, tol)
            if k is None:
                continue
            if is_valid(L[k]):
                contrib.append(L[k])

        if not contrib:
            continue

        wsum = vwsum = mwsum = 0.0
        for m, v, s in contrib:
            w = 1.0 / (s * s)
            wsum += w
            vwsum += v * w
            mwsum += m * w

        val_comb  = vwsum / wsum
        sig_comb  = 1.0 / math.sqrt(wsum)
        mean_comb = mwsum / wsum

        combined.append((mean_comb, val_comb, sig_comb))
        matched_stats.append(len(contrib))

    if want_summary:
        if matched_stats:
            srt = sorted(matched_stats)
            mn, md, mx = srt[0], srt[len(srt)//2], srt[-1]
            print(f"[COMBINE mean] {name_for_log}: {len(matched_stats)} targets, per-target contributions (min/med/max) = {mn}/{md}/{mx}")
        else:
            print(f"[COMBINE mean] {name_for_log}: no matched targets")

    return combined

# ── Index-based alignment ────────────────────────────────────────────

def combine_by_index(lists_by_file, want_summary=False, name_for_log=""):
    """
    Combine lists of (mean, value, sigma) aligned STRICTLY by index.
    For idx = 0..max_len-1, gather that index from each list (if present/valid),
    then inverse-variance weight. Skip indices with no valid contributions.
    """
    max_len = 0
    for L in lists_by_file:
        if L:
            max_len = max(max_len, len(L))

    combined = []
    matched_stats = []

    for idx in range(max_len):
        contrib = []
        for L in lists_by_file:
            if not L or idx >= len(L):
                continue
            meas = L[idx]
            if is_valid(meas):
                contrib.append(meas)

        if not contrib:
            continue

        wsum = vwsum = mwsum = 0.0
        for m, v, s in contrib:
            w = 1.0 / (s * s)
            wsum += w
            vwsum += v * w
            mwsum += m * w

        val_comb  = vwsum / wsum
        sig_comb  = 1.0 / math.sqrt(wsum)
        mean_comb = mwsum / wsum

        combined.append((mean_comb, val_comb, sig_comb))
        matched_stats.append(len(contrib))

    if want_summary:
        if matched_stats:
            srt = sorted(matched_stats)
            mn, md, mx = srt[0], srt[len(srt)//2], srt[-1]
            print(f"[COMBINE index] {name_for_log}: {len(matched_stats)} indices, per-index contributions (min/med/max) = {mn}/{md}/{mx}")
        else:
            print(f"[COMBINE index] {name_for_log}: no matched indices")

    return combined

# ── Main ─────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Combine asymmetry text files with mean- or index-based alignment.")
    ap.add_argument("--tol", type=float, default=5e-4,
                    help="tolerance for matching means across files (mean-based mode, default: 5e-4)")
    ap.add_argument("--by-index", action="store_true",
                    help="align strictly by index instead of mean matching")
    ap.add_argument("files", nargs="+",
                    help="input1.txt input2.txt ... output.txt (last arg is output)")
    args = ap.parse_args()

    if len(args.files) < 2:
        print("Usage: python combine_asymmetry_texts.py [--tol 5e-4] [--by-index] input1.txt [input2.txt ...] output.txt")
        sys.exit(1)

    *inputs, output_path = args.files
    if not inputs:
        print("Error: Need at least one input file.")
        sys.exit(1)

    # Parse first file (get order) and the rest
    first_dict, first_order = parse_file(inputs[0], return_order=True)
    parsed_rest = [parse_file(p) for p in inputs[1:]]

    # Collect union of all section names
    all_names = set(first_dict.keys())
    for d in parsed_rest:
        all_names.update(d.keys())

    # Output section order: first file order, then alphabetical for the rest
    ordered_names = []
    seen = set()
    for name in first_order:
        if name in all_names and name not in seen:
            ordered_names.append(name)
            seen.add(name)
    ordered_names.extend(sorted(n for n in all_names if n not in seen))

    # Combine per section
    combined_results = {}
    for name in ordered_names:
        lists_by_file = [first_dict.get(name)] + [d.get(name) for d in parsed_rest]

        if args.by_index:
            combo = combine_by_index(lists_by_file, want_summary=True, name_for_log=name)
        else:
            combo = combine_by_mean(lists_by_file, tol=args.tol, want_summary=True, name_for_log=name)

        if combo:
            combined_results[name] = combo

    # Write output in the same format
    with open(output_path, 'w') as out:
        for name in ordered_names:
            triples = combined_results.get(name)
            if not triples:
                continue
            out.write(f"{name} = {{")
            rows = [f"{{{m:.9f}, {v:.9f}, {s:.9f}}}" for (m, v, s) in triples]
            out.write(", ".join(rows))
            out.write("};\n")

if __name__ == "__main__":
    main()