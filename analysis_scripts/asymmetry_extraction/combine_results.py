#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine asymmetry text files by inverse-variance weighting.

SPECIAL RULE (requested):
    For each matched point:
      • Let file ordering be [Su22, Fa22, Sp23, ...] (i.e., your first, second, third inputs).
      • If the uncertainty σ in Su22 (file 0) OR in Sp23 (file 2) is smaller than the
        uncertainty σ in Fa22 (file 1), then that contribution is replaced by
        (value = 0, sigma = 1) for that file before combining.

Two alignment modes:
  1) Default (mean-based): align points by their 'mean' coordinate within a numeric
     tolerance (--tol). Robust to tiny bin-center drifts across periods.
  2) Index-based (--by-index): align strictly by list index (0..N-1). Assumes the
     same bin order/structure across inputs.

Usage:
  python combine_asymmetry_texts.py [--tol 5e-4] [--by-index] input1.txt input2.txt ... output.txt
"""

import sys
import re
import math
import argparse

# ─────────────────────────────────────────────────────────────────────

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

# ── Combination cores (with Su22/Fa22/Sp23 σ-floor rule) ────────────

def _apply_fa22_sigma_floor(matches):
    """
    matches: list with per-file matched measurements (or None).
             Each measurement is (mean, value, sigma).
    Rule:
      If file 1 (Fa22) exists and is valid, and file 0 (Su22) OR file 2 (Sp23)
      exists+valid with sigma < sigma_Fa22, then replace that file's contribution
      by (same mean, value=0.0, sigma=1.0).
    """
    if len(matches) < 2:
        return matches  # Need at least Fa22 to compare

    fa = matches[1]
    if fa is None or not is_valid(fa):
        return matches

    sigma_fa = fa[2]
    for idx in (0, 2):
        if idx < len(matches) and matches[idx] is not None and is_valid(matches[idx]):
            m, v, s = matches[idx]
            if s < sigma_fa:
                matches[idx] = (m, 0.0, 1.0)
    return matches

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
        # Collect the matched element from each file (or None)
        matches = []
        for L in lists_by_file:
            if not L:
                matches.append(None)
                continue
            k = nearest_index_by_mean(T, L, tol)
            matches.append(L[k] if k is not None else None)

        # Apply the σ-floor rule using Fa22 as reference
        matches = _apply_fa22_sigma_floor(matches)

        # Keep valid contributions only
        contrib = [meas for meas in matches if (meas is not None and is_valid(meas))]
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

def combine_by_index(lists_by_file, want_summary=False, name_for_log=""):
    """
    Combine lists of (mean, value, sigma) aligned STRICTLY by index.
    For idx = 0..max_len-1, gather that index from each list (if present/valid),
    apply the σ-floor rule, then inverse-variance weight. Skip indices with no valid contributions.
    """
    max_len = 0
    for L in lists_by_file:
        if L:
            max_len = max(max_len, len(L))

    combined = []
    matched_stats = []

    for idx in range(max_len):
        # Collect index-aligned element from each file (or None)
        matches = []
        for L in lists_by_file:
            if not L or idx >= len(L):
                matches.append(None)
            else:
                matches.append(L[idx])

        # Apply the σ-floor rule using Fa22 as reference
        matches = _apply_fa22_sigma_floor(matches)

        contrib = [meas for meas in matches if (meas is not None and is_valid(meas))]
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