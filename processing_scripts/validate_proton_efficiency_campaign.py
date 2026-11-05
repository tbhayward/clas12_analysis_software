#!/usr/bin/env python3
"""Validate proton-efficiency campaign ROOT parts, markers, and tag cutflows.

Run after loading the same ROOT environment used by the campaign.  The script
is read-only and returns nonzero if any required invariant fails.
"""

import argparse
import json
import math
import pathlib
import re
import sys

try:
    import ROOT
except ImportError as error:
    raise SystemExit("PyROOT is unavailable; load the ROOT environment first: " + str(error))

P = pathlib.Path
PERIODS = ("fa18_inb", "fa18_out", "sp18_inb", "sp18_out")
MODES = ("data", "mc")
TAG_BRANCHES = {
    "file_index", "event_sequence", "tag_sequence", "run", "event", "torus",
    "beam_energy", "e_p", "rec_e_p_gt_2", "Q2", "W", "xB", "y",
    "missing_p", "missing_theta", "missing_phi", "n_proton_eb",
    "n_proton_selected", "truth_tag_match_valid", "truth_e_p",
    "truth_e_p_gt_2", "truth_missing_p", "truth_n_proton",
}
PROTON_BRANCHES = {
    "tag_sequence", "rec_index", "selected", "vertex_pass", "fiducial_pass",
    "momentum_pass", "status", "detector", "p", "theta", "phi",
    "match_angle", "delta_p", "relative_delta_p",
}
CUTFLOW_KEYS = (
    "events", "banks", "qa", "raw_epi+pi-", "fitted_e", "fitted_pi+",
    "fitted_pi-", "fitted_epi+pi-", "combinations", "q2", "w",
    "y_upper", "y_lower", "tags",
)


def branch_names(tree):
    return {branch.GetName() for branch in tree.GetListOfBranches()}


def parse_cutflow(log_path):
    if not log_path.is_file():
        return None
    lines = [line for line in log_path.read_text(errors="replace").splitlines()
             if line.startswith("[TAG_CUTFLOW]")]
    if len(lines) != 1:
        return None
    result = {}
    for key in CUTFLOW_KEYS:
        match = re.search(r"(?:^|\s)" + re.escape(key) + r"=(\d+)(?:\s|$)", lines[0])
        if match:
            result[key] = int(match.group(1))
    return result if len(result) == len(CUTFLOW_KEYS) else None


def first_zero_transition(cutflow):
    stages = (
        ("input events", "events"),
        ("required banks", "banks"),
        ("QADB/run quality", "qa"),
        ("raw e pi+ pi- topology", "raw_epi+pi-"),
        ("fitted e pi+ pi- topology", "fitted_epi+pi-"),
        ("particle combinations", "combinations"),
        ("Q2 > 1", "q2"),
        ("W > 2", "w"),
        ("y < 0.85", "y_upper"),
        ("y > 0.19", "y_lower"),
        ("selected tags", "tags"),
    )
    for label, key in stages:
        if cutflow.get(key, 0) == 0:
            return label
    return "none"


def finite_and_range(value, low, high):
    return math.isfinite(value) and low <= value <= high


def validate_part(root_path, period, mode, errors, warnings):
    prefix = str(root_path)[:-5]
    marker_path = P(prefix + ".ok.json")
    log_path = P(prefix + ".processing.log")
    if not marker_path.is_file():
        errors.append(f"{period}/{mode}: missing marker for {root_path.name}")
        return None
    try:
        marker = json.loads(marker_path.read_text())
    except (OSError, ValueError) as error:
        errors.append(f"{period}/{mode}: unreadable marker {marker_path}: {error}")
        return None

    source = ROOT.TFile.Open(str(root_path), "READ")
    if not source or source.IsZombie():
        errors.append(f"{period}/{mode}: cannot open {root_path}")
        return None
    tags = source.Get("Tags")
    protons = source.Get("Protons")
    if not tags or not protons:
        errors.append(f"{period}/{mode}: {root_path.name} lacks Tags or Protons")
        source.Close()
        return None
    missing_tags = sorted(TAG_BRANCHES - branch_names(tags))
    missing_protons = sorted(PROTON_BRANCHES - branch_names(protons))
    if missing_tags:
        errors.append(f"{period}/{mode}: Tags missing {missing_tags}")
    if missing_protons:
        errors.append(f"{period}/{mode}: Protons missing {missing_protons}")
    if missing_tags or missing_protons:
        source.Close()
        return None

    n_tags = int(tags.GetEntries())
    n_protons = int(protons.GetEntries())
    if marker.get("selected_tags") != n_tags:
        errors.append(f"{period}/{mode}: marker/Tags mismatch {marker.get('selected_tags')} != {n_tags}")
    if marker.get("proton_rows") != n_protons:
        errors.append(f"{period}/{mode}: marker/Protons mismatch {marker.get('proton_rows')} != {n_protons}")
    if marker.get("period") != period or marker.get("mode") != mode:
        errors.append(f"{period}/{mode}: marker group identity is inconsistent")

    cutflow = parse_cutflow(log_path)
    if cutflow is None:
        errors.append(f"{period}/{mode}: missing or invalid TAG_CUTFLOW in {log_path.name}")
    elif cutflow["tags"] != n_tags:
        errors.append(f"{period}/{mode}: cutflow/Tags mismatch {cutflow['tags']} != {n_tags}")
    if n_tags == 0:
        stage = first_zero_transition(cutflow) if cutflow else "unknown"
        detail = " ".join(f"{key}={cutflow[key]}" for key in CUTFLOW_KEYS) if cutflow else "cutflow unavailable"
        warnings.append(f"{period}/{mode}: ZERO TAGS; first zero stage: {stage}; {detail}")

    tag_ids = set()
    sum_eb = 0
    sum_selected = 0
    bad_tag_values = 0
    truth_matches = 0
    truth_ep_gt_2 = 0
    missing_theta_min = math.inf
    missing_theta_max = -math.inf
    for entry in tags:
        tag_id = int(round(entry.tag_sequence))
        if tag_id in tag_ids:
            errors.append(f"{period}/{mode}: duplicate tag_sequence {tag_id} in {root_path.name}")
            break
        tag_ids.add(tag_id)
        sum_eb += int(round(entry.n_proton_eb))
        sum_selected += int(round(entry.n_proton_selected))
        if not (entry.rec_e_p_gt_2 == 1 and entry.e_p > 2.0 and entry.Q2 >= 1.0 and
                entry.W >= 2.0 and 0.19 <= entry.y <= 0.85 and
                finite_and_range(entry.missing_theta, 0.0, 180.0) and
                finite_and_range(entry.missing_phi, 0.0, 360.0) and
                math.isfinite(entry.missing_p)):
            bad_tag_values += 1
        missing_theta_min = min(missing_theta_min, float(entry.missing_theta))
        missing_theta_max = max(missing_theta_max, float(entry.missing_theta))
        if entry.truth_tag_match_valid == 1:
            truth_matches += 1
            if entry.truth_e_p_gt_2 == 1:
                truth_ep_gt_2 += 1
        elif mode == "data" and entry.truth_tag_match_valid != 0:
            bad_tag_values += 1
    if bad_tag_values:
        errors.append(f"{period}/{mode}: {bad_tag_values} Tags violate required kinematic/schema invariants")

    linked_protons = 0
    selected_protons = 0
    bad_proton_values = 0
    proton_theta_min = math.inf
    proton_theta_max = -math.inf
    for entry in protons:
        if int(round(entry.tag_sequence)) in tag_ids:
            linked_protons += 1
        if entry.selected == 1:
            selected_protons += 1
        if entry.selected not in (0, 1) or entry.vertex_pass not in (0, 1) or \
           entry.fiducial_pass not in (0, 1) or entry.momentum_pass not in (0, 1) or \
           not finite_and_range(entry.theta, 0.0, 180.0) or \
           not finite_and_range(entry.phi, 0.0, 360.0) or \
           not finite_and_range(entry.match_angle, 0.0, 180.0) or not math.isfinite(entry.p):
            bad_proton_values += 1
        proton_theta_min = min(proton_theta_min, float(entry.theta))
        proton_theta_max = max(proton_theta_max, float(entry.theta))
    if linked_protons != n_protons:
        errors.append(f"{period}/{mode}: {n_protons-linked_protons} Protons have no matching tag_sequence")
    if sum_eb != n_protons:
        errors.append(f"{period}/{mode}: sum(n_proton_eb)={sum_eb}, Protons rows={n_protons}")
    if sum_selected != selected_protons:
        errors.append(f"{period}/{mode}: selected summary={sum_selected}, candidate rows={selected_protons}")
    if bad_proton_values:
        errors.append(f"{period}/{mode}: {bad_proton_values} Protons violate required invariants")

    truth_fraction = truth_matches / n_tags if n_tags else float("nan")
    truth_ep_fraction = truth_ep_gt_2 / truth_matches if truth_matches else float("nan")
    result = {
        "period": period,
        "mode": mode,
        "file": root_path.name,
        "events": marker.get("processed_events"),
        "tags": n_tags,
        "protons": n_protons,
        "selected_protons": selected_protons,
        "cutflow": cutflow,
        "truth_matches": truth_matches,
        "truth_ep_gt_2": truth_ep_gt_2,
        "truth_match_fraction": truth_fraction,
        "truth_ep_gt_2_fraction": truth_ep_fraction,
        "missing_theta_range": None if not n_tags else [missing_theta_min, missing_theta_max],
        "proton_theta_range": None if not n_protons else [proton_theta_min, proton_theta_max],
    }
    source.Close()
    return result


def fmt_fraction(value):
    return "n/a" if not math.isfinite(value) else f"{100.0*value:.1f}%"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root")
    parser.add_argument("--json", dest="json_output")
    parser.add_argument("--allow-zero-tags", action="store_true",
                        help="Report zero-tag groups as warnings instead of a failed validation")
    args = parser.parse_args()
    root = P(args.output_root)
    errors = []
    warnings = []
    results = []

    for period in PERIODS:
        for mode in MODES:
            parts = root / period / mode / "parts"
            files = sorted(parts.glob("*.root")) if parts.is_dir() else []
            if not files:
                errors.append(f"{period}/{mode}: no ROOT parts found")
                continue
            for path in files:
                result = validate_part(path, period, mode, errors, warnings)
                if result:
                    results.append(result)

    print("period      mode files events      tags    protons selected truth-match truth-e>2")
    for period in PERIODS:
        for mode in MODES:
            group = [result for result in results if result["period"] == period and result["mode"] == mode]
            group_tags = sum(result["tags"] for result in group)
            group_truth_matches = sum(result["truth_matches"] for result in group)
            group_truth_ep_gt_2 = sum(result["truth_ep_gt_2"] for result in group)
            print(f"{period:11} {mode:4} {len(group):5d} "
                  f"{sum(result['events'] or 0 for result in group):10d} "
                  f"{group_tags:9d} "
                  f"{sum(result['protons'] for result in group):10d} "
                  f"{sum(result['selected_protons'] for result in group):8d} "
                  f"{fmt_fraction(group_truth_matches/group_tags if group_tags else float('nan')):>11} "
                  f"{fmt_fraction(group_truth_ep_gt_2/group_truth_matches if group_truth_matches else float('nan')):>9}")

    for warning in warnings:
        print("[WARNING]", warning)
    if warnings and not args.allow_zero_tags:
        errors.extend(warnings)
    for error in errors:
        print("[ERROR]", error)

    report = {"results": results, "warnings": warnings, "errors": errors}
    if args.json_output:
        P(args.json_output).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if errors:
        print(f"[VALIDATION_FAILED] errors={len(errors)} warnings={len(warnings)}")
        return 1
    print(f"[VALIDATION_OK] parts={len(results)} warnings={len(warnings)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
