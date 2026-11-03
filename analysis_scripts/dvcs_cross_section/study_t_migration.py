#!/usr/bin/env python3
"""Paired-MC |t| migration, closure, and correction-method feasibility study.

This diagnostic is deliberately separate from the production correction.  It
uses the generated and reconstructed records in the same MCEvents entry,
applies the production reconstructed-MC selection, and compares ordinary and
iteratively reweighted bin-by-bin acceptance with iterative Bayesian unfolding.
Independent half-sample closure, spectrum stresses, toy uncertainties,
correlations, topology response, and regularization scans are written as plots
and machine-readable tables.
"""

from __future__ import annotations

import argparse
import concurrent.futures as futures
import csv
import json
import math
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
try:
    import uproot
except ImportError:  # Allow --help and static checks outside the farm environment.
    uproot = None


PERIODS = {
    "sp18_inb": ("Sp18 Inb", "50nA"),
    "sp18_out": ("Sp18 Out", "45nA"),
    "fa18_inb": ("Fa18 Inb", "50nA"),
    "fa18_out": ("Fa18 Out", "50nA"),
    "sp19_inb": ("Sp19 Inb", "50nA"),
}
PERIOD_CODE = {
    "sp18_inb": "Sp18_Inb", "sp18_out": "Sp18_Out",
    "fa18_inb": "Fa18_Inb", "fa18_out": "Fa18_Out",
    "sp19_inb": "Sp19_Inb",
}
TOPOLOGIES = {(1, 1): "FD_FD", (2, 1): "CD_FD", (2, 0): "CD_FT"}
CUT_VARIABLES = ("Delta_phi", "theta", "theta_gamma_gamma", "pTmiss", "Emiss2", "Mx2", "Mx2_2")

# Extra sub-threshold bins make feed-in to the first published bin visible.
T_EDGES = np.asarray([0.0, 0.05, 0.08, 0.11, 0.15, 0.25, 0.40, 0.60, 0.80, 1.00])
T_LABELS = (
    "<0", "0-.05", ".05-.08", ".08-.11", ".11-.15",
    ".15-.25", ".25-.40", ".40-.60", ".60-.80", ".80-1.0", ">=1.0",
)
MISS_LABEL = "miss/not selected"
DELTA_EDGES = np.linspace(-0.50, 0.50, 101)
CUTFLOW_NAMES = ("all entries", "valid truth", "valid reco", "valid topology", "global cuts", "exclusivity", "selected")


def t_category(t_abs: np.ndarray) -> np.ndarray:
    return np.searchsorted(T_EDGES, t_abs, side="right").astype(np.int16)


def fd_sector(phi_rad: np.ndarray) -> np.ndarray:
    p = np.mod(np.degrees(phi_rad), 360.0)
    out = np.full(p.shape, -1, dtype=np.int16)
    out[(p >= 330.0) | (p < 30.0)] = 1
    out[(p >= 30.0) & (p < 90.0)] = 2
    out[(p >= 90.0) & (p < 150.0)] = 3
    out[(p >= 150.0) & (p < 210.0)] = 4
    out[(p >= 210.0) & (p < 270.0)] = 5
    out[(p >= 270.0) & (p < 330.0)] = 6
    return out


def cd_sector(phi_rad: np.ndarray) -> np.ndarray:
    p = np.mod(np.degrees(phi_rad), 360.0)
    out = np.full(p.shape, -1, dtype=np.int16)
    out[(p >= 272.5) | (p < 32.5)] = 1
    out[(p >= 32.5) & (p < 150.5)] = 2
    out[(p >= 150.5) & (p < 272.5)] = 3
    return out


def branch(tree_keys: set[str], logical: str, prefix: str = "rec_") -> str:
    candidate = prefix + logical
    if candidate in tree_keys:
        return candidate
    if logical in {"runnum", "detector1", "detector2"} and logical in tree_keys:
        return logical
    return candidate


def required_branches(keys: set[str], cuts: Mapping[str, object]) -> Tuple[Dict[str, str], Sequence[str]]:
    logical = {
        "runnum", "detector1", "detector2", "t1", "open_angle_ep2", "pTmiss",
        "e_phi", "p1_phi", "p2_phi", "e_theta", "p1_theta", "p2_theta",
        "p1_p", "p2_p", *CUT_VARIABLES,
    }
    names = {name: branch(keys, name) for name in logical}
    if names["Delta_phi"] not in keys:
        names["Delta_phi"] = ""  # C++ fallback uses p1_phi and p2_phi.
    physical = {"gen_kinematics_valid", "rec_kinematics_valid", "gen_t1"}
    physical.update(value for value in names.values() if value)
    missing = sorted(item for item in physical if item not in keys)
    if missing:
        raise KeyError(f"paired tree is missing required branches: {missing}")
    return names, sorted(physical)


def global_mask(a: Mapping[str, np.ndarray], n: Mapping[str, str], period: str, cfg: Mapping[str, object]) -> np.ndarray:
    if bool(cfg.get("enable_dvcsgen_ycol_cut", False)):
        raise RuntimeError("The migration study does not yet support an enabled dvcsgen ycol cut")
    d1 = np.asarray(a[n["detector1"]], dtype=np.int16)
    d2 = np.asarray(a[n["detector2"]], dtype=np.int16)
    t1 = np.asarray(a[n["t1"]], dtype=float)
    opening = np.asarray(a[n["open_angle_ep2"]], dtype=float)
    ptmiss = np.asarray(a[n["pTmiss"]], dtype=float)
    mask = np.isfinite(t1) & np.isfinite(opening) & np.isfinite(ptmiss)
    mask &= (-t1) < float(cfg["t1_abs_max"])
    mask &= opening > float(cfg["open_angle_min_deg"])
    if bool(cfg.get("enable_pTmiss_cut", False)):
        mask &= ptmiss <= float(cfg["pTmiss_max"])
    excluded = np.asarray(cfg.get("excluded_runs", []), dtype=np.int64)
    if excluded.size:
        mask &= ~np.isin(np.asarray(a[n["runnum"]], dtype=np.int64), excluded)
    if bool(cfg.get("enable_topology_filter", False)):
        mask &= (d1 == int(cfg["required_detector1"])) & (d2 == int(cfg["required_detector2"]))

    e_phi = np.asarray(a[n["e_phi"]], dtype=float)
    p_phi = np.asarray(a[n["p1_phi"]], dtype=float)
    g_phi = np.asarray(a[n["p2_phi"]], dtype=float)
    e_sec, p_fd_sec, p_cd_sec, g_sec = fd_sector(e_phi), fd_sector(p_phi), cd_sector(p_phi), fd_sector(g_phi)
    if bool(cfg.get("enable_electron_fd_sector_filter", False)):
        mask &= e_sec == int(cfg["electron_fd_sector"])
    if bool(cfg.get("enable_proton_fd_sector_filter", False)):
        mask &= (d1 == 1) & (p_fd_sec == int(cfg["proton_fd_sector"]))
    if bool(cfg.get("enable_proton_cd_sector_filter", False)):
        mask &= (d1 == 2) & (p_cd_sec == int(cfg["proton_cd_sector"]))
    if bool(cfg.get("enable_photon_fd_sector_filter", False)):
        mask &= (d2 == 1) & (g_sec == int(cfg["photon_fd_sector"]))

    if not bool(cfg.get("enable_auxiliary_fiducial_cuts", False)):
        return mask
    et = np.degrees(np.asarray(a[n["e_theta"]], dtype=float))
    pt = np.degrees(np.asarray(a[n["p1_theta"]], dtype=float))
    gt = np.degrees(np.asarray(a[n["p2_theta"]], dtype=float))
    pp = np.asarray(a[n["p1_p"]], dtype=float)
    gp = np.asarray(a[n["p2_p"]], dtype=float)
    mask &= np.isfinite(et) & np.isfinite(pt) & np.isfinite(gt) & np.isfinite(pp) & np.isfinite(gp)
    mask &= et > float(cfg["auxiliary_e_theta_min_deg"])
    mask &= (d2 != 1) | (et < float(cfg["auxiliary_e_theta_max_deg"]))
    fd_pmin = float(cfg["auxiliary_fd_proton_p_min_inb_GeV"] if period.endswith("_inb")
                    else cfg["auxiliary_fd_proton_p_min_out_GeV"])
    proton_ok = ((d1 == 1) & (pp > fd_pmin)
                 & (pt > float(cfg["auxiliary_fd_proton_theta_min_deg"]))
                 & (pt < float(cfg["auxiliary_fd_proton_theta_max_deg"])))
    proton_ok |= ((d1 == 2) & (pp > float(cfg["auxiliary_cd_proton_p_min_GeV"]))
                  & (pt > float(cfg["auxiliary_cd_proton_theta_min_deg"]))
                  & (pt < float(cfg["auxiliary_cd_proton_theta_max_deg"])))
    mask &= proton_ok
    photon_ok = ((d2 == 1)
                 & ((d1 != 1) | (gt > float(cfg["auxiliary_fd_photon_theta_min_deg"])))
                 & (gt < float(cfg["auxiliary_fd_photon_theta_max_deg"])))
    photon_ok |= (d2 == 0) & (gp > float(cfg["auxiliary_ft_photon_p_min_GeV"]))
    mask &= photon_ok
    if bool(cfg.get("auxiliary_require_distinct_fd_sectors", False)):
        mask &= e_sec >= 1
        mask &= (d1 != 1) | ((p_fd_sec >= 1) & (p_fd_sec != e_sec))
        mask &= (d2 != 1) | ((g_sec >= 1) & (g_sec != e_sec))
        mask &= ~((d1 == 1) & (d2 == 1) & (p_fd_sec == g_sec))
    return mask


def topology_codes(a: Mapping[str, np.ndarray], n: Mapping[str, str]) -> Tuple[np.ndarray, np.ndarray]:
    d1 = np.asarray(a[n["detector1"]], dtype=np.int16)
    d2 = np.asarray(a[n["detector2"]], dtype=np.int16)
    valid = np.zeros(d1.shape, dtype=bool)
    code = np.full(d1.shape, -1, dtype=np.int8)
    for index, pair in enumerate(TOPOLOGIES):
        here = (d1 == pair[0]) & (d2 == pair[1])
        valid |= here
        code[here] = index
    return valid, code


def exclusivity_mask(a: Mapping[str, np.ndarray], n: Mapping[str, str], period: str,
                     topology: np.ndarray, cuts: Mapping[str, object]) -> np.ndarray:
    result = np.zeros(topology.shape, dtype=bool)
    for topo_index, topo_name in enumerate(TOPOLOGIES.values()):
        here = topology == topo_index
        if not np.any(here):
            continue
        key = f"DVCS_{PERIOD_CODE[period]}_{topo_name}"
        if key not in cuts or "mc" not in cuts[key]:
            raise KeyError(f"missing reconstructed-MC exclusivity block {key}")
        passed = here.copy()
        for variable, stats in cuts[key]["mc"].items():
            if variable not in CUT_VARIABLES:
                raise KeyError(f"unsupported production cut variable {variable} in {key}")
            if variable == "Delta_phi" and not n["Delta_phi"]:
                delta = np.asarray(a[n["p1_phi"]], dtype=float) - np.asarray(a[n["p2_phi"]], dtype=float)
                values = np.abs((delta + np.pi) % (2.0 * np.pi) - np.pi)
            else:
                values = np.asarray(a[n[variable]], dtype=float)
            low = float(stats.get("cut_low", math.nan))
            high = float(stats.get("cut_high", math.nan))
            if not (math.isfinite(low) and math.isfinite(high) and high > low):
                mean, std = float(stats.get("mean", math.nan)), float(stats.get("std", math.nan))
                if not (math.isfinite(mean) and math.isfinite(std) and std > 0.0):
                    continue
                low, high = mean - 3.0 * std, mean + 3.0 * std
            passed &= np.isfinite(values) & (values >= low) & (values <= high)
        result |= passed
    return result


def empty_result() -> Dict[str, np.ndarray]:
    ncat = len(T_LABELS)
    return {
        "matrix": np.zeros((ncat + 1, ncat), dtype=np.int64),
        "train_matrix": np.zeros((ncat + 1, ncat), dtype=np.int64),
        "test_matrix": np.zeros((ncat + 1, ncat), dtype=np.int64),
        "delta": np.zeros((ncat, len(DELTA_EDGES) - 1), dtype=np.int64),
        "cutflow": np.zeros(len(CUTFLOW_NAMES), dtype=np.int64),
        "topology": np.zeros((len(TOPOLOGIES), ncat), dtype=np.int64),
        "topology_matrix": np.zeros((len(TOPOLOGIES), ncat, ncat), dtype=np.int64),
    }


def add_result(target: Dict[str, np.ndarray], source: Mapping[str, np.ndarray]) -> None:
    for key in target:
        target[key] += source[key]


def process_file(task: Tuple[str, str, str, dict, dict, str]) -> Tuple[str, Dict[str, np.ndarray]]:
    period, filename, step_size, cfg, cuts, _ = task
    out = empty_result()
    with uproot.open(filename) as root_file:
        if "MCEvents" not in root_file:
            raise KeyError(f"MCEvents tree missing in {filename}")
        keys = set(root_file["MCEvents"].keys())
        names, expressions = required_branches(keys, cuts)
    spec = f"{filename}:MCEvents"
    entry_offset = 0
    for a in uproot.iterate(spec, expressions=expressions, step_size=step_size, library="np"):
        size = len(a["gen_kinematics_valid"])
        # A deterministic even/odd entry split provides statistically independent
        # response-training and closure-test samples without requiring an event ID.
        is_train = ((entry_offset + np.arange(size, dtype=np.int64)) % 2) == 0
        entry_offset += size
        gt = -np.asarray(a["gen_t1"], dtype=float)
        rt = -np.asarray(a[names["t1"]], dtype=float)
        gen_valid = (np.asarray(a["gen_kinematics_valid"], dtype=float) >= 0.5) & np.isfinite(gt)
        rec_valid = np.asarray(a["rec_kinematics_valid"], dtype=float) >= 0.5
        topo_valid, topo_code = topology_codes(a, names)
        gmask = global_mask(a, names, period, cfg)
        emask = exclusivity_mask(a, names, period, topo_code, cuts)
        exclusivity_selected = gen_valid & rec_valid & topo_valid & gmask & emask
        selected = exclusivity_selected & np.isfinite(rt)
        stages = (np.ones(size, bool), gen_valid, gen_valid & rec_valid,
                  gen_valid & rec_valid & topo_valid,
                  gen_valid & rec_valid & topo_valid & gmask,
                  exclusivity_selected, selected)
        out["cutflow"] += np.asarray([np.count_nonzero(stage) for stage in stages], dtype=np.int64)

        if not np.any(gen_valid):
            continue
        gcat = t_category(gt)
        rcat = t_category(rt)
        valid_indices = np.flatnonzero(gen_valid)
        selected_valid = selected[valid_indices]
        reco_outcome = np.where(selected_valid, rcat[valid_indices], len(T_LABELS))
        np.add.at(out["matrix"], (reco_outcome, gcat[valid_indices]), 1)
        train_valid = is_train[valid_indices]
        np.add.at(out["train_matrix"],
                  (reco_outcome[train_valid], gcat[valid_indices][train_valid]), 1)
        np.add.at(out["test_matrix"],
                  (reco_outcome[~train_valid], gcat[valid_indices][~train_valid]), 1)
        chosen = np.flatnonzero(selected)
        if chosen.size:
            delta_bin = np.searchsorted(DELTA_EDGES, rt[chosen] - gt[chosen], side="right") - 1
            keep = (delta_bin >= 0) & (delta_bin < len(DELTA_EDGES) - 1)
            np.add.at(out["delta"], (gcat[chosen][keep], delta_bin[keep]), 1)
            np.add.at(out["topology"], (topo_code[chosen], gcat[chosen]), 1)
            np.add.at(out["topology_matrix"],
                      (topo_code[chosen], rcat[chosen], gcat[chosen]), 1)
    return period, out


def inventory(root: Path, max_files: int) -> Dict[str, Sequence[Path]]:
    found: Dict[str, Sequence[Path]] = {}
    for period, (_, current) in PERIODS.items():
        directory = root / "mc" / period / "dvcs_mc" / current / "parts"
        files = sorted(path for path in directory.glob("*.root")
                       if ".rec_only.root" not in path.name and ".partial" not in path.name)
        if max_files > 0:
            files = files[:max_files]
        if not files:
            raise FileNotFoundError(f"no paired DVCS MC files found in {directory}")
        found[period] = files
    return found


def metrics(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ncat = len(T_LABELS)
    truth = matrix.sum(axis=0).astype(float)
    selected_truth = matrix[:ncat].sum(axis=0).astype(float)
    reconstructed = matrix[:ncat].sum(axis=1).astype(float)
    diagonal = np.diag(matrix[:ncat, :ncat]).astype(float)
    efficiency = np.divide(selected_truth, truth, out=np.zeros_like(truth), where=truth > 0)
    stability = np.divide(diagonal, selected_truth, out=np.zeros_like(truth), where=selected_truth > 0)
    diagonal_efficiency = np.divide(diagonal, truth, out=np.zeros_like(truth), where=truth > 0)
    purity = np.divide(diagonal, reconstructed, out=np.zeros_like(reconstructed), where=reconstructed > 0)
    miss = np.divide(matrix[ncat], truth, out=np.zeros_like(truth), where=truth > 0)
    return efficiency, stability, diagonal_efficiency, purity, miss


def save_csvs(outdir: Path, label: str, result: Mapping[str, np.ndarray]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    matrix = result["matrix"]
    with (outdir / "response_counts.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["reco\\truth", *T_LABELS])
        for row_label, row in zip((*T_LABELS, MISS_LABEL), matrix):
            writer.writerow([row_label, *row.tolist()])
    efficiency, stability, diagonal_efficiency, purity, miss = metrics(matrix)
    with (outdir / "metrics.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["period", "truth_bin", "generated", "selected", "efficiency",
                         "conditional_stability", "diagonal_efficiency", "purity", "miss_fraction"])
        for i, name in enumerate(T_LABELS):
            writer.writerow([label, name, int(matrix[:, i].sum()), int(matrix[:-1, i].sum()),
                             efficiency[i], stability[i], diagonal_efficiency[i], purity[i], miss[i]])
    with (outdir / "cutflow.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["stage", "events"])
        writer.writerows(zip(CUTFLOW_NAMES, result["cutflow"].tolist()))
    with (outdir / "selected_topology_by_truth_bin.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["reconstructed_topology", *T_LABELS])
        for topology_name, row in zip(TOPOLOGIES.values(), result["topology"]):
            writer.writerow([topology_name, *row.tolist()])
    with (outdir / "topology_migration_metrics.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["period", "topology", "truth_bin", "selected_truth", "same_bin",
                         "conditional_stability", "purity"])
        for topology_name, topology_matrix in zip(TOPOLOGIES.values(), result["topology_matrix"]):
            selected_truth = topology_matrix.sum(axis=0).astype(float)
            reconstructed = topology_matrix.sum(axis=1).astype(float)
            diagonal = np.diag(topology_matrix).astype(float)
            topo_stability = np.divide(diagonal, selected_truth, out=np.zeros_like(diagonal),
                                       where=selected_truth > 0)
            topo_purity = np.divide(diagonal, reconstructed, out=np.zeros_like(diagonal),
                                    where=reconstructed > 0)
            for i, name in enumerate(T_LABELS):
                writer.writerow([label, topology_name, name, int(selected_truth[i]), int(diagonal[i]),
                                 topo_stability[i], topo_purity[i]])


def annotate_heatmap(ax, values: np.ndarray, threshold: float = 0.01) -> None:
    for iy in range(values.shape[0]):
        for ix in range(values.shape[1]):
            value = values[iy, ix]
            if np.isfinite(value) and value >= threshold:
                ax.text(ix, iy, f"{100.0 * value:.1f}", ha="center", va="center", fontsize=7,
                        color="white" if value > 0.45 else "black")


def plot_result(outdir: Path, title: str, result: Mapping[str, np.ndarray], dpi: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    matrix = result["matrix"].astype(float)
    ncat = len(T_LABELS)
    truth = matrix.sum(axis=0)
    response = np.divide(matrix, truth[None, :], out=np.zeros_like(matrix), where=truth[None, :] > 0)
    reco = matrix[:ncat].sum(axis=1)
    origin = np.divide(matrix[:ncat], reco[:, None], out=np.zeros((ncat, ncat)), where=reco[:, None] > 0)
    efficiency, stability, diagonal_efficiency, purity, miss = metrics(matrix)

    fig, ax = plt.subplots(figsize=(12.5, 8.0))
    positive = matrix[matrix > 0]
    vmin = max(1.0, positive.min() if positive.size else 1.0)
    norm = LogNorm(vmin=vmin, vmax=max(vmin * 1.01, matrix.max()))
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap="viridis", norm=norm)
    ax.set_xticks(range(ncat), T_LABELS, rotation=45, ha="right")
    ax.set_yticks(range(ncat + 1), (*T_LABELS, MISS_LABEL))
    ax.set_xlabel("generated |t| bin [GeV²]"); ax.set_ylabel("reconstructed outcome")
    ax.set_title(f"{title}: paired-MC response counts")
    fig.colorbar(image, ax=ax, label="events (log scale)")
    fig.tight_layout(); fig.savefig(outdir / "response_counts_log.png", dpi=dpi); plt.close(fig)

    fig, ax = plt.subplots(figsize=(12.5, 8.0))
    image = ax.imshow(response, origin="lower", aspect="auto", cmap="magma", vmin=0.0, vmax=1.0)
    annotate_heatmap(ax, response)
    ax.set_xticks(range(ncat), T_LABELS, rotation=45, ha="right")
    ax.set_yticks(range(ncat + 1), (*T_LABELS, MISS_LABEL))
    ax.set_xlabel("generated |t| bin [GeV²]"); ax.set_ylabel("reconstructed outcome")
    ax.set_title(f"{title}: response probability by truth bin (%)")
    fig.colorbar(image, ax=ax, label="probability")
    fig.tight_layout(); fig.savefig(outdir / "response_probability.png", dpi=dpi); plt.close(fig)

    selected_truth = matrix[:ncat].sum(axis=0)
    conditional_response = np.divide(
        matrix[:ncat], selected_truth[None, :], out=np.zeros((ncat, ncat)),
        where=selected_truth[None, :] > 0,
    )
    fig, ax = plt.subplots(figsize=(12.5, 8.0))
    image = ax.imshow(conditional_response, origin="lower", aspect="auto",
                      cmap="magma", vmin=0.0, vmax=1.0)
    annotate_heatmap(ax, conditional_response)
    ax.set_xticks(range(ncat), T_LABELS, rotation=45, ha="right")
    ax.set_yticks(range(ncat), T_LABELS)
    ax.set_xlabel("generated |t| bin [GeV²]"); ax.set_ylabel("selected reconstructed |t| bin")
    ax.set_title(f"{title}: selected-event conditional response (%)")
    fig.colorbar(image, ax=ax, label="P(reconstructed bin | selected truth bin)")
    fig.tight_layout(); fig.savefig(outdir / "conditional_response_probability.png", dpi=dpi); plt.close(fig)

    fig, ax = plt.subplots(figsize=(12.0, 7.5))
    image = ax.imshow(origin, origin="lower", aspect="auto", cmap="cividis", vmin=0.0, vmax=1.0)
    annotate_heatmap(ax, origin)
    ax.set_xticks(range(ncat), T_LABELS, rotation=45, ha="right")
    ax.set_yticks(range(ncat), T_LABELS)
    ax.set_xlabel("generated |t| bin [GeV²]"); ax.set_ylabel("selected reconstructed |t| bin")
    ax.set_title(f"{title}: truth-bin composition of reconstructed bins (%)")
    fig.colorbar(image, ax=ax, label="fraction of reconstructed bin")
    fig.tight_layout(); fig.savefig(outdir / "reco_origin_composition.png", dpi=dpi); plt.close(fig)

    x = np.arange(ncat)
    fig, ax = plt.subplots(figsize=(12.0, 6.5))
    ax.plot(x, efficiency, "o-", label="efficiency")
    ax.plot(x, stability, "s-", label="conditional stability")
    ax.plot(x, diagonal_efficiency, "d-", label="diagonal efficiency")
    ax.plot(x, purity, "^-", label="purity")
    ax.plot(x, miss, "x--", label="miss fraction")
    ax.set_xticks(x, T_LABELS, rotation=45, ha="right"); ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("fraction"); ax.set_title(f"{title}: response quality metrics")
    ax.grid(alpha=0.3); ax.legend(ncol=2, frameon=False)
    fig.tight_layout(); fig.savefig(outdir / "efficiency_purity_stability.png", dpi=dpi); plt.close(fig)

    delta = result["delta"].astype(float)
    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    image = ax.imshow(delta, origin="lower", aspect="auto", cmap="plasma",
                      extent=[DELTA_EDGES[0], DELTA_EDGES[-1], -0.5, ncat - 0.5],
                      norm=LogNorm(vmin=1.0, vmax=max(1.01, delta.max())))
    ax.set_yticks(range(ncat), T_LABELS); ax.axvline(0.0, color="white", lw=1.0)
    ax.set_xlabel(r"$|t|_{rec}-|t|_{gen}$ [GeV²]"); ax.set_ylabel("generated |t| bin")
    ax.set_title(f"{title}: selected-event |t| residual")
    fig.colorbar(image, ax=ax, label="events (log scale)")
    fig.tight_layout(); fig.savefig(outdir / "t_residual_by_truth_bin.png", dpi=dpi); plt.close(fig)

    flow = result["cutflow"].astype(float)
    fig, ax = plt.subplots(figsize=(10.5, 6.0))
    bars = ax.bar(np.arange(len(flow)), flow, color="#4477aa")
    ax.set_yscale("log"); ax.set_xticks(range(len(flow)), CUTFLOW_NAMES, rotation=30, ha="right")
    ax.set_ylabel("events"); ax.set_title(f"{title}: paired-MC selection cutflow")
    for bar, value in zip(bars, flow):
        ax.text(bar.get_x() + bar.get_width()/2, max(value, 1.0) * 1.08, f"{int(value):,}", ha="center", fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / "selection_cutflow.png", dpi=dpi); plt.close(fig)

    topology = result["topology"].astype(float)
    topology_total = topology.sum(axis=0)
    topology_fraction = np.divide(topology, topology_total[None, :],
                                  out=np.zeros_like(topology), where=topology_total[None, :] > 0)
    fig, axes = plt.subplots(2, 1, figsize=(12.0, 10.0), sharex=True)
    bottom = np.zeros(ncat)
    colors = ("#4477aa", "#ee6677", "#228833")
    for row, name, color in zip(topology, TOPOLOGIES.values(), colors):
        axes[0].bar(x, row, bottom=bottom, label=name, color=color)
        bottom += row
    bottom = np.zeros(ncat)
    for row, name, color in zip(topology_fraction, TOPOLOGIES.values(), colors):
        axes[1].bar(x, row, bottom=bottom, label=name, color=color)
        bottom += row
    axes[0].set_yscale("log"); axes[0].set_ylabel("selected events")
    axes[1].set_ylabel("selected topology fraction"); axes[1].set_ylim(0.0, 1.0)
    axes[0].legend(ncol=3, frameon=False); axes[0].grid(axis="y", alpha=0.25)
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].set_xticks(x, T_LABELS, rotation=45, ha="right")
    axes[1].set_xlabel("generated |t| bin [GeV²]")
    fig.suptitle(f"{title}: selected reconstructed topology versus truth |t|")
    fig.tight_layout(); fig.savefig(outdir / "selected_topology_composition.png", dpi=dpi); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(12.0, 9.0), sharex=True)
    for topology_name, topology_matrix in zip(TOPOLOGIES.values(), result["topology_matrix"]):
        selected_by_truth = topology_matrix.sum(axis=0).astype(float)
        selected_by_reco = topology_matrix.sum(axis=1).astype(float)
        diagonal = np.diag(topology_matrix).astype(float)
        topo_stability = np.divide(diagonal, selected_by_truth, out=np.zeros_like(diagonal),
                                   where=selected_by_truth > 0)
        topo_purity = np.divide(diagonal, selected_by_reco, out=np.zeros_like(diagonal),
                                where=selected_by_reco > 0)
        axes[0].plot(x, topo_stability, marker="o", label=topology_name)
        axes[1].plot(x, topo_purity, marker="o", label=topology_name)
    axes[0].set_ylabel("conditional stability")
    axes[1].set_ylabel("purity")
    for ax in axes:
        ax.set_ylim(0.0, 1.05); ax.grid(alpha=0.3)
    axes[0].legend(ncol=3, frameon=False)
    axes[1].set_xticks(x, T_LABELS, rotation=45, ha="right")
    axes[1].set_xlabel("generated |t| bin [GeV²]")
    fig.suptitle(f"{title}: topology-separated |t| migration")
    fig.tight_layout(); fig.savefig(outdir / "topology_migration_metrics.png", dpi=dpi); plt.close(fig)

    fig, axes = plt.subplots(1, len(TOPOLOGIES), figsize=(18.0, 5.8), sharex=True, sharey=True)
    last_image = None
    for ax, topology_name, topology_matrix in zip(axes, TOPOLOGIES.values(), result["topology_matrix"]):
        selected_by_truth = topology_matrix.sum(axis=0).astype(float)
        conditional = np.divide(
            topology_matrix, selected_by_truth[None, :], out=np.zeros_like(topology_matrix, dtype=float),
            where=selected_by_truth[None, :] > 0,
        )
        last_image = ax.imshow(conditional, origin="lower", aspect="auto",
                               cmap="magma", vmin=0.0, vmax=1.0)
        annotate_heatmap(ax, conditional, threshold=0.02)
        ax.set_xticks(range(ncat), T_LABELS, rotation=55, ha="right")
        ax.set_title(topology_name); ax.set_xlabel("generated |t|")
    axes[0].set_yticks(range(ncat), T_LABELS)
    axes[0].set_ylabel("selected reconstructed |t|")
    fig.colorbar(last_image, ax=axes, label="conditional probability", shrink=0.85)
    fig.suptitle(f"{title}: topology-separated selected-event response")
    fig.subplots_adjust(left=0.07, right=0.91, bottom=0.20, top=0.88, wspace=0.08)
    fig.savefig(outdir / "topology_conditional_response.png", dpi=dpi); plt.close(fig)

    published_indices = (4, 5)  # .11-.15 and .15-.25 reconstructed bins.
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 6.0), sharey=True)
    for ax, reco_index in zip(axes, published_indices):
        values = matrix[reco_index, :ncat]
        fractions = np.divide(values, values.sum(), out=np.zeros_like(values), where=values.sum() > 0)
        ax.bar(x, fractions, color="#66c2a5")
        ax.axvspan(-0.5, 3.5, color="#fc8d62", alpha=0.18, label="generated below 0.11")
        ax.set_title(f"reconstructed {T_LABELS[reco_index]} GeV²")
        ax.set_xticks(x, T_LABELS, rotation=55, ha="right"); ax.grid(axis="y", alpha=0.25)
        ax.set_ylim(0.0, max(1.0, 1.1 * fractions.max()))
    axes[0].set_ylabel("fraction by generated |t| origin")
    axes[0].legend(frameon=False)
    fig.suptitle(f"{title}: origin of the first two published reconstructed |t| bins")
    fig.tight_layout(); fig.savefig(outdir / "low_t_feed_in_origins.png", dpi=dpi); plt.close(fig)


def iterative_bayes_unfold(response: np.ndarray, observed: np.ndarray,
                           prior_counts: np.ndarray, iterations: int = 4) -> np.ndarray:
    """D'Agostini unfolding for a response that includes inefficiency."""
    efficiency = response.sum(axis=0)
    prior = np.asarray(prior_counts, dtype=float).copy()
    if prior.sum() <= 0.0:
        return np.zeros_like(prior)
    prior /= prior.sum()
    unfolded = np.zeros_like(prior)
    for _ in range(iterations):
        denominator = response @ prior
        posterior = np.divide(
            response * prior[None, :], denominator[:, None],
            out=np.zeros_like(response), where=denominator[:, None] > 0,
        )
        unfolded = np.divide(
            (observed[:, None] * posterior).sum(axis=0), efficiency,
            out=np.zeros_like(prior), where=efficiency > 0,
        )
        if unfolded.sum() > 0.0:
            prior = unfolded / unfolded.sum()
    return unfolded


PUBLISHED = np.arange(4, 10)
PUBLISHED_LABELS = T_LABELS[4:10]


def default_reweight_scenarios() -> Dict[str, np.ndarray]:
    return {
        "nominal": np.ones(6),
        "moderate_low_t_enhancement": np.asarray([1.30, 1.16, 1.07, 1.00, 0.94, 0.90]),
        "strong_low_t_enhancement": np.asarray([1.80, 1.40, 1.10, 0.90, 0.75, 0.65]),
        "moderate_low_t_suppression": np.asarray([0.72, 0.86, 0.95, 1.00, 1.05, 1.10]),
    }


def load_reweight_scenarios(path: Path | None) -> Dict[str, np.ndarray]:
    scenarios = default_reweight_scenarios()
    if path is None:
        return scenarios
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError("reweight JSON must be an object mapping scenario names to six weights")
    for name, values in payload.items():
        weights = np.asarray(values, dtype=float)
        if weights.shape != (6,) or not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
            raise ValueError(f"reweight scenario {name!r} must contain six finite positive weights")
        scenarios[str(name)] = weights
    return scenarios


def closure_inputs(result: Mapping[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    return (result["train_matrix"].astype(float), result["test_matrix"].astype(float))


def response_from_train(train: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    train_truth = train[:, PUBLISHED].sum(axis=0)
    selected = train[PUBLISHED, :][:, PUBLISHED]
    response = np.divide(
        selected, train_truth[None, :], out=np.zeros((6, 6)), where=train_truth[None, :] > 0,
    )
    return train_truth, selected, response


def bin_by_bin_correct(train: np.ndarray, observed: np.ndarray,
                       reconstructed_reweight_iterations: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    """Correct reco counts; optional weights are inferred only from reco-level ratios."""
    train_truth, selected, _ = response_from_train(train)
    weights = np.ones(6)
    for _ in range(reconstructed_reweight_iterations):
        predicted = selected @ weights
        scale = np.divide(predicted.sum(), observed.sum(), out=np.asarray(1.0), where=observed.sum() > 0)
        ratio = np.divide(observed * scale, predicted, out=np.ones_like(observed), where=predicted > 0)
        # Damping prevents an alternating-bin correction from overshooting.
        weights *= np.clip(ratio, 0.25, 4.0) ** 0.7
        weights /= np.exp(np.mean(np.log(weights)))
    weighted_reco = selected @ weights
    weighted_truth = train_truth * weights
    acceptance = np.divide(weighted_reco, weighted_truth, out=np.zeros_like(weighted_truth),
                           where=weighted_truth > 0)
    corrected = np.divide(observed, acceptance, out=np.zeros_like(observed), where=acceptance > 0)
    return corrected, weights


def evaluate_closure(train: np.ndarray, test: np.ndarray, weights: np.ndarray,
                     bayes_iterations: int = 4, prior: np.ndarray | None = None,
                     acceptance_reweight_iterations: int = 3) -> Dict[str, np.ndarray]:
    train_truth, _, response = response_from_train(train)
    weighted_test = test[PUBLISHED, :][:, PUBLISHED] * weights[None, :]
    truth = test[:, PUBLISHED].sum(axis=0) * weights
    observed = weighted_test.sum(axis=1)
    baseline, _ = bin_by_bin_correct(train, observed, 0)
    reweighted, inferred_weights = bin_by_bin_correct(
        train, observed, acceptance_reweight_iterations,
    )
    unfolded = iterative_bayes_unfold(
        response, observed, train_truth if prior is None else prior, iterations=bayes_iterations,
    )
    return {
        "truth": truth, "observed": observed, "bin_by_bin": baseline,
        "reweighted_acceptance": reweighted, "unfolded": unfolded,
        "inferred_weights": inferred_weights,
    }


def ratio_to_truth(values: np.ndarray, truth: np.ndarray) -> np.ndarray:
    return np.divide(values, truth, out=np.zeros_like(values), where=truth > 0)


def closure_toys(train: np.ndarray, test: np.ndarray, weights: np.ndarray, toys: int,
                 seed: int) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    samples = {name: [] for name in ("bin_by_bin", "reweighted_acceptance", "unfolded")}
    if toys <= 1:
        zeros = {name: np.zeros(6) for name in samples}
        covariances = {name: np.zeros((6, 6)) for name in samples}
        return zeros, covariances
    rng = np.random.default_rng(seed)
    for _ in range(toys):
        toy_train = rng.poisson(train).astype(float)
        toy_test = rng.poisson(test).astype(float)
        evaluated = evaluate_closure(toy_train, toy_test, weights)
        for name in samples:
            samples[name].append(ratio_to_truth(evaluated[name], evaluated["truth"]))
    arrays = {name: np.asarray(values) for name, values in samples.items()}
    errors = {name: values.std(axis=0, ddof=1) for name, values in arrays.items()}
    covariances = {name: np.cov(values, rowvar=False, ddof=1) for name, values in arrays.items()}
    return errors, covariances


def save_closure(outdir: Path, title: str, result: Mapping[str, np.ndarray], dpi: int,
                 scenarios: Mapping[str, np.ndarray], toys: int, seed: int) -> None:
    train, test = closure_inputs(result)
    evaluated = {name: evaluate_closure(train, test, weights)
                 for name, weights in scenarios.items()}
    uncertainties = {}
    covariances = {}
    for index, (name, weights) in enumerate(scenarios.items()):
        uncertainties[name], covariances[name] = closure_toys(
            train, test, weights, toys, seed + 1009 * index,
        )

    method_names = ("bin_by_bin", "reweighted_acceptance", "unfolded")
    with (outdir / "closure_metrics.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["period", "scenario", "truth_bin", "truth", "observed",
                         "method", "corrected", "corrected_over_truth", "toy_ratio_uncertainty"])
        for scenario, values in evaluated.items():
            for method in method_names:
                ratios = ratio_to_truth(values[method], values["truth"])
                for i, bin_label in enumerate(PUBLISHED_LABELS):
                    writer.writerow([title, scenario, bin_label, values["truth"][i],
                                     values["observed"][i], method, values[method][i], ratios[i],
                                     uncertainties[scenario][method][i]])

    x = np.arange(6)
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 10.5), sharex=True, sharey=True)
    styles = {
        "bin_by_bin": ("o-", "bin-by-bin acceptance"),
        "reweighted_acceptance": ("^-", "reweighted acceptance (3 iter.)"),
        "unfolded": ("s-", "Bayes unfolding (4 iter.)"),
    }
    for ax, (scenario, values) in zip(axes.ravel(), evaluated.items()):
        ax.axhline(1.0, color="black", lw=1.0)
        for method in method_names:
            style, label = styles[method]
            ax.errorbar(x, ratio_to_truth(values[method], values["truth"]),
                        yerr=uncertainties[scenario][method], fmt=style, capsize=2, label=label)
        ax.set_title(scenario.replace("_", " "))
        ax.grid(alpha=0.3); ax.set_ylim(0.94, 1.06)
        ax.set_xticks(x, PUBLISHED_LABELS, rotation=45, ha="right")
    axes[0, 0].set_ylabel("corrected / known truth")
    axes[1, 0].set_ylabel("corrected / known truth")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"{title}: independent half-sample closure with toy uncertainties")
    fig.tight_layout(); fig.savefig(outdir / "closure_comparison.png", dpi=dpi); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    for name, weights in scenarios.items():
        ax.plot(x, weights, marker="o", label=name.replace("_", " "))
    ax.axhline(1.0, color="black", lw=1.0); ax.grid(alpha=0.3)
    ax.set_xticks(x, PUBLISHED_LABELS, rotation=45, ha="right")
    ax.set_ylabel("truth-spectrum weight"); ax.set_title(f"{title}: closure spectra being tested")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(outdir / "closure_tested_spectra.png", dpi=dpi); plt.close(fig)

    stress_name = "strong_low_t_enhancement"
    if stress_name not in scenarios:
        stress_name = next(reversed(scenarios))
    stress_weights = scenarios[stress_name]
    bayes_iterations = (1, 2, 3, 4, 6, 8, 10)
    acceptance_iterations = tuple(range(0, 7))
    bayes_bias = []
    for count in bayes_iterations:
        values = evaluate_closure(train, test, stress_weights, bayes_iterations=count)
        bayes_bias.append(100.0 * (ratio_to_truth(values["unfolded"], values["truth"]) - 1.0))
    acceptance_bias = []
    for count in acceptance_iterations:
        values = evaluate_closure(train, test, stress_weights,
                                  acceptance_reweight_iterations=count)
        acceptance_bias.append(100.0 * (ratio_to_truth(values["reweighted_acceptance"], values["truth"]) - 1.0))
    bayes_bias = np.asarray(bayes_bias)
    acceptance_bias = np.asarray(acceptance_bias)
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.5))
    axes[0].plot(bayes_iterations, bayes_bias[:, 0], "o-", label="first-bin bias")
    axes[0].plot(bayes_iterations, np.max(np.abs(bayes_bias), axis=1), "s-", label="largest |bias|")
    axes[0].set_xlabel("Bayesian unfolding iterations")
    axes[1].plot(acceptance_iterations, acceptance_bias[:, 0], "o-", label="first-bin bias")
    axes[1].plot(acceptance_iterations, np.max(np.abs(acceptance_bias), axis=1), "s-", label="largest |bias|")
    axes[1].set_xlabel("reconstructed-ratio reweighting iterations")
    for ax in axes:
        ax.axhline(0.0, color="black", lw=1.0); ax.set_ylabel("closure bias [%]")
        ax.grid(alpha=0.3); ax.legend(frameon=False)
    fig.suptitle(f"{title}: regularization/convergence scan ({stress_name.replace('_', ' ')})")
    fig.tight_layout(); fig.savefig(outdir / "closure_iteration_scan.png", dpi=dpi); plt.close(fig)

    train_truth, _, _ = response_from_train(train)
    priors = {
        "MC prior": train_truth,
        "flat prior": np.ones(6),
        "low-t enhanced prior": train_truth * np.asarray([1.8, 1.4, 1.1, .9, .75, .65]),
        "low-t suppressed prior": train_truth * np.asarray([.7, .85, .95, 1., 1.05, 1.1]),
    }
    fig, ax = plt.subplots(figsize=(10.5, 6.0))
    for prior_name, prior in priors.items():
        first_bin_bias = []
        for count in bayes_iterations:
            values = evaluate_closure(train, test, stress_weights, bayes_iterations=count, prior=prior)
            first_bin_bias.append(100.0 * (values["unfolded"][0] / values["truth"][0] - 1.0))
        ax.plot(bayes_iterations, first_bin_bias, marker="o", label=prior_name)
    ax.axhline(0.0, color="black", lw=1.0); ax.grid(alpha=0.3)
    ax.set_xlabel("Bayesian unfolding iterations"); ax.set_ylabel("first-bin closure bias [%]")
    ax.set_title(f"{title}: unfolding prior dependence under strong spectrum stress")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(outdir / "closure_prior_scan.png", dpi=dpi); plt.close(fig)

    for method in method_names:
        covariance = covariances[stress_name][method]
        sigma = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
        correlation = np.divide(covariance, sigma[:, None] * sigma[None, :],
                                out=np.zeros_like(covariance),
                                where=(sigma[:, None] * sigma[None, :]) > 0)
        fig, ax = plt.subplots(figsize=(7.2, 6.2))
        image = ax.imshow(correlation, origin="lower", cmap="coolwarm", vmin=-1.0, vmax=1.0)
        for iy in range(6):
            for ix in range(6):
                ax.text(ix, iy, f"{correlation[iy, ix]:.2f}", ha="center", va="center", fontsize=8)
        ax.set_xticks(range(6), PUBLISHED_LABELS, rotation=45, ha="right")
        ax.set_yticks(range(6), PUBLISHED_LABELS)
        ax.set_title(f"{title}: {method.replace('_', ' ')} toy correlation")
        fig.colorbar(image, ax=ax, label="correlation")
        fig.tight_layout(); fig.savefig(outdir / f"closure_correlation_{method}.png", dpi=dpi); plt.close(fig)


def period_comparison(outdir: Path, results: Mapping[str, Mapping[str, np.ndarray]], dpi: int) -> None:
    ncat = len(T_LABELS)
    fig, axes = plt.subplots(3, 1, figsize=(13.0, 14.0), sharex=True)
    for period, result in results.items():
        efficiency, stability, _, purity, _ = metrics(result["matrix"])
        label = PERIODS[period][0]
        axes[0].plot(range(ncat), efficiency, marker="o", label=label)
        axes[1].plot(range(ncat), purity, marker="o", label=label)
        axes[2].plot(range(ncat), stability, marker="o", label=label)
    for ax, name in zip(axes, ("efficiency", "purity", "stability")):
        ax.set_ylim(0.0, 1.05); ax.set_ylabel(name); ax.grid(alpha=0.3)
    axes[0].legend(ncol=3, frameon=False)
    axes[-1].set_xticks(range(ncat), T_LABELS, rotation=45, ha="right")
    axes[-1].set_xlabel("generated |t| bin [GeV²]")
    fig.suptitle("Period comparison of paired-MC |t| response diagnostics")
    fig.tight_layout(); fig.savefig(outdir / "period_comparison_metrics.png", dpi=dpi); plt.close(fig)


def write_root(path: Path, results: Mapping[str, Mapping[str, np.ndarray]]) -> None:
    with uproot.recreate(path) as root_file:
        for period, result in results.items():
            matrix = result["matrix"]
            rr, tt = np.indices(matrix.shape)
            root_file[f"{period}/response_cells"] = {
                "reco_category": rr.ravel().astype(np.int32),
                "truth_category": tt.ravel().astype(np.int32),
                "count": matrix.ravel().astype(np.int64),
            }
            efficiency, stability, diagonal_efficiency, purity, miss = metrics(matrix)
            root_file[f"{period}/metrics"] = {
                "truth_category": np.arange(len(T_LABELS), dtype=np.int32),
                "efficiency": efficiency, "conditional_stability": stability,
                "diagonal_efficiency": diagonal_efficiency,
                "purity": purity, "miss_fraction": miss,
            }


def write_plot_guide(path: Path) -> None:
    path.write_text(
        """PAIRED-MC |t| MIGRATION STUDY: PLOT GUIDE

Decision plots
--------------
closure_tested_spectra.png
  Shows the exact truth-spectrum distortions used in closure. These are tests,
  not fitted statements about the data.

closure_comparison.png
  Corrected result divided by known independent test truth. One is perfect.
  Error bars are the spread of Poisson response/test toys. Compare ordinary
  acceptance, reconstructed-ratio reweighted acceptance, and Bayes unfolding.

closure_iteration_scan.png
  First-bin and worst-bin bias versus unfolding/reweighting iteration count.
  Prefer the smallest iteration count on the stable bias plateau; continued
  movement or growth indicates prior dependence or over-iteration.

closure_prior_scan.png
  First-bin unfolding bias for deliberately different starting priors. Curves
  should converge before statistical fluctuations become excessive.

closure_correlation_*.png
  Toy correlation of corrected/truth ratios. Off-diagonal structure is the
  statistical coupling introduced by each correction method.

Response and topology plots
---------------------------
conditional_response_probability.png
  Migration among selected events only; columns sum to one. This isolates
  resolution from inefficiency.

response_probability.png
  Full response including miss/not-selected. This shows absolute efficiency,
  but the selected diagonal is visually compressed by the large miss fraction.

reco_origin_composition.png
  For each reconstructed bin, where its events were generated (purity view).

topology_conditional_response.png / topology_migration_metrics.png
  Selected-event migration separated into FD_FD, CD_FD, and CD_FT. Empty or
  tiny topology bins must not be interpreted as precise measurements.

low_t_feed_in_origins.png
  Direct test of whether events generated below the published boundary feed
  the first two reconstructed bins.

Machine-readable outputs
------------------------
closure_metrics.csv contains method values, corrected/truth ratios, and toy
uncertainties. metrics.csv separates efficiency, conditional stability,
diagonal efficiency, purity, and miss fraction.
"""
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path("/volatile/clas12/thayward/dvcs_pass2_full_rebuild"))
    parser.add_argument("--global-cuts", type=Path, default=Path("output/jsons/global_cuts_config.json"))
    parser.add_argument("--exclusivity-cuts", type=Path, default=Path("output/jsons/combined_cuts.json"))
    parser.add_argument("--output", type=Path, default=Path("output/migration_study/t_response"))
    parser.add_argument("--reuse-response-dir", type=Path,
                        help="Reuse per-period response_arrays.npz files; skip all ROOT event processing")
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--step-size", default="200 MB")
    parser.add_argument("--max-files-per-period", type=int, default=0,
                        help="Development/test limit; zero processes all files")
    parser.add_argument("--reweight-json", type=Path,
                        help="Optional JSON object adding named six-bin truth-weight scenarios")
    parser.add_argument("--closure-toys", type=int, default=300,
                        help="Poisson toys used for closure uncertainties/correlations; zero disables")
    parser.add_argument("--closure-seed", type=int, default=20261102)
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    if uproot is None and args.reuse_response_dir is None:
        parser.error(
            "uproot is required for paired-file streaming. Run in the same "
            "Python environment used by plot_exclusivity_data_dvcs_pi0_mc.py."
        )
    if not 1 <= args.workers <= 10:
        parser.error("--workers must be in [1,10]")
    if args.closure_toys < 0:
        parser.error("--closure-toys must be non-negative")
    scenarios = load_reweight_scenarios(args.reweight_json)
    files: Dict[str, Sequence[Path]] = {}
    if args.reuse_response_dir is not None:
        results = {}
        for period in PERIODS:
            source = args.reuse_response_dir / period / "response_arrays.npz"
            if not source.is_file():
                parser.error(f"missing reusable response array: {source}")
            with np.load(source) as payload:
                result = {key: payload[key] for key in payload.files}
            required = set(empty_result())
            missing = sorted(required - set(result))
            if missing:
                parser.error(f"{source} predates required arrays: {missing}; rerun event processing once")
            results[period] = result
        combined = empty_result()
        for result in results.values():
            add_result(combined, result)
        print(f"Reusing response arrays from {args.reuse_response_dir}; no ROOT files will be read", flush=True)
    else:
        cfg = json.loads(args.global_cuts.read_text())
        cuts = json.loads(args.exclusivity_cuts.read_text())
        files = inventory(args.input_root, args.max_files_per_period)
        tasks = [(period, str(path), args.step_size, cfg, cuts, "")
                 for period, paths in files.items() for path in paths]
        print(f"Processing {len(tasks)} paired files with up to {args.workers} workers", flush=True)
        results = {period: empty_result() for period in PERIODS}
        with futures.ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as pool:
            pending = {pool.submit(process_file, task): task for task in tasks}
            completed = 0
            for future in futures.as_completed(pending):
                period, result = future.result()
                add_result(results[period], result)
                completed += 1
                print(f"[{completed}/{len(tasks)}] {period}: {Path(pending[future][1]).name}", flush=True)
        combined = empty_result()
        for result in results.values():
            add_result(combined, result)
    all_results = dict(results)
    all_results["combined"] = combined
    args.output.mkdir(parents=True, exist_ok=True)
    for result_index, (period, result) in enumerate(all_results.items()):
        label = "All periods combined" if period == "combined" else PERIODS[period][0]
        destination = args.output / period
        save_csvs(destination, label, result)
        plot_result(destination, label, result, args.dpi)
        save_closure(destination, label, result, args.dpi, scenarios,
                     args.closure_toys, args.closure_seed + 7919 * result_index)
        np.savez_compressed(destination / "response_arrays.npz", **result)
    period_comparison(args.output, results, args.dpi)
    if uproot is not None:
        write_root(args.output / "migration_response.root", all_results)
    else:
        print("uproot unavailable in reuse mode; skipped migration_response.root", flush=True)
    write_plot_guide(args.output / "PLOT_GUIDE.txt")
    manifest = {
        "study": "paired-MC |t| response, reweighting, and unfolding closure",
        "input_root": str(args.input_root), "global_cuts": str(args.global_cuts),
        "exclusivity_cuts": str(args.exclusivity_cuts), "workers": args.workers,
        "step_size": args.step_size, "t_edges_GeV2": T_EDGES.tolist(),
        "closure_toys": args.closure_toys, "closure_seed": args.closure_seed,
        "reweight_scenarios": {name: weights.tolist() for name, weights in scenarios.items()},
        "reused_response_dir": str(args.reuse_response_dir) if args.reuse_response_dir else None,
        "files": {period: [str(path) for path in paths] for period, paths in files.items()},
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote migration study to {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
