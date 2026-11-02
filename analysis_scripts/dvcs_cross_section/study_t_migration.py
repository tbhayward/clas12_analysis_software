#!/usr/bin/env python3
"""First-stage paired-MC |t| migration and response feasibility study.

This diagnostic is deliberately separate from the production correction.  It
uses the generated and reconstructed records in the same MCEvents entry,
applies the production reconstructed-MC selection, and reports whether |t|
migration is large enough to motivate matrix unfolding.
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


def closure_scenarios(result: Mapping[str, np.ndarray]) -> Sequence[Tuple[str, np.ndarray, np.ndarray, np.ndarray]]:
    """Return truth, bin-by-bin, and unfolded spectra for two closure stresses."""
    # Limit closure to the six published bins. Feed-in from below 0.11 and
    # above 1.0 is reported separately by the response diagnostics.
    published = np.arange(4, 10)
    train = result["train_matrix"].astype(float)
    test = result["test_matrix"].astype(float)
    train_truth = train[:, published].sum(axis=0)
    train_reco = train[published, :][:, published].sum(axis=1)
    response = np.divide(
        train[published, :][:, published], train_truth[None, :],
        out=np.zeros((published.size, published.size)), where=train_truth[None, :] > 0,
    )
    bin_by_bin_acceptance = np.divide(
        train_reco, train_truth, out=np.zeros_like(train_truth), where=train_truth > 0,
    )
    scenarios = (
        ("nominal", np.ones(published.size)),
        # Deliberately unlike the training prior: enhance the first bins and
        # progressively suppress the high-|t| tail.
        ("low_t_enhanced_stress", np.asarray([1.80, 1.40, 1.10, 0.90, 0.75, 0.65])),
    )
    output = []
    for name, weights in scenarios:
        weighted_test = test[published, :][:, published] * weights[None, :]
        truth = test[:, published].sum(axis=0) * weights
        observed = weighted_test.sum(axis=1)
        corrected = np.divide(
            observed, bin_by_bin_acceptance, out=np.zeros_like(observed),
            where=bin_by_bin_acceptance > 0,
        )
        unfolded = iterative_bayes_unfold(response, observed, train_truth, iterations=4)
        output.append((name, truth, corrected, unfolded))
    return output


def save_closure(outdir: Path, title: str, result: Mapping[str, np.ndarray], dpi: int) -> None:
    published_labels = T_LABELS[4:10]
    scenarios = closure_scenarios(result)
    with (outdir / "closure_metrics.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["period", "scenario", "truth_bin", "truth", "bin_by_bin_corrected",
                         "iterative_bayes_unfolded", "bin_by_bin_over_truth", "unfolded_over_truth"])
        for scenario, truth, corrected, unfolded in scenarios:
            corrected_ratio = np.divide(corrected, truth, out=np.zeros_like(truth), where=truth > 0)
            unfolded_ratio = np.divide(unfolded, truth, out=np.zeros_like(truth), where=truth > 0)
            for i, bin_label in enumerate(published_labels):
                writer.writerow([title, scenario, bin_label, truth[i], corrected[i], unfolded[i],
                                 corrected_ratio[i], unfolded_ratio[i]])

    fig, axes = plt.subplots(1, len(scenarios), figsize=(14.0, 5.5), sharey=True)
    x = np.arange(len(published_labels))
    for ax, (scenario, truth, corrected, unfolded) in zip(np.atleast_1d(axes), scenarios):
        corrected_ratio = np.divide(corrected, truth, out=np.zeros_like(truth), where=truth > 0)
        unfolded_ratio = np.divide(unfolded, truth, out=np.zeros_like(truth), where=truth > 0)
        ax.axhline(1.0, color="black", lw=1.0)
        ax.plot(x, corrected_ratio, "o-", label="bin-by-bin acceptance")
        ax.plot(x, unfolded_ratio, "s-", label="Bayes unfolding (4 iter.)")
        ax.set_xticks(x, published_labels, rotation=45, ha="right")
        ax.set_ylim(0.75, 1.25); ax.grid(alpha=0.3)
        ax.set_title(scenario.replace("_", " "))
    axes[0].set_ylabel("corrected / known test truth")
    axes[0].legend(frameon=False)
    fig.suptitle(f"{title}: independent half-sample closure")
    fig.tight_layout(); fig.savefig(outdir / "closure_comparison.png", dpi=dpi); plt.close(fig)


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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path("/volatile/clas12/thayward/dvcs_pass2_full_rebuild"))
    parser.add_argument("--global-cuts", type=Path, default=Path("output/jsons/global_cuts_config.json"))
    parser.add_argument("--exclusivity-cuts", type=Path, default=Path("output/jsons/combined_cuts.json"))
    parser.add_argument("--output", type=Path, default=Path("output/migration_study/t_response"))
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--step-size", default="200 MB")
    parser.add_argument("--max-files-per-period", type=int, default=0,
                        help="Development/test limit; zero processes all files")
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    if uproot is None:
        parser.error(
            "uproot is required for paired-file streaming. Run in the same "
            "Python environment used by plot_exclusivity_data_dvcs_pi0_mc.py."
        )
    if not 1 <= args.workers <= 10:
        parser.error("--workers must be in [1,10]")
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
    for period, result in all_results.items():
        label = "All periods combined" if period == "combined" else PERIODS[period][0]
        destination = args.output / period
        save_csvs(destination, label, result)
        plot_result(destination, label, result, args.dpi)
        save_closure(destination, label, result, args.dpi)
        np.savez_compressed(destination / "response_arrays.npz", **result)
    period_comparison(args.output, results, args.dpi)
    write_root(args.output / "migration_response.root", all_results)
    manifest = {
        "study": "paired-MC |t| response feasibility, milestone 1",
        "input_root": str(args.input_root), "global_cuts": str(args.global_cuts),
        "exclusivity_cuts": str(args.exclusivity_cuts), "workers": args.workers,
        "step_size": args.step_size, "t_edges_GeV2": T_EDGES.tolist(),
        "files": {period: [str(path) for path in paths] for period, paths in files.items()},
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote migration study to {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
