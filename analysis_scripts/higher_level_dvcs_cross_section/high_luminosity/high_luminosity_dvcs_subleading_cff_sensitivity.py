#!/usr/bin/env python3
"""
High-luminosity CLAS12 DVCS: local subleading-CFF separation study.

Purpose
-------
Use the same matched Pass-2 XS+BSA pseudo-data inputs as the D-term study, but
ask a different question: how well can an unpolarized proton target separate
H, Htilde, and E as luminosity increases?

No dispersion relation and no D-term constraint are used here.  In each
(xB,Q2,t) cell the selected real and imaginary CFF components are varied
independently around the KM15 pseudo-truth and propagated through Gepard to
the measured phi-dependent unpolarized cross section and beam-spin asymmetry.

Fit modes
---------
  H              : ReH, ImH
  H_Ht           : ReH, ImH, ReHt, ImHt
  H_Ht_E         : ReH, ImH, ReHt, ImHt, ReE, ImE

The comparison H_Ht -> H_Ht_E directly quantifies how much apparent Htilde
precision is lost when E is no longer held fixed.

Experimental treatment
----------------------
  * same Stage-2 Pass-2 kinematics and pseudo-data uncertainties
  * XS and BSA used jointly in every fit
  * one shared XS normalization nuisance
  * one shared BSA beam-polarization nuisance
  * statistical uncertainties scale with luminosity in the Stage-2 inputs
  * point-to-point and scale systematics remain as supplied by Stage 2

Outputs
-------
  output/stage4_cff_separation/
      tables/
      figures/

Important
---------
Relative errors can become meaningless near a KM15 zero crossing.  Absolute
uncertainties are always saved.  Relative uncertainties are saved only when
|CFF| exceeds --relative-floor.

This is a local sensitivity/separation projection, not a model-independent
global CFF extraction.
"""

from __future__ import annotations

import argparse
import contextlib
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)

FIT_MODES = {
    "H": ("ReH", "ImH"),
    "H_Ht": ("ReH", "ImH", "ReHt", "ImHt"),
    "H_Ht_E": ("ReH", "ImH", "ReHt", "ImHt", "ReE", "ImE"),
}

ALL_CFFS = ("ReH", "ImH", "ReHt", "ImHt", "ReE", "ImE")

LABELS = {
    "ReH": r"$\mathrm{Re}\,\mathcal{H}$",
    "ImH": r"$\mathrm{Im}\,\mathcal{H}$",
    "ReHt": r"$\mathrm{Re}\,\widetilde{\mathcal{H}}$",
    "ImHt": r"$\mathrm{Im}\,\widetilde{\mathcal{H}}$",
    "ReE": r"$\mathrm{Re}\,\mathcal{E}$",
    "ImE": r"$\mathrm{Im}\,\mathcal{E}$",
}

DEFAULT_MIN_NPHI = 3
DEFAULT_REL_STEP = 0.02
DEFAULT_ABS_STEP = 0.05
DEFAULT_RELATIVE_FLOOR = 0.10
SVD_RTOL = 1.0e-11


def make_point(g, row, helicity: int = 0):
    phi_trento = math.pi - math.radians(float(row.phi_deg))
    pt = g.DataPoint(
        xB=float(row.xB),
        t=-abs(float(row.t_abs)),
        Q2=float(row.Q2),
        phi=float(phi_trento),
        observable="XS",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(row.ebeam),
        in1charge=-1,
        in1polarization=int(helicity),
        in2particle="p",
    )
    pt.prepare()
    return pt


def cff_owner(th, name: str):
    candidates = [th.m] if hasattr(th, "m") else []
    candidates.append(th)
    seen = set()
    for obj in candidates:
        if obj is None or id(obj) in seen:
            continue
        seen.add(id(obj))
        if hasattr(obj, name) and callable(getattr(obj, name)):
            return obj
    raise AttributeError(f"Could not find callable Gepard CFF method {name}")


@contextlib.contextmanager
def shifted_cff(th, name: str, shift: float):
    owner = cff_owner(th, name)
    old = getattr(owner, name)
    setattr(
        owner,
        name,
        lambda pt, _old=old, _shift=float(shift): float(_old(pt)) + _shift,
    )
    try:
        yield
    finally:
        setattr(owner, name, old)


def pred(th, pt) -> float:
    return float(th.predict(pt))


def cff_value(th, name: str, pt) -> float:
    return float(getattr(cff_owner(th, name), name)(pt))


def bsa_from_points(th, plus_pt, minus_pt) -> float:
    sp = pred(th, plus_pt)
    sm = pred(th, minus_pt)
    den = sp + sm
    if not np.isfinite(den) or abs(den) < 1e-30:
        return np.nan
    return (sp - sm) / den


def finite_derivative(
    th,
    points: Sequence,
    observable: str,
    name: str,
    nominal_cff: float,
    rel_step: float,
    abs_step: float,
) -> Tuple[np.ndarray, float]:
    step = max(abs_step, rel_step * max(abs(nominal_cff), 1.0))

    def evaluate():
        if observable == "xs":
            return np.asarray([pred(th, p) for p in points], dtype=float)
        return np.asarray(
            [bsa_from_points(th, pair[0], pair[1]) for pair in points],
            dtype=float,
        )

    with shifted_cff(th, name, +step):
        plus = evaluate()
    with shifted_cff(th, name, -step):
        minus = evaluate()

    return (plus - minus) / (2.0 * step), step


def load_inputs(indir: Path) -> Dict[int, pd.DataFrame]:
    out = {}
    for factor in LUMI_FACTORS:
        path = indir / f"joint_fit_input_km15_{factor}x.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        d = pd.read_csv(path)
        d["bin"] = d["bin"].astype(int)
        out[factor] = d.sort_values(["bin", "phi_deg"]).reset_index(drop=True)

    ref = out[1][["point_id", "bin"]].astype(str)
    for factor in LUMI_FACTORS[1:]:
        if not ref.equals(out[factor][["point_id", "bin"]].astype(str)):
            raise RuntimeError("Luminosity input kinematics/order mismatch.")
    return out


def build_derivatives(df, th, g, rel_step, abs_step, min_nphi):
    point_rows, cell_rows = [], []
    groups = list(df.groupby("bin", sort=True))
    print(f"[CFF separation] calculating derivatives for {len(groups)} cells")

    for icell, (bid, cell) in enumerate(groups, 1):
        cell = cell.sort_values("phi_deg")
        if len(cell) < min_nphi:
            print(f"[skip] bin {bid}: only {len(cell)} phi points")
            continue

        xs_pts = [make_point(g, r, 0) for r in cell.itertuples(index=False)]
        bsa_pts = [
            (make_point(g, r, +1), make_point(g, r, -1))
            for r in cell.itertuples(index=False)
        ]
        pt0 = xs_pts[0]

        cff0 = {name: cff_value(th, name, pt0) for name in ALL_CFFS}
        xs0 = np.asarray([pred(th, p) for p in xs_pts], dtype=float)
        a0 = np.asarray([bsa_from_points(th, p, m) for p, m in bsa_pts], dtype=float)

        dx, da, steps = {}, {}, {}
        for name in ALL_CFFS:
            dx[name], steps[name] = finite_derivative(
                th, xs_pts, "xs", name, cff0[name], rel_step, abs_step
            )
            da[name], _ = finite_derivative(
                th, bsa_pts, "bsa", name, cff0[name], rel_step, abs_step
            )

        # A transparent observable-level diagnostic: response to a +10% change
        # in Im(Htilde).  This is not a fit result.
        delta_imht_10 = 0.10 * cff0["ImHt"]
        xs_frac_response_10 = np.full(len(cell), np.nan)
        good_xs = np.isfinite(xs0) & (np.abs(xs0) > 1e-30)
        xs_frac_response_10[good_xs] = (
            dx["ImHt"][good_xs] * delta_imht_10 / xs0[good_xs]
        )
        bsa_abs_response_10 = da["ImHt"] * delta_imht_10

        for j, row in enumerate(cell.itertuples(index=False)):
            rec = {
                "point_id": row.point_id,
                "bin": int(bid),
                "phi_deg": float(row.phi_deg),
                "xs_km15": float(xs0[j]),
                "bsa_km15": float(a0[j]),
                "xs_frac_response_to_10pct_ImHt": float(xs_frac_response_10[j]),
                "bsa_abs_response_to_10pct_ImHt": float(bsa_abs_response_10[j]),
            }
            for name in ALL_CFFS:
                rec[f"d_xs_d_{name}"] = float(dx[name][j])
                rec[f"d_bsa_d_{name}"] = float(da[name][j])
            point_rows.append(rec)

        xi = float(np.median(cell["xB"])) / (2.0 - float(np.median(cell["xB"])))
        meta = {
            "bin": int(bid),
            "xB": float(np.median(cell["xB"])),
            "xi": xi,
            "Q2": float(np.median(cell["Q2"])),
            "t_abs": float(np.median(cell["t_abs"])),
            "t_over_Q2": float(np.median(cell["t_over_Q2"])),
            "n_phi": len(cell),
            "median_abs_xs_frac_response_to_10pct_ImHt": float(
                np.nanmedian(np.abs(xs_frac_response_10))
            ),
            "median_abs_bsa_response_to_10pct_ImHt": float(
                np.nanmedian(np.abs(bsa_abs_response_10))
            ),
        }
        for name in ALL_CFFS:
            meta[f"{name}_KM15"] = cff0[name]
            meta[f"finite_step_{name}"] = steps[name]
        cell_rows.append(meta)

        if icell % 10 == 0 or icell == len(groups):
            print(f"[CFF separation] {icell:3d}/{len(groups)} cells")

    return pd.DataFrame(point_rows), pd.DataFrame(cell_rows)


def layout_for(cell_meta, mode):
    cffs = FIT_MODES[mode]
    pars = []
    for b in cell_meta["bin"].astype(int):
        for name in cffs:
            pars.append((b, name))
    pars += [(-1, "xs_scale_beta"), (-1, "bsa_scale_beta")]
    return pars


def covariance(df, deriv, cell_meta, mode):
    use = df.merge(deriv, on=["point_id", "bin", "phi_deg"], how="inner")
    layout = layout_for(cell_meta, mode)
    idx = {k: i for i, k in enumerate(layout)}
    npar = len(layout)
    cffs = FIT_MODES[mode]

    rows, sigmas = [], []

    for r in use.itertuples(index=False):
        v = np.zeros(npar)
        b = int(r.bin)
        for name in cffs:
            v[idx[(b, name)]] = getattr(r, f"d_xs_d_{name}")
        v[idx[(-1, "xs_scale_beta")]] = r.xs_scale_frac * r.xs_km15
        s = math.hypot(r.xs_stat_pseudo_abs, r.xs_ptp_sys_pseudo_abs)
        if np.isfinite(s) and s > 0:
            rows.append(v)
            sigmas.append(s)

    for r in use.itertuples(index=False):
        v = np.zeros(npar)
        b = int(r.bin)
        for name in cffs:
            v[idx[(b, name)]] = getattr(r, f"d_bsa_d_{name}")
        v[idx[(-1, "bsa_scale_beta")]] = r.bsa_scale_frac * r.bsa_km15
        s = math.hypot(r.bsa_stat_pseudo_abs, r.bsa_ptp_sys_pseudo_abs)
        if np.isfinite(s) and s > 0:
            rows.append(v)
            sigmas.append(s)

    A = np.asarray(rows, dtype=float)
    s = np.asarray(sigmas, dtype=float)
    B = A / s[:, None]

    # Unit Gaussian priors only on the two experimental scale nuisances.
    priors = []
    for key in [(-1, "xs_scale_beta"), (-1, "bsa_scale_beta")]:
        p = np.zeros(npar)
        p[idx[key]] = 1.0
        priors.append(p)
    B = np.vstack([B] + priors)

    # Explicit SVD so we never silently interpret an unconstrained null
    # direction as a zero uncertainty.
    U, sv, Vt = np.linalg.svd(B, full_matrices=False)
    tol = max(sv[0], 1.0) * SVD_RTOL
    keep = sv > tol
    rank = int(np.sum(keep))
    rank_deficit = npar - rank

    inv_s2 = np.zeros_like(sv)
    inv_s2[keep] = 1.0 / (sv[keep] ** 2)
    cov = (Vt.T * inv_s2) @ Vt

    # Fraction of each parameter basis vector lying in unresolved SVD modes.
    # This is zero for a fully identified fit.
    unresolved_fraction = np.zeros(npar)
    if rank_deficit > 0:
        null = Vt[~keep, :]
        unresolved_fraction = np.sum(null * null, axis=0)

    cond = float(sv[keep][0] / sv[keep][-1]) if rank > 1 else np.inf

    diag = {
        "fit_mode": mode,
        "n_measurement_rows": len(A),
        "n_parameters": npar,
        "rank": rank,
        "rank_deficit": rank_deficit,
        "condition_number_effective": cond,
        "smallest_kept_singular_value": float(sv[keep][-1]) if rank else np.nan,
        "largest_singular_value": float(sv[0]) if len(sv) else np.nan,
    }
    return cov, layout, unresolved_fraction, diag


def summarize(cell_meta, cov, layout, unresolved, mode, factor, relative_floor):
    idx = {k: i for i, k in enumerate(layout)}
    cffs = FIT_MODES[mode]
    rows = []

    for c in cell_meta.itertuples(index=False):
        b = int(c.bin)
        row = {
            "fit_mode": mode,
            "luminosity_factor": factor,
            "bin": b,
            "xB": c.xB,
            "xi": c.xi,
            "Q2": c.Q2,
            "t_abs": c.t_abs,
            "t_over_Q2": c.t_over_Q2,
            "n_phi": c.n_phi,
        }

        local_indices = [idx[(b, name)] for name in cffs]

        for name in ALL_CFFS:
            truth = getattr(c, f"{name}_KM15")
            row[f"{name}_KM15"] = truth
            row[f"sigma_{name}"] = np.nan
            row[f"relative_sigma_{name}"] = np.nan
            row[f"unresolved_fraction_{name}"] = np.nan

        for name in cffs:
            j = idx[(b, name)]
            sig = math.sqrt(max(cov[j, j], 0.0))
            truth = getattr(c, f"{name}_KM15")
            row[f"sigma_{name}"] = sig
            row[f"unresolved_fraction_{name}"] = unresolved[j]
            if abs(truth) >= relative_floor and unresolved[j] < 1e-8:
                row[f"relative_sigma_{name}"] = sig / abs(truth)

        # Pairwise local correlations among all floated CFFs.
        for i, n1 in enumerate(cffs):
            j1 = idx[(b, n1)]
            for n2 in cffs[i + 1 :]:
                j2 = idx[(b, n2)]
                den = math.sqrt(max(cov[j1, j1] * cov[j2, j2], 0.0))
                row[f"corr_{n1}_{n2}"] = cov[j1, j2] / den if den > 0 else np.nan

        # Maximum absolute correlation of each CFF with another floated CFF.
        for name in cffs:
            vals = []
            j1 = idx[(b, name)]
            for other in cffs:
                if other == name:
                    continue
                j2 = idx[(b, other)]
                den = math.sqrt(max(cov[j1, j1] * cov[j2, j2], 0.0))
                if den > 0:
                    vals.append(abs(cov[j1, j2] / den))
            row[f"max_abs_corr_{name}"] = max(vals) if vals else np.nan

        rows.append(row)

    return pd.DataFrame(rows)


def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def finite_rel(d, col):
    return d[np.isfinite(d[col])].copy()


def make_plots(results, cell_meta, deriv, figures):
    # 1: Im Htilde relative precision vs xB, H+Htilde fit.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for factor in LUMI_FACTORS:
        d = results[
            (results.fit_mode == "H_Ht") &
            (results.luminosity_factor == factor)
        ]
        d = finite_rel(d, "relative_sigma_ImHt")
        ax.scatter(d.xB, 100 * d.relative_sigma_ImHt, s=30, alpha=.70, label=f"{factor}x")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$ (%)")
    ax.set_title(r"Local $H+\widetilde H$ fit: luminosity dependence")
    ax.grid(alpha=.2); ax.legend(title="Luminosity")
    savefig(fig, figures / "01_ImHt_relative_uncertainty_vs_xB.png")

    # 2: same versus xi.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for factor in LUMI_FACTORS:
        d = results[
            (results.fit_mode == "H_Ht") &
            (results.luminosity_factor == factor)
        ]
        d = finite_rel(d, "relative_sigma_ImHt")
        ax.scatter(d.xi, 100 * d.relative_sigma_ImHt, s=30, alpha=.70, label=f"{factor}x")
    ax.set_xlabel(r"$\xi \simeq x_B/(2-x_B)$")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$ (%)")
    ax.set_title(r"Does $\widetilde H$ sensitivity emerge at high $\xi$?")
    ax.grid(alpha=.2); ax.legend(title="Luminosity")
    savefig(fig, figures / "02_ImHt_relative_uncertainty_vs_xi.png")

    # 3: freeing E penalty at 10x.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for mode, label in [
        ("H_Ht", r"$H+\widetilde H$"),
        ("H_Ht_E", r"$H+\widetilde H+E$"),
    ]:
        d = results[(results.fit_mode == mode) & (results.luminosity_factor == 10)]
        d = finite_rel(d, "relative_sigma_ImHt")
        ax.scatter(d.xB, 100 * d.relative_sigma_ImHt, s=32, alpha=.72, label=label)
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$ (%)")
    ax.set_title(r"CFF-separation penalty when $\mathcal{E}$ is also free at 10x")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures / "03_ImHt_HHt_vs_HHtE_10x.png")

    # 4: absolute uncertainty, robust through zero crossings.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for factor in LUMI_FACTORS:
        d = results[
            (results.fit_mode == "H_Ht_E") &
            (results.luminosity_factor == factor)
        ]
        ax.scatter(d.xB, d.sigma_ImHt, s=30, alpha=.70, label=f"{factor}x")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Absolute uncertainty on $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.set_title(r"Full local $H+\widetilde H+E$ fit")
    ax.grid(alpha=.2); ax.legend(title="Luminosity")
    savefig(fig, figures / "04_ImHt_absolute_uncertainty_full_fit.png")

    # 5: 10x kinematic map.  Marker size inversely tracks absolute uncertainty.
    d = results[(results.fit_mode == "H_Ht_E") & (results.luminosity_factor == 10)].copy()
    finite = np.isfinite(d.sigma_ImHt) & (d.sigma_ImHt > 0)
    d = d[finite]
    inv = 1.0 / d.sigma_ImHt
    sizes = 25 + 175 * (inv - inv.min()) / max(inv.max() - inv.min(), 1e-30)
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    sc = ax.scatter(d.xB, d.t_abs, c=d.sigma_ImHt, s=sizes, alpha=.78)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$\sigma[\mathrm{Im}\,\widetilde{\mathcal{H}}]$")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"$|t|$ (GeV$^2$)")
    ax.set_title(r"Where does XS+BSA constrain $\widetilde H$?  Full fit, 10x")
    ax.grid(alpha=.2)
    savefig(fig, figures / "05_ImHt_sensitivity_map_10x.png")

    # 6: maximum local CFF correlation for ImHt.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for mode, label in [
        ("H_Ht", r"$H+\widetilde H$"),
        ("H_Ht_E", r"$H+\widetilde H+E$"),
    ]:
        d = results[(results.fit_mode == mode) & (results.luminosity_factor == 10)]
        col = "max_abs_corr_ImHt"
        ax.scatter(d.xB, d[col], s=30, alpha=.72, label=label)
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Maximum $|\rho|$ involving $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.set_ylim(0, 1.03)
    ax.set_title(r"Is $\widetilde H$ precision limited by CFF separation?")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures / "06_ImHt_max_correlation_10x.png")

    # 7: observable-level response to a +10% ImHt change.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    ax.scatter(
        cell_meta.xB,
        100 * cell_meta.median_abs_xs_frac_response_to_10pct_ImHt,
        s=34, alpha=.72, label="XS: median fractional response"
    )
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Median $|\Delta\sigma/\sigma|$ for a 10% change in $\mathrm{Im}\,\widetilde{\mathcal{H}}$ (%)")
    ax.set_title(r"Direct observable sensitivity to $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.grid(alpha=.2)
    savefig(fig, figures / "07_XS_response_to_10pct_ImHt.png")

    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    ax.scatter(
        cell_meta.xB,
        cell_meta.median_abs_bsa_response_to_10pct_ImHt,
        s=34, alpha=.72
    )
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Median $|\Delta A_{LU}|$ for a 10% change in $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.set_title(r"Direct BSA sensitivity to $\mathrm{Im}\,\widetilde{\mathcal{H}}$")
    ax.grid(alpha=.2)
    savefig(fig, figures / "08_BSA_response_to_10pct_ImHt.png")

    # 9: H degradation as progressively more CFFs are released.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    for mode, label in [
        ("H", r"$H$ only"),
        ("H_Ht", r"$H+\widetilde H$"),
        ("H_Ht_E", r"$H+\widetilde H+E$"),
    ]:
        d = results[(results.fit_mode == mode) & (results.luminosity_factor == 10)]
        d = finite_rel(d, "relative_sigma_ImH")
        ax.scatter(d.xB, 100 * d.relative_sigma_ImH, s=30, alpha=.70, label=label)
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Im}\,\mathcal{H}$ (%)")
    ax.set_title(r"Cost of progressively relaxing the $H$-dominance assumption at 10x")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures / "09_ImH_degradation_when_subleading_CFFs_float.png")


def print_summary(results, diagnostics):
    print("\n" + "=" * 108)
    print("LOCAL H / Htilde / E CFF-SEPARATION PROJECTION")
    print("=" * 108)

    rows = []
    for mode in FIT_MODES:
        for factor in LUMI_FACTORS:
            d = results[
                (results.fit_mode == mode) &
                (results.luminosity_factor == factor)
            ]
            diag = diagnostics[
                (diagnostics.fit_mode == mode) &
                (diagnostics.luminosity_factor == factor)
            ].iloc[0]
            rec = {
                "mode": mode,
                "L": f"{factor}x",
                "rank_def": int(diag.rank_deficit),
                "cond": diag.condition_number_effective,
                "med_rel_ImH_%": 100 * np.nanmedian(d.relative_sigma_ImH),
                "med_rel_ImHt_%": np.nan,
                "med_abs_ImHt": np.nan,
                "med_maxcorr_ImHt": np.nan,
            }
            if "ImHt" in FIT_MODES[mode]:
                rec["med_rel_ImHt_%"] = 100 * np.nanmedian(d.relative_sigma_ImHt)
                rec["med_abs_ImHt"] = np.nanmedian(d.sigma_ImHt)
                rec["med_maxcorr_ImHt"] = np.nanmedian(d.max_abs_corr_ImHt)
            rows.append(rec)

    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # Explicitly quantify the penalty from freeing E.
    print("\n10x Htilde separation penalty from freeing E:")
    a = results[(results.fit_mode == "H_Ht") & (results.luminosity_factor == 10)]
    b = results[(results.fit_mode == "H_Ht_E") & (results.luminosity_factor == 10)]
    m = a[["bin", "sigma_ImHt"]].merge(
        b[["bin", "sigma_ImHt"]], on="bin", suffixes=("_HHt", "_HHtE")
    )
    good = np.isfinite(m.sigma_ImHt_HHt) & np.isfinite(m.sigma_ImHt_HHtE) & (m.sigma_ImHt_HHt > 0)
    ratio = m.loc[good, "sigma_ImHt_HHtE"] / m.loc[good, "sigma_ImHt_HHt"]
    if len(ratio):
        print(
            f"  median sigma(ImHt) ratio [H+Ht+E]/[H+Ht] = {np.median(ratio):.3g}"
        )
        print(
            f"  central 68% range of ratio = "
            f"{np.percentile(ratio,16):.3g} -- {np.percentile(ratio,84):.3g}"
        )


def main(argv: Optional[List[str]] = None) -> int:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input-dir",
        default=str(here / "output" / "stage2" / "tables"),
        help="Directory containing joint_fit_input_km15_{1,2,5,10}x.csv",
    )
    p.add_argument(
        "--outdir",
        default=str(here / "output" / "stage4_cff_separation"),
    )
    p.add_argument("--min-nphi", type=int, default=DEFAULT_MIN_NPHI)
    p.add_argument("--relative-step", type=float, default=DEFAULT_REL_STEP)
    p.add_argument("--absolute-step", type=float, default=DEFAULT_ABS_STEP)
    p.add_argument(
        "--relative-floor",
        type=float,
        default=DEFAULT_RELATIVE_FLOOR,
        help="Do not quote relative CFF error when |CFF_truth| is below this.",
    )
    args = p.parse_args(argv)

    indir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    tables = outdir / "tables"
    figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    inputs = load_inputs(indir)
    one = inputs[1]

    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError("Could not import Gepard/KM15.") from exc
    th = th_KM15

    # Preflight verifies the exact Gepard interface and demonstrates observable
    # response before the expensive derivative pass.
    row = one.iloc[len(one) // 2]
    xs0 = make_point(g, row, 0)
    bp = make_point(g, row, +1)
    bm = make_point(g, row, -1)

    print("[preflight] representative KM15 point")
    print(
        f"[preflight] xB={row.xB:.4g}, Q2={row.Q2:.4g}, "
        f"|t|={row.t_abs:.4g}, phi={row.phi_deg:.4g}"
    )
    print(
        f"[preflight] XS={pred(th,xs0):.8g}, "
        f"BSA={bsa_from_points(th,bp,bm):+.8g}"
    )
    for name in ALL_CFFS:
        val = cff_value(th, name, xs0)
        step = max(args.absolute_step, args.relative_step * max(abs(val), 1.0))
        base_xs = pred(th, xs0)
        base_a = bsa_from_points(th, bp, bm)
        with shifted_cff(th, name, step):
            xs1 = pred(th, xs0)
            a1 = bsa_from_points(th, bp, bm)
        print(
            f"[preflight] {name:4s}: value={val:+.6g}, step={step:.4g}, "
            f"XS response={100*(xs1/base_xs-1):+.4g}%, "
            f"BSA change={a1-base_a:+.4g}"
        )

    deriv, cell_meta = build_derivatives(
        one,
        th,
        g,
        args.relative_step,
        args.absolute_step,
        args.min_nphi,
    )
    deriv.to_csv(tables / "point_xs_bsa_allcff_derivatives.csv", index=False)
    cell_meta.to_csv(tables / "cell_km15_allcffs.csv", index=False)

    result_frames = []
    diag_rows = []

    for mode in FIT_MODES:
        for factor in LUMI_FACTORS:
            print(f"[fit] mode={mode:7s} luminosity={factor}x")
            cov, layout, unresolved, diag = covariance(
                inputs[factor], deriv, cell_meta, mode
            )
            diag["luminosity_factor"] = factor
            diag_rows.append(diag)
            result_frames.append(
                summarize(
                    cell_meta,
                    cov,
                    layout,
                    unresolved,
                    mode,
                    factor,
                    args.relative_floor,
                )
            )

    results = pd.concat(result_frames, ignore_index=True)
    diagnostics = pd.DataFrame(diag_rows)

    results.to_csv(tables / "cff_separation_cell_uncertainties.csv", index=False)
    diagnostics.to_csv(tables / "cff_separation_fit_diagnostics.csv", index=False)

    # Save a compact direct-sensitivity table for easy inspection.
    sens_cols = [
        "bin", "xB", "xi", "Q2", "t_abs", "t_over_Q2", "n_phi",
        "ReH_KM15", "ImH_KM15", "ReHt_KM15", "ImHt_KM15",
        "ReE_KM15", "ImE_KM15",
        "median_abs_xs_frac_response_to_10pct_ImHt",
        "median_abs_bsa_response_to_10pct_ImHt",
    ]
    cell_meta[sens_cols].to_csv(
        tables / "ImHt_direct_observable_sensitivity.csv", index=False
    )

    make_plots(results, cell_meta, deriv, figures)
    print_summary(results, diagnostics)

    print("\n[output]", outdir)
    print("[interpretation]")
    print("  H       : optimistic H-dominance reference.")
    print("  H_Ht    : asks whether XS+BSA can separate Htilde once H is free.")
    print("  H_Ht_E  : asks how much that conclusion survives when E is also free.")
    print("  Absolute errors remain meaningful near CFF zero crossings.")
    print("  Relative errors are suppressed below the configured truth floor.")
    print("  Any nonzero rank deficit is a warning that the corresponding fit")
    print("  contains an exactly unresolved CFF combination and must not be")
    print("  interpreted from pseudoinverse errors alone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
