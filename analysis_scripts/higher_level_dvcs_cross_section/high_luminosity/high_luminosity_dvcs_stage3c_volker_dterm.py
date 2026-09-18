#!/usr/bin/env python3
"""
Stage 3C: Volker/Burkert-style generalized-multipole D-term projection.

This stage keeps the KM15 ImH ansatz and fixed-t dispersion integral, but
replaces KM15's fixed dipole subtraction function by the functional family
used in the Burkert-Elouadrhiri-Girod Nature analysis and the later Volker
presentation:

    d1^Q(t) = d1^Q(0) * (1 - t/M^2)^(-alpha)

with d1^Q(0), M^2, and alpha all allowed to vary.

The default pseudo-truth uses the values displayed in the Volker presentation:
    d1^Q(0) = -2.04
    M^2     =  1.02 GeV^2
    alpha   =  2.76

The 2018 Nature paper used the same multipole functional family; its text
describes alpha=3 as the chosen asymptotic form.  Here alpha is deliberately
floated so that the high-luminosity projection tests how well the new data can
determine curvature rather than imposing the dipole/tripole shape.

Primary story:
  1. statistics-only: what high luminosity can buy at high |t|;
  2. baseline systematics: how much of that statistical potential survives;
  3. PTP/2: a realistic analysis-improvement target.

The existing CLAS12 |t| reach is not treated as a future-upgrade gain.

The actual observable derivatives are still calculated through Gepard, so
changes to the ImH parameters propagate into ReH through the fixed-t
dispersion relation.  Only the subtraction function is generalized.

Outputs are written to the stable directory output/stage3c/.
"""

from __future__ import annotations

import argparse
import copy
import math
import os
import time
import types
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import gepard as g
from gepard.fits import th_KM15


LUMI_FACTORS = (1, 2, 5, 10)
SCENARIOS = ("baseline", "ptp_half", "statistics_only")

# Historical CLAS6 input range used in the 2018 D-term analysis chain.
# The public dterm18 CLAS files used by Kumericki contain -t values from
# 0.11 through 0.45 GeV^2.  We show this only as historical kinematic
# context; it is NOT digitized C_H(t) data and is not used in the fit.
CLAS6_DTERM_TMIN = 0.11
CLAS6_DTERM_TMAX = 0.45

# Present pass-2 theory-controlled sample (|t|/Q^2 < 0.2).  The exact
# endpoints are recomputed from the Stage-2 input at runtime; these values are
# only fallbacks for labels before input is loaded.
CLAS12_CONTROLLED_TMIN_FALLBACK = 0.129
CLAS12_CONTROLLED_TMAX_FALLBACK = 0.912

# Parameters entering DispersionFixedPoleCFF.ImH in Gepard.
IMH_PARAMETER_CANDIDATES = (
    "Nsea", "alS", "alpS", "mS2", "rS", "bS",
    "Nv", "alv", "alpv", "mv2", "rv", "bv",
)
TARGET_PARAMETERS = ("C", "dM2", "dalpha")

# Volker-presentation central values for the generalized multipole.
VOLKER_D1_0 = -2.04
VOLKER_M2 = 1.02
VOLKER_ALPHA = 2.76
VOLKER_C = -VOLKER_D1_0 / 0.9

def _volker_subtraction(self, pt):
    """Positive Gepard subtraction S(t)=-C_H(t) for the generalized multipole."""
    p = self.parameters
    return p["C"] / (1.0 - pt.t/p["dM2"])**p["dalpha"]


def configure_volker_subtraction(model, dterm_parameters=None):
    """Install the generalized multipole on an independent Gepard model."""
    vals = {
        "C": VOLKER_C,
        "dM2": VOLKER_M2,
        "dalpha": VOLKER_ALPHA,
    }
    if dterm_parameters:
        vals.update(dterm_parameters)
    model.parameters.update(vals)
    # Keep limits local to this projection; they are used only to choose safe
    # symmetric finite steps.
    if not hasattr(model, "parameters_limits"):
        model.parameters_limits = {}
    model.parameters_limits.update({
        "C": (0.2, 8.0),
        "dM2": (0.16, 4.0),
        "dalpha": (1.0, 6.0),
    })
    model.subtraction = types.MethodType(_volker_subtraction, model)
    return model



def make_point(xB: float, Q2: float, t_abs: float, phi_deg: float,
               ebeam: float, polarization: float = 0.0):
    phi_trento = math.pi - math.radians(phi_deg)
    pt = g.DataPoint(
        xB=float(xB), t=-abs(float(t_abs)), Q2=float(Q2), phi=phi_trento,
        observable="XS", frame="trento", process="ep2epgamma",
        exptype="fixed target", in1energy=float(ebeam), in1charge=-1,
        in1polarization=float(polarization), in2particle="p",
    )
    pt.prepare()
    return pt


def predict_xs_bsa(model, row) -> Tuple[float, float]:
    ebeam = float(getattr(row, "ebeam", 10.6041))
    p0 = make_point(row.xB, row.Q2, row.t_abs, row.phi_deg, ebeam, 0.0)
    xs = float(model.predict(p0))

    pp = make_point(row.xB, row.Q2, row.t_abs, row.phi_deg, ebeam, +1.0)
    pm = make_point(row.xB, row.Q2, row.t_abs, row.phi_deg, ebeam, -1.0)
    sp = float(model.predict(pp))
    sm = float(model.predict(pm))
    bsa = (sp - sm) / (sp + sm)
    return xs, bsa


def choose_step(name: str, value: float, limits: Dict) -> float:
    """Finite step: normally 1%, with safeguards near zero and parameter limits."""
    scale = max(abs(float(value)), 0.1)
    step = 0.01 * scale

    # Shape masses benefit from a slightly larger numerical step.
    if name in ("dM2", "mC2", "mS2", "mv2"):
        step = 0.005 * max(abs(float(value)), 0.5)

    lo, hi = limits.get(name, (-np.inf, np.inf))
    if np.isfinite(hi):
        step = min(step, 0.2 * max(hi - value, 1e-8))
    if np.isfinite(lo):
        step = min(step, 0.2 * max(value - lo, 1e-8))
    return max(step, 1e-6)



def _rows_to_records(base: pd.DataFrame):
    """Small, picklable representation of the kinematics sent to workers."""
    cols = ["point_id", "bin", "xB", "Q2", "t_abs", "phi_deg"]
    if "ebeam" in base.columns:
        cols.append("ebeam")
    return base[cols].to_dict("records")


class _Row:
    """Attribute-style view used by predict_xs_bsa without pandas overhead."""
    __slots__ = ("point_id", "bin", "xB", "Q2", "t_abs", "phi_deg", "ebeam")
    def __init__(self, r):
        self.point_id = r["point_id"]
        self.bin = r["bin"]
        self.xB = r["xB"]
        self.Q2 = r["Q2"]
        self.t_abs = r["t_abs"]
        self.phi_deg = r["phi_deg"]
        self.ebeam = r.get("ebeam", 10.6041)


def _predict_records(model, records, progress_label=None, progress_every=100):
    xs, bsa = [], []
    n = len(records)
    t0 = time.time()
    for i, rec in enumerate(records, 1):
        x, a = predict_xs_bsa(model, _Row(rec))
        xs.append(x)
        bsa.append(a)
        if progress_label and (i % progress_every == 0 or i == n):
            dt = time.time() - t0
            rate = i / dt if dt > 0 else 0.0
            eta = (n-i)/rate if rate > 0 else float("nan")
            print(f"[{progress_label}] {i}/{n}  elapsed={dt:.1f}s  ETA={eta:.1f}s", flush=True)
    return np.asarray(xs), np.asarray(bsa)


def _parameter_worker(payload):
    """
    One process owns one independent KM15 copy and computes all finite-step
    information for one parameter.  Parallelizing by parameter is deliberate:
    Gepard model mutation is never shared between processes, and each worker
    gets a large enough block of work that process overhead is negligible.
    """
    par, value, h0, records, do_stability = payload
    model = configure_volker_subtraction(copy.deepcopy(th_KM15))

    results = {"parameter": par, "value": value, "step": h0}

    factors = (1.0, 0.5, 2.0) if do_stability else (1.0,)
    for factor in factors:
        h = h0 * factor

        model.parameters[par] = value + h
        px, pa = _predict_records(model, records)

        model.parameters[par] = value - h
        mx, ma = _predict_records(model, records)

        model.parameters[par] = value
        results[f"dx_{factor:g}"] = (px - mx) / (2*h)
        results[f"da_{factor:g}"] = (pa - ma) / (2*h)

    return results


def _central_worker(records):
    model = configure_volker_subtraction(copy.deepcopy(th_KM15))
    return _predict_records(model, records, progress_label="central", progress_every=100)


def finite_derivatives(
    base: pd.DataFrame,
    model,
    parameters: List[str],
    cache_csv: Path,
    workers: int = 1,
    compute_stability: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute central predictions and finite derivatives.

    Expensive parameter variations are parallelized across independent
    processes.  Each process owns its own KM15 object; no mutable Gepard state
    is shared.  The same worker also computes h/2 and 2h derivatives, avoiding
    a second complete validation pass.
    """
    steps_path = cache_csv.with_name("finite_difference_steps.csv")
    stability_path = cache_csv.with_name("finite_difference_stability.csv")

    if cache_csv.exists() and steps_path.exists():
        d = pd.read_csv(cache_csv)
        needed = {"point_id", "xs0", "bsa0"}
        for p in parameters:
            needed |= {f"d_xs_d_{p}", f"d_bsa_d_{p}"}
        if needed.issubset(d.columns) and len(d) == len(base):
            steps = pd.read_csv(steps_path)
            if stability_path.exists():
                stability = pd.read_csv(stability_path)
                print(f"[cache] using derivatives + stability from {cache_csv.parent}")
                return d, steps, stability
            if not compute_stability:
                return d, steps, pd.DataFrame()

    central = copy.deepcopy(model.parameters)
    limits = getattr(model, "parameters_limits", {})
    records = _rows_to_records(base)

    print(f"[performance] workers={workers}; parallel unit=KM15 parameter")
    print(f"[derivatives] central predictions for {len(base)} points")
    xs0, a0 = _central_worker(records)

    out = base[["point_id", "bin", "xB", "Q2", "t_abs", "phi_deg"]].copy()
    out["xs0"] = xs0
    out["bsa0"] = a0

    payloads = []
    step_info = {}
    for par in parameters:
        val = float(central[par])
        h = choose_step(par, val, limits)
        step_info[par] = (val, h)
        payloads.append((par, val, h, records, compute_stability))

    results = {}
    if workers <= 1:
        for i, payload in enumerate(payloads, 1):
            par = payload[0]
            print(f"[derivatives] {i}/{len(payloads)} {par}")
            results[par] = _parameter_worker(payload)
    else:
        # Linux ifarm: fork is both fast and compatible with the existing
        # scientific Python/Gepard environment.  Explicit context also avoids
        # Python-version-dependent defaults.
        ctx = mp.get_context("fork")
        with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
            futures = {ex.submit(_parameter_worker, x): x[0] for x in payloads}
            done = 0
            for fut in as_completed(futures):
                par = futures[fut]
                results[par] = fut.result()
                done += 1
                print(f"[derivatives] completed {done}/{len(payloads)} {par}")

    steps_rows = []
    stability_rows = []

    for par in parameters:
        r = results[par]
        val, h = step_info[par]
        dx = r["dx_1"]
        da = r["da_1"]
        out[f"d_xs_d_{par}"] = dx
        out[f"d_bsa_d_{par}"] = da

        rel_x = np.nanmedian(np.abs(dx*h) / np.maximum(np.abs(xs0), 1e-30))
        abs_a = np.nanmedian(np.abs(da*h))
        steps_rows.append({
            "parameter": par, "central_value": val, "step": h,
            "median_fractional_XS_half_response": rel_x,
            "median_absolute_BSA_half_response": abs_a,
        })

        if compute_stability:
            ref = np.concatenate([dx, da])
            ref_norm = float(np.linalg.norm(ref))
            half = np.concatenate([r["dx_0.5"], r["da_0.5"]])
            double = np.concatenate([r["dx_2"], r["da_2"]])
            dh = float(np.linalg.norm(half-ref) / max(ref_norm, 1e-300))
            dd = float(np.linalg.norm(double-ref) / max(ref_norm, 1e-300))
            stability_rows.append({
                "parameter": par,
                "central_step": h,
                "half_step_relative_vector_change": dh,
                "double_step_relative_vector_change": dd,
                "max_relative_vector_change": max(dh, dd),
            })

    steps = pd.DataFrame(steps_rows)
    stability = pd.DataFrame(stability_rows)
    out.to_csv(cache_csv, index=False)
    steps.to_csv(steps_path, index=False)
    if compute_stability:
        stability.to_csv(stability_path, index=False)
    return out, steps, stability


def classify_active_parameters(
    deriv: pd.DataFrame,
    steps: pd.DataFrame,
    requested: List[str],
    mandatory=("C", "dM2", "dalpha"),
    response_threshold: float = 1e-8,
    boundary_step_fraction: float = 1e-4,
) -> Tuple[List[str], pd.DataFrame]:
    """
    Keep only genuine local directions of the installed KM15 solution.

    A candidate ImH parameter is excluded when:
      1) its finite change produces essentially zero XS and BSA response, or
      2) the available symmetric step has collapsed against a parameter
         boundary (the bv situation in KM15).

    C, dM2, and dalpha are mandatory because they define the subtraction function.
    """
    step_map = steps.set_index("parameter").to_dict("index")
    rows = []
    active = []

    for p in requested:
        info = step_map.get(p, {})
        value = float(info.get("central_value", np.nan))
        step = float(info.get("step", np.nan))
        scale = max(abs(value), 0.1) if np.isfinite(value) else 1.0
        step_fraction = step / scale if np.isfinite(step) else np.nan

        dx = np.asarray(deriv[f"d_xs_d_{p}"], dtype=float)
        da = np.asarray(deriv[f"d_bsa_d_{p}"], dtype=float)
        xs0 = np.maximum(np.abs(np.asarray(deriv["xs0"], dtype=float)), 1e-30)

        # Translate derivative back to the actual response produced by the
        # finite step, and use the 90th percentile so parameters that matter
        # only in part of phase space are not discarded by a zero median.
        rel_xs = np.abs(dx * step) / xs0
        abs_bsa = np.abs(da * step)
        p90_xs = float(np.nanpercentile(rel_xs, 90))
        p90_bsa = float(np.nanpercentile(abs_bsa, 90))
        response = max(p90_xs, p90_bsa)

        boundary_limited = (
            p not in mandatory
            and np.isfinite(step_fraction)
            and step_fraction < boundary_step_fraction
        )
        zero_response = p not in mandatory and response < response_threshold
        keep = (p in mandatory) or (not boundary_limited and not zero_response)

        if keep:
            active.append(p)

        if p in mandatory:
            reason = "mandatory subtraction parameter"
        elif boundary_limited:
            reason = "excluded: symmetric step collapsed at/near parameter boundary"
        elif zero_response:
            reason = "excluded: no observable response around installed KM15 point"
        else:
            reason = "active"

        rows.append({
            "parameter": p,
            "central_value": value,
            "step": step,
            "step_fraction_of_parameter_scale": step_fraction,
            "p90_fractional_XS_step_response": p90_xs,
            "p90_absolute_BSA_step_response": p90_bsa,
            "activity_response": response,
            "active": keep,
            "reason": reason,
        })

    return active, pd.DataFrame(rows)




def remove_exact_response_degeneracies(
    deriv: pd.DataFrame,
    parameters: List[str],
    mandatory=("C", "dM2", "dalpha"),
    collinearity_tol: float = 1e-10,
) -> Tuple[List[str], pd.DataFrame]:
    """
    Remove parameter directions whose observable derivative is exactly
    collinear with one already retained.

    This is not an arbitrary numerical regularization.  For the Gepard
    DispersionFixedPoleCFF ansatz, for example, the valence normalization
    enters ImH as Nv*rv, so Nv and rv cannot both be independently determined
    by these observables.  We retain one representative normalization
    parameter and document the removed direction.
    """
    kept = []
    rows = []

    # Prefer mandatory subtraction parameters, then preserve original order.
    ordered = list(mandatory) + [p for p in parameters if p not in mandatory]

    for p in ordered:
        v = np.concatenate([
            np.asarray(deriv[f"d_xs_d_{p}"], dtype=float),
            np.asarray(deriv[f"d_bsa_d_{p}"], dtype=float),
        ])
        nv = np.linalg.norm(v)
        redundant_with = None
        cosine = np.nan

        for q in kept:
            w = np.concatenate([
                np.asarray(deriv[f"d_xs_d_{q}"], dtype=float),
                np.asarray(deriv[f"d_bsa_d_{q}"], dtype=float),
            ])
            nw = np.linalg.norm(w)
            if nv == 0 or nw == 0:
                continue
            c = float(np.dot(v, w)/(nv*nw))
            # Exact same or opposite response direction.
            if abs(abs(c)-1.0) < collinearity_tol:
                residual = np.linalg.norm(v - (np.dot(v,w)/np.dot(w,w))*w) / nv
                if residual < collinearity_tol:
                    redundant_with = q
                    cosine = c
                    break

        if redundant_with is None:
            kept.append(p)
            rows.append({
                "parameter": p, "kept": True,
                "redundant_with": "", "cosine": np.nan,
                "reason": "independent observable-response direction",
            })
        else:
            rows.append({
                "parameter": p, "kept": False,
                "redundant_with": redundant_with, "cosine": cosine,
                "reason": f"exactly redundant observable response with {redundant_with}",
            })

    return kept, pd.DataFrame(rows)


def scenario_input(d: pd.DataFrame, scenario: str) -> pd.DataFrame:
    q = d.copy()
    if scenario == "baseline":
        return q
    if scenario == "ptp_half":
        q["xs_ptp_sys_pseudo_abs"] *= 0.5
        q["bsa_ptp_sys_pseudo_abs"] *= 0.5
        return q
    if scenario == "statistics_only":
        q["xs_ptp_sys_pseudo_abs"] = 0.0
        q["bsa_ptp_sys_pseudo_abs"] = 0.0
        q["xs_scale_frac"] = 0.0
        q["bsa_scale_frac"] = 0.0
        return q
    raise ValueError(scenario)


def covariance(data: pd.DataFrame, deriv: pd.DataFrame,
               parameters: List[str], central: Dict[str, float]):
    use = data.merge(
        deriv, on=["point_id", "bin", "xB", "Q2", "t_abs", "phi_deg"],
        how="inner", validate="one_to_one"
    )

    layout = list(parameters) + ["beta_xs_norm", "beta_bsa_pol"]
    idx = {p: i for i, p in enumerate(layout)}
    rows, sigmas = [], []

    for r in use.itertuples(index=False):
        # XS
        v = np.zeros(len(layout))
        for p in parameters:
            v[idx[p]] = float(getattr(r, f"d_xs_d_{p}"))
        v[idx["beta_xs_norm"]] = float(r.xs_scale_frac) * float(r.xs0)
        sx = math.hypot(float(r.xs_stat_pseudo_abs),
                        float(r.xs_ptp_sys_pseudo_abs))
        if sx > 0 and np.isfinite(sx):
            rows.append(v); sigmas.append(sx)

        # BSA
        v = np.zeros(len(layout))
        for p in parameters:
            v[idx[p]] = float(getattr(r, f"d_bsa_d_{p}"))
        v[idx["beta_bsa_pol"]] = float(r.bsa_scale_frac) * float(r.bsa0)
        sa = math.hypot(float(r.bsa_stat_pseudo_abs),
                        float(r.bsa_ptp_sys_pseudo_abs))
        if sa > 0 and np.isfinite(sa):
            rows.append(v); sigmas.append(sa)

    A = np.asarray(rows)
    sig = np.asarray(sigmas)
    B = A / sig[:, None]

    # Unit priors only on experimental scale nuisances.
    for nuisance in ("beta_xs_norm", "beta_bsa_pol"):
        q = np.zeros(len(layout))
        q[idx[nuisance]] = 1.0
        B = np.vstack([B, q])

    info = B.T @ B
    cov = np.linalg.pinv(info, rcond=1e-12)
    rank = int(np.linalg.matrix_rank(B))
    sv = np.linalg.svd(B, compute_uv=False)
    good = sv[sv > max(sv[0], 1.0)*1e-12]
    cond = float(good[0]/good[-1]) if len(good) > 1 else np.inf

    diag = {
        "n_measurement_rows": len(A),
        "n_parameters": len(layout),
        "rank": rank,
        "rank_deficit": len(layout)-rank,
        "condition_number_effective": cond,
        "xs_norm_beta_sigma": math.sqrt(max(cov[idx["beta_xs_norm"], idx["beta_xs_norm"]], 0)),
        "bsa_pol_beta_sigma": math.sqrt(max(cov[idx["beta_bsa_pol"], idx["beta_bsa_pol"]], 0)),
    }
    return cov, layout, diag


def ch_and_gradient(t_abs: np.ndarray, C: float, dM2: float, dalpha: float):
    """
    C_H(t) = -C * (1 + |t|/M^2)^(-alpha), plus analytic derivatives.
    """
    u = 1.0 + t_abs/dM2
    ch = -C * u**(-dalpha)
    dC = -u**(-dalpha)
    dM = -C * dalpha * t_abs/(dM2**2) * u**(-dalpha-1.0)
    da = C * np.log(u) * u**(-dalpha)
    return ch, dC, dM, da


def propagate_curve(cov, layout, central, tgrid, scenario, lumi):
    idx = {p:i for i,p in enumerate(layout)}
    C = central["C"]; M2 = central["dM2"]; alpha = central["dalpha"]
    ch, dC, dM, da = ch_and_gradient(tgrid, C, M2, alpha)
    rows = []
    for t, y, gC, gM, ga in zip(tgrid, ch, dC, dM, da):
        grad = np.zeros(len(layout))
        grad[idx["C"]] = gC
        grad[idx["dM2"]] = gM
        grad[idx["dalpha"]] = ga
        var = float(grad @ cov @ grad)
        sch = math.sqrt(max(var, 0))
        rows.append({
            "scenario": scenario, "luminosity_factor": lumi,
            "t_abs": float(t), "CH": float(y), "sigma_CH": sch,
            "d1Q": 0.9*float(y), "sigma_d1Q": 0.9*sch,
        })
    return pd.DataFrame(rows)


def savefig(fig, path):
    # Leave room for any annotations intentionally drawn above the axes.
    fig.tight_layout(rect=(0, 0, 1, 0.88) if len(fig.axes) and
                     any(t.get_position()[1] > 1.0 for t in fig.axes[0].texts)
                     else None)
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def add_t_range_context(ax, clas12_tmin, clas12_tmax):
    """Small in-axis range markers; context only, not part of the upgrade gain."""
    trans = ax.get_xaxis_transform()
    def bracket(x0, x1, y, text):
        ax.annotate(
            "", xy=(x1, y), xytext=(x0, y),
            xycoords=trans, textcoords=trans,
            arrowprops=dict(arrowstyle="|-|", lw=0.9, alpha=0.65),
        )
        ax.text(
            0.5*(x0+x1), y+0.012, text, transform=trans,
            ha="center", va="bottom", fontsize=8, alpha=0.72,
        )
    bracket(CLAS6_DTERM_TMIN, CLAS6_DTERM_TMAX, 0.08, "CLAS6 D-term range")
    bracket(clas12_tmin, clas12_tmax, 0.15, "CLAS12 controlled range")


def _validated_covariance(q, deriv, parameters, central, label="fit"):
    cov, layout, diag = covariance(q, deriv, parameters, central)
    if diag["rank_deficit"] != 0:
        raise RuntimeError(
            f"Rank-deficient {label}: {diag['rank']}/{diag['n_parameters']}. "
            "Do not regularize this projection with a pseudoinverse."
        )
    return cov, layout, diag



def make_plots(curves, pars, figures, clas12_tmin, clas12_tmax):
    labels = {
        "baseline":"Current point-to-point systematics",
        "ptp_half":"Point-to-point systematics / 2",
        "statistics_only":"Statistics only",
    }

    # 01: lead with the statistical potential of the luminosity upgrade.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    add_t_range_context(ax, clas12_tmin, clas12_tmax)
    for L in LUMI_FACTORS:
        d = curves[(curves.scenario=="statistics_only") &
                   (curves.luminosity_factor==L)].sort_values("t_abs")
        ax.fill_between(d.t_abs, d.d1Q-d.sigma_d1Q, d.d1Q+d.sigma_d1Q, alpha=0.16)
        ax.plot(d.t_abs, d.d1Q, label=f"{L}x")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$d_1^Q(t)$")
    ax.set_title("Statistical potential of increased luminosity")
    ax.grid(alpha=.2); ax.legend(title="Pass-2 exposure")
    savefig(fig, figures/"01_d1_statistics_only_luminosity.png")

    # 02: direct precision on both free shape parameters.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    d = pars[pars.scenario=="statistics_only"].sort_values("luminosity_factor")
    ax.plot(d.luminosity_factor, 100*d.sigma_dM2/abs(d.dM2),
            marker="o", label=r"$M^2$")
    ax.plot(d.luminosity_factor, 100*d.sigma_dalpha/abs(d.dalpha),
            marker="o", label=r"$\alpha$")
    ax.set_xscale("log")
    ax.set_xticks(LUMI_FACTORS, [f"{L}x" for L in LUMI_FACTORS])
    ax.set_xlabel("Luminosity relative to pass-2 exposure")
    ax.set_ylabel("Relative parameter uncertainty (%)")
    ax.set_title("Statistical precision on the D-term shape")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures/"02_shape_parameters_statistics_only.png")

    # 03: then introduce the realistic systematic limitation at 10x.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    add_t_range_context(ax, clas12_tmin, clas12_tmax)
    for scenario in SCENARIOS:
        d = curves[(curves.scenario==scenario) &
                   (curves.luminosity_factor==10)].sort_values("t_abs")
        ax.fill_between(d.t_abs, d.d1Q-d.sigma_d1Q, d.d1Q+d.sigma_d1Q, alpha=.16)
        ax.plot(d.t_abs, d.d1Q, label=labels[scenario])
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$d_1^Q(t)$")
    ax.set_title("What limits the 10x D-term projection?")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures/"03_d1_systematics_10x.png")

    # 04: show how systematics alter the luminosity progression of alpha.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    for scenario, label in labels.items():
        d = pars[pars.scenario==scenario].sort_values("luminosity_factor")
        ax.plot(d.luminosity_factor, 100*d.sigma_dalpha/abs(d.dalpha),
                marker="o", label=label)
    ax.set_xscale("log")
    ax.set_xticks(LUMI_FACTORS, [f"{L}x" for L in LUMI_FACTORS])
    ax.set_xlabel("Luminosity relative to pass-2 exposure")
    ax.set_ylabel(r"Relative uncertainty on $\alpha$ (%)")
    ax.set_title(r"Luminosity and systematic control of D-term curvature")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures/"04_alpha_precision_systematics.png")

    # 05: baseline luminosity bands, useful after the systematics caveat.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    add_t_range_context(ax, clas12_tmin, clas12_tmax)
    for L in LUMI_FACTORS:
        d = curves[(curves.scenario=="baseline") &
                   (curves.luminosity_factor==L)].sort_values("t_abs")
        ax.fill_between(d.t_abs, d.d1Q-d.sigma_d1Q, d.d1Q+d.sigma_d1Q, alpha=.16)
        ax.plot(d.t_abs, d.d1Q, label=f"{L}x")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$d_1^Q(t)$")
    ax.set_title("D-term projection with current point-to-point systematics")
    ax.grid(alpha=.2); ax.legend(title="Pass-2 exposure")
    savefig(fig, figures/"05_d1_baseline_luminosity.png")

def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default=str(here/"output"/"stage2"/"tables"))
    ap.add_argument("--outdir", default=str(here/"output"/"stage3c"))
    ap.add_argument("--force-derivatives", action="store_true")
    ap.add_argument(
        "--workers", type=int, default=min(8, os.cpu_count() or 1),
        help="Parallel KM15 parameter workers (default: min(8, available CPUs))."
    )
    args = ap.parse_args()

    inp = Path(args.input_dir).resolve()
    # Transitional compatibility with the pre-standardization Stage-2 output.
    # New runs should use output/stage2; this fallback can be removed once the
    # Stage-2 producer is updated.
    if not inp.exists():
        legacy = here/"output_stage2_pass2"/"tables"
        if legacy.exists():
            print(f"[input] standardized Stage-2 directory not found; using legacy {legacy}")
            inp = legacy.resolve()
    out = Path(args.outdir).resolve()
    tables = out/"tables"; figures = out/"figures"
    tables.mkdir(parents=True, exist_ok=True); figures.mkdir(parents=True, exist_ok=True)

    base1 = pd.read_csv(inp/"joint_fit_input_km15_1x.csv")
    base1["bin"] = base1["bin"].astype(int)
    clas12_tmin = float(base1["t_abs"].min())
    clas12_tmax = float(base1["t_abs"].max())

    # The installed KM15 model is the authority for central parameter values.
    projection_model = configure_volker_subtraction(copy.deepcopy(th_KM15))
    central = copy.deepcopy(projection_model.parameters)
    candidate_parameters = [
        p for p in TARGET_PARAMETERS + IMH_PARAMETER_CANDIDATES if p in central
    ]

    print("="*96)
    print("STAGE 3C: VOLKER/BURKERT GENERALIZED-MULTIPOLE D-TERM PROJECTION")
    print("="*96)
    print("Gepard relation : ReH = DR[ImH] - S(t)")
    print("subtraction     : S(t)=C*(1-t/M^2)^(-alpha)")
    print(f"d1^Q(0) truth   : {-0.9*central['C']:.10g}")
    print(f"C_H(0) truth    : {-central['C']:.10g}")
    print(f"M^2 truth       : {central['dM2']:.10g} GeV^2")
    print(f"alpha truth     : {central['dalpha']:.10g}")
    print("candidate pars  :", ", ".join(candidate_parameters))
    print("important       : ImH parameters propagate into ReH through Gepard's DR")
    print("pseudo-data     : KM15 ImH + Volker-style D-term truth; unpublished pass-2 central values are not used")
    print(f"CLAS6 context   : historical input range |t|={CLAS6_DTERM_TMIN:.2f}"
          f"--{CLAS6_DTERM_TMAX:.2f} GeV^2 (plot context only)")
    print(f"CLAS12 context  : pass-2 controlled range |t|={clas12_tmin:.3f}"
          f"--{clas12_tmax:.3f} GeV^2")

    cache = tables/"volker_observable_parameter_derivatives.csv"
    if args.force_derivatives and cache.exists():
        cache.unlink()
    deriv, steps, stability = finite_derivatives(
        base1, projection_model, candidate_parameters, cache,
        workers=max(1, args.workers), compute_stability=True
    )

    parameters, activity = classify_active_parameters(
        deriv, steps, candidate_parameters
    )
    activity.to_csv(tables/"parameter_activity.csv", index=False)

    print("\nParameter activity:")
    print(activity[[
        "parameter","active","activity_response",
        "step_fraction_of_parameter_scale","reason"
    ]].to_string(index=False, float_format=lambda x:f"{x:.3g}"))
    print("\nACTIVE PARAMETERS BEFORE DEGENERACY CHECK:", ", ".join(parameters))

    parameters, degeneracy = remove_exact_response_degeneracies(
        deriv, parameters
    )
    degeneracy.to_csv(tables/"parameter_degeneracies.csv", index=False)
    print("\nExact observable-response degeneracy check:")
    print(degeneracy.to_string(index=False))
    print("\nFINAL INDEPENDENT FIT PARAMETERS:", ", ".join(parameters))

    # The parallel derivative pass already computed h/2 and 2h.  Report only
    # the parameters that survive activity + exact-degeneracy filtering.
    stability = stability[stability["parameter"].isin(parameters)].copy()
    stability.to_csv(tables/"finite_difference_stability.csv", index=False)
    print("\nFinite-step stability (active parameters):")
    print(stability.to_string(index=False, float_format=lambda x:f"{x:.3g}"))

    all_pars, all_curves, all_diag = [], [], []
    tgrid = np.linspace(0.0, 0.95, 191)

    for L in LUMI_FACTORS:
        data = pd.read_csv(inp/f"joint_fit_input_km15_{L}x.csv")
        data["bin"] = data["bin"].astype(int)
        for scenario in SCENARIOS:
            q = scenario_input(data, scenario)
            cov, layout, diag = _validated_covariance(
                q, deriv, parameters, central,
                label=f"{scenario} {L}x"
            )
            idx = {p:i for i,p in enumerate(layout)}

            C = central["C"]; m = central["dM2"]; alpha = central["dalpha"]
            sC = math.sqrt(max(cov[idx["C"],idx["C"]],0))
            sm = math.sqrt(max(cov[idx["dM2"],idx["dM2"]],0))
            sa = math.sqrt(max(cov[idx["dalpha"],idx["dalpha"]],0))
            corr_C_M = cov[idx["C"],idx["dM2"]] / max(sC*sm,1e-300)
            corr_C_a = cov[idx["C"],idx["dalpha"]] / max(sC*sa,1e-300)
            corr_M_a = cov[idx["dM2"],idx["dalpha"]] / max(sm*sa,1e-300)

            all_pars.append({
                "scenario":scenario, "luminosity_factor":L,
                "C":C, "sigma_C":sC,
                "CH0":-C, "sigma_CH0":sC,
                "d1Q0":-0.9*C, "sigma_d1Q0":0.9*sC,
                "dM2":m, "sigma_dM2":sm,
                "dalpha":alpha, "sigma_dalpha":sa,
                "corr_C_dM2":corr_C_M,
                "corr_C_dalpha":corr_C_a,
                "corr_dM2_dalpha":corr_M_a,
            })
            curve = propagate_curve(cov, layout, central, tgrid, scenario, L)
            all_curves.append(curve)
            diag.update({"scenario":scenario,"luminosity_factor":L})
            all_diag.append(diag)

    pars = pd.DataFrame(all_pars)
    curves = pd.concat(all_curves, ignore_index=True)
    diags = pd.DataFrame(all_diag)
    pars.to_csv(tables/"dr_parameter_uncertainties.csv", index=False)
    curves.to_csv(tables/"ch_d1_bands.csv", index=False)
    diags.to_csv(tables/"global_fit_diagnostics.csv", index=False)
    make_plots(curves, pars, figures, clas12_tmin, clas12_tmax)

    print("\nC_H(0), d1^Q(0), and t-shape precision:")
    show = pars[["scenario","luminosity_factor","CH0","sigma_CH0",
                  "d1Q0","sigma_d1Q0","dM2","sigma_dM2","dalpha","sigma_dalpha",
                 "corr_C_dM2","corr_C_dalpha","corr_dM2_dalpha"]]
    print(show.to_string(index=False, float_format=lambda x:f"{x:.4g}"))

    # Explicit high-|t| shape uncertainty at the existing clean endpoint.
    target_t = 0.900991
    high = curves.iloc[(curves.t_abs-target_t).abs().argsort()].copy()
    high = high.groupby(["scenario","luminosity_factor"],as_index=False).first()
    print("\nC_H(t) uncertainty near the existing |t|=0.901 GeV^2 endpoint:")
    print(high[["scenario","luminosity_factor","t_abs","CH","sigma_CH"]]
          .sort_values(["scenario","luminosity_factor"])
          .to_string(index=False,float_format=lambda x:f"{x:.4g}"))

    print("\nValidated global-fit diagnostics:")
    print(diags[[
        "scenario","luminosity_factor","n_parameters","rank",
        "rank_deficit","condition_number_effective"
    ]].to_string(index=False, float_format=lambda x:f"{x:.4g}"))

    print("\n[output]", out)
    print("[next] Inspect the statistical-only luminosity progression first,")
    print("       then quantify how much current and improved PTP systematics")
    print("       reduce that high-luminosity potential.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
