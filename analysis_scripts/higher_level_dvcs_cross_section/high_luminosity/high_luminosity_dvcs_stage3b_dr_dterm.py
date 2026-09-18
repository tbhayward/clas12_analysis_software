#!/usr/bin/env python3
"""
Stage 3B v3: validated KM15 fixed-t dispersion-relation projection.

This is the first stage in which the subtraction function is fitted *inside*
the same dispersion-relation model that generates ImH and ReH.

KM15/Gepard convention
----------------------
Gepard evaluates

    ReH(xi,t) = DR[ImH](xi,t) - S(t)

with

    S(t) = C / (1 - t/mC2)^2.

Therefore, in the Volker/Burkert-Elouadrhiri-Girod convention used for the
mechanical-structure presentation,

    C_H(t) = -S(t) = -C / (1 - t/mC2)^2
    d1^Q(t) = (9/10) C_H(t).

For the installed KM15 values supplied by the user:
    C   = 2.7678681812890016
    mC2 = 1.4497411858248308 GeV^2

so KM15 corresponds to C_H(0)=-2.767868... and d1^Q(0)=-2.491081....

What is fitted
--------------
The fit parameters are:
  * C and mC2: normalization and t-shape of the subtraction function;
  * all available KM15 parameters that enter its ImH ansatz:
      Nsea, alS, alpS, mS2, rS, bS,
      Nv,   alv, alpv, mv2, rv,  bv
    (only names actually present in the installed KM15 are used);
  * one common XS normalization nuisance;
  * one common BSA beam-polarization nuisance.

This is crucial: when an ImH parameter is varied, Gepard recomputes ReH through
the principal-value dispersion integral.  ReH is NOT floated independently.

The covariance is obtained from finite changes of the actual XS and BSA
predictions around KM15.  Pseudo-data equal KM15, so this is a projected
sensitivity/covariance study, not a fit to unpublished pass-2 central values.

Scenarios:
  baseline        current PTP + 10% XS normalization + 4% BSA polarization
  ptp_half        XS and BSA point-to-point systematics divided by two
  statistics_only statistical errors only

Outputs:
  output_stage3b_dr_dterm_v3_parallel/
    tables/dr_parameter_uncertainties.csv
    tables/ch_d1_bands.csv
    tables/global_fit_diagnostics.csv
    tables/finite_difference_steps.csv
    figures/01_CH_t_luminosity.png
    figures/02_d1_t_luminosity.png
    figures/03_CH_t_systematics_10x.png
    figures/04_mC2_shape_precision.png
"""

from __future__ import annotations

import argparse
import copy
import math
import os
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

# Parameters entering DispersionFixedPoleCFF.ImH in Gepard.
IMH_PARAMETER_CANDIDATES = (
    "Nsea", "alS", "alpS", "mS2", "rS", "bS",
    "Nv", "alv", "alpv", "mv2", "rv", "bv",
)
TARGET_PARAMETERS = ("C", "mC2")


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
    if name in ("mC2", "mS2", "mv2"):
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


def _predict_records(model, records):
    xs, bsa = [], []
    for rec in records:
        x, a = predict_xs_bsa(model, _Row(rec))
        xs.append(x)
        bsa.append(a)
    return np.asarray(xs), np.asarray(bsa)


def _parameter_worker(payload):
    """
    One process owns one independent KM15 copy and computes all finite-step
    information for one parameter.  Parallelizing by parameter is deliberate:
    Gepard model mutation is never shared between processes, and each worker
    gets a large enough block of work that process overhead is negligible.
    """
    par, value, h0, records, do_stability = payload
    model = copy.deepcopy(th_KM15)

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
    model = copy.deepcopy(th_KM15)
    return _predict_records(model, records)


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
    mandatory=("C", "mC2"),
    response_threshold: float = 1e-8,
    boundary_step_fraction: float = 1e-4,
) -> Tuple[List[str], pd.DataFrame]:
    """
    Keep only genuine local directions of the installed KM15 solution.

    A candidate ImH parameter is excluded when:
      1) its finite change produces essentially zero XS and BSA response, or
      2) the available symmetric step has collapsed against a parameter
         boundary (the bv situation in KM15).

    C and mC2 are mandatory because they define the subtraction function.
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


def ch_and_gradient(t_abs: np.ndarray, C: float, mC2: float):
    """C_H(t)=-C/(1+|t|/mC2)^2 and derivatives wrt C,mC2."""
    u = 1.0 + t_abs/mC2
    ch = -C/u**2
    dC = -1.0/u**2
    dm = -C * 2.0 * t_abs / (mC2**2 * u**3)
    return ch, dC, dm


def propagate_curve(cov, layout, central, tgrid, scenario, lumi):
    idx = {p:i for i,p in enumerate(layout)}
    C = central["C"]; m = central["mC2"]
    ch, dC, dm = ch_and_gradient(tgrid, C, m)
    rows = []
    for t, y, gC, gm in zip(tgrid, ch, dC, dm):
        grad = np.zeros(len(layout))
        grad[idx["C"]] = gC
        grad[idx["mC2"]] = gm
        var = float(grad @ cov @ grad)
        sch = math.sqrt(max(var, 0))
        rows.append({
            "scenario": scenario, "luminosity_factor": lumi,
            "t_abs": float(t), "CH": float(y), "sigma_CH": sch,
            "d1Q": 0.9*float(y), "sigma_d1Q": 0.9*sch,
        })
    return pd.DataFrame(rows)


def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def add_clas6_context(ax):
    ax.axvspan(
        CLAS6_DTERM_TMIN, CLAS6_DTERM_TMAX,
        alpha=0.10, zorder=0,
        label=r"CLAS6 D-term analysis input $|t|$ range"
    )


def make_plots(curves, pars, figures):
    # CH(t), baseline luminosity evolution.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    add_clas6_context(ax)
    for L in LUMI_FACTORS:
        d = curves[(curves.scenario=="baseline") &
                   (curves.luminosity_factor==L)].sort_values("t_abs")
        ax.fill_between(d.t_abs, d.CH-d.sigma_CH, d.CH+d.sigma_CH, alpha=0.18)
        ax.plot(d.t_abs, d.CH, label=f"{L}x")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$C_H(t)$")
    ax.set_title(r"DR-constrained subtraction function: luminosity projection")
    ax.grid(alpha=.2); ax.legend(title="Luminosity")
    savefig(fig, figures/"01_CH_t_luminosity.png")

    # d1(t) in Volker convention.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    add_clas6_context(ax)
    for L in LUMI_FACTORS:
        d = curves[(curves.scenario=="baseline") &
                   (curves.luminosity_factor==L)].sort_values("t_abs")
        ax.fill_between(d.t_abs, d.d1Q-d.sigma_d1Q, d.d1Q+d.sigma_d1Q, alpha=0.18)
        ax.plot(d.t_abs, d.d1Q, label=f"{L}x")
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$d_1^Q(t)=0.9\,C_H(t)$")
    ax.set_title(r"Projected quark D-term form factor")
    ax.grid(alpha=.2); ax.legend(title="Luminosity")
    savefig(fig, figures/"02_d1_t_luminosity.png")

    # 10x systematic comparison.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    add_clas6_context(ax)
    labels = {
        "baseline":"Current point-to-point systematics",
        "ptp_half":"Point-to-point systematics / 2",
        "statistics_only":"Statistics only",
    }
    for scenario in SCENARIOS:
        d = curves[(curves.scenario==scenario) &
                   (curves.luminosity_factor==10)].sort_values("t_abs")
        ax.fill_between(d.t_abs, d.CH-d.sigma_CH, d.CH+d.sigma_CH, alpha=.18)
        ax.plot(d.t_abs, d.CH, label=labels[scenario])
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$C_H(t)$")
    ax.set_title(r"Systematic limitation on $C_H(t)$ at 10x")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures/"03_CH_t_systematics_10x.png")

    # mC2 is the direct shape/large-|t| lever-arm parameter.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    for scenario, label in labels.items():
        d = pars[pars.scenario==scenario].sort_values("luminosity_factor")
        ax.plot(d.luminosity_factor, 100*d.sigma_mC2/abs(d.mC2),
                marker="o", label=label)
    ax.set_xscale("log")
    ax.set_xticks(LUMI_FACTORS, [f"{L}x" for L in LUMI_FACTORS])
    ax.set_xlabel("Luminosity relative to pass-2 exposure")
    ax.set_ylabel(r"Relative uncertainty on $m_C^2$ (%)")
    ax.set_title(r"How luminosity constrains the $t$-shape of the D-term")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures/"04_mC2_shape_precision.png")


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default=str(here/"output_stage2_pass2"/"tables"))
    ap.add_argument("--outdir", default=str(here/"output_stage3b_dr_dterm_v3_parallel"))
    ap.add_argument("--force-derivatives", action="store_true")
    ap.add_argument(
        "--workers", type=int, default=min(8, os.cpu_count() or 1),
        help="Parallel KM15 parameter workers (default: min(8, available CPUs))."
    )
    args = ap.parse_args()

    inp = Path(args.input_dir).resolve()
    out = Path(args.outdir).resolve()
    tables = out/"tables"; figures = out/"figures"
    tables.mkdir(parents=True, exist_ok=True); figures.mkdir(parents=True, exist_ok=True)

    base1 = pd.read_csv(inp/"joint_fit_input_km15_1x.csv")
    base1["bin"] = base1["bin"].astype(int)

    # The installed KM15 model is the authority for central parameter values.
    central = copy.deepcopy(th_KM15.parameters)
    candidate_parameters = [
        p for p in TARGET_PARAMETERS + IMH_PARAMETER_CANDIDATES if p in central
    ]

    print("="*96)
    print("STAGE 3B v3: VALIDATED KM15 FIXED-t DISPERSION-RELATION PROJECTION")
    print("="*96)
    print("Gepard relation : ReH = DR[ImH] - C/(1-t/mC2)^2")
    print(f"KM15 C          : {central['C']:.10g}")
    print(f"KM15 mC2        : {central['mC2']:.10g} GeV^2")
    print(f"Volker C_H(0)   : {-central['C']:.10g}")
    print(f"Volker d1^Q(0)  : {-0.9*central['C']:.10g}")
    print("candidate pars  :", ", ".join(candidate_parameters))
    print("important       : ImH parameters propagate into ReH through Gepard's DR")
    print("pseudo-data     : KM15 only; unpublished pass-2 central values are not used")
    print(f"CLAS6 context   : historical input range |t|={CLAS6_DTERM_TMIN:.2f}"
          f"--{CLAS6_DTERM_TMAX:.2f} GeV^2 (plot context only)")

    cache = tables/"km15_observable_parameter_derivatives.csv"
    if args.force_derivatives and cache.exists():
        cache.unlink()
    deriv, steps, stability = finite_derivatives(
        base1, th_KM15, candidate_parameters, cache,
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
    print("\nACTIVE FIT PARAMETERS:", ", ".join(parameters))

    # The parallel derivative pass already computed h/2 and 2h.  Report only
    # the parameters that survive the activity filter.
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
            cov, layout, diag = covariance(q, deriv, parameters, central)
            if diag["rank_deficit"] != 0:
                raise RuntimeError(
                    f"Rank-deficient validated fit for {scenario} {L}x: "
                    f"{diag['rank']}/{diag['n_parameters']}. "
                    "Do not interpret D-term uncertainties until resolved."
                )
            idx = {p:i for i,p in enumerate(layout)}

            C = central["C"]; m = central["mC2"]
            sC = math.sqrt(max(cov[idx["C"],idx["C"]],0))
            sm = math.sqrt(max(cov[idx["mC2"],idx["mC2"]],0))
            corr = cov[idx["C"],idx["mC2"]] / max(sC*sm,1e-300)

            all_pars.append({
                "scenario":scenario, "luminosity_factor":L,
                "C":C, "sigma_C":sC,
                "CH0":-C, "sigma_CH0":sC,
                "d1Q0":-0.9*C, "sigma_d1Q0":0.9*sC,
                "mC2":m, "sigma_mC2":sm,
                "corr_C_mC2":corr,
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
    make_plots(curves, pars, figures)

    print("\nC_H(0), d1^Q(0), and t-shape precision:")
    show = pars[["scenario","luminosity_factor","CH0","sigma_CH0",
                 "d1Q0","sigma_d1Q0","mC2","sigma_mC2","corr_C_mC2"]]
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
    print("[next] Inspect rank/conditioning and finite-step responses before")
    print("       interpreting pressure/shear.  If stable, propagate the fitted")
    print("       d1^Q(t) covariance through the spatial Fourier transform.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
