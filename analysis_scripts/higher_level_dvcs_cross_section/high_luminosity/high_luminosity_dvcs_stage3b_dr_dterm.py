#!/usr/bin/env python3
"""
Stage 3B: genuine KM15 fixed-t dispersion-relation projection.

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
  output_stage3b_dr_dterm/
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
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import gepard as g
from gepard.fits import th_KM15


LUMI_FACTORS = (1, 2, 5, 10)
SCENARIOS = ("baseline", "ptp_half", "statistics_only")

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


def finite_derivatives(base: pd.DataFrame, model, parameters: List[str],
                       cache_csv: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if cache_csv.exists():
        d = pd.read_csv(cache_csv)
        needed = {"point_id", "xs0", "bsa0"}
        for p in parameters:
            needed |= {f"d_xs_d_{p}", f"d_bsa_d_{p}"}
        if needed.issubset(d.columns) and len(d) == len(base):
            steps_path = cache_csv.with_name("finite_difference_steps.csv")
            steps = pd.read_csv(steps_path) if steps_path.exists() else pd.DataFrame()
            print(f"[cache] using {cache_csv}")
            return d, steps

    central = copy.deepcopy(model.parameters)
    limits = getattr(model, "parameters_limits", {})
    steps_rows = []
    out = base[["point_id", "bin", "xB", "Q2", "t_abs", "phi_deg"]].copy()

    print(f"[derivatives] central predictions for {len(base)} points")
    xs0, a0 = [], []
    for i, row in enumerate(base.itertuples(index=False), 1):
        x, a = predict_xs_bsa(model, row)
        xs0.append(x); a0.append(a)
        if i % 100 == 0:
            print(f"  central {i}/{len(base)}")
    out["xs0"] = xs0
    out["bsa0"] = a0

    for ip, par in enumerate(parameters, 1):
        val = float(central[par])
        h = choose_step(par, val, limits)
        print(f"[derivatives] {ip}/{len(parameters)} {par}: value={val:.8g}, step={h:.4g}")

        plus_xs, plus_a, minus_xs, minus_a = [], [], [], []

        model.parameters[par] = val + h
        for row in base.itertuples(index=False):
            x, a = predict_xs_bsa(model, row)
            plus_xs.append(x); plus_a.append(a)

        model.parameters[par] = val - h
        for row in base.itertuples(index=False):
            x, a = predict_xs_bsa(model, row)
            minus_xs.append(x); minus_a.append(a)

        model.parameters[par] = val
        out[f"d_xs_d_{par}"] = (
            np.asarray(plus_xs) - np.asarray(minus_xs)
        ) / (2*h)
        out[f"d_bsa_d_{par}"] = (
            np.asarray(plus_a) - np.asarray(minus_a)
        ) / (2*h)

        # A simple finite-step stability indicator: response relative to central.
        rel_x = np.nanmedian(
            np.abs(np.asarray(plus_xs) - np.asarray(minus_xs)) /
            np.maximum(2*np.abs(np.asarray(xs0)), 1e-30)
        )
        abs_a = np.nanmedian(
            np.abs(np.asarray(plus_a) - np.asarray(minus_a)) / 2
        )
        steps_rows.append({
            "parameter": par, "central_value": val, "step": h,
            "median_fractional_XS_half_response": rel_x,
            "median_absolute_BSA_half_response": abs_a,
        })

    model.parameters.update(central)
    out.to_csv(cache_csv, index=False)
    steps = pd.DataFrame(steps_rows)
    steps.to_csv(cache_csv.with_name("finite_difference_steps.csv"), index=False)
    return out, steps


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


def make_plots(curves, pars, figures):
    # CH(t), baseline luminosity evolution.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
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
    ap.add_argument("--outdir", default=str(here/"output_stage3b_dr_dterm"))
    ap.add_argument("--force-derivatives", action="store_true")
    args = ap.parse_args()

    inp = Path(args.input_dir).resolve()
    out = Path(args.outdir).resolve()
    tables = out/"tables"; figures = out/"figures"
    tables.mkdir(parents=True, exist_ok=True); figures.mkdir(parents=True, exist_ok=True)

    base1 = pd.read_csv(inp/"joint_fit_input_km15_1x.csv")
    base1["bin"] = base1["bin"].astype(int)

    # The installed KM15 model is the authority for central parameter values.
    central = copy.deepcopy(th_KM15.parameters)
    parameters = [p for p in TARGET_PARAMETERS + IMH_PARAMETER_CANDIDATES if p in central]

    print("="*96)
    print("STAGE 3B: GENUINE KM15 FIXED-t DISPERSION-RELATION PROJECTION")
    print("="*96)
    print("Gepard relation : ReH = DR[ImH] - C/(1-t/mC2)^2")
    print(f"KM15 C          : {central['C']:.10g}")
    print(f"KM15 mC2        : {central['mC2']:.10g} GeV^2")
    print(f"Volker C_H(0)   : {-central['C']:.10g}")
    print(f"Volker d1^Q(0)  : {-0.9*central['C']:.10g}")
    print("fit parameters  :", ", ".join(parameters))
    print("important       : ImH parameters propagate into ReH through Gepard's DR")
    print("pseudo-data     : KM15 only; unpublished pass-2 central values are not used")

    cache = tables/"km15_observable_parameter_derivatives.csv"
    if args.force_derivatives and cache.exists():
        cache.unlink()
    deriv, steps = finite_derivatives(base1, th_KM15, parameters, cache)

    all_pars, all_curves, all_diag = [], [], []
    tgrid = np.linspace(0.0, 0.95, 191)

    for L in LUMI_FACTORS:
        data = pd.read_csv(inp/f"joint_fit_input_km15_{L}x.csv")
        data["bin"] = data["bin"].astype(int)
        for scenario in SCENARIOS:
            q = scenario_input(data, scenario)
            cov, layout, diag = covariance(q, deriv, parameters, central)
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

    print("\n[output]", out)
    print("[next] Inspect rank/conditioning and finite-step responses before")
    print("       interpreting pressure/shear.  If stable, propagate the fitted")
    print("       d1^Q(t) covariance through the spatial Fourier transform.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
