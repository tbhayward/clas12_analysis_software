#!/usr/bin/env python3
"""
Validate the AAOgen/VPK pi0 model against published Hall-A Dlamini et al. data.

This is a literal Python transcription of the March-2021 Valery Kubarovsky
parameterization in AAOgen aao_norad/dvmpx.F for pi0 production.

Run from:
  /u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/external_pi0_normalization/

Typical:
  python3 validate_aao_pi0_model.py

Input:
  import/dlamini2021_halla_pi0_structure_functions.csv

Outputs:
  output/dlamini2021_vpk_point_comparison.csv
  output/dlamini2021_vpk_setting_summary.csv
  output/png/dlamini2021_sigmaU_vs_tprime.png
  output/png/dlamini2021_data_over_vpk_vs_Q2.png
  output/png/dlamini2021_data_over_vpk_all_points.png
  output/png/dlamini2021_interference_checks.png
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MP = 0.938272
MPI0 = 0.134976
ALPHA = 0.00729927
HC2 = 389379.36  # GeV^-2 -> nb

# Exact pi0 parameters from AAOgen aao_norad/dvmpx.F XSINIT, March 2021.
P = np.array([
    6.2054, 0.8020, -0.1066, 1.8364, 0.0,
    94.6474, 3.4426, -1.9769, -0.1324, 0.0,
    17.0001, 2.1228
], dtype=float)


def epsilon(x: float, q2: float, ebeam: float) -> float:
    y = q2 / (2.0 * MP * x * ebeam)
    e1 = (y * x * MP) ** 2 / q2
    return (1.0 - y - e1) / (1.0 - y + y * y / 2.0 + e1)


def tminq(q2: float, x: float) -> float:
    """Literal AAOgen tminq: returns -t_min = |t_min| (>0), in GeV^2."""
    if not (0.0 < x < 1.0):
        return float("nan")
    w2 = q2 * (1.0 / x - 1.0) + MP * MP
    if w2 <= 0:
        return float("nan")
    w = math.sqrt(w2)
    if w < MP + MPI0:
        return float("nan")

    e1cm = (w2 + q2 + MP * MP) / (2.0 * w)
    p1sq = e1cm * e1cm - MP * MP
    e3cm = (w2 - MPI0 * MPI0 + MP * MP) / (2.0 * w)
    p3sq = e3cm * e3cm - MP * MP
    if p1sq < 0 or p3sq < 0:
        return float("nan")
    p1cm = math.sqrt(p1sq)
    p3cm = math.sqrt(p3sq)
    return -(((q2 + MPI0 * MPI0) ** 2) / (4.0 * w2) - (p1cm - p3cm) ** 2)


def phase(q2: float, x: float) -> float:
    if x >= 1.0:
        return 0.0
    w2 = q2 * (1.0 / x - 1.0) + MP * MP
    w4 = w2 * w2
    mp2 = MP * MP
    mp4 = mp2 * mp2
    q4 = q2 * q2
    lam = w4 + q4 + mp4 + 2.0*w2*q2 - 2.0*w2*mp2 + 2.0*q2*mp2
    if lam <= 0:
        return 0.0
    return 16.0 * math.pi * (w2 - mp2) * math.sqrt(lam)


def ht(t: float, x: float, q2: float) -> float:
    T = -t
    slope = P[1] + P[2] * (math.log(x) - math.log(0.15))
    return P[0] * math.exp(-slope * T) * q2 ** (P[3] / 2.0)


def et(t: float, x: float, q2: float) -> float:
    T = -t
    slope = P[6] + P[7] * (math.log(x) - math.log(0.15))
    return P[5] * math.exp(-slope * T) * q2 ** (P[8] / 2.0)


def htebar(t: float, x: float, q2: float) -> float:
    T = -t
    return P[10] * math.exp(-P[11] * T)


def xcheck_kine(t: float, x: float, q2: float, ebeam: float) -> bool:
    # Literal physics checks from XCHECK_KINE, simplified algebraically only
    # where no numerical convention changes.
    if q2 <= 0 or ebeam <= 0:
        return False
    xmin = q2 / (2.0 * MP * ebeam)
    if x <= xmin or x > 1.0:
        return False
    w2 = MP*MP + q2*(1.0/x - 1.0)
    if w2 <= (MP + MPI0)**2:
        return False

    nu = q2 / (2.0 * MP * x)
    w = math.sqrt(w2)
    qmod = math.sqrt(nu*nu + q2)
    e1cm = MP*(MP + nu)/w
    p1cm = MP*qmod/w
    e2cm = (w2 + MP*MP - MPI0*MPI0)/(2.0*w)
    if e2cm <= MP:
        return False
    p2cm = math.sqrt(max(0.0, e2cm*e2cm - MP*MP))
    tmax_neg = 2.0*(MP*MP - e1cm*e2cm - p1cm*p2cm)
    tmin_neg = 2.0*(MP*MP - e1cm*e2cm + p1cm*p2cm)
    if t >= tmin_neg or t <= tmax_neg:
        return False

    eps = epsilon(x, q2, ebeam)
    if not (0.0 < eps < 1.0):
        return False
    ksi = x/(2.0-x)*(1.0 + MP*MP/q2)
    return ksi < 1.0


def xsigma_t(t: float, x: float, q2: float, ebeam: float) -> float:
    if not xcheck_kine(t, x, q2, ebeam):
        return 0.0
    t0 = tminq(q2, x)
    T = -t
    ksi = x/(2.0-x)*(1.0 + MP*MP/q2)
    val = (4.0*math.pi*ALPHA/2.0/phase(q2,x)/q2**2 *
           ((1.0-ksi**2)*ht(t,x,q2)**2 +
            (T-t0)/8.0/MP**2*et(t,x,q2)**2))
    return val * HC2


def xsigma_tt(t: float, x: float, q2: float, ebeam: float) -> float:
    if not xcheck_kine(t, x, q2, ebeam):
        return 0.0
    t0 = tminq(q2, x)
    T = -t
    val = (-4.0*math.pi*ALPHA/2.0/phase(q2,x)/q2**2 *
           (T-t0)/8.0/MP**2*et(t,x,q2)**2)
    return val * HC2


def xsigma_l(t: float, x: float, q2: float, ebeam: float) -> float:
    return 0.0


def xsigma_lt(t: float, x: float, q2: float, ebeam: float) -> float:
    if not xcheck_kine(t, x, q2, ebeam):
        return 0.0
    t0 = tminq(q2, x)
    T = -t
    ksi = x/(2.0-x)*(1.0 + MP*MP/q2)
    rad = T - t0
    if rad < 0:
        return 0.0
    val = (4.0*math.pi*ALPHA/math.sqrt(2.0)/phase(q2,x)/q2**1.5 *
           ksi*math.sqrt(max(0.0, 1.0-ksi**2)) *
           math.sqrt(rad)/2.0/MP * htebar(t,x,q2)**2)
    return val * HC2


def evaluate_row(row: pd.Series) -> dict:
    x = float(row["xB_mean"])
    q2 = float(row["Q2_GeV2"])
    e = float(row["Ebeam_GeV"])
    tp = float(row["tprime_mean_GeV2"])

    tm = tminq(q2, x)
    # AAOgen tminq returns -t_min = |t_min|. Dlamini t' = t_min - t,
    # so -t = (-t_min) + t' and the signed AAOgen del2 is:
    t = -(tm + tp)
    eps = epsilon(x, q2, e)
    st = xsigma_t(t, x, q2, e)
    sl = xsigma_l(t, x, q2, e)
    su = st + eps*sl
    stt = xsigma_tt(t, x, q2, e)
    slt = xsigma_lt(t, x, q2, e)

    return {
        "vpk_tmin_abs_GeV2": tm,
        "vpk_minus_tmin_GeV2": tm,
        "vpk_t_GeV2": t,
        "vpk_minus_t_GeV2": -t,
        "vpk_epsilon": eps,
        "vpk_sigmaT_nb_per_GeV2": st,
        "vpk_sigmaL_nb_per_GeV2": sl,
        "vpk_sigmaU_nb_per_GeV2": su,
        "vpk_sigmaTT_nb_per_GeV2": stt,
        "vpk_sigmaLT_nb_per_GeV2": slt,
        "vpk_sigmaLTprime_nb_per_GeV2": 0.0,
        "vpk_kine_valid": xcheck_kine(t, x, q2, e),
    }


def add_ratios(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    pairs = [
        ("U", "sigmaU", "vpk_sigmaU"),
        ("TT", "sigmaTT", "vpk_sigmaTT"),
        ("LT", "sigmaLT", "vpk_sigmaLT"),
    ]
    for tag, data, model in pairs:
        d = out[f"{data}_nb_per_GeV2"].astype(float)
        stat = out[f"{data}_stat_nb_per_GeV2"].astype(float)
        syst = out[f"{data}_syst_nb_per_GeV2"].astype(float)
        m = out[f"{model}_nb_per_GeV2"].astype(float)
        out[f"sigma{tag}_totalerr_nb_per_GeV2"] = np.hypot(stat, syst)
        out[f"data_over_vpk_{tag}"] = np.where(m != 0, d/m, np.nan)
        out[f"data_over_vpk_{tag}_err"] = np.where(
            m != 0, out[f"sigma{tag}_totalerr_nb_per_GeV2"]/np.abs(m), np.nan)
    return out


def weighted_mean(values, errors):
    values = np.asarray(values, float)
    errors = np.asarray(errors, float)
    mask = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    if not np.any(mask):
        return np.nan, np.nan, 0
    w = 1.0/errors[mask]**2
    return np.sum(w*values[mask])/np.sum(w), math.sqrt(1.0/np.sum(w)), int(mask.sum())


def make_setting_summary(df):
    rows = []
    keys = ["xB_mean", "Q2_GeV2", "Ebeam_GeV"]
    for key, g in df.groupby(keys, sort=True):
        r, er, n = weighted_mean(g["data_over_vpk_U"], g["data_over_vpk_U_err"])
        rows.append({
            "xB_mean": key[0], "Q2_GeV2": key[1], "Ebeam_GeV": key[2],
            "n_tprime_points": len(g),
            "weighted_data_over_vpk_U": r,
            "weighted_data_over_vpk_U_err": er,
            "n_used": n,
            "mean_vpk_epsilon": g["vpk_epsilon"].mean(),
            "paper_epsilon": g["epsilon"].mean(),
        })
    return pd.DataFrame(rows)


def plot_sigma_u(df, outdir):
    fig, ax = plt.subplots(figsize=(9, 6))
    for (x, q2), g in df.groupby(["xB_mean", "Q2_GeV2"], sort=True):
        g = g.sort_values("tprime_mean_GeV2")
        err = np.hypot(g["sigmaU_stat_nb_per_GeV2"], g["sigmaU_syst_nb_per_GeV2"])
        ax.errorbar(g["tprime_mean_GeV2"], g["sigmaU_nb_per_GeV2"], yerr=err,
                    fmt="o", capsize=2, label=f"data xB={x:.2f}, Q²={q2:.2f}")
        ax.plot(g["tprime_mean_GeV2"], g["vpk_sigmaU_nb_per_GeV2"], marker=".",
                linestyle="--", label=f"VPK xB={x:.2f}, Q²={q2:.2f}")
    ax.set_xlabel("t' = tmin - t (GeV²)")
    ax.set_ylabel("dσU/dt (nb/GeV²)")
    ax.set_yscale("log")
    ax.set_title("Dlamini 2021 vs AAOgen/VPK: unseparated π⁰ cross section")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir/"dlamini2021_sigmaU_vs_tprime.png", dpi=220)
    plt.close(fig)


def plot_setting_ratio(summary, outdir):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for x, g in summary.groupby("xB_mean", sort=True):
        ax.errorbar(g["Q2_GeV2"], g["weighted_data_over_vpk_U"],
                    yerr=g["weighted_data_over_vpk_U_err"],
                    fmt="o-", capsize=3, label=f"xB≈{x:.2f}")
    ax.axhline(1.0, linewidth=1)
    ax.set_xlabel("Q² (GeV²)")
    ax.set_ylabel("Dlamini / VPK for dσU/dt")
    ax.set_title("AAOgen/VPK normalization check by Hall-A setting")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir/"dlamini2021_data_over_vpk_vs_Q2.png", dpi=220)
    plt.close(fig)


def plot_all_ratios(df, outdir):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    order = df.sort_values(["xB_mean","Q2_GeV2","tprime_mean_GeV2"]).reset_index(drop=True)
    xidx = np.arange(len(order))
    ax.errorbar(xidx, order["data_over_vpk_U"], yerr=order["data_over_vpk_U_err"],
                fmt="o", capsize=2)
    ax.axhline(1.0, linewidth=1)
    ax.set_xlabel("Published point index (sorted by xB, Q², t')")
    ax.set_ylabel("Dlamini / VPK for dσU/dt")
    ax.set_title("Point-by-point AAOgen/VPK normalization residual")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir/"dlamini2021_data_over_vpk_all_points.png", dpi=220)
    plt.close(fig)


def plot_interference(df, outdir):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    order = df.sort_values(["xB_mean","Q2_GeV2","tprime_mean_GeV2"]).reset_index(drop=True)
    xidx = np.arange(len(order))
    ax.plot(xidx, order["sigmaTT_nb_per_GeV2"], "o", label="Dlamini TT")
    ax.plot(xidx, order["vpk_sigmaTT_nb_per_GeV2"], ".", label="VPK TT")
    ax.plot(xidx, order["sigmaLT_nb_per_GeV2"], "s", label="Dlamini LT")
    ax.plot(xidx, order["vpk_sigmaLT_nb_per_GeV2"], ".", label="VPK LT")
    ax.axhline(0.0, linewidth=1)
    ax.set_xlabel("Published point index (sorted by xB, Q², t')")
    ax.set_ylabel("Structure function (nb/GeV²)")
    ax.set_title("Dlamini interference structure functions vs AAOgen/VPK")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir/"dlamini2021_interference_checks.png", dpi=220)
    plt.close(fig)


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path,
                    default=here/"import"/"dlamini2021_halla_pi0_structure_functions.csv")
    ap.add_argument("--output-dir", type=Path, default=here/"output")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    required = {
        "xB_mean","Q2_GeV2","Ebeam_GeV","minus_tmin_GeV2","epsilon",
        "tprime_mean_GeV2",
        "sigmaU_nb_per_GeV2","sigmaU_stat_nb_per_GeV2","sigmaU_syst_nb_per_GeV2",
        "sigmaTT_nb_per_GeV2","sigmaTT_stat_nb_per_GeV2","sigmaTT_syst_nb_per_GeV2",
        "sigmaLT_nb_per_GeV2","sigmaLT_stat_nb_per_GeV2","sigmaLT_syst_nb_per_GeV2",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit("Missing required columns: " + ", ".join(missing))

    calc = pd.DataFrame([evaluate_row(r) for _, r in df.iterrows()])
    out = pd.concat([df.reset_index(drop=True), calc], axis=1)
    out = add_ratios(out)

    # Explicit convention/implementation checks.
    out["delta_minus_tmin_GeV2"] = out["vpk_minus_tmin_GeV2"] - out["minus_tmin_GeV2"]
    out["delta_epsilon"] = out["vpk_epsilon"] - out["epsilon"]

    outdir = args.output_dir
    pngdir = outdir/"png"
    pngdir.mkdir(parents=True, exist_ok=True)

    summary = make_setting_summary(out)
    out.to_csv(outdir/"dlamini2021_vpk_point_comparison.csv", index=False)
    summary.to_csv(outdir/"dlamini2021_vpk_setting_summary.csv", index=False)

    plot_sigma_u(out, pngdir)
    plot_setting_ratio(summary, pngdir)
    plot_all_ratios(out, pngdir)
    plot_interference(out, pngdir)

    print("\n=== Dlamini 2021 / AAOgen-VPK validation ===")
    print(f"Input points: {len(out)}")
    print(f"Valid VPK kinematics: {int(out['vpk_kine_valid'].sum())}/{len(out)}")
    print(f"max |AAOgen - paper epsilon| = {out['delta_epsilon'].abs().max():.5f}")
    print(f"max |AAOgen - paper (-tmin)| = {out['delta_minus_tmin_GeV2'].abs().max():.5f} GeV^2")
    print("\nPer-setting weighted Dlamini/VPK ratio for sigma_U:")
    print(summary[["xB_mean","Q2_GeV2","Ebeam_GeV",
                   "weighted_data_over_vpk_U","weighted_data_over_vpk_U_err",
                   "n_tprime_points"]].to_string(index=False))
    print(f"\nWrote: {outdir/'dlamini2021_vpk_point_comparison.csv'}")
    print(f"Wrote: {outdir/'dlamini2021_vpk_setting_summary.csv'}")
    print(f"Wrote PNGs under: {pngdir}")


if __name__ == "__main__":
    main()
