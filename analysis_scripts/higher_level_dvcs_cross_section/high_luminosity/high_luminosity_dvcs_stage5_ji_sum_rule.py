#!/usr/bin/env python3
"""
high_luminosity_dvcs_stage5_ji_sum_rule.py

Stage 5, v1: CLAS12 high-luminosity flavor separation / Ji-sum-rule framework.

This first version deliberately does three things before attempting a global fit:
  1. Loads and audits the preliminary RGB neutron unpolarized cross-section table.
  2. Builds 1x/2x/5x/10x RGB statistical projections while holding quoted systematics fixed.
  3. Implements the DFJK-inspired zero-skewness valence E^u_v,E^d_v model and computes
     B20^{u_v}, B20^{d_v}; it also provides the bookkeeping needed for A20 and J.

The next version will connect these parameters to finite-skewness proton/neutron DVCS
observables through Gepard and add the published RGB neutron BSA.

Physics scope:
  * workshop projection, not a precision global GPD extraction;
  * valence-dominated H/E model;
  * E zeroth moments fixed by proton/neutron anomalous magnetic moments;
  * finite-skewness double-distribution profile will be fixed, not floated, in v2;
  * RGA and RGB unpolarized luminosities scale together in the high-L scenarios;
  * polarized-target RGC/RGH precision is NOT scaled with RGA/RGB luminosity.

Expected input beside this script:
  import/ndvcs_clas12_preliminary_unpolarized.txt

Canonical output:
  output/stage5_ji/
"""

from __future__ import annotations
import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

LUMI = (1, 2, 5, 10)
KAPPA_P = 1.79284734463
KAPPA_N = -1.91304273
# Flavor Pauli zeroth moments, neglecting strange:
# kappa_p = 2/3 kappa_u - 1/3 kappa_d
# kappa_n = 2/3 kappa_d - 1/3 kappa_u
KAPPA_U = 2.0 * KAPPA_P + KAPPA_N
KAPPA_D = KAPPA_P + 2.0 * KAPPA_N

def beta_fn(a: float, b: float) -> float:
    return math.exp(math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))

def normalized_e_valence(x, kappa_q, alpha_q, beta_q):
    """e_v^q(x)=N x^{-alpha}(1-x)^beta with integral fixed to kappa_q."""
    x = np.asarray(x, dtype=float)
    norm = kappa_q / beta_fn(1.0-alpha_q, 1.0+beta_q)
    return norm * np.power(x, -alpha_q) * np.power(1.0-x, beta_q)

def b20_valence(kappa_q: float, alpha_q: float, beta_q: float) -> float:
    """Integral_0^1 dx x E_v^q(x,0,0), analytic for the chosen forward limit."""
    return kappa_q * beta_fn(2.0-alpha_q, 1.0+beta_q) / beta_fn(1.0-alpha_q, 1.0+beta_q)

def load_ndvcs_xs(path: Path) -> pd.DataFrame:
    names = ["kin_bin","phi_deg","Q2_GeV2","xB","t_abs_GeV2",
             "xs_pb_GeV4","stat_pb_GeV4","sys_pb_GeV4",
             "xs_nb_GeV4","stat_nb_GeV4","sys_nb_GeV4"]
    df = pd.read_csv(path, comment="#", sep=r"\s+", names=names, engine="python")
    for c in names:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna().reset_index(drop=True)
    if len(df) == 0:
        raise RuntimeError(f"No neutron cross-section rows parsed from {path}")
    return df

def make_rgb_projection_table(df):
    out = []
    for L in LUMI:
        d = df.copy()
        d["luminosity_factor"] = L
        d["stat_projected_nb_GeV4"] = d["stat_nb_GeV4"] / np.sqrt(L)
        d["sys_held_fixed_nb_GeV4"] = d["sys_nb_GeV4"]
        d["total_projected_nb_GeV4"] = np.hypot(
            d["stat_projected_nb_GeV4"], d["sys_held_fixed_nb_GeV4"])
        d["rel_stat_projected"] = d["stat_projected_nb_GeV4"] / d["xs_nb_GeV4"].abs()
        d["rel_total_projected"] = d["total_projected_nb_GeV4"] / d["xs_nb_GeV4"].abs()
        out.append(d)
    return pd.concat(out, ignore_index=True)

def summarize_rgb(proj):
    rows = []
    for L in LUMI:
        d = proj[proj.luminosity_factor == L]
        rows.append(dict(
            luminosity=f"{L}x",
            N=len(d),
            N_kin_bins=d.kin_bin.nunique(),
            median_rel_stat_pct=100*np.median(d.rel_stat_projected),
            median_rel_total_pct=100*np.median(d.rel_total_projected),
            p90_rel_stat_pct=100*np.quantile(d.rel_stat_projected, .90),
            p90_rel_total_pct=100*np.quantile(d.rel_total_projected, .90),
            stat_dominated_fraction=float(np.mean(
                d.stat_projected_nb_GeV4 > d.sys_held_fixed_nb_GeV4)),
        ))
    return pd.DataFrame(rows)

def plot_rgb_precision(summary, figdir):
    L = np.array([1,2,5,10], float)
    fig, ax = plt.subplots(figsize=(7.4,5.2))
    ax.plot(L, summary.median_rel_stat_pct, "o-", label="statistical only")
    ax.plot(L, summary.median_rel_total_pct, "s-", label="quoted systematics held fixed")
    ax.set_xscale("log")
    ax.set_xticks(L, [f"{int(v)}x" for v in L])
    ax.set_xlabel("RGB unpolarized luminosity")
    ax.set_ylabel("Median relative cross-section uncertainty (%)")
    ax.set_title("RGB neutron DVCS: where additional luminosity stops helping")
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figdir/"rgb_precision_vs_luminosity.png", dpi=220)
    plt.close(fig)

def plot_rgb_kinematic_map(df, figdir):
    # One point per preliminary kin_bin, using means over phi rows.
    g = df.groupby("kin_bin", as_index=False).agg(
        xB=("xB","mean"), Q2=("Q2_GeV2","mean"), tabs=("t_abs_GeV2","mean"),
        med_rel_stat=("stat_nb_GeV4", lambda x: np.nan))
    # compute relative stat correctly from original groups
    vals=[]
    for k, d in df.groupby("kin_bin"):
        vals.append((k, np.median(d.stat_nb_GeV4/d.xs_nb_GeV4.abs())))
    rel=dict(vals)
    g["med_rel_stat"]=g.kin_bin.map(rel)
    fig, ax = plt.subplots(figsize=(7.2,5.4))
    sc=ax.scatter(g.xB, g.Q2, s=80+350*g.tabs, c=100*g.med_rel_stat)
    for _,r in g.iterrows():
        ax.annotate(f"|t|={r.tabs:.2f}", (r.xB,r.Q2), xytext=(5,5),
                    textcoords="offset points", fontsize=8)
    cb=fig.colorbar(sc, ax=ax)
    cb.set_label("Median current statistical uncertainty (%)")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
    ax.set_title("Preliminary RGB neutron-DVCS kinematic leverage")
    ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(figdir/"rgb_kinematic_leverage.png", dpi=220)
    plt.close(fig)

def make_dfjk_seed_scan():
    """
    Workshop-level seed scan, not a fit.
    alpha_E is held common here only to visualize how beta_u,beta_d move B20.
    v2 will use the actual observable likelihood.
    """
    alpha = 0.55
    bu = np.linspace(2.0, 8.0, 121)
    bd = np.linspace(2.0, 10.0, 161)
    rows=[]
    for u in bu:
        for d in bd:
            rows.append((u,d,b20_valence(KAPPA_U,alpha,u),
                         b20_valence(KAPPA_D,alpha,d)))
    return pd.DataFrame(rows, columns=["beta_u","beta_d","B20_uv","B20_dv"])

def plot_dfjk_map(scan, figdir):
    fig, ax = plt.subplots(figsize=(7.0,5.4))
    sc=ax.scatter(scan.B20_uv, scan.B20_dv, c=scan.beta_u, s=5, alpha=.25)
    cb=fig.colorbar(sc, ax=ax)
    cb.set_label(r"$\beta_u^E$ (seed scan)")
    ax.set_xlabel(r"$B_{20}^{u_v}(0)=\int dx\,xE_v^u$")
    ax.set_ylabel(r"$B_{20}^{d_v}(0)=\int dx\,xE_v^d$")
    ax.set_title("DFJK-inspired E-sector parameter space (prior map, not data constraint)")
    ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(figdir/"dfjk_B20_parameter_map_prior_only.png", dpi=220)
    plt.close(fig)

def check_gepard():
    try:
        import gepard
        return True, getattr(gepard, "__file__", "available")
    except Exception as e:
        return False, str(e)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ndvcs-xs", default="import/ndvcs_clas12_preliminary_unpolarized.txt")
    ap.add_argument("--output", default="output/stage5_ji")
    args=ap.parse_args()

    out=Path(args.output)
    tab=out/"tables"; fig=out/"figures"
    tab.mkdir(parents=True, exist_ok=True); fig.mkdir(parents=True, exist_ok=True)

    xs_path=Path(args.ndvcs_xs)
    if not xs_path.exists():
        raise FileNotFoundError(
            f"{xs_path} not found. Put ndvcs_clas12_preliminary_unpolarized.txt "
            "in the import directory beside the high-luminosity scripts."
        )

    df=load_ndvcs_xs(xs_path)
    proj=make_rgb_projection_table(df)
    summ=summarize_rgb(proj)
    scan=make_dfjk_seed_scan()

    df.to_csv(tab/"rgb_ndvcs_xs_input_audit.csv", index=False)
    proj.to_csv(tab/"rgb_ndvcs_xs_luminosity_projections.csv", index=False)
    summ.to_csv(tab/"rgb_ndvcs_xs_projection_summary.csv", index=False)
    scan.to_csv(tab/"dfjk_B20_seed_scan_prior_only.csv", index=False)

    plot_rgb_precision(summ, fig)
    plot_rgb_kinematic_map(df, fig)
    plot_dfjk_map(scan, fig)

    gp_ok,gp_msg=check_gepard()

    print("="*100)
    print("STAGE 5 v1 — RGB + DFJK/Ji FRAMEWORK")
    print("="*100)
    print(f"RGB neutron XS: {len(df)} phi points in {df.kin_bin.nunique()} kinematic bins")
    print(f"xB range      : {df.xB.min():.3f} -- {df.xB.max():.3f}")
    print(f"Q2 range      : {df.Q2_GeV2.min():.3f} -- {df.Q2_GeV2.max():.3f} GeV^2")
    print(f"|t| range     : {df.t_abs_GeV2.min():.3f} -- {df.t_abs_GeV2.max():.3f} GeV^2")
    print()
    print(summ.to_string(index=False, float_format=lambda x:f"{x:.4g}"))
    print()
    print("DFJK-inspired E normalization from anomalous magnetic moments:")
    print(f"  kappa_u = 2*kappa_p + kappa_n = {KAPPA_U:+.6f}")
    print(f"  kappa_d = kappa_p + 2*kappa_n = {KAPPA_D:+.6f}")
    print("  B20 is computed analytically from the x-weighted E_v moment.")
    print()
    print(f"Gepard available: {gp_ok} ({gp_msg})")
    print()
    print("NEXT IMPLEMENTATION STEP (v2):")
    print("  1) add published RGB nDVCS BSA input;")
    print("  2) implement fixed-profile double-distribution skewness for DFJK H/E;")
    print("  3) predict p/n XS and BSA through Gepard;")
    print("  4) fit p(1x) -> p(1x)+n(1x) -> p(10x)+n(10x);")
    print("  5) propagate replicas to B20_uv,B20_dv and then Ju,Jd.")
    print()
    print(f"Wrote: {out}")

if __name__ == "__main__":
    raise SystemExit(main())
