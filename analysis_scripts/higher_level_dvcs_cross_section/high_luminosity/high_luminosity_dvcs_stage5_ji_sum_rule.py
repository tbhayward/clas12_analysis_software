#!/usr/bin/env python3
"""
Stage 5 — neutron DVCS / flavor separation / Ji-sum-rule preparation.

This is a standalone Stage-5 script.  The present revision does three things:

  1. parses the preliminary RGB neutron unpolarized cross-section table;
  2. treats the preliminary systematic uncertainty honestly as an *unresolved
     covariance problem*, rather than automatically converting the quoted total
     into independent point-to-point noise;
  3. keeps the DFJK-inspired E_v moment map that will become the parameter space
     for the proton+neutron fit.

The preliminary RGB note (v0.0, 22 Apr 2026) constructs the quoted systematic
uncertainty from analysis variations.  Several sources are intrinsically
coherent across phi bins (and, in some cases, across broader kinematics).
Therefore the table's single total systematic error per point does not determine
a covariance matrix.

For projections we consequently show three *brackets*, not three claims:
  - statistics only;
  - quoted total systematic treated as diagonal (maximally conservative for
    shape information);
  - quoted total systematic treated as fully correlated within each
    (Q2,xB,t) cell (optimistic shape-information bracket).

Neither systematic bracket is the production covariance model.  A production
model requires the individual signed shifts from the seven systematic studies.

Published RGB nDVCS BSA:
  If import/ndvcs_rgb_published_bsa.txt is present, it is detected and audited.
  The next fit revision will consume it.  We deliberately do not digitize the
  PRL figure in this script.

Output:
  output/stage5_ji/
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LUMI_FACTORS = (1, 2, 5, 10)

KAPPA_P = 1.79284734463
KAPPA_N = -1.91304273
KAPPA_U = 2.0 * KAPPA_P + KAPPA_N
KAPPA_D = KAPPA_P + 2.0 * KAPPA_N

# Simple beta-family seed used only to map the E_v moment freedom.
# e_v^q(x) = N_q x^{-alpha} (1-x)^beta, with N_q fixed by kappa_q.
DEFAULT_ALPHA_U = 0.55
DEFAULT_ALPHA_D = 0.55

RGB_XS_DEFAULT = Path("import/ndvcs_clas12_preliminary_unpolarized.txt")
RGB_BSA_DEFAULT = Path("import/ndvcs_rgb_published_bsa_digitized_t_projection.csv")
OUT_DEFAULT = Path("output/stage5_ji")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def beta_fn(a: float, b: float) -> float:
    return math.gamma(a) * math.gamma(b) / math.gamma(a + b)


def b20_from_beta(kappa: float, alpha: float, beta: float) -> float:
    """B20 = integral dx x E_v(x), with integral E_v = kappa."""
    den = beta_fn(1.0 - alpha, beta + 1.0)
    num = beta_fn(2.0 - alpha, beta + 1.0)
    return kappa * num / den


def _float_tokens(line: str):
    vals = []
    for tok in re.split(r"[\s,]+", line.strip()):
        try:
            vals.append(float(tok))
        except ValueError:
            pass
    return vals


def load_rgb_xs(path: Path) -> pd.DataFrame:
    """
    Parse the Stage-5 RGB text table.

    Preferred format:
      kin_bin phi_deg Q2_GeV2 xB t_abs_GeV2 xs_pb_GeV4 stat_pb_GeV4 sys_pb_GeV4

    Comment/header lines are ignored.  A whitespace table with >=8 numeric
    columns is also accepted.
    """
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")

    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            vals = _float_tokens(s)
            if len(vals) < 8:
                continue
            rows.append(vals[:8])

    if not rows:
        raise RuntimeError(f"No 8-column RGB cross-section rows found in {path}")

    df = pd.DataFrame(
        rows,
        columns=[
            "kin_bin", "phi_deg", "Q2_GeV2", "xB", "t_abs_GeV2",
            "xs_pb_GeV4", "stat_pb_GeV4", "sys_pb_GeV4",
        ],
    )
    df["kin_bin"] = df["kin_bin"].astype(int)
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna().reset_index(drop=True)
    df["rel_stat"] = np.abs(df["stat_pb_GeV4"] / df["xs_pb_GeV4"])
    df["rel_sys"] = np.abs(df["sys_pb_GeV4"] / df["xs_pb_GeV4"])
    return df


def load_optional_bsa(path: Path):
    """
    Audit an optional user-supplied published RGB BSA table.

    Accepted rows must contain at least:
      phi  Q2  xB  |t|  ALU  stat
    with an optional seventh column for systematic uncertainty.

    This loader is intentionally permissive because the exact CLAS-database
    export format may differ.  It does not silently fabricate missing errors.
    """
    if not path.exists():
        return None

    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
        required = {"phi_deg", "Q2_GeV2", "xB", "t_abs_GeV2", "stat_digitized"}
        missing = required - set(df.columns)
        if missing:
            raise RuntimeError(f"{path}: missing required BSA columns {sorted(missing)}")
        return df

    rows = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            vals = _float_tokens(s)
            if len(vals) >= 6:
                rows.append(vals[:7])

    if not rows:
        raise RuntimeError(f"{path} exists but no >=6-column numeric BSA rows were found")

    ncol = max(len(r) for r in rows)
    rows = [r + [np.nan] * (ncol - len(r)) for r in rows]
    cols = ["phi_deg", "Q2_GeV2", "xB", "t_abs_GeV2", "ALU", "stat"]
    if ncol >= 7:
        cols.append("sys")
    df = pd.DataFrame(rows, columns=cols)
    return df


# ---------------------------------------------------------------------------
# RGB systematic-covariance audit
# ---------------------------------------------------------------------------

def make_rgb_projection_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for L in LUMI_FACTORS:
        stat = df["stat_pb_GeV4"].to_numpy() / math.sqrt(L)
        sys = df["sys_pb_GeV4"].to_numpy()
        xs = np.abs(df["xs_pb_GeV4"].to_numpy())

        rel_stat = np.abs(stat / xs)
        rel_diag = np.sqrt(stat**2 + sys**2) / xs

        rows.append({
            "luminosity": f"{L}x",
            "N": len(df),
            "N_kin_bins": df["kin_bin"].nunique(),
            "median_rel_stat_pct": 100*np.median(rel_stat),
            "median_rel_total_diag_pct": 100*np.median(rel_diag),
            "p90_rel_stat_pct": 100*np.quantile(rel_stat, 0.90),
            "p90_rel_total_diag_pct": 100*np.quantile(rel_diag, 0.90),
            "stat_dominated_fraction_diag": float(np.mean(stat > sys)),
        })
    return pd.DataFrame(rows)


def cell_shape_precision(df: pd.DataFrame, L: float, mode: str) -> np.ndarray:
    """
    Return a simple per-cell fractional *shape* precision diagnostic.

    We remove one arbitrary normalization direction in each phi distribution.
    This is not a physics fit; it only illustrates how covariance assumptions
    change the amount of relative-phi information.

    mode:
      stat          -> diagonal statistical covariance
      sys_diagonal  -> stat + quoted total systematic diagonal
      sys_cellcorr  -> stat diagonal + quoted total systematic represented as
                       a single coherent cell-scale direction.

    For the correlated bracket, the point-dependent quoted sys magnitudes are
    used to construct one rank-1 shift vector s_i s_j.  This is deliberately an
    optimistic bracket, not the production covariance.
    """
    out = []
    for _, g in df.groupby("kin_bin"):
        y = g["xs_pb_GeV4"].to_numpy(float)
        st = g["stat_pb_GeV4"].to_numpy(float) / math.sqrt(L)
        sy = g["sys_pb_GeV4"].to_numpy(float)
        n = len(g)
        if n < 3:
            continue

        if mode == "stat":
            C = np.diag(st**2)
        elif mode == "sys_diagonal":
            C = np.diag(st**2 + sy**2)
        elif mode == "sys_cellcorr":
            C = np.diag(st**2) + np.outer(sy, sy)
        else:
            raise ValueError(mode)

        # Remove one normalization-like direction.  Work in fractional units.
        scale = np.maximum(np.abs(y), 1e-30)
        Cf = C / np.outer(scale, scale)
        one = np.ones(n)
        try:
            W = np.linalg.pinv(Cf, rcond=1e-12)
            denom = one @ W @ one
            P = W - np.outer(W @ one, W @ one) / denom if denom > 0 else W
            # Effective average uncertainty of the n-1 relative-shape directions.
            evals = np.linalg.eigvalsh((P + P.T)/2)
            pos = evals[evals > max(evals.max(initial=0)*1e-10, 1e-14)]
            if len(pos):
                out.append(math.sqrt(np.mean(1.0/pos)))
        except np.linalg.LinAlgError:
            pass
    return np.asarray(out)


def plot_rgb_precision_brackets(df: pd.DataFrame, outdir: Path):
    xs = np.arange(len(LUMI_FACTORS))
    med_stat, med_diag, med_corr = [], [], []

    for L in LUMI_FACTORS:
        a = cell_shape_precision(df, L, "stat")
        b = cell_shape_precision(df, L, "sys_diagonal")
        c = cell_shape_precision(df, L, "sys_cellcorr")
        med_stat.append(100*np.median(a))
        med_diag.append(100*np.median(b))
        med_corr.append(100*np.median(c))

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.plot(xs, med_stat, "o-", label="statistics only")
    ax.plot(xs, med_diag, "s-", label="quoted total sys treated point-to-point")
    ax.plot(xs, med_corr, "^-", label="quoted total sys correlated within cell")
    ax.set_xticks(xs, [f"{L}x" for L in LUMI_FACTORS])
    ax.set_xlabel("RGB unpolarized luminosity")
    ax.set_ylabel("Median relative-phi shape uncertainty (%)")
    ax.set_title("RGB neutron DVCS: covariance assumption controls the luminosity projection")
    ax.grid(alpha=0.25)
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(outdir / "rgb_shape_precision_covariance_brackets.png", dpi=180)
    plt.close(fig)

    return pd.DataFrame({
        "luminosity": [f"{L}x" for L in LUMI_FACTORS],
        "shape_stat_only_pct": med_stat,
        "shape_sys_diagonal_pct": med_diag,
        "shape_sys_cellcorr_pct": med_corr,
    })


def plot_rgb_kinematics(df: pd.DataFrame, outdir: Path):
    # One marker per kinematic cell; labels deliberately omitted to avoid the
    # clutter in v1.  Color = current median statistical uncertainty in that cell.
    cells = []
    for k, g in df.groupby("kin_bin"):
        cells.append({
            "kin_bin": k,
            "xB": g["xB"].mean(),
            "Q2": g["Q2_GeV2"].mean(),
            "t": g["t_abs_GeV2"].mean(),
            "stat": 100*np.median(g["rel_stat"]),
        })
    c = pd.DataFrame(cells)

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    sc = ax.scatter(c["xB"], c["Q2"], c=c["stat"], s=55 + 110*c["t"])
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Median current statistical uncertainty (%)")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
    ax.set_title("Preliminary RGB neutron-DVCS kinematic leverage")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "rgb_kinematic_leverage.png", dpi=180)
    plt.close(fig)


# ---------------------------------------------------------------------------
# DFJK-inspired moment map
# ---------------------------------------------------------------------------

def make_b20_prior_map(figdir: Path, tabdir: Path):
    beta_u = np.linspace(2.0, 8.0, 121)
    beta_d = np.linspace(2.0, 8.0, 121)
    rows = []
    for bu in beta_u:
        B_u = b20_from_beta(KAPPA_U, DEFAULT_ALPHA_U, bu)
        for bd in beta_d:
            B_d = b20_from_beta(KAPPA_D, DEFAULT_ALPHA_D, bd)
            rows.append((bu, bd, B_u, B_d))
    df = pd.DataFrame(rows, columns=["beta_u_E", "beta_d_E", "B20_uv", "B20_dv"])
    df.to_csv(tabdir / "dfjk_B20_prior_map.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    sc = ax.scatter(df["B20_uv"], df["B20_dv"], c=df["beta_u_E"], s=7, alpha=0.65)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$\beta_u^E$ (seed scan)")
    ax.set_xlabel(r"$B_{20}^{u_v}(0)=\int dx\,xE_v^u$")
    ax.set_ylabel(r"$B_{20}^{d_v}(0)=\int dx\,xE_v^d$")
    ax.set_title("DFJK-inspired E-sector parameter space (prior map, not data constraint)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(figdir / "dfjk_B20_parameter_map_prior_only.png", dpi=180)
    plt.close(fig)
    return df



def bsa_4d_feasibility(bsa: pd.DataFrame, figdir: Path, tabdir: Path):
    """
    Workshop-level feasibility diagnostic for replacing the published 1D neutron
    BSA projections by multidimensional binning.

    If a published bin is subdivided into S equally populated multidimensional
    cells, the statistical uncertainty scales approximately as sqrt(S/L).
    We therefore propagate the *measured published statistical error bars* as

        sigma_4D = sigma_published * sqrt(S / L).

    S is a split factor, not a proposed final binning.  This deliberately avoids
    claiming that 10x luminosity guarantees a particular 4D grid.
    """
    if "stat_digitized" not in bsa.columns:
        return None

    split_factors = [1, 2, 4, 6, 8, 10]
    rows = []
    stat0 = bsa["stat_digitized"].to_numpy(float)
    for L in LUMI_FACTORS:
        for S in split_factors:
            e = stat0 * np.sqrt(float(S) / float(L))
            rows.append({
                "luminosity": f"{L}x",
                "subcells_per_published_bin": S,
                "median_sigma_ALU": np.median(e),
                "p90_sigma_ALU": np.quantile(e, 0.90),
                "fraction_sigma_below_0p03": np.mean(e < 0.03),
                "fraction_sigma_below_0p05": np.mean(e < 0.05),
                "fraction_sigma_below_0p10": np.mean(e < 0.10),
            })
    out = pd.DataFrame(rows)
    out.to_csv(tabdir / "rgb_bsa_multidimensional_feasibility.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    for L in LUMI_FACTORS:
        q = out[out["luminosity"] == f"{L}x"]
        ax.plot(q["subcells_per_published_bin"], q["median_sigma_ALU"],
                marker="o", label=f"{L}x")
    ax.set_xlabel("Equal-statistics subcells per published 1D bin")
    ax.set_ylabel(r"Median projected statistical $\sigma(A_{LU})$")
    ax.set_title("RGB neutron BSA: statistical cost of multidimensional binning")
    ax.grid(alpha=0.25)
    ax.legend(title="Luminosity")
    fig.tight_layout()
    fig.savefig(figdir / "rgb_bsa_multidimensional_feasibility.png", dpi=180)
    plt.close(fig)
    return out

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rgb-xs", type=Path, default=RGB_XS_DEFAULT)
    ap.add_argument("--rgb-bsa", type=Path, default=RGB_BSA_DEFAULT)
    ap.add_argument("--output", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args()

    figdir = args.output / "figures"
    tabdir = args.output / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    # Remove obsolete v2 root-level products so output/stage5_ji has only
    # figures/ and tables/ after a successful run.
    for stale in args.output.glob("*"):
        if stale.is_file():
            stale.unlink()

    xs = load_rgb_xs(args.rgb_xs)
    proj = make_rgb_projection_table(xs)
    proj.to_csv(tabdir / "rgb_xs_luminosity_projection.csv", index=False)

    brackets = plot_rgb_precision_brackets(xs, figdir)
    brackets.to_csv(tabdir / "rgb_xs_covariance_brackets.csv", index=False)

    plot_rgb_kinematics(xs, figdir)
    make_b20_prior_map(figdir, tabdir)

    bsa = load_optional_bsa(args.rgb_bsa)
    if bsa is not None:
        bsa.to_csv(tabdir / "rgb_published_bsa_audit.csv", index=False)
        feasibility = bsa_4d_feasibility(bsa, figdir, tabdir)
    else:
        feasibility = None

    print("=" * 100)
    print("STAGE 5 v3 — RGB COVARIANCE + MULTIDIMENSIONAL-BSA FEASIBILITY + DFJK/Ji FRAMEWORK")
    print("=" * 100)
    print(f"RGB neutron XS: {len(xs)} phi points in {xs['kin_bin'].nunique()} kinematic bins")
    print(f"xB range      : {xs.xB.min():.3f} -- {xs.xB.max():.3f}")
    print(f"Q2 range      : {xs.Q2_GeV2.min():.3f} -- {xs.Q2_GeV2.max():.3f} GeV^2")
    print(f"|t| range     : {xs.t_abs_GeV2.min():.3f} -- {xs.t_abs_GeV2.max():.3f} GeV^2")
    print()
    print(proj.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print()
    print("IMPORTANT SYSTEMATICS INTERPRETATION:")
    print("  The preliminary note's total systematic error is NOT a measured point-to-point covariance.")
    print("  It is the quadrature sum of analysis-variation sources.")
    print("  v2 therefore keeps diagonal-total and cell-correlated treatments only as projection brackets.")
    print("  The production p+n GPD fit should use source-by-source signed shifts/covariances if available.")
    print()
    print("RGB covariance-bracket shape diagnostic:")
    print(brackets.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print()
    print("DFJK-inspired E normalization from anomalous magnetic moments:")
    print(f"  kappa_u = 2*kappa_p + kappa_n = {KAPPA_U:+.6f}")
    print(f"  kappa_d = kappa_p + 2*kappa_n = {KAPPA_D:+.6f}")
    print()
    if bsa is None:
        print(f"Published RGB BSA: not yet present ({args.rgb_bsa})")
        print("  Put the CLAS-database export there; v2 will audit it automatically.")
    else:
        print(f"Published RGB BSA template: loaded {len(bsa)} rows from {args.rgb_bsa}")
        print("  Columns:", ", ".join(bsa.columns))
        if feasibility is not None:
            print()
            print("Multidimensional-BSA feasibility (selected rows):")
            sel = feasibility[
                feasibility["subcells_per_published_bin"].isin([1, 4, 8])
            ]
            print(sel.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print()
    print("NEXT FIT STEP:")
    print("  Use the published neutron BSA + proton XS/BSA + preliminary neutron XS in a")
    print("  finite-skewness H/E model, then compare p(1x), p(1x)+n(1x), p(10x),")
    print("  and p(10x)+n(10x) directly in the B20_uv-B20_dv plane.")
    print()
    print()
    print("OUTPUT LAYOUT:")
    print(f"  figures -> {figdir}")
    print(f"  tables  -> {tabdir}")
    print(f"Wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
