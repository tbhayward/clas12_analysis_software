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
import re
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
RGB_BSA_DIR_DEFAULT = None  # resolved relative to this script in main()
RGA_PASS2_DEFAULT = None  # resolved relative to this script below
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


def load_actual_rgb_bsa_directory(path: Path) -> pd.DataFrame:
    """Load the actual published nDVCS per-phi BSA points from Adam/Silvia.

    The three per-phi files are alternative 1D projections of the same event
    sample (xB, Q2, and t). They are retained with a `projection` label and
    MUST NOT be combined as statistically independent measurements.
    """
    files = {
        "xB": path / "ndvcs_xbbins_centralkinematics_BSA_per_phi_bin.dat",
        "Q2": path / "ndvcs_q2bins_centralkinematics_BSA_per_phi_bin.dat",
        "t":  path / "ndvcs_tbins_centralkinematics_BSA_per_phi_bin.dat",
    }
    missing = [str(p) for p in files.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing actual nDVCS BSA per-phi files:\n  " + "\n  ".join(missing)
        )

    pat = re.compile(
        r"<Q2>\s*([-+0-9.eE]+)\s+"
        r"<xb>\s*([-+0-9.eE]+)\s+"
        r"<t>\s*([-+0-9.eE]+)\s+"
        r"<phi>\s*([-+0-9.eE]+)\s+"
        r"BSA\s*([-+0-9.eE]+)\s+"
        r"BSA_error_stat\s*([-+0-9.eE]+)\s+"
        r"BSA_error_syst\s*([-+0-9.eE]+)"
    )
    frames = []
    for projection, fpath in files.items():
        rows = []
        for line in fpath.read_text().splitlines():
            m = pat.search(line)
            if not m:
                continue
            Q2, xb, t, phi, alu, estat, esys = map(float, m.groups())
            rows.append({
                "projection": projection,
                "Q2_GeV2": Q2,
                "xB": xb,
                "t_abs_GeV2": abs(t),
                "phi_deg": phi,
                "ALU": alu,
                "stat": estat,
                "sys": esys,
                "source_file": fpath.name,
            })
        if not rows:
            raise RuntimeError(f"No nDVCS BSA rows parsed from {fpath}")
        df = pd.DataFrame(rows)
        # Adam/Silvia files contain 3 projected bins x 7 phi bins = 21 rows.
        df["projected_bin"] = np.repeat(
            np.arange(1, 1 + int(np.ceil(len(df) / 7))), 7
        )[:len(df)]
        frames.append(df)

    return pd.concat(frames, ignore_index=True)


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



def _tuple_component_stage5(series: pd.Series, index: int) -> np.ndarray:
    """Extract a numeric component from tuple-valued pass-2 CSV cells."""
    import ast
    out = np.full(len(series), np.nan, dtype=float)
    for i, value in enumerate(series):
        if pd.isna(value):
            continue
        try:
            parsed = ast.literal_eval(str(value).strip())
            if isinstance(parsed, (tuple, list)) and len(parsed) > index:
                out[i] = float(parsed[index])
        except Exception:
            continue
    return out


def load_rga_pass2_bsa(path: Path):
    """Load current pass-2 RGA 4D BSA and its statistical uncertainty."""
    if not path.exists():
        return None

    raw = pd.read_csv(path, low_memory=False)
    bsa_col = "BSA, counts, 10.6 GeV"
    if bsa_col not in raw.columns:
        raise KeyError(f"{path}: missing required column {bsa_col!r}")

    bsa = _tuple_component_stage5(raw[bsa_col], 0)
    stat = _tuple_component_stage5(raw[bsa_col], 1)

    out = pd.DataFrame(index=raw.index)
    for old, new in [
        ("Bin Name", "bin"),
        ("xBavg, 10.6 GeV", "xB"),
        ("Q2avg, 10.6 GeV", "Q2"),
        ("t_abs_avg, 10.6 GeV", "t_abs"),
        ("phiavg, 10.6 GeV", "phi_deg"),
    ]:
        if old in raw.columns:
            out[new] = pd.to_numeric(raw[old], errors="coerce")
    out["BSA_pass2"] = bsa
    out["stat_error"] = stat
    good = np.isfinite(out["BSA_pass2"]) & np.isfinite(out["stat_error"]) & (out["stat_error"] > 0)
    out = out.loc[good].reset_index(drop=True)
    out.attrs["bsa_column"] = bsa_col
    out.attrs["stat_column"] = f"{bsa_col} tuple component 1"
    return out


def rga_binning_transfer_projection(rgb_bsa: pd.DataFrame, rga_bsa: pd.DataFrame,
                                    figdir: Path, tabdir: Path,
                                    reference_projection: str = "t"):
    """Project RGB BSA precision into the *actual RGA pass-2 4D bin pattern*.

    Workshop assumption requested by the user:
      neutron and proton accepted phase-space populations have the same shape.

    We use the RGA point-by-point statistical-error pattern as the empirical
    relative occupancy template.  Since statistical information scales
    approximately as 1/sigma^2, define RGA weights w_i proportional to
    1/sigma_RGA,i^2.

    The absolute RGB information budget is obtained from ONE published neutron
    projection only (default: t), because xB/Q2/t projection files are three
    views of the same events and cannot be added as independent statistics.

    For luminosity factor L:
        I_RGB(L) = L * sum_j 1/sigma_RGB,j^2
        I_i      = I_RGB(L) * w_i
        sigma_i  = 1/sqrt(I_i)

    This deliberately assumes identical p/n phase-space shape and ignores
    possible differences in analyzing power, beam polarization, efficiency,
    backgrounds, and migration.  It is a workshop-level binning/statistics
    projection, not a detector simulation.
    """
    nref = rgb_bsa[rgb_bsa["projection"] == reference_projection].copy()
    if nref.empty:
        raise RuntimeError(f"No RGB BSA rows for reference projection {reference_projection!r}")

    nstat = nref["stat"].to_numpy(float)
    pstat = rga_bsa["stat_error"].to_numpy(float)
    nstat = nstat[np.isfinite(nstat) & (nstat > 0)]
    pstat = pstat[np.isfinite(pstat) & (pstat > 0)]

    I_rgb_1x = float(np.sum(1.0 / nstat**2))
    raw_w = 1.0 / pstat**2
    w = raw_w / np.sum(raw_w)

    summary_rows = []
    point_tables = []
    for L in [1, 2, 5, 10]:
        sigma = 1.0 / np.sqrt((L * I_rgb_1x) * w)
        pt = rga_bsa.copy()
        pt["projection_luminosity"] = L
        pt["rgb_projected_stat_error"] = sigma
        pt["ratio_rgb_projected_to_rga_stat"] = sigma / pt["stat_error"].to_numpy(float)
        point_tables.append(pt)

        summary_rows.append({
            "luminosity": f"{L}x",
            "N_RGA_like_4D_points": len(sigma),
            "median_projected_sigma_ALU": float(np.median(sigma)),
            "p16_projected_sigma_ALU": float(np.quantile(sigma, 0.16)),
            "p84_projected_sigma_ALU": float(np.quantile(sigma, 0.84)),
            "median_current_RGA_sigma_ALU": float(np.median(pstat)),
            "median_ratio_projected_RGB_to_current_RGA": float(np.median(sigma) / np.median(pstat)),
            "fraction_projected_better_than_current_RGA_point": float(np.mean(sigma < pstat)),
            "RGB_reference_projection": reference_projection,
            "RGB_reference_points": len(nstat),
            "RGB_reference_median_sigma_ALU": float(np.median(nstat)),
        })

    summary = pd.DataFrame(summary_rows)
    points = pd.concat(point_tables, ignore_index=True)
    summary.to_csv(tabdir / "rgb_on_rga_4d_binning_summary.csv", index=False)
    points.to_csv(tabdir / "rgb_on_rga_4d_binning_points.csv", index=False)

    # Clean main-talk plot: distributions for current RGA and projected RGB at 10x.
    rgb10 = points[points["projection_luminosity"] == 10]["rgb_projected_stat_error"].to_numpy(float)
    bins = np.linspace(
        0.0,
        max(np.quantile(pstat, 0.95), np.quantile(rgb10, 0.95)) * 1.08,
        34,
    )
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.hist(pstat, bins=bins, histtype="step", linewidth=2.2,
            label=f"Current RGA proton 4D (N={len(pstat)})")
    ax.hist(rgb10, bins=bins, histtype="step", linewidth=2.2,
            label=f"RGB neutron projected 10x on same 4D pattern (N={len(rgb10)})")
    ax.axvline(np.median(pstat), linestyle=":", linewidth=2.0,
               label=f"RGA median = {np.median(pstat):.3f}")
    ax.axvline(np.median(rgb10), linestyle="--", linewidth=2.0,
               label=f"RGB 10x median = {np.median(rgb10):.3f}")
    ax.set_xlabel(r"Statistical $\sigma(A_{LU})$")
    ax.set_ylabel("Number of 4D points")
    ax.set_title("Neutron BSA projected onto the actual RGA pass-2 4D bin pattern")
    ax.grid(alpha=0.20)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(figdir / "rgb_10x_on_rga_4d_binning_precision.png", dpi=180)
    plt.close(fig)

    # Luminosity progression in the same fixed RGA 4D binning.
    #
    # Once the bin pattern is fixed, this is pure counting-statistics scaling:
    # sigma(L) = sigma(1x) / sqrt(L).  Plot it continuously to 30x and solve
    # explicitly for the luminosity at which the RGB and current-RGA medians
    # are equal.
    median_rgb_1x = float(summary.loc[summary["luminosity"] == "1x",
                                      "median_projected_sigma_ALU"].iloc[0])
    median_rga = float(np.median(pstat))
    L_cross = (median_rgb_1x / median_rga) ** 2

    L_curve = np.linspace(1.0, 30.0, 400)
    sigma_curve = median_rgb_1x / np.sqrt(L_curve)

    fig, ax = plt.subplots(figsize=(7.7, 5.2))
    ax.plot(L_curve, sigma_curve, linewidth=2.2,
            label="Projected RGB neutron")
    ax.scatter([1, 2, 5, 10],
               summary["median_projected_sigma_ALU"],
               s=42, zorder=3, label="Calculated luminosity points")
    ax.axhline(median_rga, linestyle=":", linewidth=2.0,
               label=f"Current RGA proton median = {median_rga:.3f}")

    if 1.0 <= L_cross <= 30.0:
        ax.scatter([L_cross], [median_rga], s=70, zorder=5)
        ax.axvline(L_cross, linestyle="--", linewidth=1.4, alpha=0.7)
        ax.annotate(
            f"Equal median precision\nL = {L_cross:.1f}x",
            xy=(L_cross, median_rga),
            xytext=(L_cross - 7.5, median_rga + 0.055),
            arrowprops=dict(arrowstyle="->", lw=1.2),
            ha="center",
            va="bottom",
        )

    ax.set_xlim(1, 30)
    ax.set_xticks([1, 2, 5, 10, 15, 20, 25, 30])
    ax.set_xlabel("RGB luminosity factor")
    ax.set_ylabel(r"Median statistical $\sigma(A_{LU})$")
    ax.set_title("Precision at fixed RGA-like 4D granularity")
    ax.grid(alpha=0.20)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(figdir / "rgb_rga_4d_binning_luminosity_progression.png", dpi=180)
    plt.close(fig)

    # Save the crossing as an explicit numerical result for audit/use in slides.
    crossing = pd.DataFrame([{
        "median_RGB_1x_on_RGA_binning": median_rgb_1x,
        "median_current_RGA": median_rga,
        "equal_median_precision_luminosity_factor": L_cross,
    }])
    crossing.to_csv(tabdir / "rgb_rga_4d_binning_precision_crossing.csv", index=False)

    return summary, points

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rgb-xs", type=Path, default=RGB_XS_DEFAULT)
    ap.add_argument(
        "--rgb-bsa-dir",
        type=Path,
        default=None,
        help="Actual published nDVCS BSA directory; default is ../import/nDVCS_BSA relative to this script.",
    )
    ap.add_argument(
        "--rga-pass2",
        type=Path,
        default=None,
        help="Pass-2 CSV; default: ../higher_level_dvcs_cross_section/import/dvcs_pass2_analysis.csv relative to this script.",
    )
    ap.add_argument("--output", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args()
    here = Path(__file__).resolve().parent
    if args.rga_pass2 is None:
        args.rga_pass2 = (
            here.parent
            / "import"
            / "dvcs_pass2_analysis.csv"
        ).resolve()
    else:
        args.rga_pass2 = args.rga_pass2.expanduser().resolve()

    print(f"[Stage5] RGA pass-2 CSV path: {args.rga_pass2}")
    print(f"[Stage5] RGA pass-2 CSV exists: {args.rga_pass2.exists()}")

    if args.rgb_bsa_dir is None:
        args.rgb_bsa_dir = (here.parent / "import" / "nDVCS_BSA").resolve()
    else:
        args.rgb_bsa_dir = args.rgb_bsa_dir.expanduser().resolve()

    print(f"[Stage5] RGB BSA directory: {args.rgb_bsa_dir}")
    print(f"[Stage5] RGB BSA directory exists: {args.rgb_bsa_dir.exists()}")

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

    bsa = load_actual_rgb_bsa_directory(args.rgb_bsa_dir)
    bsa.to_csv(tabdir / "rgb_published_bsa_actual_points.csv", index=False)

    rga_bsa = load_rga_pass2_bsa(args.rga_pass2)
    if rga_bsa is not None:
        rga_bsa.to_csv(tabdir / "rga_pass2_bsa_audit.csv", index=False)
        transfer_summary, transfer_points = rga_binning_transfer_projection(
            bsa, rga_bsa, figdir, tabdir, reference_projection="t"
        )
    else:
        transfer_summary = None
        transfer_points = None

    print("=" * 100)
    print("STAGE 5 v12 — ACTUAL RGB BSA + RGA-LIKE 4D BINNING TRANSFER + DFJK/Ji FRAMEWORK")
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
    print(f"Actual published RGB BSA: loaded {len(bsa)} rows from {args.rgb_bsa_dir}")
    for projection in ["xB", "Q2", "t"]:
        q = bsa[bsa["projection"] == projection]
        p16, p84 = np.quantile(q["stat"], [0.16, 0.84])
        print(
            f"  {projection:>2s} projection: N={len(q):2d}, "
            f"median sigma(A_LU)={np.median(q['stat']):.5f}, "
            f"central 68%={p16:.5f}--{p84:.5f}"
        )
    print("  NOTE: xB, Q2, and t files are alternative projections of the same event sample;")
    print("        they are NOT summed as independent statistics.")
    print()

    if transfer_summary is None:
        print(f"RGA proton pass-2 BSA benchmark: not loaded ({args.rga_pass2})")
    else:
        q10 = transfer_summary[transfer_summary["luminosity"] == "10x"].iloc[0]
        print("RGA-like 4D neutron-BSA projection:")
        print("  Workshop assumption: neutron and proton accepted phase-space populations")
        print("                       have the same shape.")
        print(f"  RGA template: {len(rga_bsa)} actual pass-2 4D phi points")
        print("  RGB normalization: actual published t-projection only (21 disjoint t/phi points)")
        print(f"  Current RGA median sigma(A_LU)       = {q10.median_current_RGA_sigma_ALU:.5f}")
        print(f"  Projected RGB 10x median on RGA bins = {q10.median_projected_sigma_ALU:.5f}")
        print(f"  Ratio RGB(10x)/RGA(current) median   = {q10.median_ratio_projected_RGB_to_current_RGA:.3f}")
        median_rgb_1x = float(
            transfer_summary.loc[transfer_summary["luminosity"] == "1x",
                                 "median_projected_sigma_ALU"].iloc[0]
        )
        L_cross = (median_rgb_1x / float(q10.median_current_RGA_sigma_ALU)) ** 2
        print(f"  Equal-median-precision luminosity    = {L_cross:.2f}x")
        print(f"  Fraction of RGA-like bins where projected RGB 10x has smaller stat error")
        print(f"                                        = {q10.fraction_projected_better_than_current_RGA_point:.3f}")
        print()
        print("Luminosity progression at fixed RGA-like 4D granularity:")
        print(transfer_summary.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
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
