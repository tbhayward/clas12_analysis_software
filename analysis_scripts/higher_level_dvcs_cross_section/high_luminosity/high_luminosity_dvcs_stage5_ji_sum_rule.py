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

World-comparison sources used by the Ji comparison figure
---------------------------------------------------------
The numerical literature inputs below are intentionally documented here so a
future rerun can recover exactly what was plotted and why.

HERMES proton DVCS transverse-target constraint (full-flavor Ju,Jd, NOT
valence-only):
  HERMES Collaboration / Airapetian et al., 2008 result as summarized in
  Aidala et al., Rev. Mod. Phys. 85, 655 (2013):
      Ju + Jd/2.8 = 0.49 +/- 0.17 (experimental)
  https://cds.cern.ch/record/1478050/files/RevModPhys.85.655.pdf
  The HERMES 2002--04 preliminary form was Ju + Jd/2.9 = 0.42 +/- 0.21
  (exp_tot) +/- 0.06 model; see the HERMES/DESY transparency archive.
  IMPORTANT: the HERMES archive explicitly warns that its old 2008 Ju-vs-Jd
  conference plot should no longer be shown because the Dual-model curve had a
  theory error and the remaining VGG variants did not describe all HERMES data.
  We therefore reconstruct ONLY the published one-dimensional VGG/DD constraint
  above; we do not reproduce the deprecated HERMES plot.
  Archive: https://www.desy.de/~w3hermes/trans-public-author.html

Diehl--Kroll elastic-form-factor GPD extraction (VALENCE):
  M. Diehl and P. Kroll, Eur. Phys. J. C 73, 2397 (2013), arXiv:1302.4604
  https://arxiv.org/abs/1302.4604
      J_uv = 0.230 +0.009 -0.024
      J_dv = -0.004 +0.010 -0.016
  quoted at mu = 2 GeV.

Cichy--Constantinou--Sznajder--Wagner elastic+lattice extraction (VALENCE):
  K. Cichy et al., Phys. Rev. D 110, 114025 (2024)
  https://doi.org/10.1103/PhysRevD.110.114025
      J_uv = 0.195 +/- 0.010
      J_dv = 0.0173 +/- 0.0046
  Values are also explicitly tabulated in P. Sznajder's BNL 2024 GPD slides:
  https://indico.bnl.gov/event/24891/contributions/96844/attachments/58671/100775/BNL24_SB.pdf

Hall-A neutron DVCS historical constraint (full-flavor Ju,Jd, VGG model):
  Mazouz et al., Phys. Rev. Lett. 99, 242501 (2007), as summarized together
  with HERMES in Aidala et al., Rev. Mod. Phys. 85, 655 (2013):
      Jd + Ju/5.0 = 0.18 +/- 0.14 (experimental)
  https://cds.cern.ch/record/1478050/files/RevModPhys.85.655.pdf
  This is a one-dimensional model-dependent combination, not separate Ju/Jd
  marginal uncertainties.

Future EIC/EicC DVCS projection (full-flavor Ju,Jd, VGG-style model):
  X. Ji et al., ``Deeply virtual compton scattering at future electron-ion
  colliders'', Eur. Phys. J. C 83 (2023), DOI 10.1140/epjc/s10052-023-12065-x
  https://link.springer.com/article/10.1140/epjc/s10052-023-12065-x
  In the HERMES-like kinematic region, including statistical + assumed
  experimental systematic errors (their Eqs. 49--50):
      EicC: Ju + Jd/2.9 = 0.41 +/- 0.08
      EIC : Ju + Jd/3.0 = 0.39 +/- 0.06
  Statistics-only values are 0.06 (EicC) and 0.04 (EIC), but the production
  comparison uses the stat+syst values to avoid giving the collider projection
  an artificially favorable treatment.  Their small-x projections reach
  roughly 0.04--0.05 on related Ju+c*Jd combinations; these are documented but
  not mixed into the HERMES-region benchmark because the kinematics differ.

Hashamipour et al. valence extraction -- documented, not numerically plotted:
  H. Hashamipour et al., Phys. Rev. D 105, 054002 (2022),
  DOI 10.1103/PhysRevD.105.054002.
  https://journals.aps.org/prd/abstract/10.1103/PhysRevD.105.054002
  This paper extracts valence H_v and E_v from proton/neutron elastic FF and
  radius data and presents J_uv,J_dv with uncertainties at mu=2 GeV in Fig. 11.
  The article does not tabulate compact numerical Ji uncertainties.  We do NOT
  digitize the figure or invent numbers here; add it to the numerical comparison
  only if an author table/data file or a documented digitization is adopted.

Comparability warning:
  The CLAS12 projection, Diehl--Kroll, and Cichy et al. entries are valence
  quantities. HERMES, Hall-A, EIC and EicC entries are historical/projected
  full-flavor one-combination constraints. The uncertainty-only figure therefore
  separates these into two groups rather than treating all markers as the same
  observable.
"""

from __future__ import annotations

import argparse
import re
import math
import contextlib
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

# ---------------------------------------------------------------------------
# Stage-5 Step 1: minimal DFJK E_v forward-limit model
# ---------------------------------------------------------------------------
#
# We deliberately start ONLY at xi=0, t=0:
#
#   E_v^q(x,0,0) = e_v^q(x)
#                 = N_q x^{-alpha} (1-x)^{beta_q}.
#
# This is the forward-limit ansatz used by Diehl, Feldmann, Jakob and Kroll
# (EPJC 39 (2005) 1, hep-ph/0408173).  Their normalization is fixed by the
# proton/neutron anomalous magnetic moments:
#
#   integral_0^1 dx e_v^u = kappa_u,
#   integral_0^1 dx e_v^d = kappa_d.
#
# For this workshop projection we fix alpha=0.55, the value used in the DFJK
# E_v study, and allow only beta_u and beta_d to control the x-shape.
#
# IMPORTANT: this step does NOT yet add t dependence, skewness xi, a
# double-distribution profile, CFFs, or DVCS observables.  Those belong to
# later steps and should only be added after this simple map is understood.
DFJK_ALPHA_E = 0.55

# Transparent workshop reference point, chosen near the middle of the
# beta ranges explored by DFJK.  It is NOT advertised as a new fit or as the
# exact DFJK best-fit parameter set.
WORKSHOP_BETA_U_E = 4.0
WORKSHOP_BETA_D_E = 6.0

# Scan ranges chosen to cover the beta variations explicitly discussed/shown
# in the DFJK E_v study (u roughly 3.5--5.5, d roughly 5--9).
BETA_U_E_SCAN = (3.5, 5.5)
BETA_D_E_SCAN = (5.0, 9.0)

# ---------------------------------------------------------------------------
# Stage-5 Step 2A: add ONLY t dependence
# ---------------------------------------------------------------------------
#
# Before introducing skewness xi, use the simplest possible factorized
# extension away from t=0:
#
#   E_v^q(x, xi=0, t) = E_v^q(x,0,0) * exp(B_E * t).
#
# Physical DVCS kinematics have t < 0, so positive B_E suppresses the GPD as
# |t| grows.  A single common slope is used for u and d on purpose: at this
# stage it is only a transparent bridge from the forward limit to nonzero t,
# not a precision GPD parameterization.
WORKSHOP_E_T_SLOPE = 1.0  # GeV^-2; illustrative fixed workshop value

RGB_XS_DEFAULT = Path("import/ndvcs_clas12_preliminary_unpolarized.txt")
RGB_BSA_DIR_DEFAULT = None  # resolved relative to this script in main()
RGA_PASS2_DEFAULT = None  # resolved relative to this script below
OUT_DEFAULT = Path("output/stage5_ji")
A20_Q2_GEV2 = 3.0
A20_SUMMARY_DEFAULT = Path("output/stage5_ji/tables/nnpdf40_valence_A20_summary_Q2_3.csv")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def beta_fn(a: float, b: float) -> float:
    return math.gamma(a) * math.gamma(b) / math.gamma(a + b)


def ev_normalization(kappa: float, alpha: float, beta: float) -> float:
    """N_q such that integral_0^1 e_v^q(x) dx = kappa_q."""
    return kappa / beta_fn(1.0 - alpha, beta + 1.0)


def ev_forward(x, kappa: float, alpha: float, beta: float):
    """Minimal DFJK forward limit e_v^q(x) = E_v^q(x, xi=0, t=0)."""
    x = np.asarray(x, dtype=float)
    Nq = ev_normalization(kappa, alpha, beta)
    return Nq * np.power(x, -alpha) * np.power(1.0 - x, beta)


def ev_zero_skewness_t(x, t: float, kappa: float, alpha: float, beta: float,
                       slope: float = WORKSHOP_E_T_SLOPE):
    """Step 2A: simplest nonzero-t extension at xi=0.

    E_v^q(x,0,t) = E_v^q(x,0,0) exp(slope*t).

    t is in GeV^2 and is negative in the physical region.
    """
    return ev_forward(x, kappa, alpha, beta) * np.exp(slope * t)


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
    """Stage-5 Step 1: show directly what beta_q does to E_v^q(x).

    This is intentionally NOT a DVCS fit.  Each curve has the same zeroth
    moment (fixed by kappa_q); changing beta_q only redistributes that fixed
    normalization in x, thereby changing the x-weighted B20 moment.
    """
    x = np.linspace(0.01, 0.99, 600)

    beta_u_values = [3.5, WORKSHOP_BETA_U_E, 5.5]
    beta_d_values = [5.0, WORKSHOP_BETA_D_E, 9.0]

    rows = []
    for flavor, kappa, beta_values in [
        ("u_v", KAPPA_U, beta_u_values),
        ("d_v", KAPPA_D, beta_d_values),
    ]:
        for beta in beta_values:
            b20 = b20_from_beta(kappa, DFJK_ALPHA_E, beta)
            norm = ev_normalization(kappa, DFJK_ALPHA_E, beta)
            rows.append({
                "flavor": flavor,
                "alpha_E": DFJK_ALPHA_E,
                "beta_E": beta,
                "N_q": norm,
                "zeroth_moment_kappa": kappa,
                "B20_qv": b20,
            })

    summary = pd.DataFrame(rows)
    summary.to_csv(tabdir / "dfjk_step1_Ev_shape_summary.csv", index=False)

    assumptions = pd.DataFrame([
        ["alpha_E", DFJK_ALPHA_E,
         "fixed; DFJK E_v forward-limit choice"],
        ["beta_u_E_reference", WORKSHOP_BETA_U_E,
         "workshop reference; not fitted yet"],
        ["beta_d_E_reference", WORKSHOP_BETA_D_E,
         "workshop reference; not fitted yet"],
        ["kappa_u", KAPPA_U,
         "fixes integral E_v^u dx"],
        ["kappa_d", KAPPA_D,
         "fixes integral E_v^d dx"],
    ], columns=["quantity", "value", "role"])
    assumptions.to_csv(tabdir / "dfjk_step1_model_assumptions.csv", index=False)

    # Two separate figures rather than subplots: each has one simple message.
    for flavor, latex_flavor, kappa, beta_values in [
        ("u", "u", KAPPA_U, beta_u_values),
        ("d", "d", KAPPA_D, beta_d_values),
    ]:
        fig, ax = plt.subplots(figsize=(8.0, 5.4))

        for beta in beta_values:
            y = ev_forward(x, kappa, DFJK_ALPHA_E, beta)
            b20 = b20_from_beta(kappa, DFJK_ALPHA_E, beta)
            is_ref = (
                (flavor == "u" and np.isclose(beta, WORKSHOP_BETA_U_E))
                or
                (flavor == "d" and np.isclose(beta, WORKSHOP_BETA_D_E))
            )
            label = (
                rf"$\beta_{latex_flavor}^E={beta:.1f}$"
                + rf"  [$B_{{20}}^{{{latex_flavor}_v}}={b20:.3f}$]"
            )
            ax.plot(
                x, y,
                linewidth=3.0 if is_ref else 2.0,
                linestyle="-" if is_ref else "--",
                label=label,
            )

        ax.axhline(0.0, linewidth=0.8, color="black")
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(rf"$E_v^{latex_flavor}(x,0,0)$")
        ax.set_title(
            rf"Step 1: changing $\beta_{latex_flavor}^E$ changes the "
            rf"$x$-shape of $E_v^{latex_flavor}$"
        )
        ax.text(
            0.98, 0.96,
            rf"All curves satisfy $\int_0^1 E_v^{latex_flavor}(x)\,dx"
            rf"={kappa:+.3f}$",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
        )
        ax.grid(alpha=0.22)
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(
            figdir / f"dfjk_step1_Ev_{flavor}_shape.png",
            dpi=180,
        )
        plt.close(fig)

    return summary


def make_step2a_t_dependence(figdir: Path, tabdir: Path):
    """Step 2A diagnostic: add t dependence, but still keep xi=0.

    The goal is pedagogical: show exactly what the new variable t does before
    introducing skewness or calculating a DVCS observable.
    """
    x = np.linspace(0.01, 0.99, 600)
    t_values = [0.0, -0.3, -0.6, -0.9]

    rows = []
    for flavor, kappa, beta in [
        ("u_v", KAPPA_U, WORKSHOP_BETA_U_E),
        ("d_v", KAPPA_D, WORKSHOP_BETA_D_E),
    ]:
        b20_0 = b20_from_beta(kappa, DFJK_ALPHA_E, beta)
        for t in t_values:
            scale = float(np.exp(WORKSHOP_E_T_SLOPE * t))
            rows.append({
                "flavor": flavor,
                "t_GeV2": t,
                "beta_E": beta,
                "B_E_GeV_minus2": WORKSHOP_E_T_SLOPE,
                "exp_Bt": scale,
                "integral_Ev_dx": kappa * scale,
                "B20_qv_t": b20_0 * scale,
            })

    pd.DataFrame(rows).to_csv(
        tabdir / "dfjk_step2a_t_dependence_summary.csv", index=False
    )

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2))

    for ax, flavor, latex_flavor, kappa, beta in [
        (axes[0], "u", "u", KAPPA_U, WORKSHOP_BETA_U_E),
        (axes[1], "d", "d", KAPPA_D, WORKSHOP_BETA_D_E),
    ]:
        for t in t_values:
            y = ev_zero_skewness_t(
                x, t, kappa, DFJK_ALPHA_E, beta, WORKSHOP_E_T_SLOPE
            )
            ax.plot(x, y, linewidth=2.2, label=rf"$t={t:.1f}\ \mathrm{{GeV}}^2$")

        ax.axhline(0.0, linewidth=0.8, color="black")
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(rf"$E_v^{latex_flavor}(x,\xi=0,t)$")
        ax.set_title(rf"$E_v^{latex_flavor}$")
        ax.grid(alpha=0.22)
        ax.legend(fontsize=8)

    fig.suptitle(
        rf"Step 2A: adding $t$ dependence at $\xi=0$: "
        rf"$E_v(x,0,t)=E_v(x,0,0)e^{{B_Et}}$, "
        rf"$B_E={WORKSHOP_E_T_SLOPE:.1f}\ \mathrm{{GeV}}^{{-2}}$",
        fontsize=15,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(figdir / "dfjk_step2a_Ev_t_dependence.png", dpi=180)
    plt.close(fig)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Stage-5 Step 2B: introduce skewness as a redistribution variable
# ---------------------------------------------------------------------------
#
# A Radyushkin-style double distribution introduces a second variable alpha.
# For a fixed forward momentum fraction beta, alpha describes how longitudinal
# momentum transfer is shared.  We start by plotting ONLY the normalized
# profile pi_b(beta, alpha), before integrating it into E(x,xi,t).
#
# Standard profile:
#   pi_b(beta,alpha) =
#     Gamma(2b+2) / [2^(2b+1) Gamma(b+1)^2]
#     * [((1-beta)^2-alpha^2)^b / (1-beta)^(2b+1)]
#
# for 0 <= beta <= 1 and |alpha| <= 1-beta.
#
# Its key property is integral d alpha pi_b = 1.  Thus it redistributes a
# given forward E_v(beta) in the new alpha direction without changing its
# total weight.  Only in the NEXT step will xi connect beta and alpha through
# x = beta + xi*alpha.
WORKSHOP_DD_PROFILE_B = 1.0


def dd_profile(beta, alpha, b=WORKSHOP_DD_PROFILE_B):
    """Normalized Radyushkin DD profile for valence beta in [0,1]."""
    beta = float(beta)
    alpha = np.asarray(alpha, dtype=float)
    one_minus = 1.0 - beta
    pref = (
        math.gamma(2.0 * b + 2.0)
        / (2.0 ** (2.0 * b + 1.0) * math.gamma(b + 1.0) ** 2)
    )
    out = np.zeros_like(alpha)
    mask = np.abs(alpha) <= one_minus
    core = np.maximum(one_minus**2 - alpha[mask]**2, 0.0)
    out[mask] = pref * core**b / one_minus**(2.0 * b + 1.0)
    return out


def make_step2b_skewness_profile(figdir: Path, tabdir: Path):
    """Pedagogical Step 2B: show the normalized DD profile only.

    We intentionally stop before constructing E(x,xi,t).  This isolates the
    one new idea: skewness requires a second longitudinal-momentum variable.
    """
    beta_values = [0.1, 0.3, 0.5]
    rows = []

    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    for beta in beta_values:
        amax = 1.0 - beta
        alpha = np.linspace(-amax, amax, 600)
        profile = dd_profile(beta, alpha, WORKSHOP_DD_PROFILE_B)
        integral = np.trapezoid(profile, alpha)

        ax.plot(
            alpha, profile, linewidth=2.2,
            label=rf"$\beta={beta:.1f}$  [$\int d\alpha\,\pi={integral:.3f}$]"
        )

        for a, p in zip(alpha, profile):
            rows.append({
                "beta": beta,
                "alpha": a,
                "profile_pi": p,
                "profile_b": WORKSHOP_DD_PROFILE_B,
            })

    ax.set_xlabel(r"$\alpha$ (longitudinal momentum-transfer sharing)")
    ax.set_ylabel(r"$\pi_b(\beta,\alpha)$")
    ax.set_title(
        r"Step 2B: skewness starts by spreading each forward $\beta$ "
        r"over a second variable $\alpha$"
    )
    ax.grid(alpha=0.22)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(figdir / "dfjk_step2b_skewness_profile.png", dpi=180)
    plt.close(fig)

    pd.DataFrame(rows).to_csv(
        tabdir / "dfjk_step2b_skewness_profile.csv", index=False
    )


def make_step2c_beta_to_x_mapping(figdir: Path, tabdir: Path):
    """Step 2C: visualize x = beta + xi*alpha for one beta.

    This still does NOT construct the full GPD E(x,xi,t).  It follows one
    forward beta contribution through the DD profile and shows how finite xi
    spreads that one contribution over x.
    """
    beta = 0.30
    xi = 0.20
    b = WORKSHOP_DD_PROFILE_B

    amax = 1.0 - beta
    alpha = np.linspace(-amax, amax, 800)
    pi_alpha = dd_profile(beta, alpha, b)

    # Mapping from the DD variable alpha to the GPD integration variable x.
    x = beta + xi * alpha

    # Probability-density transformation:
    # p_x(x) dx = pi(alpha) d alpha, with dx = xi d alpha.
    p_x = pi_alpha / abs(xi)

    int_alpha = np.trapezoid(pi_alpha, alpha)
    int_x = np.trapezoid(p_x, x)

    # Save the exact pedagogical mapping used in the figure.
    pd.DataFrame({
        "beta": beta,
        "xi": xi,
        "alpha": alpha,
        "x_equals_beta_plus_xi_alpha": x,
        "pi_beta_alpha": pi_alpha,
        "mapped_density_in_x": p_x,
    }).to_csv(tabdir / "dfjk_step2c_beta_to_x_mapping.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2))

    axes[0].plot(alpha, pi_alpha, linewidth=2.5)
    axes[0].axvline(0.0, linewidth=0.9, color="black", alpha=0.7)
    axes[0].set_xlabel(r"$\alpha$")
    axes[0].set_ylabel(r"$\pi_b(\beta,\alpha)$")
    axes[0].set_title(
        rf"Start with one forward contribution: $\beta={beta:.2f}$"
    )
    axes[0].text(
        0.04, 0.94,
        rf"$\int d\alpha\,\pi={int_alpha:.3f}$",
        transform=axes[0].transAxes, ha="left", va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )
    axes[0].grid(alpha=0.22)

    axes[1].plot(x, p_x, linewidth=2.5)
    axes[1].axvline(beta, linewidth=1.2, color="black", linestyle="--",
                    label=rf"center $x=\beta={beta:.2f}$")
    axes[1].set_xlabel(r"$x=\beta+\xi\alpha$")
    axes[1].set_ylabel(r"mapped weight density in $x$")
    axes[1].set_title(
        rf"Turn on $\xi={xi:.2f}$: that contribution spreads over $x$"
    )
    axes[1].text(
        0.04, 0.94,
        rf"$x\in[{x.min():.2f},{x.max():.2f}]$" + "\n"
        + rf"$\int dx\,p_x={int_x:.3f}$",
        transform=axes[1].transAxes, ha="left", va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.22)

    fig.suptitle(
        r"Step 2C: finite skewness maps the DD variables through "
        r"$x=\beta+\xi\alpha$",
        fontsize=15,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(figdir / "dfjk_step2c_beta_to_x_mapping.png", dpi=180)
    plt.close(fig)

    return {
        "beta": beta,
        "xi": xi,
        "x_min": float(x.min()),
        "x_max": float(x.max()),
        "integral_alpha": float(int_alpha),
        "integral_x": float(int_x),
    }


# ---------------------------------------------------------------------------
# Production finite-skewness E_v model
# ---------------------------------------------------------------------------
#
# The pedagogical alpha/beta plots above are retained as helper functions but
# are no longer produced by main().  From here on Stage 5 uses the complete
# DD integral internally:
#
#   E_v^q(x,xi,t) = integral d beta E_v^q(beta,0,t)
#                    pi_b(beta, alpha=(x-beta)/xi) / |xi|
#
# over the support |alpha| <= 1-beta.
#
# This is still a deliberately simplified workshop model:
#   * valence E only;
#   * common fixed factorized t slope;
#   * fixed DD profile b=1;
#   * beta_u and beta_d are the shape parameters of interest.
#
# The model is used internally; alpha/beta are not conference observables.


def _positive_unit_grid(n: int = 1800):
    """Grid on 0<z<1 with strong resolution near z=0.

    E_v ~ z^{-alpha_E}, so a uniform grid badly under-resolves the small-z
    contribution to its moments.  This composite grid fixes that numerical
    problem without changing the model.
    """
    n_log = max(200, n // 2)
    n_lin = max(200, n - n_log)
    return np.unique(np.concatenate([
        np.geomspace(1.0e-9, 0.03, n_log),
        np.linspace(0.03, 1.0 - 1.0e-8, n_lin),
    ]))


def _dd_profile_vectorized(beta, alpha, b=WORKSHOP_DD_PROFILE_B):
    """Vectorized equivalent of dd_profile for numerical DD integration."""
    beta = np.asarray(beta, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    one_minus = 1.0 - beta
    pref = (
        math.gamma(2.0 * b + 2.0)
        / (2.0 ** (2.0 * b + 1.0) * math.gamma(b + 1.0) ** 2)
    )
    core = one_minus**2 - alpha**2
    out = np.zeros_like(beta)
    mask = np.abs(alpha) <= one_minus
    out[mask] = (
        pref
        * np.maximum(core[mask], 0.0) ** b
        / one_minus[mask] ** (2.0 * b + 1.0)
    )
    return out


def ev_finite_xi(x, xi: float, t: float, kappa: float, alpha_e: float,
                 beta_e: float, slope: float = WORKSHOP_E_T_SLOPE,
                 profile_b: float = WORKSHOP_DD_PROFILE_B,
                 n_beta: int = 1800):
    """Finite-skewness valence E from the fixed-profile DD construction."""
    xarr = np.atleast_1d(np.asarray(x, dtype=float))

    # xi -> 0 must reduce to the zero-skewness input.
    if abs(xi) < 1.0e-8:
        out = np.zeros_like(xarr)
        mask = (xarr > 0.0) & (xarr < 1.0)
        out[mask] = ev_zero_skewness_t(
            xarr[mask], t, kappa, alpha_e, beta_e, slope
        )
        return float(out[0]) if np.ndim(x) == 0 else out

    beta_grid = _positive_unit_grid(n_beta)
    forward = ev_zero_skewness_t(
        beta_grid, t, kappa, alpha_e, beta_e, slope
    )

    out = np.zeros_like(xarr)
    for ix, xv in enumerate(xarr):
        alpha_dd = (xv - beta_grid) / xi
        prof = _dd_profile_vectorized(beta_grid, alpha_dd, profile_b)
        out[ix] = np.trapezoid(
            forward * prof / abs(xi), beta_grid
        )

    return float(out[0]) if np.ndim(x) == 0 else out


def validate_finite_xi_model(tabdir: Path):
    """Numerically verify the DD zeroth-moment sum rule.

    Two corrections relative to v21 are essential:
      1. resolve the integrable x^{-alpha_E} behavior near x=0;
      2. at finite xi integrate the full DD support -xi <= x <= 1,
         including the negative-x part of the ERBL region.

    No conference figure is produced; this is an internal numerical check.
    """
    rows = []

    # Dense near x=0 so the forward singular behavior is integrated correctly.
    x_forward = _positive_unit_grid(4000)

    for flavor, kappa, beta_e in [
        ("u_v", KAPPA_U, WORKSHOP_BETA_U_E),
        ("d_v", KAPPA_D, WORKSHOP_BETA_D_E),
    ]:
        for t in [0.0, -0.3, -0.6]:
            expected_zeroth = kappa * np.exp(WORKSHOP_E_T_SLOPE * t)

            e0 = ev_zero_skewness_t(
                x_forward, t, kappa, DFJK_ALPHA_E, beta_e,
                WORKSHOP_E_T_SLOPE
            )
            zeroth_xi0 = np.trapezoid(e0, x_forward)

            for xi in [0.05, 0.10, 0.20, 0.30]:
                # For beta in [0,1] and |alpha|<=1-beta, the finite-xi
                # construction has support down to x=-xi and up to x=1.
                x_negative = np.linspace(-xi, 0.0, 500, endpoint=False)
                x_positive = np.unique(np.concatenate([
                    np.geomspace(1.0e-8, 0.03, 700),
                    np.linspace(0.03, 1.0, 900),
                ]))
                x_full = np.concatenate([x_negative, x_positive])

                exi = ev_finite_xi(
                    x_full, xi, t, kappa, DFJK_ALPHA_E, beta_e,
                    WORKSHOP_E_T_SLOPE, WORKSHOP_DD_PROFILE_B,
                    n_beta=1800
                )
                zeroth = np.trapezoid(exi, x_full)

                rows.append({
                    "flavor": flavor,
                    "t_GeV2": t,
                    "xi": xi,
                    "expected_zeroth_moment": expected_zeroth,
                    "numerical_xi0_moment": zeroth_xi0,
                    "xi0_over_expected": zeroth_xi0 / expected_zeroth,
                    "numerical_finite_xi_moment": zeroth,
                    "finite_xi_over_expected": zeroth / expected_zeroth,
                    "finite_xi_fractional_deviation":
                        zeroth / expected_zeroth - 1.0,
                })

    out = pd.DataFrame(rows)
    out.to_csv(tabdir / "finite_xi_model_validation.csv", index=False)

    max_forward_dev = float(
        np.max(np.abs(out["xi0_over_expected"] - 1.0))
    )
    max_finite_dev = float(
        np.max(np.abs(out["finite_xi_over_expected"] - 1.0))
    )

    # This should be a numerical identity at the sub-percent level.  Fail
    # loudly rather than allowing later observable fits to use a broken grid.
    tolerance = 2.0e-3
    if max_forward_dev > tolerance or max_finite_dev > tolerance:
        raise RuntimeError(
            "Finite-xi DD validation failed: "
            f"max forward deviation={max_forward_dev:.4g}, "
            f"max finite-xi deviation={max_finite_dev:.4g}, "
            f"tolerance={tolerance:.4g}"
        )

    return out



def write_stage5_projection_plan(tabdir: Path):
    """Machine-readable statement of the only conference comparison we want.

    No fake contour is produced here.  The contour must come from an
    observable-level proton/neutron DVCS fit, not from assigning an arbitrary
    conversion between sigma(A_LU) and beta_q.
    """
    plan = pd.DataFrame([
        ["p_1x", "proton", 1.0,
         "current proton XS+BSA statistical precision"],
        ["p1x_plus_n1x", "proton+neutron", 1.0,
         "current proton plus current/projected neutron information"],
        ["p10x_plus_n10x", "proton+neutron", 10.0,
         "high-luminosity proton+neutron projection"],
    ], columns=["scenario", "targets", "luminosity_factor", "meaning"])
    plan.to_csv(tabdir / "b20_conference_projection_scenarios.csv", index=False)
    return plan



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

    # Main-talk distribution plot:
    #   - current RGA at its full 1444-point 4D granularity;
    #   - RGB at 2x luminosity with FOUR RGA-like points combined per projected
    #     neutron bin (~1/4 as many bins);
    #   - RGB at 10x luminosity at the full RGA-like 1444-point granularity.
    #
    # For the 2x/coarser case, combine consecutive groups of four entries after
    # sorting by the actual RGA kinematic coordinates. Statistical information
    # adds, so sigma_group = 1/sqrt(sum_k 1/sigma_k^2). This preserves the total
    # 2x RGB information budget while spending it on ~361 rather than 1444 bins.
    rgb2_full = points[
        points["projection_luminosity"] == 2
    ].copy()
    rgb10 = points[
        points["projection_luminosity"] == 10
    ]["rgb_projected_stat_error"].to_numpy(float)

    sort_cols = [c for c in ["bin", "xB", "Q2", "t_abs", "phi_deg"] if c in rgb2_full.columns]
    if sort_cols:
        rgb2_full = rgb2_full.sort_values(sort_cols).reset_index(drop=True)

    rgb2_sigma_full = rgb2_full["rgb_projected_stat_error"].to_numpy(float)
    rgb2_coarse = []
    for i in range(0, len(rgb2_sigma_full), 4):
        group = rgb2_sigma_full[i:i + 4]
        group = group[np.isfinite(group) & (group > 0)]
        if len(group):
            rgb2_coarse.append(1.0 / np.sqrt(np.sum(1.0 / group**2)))
    rgb2_coarse = np.asarray(rgb2_coarse, dtype=float)

    # Common bin edges for a direct visual comparison.
    bins = np.linspace(
        0.0,
        max(
            np.quantile(pstat, 0.99),
            np.quantile(rgb2_coarse, 0.99),
            np.quantile(rgb10, 0.99),
        ) * 1.08,
        36,
    )

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.hist(
        pstat, bins=bins, histtype="step", linewidth=2.6, color="black",
        label=f"Current RGA proton 4D (N={len(pstat)})"
    )
    ax.hist(
        rgb2_coarse, bins=bins, histtype="step", linewidth=2.2,
        label=f"RGB neutron 2x, 1/4 RGA granularity (N={len(rgb2_coarse)})"
    )
    ax.hist(
        rgb10, bins=bins, histtype="step", linewidth=2.2,
        label=f"RGB neutron 10x, full RGA granularity (N={len(rgb10)})"
    )

    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax.axvline(
        np.median(pstat), color="black", linestyle=":", linewidth=2.0,
        label=f"RGA median = {np.median(pstat):.3f}"
    )
    ax.axvline(
        np.median(rgb2_coarse), color=cycle[0], linestyle="--", linewidth=1.8,
        label=f"RGB 2x coarse median = {np.median(rgb2_coarse):.3f}"
    )
    ax.axvline(
        np.median(rgb10), color=cycle[1], linestyle="--", linewidth=1.8,
        label=f"RGB 10x median = {np.median(rgb10):.3f}"
    )

    ax.set_xlabel(r"Statistical $\sigma(A_{LU})$")
    ax.set_ylabel("Number of bins")
    ax.set_title("Neutron BSA projected onto RGA-like multidimensional binning")
    ax.grid(alpha=0.20)
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    fig.savefig(figdir / "rgb_10x_on_rga_4d_binning_precision.png", dpi=180)
    plt.close(fig)

    # Audit the coarser 2x scenario explicitly.
    pd.DataFrame({
        "rgb_2x_coarse_stat_error": rgb2_coarse
    }).to_csv(tabdir / "rgb_2x_quarter_rga_granularity.csv", index=False)

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
# Stage-5 production p+n H+E projection
# ---------------------------------------------------------------------------
#
# This is the first conference-facing fit.  beta_u and beta_d are global E_v
# shape parameters.  ReH and ImH are NOT fixed: every proton/neutron kinematic
# cell receives its own local H nuisance pair, constrained by the same XS+BSA
# rows used in the fit.  Thus the 1x -> 10x improvement of H is propagated
# automatically rather than imposed as an external percentage prior.
#
# E is mapped to ReE/ImE at leading order from the finite-xi DD model.  Gepard
# is then used only for the observable response (BH+DVCS+interference) to the
# resulting CFF changes.  This is intentionally a workshop projection, not a
# global GPD extraction.

STAGE5_BETA_STEP_U = 0.08
STAGE5_BETA_STEP_D = 0.10
STAGE5_MIN_NPHI = 3


def _cff_owner_stage5(th, name):
    candidates = [th.m] if hasattr(th, "m") else []
    candidates.append(th)
    for obj in candidates:
        if hasattr(obj, name) and callable(getattr(obj, name)):
            return obj
    raise AttributeError(f"Could not find callable Gepard CFF method {name}")


@contextlib.contextmanager
def _shifted_cffs_stage5(th, shifts):
    saved = []
    try:
        for name, shift in shifts.items():
            owner = _cff_owner_stage5(th, name)
            old = getattr(owner, name)
            setattr(owner, name, lambda pt, _old=old, _s=float(shift):
                    float(_old(pt)) + _s)
            saved.append((owner, name, old))
        yield
    finally:
        for owner, name, old in reversed(saved):
            setattr(owner, name, old)


def _make_gepard_point_stage5(g, xB, Q2, t_abs, phi_deg, ebeam,
                              target="p", helicity=0):
    pt = g.DataPoint(
        xB=float(xB), t=-abs(float(t_abs)), Q2=float(Q2),
        phi=math.pi - math.radians(float(phi_deg)),
        observable="XS", frame="trento", process="ep2epgamma",
        exptype="fixed target", in1energy=float(ebeam), in1charge=-1,
        in1polarization=int(helicity), in2particle=str(target),
    )
    pt.prepare()
    return pt


def _bsa_gepard_stage5(th, pp, pm):
    sp, sm = float(th.predict(pp)), float(th.predict(pm))
    return (sp-sm)/(sp+sm) if np.isfinite(sp+sm) and abs(sp+sm)>1e-30 else np.nan


def _lo_flavor_E_cff(beta_u, beta_d, xB, t_abs, target):
    """Workshop LO CFF E from the finite-xi valence DD model.

    Convention used internally:
      Re E = PV int_-1^1 dx E(x,xi,t)[1/(xi-x)-1/(xi+x)]
      Im E = pi [E(xi,xi,t)-E(-xi,xi,t)].

    Proton charge weights are 4/9 u + 1/9 d; neutron interchanges u,d by
    isospin.  Only *changes* around the reference beta point drive the
    projected covariance, which reduces sensitivity to an overall CFF sign
    convention.
    """
    xi = float(xB)/(2.0-float(xB))
    t = -abs(float(t_abs))
    # Dense full support; avoid evaluating exactly at +/-xi in the subtracted
    # principal-value integrands.
    xneg = np.linspace(-xi, -1.0e-6, 220, endpoint=True)
    xpos = np.unique(np.concatenate([
        np.geomspace(1.0e-6, 0.03, 260),
        np.linspace(0.03, 1.0-1.0e-6, 520),
    ]))
    x = np.unique(np.concatenate([xneg, xpos]))

    def one(kappa, beta_e):
        E = ev_finite_xi(x, xi, t, kappa, DFJK_ALPHA_E, beta_e,
                         WORKSHOP_E_T_SLOPE, WORKSHOP_DD_PROFILE_B,
                         n_beta=1100)
        Ep = float(ev_finite_xi(xi, xi, t, kappa, DFJK_ALPHA_E, beta_e,
                                WORKSHOP_E_T_SLOPE, WORKSHOP_DD_PROFILE_B,
                                n_beta=1400))
        Em = float(ev_finite_xi(-xi, xi, t, kappa, DFJK_ALPHA_E, beta_e,
                                WORKSHOP_E_T_SLOPE, WORKSHOP_DD_PROFILE_B,
                                n_beta=1400))
        # Subtracted PV forms are finite at x=+/-xi.
        d1 = xi-x
        d2 = xi+x
        f1 = np.where(np.abs(d1)>2e-6, (E-Ep)/d1, 0.0)
        f2 = np.zeros_like(E, dtype=float)
        np.divide(E - Em, d2, out=f2, where=np.abs(d2) > 2e-6)
        logpv = math.log((1.0+xi)/(1.0-xi))
        re = np.trapezoid(f1, x) + Ep*logpv \
             - np.trapezoid(f2, x) - Em*logpv
        im = math.pi*(Ep-Em)
        return float(re), float(im)

    ru, iu = one(KAPPA_U, beta_u)
    rd, id_ = one(KAPPA_D, beta_d)
    if str(target).lower().startswith("n"):
        return (4.0*rd+ru)/9.0, (4.0*id_+iu)/9.0
    return (4.0*ru+rd)/9.0, (4.0*iu+id_)/9.0


def _load_stage2_proton_inputs(stage2_dir: Path, factor: int):
    p = stage2_dir / f"joint_fit_input_km15_{factor}x.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"Stage-5 H+E fit needs the existing Stage-2 input {p}. "
            "Run high_luminosity_dvcs_stage2_pseudodata_pass2.py first."
        )
    d = pd.read_csv(p)
    d["bin"] = d["bin"].astype(int)
    return d.sort_values(["bin","phi_deg"]).reset_index(drop=True)


def _build_he_observable_rows(th, g, proton, rgb_xs, rgb_bsa, factor):
    """Build linearized p+n observable response rows around the reference model."""
    obs = []
    cell_cff_cache = {}

    def cff_triplet(xB,Q2,t_abs,target):
        key=(target,round(float(xB),6),round(float(Q2),6),round(float(t_abs),6))
        if key not in cell_cff_cache:
            ref = _lo_flavor_E_cff(WORKSHOP_BETA_U_E, WORKSHOP_BETA_D_E,
                                   xB,t_abs,target)
            up_p = _lo_flavor_E_cff(WORKSHOP_BETA_U_E+STAGE5_BETA_STEP_U,
                                     WORKSHOP_BETA_D_E,xB,t_abs,target)
            up_m = _lo_flavor_E_cff(WORKSHOP_BETA_U_E-STAGE5_BETA_STEP_U,
                                     WORKSHOP_BETA_D_E,xB,t_abs,target)
            dn_p = _lo_flavor_E_cff(WORKSHOP_BETA_U_E,
                                     WORKSHOP_BETA_D_E+STAGE5_BETA_STEP_D,xB,t_abs,target)
            dn_m = _lo_flavor_E_cff(WORKSHOP_BETA_U_E,
                                     WORKSHOP_BETA_D_E-STAGE5_BETA_STEP_D,xB,t_abs,target)
            cell_cff_cache[key]=(ref,up_p,up_m,dn_p,dn_m)
        return cell_cff_cache[key]

    def prediction(xB,Q2,t_abs,phi,ebeam,target,helicity,beta_u,beta_d):
        p0=_make_gepard_point_stage5(g,xB,Q2,t_abs,phi,ebeam,target,helicity)
        re0=float(getattr(_cff_owner_stage5(th,"ReE"),"ReE")(p0))
        im0=float(getattr(_cff_owner_stage5(th,"ImE"),"ImE")(p0))
        rem,imm=_lo_flavor_E_cff(beta_u,beta_d,xB,t_abs,target)
        with _shifted_cffs_stage5(th,{"ReE":rem-re0,"ImE":imm-im0}):
            return float(th.predict(p0))

    def bsa_prediction(xB,Q2,t_abs,phi,ebeam,target,beta_u,beta_d):
        pp=_make_gepard_point_stage5(g,xB,Q2,t_abs,phi,ebeam,target,+1)
        pm=_make_gepard_point_stage5(g,xB,Q2,t_abs,phi,ebeam,target,-1)
        re0=float(getattr(_cff_owner_stage5(th,"ReE"),"ReE")(pp))
        im0=float(getattr(_cff_owner_stage5(th,"ImE"),"ImE")(pp))
        rem,imm=_lo_flavor_E_cff(beta_u,beta_d,xB,t_abs,target)
        with _shifted_cffs_stage5(th,{"ReE":rem-re0,"ImE":imm-im0}):
            return _bsa_gepard_stage5(th,pp,pm)

    def derivatives(xB,Q2,t_abs,phi,ebeam,target,kind):
        if kind=="xs":
            fun=lambda bu,bd: prediction(xB,Q2,t_abs,phi,ebeam,target,0,bu,bd)
        else:
            fun=lambda bu,bd: bsa_prediction(xB,Q2,t_abs,phi,ebeam,target,bu,bd)
        y0=fun(WORKSHOP_BETA_U_E,WORKSHOP_BETA_D_E)
        du=(fun(WORKSHOP_BETA_U_E+STAGE5_BETA_STEP_U,WORKSHOP_BETA_D_E)-
            fun(WORKSHOP_BETA_U_E-STAGE5_BETA_STEP_U,WORKSHOP_BETA_D_E))/(2*STAGE5_BETA_STEP_U)
        dd=(fun(WORKSHOP_BETA_U_E,WORKSHOP_BETA_D_E+STAGE5_BETA_STEP_D)-
            fun(WORKSHOP_BETA_U_E,WORKSHOP_BETA_D_E-STAGE5_BETA_STEP_D))/(2*STAGE5_BETA_STEP_D)
        # Local H derivatives around the same E-reference pseudo-truth.
        p0=_make_gepard_point_stage5(g,xB,Q2,t_abs,phi,ebeam,target,0)
        reh=float(getattr(_cff_owner_stage5(th,"ReH"),"ReH")(p0))
        imh=float(getattr(_cff_owner_stage5(th,"ImH"),"ImH")(p0))
        hs_re=max(0.05,0.02*max(abs(reh),1.0)); hs_im=max(0.05,0.02*max(abs(imh),1.0))
        def h_eval(name,step):
            if kind=="xs":
                # Preserve the reference E replacement while shifting H.
                re0=float(getattr(_cff_owner_stage5(th,"ReE"),"ReE")(p0)); im0=float(getattr(_cff_owner_stage5(th,"ImE"),"ImE")(p0))
                rem,imm=_lo_flavor_E_cff(WORKSHOP_BETA_U_E,WORKSHOP_BETA_D_E,xB,t_abs,target)
                with _shifted_cffs_stage5(th,{"ReE":rem-re0,"ImE":imm-im0,name:step}):
                    return float(th.predict(p0))
            pp=_make_gepard_point_stage5(g,xB,Q2,t_abs,phi,ebeam,target,+1)
            pm=_make_gepard_point_stage5(g,xB,Q2,t_abs,phi,ebeam,target,-1)
            re0=float(getattr(_cff_owner_stage5(th,"ReE"),"ReE")(pp)); im0=float(getattr(_cff_owner_stage5(th,"ImE"),"ImE")(pp))
            rem,imm=_lo_flavor_E_cff(WORKSHOP_BETA_U_E,WORKSHOP_BETA_D_E,xB,t_abs,target)
            with _shifted_cffs_stage5(th,{"ReE":rem-re0,"ImE":imm-im0,name:step}):
                return _bsa_gepard_stage5(th,pp,pm)
        dre=(h_eval("ReH",+hs_re)-h_eval("ReH",-hs_re))/(2*hs_re)
        dim=(h_eval("ImH",+hs_im)-h_eval("ImH",-hs_im))/(2*hs_im)
        return y0,du,dd,dre,dim

    # Proton: actual Stage-2 matched 4D XS+BSA rows.  H is refit locally here;
    # no Stage-2 H prior is added, avoiding double counting the same data.
    for row in proton.to_dict("records"):
        cell = f"p:{int(row['bin'])}"
        for kind in ("xs", "bsa"):
            y0, du, dd, dhre, dhim = derivatives(
                row["xB"], row["Q2"], row["t_abs"], row["phi_deg"],
                row["ebeam"], "p", kind
            )
            if kind == "xs":
                # Stage-2's joint_fit_input file contains the experimental
                # absolute XS uncertainties, but not the KM15 prediction.
                # Keep those uncertainties in their native XS units, exactly
                # as Stage-2 CFF sensitivity does.  Do not reconstruct a
                # relative uncertainty using an xs_km15 column that is only
                # created later inside the Stage-2 CFF-derivative script.
                sig = math.hypot(
                    float(row["xs_stat_pseudo_abs"]),
                    float(row["xs_ptp_sys_pseudo_abs"])
                )
            else:
                sig = math.hypot(
                    float(row["bsa_stat_pseudo_abs"]),
                    float(row["bsa_ptp_sys_pseudo_abs"])
                )
            if np.isfinite(sig) and sig > 0:
                obs.append(dict(
                    target="p", cell=cell, kind=kind, y0=y0, sigma=sig,
                    d_bu=du, d_bd=dd, d_ReH=dhre, d_ImH=dhim
                ))

    # Neutron XS: use the measured relative statistical precision, scaled with L.
    # The unresolved quoted total systematic is not silently made diagonal here.
    for r in rgb_xs.itertuples(index=False):
        cell=f"nxs:{int(r.kin_bin)}"
        y0,du,dd,dhre,dhim=derivatives(r.xB,r.Q2_GeV2,r.t_abs_GeV2,r.phi_deg,10.45,"n","xs")
        sig=abs(y0)*abs(r.rel_stat)/math.sqrt(factor)
        if np.isfinite(sig) and sig>0:
            obs.append(dict(target="n",cell=cell,kind="xs",y0=y0,sigma=sig,
                            d_bu=du,d_bd=dd,d_ReH=dhre,d_ImH=dhim))

    # Neutron BSA: use ONE published projection only (t), never all 63 points.
    nb=rgb_bsa[rgb_bsa["projection"]=="t"]
    for r in nb.itertuples(index=False):
        cell=f"nbsa:{int(r.projected_bin)}"
        y0,du,dd,dhre,dhim=derivatives(r.xB,r.Q2_GeV2,r.t_abs_GeV2,r.phi_deg,10.45,"n","bsa")
        sig=abs(r.stat)/math.sqrt(factor)
        if np.isfinite(sig) and sig>0:
            obs.append(dict(target="n",cell=cell,kind="bsa",y0=y0,sigma=sig,
                            d_bu=du,d_bd=dd,d_ReH=dhre,d_ImH=dhim))
    return pd.DataFrame(obs)



def _assign_rgb_rga_ptp_fraction(rgb_xs, proton_1x, obs_1x, scale=1.25):
    """Assign each RGB XS point a conservative RGA-like relative PTP systematic.

    The RGA relative PTP pattern is taken from the Stage-2 proton XS input and
    matched to RGB by nearest (xB,Q2,|t|,phi) kinematics.  The assigned
    fractional uncertainty is then multiplied by `scale` (default 1.25).
    """
    pobs = obs_1x[(obs_1x.target == "p") & (obs_1x.kind == "xs")].reset_index(drop=True)
    pdat = proton_1x.reset_index(drop=True)
    if len(pobs) != len(pdat):
        raise RuntimeError(f"RGA PTP mapping mismatch: {len(pobs)} XS derivative rows vs {len(pdat)} Stage-2 rows")
    rel = np.abs(pdat["xs_ptp_sys_pseudo_abs"].to_numpy(float)) / np.maximum(np.abs(pobs["y0"].to_numpy(float)), 1e-12)
    coords = pdat[["xB","Q2","t_abs","phi_deg"]].to_numpy(float)
    # Dimension scales only define the nearest-neighbour metric; they do not
    # alter the RGA systematic values themselves.
    scales = np.array([0.08, 0.8, 0.15, 30.0], dtype=float)
    out=[]
    for r in rgb_xs.itertuples(index=False):
        q=np.array([r.xB,r.Q2_GeV2,r.t_abs_GeV2,r.phi_deg],dtype=float)
        d=coords-q[None,:]
        # phi is periodic.
        d[:,3]=((d[:,3]+180.0)%360.0)-180.0
        j=int(np.argmin(np.sum((d/scales[None,:])**2,axis=1)))
        out.append(scale*rel[j])
    return np.asarray(out,dtype=float)


def _set_observable_uncertainties(obs_template, factor, proton_factor, rgb_xs, rgb_bsa, rgb_rel_ptp, include_ptp=True):
    """Reuse luminosity-independent derivatives and update only uncertainties.

    include_ptp=True is the conference baseline: RGA uses its measured PTP
    systematics and RGB XS receives the matched RGA PTP pattern scaled by 1.25.
    include_ptp=False is the strict statistics-only backup: ALL RGA XS/BSA PTP
    terms and the assigned RGB XS PTP term are removed.
    """
    obs=obs_template.copy()
    # Proton rows: exact Stage-2 uncertainty prescription at this luminosity.
    p_xs_idx=obs.index[(obs.target=="p")&(obs.kind=="xs")].to_numpy()
    p_bsa_idx=obs.index[(obs.target=="p")&(obs.kind=="bsa")].to_numpy()
    if len(p_xs_idx)!=len(proton_factor) or len(p_bsa_idx)!=len(proton_factor):
        raise RuntimeError("Proton derivative/template row count does not match Stage-2 input")
    if include_ptp:
        obs.loc[p_xs_idx,"sigma"]=[math.hypot(float(r["xs_stat_pseudo_abs"]),float(r["xs_ptp_sys_pseudo_abs"])) for r in proton_factor.to_dict("records")]
        obs.loc[p_bsa_idx,"sigma"]=[math.hypot(float(r["bsa_stat_pseudo_abs"]),float(r["bsa_ptp_sys_pseudo_abs"])) for r in proton_factor.to_dict("records")]
    else:
        obs.loc[p_xs_idx,"sigma"]=[abs(float(r["xs_stat_pseudo_abs"])) for r in proton_factor.to_dict("records")]
        obs.loc[p_bsa_idx,"sigma"]=[abs(float(r["bsa_stat_pseudo_abs"])) for r in proton_factor.to_dict("records")]

    # RGB XS: statistical error scales as 1/sqrt(L); the assigned RGA-like PTP
    # systematic is fixed with luminosity and is 1.25x the matched RGA fraction.
    n_xs_idx=obs.index[(obs.target=="n")&(obs.kind=="xs")].to_numpy()
    if len(n_xs_idx)!=len(rgb_xs):
        raise RuntimeError("RGB XS derivative/template row count mismatch")
    y=np.abs(obs.loc[n_xs_idx,"y0"].to_numpy(float))
    relstat=rgb_xs["rel_stat"].to_numpy(float)/math.sqrt(factor)
    obs.loc[n_xs_idx,"sigma"]=y*(np.hypot(relstat,rgb_rel_ptp) if include_ptp else relstat)

    # RGB BSA: retain the published statistical projection used previously.
    nb=rgb_bsa[rgb_bsa["projection"]=="t"].reset_index(drop=True)
    n_bsa_idx=obs.index[(obs.target=="n")&(obs.kind=="bsa")].to_numpy()
    if len(n_bsa_idx)!=len(nb):
        raise RuntimeError("RGB BSA derivative/template row count mismatch")
    obs.loc[n_bsa_idx,"sigma"]=np.abs(nb["stat"].to_numpy(float))/math.sqrt(factor)
    return obs

def _he_global_covariance(obs, include_neutron=False, neutron_kinds=None):
    """Profile local ReH/ImH while selecting the requested target/observable set."""
    if not include_neutron:
        use = obs[obs.target == "p"].copy()
    else:
        use = obs.copy()
        if neutron_kinds is not None:
            neutron_kinds = set(neutron_kinds)
            use = use[(use.target == "p") |
                      ((use.target == "n") & use.kind.isin(neutron_kinds))].copy()

    cells = sorted(use.cell.unique())
    layout = ["beta_u", "beta_d"] + [(c, "ReH") for c in cells] + [(c, "ImH") for c in cells]
    idx = {k:i for i,k in enumerate(layout)}
    A=[]; sig=[]
    for r in use.itertuples(index=False):
        v=np.zeros(len(layout))
        v[idx["beta_u"]]=r.d_bu
        v[idx["beta_d"]]=r.d_bd
        v[idx[(r.cell,"ReH")]]=r.d_ReH
        v[idx[(r.cell,"ImH")]]=r.d_ImH
        A.append(v); sig.append(r.sigma)
    B=np.asarray(A)/np.asarray(sig)[:,None]
    info=B.T@B
    cov=np.linalg.pinv(info,rcond=1e-11)
    rank=int(np.linalg.matrix_rank(B))
    return cov,layout,rank,len(layout),len(B)


def _b20_covariance_from_beta(cov,layout):
    idx={k:i for i,k in enumerate(layout)}
    cb=cov[np.ix_([idx["beta_u"],idx["beta_d"]],[idx["beta_u"],idx["beta_d"]])]
    h=1.0e-4
    du=(b20_from_beta(KAPPA_U,DFJK_ALPHA_E,WORKSHOP_BETA_U_E+h)-b20_from_beta(KAPPA_U,DFJK_ALPHA_E,WORKSHOP_BETA_U_E-h))/(2*h)
    dd=(b20_from_beta(KAPPA_D,DFJK_ALPHA_E,WORKSHOP_BETA_D_E+h)-b20_from_beta(KAPPA_D,DFJK_ALPHA_E,WORKSHOP_BETA_D_E-h))/(2*h)
    J=np.diag([du,dd])
    return J@cb@J.T,cb


def _ellipse_points(mean,cov,delta_chi2=2.30,n=300):
    vals,vecs=np.linalg.eigh(cov); vals=np.maximum(vals,0.0)
    ang=np.linspace(0,2*np.pi,n)
    circ=np.vstack([np.cos(ang),np.sin(ang)])
    pts=np.asarray(mean)[:,None]+vecs@np.diag(np.sqrt(delta_chi2*vals))@circ
    return pts.T




def load_nnpdf_a20_summary(path: Path):
    """Read the one-time NNPDF4.0 valence A20 extraction.

    Returns the mean vector (u_v,d_v) and its 2x2 replica covariance.
    The Stage-5 workshop projection treats this external PDF covariance as
    independent of the projected DVCS B20 covariance.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"NNPDF A20 summary not found: {path}\n"
            "Run extract_nnpdf40_valence_moments.py first."
        )
    df = pd.read_csv(path)
    if not {"quantity", "value"}.issubset(df.columns):
        raise RuntimeError(f"Unexpected A20 summary format in {path}")
    vals = dict(zip(df["quantity"].astype(str), df["value"]))
    required = ["Q2_GeV2", "A20_uv_mean", "A20_uv_std",
                "A20_dv_mean", "A20_dv_std", "cov_uv_dv"]
    missing = [k for k in required if k not in vals]
    if missing:
        raise RuntimeError(f"A20 summary is missing: {missing}")
    q2 = float(vals["Q2_GeV2"])
    if not np.isclose(q2, A20_Q2_GEV2, rtol=0.0, atol=1e-9):
        raise RuntimeError(
            f"A20 summary has Q^2={q2:g} GeV^2; Stage 5 expects "
            f"Q^2={A20_Q2_GEV2:g} GeV^2."
        )
    mean = np.array([float(vals["A20_uv_mean"]), float(vals["A20_dv_mean"])])
    su = float(vals["A20_uv_std"]); sd = float(vals["A20_dv_std"])
    covud = float(vals["cov_uv_dv"])
    cov = np.array([[su*su, covud], [covud, sd*sd]], dtype=float)
    return mean, cov, vals


def _corr_from_cov(cov):
    den = math.sqrt(max(float(cov[0,0]*cov[1,1]), 1e-300))
    return float(cov[0,1]/den)


def make_ji_projection(contour_data, a20_mean, a20_cov, figdir, tabdir, tag="baseline"):
    """Combine external A20 with projected B20 using the valence Ji sum rule.

      J_qv = 1/2 (A20_qv + B20_qv)
      C_J  = 1/4 (C_A20 + C_B20)

    C_A20 and C_B20 are taken independent for this workshop projection.
    """
    rows=[]; ji_data=[]
    for name, bmean, bcov, betacov in contour_data:
        jmean = 0.5*(a20_mean + bmean)
        jcov = 0.25*(a20_cov + bcov)
        ji_data.append((name, jmean, jcov))
        rows.append(dict(
            scenario=name,
            Q2_GeV2=A20_Q2_GEV2,
            A20_uv=float(a20_mean[0]), A20_dv=float(a20_mean[1]),
            B20_uv=float(bmean[0]), B20_dv=float(bmean[1]),
            J_uv=float(jmean[0]), J_dv=float(jmean[1]),
            sigma_A20_uv=math.sqrt(max(a20_cov[0,0],0.0)),
            sigma_A20_dv=math.sqrt(max(a20_cov[1,1],0.0)),
            sigma_B20_uv=math.sqrt(max(bcov[0,0],0.0)),
            sigma_B20_dv=math.sqrt(max(bcov[1,1],0.0)),
            sigma_J_uv=math.sqrt(max(jcov[0,0],0.0)),
            sigma_J_dv=math.sqrt(max(jcov[1,1],0.0)),
            corr_A20=_corr_from_cov(a20_cov),
            corr_B20=_corr_from_cov(bcov),
            corr_J=_corr_from_cov(jcov),
        ))
    out=pd.DataFrame(rows)
    out.to_csv(tabdir/("ji_valence_projection_summary.csv" if tag=="baseline" else f"ji_valence_projection_summary_{tag}.csv"), index=False)

    labels={
        "p1x_plus_n1x": r"$p+n$ 1x: 68% contour",
        "p5x_plus_n5x": r"$p+n$ 5x: 68% contour",
        "p10x_plus_n10x": r"$p+n$ 10x: 68% contour",
    }
    fig,ax=plt.subplots(figsize=(7.4,6.2))
    pts_all=[]
    for name,mean,cov in ji_data:
        pts=_ellipse_points(mean,cov); pts_all.append(pts)
        ax.plot(pts[:,0],pts[:,1],linewidth=2.6,label=labels[name])
    allpn=np.vstack(pts_all)
    xmin,xmax=float(allpn[:,0].min()),float(allpn[:,0].max())
    ymin,ymax=float(allpn[:,1].min()),float(allpn[:,1].max())
    dx=max(xmax-xmin,0.04); dy=max(ymax-ymin,0.04)
    xmin-=0.25*dx; xmax+=0.25*dx; ymin-=0.25*dy; ymax+=0.25*dy
    ref=ji_data[0][1]
    ax.scatter([ref[0]],[ref[1]],marker="*",s=120,label="reference model",zorder=5)
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin,ymax)
    ax.set_xlabel(r"$J_{u_v}$"); ax.set_ylabel(r"$J_{d_v}$")
    ax.set_title(r"Projected valence Ji sum-rule flavor separation at $Q^2=3\,\mathrm{GeV}^2$")
    ax.grid(alpha=.22); ax.legend(fontsize=8)
    fig.tight_layout()
    fname = "ji_uv_dv_projected_contours.png" if tag=="baseline" else f"ji_uv_dv_projected_contours_{tag}.png"
    fig.savefig(figdir/fname,dpi=200); plt.close(fig)
    return out


def make_ji_world_comparison(ji_summary, figdir, tabdir):
    """World-context plot for the Stage-5 valence Ji projection.

    The CLAS12 contours and the two modern elastic/lattice points are valence
    quantities.  The HERMES band is *not*: it is the historical full-flavor
    proton-DVCS model constraint Ju + Jd/2.8 = 0.49 +/- 0.17.  It is included
    only to show the scale/direction of the early direct-DVCS constraint.

    Literature numbers and URLs are documented in the module docstring.
    No covariance between J_uv and J_dv is published with the compact Cichy or
    Diehl--Kroll numbers used here, so they are shown as marginal x/y error bars,
    NOT reconstructed 2D confidence ellipses.
    """
    rows = [
        dict(source="HERMES 2008 proton DVCS (full flavor)", kind="full_flavor_band",
             Ju=np.nan, Jd=np.nan, Ju_lo=np.nan, Ju_hi=np.nan, Jd_lo=np.nan, Jd_hi=np.nan,
             note="Ju + Jd/2.8 = 0.49 +/- 0.17 exp; historical full-flavor VGG/DD constraint"),
        dict(source="Diehl-Kroll 2013 elastic FF (valence)", kind="valence_point",
             Ju=0.230, Jd=-0.004, Ju_lo=0.024, Ju_hi=0.009, Jd_lo=0.016, Jd_hi=0.010,
             note="EPJC 73, 2397 (2013), arXiv:1302.4604; mu=2 GeV"),
        dict(source="Cichy et al. 2024 elastic+lattice (valence)", kind="valence_point",
             Ju=0.195, Jd=0.0173, Ju_lo=0.010, Ju_hi=0.010, Jd_lo=0.0046, Jd_hi=0.0046,
             note="PRD 110, 114025 (2024)"),
    ]
    pd.DataFrame(rows).to_csv(tabdir/"ji_world_comparison_literature_inputs.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.0,6.4))

    # Historical HERMES direct-DVCS band.  Axes are Ju horizontal, Jd vertical.
    # Central relation: Jd = 2.8*(0.49-Ju); +/-0.17 widens the band in the
    # measured linear combination.  This is full flavor, not valence.
    xu = np.linspace(0.05, 0.36, 400)
    yd_c = 2.8*(0.49-xu)
    yd_lo = 2.8*((0.49-0.17)-xu)
    yd_hi = 2.8*((0.49+0.17)-xu)
    ax.fill_between(xu, yd_lo, yd_hi, alpha=0.13,
                    label=r"HERMES DVCS: $J_u+J_d/2.8=0.49\pm0.17$ (full flavor)")
    ax.plot(xu, yd_c, linewidth=1.4, alpha=0.65)

    # Current Stage-5 projected valence contours.
    labels = {1:r"CLAS12 $p+n$ 1x (valence)", 5:r"CLAS12 $p+n$ 5x (valence)",
              10:r"CLAS12 $p+n$ 10x (valence)"}
    for L in (1,5,10):
        r = ji_summary[ji_summary["scenario"]==f"p{L}x_plus_n{L}x"].iloc[0]
        mean=np.array([r.J_uv,r.J_dv],float)
        su,sd=float(r.sigma_J_uv),float(r.sigma_J_dv); rho=float(r.corr_J)
        cov=np.array([[su*su,rho*su*sd],[rho*su*sd,sd*sd]])
        pts=_ellipse_points(mean,cov)
        ax.plot(pts[:,0],pts[:,1],linewidth=2.4,label=labels[L])

    # Directly comparable valence determinations.  These are marginal error bars
    # because a 2D covariance is not supplied by the compact published numbers.
    ax.errorbar(0.230,-0.004,xerr=np.array([[0.024],[0.009]]),
                yerr=np.array([[0.016],[0.010]]),fmt="s",capsize=3,markersize=6,
                label="Diehl-Kroll 2013 elastic FF (valence)")
    ax.errorbar(0.195,0.0173,xerr=0.010,yerr=0.0046,fmt="D",capsize=3,markersize=6,
                label="Cichy et al. 2024 elastic+lattice (valence)")

    ax.axhline(0,linewidth=0.8,alpha=0.25); ax.axvline(0,linewidth=0.8,alpha=0.25)
    ax.set_xlim(0.05,0.36); ax.set_ylim(-0.22,0.42)
    ax.set_xlabel(r"$J_u$ or $J_{u_v}$")
    ax.set_ylabel(r"$J_d$ or $J_{d_v}$")
    ax.set_title("Projected CLAS12 valence Ji constraints in world context")
    ax.text(0.02,0.02,"HERMES band is full flavor; all other entries shown are valence.\n"
            "Literature point error bars are marginal uncertainties, not 2D covariances.",
            transform=ax.transAxes,fontsize=8,va="bottom")
    ax.grid(alpha=.20); ax.legend(fontsize=7.4,loc="upper right")
    fig.tight_layout(); fig.savefig(figdir/"ji_world_comparison.png",dpi=200); plt.close(fig)



def make_ji_world_uncertainty_comparison(ji_summary, ji_summary_stat, figdir, tabdir):
    """Uncertainty-only world context, separated by what is actually constrained.

    The upper group contains direct/projected DVCS constraints on ONE full-flavor
    linear combination of Ju and Jd.  Those numbers are not sigma(J_uv) or
    sigma(J_dv), so they are deliberately not mixed with the lower group.

    The lower group contains genuinely flavor-separated VALENCE marginal
    uncertainties sigma(J_uv), sigma(J_dv).  This is where the CLAS12 RGA+RGB
    projection can be compared directly with Diehl--Kroll and Cichy et al.

    EIC/EicC use the stat+syst HERMES-region projections from Eqs. (49--50) of
    EPJC 83 (2023), not their smaller statistics-only numbers.  Hashamipour et
    al. PRD 105, 054002 (2022) is documented in the module header but omitted
    numerically because its J_uv/J_dv uncertainties are figure-only in the
    publication; no undocumented digitization is used here.
    """
    rows=[]
    # Direct DVCS one-combination benchmarks.  coefficient means
    # combination = first_flavor + second_flavor/coefficient.
    direct=[
        dict(label="HERMES p-DVCS [Airapetian et al. 2008]", status="measured", sigma_combination=0.21,
             combination="Ju + Jd/2.9", flavor_scope="full flavor",
             note="HERMES result quoted in Xie et al., EPJC 83, 900 (2023), Eq. 51; uncertainty used here matches the EIC/EicC comparison"),
        dict(label="Hall-A n-DVCS [Mazouz et al. 2007]", status="measured", sigma_combination=0.14,
             combination="Jd + Ju/5.0", flavor_scope="full flavor",
             note="Mazouz et al. PRL 99, 242501 (2007); VGG interpretation"),
        dict(label="EicC projection [Xie et al. 2023]", status="projection stat+syst", sigma_combination=0.08,
             combination="Ju + Jd/2.9", flavor_scope="full flavor",
             note="EPJC 83 (2023), Eq. 49, HERMES-like region"),
        dict(label="EIC projection [Xie et al. 2023]", status="projection stat+syst", sigma_combination=0.06,
             combination="Ju + Jd/3.0", flavor_scope="full flavor",
             note="EPJC 83 (2023), Eq. 50, HERMES-like region"),
    ]
    rows.extend(direct)

    for L in (1,5,10):
        r=ji_summary[ji_summary["scenario"]==f"p{L}x_plus_n{L}x"].iloc[0]
        rows.append(dict(label=f"CLAS12 p+n {L}x", status="projected valence baseline",
                         sigma_Juv=float(r.sigma_J_uv), sigma_Jdv=float(r.sigma_J_dv),
                         flavor_scope="valence",
                         note="RGA+RGB only; baseline PTP systematics; RGH AUT not included"))
    r=ji_summary_stat[ji_summary_stat["scenario"]=="p10x_plus_n10x"].iloc[0]
    rows.append(dict(label="CLAS12 p+n 10x stat-only", status="projected valence stat-only",
                     sigma_Juv=float(r.sigma_J_uv), sigma_Jdv=float(r.sigma_J_dv),
                     flavor_scope="valence",
                     note="all RGA/RGB experimental PTP removed; RGH AUT not included"))
    rows += [
        dict(label="Diehl-Kroll 2013 [EPJC 73, 2397]", status="published valence asymmetric",
             sigma_Juv_low=0.024, sigma_Juv_high=0.009,
             sigma_Jdv_low=0.016, sigma_Jdv_high=0.010,
             flavor_scope="valence", note="elastic FF + GPD model; mu=2 GeV"),
        dict(label="Cichy et al. 2024 [PRD 110, 114025]", status="published valence",
             sigma_Juv=0.010, sigma_Jdv=0.0046, flavor_scope="valence",
             note="elastic FF + lattice; mu=2 GeV"),
        dict(label="Hashamipour et al. 2022", status="documented_not_plotted",
             flavor_scope="valence",
             note="Juv/Jdv shown with uncertainties in Fig. 11 but compact numerical errors not tabulated; no digitization used"),
    ]
    pd.DataFrame(rows).to_csv(tabdir/"ji_world_uncertainty_comparison_inputs.csv", index=False)

    # One figure, two panels: do not make full-flavor combination widths look
    # like valence marginal errors simply because all numbers are dimensionless.
    fig,(ax0,ax1)=plt.subplots(2,1,figsize=(8.6,7.4),gridspec_kw={"height_ratios":[0.85,1.35]})

    # --- direct DVCS combination constraints ---
    labels0=[d["label"] for d in direct]
    y0=np.arange(len(labels0))[::-1]
    for yy,d in zip(y0,direct):
        marker="o"
        projected = d["status"] != "measured"
        ax0.scatter(d["sigma_combination"], yy, marker=marker, s=65, zorder=4,
                    facecolors="none" if projected else "black", edgecolors="black", linewidths=1.5)
        ax0.text(d["sigma_combination"]+0.004,yy,d["combination"],va="center",fontsize=8)
    ax0.set_yticks(y0); ax0.set_yticklabels(labels0)
    ax0.set_xlim(0,0.23)
    ax0.set_xlabel("Uncertainty on quoted $J_u$--$J_d$ linear combination")
    ax0.set_title("Direct DVCS constraints and collider projections (full flavor)",fontsize=11)
    ax0.grid(axis="x",alpha=.22)

    # --- flavor-separated valence marginal uncertainties ---
    labels1=["CLAS12 p+n 1x","CLAS12 p+n 5x","CLAS12 p+n 10x",
             "CLAS12 p+n 10x stat-only","Diehl-Kroll 2013 [EPJC 73, 2397]","Cichy et al. 2024 [PRD 110, 114025]"]
    y1=np.arange(len(labels1))[::-1]
    for i,L in enumerate((1,5,10)):
        r=ji_summary[ji_summary["scenario"]==f"p{L}x_plus_n{L}x"].iloc[0]
        ax1.scatter(float(r.sigma_J_uv),y1[i]+0.10,marker="o",s=55,zorder=4,color="black",
                    label=r"$\sigma(J_{u_v})$" if i==0 else None)
        ax1.scatter(float(r.sigma_J_dv),y1[i]-0.10,marker="o",s=55,zorder=4,facecolors="none",edgecolors="black",linewidths=1.5,
                    label=r"$\sigma(J_{d_v})$" if i==0 else None)
    r=ji_summary_stat[ji_summary_stat["scenario"]=="p10x_plus_n10x"].iloc[0]
    ax1.scatter(float(r.sigma_J_uv),y1[3]+0.10,marker="o",s=55,zorder=4,color="black")
    ax1.scatter(float(r.sigma_J_dv),y1[3]-0.10,marker="o",s=55,zorder=4,facecolors="none",edgecolors="black",linewidths=1.5)

    # DK13 quotes asymmetric errors on J itself.  For this uncertainty-scale
    # comparison use the arithmetic mean of the magnitudes of the two sides:
    #   u_v: (0.024 + 0.009)/2 = 0.0165
    #   d_v: (0.016 + 0.010)/2 = 0.0130
    # This plotting convention is stated explicitly below the y-axis labels.
    dk_u = 0.5*(0.024 + 0.009)
    dk_d = 0.5*(0.016 + 0.010)
    ax1.scatter(dk_u,y1[4]+0.10,marker="o",s=55,zorder=4,color="black")
    ax1.scatter(dk_d,y1[4]-0.10,marker="o",s=55,zorder=4,facecolors="none",edgecolors="black",linewidths=1.5)
    ax1.scatter(0.010,y1[5]+0.10,marker="o",s=55,zorder=4,color="black")
    ax1.scatter(0.0046,y1[5]-0.10,marker="o",s=55,zorder=4,facecolors="none",edgecolors="black",linewidths=1.5)

    ax1.set_yticks(y1); ax1.set_yticklabels(labels1)
    ax1.set_xlim(0,0.085)
    ax1.set_xlabel("Marginal uncertainty magnitude")
    ax1.set_title(r"Flavor-separated valence $J_{u_v},J_{d_v}$ constraints",fontsize=11)
    ax1.grid(axis="x",alpha=.22); ax1.legend(fontsize=8,loc="lower right",frameon=False)
    ax1.text(-0.01,-0.17,"Diehl--Kroll: average of asymmetric errors",
             transform=ax1.transAxes,ha="left",va="top",fontsize=7.2,clip_on=False)
    ax1.text(0.99,-0.17,
             "CLAS12 = RGA+RGB only; projected RGH $A_{UT}$ sensitivity to $E$ is not included.",
             transform=ax1.transAxes,ha="right",va="top",fontsize=7.2,clip_on=False)

    fig.suptitle("Ji-sum-rule uncertainty scales in world context",fontsize=14)
    fig.tight_layout(rect=(0,0.065,1,.965))
    fig.savefig(figdir/"ji_world_uncertainty_comparison.png",dpi=200)
    plt.close(fig)

def run_stage5_he_projection(stage2_dir, rgb_xs, rgb_bsa, figdir, tabdir, a20_mean, a20_cov):
    """Run joint H+E projection and target/observable ablations.

    Proton-only local covariance contours are retained only as diagnostics in
    the tables; the conference figure does not draw them as literal 68% B20
    regions because the corresponding beta uncertainties leave the intended
    model domain.
    """
    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError("Stage-5 H+E projection requires Gepard/KM15 on ifarm.") from exc

    th = th_KM15
    scenarios = []
    contour_data = []
    obs_by_factor = {}
    obs_by_factor_stat_only = {}

    # The observable derivatives themselves do not depend on luminosity.  Build
    # them once (or load the existing 1x cache), then change only the uncertainty
    # model for 1x/5x/10x.  This avoids repeated slow Gepard derivative calls.
    cache_file = tabdir / "he_observable_derivatives_1x.csv"
    if cache_file.exists():
        print("[Stage5 H+E] loading cached observable derivatives ...")
        obs_template = pd.read_csv(cache_file)
        print(f"[Stage5 H+E] using cached derivatives: {cache_file}")
    else:
        print("[Stage5 H+E] building observable derivatives once ...")
        p1 = _load_stage2_proton_inputs(stage2_dir, 1)
        obs_template = _build_he_observable_rows(th, g, p1, rgb_xs, rgb_bsa, 1)
        obs_template.to_csv(cache_file, index=False)
        print(f"[Stage5 H+E] cached observable derivatives: {cache_file}")

    p1 = _load_stage2_proton_inputs(stage2_dir, 1)
    rgb_rel_ptp = _assign_rgb_rga_ptp_fraction(rgb_xs, p1, obs_template, scale=1.25)
    pd.DataFrame({
        "kin_bin": rgb_xs["kin_bin"].to_numpy(),
        "phi_deg": rgb_xs["phi_deg"].to_numpy(),
        "xB": rgb_xs["xB"].to_numpy(),
        "Q2_GeV2": rgb_xs["Q2_GeV2"].to_numpy(),
        "t_abs_GeV2": rgb_xs["t_abs_GeV2"].to_numpy(),
        "assigned_RGA_ptp_frac_times_1p25": rgb_rel_ptp,
    }).to_csv(tabdir/"rgb_assigned_rga_ptp_systematics.csv",index=False)
    print(f"[Stage5 H+E] RGB XS assigned RGA-like PTP: median={100*np.median(rgb_rel_ptp):.2f}% (includes 1.25x scale)")

    for factor in (1, 5, 10):
        pf = _load_stage2_proton_inputs(stage2_dir, factor)
        obs_by_factor[factor] = _set_observable_uncertainties(
            obs_template, factor, pf, rgb_xs, rgb_bsa, rgb_rel_ptp, include_ptp=True
        )
        obs_by_factor_stat_only[factor] = _set_observable_uncertainties(
            obs_template, factor, pf, rgb_xs, rgb_bsa, rgb_rel_ptp, include_ptp=False
        )

    # Four conference-facing scenarios.  p10x is deliberately included so the
    # plot separates "more proton statistics" from "new neutron flavor info".
    configs = [
        ("p1x_plus_n1x",       1,  True,  ("xs","bsa")),
        ("p5x_plus_n5x",       5,  True,  ("xs","bsa")),
        ("p10x_plus_n10x",    10,  True,  ("xs","bsa")),
    ]

    mean_beta = np.array([WORKSHOP_BETA_U_E, WORKSHOP_BETA_D_E], dtype=float)
    mean_b20 = np.array([
        b20_from_beta(KAPPA_U, DFJK_ALPHA_E, WORKSHOP_BETA_U_E),
        b20_from_beta(KAPPA_D, DFJK_ALPHA_E, WORKSHOP_BETA_D_E),
    ])

    for name, factor, incn, nkinds in configs:
        cov, layout, rank, npar, nrow = _he_global_covariance(
            obs_by_factor[factor], incn, nkinds
        )
        bcov, betacov = _b20_covariance_from_beta(cov, layout)
        contour_data.append((name, mean_b20.copy(), bcov, betacov))


        scenarios.append(dict(
            scenario=name, luminosity_factor=factor, include_neutron=incn,
            neutron_observables="none" if not incn else "+".join(nkinds),
            n_rows=nrow, n_parameters=npar, rank=rank, rank_deficit=npar-rank,
            sigma_beta_u=math.sqrt(max(betacov[0,0],0)),
            sigma_beta_d=math.sqrt(max(betacov[1,1],0)),
            sigma_B20_uv=math.sqrt(max(bcov[0,0],0)),
            sigma_B20_dv=math.sqrt(max(bcov[1,1],0)),
            corr_B20=float(bcov[0,1]/math.sqrt(max(bcov[0,0]*bcov[1,1],1e-300)))
        ))

    # Strict statistics-only backup: remove PTP terms from BOTH RGA and RGB.
    stat_contour_data=[]
    stat_rows=[]
    for name, factor, incn, nkinds in configs:
        cov, layout, rank, npar, nrow = _he_global_covariance(
            obs_by_factor_stat_only[factor], incn, nkinds
        )
        bcov, betacov = _b20_covariance_from_beta(cov, layout)
        stat_contour_data.append((name, mean_b20.copy(), bcov, betacov))
        stat_rows.append(dict(
            scenario=name, luminosity_factor=factor, n_rows=nrow, n_parameters=npar, rank=rank,
            sigma_B20_uv=math.sqrt(max(bcov[0,0],0)),
            sigma_B20_dv=math.sqrt(max(bcov[1,1],0)),
            corr_B20=_corr_from_cov(bcov)))
    pd.DataFrame(stat_rows).to_csv(tabdir/"b20_he_projection_summary_stat_only.csv",index=False)

    summary = pd.DataFrame(scenarios)
    summary.to_csv(tabdir/"b20_he_projection_summary.csv", index=False)

    # Neutron ablation: diagnostic only, not extra conference contours.
    ablations = []
    for factor in (1, 5, 10):
        for label, kinds in [
            ("proton_only", None),
            ("plus_neutron_XS_only", ("xs",)),
            ("plus_neutron_BSA_only", ("bsa",)),
            ("plus_neutron_XS_and_BSA", ("xs","bsa")),
        ]:
            incn = kinds is not None
            cov, layout, rank, npar, nrow = _he_global_covariance(
                obs_by_factor[factor], incn, kinds
            )
            bcov, betacov = _b20_covariance_from_beta(cov, layout)
            ablations.append(dict(
                luminosity_factor=factor, configuration=label,
                n_rows=nrow, n_parameters=npar, rank=rank,
                sigma_B20_uv=math.sqrt(max(bcov[0,0],0)),
                sigma_B20_dv=math.sqrt(max(bcov[1,1],0)),
                corr_B20=float(bcov[0,1]/math.sqrt(max(bcov[0,0]*bcov[1,1],1e-300)))
            ))
    ablation_df = pd.DataFrame(ablations)
    ablation_df.to_csv(tabdir/"b20_neutron_observable_ablation.csv", index=False)

    # Conference-facing B20 figure: only the quantitative p+n contours.
    fig, ax = plt.subplots(figsize=(7.4,6.2))
    labels = {
        "p1x_plus_n1x": r"$p+n$ 1x: 68% contour",
        "p5x_plus_n5x": r"$p+n$ 5x: 68% contour",
        "p10x_plus_n10x": r"$p+n$ 10x: 68% contour",
    }
    allpts=[]
    for name, mean, cov, betacov in contour_data:
        pts=_ellipse_points(mean,cov)
        allpts.append(pts)
        ax.plot(pts[:,0],pts[:,1],linewidth=2.6,label=labels[name])
    allpn=np.vstack(allpts)
    xmin,xmax=float(allpn[:,0].min()),float(allpn[:,0].max())
    ymin,ymax=float(allpn[:,1].min()),float(allpn[:,1].max())
    dx=max(xmax-xmin,0.08); dy=max(ymax-ymin,0.08)
    xmin-=0.25*dx; xmax+=0.25*dx; ymin-=0.25*dy; ymax+=0.25*dy
    ax.scatter([mean_b20[0]],[mean_b20[1]],marker="*",s=120,label="reference model",zorder=5)
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin,ymax)
    ax.set_xlabel(r"$B_{20}^{u_v}(0)$"); ax.set_ylabel(r"$B_{20}^{d_v}(0)$")
    ax.set_title(r"Projected valence flavor separation allowing $\mathcal{H}$ to vary")
    ax.grid(alpha=.22); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(figdir/"b20_uv_dv_he_marginalized_contours.png",dpi=200); plt.close(fig)

    # Statistics-only B20 backup figure.
    fig, ax = plt.subplots(figsize=(7.4,6.2))
    allpts=[]
    for name, mean, cov, betacov in stat_contour_data:
        pts=_ellipse_points(mean,cov); allpts.append(pts)
        ax.plot(pts[:,0],pts[:,1],linewidth=2.6,label=labels[name])
    allpn=np.vstack(allpts); xmin,xmax=allpn[:,0].min(),allpn[:,0].max(); ymin,ymax=allpn[:,1].min(),allpn[:,1].max()
    dx=max(xmax-xmin,0.08); dy=max(ymax-ymin,0.08)
    ax.set_xlim(xmin-.25*dx,xmax+.25*dx); ax.set_ylim(ymin-.25*dy,ymax+.25*dy)
    ax.scatter([mean_b20[0]],[mean_b20[1]],marker="*",s=120,label="reference model",zorder=5)
    ax.set_xlabel(r"$B_{20}^{u_v}(0)$"); ax.set_ylabel(r"$B_{20}^{d_v}(0)$")
    ax.set_title(r"Projected valence flavor separation — statistics only")
    ax.grid(alpha=.22); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(figdir/"b20_uv_dv_he_marginalized_contours_stat_only.png",dpi=200); plt.close(fig)

    # Remove stale v26 nonlinear-validation products so the output directory
    # reflects the current production analysis after an in-place rerun.
    for stale in (figdir/"b20_exact_beta_mapping_validation.png",
                  tabdir/"b20_exact_beta_mapping_validation.csv"):
        if stale.exists():
            stale.unlink()

    ji_summary = make_ji_projection(contour_data, a20_mean, a20_cov, figdir, tabdir, tag="baseline")
    ji_summary_stat = make_ji_projection(stat_contour_data, a20_mean, a20_cov, figdir, tabdir, tag="stat_only")
    make_ji_world_comparison(ji_summary, figdir, tabdir)
    make_ji_world_uncertainty_comparison(ji_summary, ji_summary_stat, figdir, tabdir)
    print("\n[Stage5 H+E] strict statistics-only Ji projection (all RGA/RGB PTP removed):")
    print(ji_summary_stat.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\n[Stage5 H+E] neutron observable ablation:")
    print(ablation_df.to_string(index=False))

    return summary, ji_summary

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
    ap.add_argument("--a20-summary", type=Path, default=None,
                    help="NNPDF4.0 A20 summary; default output/stage5_ji/tables/nnpdf40_valence_A20_summary_Q2_3.csv next to this script.")
    ap.add_argument("--stage2-dir", type=Path, default=None,
                    help="Stage-2 pseudo-data tables; default output/stage2/tables next to this script.")
    args = ap.parse_args()
    here = Path(__file__).resolve().parent
    if args.a20_summary is None:
        args.a20_summary = (here / A20_SUMMARY_DEFAULT).resolve()
    else:
        args.a20_summary = args.a20_summary.expanduser().resolve()
    if args.stage2_dir is None:
        args.stage2_dir = (here / "output" / "stage2" / "tables").resolve()
    else:
        args.stage2_dir = args.stage2_dir.expanduser().resolve()
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
    make_step2a_t_dependence(figdir, tabdir)
    finite_xi_validation = validate_finite_xi_model(tabdir)
    conference_scenarios = write_stage5_projection_plan(tabdir)

    bsa = load_actual_rgb_bsa_directory(args.rgb_bsa_dir)
    bsa.to_csv(tabdir / "rgb_published_bsa_actual_points.csv", index=False)

    a20_mean, a20_cov, a20_meta = load_nnpdf_a20_summary(args.a20_summary)
    he_summary, ji_summary = run_stage5_he_projection(
        args.stage2_dir, xs, bsa, figdir, tabdir, a20_mean, a20_cov
    )

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
    print("STAGE 5 v29 — VALENCE p+n H+E → Ji PROJECTION")
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
    print("Stage-5 production GPD model:")
    print("  valence E_v^u and E_v^d with beta_u,beta_d shape parameters")
    print("  anomalous magnetic moments fix the zeroth moments")
    print(f"  fixed alpha_E={DFJK_ALPHA_E:.2f}, B_E={WORKSHOP_E_T_SLOPE:.1f} GeV^-2, DD profile b={WORKSHOP_DD_PROFILE_B:.1f}")
    print("  finite skewness is now constructed internally with the full DD integral")
    max_forward_dev = np.max(
        np.abs(finite_xi_validation["xi0_over_expected"] - 1.0)
    )
    max_finite_dev = np.max(
        np.abs(finite_xi_validation["finite_xi_over_expected"] - 1.0)
    )
    print(f"  forward zeroth-moment validation:  max fractional deviation = {max_forward_dev:.4g}")
    print(f"  finite-xi zeroth-moment validation: max fractional deviation = {max_finite_dev:.4g}")
    print("  validation includes full finite-xi support -xi <= x <= 1")
    print("  pedagogical alpha/beta/skewness figures are no longer produced")
    print()
    print("Conference target:")
    print("  one B20_uv vs B20_dv figure separating luminosity from target complementarity:")
    print("    p-only: local flavor-degeneracy directions (not literal closed 68% regions)")
    print("    p+n: quantitative 68% contours at 1x and 10x")
    print("  IMPORTANT: no contour is fabricated from raw BSA errors.")
    print("  It must come from the observable-level DVCS response to beta_u,beta_d.")
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
    print("JOINT H+E B20 PROJECTION:")
    print("  beta_u,beta_d are global E_v parameters of interest")
    print("  ReH,ImH are free local nuisance parameters in every p/n kinematic cell")
    print("  therefore H is constrained by the same RGA/RGB observables, not held fixed")
    print("  neutron fit uses preliminary XS plus ONE published BSA projection (t)")
    print("  neutron XS uses an RGA point-to-point systematic pattern, matched in kinematics and scaled by 1.25x")
    print("  p-only closed B20 ellipses are not plotted: their local beta covariance leaves the intended model domain")
    print("  p-only dashed lines show only the locally constrained flavor direction")
    print("  p+n contours are the quantitative 68% regions used for the workshop projection")
    print(he_summary.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print()
    print("VALENCE Ji SUM-RULE PROJECTION:")
    print(f"  NNPDF4.0 A20 input: {args.a20_summary}")
    print(f"  fixed quoted scale: Q^2={A20_Q2_GEV2:g} GeV^2")
    print("  J_qv = 0.5 * (A20_qv + B20_qv)")
    print("  external PDF A20 covariance and projected DVCS B20 covariance are treated as independent")
    print("  B20 model is interpreted at the same fixed workshop scale; no QCD evolution across CLAS12 bins is implemented")
    print(ji_summary.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print()
    print()
    print("OUTPUT LAYOUT:")
    print(f"  figures -> {figdir}")
    print(f"  tables  -> {tabdir}")
    print(f"Wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
