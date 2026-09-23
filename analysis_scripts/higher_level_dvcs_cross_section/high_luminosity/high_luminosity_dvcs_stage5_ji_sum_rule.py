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
        integral = np.trapz(profile, alpha)

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
    make_step2a_t_dependence(figdir, tabdir)
    make_step2b_skewness_profile(figdir, tabdir)

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
    print("STAGE 5 v18 — STEP 2B SKEWNESS PROFILE + RGB/RGA PROJECTIONS")
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
    print("Stage-5 Step 1 — minimal DFJK E_v forward-limit model:")
    print("  E_v^q(x, xi=0, t=0) = N_q x^(-alpha) (1-x)^(beta_q)")
    print(f"  fixed alpha = {DFJK_ALPHA_E:.2f}")
    print(f"  kappa_u = 2*kappa_p + kappa_n = {KAPPA_U:+.6f}")
    print(f"  kappa_d = kappa_p + 2*kappa_n = {KAPPA_D:+.6f}")
    print("  N_u and N_d are fixed so integral E_v^q dx = kappa_q")
    print(f"  workshop reference beta_u = {WORKSHOP_BETA_U_E:.1f}")
    print(f"  workshop reference beta_d = {WORKSHOP_BETA_D_E:.1f}")
    print(f"  beta_u scan = {BETA_U_E_SCAN[0]:.1f} -- {BETA_U_E_SCAN[1]:.1f}")
    print(f"  beta_d scan = {BETA_D_E_SCAN[0]:.1f} -- {BETA_D_E_SCAN[1]:.1f}")
    print("  Step 2A adds ONLY nonzero t at xi=0:")
    print("    E_v^q(x,0,t) = E_v^q(x,0,0) exp(B_E t)")
    print(f"    fixed illustrative B_E = {WORKSHOP_E_T_SLOPE:.1f} GeV^-2")
    print("  Step 2B introduces ONLY the normalized double-distribution profile:")
    print("    beta = forward parton momentum fraction")
    print("    alpha = how longitudinal momentum transfer is shared")
    print(f"    fixed profile parameter b = {WORKSHOP_DD_PROFILE_B:.1f}")
    print("    integral d alpha pi_b(beta,alpha) = 1")
    print("  NOT included yet: finite-xi E(x,xi,t), CFFs, or a DVCS fit.")
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
