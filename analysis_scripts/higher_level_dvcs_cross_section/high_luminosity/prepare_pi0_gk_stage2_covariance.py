#!/usr/bin/env python3
"""
Stage 2: blinded pi0 high-luminosity projection preparation.

Purpose
-------
Build the experimental uncertainty/covariance model and the common-kinematics
Rosenbluth model-query grid. This script deliberately does NOT use preliminary
CLAS12 reduced-cross-section central values.

Inputs
------
  output/pi0_gk_stage1/{01_blinded_rga_bins.csv,02_blinded_rgk_bins.csv,
                        06_exact_common_cells.csv}
  import/fa18_rosenbluth_inputs_20260924T165948Z.tar.gz

Optional model input
--------------------
  --gk-results CSV containing model structure functions at the Stage-2 common
  points. Model-result parsing/pseudo-data construction is intentionally deferred
  until the actual PARTONS output format is known.

Outputs
-------
  output/pi0_gk_stage2/
      tables/
      model_queries/
      figures/
      summary.txt

Blinding
--------
No measured CLAS12 cross-section central values or absolute cross-section
uncertainties are read from the CSVs. The NPZ covariance is used only to obtain
within-cell correlation matrices; absolute covariance magnitudes are discarded.
The supplied relative_uncertainty column sets the per-bin fractional precision.
"""

from __future__ import annotations
import argparse, hashlib, math, tarfile, tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TAR_NAME = "fa18_rosenbluth_inputs_20260924T165948Z.tar.gz"

# Exposure factors relative to the supplied files.  These are configuration,
# not physics assumptions hidden in the code.
DEFAULT_RGA_CURRENT_FACTOR = 1.5
DEFAULT_RGK_FINAL_FACTOR = 8.0

FORBIDDEN = (
    "reduced_cross_section_nb_per_GeV2_rad",
    "propagated_statistical_and_finite_MC_uncertainty_nb_per_GeV2_rad",
)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def parse_args():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--stage1", type=Path, default=here/"output"/"pi0_gk_stage1")
    p.add_argument("--input", type=Path, default=here/"import"/TAR_NAME)
    p.add_argument("--output", type=Path, default=here/"output"/"pi0_gk_stage2")
    p.add_argument("--rga-current-factor", type=float, default=DEFAULT_RGA_CURRENT_FACTOR)
    p.add_argument("--rgk-final-factor", type=float, default=DEFAULT_RGK_FINAL_FACTOR)
    p.add_argument("--gk-results", type=Path, default=None,
                   help="Optional future PARTONS/GK result CSV; not required for covariance stage.")
    return p.parse_args()

def assert_blinded(df: pd.DataFrame, name: str):
    bad = [c for c in FORBIDDEN if c in df.columns]
    if bad:
        raise RuntimeError(f"{name} contains forbidden preliminary central/absolute columns: {bad}")

def locate_package(inp: Path):
    if inp.is_dir():
        return inp, None
    if not inp.exists():
        raise FileNotFoundError(inp)
    if not tarfile.is_tarfile(inp):
        raise RuntimeError(f"Not a tar archive: {inp}")
    tmp = tempfile.TemporaryDirectory(prefix="pi0_stage2_")
    root = Path(tmp.name)
    with tarfile.open(inp, "r:*") as tf:
        tf.extractall(root, filter="data")
    dirs = [p for p in root.iterdir() if p.is_dir()]
    if len(dirs) == 1:
        root = dirs[0]
    return root, tmp

def find_one(root: Path, pattern: str) -> Path:
    hits = list(root.rglob(pattern))
    if not hits:
        raise FileNotFoundError(f"Could not find {pattern} under {root}")
    # Prefer generic non-duplicated export name where possible.
    hits.sort(key=lambda p: (len(p.name), str(p)))
    return hits[0]

def covariance_sources(root: Path):
    rga = find_one(root, "rga_fa18_reduced_cross_section_phi_covariance.npz")
    rgk = find_one(root, "rgk_fa18_reduced_cross_section_phi_covariance.npz")
    return rga, rgk

def safe_corr(cov: np.ndarray) -> np.ndarray:
    d = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    den = np.outer(d, d)
    out = np.zeros_like(cov, dtype=float)
    np.divide(cov, den, out=out, where=den > 0)
    np.fill_diagonal(out, np.where(d > 0, 1.0, 0.0))
    return out

def extract_cell_correlations(df: pd.DataFrame, npz_path: Path, campaign: str):
    z = np.load(npz_path, allow_pickle=True)
    if campaign == "rga":
        cov4 = z["combined_covariance_phi"]
        mask4 = z["combined_final_validity_mask"]
    else:
        cov4 = z["rgk6535_covariance_phi"]
        mask4 = z["rgk6535_final_validity_mask"]

    records = []
    matrices = {}
    for (iq, ix, it), g in df.groupby(["iq2", "ixb", "it"], sort=True):
        iq, ix, it = int(iq), int(ix), int(it)
        valid_npz = np.flatnonzero(mask4[iq, ix, it])
        phi_ids = np.sort(g["iphi"].astype(int).unique())
        # CSV validity and NPZ validity must describe the same selected phi bins.
        if not np.array_equal(valid_npz, phi_ids):
            raise RuntimeError(
                f"{campaign} validity mismatch at ({iq},{ix},{it}): "
                f"CSV={phi_ids.tolist()} NPZ={valid_npz.tolist()}"
            )
        full = cov4[iq, ix, it]
        sub = full[np.ix_(phi_ids, phi_ids)]
        corr = safe_corr(sub)
        key = f"{campaign}_{iq}_{ix}_{it}"
        matrices[key] = corr

        if len(phi_ids) > 1:
            off = corr[np.triu_indices(len(phi_ids), 1)]
            medabs = float(np.median(np.abs(off)))
            maxabs = float(np.max(np.abs(off)))
        else:
            medabs = maxabs = 0.0
        eig = np.linalg.eigvalsh((corr + corr.T) / 2)
        records.append(dict(
            campaign=campaign, iq2=iq, ixb=ix, it=it, nphi=len(phi_ids),
            median_abs_offdiag_corr=medabs,
            max_abs_offdiag_corr=maxabs,
            min_corr_eigenvalue=float(eig.min()),
            max_corr_eigenvalue=float(eig.max()),
        ))
    return pd.DataFrame(records), matrices

def common_reference_points(rga: pd.DataFrame, rgk: pd.DataFrame):
    """
    One common physics point per exact common nominal (Q2,xB,t) cell.
    We use the midpoint of the shared nominal bin edges, NOT either campaign's
    flux coordinate. This prevents an RGA/RGK model difference from masquerading
    as epsilon dependence in the Rosenbluth projection.
    """
    edge = ["Q2_low_GeV2","Q2_high_GeV2","xB_low","xB_high",
            "minus_t_low_GeV2","minus_t_high_GeV2"]
    ag = rga.groupby(edge, as_index=False).agg(
        iq2_rga=("iq2","first"), ixb_rga=("ixb","first"), it_rga=("it","first"),
        epsilon_rga=("virtual_photon_epsilon","median"),
        nphi_rga=("iphi","nunique"),
        median_rel_unc_rga=("relative_uncertainty","median"))
    bg = rgk.groupby(edge, as_index=False).agg(
        iq2_rgk=("iq2","first"), ixb_rgk=("ixb","first"), it_rgk=("it","first"),
        epsilon_rgk=("virtual_photon_epsilon","median"),
        nphi_rgk=("iphi","nunique"),
        median_rel_unc_rgk=("relative_uncertainty","median"))
    c = ag.merge(bg, on=edge, how="inner", validate="one_to_one")
    c["Q2_common_GeV2"] = 0.5*(c["Q2_low_GeV2"]+c["Q2_high_GeV2"])
    c["xB_common"] = 0.5*(c["xB_low"]+c["xB_high"])
    c["minus_t_common_GeV2"] = 0.5*(c["minus_t_low_GeV2"]+c["minus_t_high_GeV2"])
    c["xi_common"] = c["xB_common"]/(2.0-c["xB_common"])
    c["delta_epsilon"] = c["epsilon_rga"]-c["epsilon_rgk"]
    c["both_nphi_ge_8"] = (c["nphi_rga"] >= 8) & (c["nphi_rgk"] >= 8)
    return c.sort_values(["Q2_low_GeV2","xB_low","minus_t_low_GeV2"]).reset_index(drop=True)

def projected_fractional_table(df: pd.DataFrame, campaign: str, exposure_factor: float):
    out = df[["iq2","ixb","it","iphi","Q2_flux_coordinate_GeV2",
              "xB_flux_coordinate","minus_t_center_GeV2","phi_center_deg",
              "virtual_photon_epsilon","relative_uncertainty"]].copy()
    out["campaign"] = campaign
    out["exposure_factor_relative_to_supplied"] = exposure_factor
    # First sensitivity approximation. The supplied uncertainty contains data
    # statistics + finite response/radiative MC statistics, so this is explicitly
    # tagged as an all-components-scaling approximation until those components
    # can be separated.
    out["projected_relative_uncertainty_all_components_scale"] = (
        out["relative_uncertainty"] / math.sqrt(exposure_factor)
    )
    return out

def write_partons_common_kinematics(common: pd.DataFrame, path: Path):
    # Stage-1 convention retained for model-query compatibility.
    with path.open("w") as f:
        for r in common.itertuples(index=False):
            q2 = r.Q2_common_GeV2
            f.write(f"{r.xi_common:.12g}|{-r.minus_t_common_GeV2:.12g}|"
                    f"{q2:.12g}|{q2:.12g}|{q2:.12g}\n")

def make_figures(rga, rgk, common, corr_summary, outdir, rga_factor, rgk_factor):
    outdir.mkdir(parents=True, exist_ok=True)

    # Correlation strength from covariance NPZ.
    fig, ax = plt.subplots(figsize=(7.2,5.4))
    for camp, g in corr_summary.groupby("campaign"):
        ax.hist(g["median_abs_offdiag_corr"], bins=30, histtype="step",
                linewidth=1.6, label=camp.upper())
    ax.set_xlabel("Median |within-cell phi correlation|")
    ax.set_ylabel("(Q², xB, t) cells")
    ax.set_title("Within-cell covariance structure")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir/"01_phi_correlation_strength.png", dpi=180)
    plt.close(fig)

    # Supplied vs simple exposure-scaled fractional precision.
    fig, ax = plt.subplots(figsize=(7.2,5.4))
    cap = 1.0
    for vals, lab in [
        (rga["relative_uncertainty"], "RGA supplied Fa18"),
        (rga["relative_uncertainty"]/np.sqrt(rga_factor), f"RGA current ({rga_factor:g}x approx.)"),
        (rgk["relative_uncertainty"], "RGK supplied subset"),
        (rgk["relative_uncertainty"]/np.sqrt(rgk_factor), f"RGK final ({rgk_factor:g}x approx.)")]:
        ax.hist(np.clip(vals,0,cap), bins=35, histtype="step", linewidth=1.4, label=lab)
    ax.set_xlabel("Projected relative uncertainty")
    ax.set_ylabel("4D bins")
    ax.set_title("Exposure scaling diagnostic")
    ax.text(0.98,0.95,"Values >1 shown in final bin",transform=ax.transAxes,
            ha="right",va="top",fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir/"02_exposure_scaled_relative_uncertainties.png", dpi=180)
    plt.close(fig)

    # Rosenbluth leverage vs Q2 at common points.
    fig, ax = plt.subplots(figsize=(7.2,5.4))
    sc = ax.scatter(common["Q2_common_GeV2"], np.abs(common["delta_epsilon"]),
                    c=common["minus_t_common_GeV2"], s=28)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$-t$ (GeV$^2$)")
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax.set_ylabel(r"$|\Delta\epsilon|$")
    ax.set_title("Common-cell Rosenbluth lever arm")
    fig.tight_layout()
    fig.savefig(outdir/"03_rosenbluth_lever_arm_vs_q2.png", dpi=180)
    plt.close(fig)

    # Final-RGK simple projected precision vs Q2 for common cells.
    tmp = common.copy()
    tmp["rgk_final_median_rel_approx"] = tmp["median_rel_unc_rgk"]/np.sqrt(rgk_factor)
    fig, ax = plt.subplots(figsize=(7.2,5.4))
    ax.scatter(tmp["Q2_common_GeV2"], tmp["rgk_final_median_rel_approx"], s=28)
    ax.set_xlabel(r"$Q^2$ (GeV$^2$)")
    ax.set_ylabel("Approx. final RGK median relative uncertainty")
    ax.set_title("RGK precision at common Rosenbluth cells")
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(outdir/"04_rgk_final_precision_vs_q2.png", dpi=180)
    plt.close(fig)

def main():
    a = parse_args()
    if a.rga_current_factor <= 0 or a.rgk_final_factor <= 0:
        raise ValueError("Exposure factors must be positive.")

    s1 = a.stage1.resolve()
    out = a.output.resolve()
    tables = out/"tables"
    queries = out/"model_queries"
    figures = out/"figures"
    for d in (out,tables,queries,figures):
        d.mkdir(parents=True, exist_ok=True)

    rga = pd.read_csv(s1/"01_blinded_rga_bins.csv")
    rgk = pd.read_csv(s1/"02_blinded_rgk_bins.csv")
    assert_blinded(rga, "RGA Stage-1")
    assert_blinded(rgk, "RGK Stage-1")

    root, tmp = locate_package(a.input.resolve())
    try:
        rga_npz, rgk_npz = covariance_sources(root)
        rga_corr, rga_mats = extract_cell_correlations(rga, rga_npz, "rga")
        rgk_corr, rgk_mats = extract_cell_correlations(rgk, rgk_npz, "rgk")
        corr = pd.concat([rga_corr,rgk_corr], ignore_index=True)
        corr.to_csv(tables/"01_covariance_correlation_summary.csv", index=False)

        # Store only dimensionless correlation matrices, never absolute covariance.
        np.savez_compressed(tables/"02_within_cell_phi_correlations.npz",
                            **rga_mats, **rgk_mats)

        common = common_reference_points(rga, rgk)
        common.to_csv(tables/"03_common_rosenbluth_model_points.csv", index=False)
        write_partons_common_kinematics(
            common, queries/"01_partons_common_rosenbluth_kinematics.csv")

        prga = projected_fractional_table(rga, "RGA", a.rga_current_factor)
        prgk = projected_fractional_table(rgk, "RGK", a.rgk_final_factor)
        pd.concat([prga,prgk], ignore_index=True).to_csv(
            tables/"04_exposure_scaled_fractional_uncertainties.csv", index=False)

        make_figures(rga,rgk,common,corr,figures,
                     a.rga_current_factor,a.rgk_final_factor)

        summary = [
            "CLAS12 pi0 blinded projection -- stage 2",
            "="*45,
            "",
            "No preliminary CLAS12 cross-section central values were read.",
            "Absolute NPZ covariance magnitudes were discarded; only dimensionless",
            "within-cell phi correlation matrices were retained.",
            "",
            f"RGA covariance cells: {len(rga_corr)}",
            f"RGK covariance cells: {len(rgk_corr)}",
            f"Common Rosenbluth model points: {len(common)}",
            f"Common points with >=8 phi bins at both energies: {int(common.both_nphi_ge_8.sum())}",
            f"Median |Delta epsilon|: {np.abs(common.delta_epsilon).median():.6f}",
            "",
            f"Configured current-RGA exposure factor relative to supplied Fa18: {a.rga_current_factor:g}",
            f"Configured final-RGK exposure factor relative to supplied subset: {a.rgk_final_factor:g}",
            "",
            "IMPORTANT uncertainty caveat:",
            "The simple 1/sqrt(exposure) columns scale the entire supplied relative",
            "uncertainty. The source package states that this uncertainty includes data",
            "counting plus finite response-GEMC and radiative-GEMC statistics. The NPZ",
            "does not decompose those components, so this is a sensitivity approximation,",
            "not yet the final luminosity-error model.",
            "",
            "Common Rosenbluth model coordinates use nominal-bin midpoints, not either",
            "campaign's flux coordinates. This prevents model bin-centering differences",
            "from being interpreted as epsilon dependence.",
            "",
            "Stage 2 does not invent GK structure functions. The common PARTONS query",
            "grid is prepared for the next model-evaluation step.",
            "",
            f"Input tar SHA256: {sha256(a.input.resolve())}",
        ]
        (out/"summary.txt").write_text("\n".join(summary)+"\n")
        print("\n".join(summary))
        print(f"\nWrote stage-2 outputs to {out}")
    finally:
        if tmp is not None:
            tmp.cleanup()

if __name__ == "__main__":
    main()
