#!/usr/bin/env python3
import argparse
import hashlib
import tarfile
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RGA_REL = Path("rga_10604/combined_reduced_cross_sections.csv")
RGK_REL = Path("rgk_6535/rgk6535_reduced_cross_sections.csv")

BLINDED_COLUMNS = [
    "campaign_key","campaign_label","iq2","ixb","it","iphi",
    "Q2_low_GeV2","Q2_high_GeV2","Q2_flux_coordinate_GeV2",
    "xB_low","xB_high","xB_flux_coordinate","beam_energy_GeV",
    "virtual_photon_epsilon","minus_t_low_GeV2","minus_t_high_GeV2",
    "minus_t_center_GeV2","phi_low_deg","phi_high_deg","phi_center_deg",
    "relative_uncertainty",
]
FORBIDDEN = {
    "reduced_cross_section_nb_per_GeV2_rad",
    "propagated_statistical_and_finite_MC_uncertainty_nb_per_GeV2_rad",
}

def cli():
    p = argparse.ArgumentParser(description="Prepare blinded CLAS12 pi0 kinematics and PARTONS/GK inputs.")
    p.add_argument("--input", type=Path, default=None, help="Input tar.gz or extracted package directory; default: ../import/fa18_rosenbluth_inputs_20260924T165948Z.tar.gz relative to this script")
    p.add_argument("--outdir", type=Path, default=Path("output/pi0_gk_stage1"))
    p.add_argument("--partons-kinematics-path", default=None)
    p.add_argument("--n-warmups", type=int, default=10000)
    p.add_argument("--n-calls", type=int, default=100000)
    p.add_argument("--chi2-limit", type=float, default=0.8)
    return p.parse_args()

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()

def locate(inp):
    inp = inp.expanduser().resolve()
    if inp.is_dir():
        return inp, None
    if not tarfile.is_tarfile(inp):
        raise ValueError("Input must be the package directory or tar archive.")
    tmp = tempfile.TemporaryDirectory(prefix="pi0_gk_stage1_")
    root = Path(tmp.name)
    with tarfile.open(inp, "r:*") as tf:
        for m in tf.getmembers():
            q = Path(m.name)
            if q.is_absolute() or ".." in q.parts:
                raise ValueError(f"Unsafe archive member: {m.name}")
        tf.extractall(root, filter="data")
    for c in [root] + list(root.iterdir()):
        if c.is_dir() and (c/RGA_REL).exists() and (c/RGK_REL).exists():
            return c, tmp
    raise FileNotFoundError("Expected RGA/RGK files not found in package.")

def load_blinded(path, label):
    src = pd.read_csv(path)
    missing = [c for c in BLINDED_COLUMNS if c not in src.columns]
    if missing:
        raise ValueError(f"{label}: missing columns {missing}")
    if not FORBIDDEN.issubset(src.columns):
        raise ValueError(f"{label}: source schema is not the expected provisional export.")
    # Hard privacy boundary: only this explicit allow-list survives.
    out = src.loc[:, BLINDED_COLUMNS].copy()
    out.insert(0, "sample", label)
    return out

def add_gk(df):
    out = df.copy()
    xb = out["xB_flux_coordinate"].astype(float)
    q2 = out["Q2_flux_coordinate_GeV2"].astype(float)
    out["gk_xi"] = xb/(2.0-xb)
    out["gk_t_GeV2"] = -out["minus_t_center_GeV2"].astype(float)
    out["gk_MuF2_GeV2"] = q2
    out["gk_MuR2_GeV2"] = q2
    return out

def make_unique_points(rga, rgk):
    z = pd.concat([rga, rgk], ignore_index=True)
    cols = ["Q2_flux_coordinate_GeV2","xB_flux_coordinate","minus_t_center_GeV2",
            "gk_xi","gk_t_GeV2","gk_MuF2_GeV2","gk_MuR2_GeV2"]
    u = z[cols].drop_duplicates().sort_values(cols[:3]).reset_index(drop=True)
    u.insert(0, "gk_point_id", np.arange(len(u), dtype=int))
    return u

def make_common_cells(rga, rgk):
    k = ["iq2","ixb","it"]
    a = rga.groupby(k, as_index=False).agg(
        rga_nphi=("iphi","nunique"),
        rga_Q2=("Q2_flux_coordinate_GeV2","median"),
        rga_xB=("xB_flux_coordinate","median"),
        rga_minus_t=("minus_t_center_GeV2","median"),
        rga_epsilon=("virtual_photon_epsilon","median"),
        rga_median_rel_unc=("relative_uncertainty","median"))
    b = rgk.groupby(k, as_index=False).agg(
        rgk_nphi=("iphi","nunique"),
        rgk_Q2=("Q2_flux_coordinate_GeV2","median"),
        rgk_xB=("xB_flux_coordinate","median"),
        rgk_minus_t=("minus_t_center_GeV2","median"),
        rgk_epsilon=("virtual_photon_epsilon","median"),
        rgk_median_rel_unc=("relative_uncertainty","median"))
    m = a.merge(b, on=k, how="inner")
    m["delta_epsilon_rga_minus_rgk"] = m["rga_epsilon"]-m["rgk_epsilon"]
    m["both_nphi_ge_8"] = (m["rga_nphi"] >= 8) & (m["rgk_nphi"] >= 8)
    return m.sort_values(k).reset_index(drop=True)

def write_partons_kinematics(points, path):
    # Official PARTONS DVMP computeManyKinematic format:
    # xi|t|MuF2|MuR2|Q2
    with open(path, "w") as f:
        for r in points.itertuples(index=False):
            f.write(f"{r.gk_xi:.12g}|{r.gk_t_GeV2:.12g}|{r.gk_MuF2_GeV2:.12g}|"
                    f"{r.gk_MuR2_GeV2:.12g}|{r.Q2_flux_coordinate_GeV2:.12g}\n")

def write_xml(path, kin_path, nw, nc, chi2):
    lines = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>',
        '<scenario date="2026-09-24" description="CLAS12 blinded pi0 projection: GK DVMP evaluation">',
        '  <task service="DVMPConvolCoeffFunctionService" method="computeManyKinematic" storeInDB="0">',
        '    <kinematics type="DVMPConvolCoeffFunctionKinematic">',
        f'      <param name="file" value="{kin_path}" />',
        '    </kinematics>',
        '    <computation_configuration>',
        '      <module type="DVMPConvolCoeffFunctionModule" name="DVMPCFFGK06">',
        f'        <param name="nWarmUps" value="{nw}" />',
        f'        <param name="nCalls" value="{nc}" />',
        f'        <param name="chi2Limit" value="{chi2}" />',
        '        <param name="qcd_order_type" value="LO" />',
        '        <module type="GPDModule" name="GPDGK19"></module>',
        '      </module>',
        '    </computation_configuration>',
        '  </task>',
        '  <task service="DVMPConvolCoeffFunctionService" method="printResults"></task>',
        '</scenario>',
    ]
    path.write_text("\n".join(lines)+"\n")


def make_diagnostic_plots(rga, rgk, common, plotdir):
    plotdir.mkdir(parents=True, exist_ok=True)

    # One point per nominal Q2-xB cell so phi/t multiplicity does not dominate.
    rga_qx = rga.groupby(["iq2","ixb"], as_index=False).agg(
        Q2=("Q2_flux_coordinate_GeV2","median"),
        xB=("xB_flux_coordinate","median"))
    rgk_qx = rgk.groupby(["iq2","ixb"], as_index=False).agg(
        Q2=("Q2_flux_coordinate_GeV2","median"),
        xB=("xB_flux_coordinate","median"))

    fig, ax = plt.subplots(figsize=(7.2,5.4))
    ax.scatter(rga_qx["xB"], rga_qx["Q2"], marker="o", label="RGA 10.604 GeV")
    ax.scatter(rgk_qx["xB"], rgk_qx["Q2"], marker="x", label="RGK 6.535 GeV")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
    ax.set_title(r"Blinded $\pi^0$ kinematic coverage")
    ax.legend()
    fig.tight_layout()
    fig.savefig(plotdir/"01_q2_xb_coverage.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2,5.4))
    ax.hist(np.abs(common["delta_epsilon_rga_minus_rgk"]), bins=18)
    ax.axvline(np.abs(common["delta_epsilon_rga_minus_rgk"]).median(),
               linestyle="--", label="Median")
    ax.set_xlabel(r"$|\Delta\epsilon|$")
    ax.set_ylabel("Common (Q², xB, t) cells")
    ax.set_title("RGA–RGK Rosenbluth lever arm")
    ax.legend()
    fig.tight_layout()
    fig.savefig(plotdir/"02_delta_epsilon_distribution.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2,5.4))
    cap = 1.5
    ax.hist(np.clip(rga["relative_uncertainty"], 0, cap), bins=35,
            histtype="step", linewidth=1.6, label="RGA supplied sample")
    ax.hist(np.clip(rgk["relative_uncertainty"], 0, cap), bins=35,
            histtype="step", linewidth=1.6, label="RGK supplied subset")
    ax.set_xlabel("Relative uncertainty")
    ax.set_ylabel("4D bins")
    ax.set_title("Supplied per-bin relative uncertainties")
    ax.legend()
    ax.text(0.98, 0.95, "Values >1.5 shown in final bin",
            transform=ax.transAxes, ha="right", va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(plotdir/"03_relative_uncertainty_distribution.png", dpi=180)
    plt.close(fig)

    q2_edges = rga[["iq2","Q2_low_GeV2","Q2_high_GeV2"]].drop_duplicates("iq2")
    counts = common.groupby("iq2").size().rename("n_common").reset_index()
    counts8 = common[common["both_nphi_ge_8"]].groupby("iq2").size().rename("n_ge8").reset_index()
    qc = q2_edges.merge(counts, on="iq2", how="inner").merge(counts8, on="iq2", how="left").fillna({"n_ge8":0})
    qc["label"] = qc.apply(lambda r: f'{r["Q2_low_GeV2"]:g}–{r["Q2_high_GeV2"]:g}', axis=1)

    fig, ax = plt.subplots(figsize=(8.2,5.4))
    x = np.arange(len(qc))
    ax.bar(x, qc["n_common"], label="All common cells")
    ax.plot(x, qc["n_ge8"], marker="o", label=r"$\geq8$ $\phi$ bins at both energies")
    ax.set_xticks(x)
    ax.set_xticklabels(qc["label"], rotation=35, ha="right")
    ax.set_xlabel(r"$Q^2$ bin (GeV$^2$)")
    ax.set_ylabel("Common (Q², xB, t) cells")
    ax.set_title("RGA–RGK common-cell reach versus Q²")
    ax.legend()
    fig.tight_layout()
    fig.savefig(plotdir/"04_common_cells_vs_q2.png", dpi=180)
    plt.close(fig)

def main():
    a = cli()
    if a.input is None:
        a.input = Path(__file__).resolve().parent / "import" / "fa18_rosenbluth_inputs_20260924T165948Z.tar.gz"
    out = a.outdir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    root, tmp = locate(a.input)
    try:
        rp, kp = root/RGA_REL, root/RGK_REL
        rga = add_gk(load_blinded(rp, "RGA_10.604"))
        rgk = add_gk(load_blinded(kp, "RGK_6.535"))

        rga.to_csv(out/"01_blinded_rga_bins.csv", index=False)
        rgk.to_csv(out/"02_blinded_rgk_bins.csv", index=False)

        points = make_unique_points(rga, rgk)
        points.to_csv(out/"03_unique_gk_points.csv", index=False)

        kin = out/"04_partons_dvmp_kinematics.csv"
        write_partons_kinematics(points, kin)

        xml = out/"05_partons_gk_pi0_scenario.xml"
        write_xml(xml, a.partons_kinematics_path or str(kin),
                  a.n_warmups, a.n_calls, a.chi2_limit)

        common = make_common_cells(rga, rgk)
        common.to_csv(out/"06_exact_common_cells.csv", index=False)

        plotdir = out/"figures"
        make_diagnostic_plots(rga, rgk, common, plotdir)

        # Defensive privacy check on every CSV emitted here.
        for p in out.glob("*.csv"):
            bad = FORBIDDEN.intersection(pd.read_csv(p, nrows=0).columns)
            if bad:
                raise RuntimeError(f"Privacy guard failed for {p}: {sorted(bad)}")

        de = common["delta_epsilon_rga_minus_rgk"].abs()
        summary = [
            "CLAS12 pi0 blinded GK preparation -- stage 1",
            "================================================",
            "",
            "Measured CLAS12 reduced-cross-section central values are NOT exported.",
            "Outputs retain kinematics, epsilon, bin IDs, and supplied relative uncertainties.",
            "",
            f"RGA blinded 4D bins: {len(rga)}",
            f"RGK blinded 4D bins: {len(rgk)}",
            f"Unique GK query points: {len(points)}",
            f"Exact common nominal (Q2,xB,t) cells: {len(common)}",
            f"Common cells with >=8 phi bins in both: {int(common['both_nphi_ge_8'].sum())}",
            f"Median |Delta epsilon|: {de.median():.6f}",
            f"Maximum |Delta epsilon|: {de.max():.6f}",
            "",
            "PARTONS query convention: xi=xB/(2-xB), t=-|t|, MuF2=MuR2=Q2.",
            "Scenario: DVMPCFFGK06 + GPDGK19 at LO.",
            "",
            "This stage prepares the model query only. It does not run PARTONS,",
            "parse model partial cross sections, or construct pseudo-data.",
            "",
            f"RGA source SHA256: {sha256(rp)}",
            f"RGK source SHA256: {sha256(kp)}",
        ]
        (out/"07_summary.txt").write_text("\n".join(summary)+"\n")
        print("\n".join(summary))
        print(f"\nWrote stage-1 outputs to {out}")
    finally:
        if tmp is not None:
            tmp.cleanup()

if __name__ == "__main__":
    main()
