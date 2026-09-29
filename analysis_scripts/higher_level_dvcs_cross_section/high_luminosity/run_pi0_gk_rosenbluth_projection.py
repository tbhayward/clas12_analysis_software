#!/usr/bin/env python3
"""
Bridge the validated PARTONS/GK Stage-3 physical structure functions into the
existing CLAS12 pi0 RGA+RGK Rosenbluth projection.

This script does NOT invoke PARTONS.  It:
  1. reads the validated shared-point GK physical structure-function grid;
  2. maps the explicit physical column names onto the Stage-3 projection schema;
  3. writes a durable projection-input CSV;
  4. optionally runs prepare_pi0_gk_stage3_projection.py with that CSV.

All four structure functions are in nb/GeV^2.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument(
        "--gk-grid",
        type=Path,
        default=here / "output" / "pi0_gk_stage3" / "partons_gk" /
                "production_grid" / "gk_pi0_shared_physical_structure_functions.csv",
        help="Validated PARTONS/GK physical structure-function grid.",
    )
    p.add_argument(
        "--stage3",
        type=Path,
        default=here / "output" / "pi0_gk_stage3",
        help="Existing Stage-3 projection directory.",
    )
    p.add_argument(
        "--projection-script",
        type=Path,
        default=here / "prepare_pi0_gk_stage3_projection.py",
        help="Existing Stage-3 projection script.",
    )
    p.add_argument(
        "--no-run",
        action="store_true",
        help="Only write the projection-input CSV; do not launch the projection.",
    )
    return p.parse_args()


def main():
    a = parse_args()

    grid = a.gk_grid.expanduser().resolve()
    stage3 = a.stage3.expanduser().resolve()
    projection_script = a.projection_script.expanduser().resolve()

    if not grid.exists():
        raise FileNotFoundError(f"Validated GK grid not found: {grid}")

    d = pd.read_csv(grid)

    required = [
        "point_id",
        "Q2_shared_GeV2",
        "xB_shared",
        "minus_t_shared_GeV2",
        "epsilon_rga",
        "epsilon_rgk",
        "dsigma_T_dt_nb_per_GeV2",
        "dsigma_L_dt_nb_per_GeV2",
        "dsigma_TT_dt_nb_per_GeV2",
        "dsigma_LT_dt_nb_per_GeV2",
    ]
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise RuntimeError(f"{grid} is missing required columns: {missing}")

    out = pd.DataFrame({
        "point_id": d["point_id"].astype(str),
        "sigma_T": d["dsigma_T_dt_nb_per_GeV2"].astype(float),
        "sigma_L": d["dsigma_L_dt_nb_per_GeV2"].astype(float),
        "sigma_TT": d["dsigma_TT_dt_nb_per_GeV2"].astype(float),
        "sigma_LT": d["dsigma_LT_dt_nb_per_GeV2"].astype(float),
        "Q2_shared_GeV2": d["Q2_shared_GeV2"].astype(float),
        "xB_shared": d["xB_shared"].astype(float),
        "minus_t_shared_GeV2": d["minus_t_shared_GeV2"].astype(float),
        "epsilon_rga": d["epsilon_rga"].astype(float),
        "epsilon_rgk": d["epsilon_rgk"].astype(float),
    })

    physical = np.isfinite(
        out[["sigma_T", "sigma_L", "sigma_TT", "sigma_LT"]].to_numpy(dtype=float)
    ).all(axis=1)
    if not physical.all():
        bad = out.loc[~physical, "point_id"].tolist()
        raise RuntimeError(
            f"Validated grid contains {len(bad)} non-finite structure-function rows: "
            + ", ".join(bad[:10])
        )

    # Independent physical sanity checks on the exported GK responses.
    if not (out["sigma_T"] > 0).all():
        raise RuntimeError("Non-positive sigma_T in validated GK grid")
    if not (np.abs(out["sigma_TT"]) <= out["sigma_T"] + 1e-12).all():
        bad=out.loc[np.abs(out.sigma_TT)>out.sigma_T,"point_id"].tolist()
        raise RuntimeError(f"|sigma_TT| > sigma_T for {len(bad)} rows: {bad[:10]}")
    phi=np.linspace(0.0,2.0*np.pi,361)[:,None]
    for epscol in ("epsilon_rga","epsilon_rgk"):
        eps=out[epscol].to_numpy(float)[None,:]
        sig=(out.sigma_T.to_numpy(float)[None,:] + eps*out.sigma_L.to_numpy(float)[None,:]
             + eps*np.cos(2*phi)*out.sigma_TT.to_numpy(float)[None,:]
             + np.sqrt(2*eps*(1+eps))*np.cos(phi)*out.sigma_LT.to_numpy(float)[None,:])
        if np.nanmin(sig) < -1e-10:
            j=np.unravel_index(np.nanargmin(sig),sig.shape)[1]
            raise RuntimeError(f"Negative GK phi-differential bracket for {out.point_id.iloc[j]} at {epscol}: {np.nanmin(sig):.6g}")

    if out["point_id"].duplicated().any():
        dup = out.loc[out["point_id"].duplicated(False), "point_id"].unique().tolist()
        raise RuntimeError("Duplicate point_id values: " + ", ".join(dup[:10]))

    # The projection itself is homogeneous in cross-section units, but record
    # the now-verified physical unit explicitly for downstream provenance.
    out["cross_section_unit"] = "nb/GeV^2"
    out["model"] = "PARTONS_GK06_GPDGK19_LO"

    qdir = stage3 / "model_queries"
    qdir.mkdir(parents=True, exist_ok=True)
    outfile = qdir / "02_gk_structure_function_results.csv"
    out.to_csv(outfile, index=False)

    print("[GK -> Rosenbluth projection bridge]")
    print(f"  input physical grid : {grid}")
    print(f"  valid shared points : {len(out)}")
    print(f"  output projection   : {outfile}")
    print("  units               : nb/GeV^2")
    print(
        "  sigma_T range       : "
        f"{out.sigma_T.min():.6g} -- {out.sigma_T.max():.6g} nb/GeV^2"
    )
    print(
        "  sigma_L range       : "
        f"{out.sigma_L.min():.6g} -- {out.sigma_L.max():.6g} nb/GeV^2"
    )

    if a.no_run:
        print("  projection launch   : skipped (--no-run)")
        return

    if not projection_script.exists():
        raise FileNotFoundError(
            f"Projection script not found: {projection_script}\n"
            "The bridge CSV was written successfully; pass the correct script with "
            "--projection-script."
        )

    cmd = [
        sys.executable,
        str(projection_script),
        "--gk-results",
        str(outfile),
    ]
    print("\nLaunching existing Stage-3 projection:")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
