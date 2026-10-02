#!/usr/bin/env python3
"""
Step-6 audit for the pi0 Rosenbluth luminosity-reach study.

Checks:
  1) joint-fit sigma_L against the independent two-epsilon U separation;
  2) statistical sigma_L uncertainty against the analytic two-point formula;
  3) exact exposure meaning of every future-luminosity multiplier;
  4) the covariance algebra used for L/T in the model projection.

This is a diagnostic only. It does not modify production outputs.
"""

import argparse
import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--study-script", type=Path,
                   default=here/"pi0_LT_measurement_reach_study.py")
    p.add_argument("--future-lumi-scan", default="0,1,3,5,10")
    p.add_argument("--rga-recorded", type=float, default=1.5)
    p.add_argument("--rga-remaining", type=float, default=1.5)
    p.add_argument("--rgk-recorded", type=float, default=9.0)
    p.add_argument("--rgk-remaining", type=float, default=9.0)
    p.add_argument("--output", type=Path,
                   default=here/"output"/"pi0_gk_stage3"/"diagnostics"/
                           "rosenbluth_step6_audit.csv")
    a = p.parse_args()

    m = load_module(a.study_script.resolve(), "pi0_reach_study")
    # Recreate the production default paths without altering the production script.
    common_file = here/"output"/"pi0_gk_stage2"/"tables"/"03_common_rosenbluth_model_points.csv"
    rga_file = here/"import"/"fa18_rosenbluth_inputs_20260924T165948Z"/\
                    "rga_10604"/"combined_reduced_cross_sections.csv"
    rgk_file = here/"import"/"fa18_rosenbluth_inputs_20260924T165948Z"/\
                    "rgk_6535"/"rgk6535_reduced_cross_sections.csv"
    corrections = here/"output"/"pi0_gk_stage3"/"partons_gk"/\
                       "06_gk_native_to_shared_phi_corrections.csv"
    covariance = here/"output"/"pi0_gk_stage2"/"tables"/\
                      "02_within_cell_phi_correlations.npz"
    exclusions = here/"output"/"pi0_gk_stage3"/"tables"/"00_gk_model_exclusions.csv"

    common = pd.read_csv(common_file)
    data = m.build_internal_data_extraction(
        common, rga_file, rgk_file, corrections, covariance, exclusions,
        min_delta_epsilon=0.05, relative_norm_unc=0.0
    )

    # With relative_norm_unc=0, the independent-U uncertainty is purely statistical.
    rows = []
    for r in data.itertuples(index=False):
        de = float(r.delta_epsilon)
        Ur, Uk = float(r.sigma_U_rga), float(r.sigma_U_rgk)
        dUr, dUk = float(r.delta_sigma_U_rga), float(r.delta_sigma_U_rgk)

        L_analytic = (Ur-Uk)/de
        varL_analytic = (dUr*dUr+dUk*dUk)/(de*de)
        dL_analytic = math.sqrt(varL_analytic)

        # Analytic T/L covariance if U_R and U_K are independent and the
        # separation is performed from the two independently fitted U values.
        er, ek = float(r.epsilon_rga), float(r.epsilon_rgk)
        T_analytic = (er*Uk-ek*Ur)/de
        varT_analytic = (ek*ek*dUr*dUr + er*er*dUk*dUk)/(de*de)
        covTL_analytic = -(ek*dUr*dUr + er*dUk*dUk)/(de*de)

        if T_analytic != 0.0:
            R_analytic = L_analytic/T_analytic
            grad = np.array([-L_analytic/T_analytic**2, 1.0/T_analytic])
            Ctl = np.array([[varT_analytic, covTL_analytic],
                            [covTL_analytic, varL_analytic]])
            dR_analytic = math.sqrt(max(float(grad @ Ctl @ grad), 0.0))
            # Deliberately omit covariance to quantify its importance.
            Cdiag = np.diag([varT_analytic, varL_analytic])
            dR_no_cov = math.sqrt(max(float(grad @ Cdiag @ grad), 0.0))
        else:
            R_analytic = dR_analytic = dR_no_cov = np.nan

        rows.append(dict(
            point_id=r.point_id,
            delta_epsilon=de,
            sigma_T_joint=float(r.sigma_T),
            sigma_T_from_U=T_analytic,
            delta_T_joint_minus_U=float(r.sigma_T)-T_analytic,
            sigma_L_joint=float(r.sigma_L),
            sigma_L_from_U=L_analytic,
            delta_L_joint_minus_U=float(r.sigma_L)-L_analytic,
            delta_sigma_L_joint=float(r.delta_sigma_L),
            delta_sigma_L_from_U=dL_analytic,
            ratio_delta_sigma_L_joint_over_U=(
                float(r.delta_sigma_L)/dL_analytic if dL_analytic > 0 else np.nan
            ),
            cov_TL_from_U=covTL_analytic,
            corr_TL_from_U=(
                covTL_analytic/math.sqrt(varT_analytic*varL_analytic)
                if varT_analytic > 0 and varL_analytic > 0 else np.nan
            ),
            R_L_over_T_from_U=R_analytic,
            delta_R_with_TL_cov_from_U=dR_analytic,
            delta_R_without_TL_cov_from_U=dR_no_cov,
            ratio_delta_R_no_cov_over_with_cov=(
                dR_no_cov/dR_analytic if np.isfinite(dR_analytic) and dR_analytic > 0
                else np.nan
            ),
            chi2_ndf=float(r.chi2_ndf),
            joint_condition_number=float(r.condition_number),
        ))

    out = pd.DataFrame(rows)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.output, index=False)

    def stats(x):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        return np.median(x), np.percentile(x, 95), np.max(x)

    print("\n[Step 6: two-epsilon Rosenbluth audit]")
    print(f"cells audited: {len(out)}")
    print(f"max |L_joint - L_from_U|: "
          f"{np.nanmax(np.abs(out.delta_L_joint_minus_U)):.6g} nb/GeV^2")
    print(f"max |T_joint - T_from_U|: "
          f"{np.nanmax(np.abs(out.delta_T_joint_minus_U)):.6g} nb/GeV^2")

    med, p95, mx = stats(out.ratio_delta_sigma_L_joint_over_U)
    print("deltaL_joint / deltaL_independent-U:")
    print(f"  median={med:.6g}, p95={p95:.6g}, max={mx:.6g}")

    med, p95, mx = stats(np.abs(out.corr_TL_from_U))
    print("|corr(T,L)| from analytic independent-U separation:")
    print(f"  median={med:.6g}, p95={p95:.6g}, max={mx:.6g}")

    med, p95, mx = stats(out.ratio_delta_R_no_cov_over_with_cov)
    print("delta(L/T) if T-L covariance were OMITTED / correct analytic value:")
    print(f"  median={med:.6g}, p95={p95:.6g}, max={mx:.6g}")

    print("\n[Exposure interpretation]")
    fs = [float(x) for x in a.future_lumi_scan.split(",")]
    exposure_rows = []
    for f in fs:
        rga = a.rga_recorded + f*a.rga_remaining
        rgk = a.rgk_recorded + f*a.rgk_remaining
        exposure_rows.append(dict(
            future_multiplier=f,
            rga_factor_vs_supplied=rga,
            rgk_factor_vs_supplied=rgk,
            rga_factor_vs_recorded=rga/a.rga_recorded,
            rgk_factor_vs_recorded=rgk/a.rgk_recorded,
        ))
    exp = pd.DataFrame(exposure_rows)
    print(exp.to_string(index=False))

    print("\nImportant: with the current defaults, recorded == remaining for BOTH campaigns.")
    print("Therefore f=1 means twice the currently recorded statistics,")
    print("f=3 means four times recorded, f=5 means six times recorded,")
    print("and f=10 means eleven times recorded.  The absolute factors versus")
    print("the supplied source samples remain different (RGA 1.5 vs RGK 9 at f=0).")

    print(f"\nWrote: {a.output}")


if __name__ == "__main__":
    main()
