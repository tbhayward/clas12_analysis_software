#!/usr/bin/env python3
"""Step-6b audit: joint Rosenbluth covariance, L/T propagation, and pure-stat scaling."""

import argparse
import importlib.util
import math
from pathlib import Path
import numpy as np
import pandas as pd


def load_module(path):
    spec = importlib.util.spec_from_file_location("reach", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def blockdiag(a, b):
    out = np.zeros((len(a)+len(b), len(a)+len(b)))
    out[:len(a), :len(a)] = a
    out[len(a):, len(a):] = b
    return out


def independent_gls(A, y, C):
    Ci = np.linalg.pinv(0.5*(C+C.T), rcond=1e-12)
    N = A.T @ Ci @ A
    cov = np.linalg.inv(N)
    beta = cov @ A.T @ Ci @ y
    return beta, cov


def ratio_unc(beta, cov):
    T, L = float(beta[0]), float(beta[1])
    if T == 0:
        return np.nan
    g = np.array([-L/T**2, 1/T, 0.0, 0.0])
    return math.sqrt(max(float(g @ cov @ g), 0.0))


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--study-script", type=Path, default=here/"pi0_LT_measurement_reach_study.py")
    ap.add_argument("--future-scan", default="0,1,3,5,10")
    ap.add_argument("--output", type=Path,
                    default=here/"output/pi0_gk_stage3/diagnostics/joint_rosenbluth_covariance_audit.csv")
    a = ap.parse_args()
    m = load_module(a.study_script.resolve())

    common_file = here/"output/pi0_gk_stage2/tables/03_common_rosenbluth_model_points.csv"
    rga_file = here/"import/fa18_rosenbluth_inputs_20260924T165948Z/rga_10604/combined_reduced_cross_sections.csv"
    rgk_file = here/"import/fa18_rosenbluth_inputs_20260924T165948Z/rgk_6535/rgk6535_reduced_cross_sections.csv"
    corrections = here/"output/pi0_gk_stage3/partons_gk/06_gk_native_to_shared_phi_corrections.csv"
    covariance = here/"output/pi0_gk_stage2/tables/02_within_cell_phi_correlations.npz"
    exclusions = here/"output/pi0_gk_stage3/tables/00_gk_model_exclusions.csv"

    # Production results with normalization nuisance disabled: pure statistical closure target.
    common = pd.read_csv(common_file)
    prod = m.build_internal_data_extraction(
        common, rga_file, rgk_file, corrections, covariance, exclusions,
        min_delta_epsilon=0.05, relative_norm_unc=0.0
    ).set_index("point_id")

    # Reproduce the production preprocessing so that the independent calculation starts
    # from the same corrected measured phi points, but does its own matrix algebra.
    rga = m._standardize_internal_cross_sections(rga_file, "RGA")
    rgk = m._standardize_internal_cross_sections(rgk_file, "RGK")
    corr = pd.read_csv(corrections)
    covz = np.load(covariance, allow_pickle=False)
    excluded = set(pd.read_csv(exclusions).point_id.astype(str))

    rows = []
    for row in common.itertuples(index=False):
        pid = str(row.point_id)
        if pid in excluded or pid not in prod.index:
            continue

        q2, xb, mt = float(row.Q2_GeV2), float(row.xB), float(row.minus_t_GeV2)
        ir = (int(row.iq2_rga), int(row.ixb_rga), int(row.it_rga))
        ik = (int(row.iq2_rgk), int(row.ixb_rgk), int(row.it_rgk))
        ga = m._select_internal_group(rga, row, "rga")
        gk = m._select_internal_group(rgk, row, "rgk")

        # Apply exactly the production reduced-GK continuous-phi correction.
        for camp, g in (("rga", ga), ("rgk", gk)):
            cc = corr[(corr.point_id.astype(str)==pid) & (corr.campaign.str.lower()==camp)].copy()
            cc = cc.sort_values("phi_deg")
            ph = np.deg2rad(cc.phi_deg.to_numpy(float))
            H = np.column_stack([np.ones(len(cc)), np.cos(ph), np.cos(2*ph)])
            nc = np.linalg.solve(H, cc.reduced_value_native.to_numpy(float))
            sc = np.linalg.solve(H, cc.reduced_value_shared.to_numpy(float))
            pm = np.deg2rad(g.phi_deg.to_numpy(float))
            Hm = np.column_stack([np.ones(len(g)), np.cos(pm), np.cos(2*pm)])
            fac = (Hm@sc)/(Hm@nc)
            g = g.copy()
            g["sigma"] *= fac
            g["delta_sigma"] *= np.abs(fac)
            if camp == "rga":
                ga = g
            else:
                gk = g
        #endfor

        er = m._epsilon_from_q2_xb_E(q2, xb, 10.604)
        ek = m._epsilon_from_q2_xb_E(q2, xb, 6.535)
        ga["epsilon"] = er
        gk["epsilon"] = ek

        ga, Ca = m._cell_covariance(ga, covz, "rga", *ir, rel_norm=0.0)
        gk, Ck = m._cell_covariance(gk, covz, "rgk", *ik, rel_norm=0.0)

        As, ys = [], []
        for g in (ga, gk):
            ph = np.deg2rad(g.phi_deg.to_numpy(float))
            ep = g.epsilon.to_numpy(float)
            As.append(np.column_stack([
                np.ones(len(g)), ep,
                np.sqrt(2*ep*(1+ep))*np.cos(ph),
                ep*np.cos(2*ph)
            ]))
            ys.append(g.sigma.to_numpy(float))
        #endfor
        A = np.vstack(As)
        y = np.concatenate(ys)
        C = blockdiag(Ca, Ck)

        beta, V = independent_gls(A, y, C)
        beta_p, V_p, _, _, _ = m._joint_rosenbluth_fit(ga, gk, Ca, Ck)

        dR = ratio_unc(beta, V)
        V_no = V.copy()
        V_no[0,1] = V_no[1,0] = 0.0
        dR_no = ratio_unc(beta, V_no)
        rho = V[0,1]/math.sqrt(V[0,0]*V[1,1])

        pr = prod.loc[pid]
        max_beta_diff = float(np.max(np.abs(beta-beta_p)))
        max_cov_diff = float(np.max(np.abs(V-V_p)))
        dR_prod = float(pr.delta_R_L_over_T)

        scaling_residual = 0.0
        scaling_ratios = []
        for f in [float(x) for x in a.future_scan.split(",")]:
            # Here remaining == recorded for both campaigns, so pure-stat covariance
            # at final exposure is C/(1+f).
            sf = 1.0/(1.0+f)
            bf, Vf = independent_gls(A, y, C*sf)
            dRf = ratio_unc(bf, Vf)
            expected = dR/math.sqrt(1.0+f)
            rel = dRf/expected - 1.0 if expected > 0 else np.nan
            if np.isfinite(rel):
                scaling_residual = max(scaling_residual, abs(rel))
            scaling_ratios.append((f, dRf/dR if dR > 0 else np.nan))
        #endfor

        rows.append(dict(
            point_id=pid,
            max_abs_beta_independent_minus_production=max_beta_diff,
            max_abs_cov_independent_minus_production=max_cov_diff,
            delta_R_independent=dR,
            delta_R_production=dR_prod,
            abs_delta_R_difference=abs(dR-dR_prod),
            rho_TL=rho,
            delta_R_without_TL_cov=dR_no,
            ratio_delta_R_no_cov_over_full=(dR_no/dR if dR > 0 else np.nan),
            max_pure_stat_scaling_relative_residual=scaling_residual,
            condition_number=float(pr.condition_number),
        ))
    #endfor

    out = pd.DataFrame(rows)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.output, index=False)

    def report(name, x):
        z = np.asarray(x, float)
        z = z[np.isfinite(z)]
        print(f"{name}: median={np.median(z):.6g}, p95={np.percentile(z,95):.6g}, max={np.max(z):.6g}")

    print("\n[Step 6b: production joint-GLS closure]")
    print(f"cells audited: {len(out)}")
    print(f"max |independent beta - production beta|: {out.max_abs_beta_independent_minus_production.max():.6e}")
    print(f"max |independent covariance - production covariance|: {out.max_abs_cov_independent_minus_production.max():.6e}")
    print(f"max |independent delta(L/T) - production|: {out.abs_delta_R_difference.max():.6e}")

    print("\n[Actual production T-L covariance]")
    report("|rho_TL|", np.abs(out.rho_TL))
    report("deltaR(no TL cov) / deltaR(full)", out.ratio_delta_R_no_cov_over_full)

    print("\n[Pure-stat effective-statistics scaling closure]")
    print(f"max relative residual from 1/sqrt(1+f): {out.max_pure_stat_scaling_relative_residual.max():.6e}")
    print("Expected uncertainty factors relative to recorded data:")
    for f in [float(x) for x in a.future_scan.split(",")]:
        print(f"  remaining-data multiplier f={f:g}: final effective statistics={1+f:g}x recorded, "
              f"stat uncertainty factor={1/math.sqrt(1+f):.6f}")
    #endfor

    print(f"\nWrote: {a.output}")


if __name__ == "__main__":
    main()
