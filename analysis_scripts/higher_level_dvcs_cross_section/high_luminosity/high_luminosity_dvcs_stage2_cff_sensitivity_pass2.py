#!/usr/bin/env python3
"""
Stage 2B (pass-2 baseline): joint XS+BSA local CFF-H sensitivity projection.

Inputs:
  output_stage2_pass2/tables/joint_fit_input_km15_{1,2,5,10}x.csv

Pseudo-truth:
  KM15 evaluated at the exact pass-2 CLAS12 kinematics.

Fit modes:
  ReH_only
  ReH_ImH_XS
  ReH_ImH_XS_BSA

The last two modes provide the central workshop comparison: how much does
matched beam-spin-asymmetry information break the ReH/ImH degeneracy that is
present in the unpolarized cross section alone?

CFF changes are propagated by finite changes through Gepard.  For BSA the
script calculates the beam-helicity asymmetry from the +1 and -1 helicity
cross sections:
    A_LU = (sigma+ - sigma-) / (sigma+ + sigma-)
A global sign convention does not affect the projected covariance because it
multiplies both BSA derivatives by the same sign.

Nuisances:
  * one shared XS scale nuisance, with the per-point fractional scale prior
    supplied by the pass-2 CSV;
  * one shared 4% BSA beam-polarization scale nuisance.

The pseudo-data equal KM15, so all best-fit CFF shifts are zero.  This stage
projects uncertainties and correlations only.
"""

from __future__ import annotations

import argparse
import contextlib
import math
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)
FIT_MODES = ("ReH_only", "ReH_ImH_XS", "ReH_ImH_XS_BSA")
DEFAULT_MIN_NPHI = 3
DEFAULT_REL_STEP = 0.02
DEFAULT_ABS_STEP = 0.05


def make_point(g, row, helicity: int = 0):
    phi_trento = math.pi - math.radians(float(row.phi_deg))
    pt = g.DataPoint(
        xB=float(row.xB),
        t=-abs(float(row.t_abs)),
        Q2=float(row.Q2),
        phi=float(phi_trento),
        observable="XS",
        frame="trento",
        process="ep2epgamma",
        exptype="fixed target",
        in1energy=float(row.ebeam),
        in1charge=-1,
        in1polarization=int(helicity),
        in2particle="p",
    )
    pt.prepare()
    return pt


def cff_owner(th, name: str):
    candidates = [th.m] if hasattr(th, "m") else []
    candidates.append(th)
    seen = set()
    for obj in candidates:
        if id(obj) in seen:
            continue
        seen.add(id(obj))
        if hasattr(obj, name) and callable(getattr(obj, name)):
            return obj
    raise AttributeError(f"Could not find callable Gepard CFF method {name}")


@contextlib.contextmanager
def shifted_cff(th, name: str, shift: float):
    owner = cff_owner(th, name)
    old = getattr(owner, name)
    setattr(owner, name, lambda pt, _old=old, _shift=float(shift):
            float(_old(pt)) + _shift)
    try:
        yield
    finally:
        setattr(owner, name, old)


def pred(th, pt) -> float:
    return float(th.predict(pt))


def cff_value(th, name: str, pt) -> float:
    return float(getattr(cff_owner(th, name), name)(pt))


def bsa_from_points(th, plus_pt, minus_pt) -> float:
    sp = pred(th, plus_pt)
    sm = pred(th, minus_pt)
    den = sp + sm
    if not np.isfinite(den) or abs(den) < 1e-30:
        return np.nan
    return (sp - sm) / den


def finite_derivative(th, points: Sequence, observable: str, name: str,
                      nominal_cff: float, rel_step: float,
                      abs_step: float) -> Tuple[np.ndarray, float]:
    step = max(abs_step, rel_step * max(abs(nominal_cff), 1.0))

    def evaluate():
        if observable == "xs":
            return np.asarray([pred(th, p) for p in points], dtype=float)
        return np.asarray(
            [bsa_from_points(th, pair[0], pair[1]) for pair in points],
            dtype=float,
        )

    with shifted_cff(th, name, +step):
        plus = evaluate()
    with shifted_cff(th, name, -step):
        minus = evaluate()
    return (plus - minus) / (2.0 * step), step


def load_inputs(indir: Path):
    out = {}
    for factor in LUMI_FACTORS:
        path = indir / f"joint_fit_input_km15_{factor}x.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        d = pd.read_csv(path)
        d["bin"] = d["bin"].astype(int)
        out[factor] = d.sort_values(["bin", "phi_deg"]).reset_index(drop=True)
    ref = out[1][["point_id", "bin"]].astype(str)
    for factor in LUMI_FACTORS[1:]:
        if not ref.equals(out[factor][["point_id", "bin"]].astype(str)):
            raise RuntimeError("Luminosity input kinematics/order mismatch.")
    return out


def build_derivatives(df, th, g, rel_step, abs_step, min_nphi):
    point_rows, cell_rows = [], []
    groups = list(df.groupby("bin", sort=True))
    print(f"[CFF] calculating XS+BSA sensitivities for {len(groups)} cells")

    for icell, (bid, cell) in enumerate(groups, 1):
        cell = cell.sort_values("phi_deg")
        if len(cell) < min_nphi:
            print(f"[CFF] bin {bid}: skipped, only {len(cell)} phi points")
            continue

        xs_pts = [make_point(g, r, 0) for r in cell.itertuples(index=False)]
        bsa_pts = [
            (make_point(g, r, +1), make_point(g, r, -1))
            for r in cell.itertuples(index=False)
        ]
        pt0 = xs_pts[0]
        reh = cff_value(th, "ReH", pt0)
        imh = cff_value(th, "ImH", pt0)

        dx_re, step_re = finite_derivative(
            th, xs_pts, "xs", "ReH", reh, rel_step, abs_step)
        dx_im, step_im = finite_derivative(
            th, xs_pts, "xs", "ImH", imh, rel_step, abs_step)
        da_re, _ = finite_derivative(
            th, bsa_pts, "bsa", "ReH", reh, rel_step, abs_step)
        da_im, _ = finite_derivative(
            th, bsa_pts, "bsa", "ImH", imh, rel_step, abs_step)

        xs0 = np.asarray([pred(th, p) for p in xs_pts])
        a0 = np.asarray([bsa_from_points(th, p, m) for p, m in bsa_pts])

        for j, row in enumerate(cell.itertuples(index=False)):
            point_rows.append({
                "point_id": row.point_id, "bin": int(bid),
                "phi_deg": float(row.phi_deg),
                "xs_km15": float(xs0[j]), "bsa_km15": float(a0[j]),
                "d_xs_d_ReH": float(dx_re[j]),
                "d_xs_d_ImH": float(dx_im[j]),
                "d_bsa_d_ReH": float(da_re[j]),
                "d_bsa_d_ImH": float(da_im[j]),
            })

        cell_rows.append({
            "bin": int(bid),
            "xB": float(np.median(cell["xB"])),
            "Q2": float(np.median(cell["Q2"])),
            "t_abs": float(np.median(cell["t_abs"])),
            "t_over_Q2": float(np.median(cell["t_over_Q2"])),
            "n_phi": len(cell),
            "ReH_KM15": reh, "ImH_KM15": imh,
            "finite_step_ReH": step_re,
            "finite_step_ImH": step_im,
            "median_abs_dBSA_dImH": float(np.nanmedian(np.abs(da_im))),
            "median_abs_dBSA_dReH": float(np.nanmedian(np.abs(da_re))),
        })

        if icell % 10 == 0 or icell == len(groups):
            print(f"[CFF] {icell:3d}/{len(groups)} cells")

    return pd.DataFrame(point_rows), pd.DataFrame(cell_rows)


def layout_for(cell_meta, mode):
    pars = []
    for b in cell_meta["bin"].astype(int):
        pars.append((b, "ReH"))
        if mode != "ReH_only":
            pars.append((b, "ImH"))
    pars.append((-1, "xs_scale_beta"))
    if mode == "ReH_ImH_XS_BSA":
        pars.append((-1, "bsa_scale_beta"))
    return pars


def covariance(df, deriv, cell_meta, mode):
    use = df.merge(deriv, on=["point_id", "bin", "phi_deg"], how="inner")
    layout = layout_for(cell_meta, mode)
    idx = {k: i for i, k in enumerate(layout)}
    npar = len(layout)

    rows, sigmas = [], []

    # XS rows
    for r in use.itertuples(index=False):
        v = np.zeros(npar)
        b = int(r.bin)
        v[idx[(b, "ReH")]] = r.d_xs_d_ReH
        if mode != "ReH_only":
            v[idx[(b, "ImH")]] = r.d_xs_d_ImH
        # beta=1 corresponds to the quoted per-point correlated scale shift.
        v[idx[(-1, "xs_scale_beta")]] = r.xs_scale_frac * r.xs_km15
        s = math.hypot(r.xs_stat_pseudo_abs, r.xs_ptp_sys_pseudo_abs)
        if np.isfinite(s) and s > 0:
            rows.append(v)
            sigmas.append(s)

    # Matched BSA rows only in the joint mode.
    if mode == "ReH_ImH_XS_BSA":
        for r in use.itertuples(index=False):
            v = np.zeros(npar)
            b = int(r.bin)
            v[idx[(b, "ReH")]] = r.d_bsa_d_ReH
            v[idx[(b, "ImH")]] = r.d_bsa_d_ImH
            v[idx[(-1, "bsa_scale_beta")]] = r.bsa_scale_frac * r.bsa_km15
            s = math.hypot(r.bsa_stat_pseudo_abs, r.bsa_ptp_sys_pseudo_abs)
            if np.isfinite(s) and s > 0:
                rows.append(v)
                sigmas.append(s)

    A = np.asarray(rows)
    s = np.asarray(sigmas)
    B = A / s[:, None]

    # Unit Gaussian priors on nuisance betas.
    priors = []
    for key in [(-1, "xs_scale_beta"), (-1, "bsa_scale_beta")]:
        if key in idx:
            p = np.zeros(npar)
            p[idx[key]] = 1.0
            priors.append(p)
    B = np.vstack([B] + priors)

    info = B.T @ B
    rank = int(np.linalg.matrix_rank(B))
    sv = np.linalg.svd(B, compute_uv=False)
    pos = sv[sv > max(sv[0], 1.0)*1e-12]
    cond = float(pos[0]/pos[-1]) if len(pos) > 1 else np.inf
    cov = np.linalg.pinv(info, rcond=1e-12)

    diag = {
        "fit_mode": mode,
        "n_measurement_rows": len(A),
        "n_parameters": npar,
        "rank": rank,
        "rank_deficit": npar-rank,
        "condition_number_effective": cond,
        "xs_scale_beta_sigma": math.sqrt(max(cov[idx[(-1,"xs_scale_beta")],
                                                   idx[(-1,"xs_scale_beta")]], 0)),
        "bsa_scale_beta_sigma": np.nan,
    }
    if (-1, "bsa_scale_beta") in idx:
        j = idx[(-1, "bsa_scale_beta")]
        diag["bsa_scale_beta_sigma"] = math.sqrt(max(cov[j,j],0))
    return cov, layout, diag


def summarize(cell_meta, cov, layout, mode, factor):
    idx = {k:i for i,k in enumerate(layout)}
    rows = []
    for c in cell_meta.itertuples(index=False):
        b = int(c.bin)
        ir = idx[(b,"ReH")]
        sr = math.sqrt(max(cov[ir,ir],0))
        row = {
            "fit_mode": mode, "luminosity_factor": factor, "bin": b,
            "xB": c.xB, "Q2": c.Q2, "t_abs": c.t_abs,
            "t_over_Q2": c.t_over_Q2, "n_phi": c.n_phi,
            "ReH_KM15": c.ReH_KM15, "ImH_KM15": c.ImH_KM15,
            "sigma_ReH": sr,
            "relative_sigma_ReH": sr/abs(c.ReH_KM15)
                if abs(c.ReH_KM15)>1e-12 else np.nan,
            "sigma_ImH": np.nan, "relative_sigma_ImH": np.nan,
            "corr_ReH_ImH": np.nan,
        }
        if mode != "ReH_only":
            ii = idx[(b,"ImH")]
            si = math.sqrt(max(cov[ii,ii],0))
            row["sigma_ImH"] = si
            row["relative_sigma_ImH"] = si/abs(c.ImH_KM15) \
                if abs(c.ImH_KM15)>1e-12 else np.nan
            den = math.sqrt(max(cov[ir,ir]*cov[ii,ii],0))
            if den > 0:
                row["corr_ReH_ImH"] = cov[ir,ii]/den
        rows.append(row)
    return pd.DataFrame(rows)


def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def make_plots(results, figures):
    # Central comparison: XS-only vs matched XS+BSA for the two-CFF fit.
    fig, ax = plt.subplots(figsize=(8.8,6.0))
    for mode, label in [
        ("ReH_ImH_XS", "XS only"),
        ("ReH_ImH_XS_BSA", "XS + BSA"),
    ]:
        d = results[(results.fit_mode==mode) &
                    (results.luminosity_factor==10)]
        ax.scatter(d.t_abs, 100*d.relative_sigma_ReH, s=30, alpha=.72,
                   label=label)
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"Relative uncertainty on $\mathrm{Re}\,\mathcal{H}$ (%)")
    ax.set_title(r"Effect of matched BSA information at 10$\times$ luminosity")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures/"01_ReH_XS_vs_XS_BSA_10x.png")

    fig, ax = plt.subplots(figsize=(8.8,6.0))
    for mode, label in [
        ("ReH_ImH_XS", "XS only"),
        ("ReH_ImH_XS_BSA", "XS + BSA"),
    ]:
        d = results[(results.fit_mode==mode) &
                    (results.luminosity_factor==10)]
        ax.scatter(d.t_abs, np.abs(d.corr_ReH_ImH), s=30, alpha=.72,
                   label=label)
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$|\rho(\mathrm{Re}\,\mathcal{H},\mathrm{Im}\,\mathcal{H})|$")
    ax.set_ylim(0,1.05)
    ax.set_title(r"Does matched BSA break the $\mathrm{Re}H$--$\mathrm{Im}H$ degeneracy?")
    ax.grid(alpha=.2); ax.legend()
    savefig(fig, figures/"02_ReH_ImH_correlation_XS_vs_XS_BSA_10x.png")

    for quantity, ylabel, fname in [
        ("relative_sigma_ReH", r"Relative uncertainty on $\mathrm{Re}\,\mathcal{H}$ (%)",
         "03_joint_ReH_vs_t.png"),
        ("relative_sigma_ImH", r"Relative uncertainty on $\mathrm{Im}\,\mathcal{H}$ (%)",
         "04_joint_ImH_vs_t.png"),
    ]:
        fig, ax = plt.subplots(figsize=(8.8,6.0))
        for factor in LUMI_FACTORS:
            d = results[(results.fit_mode=="ReH_ImH_XS_BSA") &
                        (results.luminosity_factor==factor)]
            ax.scatter(d.t_abs, 100*d[quantity], s=27, alpha=.68,
                       label=f"{factor}x")
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(ylabel)
        ax.set_title("Matched pass-2 XS+BSA projection")
        ax.grid(alpha=.2); ax.legend(title="Luminosity")
        savefig(fig, figures/fname)


def main(argv: Optional[List[str]] = None) -> int:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir",
                   default=str(here/"output_stage2_pass2"/"tables"))
    p.add_argument("--outdir",
                   default=str(here/"output_stage2_pass2_cff"))
    p.add_argument("--min-nphi", type=int, default=DEFAULT_MIN_NPHI)
    p.add_argument("--relative-step", type=float, default=DEFAULT_REL_STEP)
    p.add_argument("--absolute-step", type=float, default=DEFAULT_ABS_STEP)
    args = p.parse_args(argv)

    indir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    tables = outdir/"tables"; figures=outdir/"figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    inputs = load_inputs(indir)
    one = inputs[1]

    try:
        import gepard as g
        from gepard.fits import th_KM15
    except Exception as exc:
        raise RuntimeError("Could not import Gepard/KM15.") from exc
    th = th_KM15

    # Preflight at a representative point.
    row = one.iloc[len(one)//2]
    xs0 = make_point(g, row, 0)
    bp = make_point(g, row, +1); bm = make_point(g, row, -1)
    print("[preflight] representative KM15 point")
    print(f"[preflight] XS={pred(th,xs0):.8g}, BSA={bsa_from_points(th,bp,bm):+.8g}")
    for name in ("ReH","ImH"):
        val = cff_value(th,name,xs0)
        step = max(args.absolute_step,args.relative_step*max(abs(val),1.0))
        with shifted_cff(th,name,step):
            xs1=pred(th,xs0); a1=bsa_from_points(th,bp,bm)
        print(
            f"[preflight] {name}: value={val:.6g}, step={step:.6g}, "
            f"XS response={100*(xs1/pred(th,xs0)-1):+.5g}%, "
            f"BSA change={a1-bsa_from_points(th,bp,bm):+.6g}"
        )

    deriv, cell_meta = build_derivatives(
        one, th, g, args.relative_step, args.absolute_step, args.min_nphi)
    deriv.to_csv(tables/"point_xs_bsa_cff_derivatives.csv", index=False)
    cell_meta.to_csv(tables/"cell_km15_cffs.csv", index=False)

    result_frames=[]; diag_rows=[]
    for mode in FIT_MODES:
        for factor in LUMI_FACTORS:
            print(f"[fit] mode={mode:16s} luminosity={factor}x")
            cov, layout, diag = covariance(
                inputs[factor], deriv, cell_meta, mode)
            diag["luminosity_factor"]=factor
            diag_rows.append(diag)
            result_frames.append(summarize(
                cell_meta,cov,layout,mode,factor))

    results=pd.concat(result_frames,ignore_index=True)
    diagnostics=pd.DataFrame(diag_rows)
    results.to_csv(tables/"cff_cell_uncertainties.csv",index=False)
    diagnostics.to_csv(tables/"global_fit_diagnostics.csv",index=False)
    make_plots(results,figures)

    print("\n"+"="*92)
    print("PASS-2 MATCHED XS+BSA CFF SENSITIVITY SUMMARY")
    print("="*92)
    print(f"usable cells: {len(cell_meta)}")
    rows=[]
    for mode in FIT_MODES:
        for factor in LUMI_FACTORS:
            d=results[(results.fit_mode==mode)&
                      (results.luminosity_factor==factor)]
            rows.append({
                "mode":mode,"L":f"{factor}x",
                "median_rel_ReH_pct":100*np.nanmedian(d.relative_sigma_ReH),
                "median_rel_ImH_pct":100*np.nanmedian(d.relative_sigma_ImH)
                    if np.any(np.isfinite(d.relative_sigma_ImH)) else np.nan,
                "median_abs_corr":np.nanmedian(np.abs(d.corr_ReH_ImH))
                    if np.any(np.isfinite(d.corr_ReH_ImH)) else np.nan,
                "ReH_lt25pct":int(np.sum(d.relative_sigma_ReH<.25)),
                "ReH_lt50pct":int(np.sum(d.relative_sigma_ReH<.50)),
            })
    print(pd.DataFrame(rows).to_string(
        index=False,float_format=lambda x:f"{x:.3g}"))
    print("\n[output]",outdir)
    print("[next] compare ReH_ImH_XS directly with ReH_ImH_XS_BSA.")
    print("[next] if BSA breaks the degeneracy, proceed to a controlled D(t) layer.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
