#!/usr/bin/env python3
"""
Stage 2A (pass-2 baseline): prepare matched CLAS12 XS+BSA Asimov pseudo-data.

The workshop projection intentionally does NOT use unpublished pass-2 central
values as pseudo-truth.  It uses:
  * pass-2 kinematics and measured uncertainties/coverage;
  * KM15 as the smooth pseudo-truth;
  * the preferred |t|/Q^2 < 0.20 theory-control cut;
  * luminosity factors 1, 2, 5, 10 relative to the existing pass-2 exposure.

Input:
  ../import/dvcs_pass2_analysis.csv

Output:
  output_stage2_pass2/tables/joint_fit_input_km15_{1,2,5,10}x.csv

The tuple-valued CSV columns are parsed as (value, statistical error, ...).
For the XS projection:
  stat uncertainty      -> scales as 1/sqrt(L)
  point-to-point syst   -> fixed
  correlated scale      -> fixed; per-point fractional scale column retained

For the BSA projection:
  stat uncertainty      -> scales as 1/sqrt(L)
  point-to-point syst   -> fixed, using quadrature of cut and pi0-subtraction
  beam-polarization     -> fixed 4% correlated multiplicative nuisance

The BSA central values from pass-2 are read only to obtain the experimental
uncertainties/systematics and for optional diagnostics.  The fit pseudo-truth
is generated later directly from KM15.
"""

from __future__ import annotations

import argparse
import ast
import math
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)
DEFAULT_CLEAN_RATIO = 0.20
DEFAULT_EBEAM = 10.604
DEFAULT_BSA_SCALE_FRAC = 0.04


def tuple_component(series: pd.Series, index: int) -> np.ndarray:
    out = np.full(len(series), np.nan, dtype=float)
    for i, value in enumerate(series):
        if pd.isna(value):
            continue
        try:
            parsed = ast.literal_eval(str(value).strip())
            if isinstance(parsed, (tuple, list)) and len(parsed) > index:
                out[i] = float(parsed[index])
        except Exception:
            pass
    return out


def num(df: pd.DataFrame, name: str, default=np.nan) -> np.ndarray:
    if name not in df.columns:
        return np.full(len(df), float(default))
    return pd.to_numeric(df[name], errors="coerce").to_numpy(float)


def load_pass2(path: Path, clean_ratio: float, ebeam: float) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    out = pd.DataFrame(index=raw.index)

    out["row_id"] = np.arange(len(raw), dtype=int)
    out["bin"] = num(raw, "bin")
    out["xB"] = num(raw, "xBavg, 10.6 GeV")
    out["Q2"] = num(raw, "Q2avg, 10.6 GeV")
    out["t_abs"] = num(raw, "t_abs_avg, 10.6 GeV")
    out["phi_deg"] = num(raw, "phiavg, 10.6 GeV")
    out["ebeam"] = float(ebeam)

    xs_col = "normed cross sections, ep->epg, exp, 10.6 GeV, unpol"
    bsa_col = "BSA, counts, 10.6 GeV"
    if xs_col not in raw or bsa_col not in raw:
        raise KeyError("Required pass-2 XS/BSA columns are missing from the CSV.")

    out["xs_data"] = tuple_component(raw[xs_col], 0)
    out["xs_stat_abs"] = tuple_component(raw[xs_col], 1)
    out["bsa_data"] = tuple_component(raw[bsa_col], 0)
    out["bsa_stat_abs"] = tuple_component(raw[bsa_col], 1)

    # Final pass-2 point-to-point XS systematic is already an absolute error.
    out["xs_ptp_sys_abs"] = num(raw, "Syst. err (point-to-point total)")

    # The combined 10.6-GeV scale uncertainty is stored as a fractional
    # uncertainty and can vary slightly with the run-period combination.
    out["xs_scale_frac"] = num(
        raw,
        "normed cross sections, ep->epg, exp, 10.6 GeV, unpol, total scale sys",
    )

    # For BSA, the cut systematic combines exclusivity+fiducial variations.
    # pi0 subtraction is a separate point-to-point contribution.
    bsa_cut = num(raw, "BSA, counts, 10.6 GeV, cut sys", default=0.0)
    bsa_pi0 = num(raw, "BSA, counts, 10.6 GeV, pi0 subtraction sys", default=0.0)
    out["bsa_ptp_sys_abs"] = np.hypot(
        np.nan_to_num(bsa_cut, nan=0.0),
        np.nan_to_num(bsa_pi0, nan=0.0),
    )

    # The CSV confirms the beam-polarization contribution is 4% of |A_LU|
    # for nonzero measured asymmetries.  Store the fractional nuisance rather
    # than an absolute error so it can be applied consistently to KM15 truth.
    out["bsa_scale_frac"] = float(DEFAULT_BSA_SCALE_FRAC)

    out["t_over_Q2"] = out["t_abs"] / out["Q2"]

    finite_kin = np.all(
        np.isfinite(out[["xB", "Q2", "t_abs", "phi_deg"]].to_numpy(float)),
        axis=1,
    )
    has_xs = np.isfinite(out["xs_stat_abs"]) & (out["xs_stat_abs"] > 0)
    has_bsa = np.isfinite(out["bsa_stat_abs"]) & (out["bsa_stat_abs"] > 0)
    clean = np.isfinite(out["t_over_Q2"]) & (out["t_over_Q2"] < clean_ratio)

    out = out[finite_kin & has_xs & has_bsa & clean].copy()
    out["bin"] = out["bin"].astype(int)
    out["point_id"] = [
        f"pass2_{int(b)}_{i:04d}" for i, b in enumerate(out["bin"].to_numpy())
    ]
    return out.reset_index(drop=True)


def main(argv: Optional[list[str]] = None) -> int:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument(
        "--pass2-csv",
        default=str(here.parent / "import" / "dvcs_pass2_analysis.csv"),
    )
    p.add_argument(
        "--outdir",
        default=str(here / "output_stage2_pass2"),
    )
    p.add_argument("--clean-ratio", type=float, default=DEFAULT_CLEAN_RATIO)
    p.add_argument("--ebeam", type=float, default=DEFAULT_EBEAM)
    args = p.parse_args(argv)

    source = Path(args.pass2_csv).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    tables = outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    d = load_pass2(source, args.clean_ratio, args.ebeam)

    print("=" * 82)
    print("PASS-2 MATCHED XS+BSA PROJECTION BASELINE")
    print("=" * 82)
    print(f"source                         : {source}")
    print(f"matched clean 4D points        : {len(d)}")
    print(f"matched clean (xB,Q2,t) cells  : {d['bin'].nunique()}")
    print(
        f"|t| range (GeV^2)              : "
        f"{d['t_abs'].min():.3f} -- {d['t_abs'].max():.3f}"
    )
    print(
        f"Q2 range (GeV^2)               : "
        f"{d['Q2'].min():.3f} -- {d['Q2'].max():.3f}"
    )
    print(
        f"median XS stat/data             : "
        f"{100*np.nanmedian(d['xs_stat_abs']/np.abs(d['xs_data'])):.2f}%"
    )
    print(
        f"median BSA statistical error    : "
        f"{np.nanmedian(d['bsa_stat_abs']):.4f}"
    )
    print(
        f"median BSA ptp systematic       : "
        f"{np.nanmedian(d['bsa_ptp_sys_abs']):.4f}"
    )
    print(
        f"median XS scale prior           : "
        f"{100*np.nanmedian(d['xs_scale_frac']):.2f}%"
    )
    print(f"BSA beam-polarization prior     : {100*DEFAULT_BSA_SCALE_FRAC:.2f}%")

    d.to_csv(tables / "pass2_matched_clean_experimental_precision.csv", index=False)

    for factor in LUMI_FACTORS:
        q = d.copy()
        q["luminosity_factor"] = factor
        q["xs_stat_pseudo_abs"] = q["xs_stat_abs"] / math.sqrt(factor)
        q["bsa_stat_pseudo_abs"] = q["bsa_stat_abs"] / math.sqrt(factor)

        # Systematics are deliberately held fixed in the first luminosity-only
        # projection.  Later we can add explicit detector/systematic-improvement
        # scenarios without changing the fit machinery.
        q["xs_ptp_sys_pseudo_abs"] = q["xs_ptp_sys_abs"]
        q["bsa_ptp_sys_pseudo_abs"] = q["bsa_ptp_sys_abs"]

        path = tables / f"joint_fit_input_km15_{factor}x.csv"
        q.to_csv(path, index=False)
        print(f"[write] {path}")

    print("\nProjection convention:")
    print("  1x = existing pass-2 RGA exposure represented by this CSV")
    print("  2x/5x/10x = same exposure time at higher luminosity")
    print("  statistical errors scale as 1/sqrt(L)")
    print("  point-to-point and correlated systematics are fixed")
    print("  unpublished pass-2 central values are NOT used as pseudo-truth")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
