#!/usr/bin/env python3
"""
Compare published unpolarized proton DVCS/BH cross-section measurements.

This is the standalone analysis driver for the external-data/model-comparison
chapter of the CLAS12 RGA pass-2 DVCS analysis note.  It deliberately does NOT
load the pass-2 result.  The first objective is to establish the mutual
consistency of the published world datasets before the new pass-2 measurement
is introduced.

Published datasets included
---------------------------
  * CLAS6  Jo et al.      2015
  * Hall A Defurne et al. 2015
  * Hall A Defurne et al. 2017
  * CLAS6  Saylor et al.  2018
  * Hall A Georges et al. 2022
  * CLAS12 Lee et al.     2026 (pass-1)

Models used
-----------
  * KM15, evaluated with Gepard.
  * GK16, evaluated with the existing PARTONS driver
        evaluate_bh_model_selection.py
    using GPDGK16 + DVCSCFFStandard(LO) + DVCSProcessGV08.
  * Pure Bethe-Heitler (BH), evaluated through the same Gepard KM15 machinery;
    only the elastic BH term is retained.

VGG99 is intentionally not used in any scientific comparison produced by this
script.  The shared PARTONS driver currently also evaluates VGG99 internally,
but those columns are ignored here.

Analysis philosophy
-------------------
1. Keep every experimental point at its native beam energy and native mean
   kinematics for the primary model/data consistency tests.
2. Treat published correlated normalization uncertainties as one nuisance per
   experiment rather than adding them independently to every point.
3. Match different experiments only when their measured kinematics are nearby.
4. For a matched A -> B comparison, transport A to B's exact kinematics with a
   LOCAL model ratio:

       sigma_A_to_B^M = sigma_A * sigma_M(k_B) / sigma_M(k_A).

5. For common-energy presentation, transport each point to Ebeam=10.6 GeV at
   fixed (xB,Q2,t,phi):

       sigma_10p6^M = sigma_data * sigma_M(10.6) / sigma_M(E_native).

6. KM15 is nominal.  The absolute KM15-vs-GK16 difference in the transported
   cross section is retained as an explicit transport-model uncertainty.
   NO point is removed merely because this uncertainty is large.

Uncertainty conventions used here
---------------------------------
CLAS6 Jo 2015:
    5% correlated elastic-normalization uncertainty.  It is removed in
    quadrature from the published total point uncertainty by the validated EMFF
    loader and then carried here as a single correlated scale nuisance.

Hall A Defurne 2015:
    Existing validated loader prescription; 2.8% correlated scale.

Hall A Defurne 2017:
    3.2% point-to-point systematic plus correlated
    sqrt(1.0%^2 + 2.0%^2 + 0.5%^2) = 2.291...% scale.

CLAS6 Saylor 2018:
    4% correlated elastic normalization, removed in quadrature from the
    published total systematic by the validated loader.

Hall A Georges 2022:
    Published symmetricized pointwise systematic is retained, and for this
    comparison study a pragmatic 5% correlated normalization prior is allowed,
    as requested for the analysis-note world-data comparison.

CLAS12 Lee 2026:
    31% correlated normalization.  The validated authoritative E214M1 loader
    removes it in quadrature from the released total systematic to recover the
    point-to-point component.

Typical use
-----------
From dvcs_cross_section/external_scripts:

  python compare_world_dvcs_cross_sections.py --workers 8

This study uses KM15 as the sole model-assisted transport prescription for
mapping nearby measurements to one another and for presentation at 10.6 GeV.
No PARTONS/GK16 transport calculation is performed in this script.

The script is intended to live beside extract_emff_from_dvcs_bh.py in
external_scripts/.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import subprocess
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =============================================================================
# Configuration
# =============================================================================

TARGET_EBEAM_GEV = 10.6

# Georges is intentionally different from the production EMFF extraction:
# for this external-data comparison we allow the ~5% normalization freedom
# indicated by the recent world-data EMFF fits.
GEORGES_COMPARISON_NORM_FRAC = 0.05

# Matching cuts are deliberately explicit and configurable.  They define only
# candidate overlap; the final score is continuous within these windows.
DEFAULT_MATCH_DXB = 0.035
DEFAULT_MATCH_DQ2 = 0.60       # GeV^2
DEFAULT_MATCH_DT = 0.10        # GeV^2
DEFAULT_MATCH_DPHI = 20.0      # degrees

DATASET_ORDER = [
    "jo2015",
    "defurne2015",
    "defurne2017",
    "saylor2018",
    "georges2022",
    "lee2026",
]

DATASET_LABELS = {
    "jo2015": "CLAS6 Jo 2015",
    "defurne2015": "Hall A Defurne 2015",
    "defurne2017": "Hall A Defurne 2017",
    "saylor2018": "CLAS6 Saylor 2018",
    "georges2022": "Hall A Georges 2022",
    "lee2026": "CLAS12 Lee 2026",
}

# Source convention needed for KM15 evaluation and for the existing PARTONS
# phi-mapping logic.
GEPARD_BMK_DATASETS = {"jo2015", "defurne2015", "defurne2017"}
DIRECT_PHI_DATASETS = {"saylor2018", "georges2022", "lee2026"}

# evaluate_bh_model_selection.py expects the dataset keys used by the EMFF
# suite.  Keep this translation in exactly one place.
PARTONS_DATASET_KEYS = {
    "jo2015": "jo2015",
    "defurne2015": "halla_defurne2015",
    "defurne2017": "halla_defurne2017",
    "saylor2018": "saylor2018",
    "georges2022": "halla_georges2022",
    "lee2026": "pass1",
}


@dataclass(frozen=True)
class MatchConfig:
    dxb: float = DEFAULT_MATCH_DXB
    dq2: float = DEFAULT_MATCH_DQ2
    dt: float = DEFAULT_MATCH_DT
    dphi: float = DEFAULT_MATCH_DPHI


# =============================================================================
# Small utilities
# =============================================================================


def circular_phi_difference_deg(a, b):
    """Smallest absolute azimuthal separation in degrees, in [0,180]."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.abs((a - b + 180.0) % 360.0 - 180.0)
#enddef


def finite_positive(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return np.isfinite(arr) & (arr > 0.0)
#enddef


def safe_fraction(num, den) -> np.ndarray:
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full(np.broadcast(num, den).shape, np.nan, dtype=float)
    good = np.isfinite(num) & np.isfinite(den) & (den != 0.0)
    out[good] = num[good] / den[good]
    return out
#enddef


def resolve_existing(path: Path, alternatives: Sequence[Path]) -> Path:
    """Resolve one required input while giving useful repository-relative fallbacks."""
    candidates = [Path(path)] + [Path(p) for p in alternatives]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
        #endif
    #endfor
    attempted = "\n".join(f"  {p}" for p in candidates)
    raise FileNotFoundError(f"Required input not found. Tried:\n{attempted}")
#enddef


def load_python_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load Python module from {path}")
    #endif
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
#enddef


# =============================================================================
# Canonicalization of the validated EMFF loaders
# =============================================================================


def _base_canonical(df: pd.DataFrame, key: str, norm_frac: float) -> pd.DataFrame:
    out = pd.DataFrame({
        "dataset": key,
        "dataset_label": DATASET_LABELS[key],
        "source_row": np.arange(len(df), dtype=int),
        "xB": pd.to_numeric(df["xB"], errors="coerce"),
        "Q2": pd.to_numeric(df["Q2"], errors="coerce"),
        "t_abs": pd.to_numeric(df["t_abs"], errors="coerce"),
        "phi_deg": np.mod(pd.to_numeric(df["phi_deg"], errors="coerce"), 360.0),
        "ebeam": pd.to_numeric(df["ebeam"], errors="coerce"),
        "xs": pd.to_numeric(df["xs"], errors="coerce"),
        "norm_frac": float(norm_frac),
    })
    return out
#enddef


def canonicalize_jo(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "jo2015", emff.JO_GLOBAL_SCALE_FRAC)

    stat = pd.to_numeric(df.get("clas6_err_stat", np.nan), errors="coerce").to_numpy(float)
    point_total = pd.to_numeric(df["clas6_err_pointwise"], errors="coerce").to_numpy(float)

    # The pointwise total from the validated loader is the published total with
    # the 5% correlated normalization removed.  Recover the residual systematic
    # only for reporting; the fits use point_unc_abs directly.
    stat_for_quad = np.where(np.isfinite(stat) & (stat >= 0.0), stat, 0.0)
    ptp_sys = np.sqrt(np.maximum(point_total**2 - stat_for_quad**2, 0.0))

    out["stat_abs"] = stat
    out["ptp_sys_abs"] = ptp_sys
    out["point_unc_abs"] = point_total
    out["phi_bmk_rad"] = pd.to_numeric(df["phi_bmk_rad"], errors="coerce")
    out["gepard_dataset_id"] = pd.to_numeric(df["gepard_dataset_id"], errors="coerce")
    return clean_canonical(out)
#enddef


def canonicalize_defurne2015(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "defurne2015", emff.HALLA_DEFURNE_GLOBAL_SCALE_FRAC)
    out["stat_abs"] = pd.to_numeric(df["halla_err_stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["halla_err_syst"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    out["phi_bmk_rad"] = pd.to_numeric(df["phi_bmk_rad"], errors="coerce")
    out["gepard_dataset_id"] = pd.to_numeric(df["gepard_dataset_id"], errors="coerce")
    return clean_canonical(out)
#enddef


def canonicalize_defurne2017(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "defurne2017", emff.HALLA_DEFURNE2017_GLOBAL_SCALE_FRAC)
    out["stat_abs"] = pd.to_numeric(df["halla_err_stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["halla_err_syst"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    out["phi_bmk_rad"] = pd.to_numeric(df["phi_bmk_rad"], errors="coerce")
    out["gepard_dataset_id"] = pd.to_numeric(df["gepard_dataset_id"], errors="coerce")
    return clean_canonical(out)
#enddef


def canonicalize_saylor(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "saylor2018", emff.SAYLOR_GLOBAL_SCALE_FRAC)
    out["stat_abs"] = pd.to_numeric(df["stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["ptp_sys"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    return clean_canonical(out)
#enddef


def canonicalize_georges(df: pd.DataFrame) -> pd.DataFrame:
    out = _base_canonical(df, "georges2022", GEORGES_COMPARISON_NORM_FRAC)
    out["stat_abs"] = pd.to_numeric(df["stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["ptp_sys_abs"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    out["setting"] = df["setting"].astype(str).to_numpy()
    out["t_minus_tmin"] = pd.to_numeric(df["t_minus_tmin"], errors="coerce")
    out["t_min"] = pd.to_numeric(df["t_min"], errors="coerce")
    return clean_canonical(out)
#enddef


def canonicalize_lee(df: pd.DataFrame, emff) -> pd.DataFrame:
    out = _base_canonical(df, "lee2026", emff.PASS1_GLOBAL_SCALE_FRAC)
    out["stat_abs"] = pd.to_numeric(df["xs_stat"], errors="coerce")
    out["ptp_sys_abs"] = pd.to_numeric(df["ptp_sys_abs"], errors="coerce")
    out["point_unc_abs"] = np.hypot(out["stat_abs"], out["ptp_sys_abs"])
    if "bin" in df.columns:
        out["published_bin"] = pd.to_numeric(df["bin"], errors="coerce")
    #endif
    if "released_syst_below_nominal_scale" in df.columns:
        out["released_syst_below_nominal_scale"] = df[
            "released_syst_below_nominal_scale"
        ].astype(bool).to_numpy()
    #endif
    return clean_canonical(out)
#enddef


def clean_canonical(df: pd.DataFrame) -> pd.DataFrame:
    required = [
        "xB", "Q2", "t_abs", "phi_deg", "ebeam", "xs",
        "stat_abs", "ptp_sys_abs", "point_unc_abs", "norm_frac",
    ]
    finite = np.ones(len(df), dtype=bool)
    for col in required:
        finite &= np.isfinite(pd.to_numeric(df[col], errors="coerce").to_numpy(float))
    #endfor
    finite &= df["xs"].to_numpy(float) > 0.0
    finite &= df["point_unc_abs"].to_numpy(float) > 0.0
    finite &= df["Q2"].to_numpy(float) > 0.0
    finite &= df["t_abs"].to_numpy(float) > 0.0

    nbad = int((~finite).sum())
    if nbad:
        warnings.warn(f"Dropping {nbad} nonfinite/nonphysical canonical rows")
    #endif

    out = df.loc[finite].copy().reset_index(drop=True)
    out["source_row"] = np.arange(len(out), dtype=int)
    out["point_id"] = [
        f"{dataset}:{row}"
        for dataset, row in zip(out["dataset"], out["source_row"])
    ]
    out["stat_frac"] = out["stat_abs"] / out["xs"]
    out["ptp_sys_frac"] = out["ptp_sys_abs"] / out["xs"]
    out["point_unc_frac"] = out["point_unc_abs"] / out["xs"]
    return out
#enddef


def load_world_data(args, emff) -> pd.DataFrame:
    """Load all six published datasets through the validated EMFF loaders."""
    print("\n" + "=" * 80)
    print("LOADING PUBLISHED WORLD DATA")
    print("=" * 80)

    jo = canonicalize_jo(emff.load_clas6_gepard_dataset(), emff)
    d15 = canonicalize_defurne2015(emff.load_halla_defurne_gepard_datasets(), emff)
    d17 = canonicalize_defurne2017(emff.load_halla_defurne2017_gepard_datasets(), emff)

    saylor_path = resolve_existing(
        Path(args.saylor_file),
        [args.script_dir / "import" / "saylor_CLAS6.txt"],
    )
    saylor = canonicalize_saylor(emff.load_saylor_supplement(str(saylor_path)), emff)

    georges_path = resolve_existing(
        Path(args.georges_file),
        [args.script_dir / "import" / "E12-06-114.xlsx"],
    )
    georges = canonicalize_georges(emff.load_georges_supplement(str(georges_path)))

    lee_path = resolve_existing(
        Path(args.lee_file),
        [
            args.script_dir / "import" / "clasdb_E214M1.txt",
            args.script_dir.parent / "imports" / "clasdb_E214M1.txt",
        ],
    )
    lee = canonicalize_lee(emff.load_clas12_pass1_csv(lee_path), emff)

    world = pd.concat([jo, d15, d17, saylor, georges, lee], ignore_index=True, sort=False)
    world["dataset"] = pd.Categorical(world["dataset"], DATASET_ORDER, ordered=True)
    world = world.sort_values(["dataset", "source_row"]).reset_index(drop=True)
    world["dataset"] = world["dataset"].astype(str)

    if world["point_id"].duplicated().any():
        raise RuntimeError("Canonical point_id collision detected")
    #endif

    print(f"[WORLD] canonicalized {len(world):,} points across 6 measurements")
    return world
#enddef


# =============================================================================
# KM15/BH model evaluation at native and common-energy kinematics
# =============================================================================


def evaluate_one_km15(emff, row, ebeam: float) -> Dict[str, float]:
    """Use the already validated dataset-specific phi convention."""
    key = str(row.dataset)
    if key in GEPARD_BMK_DATASETS:
        task = (
            0,
            float(row.xB),
            float(row.Q2),
            float(row.t_abs),
            float(row.phi_bmk_rad),
            float(ebeam),
        )
        return emff.evaluate_km15_point_bmk(task)
    #endif

    if key in DIRECT_PHI_DATASETS:
        task = (
            0,
            float(row.xB),
            float(row.Q2),
            float(row.t_abs),
            float(row.phi_deg),
            float(ebeam),
            "identity",
        )
        return emff.evaluate_km15_point(task)
    #endif

    raise KeyError(f"No KM15 phi convention configured for dataset {key}")
#enddef


def evaluate_km15_world(
        world: pd.DataFrame,
        emff,
        cache_path: Path,
        target_ebeam: float,
        force: bool = False) -> pd.DataFrame:
    """
    Evaluate KM15 total and pure BH at native E and at target E for every point.

    The calculation is intentionally cached because it is deterministic and the
    full world dataset contains several thousand points.
    """
    expected_cols = [
        "point_id",
        "km15_native", "bh_native",
        "km15_target", "bh_target",
    ]

    if cache_path.exists() and not force:
        cache = pd.read_csv(cache_path)
        if (
            len(cache) == len(world)
            and set(expected_cols).issubset(cache.columns)
            and cache["point_id"].astype(str).tolist() == world["point_id"].astype(str).tolist()
        ):
            print(f"[KM15] reusing cache: {cache_path}")
            out = world.copy()
            for col in expected_cols[1:]:
                out[col] = pd.to_numeric(cache[col], errors="coerce").to_numpy(float)
            #endfor
            return finalize_km15_columns(out)
        #endif
    #endif

    rows = []
    total = len(world)
    print(f"[KM15] evaluating native + E={target_ebeam:.3f} GeV for {total:,} points")

    for i, row in enumerate(world.itertuples(index=False), start=1):
        native = evaluate_one_km15(emff, row, float(row.ebeam))
        target = evaluate_one_km15(emff, row, float(target_ebeam))
        rows.append({
            "point_id": str(row.point_id),
            "km15_native": float(native["km15_ep"]),
            "bh_native": float(native["km15_bh"]),
            "km15_target": float(target["km15_ep"]),
            "bh_target": float(target["km15_bh"]),
        })
        if i % 250 == 0 or i == total:
            print(f"[KM15] {i:5d}/{total}")
        #endif
    #endfor

    cache = pd.DataFrame(rows)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache.to_csv(cache_path, index=False)
    print(f"[KM15] cache -> {cache_path}")

    out = world.copy()
    for col in expected_cols[1:]:
        out[col] = cache[col].to_numpy(float)
    #endfor
    return finalize_km15_columns(out)
#enddef


def finalize_km15_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["data_over_km15"] = out["xs"] / out["km15_native"]
    out["data_over_bh"] = out["xs"] / out["bh_native"]
    out["bh_fraction_km15"] = out["bh_native"] / out["km15_native"]
    out["km15_transport_factor"] = out["km15_target"] / out["km15_native"]
    out["bh_transport_factor"] = out["bh_target"] / out["bh_native"]
    out["xs_10p6_km15"] = out["xs"] * out["km15_transport_factor"]
    out["stat_10p6_km15"] = out["stat_abs"] * np.abs(out["km15_transport_factor"])
    out["ptp_sys_10p6_km15"] = out["ptp_sys_abs"] * np.abs(out["km15_transport_factor"])
    out["point_unc_10p6_km15"] = out["point_unc_abs"] * np.abs(out["km15_transport_factor"])
    return out
#enddef


# =============================================================================
# PARTONS/GK16 bridge
# =============================================================================


def write_partons_kinematics(world: pd.DataFrame, path: Path, target_ebeam: float) -> pd.DataFrame:
    """
    Export TWO rows per experimental point for the existing PARTONS driver:
    native E and target E.  This is enough to calculate the exact GK16
    transport factor without changing the PARTONS machinery.
    """
    rows = []
    source_counter: Dict[str, int] = {k: 0 for k in DATASET_ORDER}

    for row in world.itertuples(index=False):
        key = str(row.dataset)
        pkey = PARTONS_DATASET_KEYS[key]

        for state, ebeam, km15_ep, km15_bh in [
            ("native", float(row.ebeam), float(row.km15_native), float(row.bh_native)),
            ("target", float(target_ebeam), float(row.km15_target), float(row.bh_target)),
        ]:
            source_row = source_counter[key]
            source_counter[key] += 1
            rows.append({
                "point_id": f"{row.point_id}:{state}",
                "dataset": pkey,
                "source_row": source_row,
                "xB": float(row.xB),
                "Q2": float(row.Q2),
                "t_abs": float(row.t_abs),
                "phi_deg": float(row.phi_deg),
                "ebeam": float(ebeam),
                "km15_ep": float(km15_ep),
                "km15_bh": float(km15_bh),
                "km15_dvcs": np.nan,
                "km15_int": np.nan,
            })
        #endfor
    #endfor

    table = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    print(f"[PARTONS] exported {len(table):,} native/target kinematic rows -> {path}")
    return table
#enddef


def run_partons_driver(args, kinematics_path: Path, partons_outdir: Path) -> None:
    evaluator = resolve_existing(
        Path(args.partons_evaluator),
        [args.script_dir / "evaluate_bh_model_selection.py"],
    )

    cmd = [
        sys.executable,
        str(evaluator),
        "--kinematics-cache", str(kinematics_path),
        "--outdir", str(partons_outdir),
        "--run-partons",
        "--workers", str(args.partons_workers),
        "--chunk-size", str(args.partons_chunk_size),
    ]
    if args.force_partons:
        cmd.append("--force-partons")
    #endif

    print("[PARTONS] invoking existing validated evaluator:")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(args.script_dir))
#enddef


def merge_gk16_predictions(
        world: pd.DataFrame,
        partons_table_path: Path) -> Tuple[pd.DataFrame, bool]:
    """Merge native and target GK16 predictions by deterministic point_id."""
    if not partons_table_path.exists():
        print("[GK16] completed PARTONS table not present; KM15/BH-only outputs will be written")
        return world.copy(), False
    #endif

    p = pd.read_csv(partons_table_path)
    required = {"point_id", "partons_ep_gk16"}
    missing = sorted(required.difference(p.columns))
    if missing:
        raise KeyError(
            f"GK16 table {partons_table_path} is missing columns: {missing}"
        )
    #endif

    pred = dict(zip(p["point_id"].astype(str), pd.to_numeric(p["partons_ep_gk16"], errors="coerce")))
    out = world.copy()
    out["gk16_native"] = [pred.get(f"{pid}:native", np.nan) for pid in out["point_id"]]
    out["gk16_target"] = [pred.get(f"{pid}:target", np.nan) for pid in out["point_id"]]

    good = finite_positive(out["gk16_native"]) & finite_positive(out["gk16_target"])
    if not np.all(good):
        warnings.warn(
            f"GK16 predictions missing/nonpositive for {int((~good).sum())} of {len(out)} points"
        )
    #endif

    out["gk16_transport_factor"] = out["gk16_target"] / out["gk16_native"]
    out["xs_10p6_gk16"] = out["xs"] * out["gk16_transport_factor"]

    # Nominal common-energy value remains KM15.  The full model excursion from
    # KM15 to GK16 is retained as an explicit transport uncertainty.
    out["transport_model_unc_abs"] = np.abs(out["xs_10p6_gk16"] - out["xs_10p6_km15"])
    out["transport_model_unc_frac"] = out["transport_model_unc_abs"] / np.abs(out["xs_10p6_km15"])
    out["data_over_gk16"] = out["xs"] / out["gk16_native"]

    print(f"[GK16] merged completed predictions from {partons_table_path}")
    return out, True
#enddef


# =============================================================================
# Dataset summaries and normalization-aware model scores
# =============================================================================


def make_dataset_summary(world: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key in DATASET_ORDER:
        d = world.loc[world["dataset"] == key]
        if d.empty:
            continue
        #endif
        rows.append({
            "dataset": key,
            "dataset_label": DATASET_LABELS[key],
            "N": len(d),
            "ebeam_min_GeV": d["ebeam"].min(),
            "ebeam_max_GeV": d["ebeam"].max(),
            "xB_min": d["xB"].min(),
            "xB_max": d["xB"].max(),
            "Q2_min_GeV2": d["Q2"].min(),
            "Q2_max_GeV2": d["Q2"].max(),
            "t_abs_min_GeV2": d["t_abs"].min(),
            "t_abs_max_GeV2": d["t_abs"].max(),
            "median_stat_pct": 100.0 * np.nanmedian(d["stat_frac"]),
            "median_ptp_sys_pct": 100.0 * np.nanmedian(d["ptp_sys_frac"]),
            "median_point_unc_pct": 100.0 * np.nanmedian(d["point_unc_frac"]),
            "correlated_norm_pct": 100.0 * float(d["norm_frac"].iloc[0]),
        })
    #endfor
    return pd.DataFrame(rows)
#enddef


def fit_one_normalization_nuisance(
        data: np.ndarray,
        model: np.ndarray,
        sigma: np.ndarray,
        norm_frac: float) -> Dict[str, float]:
    """
    Fit one correlated multiplicative scale nuisance beta.

    Convention:
        shifted_data = data * (1 + beta * norm_frac)

    and minimize
        sum[(shifted_data-model)^2/sigma^2] + beta^2.

    beta therefore has units of published normalization standard deviations.
    """
    y = np.asarray(data, dtype=float)
    m = np.asarray(model, dtype=float)
    s = np.asarray(sigma, dtype=float)
    good = np.isfinite(y) & np.isfinite(m) & np.isfinite(s) & (s > 0.0) & (y > 0.0) & (m > 0.0)
    y, m, s = y[good], m[good], s[good]

    if len(y) == 0:
        return {"N": 0, "beta": np.nan, "scale": np.nan, "chi2": np.nan, "chi2_per_point": np.nan}
    #endif

    if norm_frac <= 0.0:
        beta = 0.0
    else:
        x = norm_frac * y
        w = 1.0 / (s * s)
        beta = -np.sum(w * x * (y - m)) / (1.0 + np.sum(w * x * x))
    #endif

    scale = 1.0 + beta * norm_frac
    residual = scale * y - m
    chi2_points = np.sum((residual / s)**2)
    chi2 = chi2_points + beta**2

    return {
        "N": len(y),
        "beta": float(beta),
        "scale": float(scale),
        "chi2": float(chi2),
        "chi2_per_point": float(chi2 / len(y)),
        "median_abs_fractional_residual": float(np.nanmedian(np.abs(residual / m))),
        "rms_pull": float(np.sqrt(np.nanmean((residual / s)**2))),
    }
#enddef


def make_native_model_scores(world: pd.DataFrame, have_gk16: bool) -> pd.DataFrame:
    models = [("km15", "km15_native"), ("bh", "bh_native")]
    if have_gk16:
        models.insert(1, ("gk16", "gk16_native"))
    #endif

    rows = []
    for key in DATASET_ORDER:
        d = world.loc[world["dataset"] == key]
        if d.empty:
            continue
        #endif
        norm = float(d["norm_frac"].iloc[0])
        for model_name, col in models:
            result = fit_one_normalization_nuisance(
                d["xs"].to_numpy(float),
                d[col].to_numpy(float),
                d["point_unc_abs"].to_numpy(float),
                norm,
            )
            result.update({
                "dataset": key,
                "dataset_label": DATASET_LABELS[key],
                "model": model_name,
                "normalization_fraction": norm,
            })
            rows.append(result)
        #endfor
    #endfor
    return pd.DataFrame(rows)
#enddef


# =============================================================================
# Pairwise kinematic matching and direct dataset consistency
# =============================================================================


def candidate_matches(a: pd.DataFrame, b: pd.DataFrame, cfg: MatchConfig) -> pd.DataFrame:
    """Construct all candidate A/B pairs inside explicit 4D matching windows."""
    rows = []

    # The datasets are only a few thousand points; grouping by coarse xB first
    # avoids a full O(N^2) Cartesian product while retaining transparent cuts.
    for ia, ra in a.iterrows():
        xbmask = np.abs(b["xB"].to_numpy(float) - float(ra["xB"])) <= cfg.dxb
        if not np.any(xbmask):
            continue
        #endif

        sub = b.loc[xbmask]
        dx = np.abs(sub["xB"].to_numpy(float) - float(ra["xB"]))
        dq = np.abs(sub["Q2"].to_numpy(float) - float(ra["Q2"]))
        dt = np.abs(sub["t_abs"].to_numpy(float) - float(ra["t_abs"]))
        dp = circular_phi_difference_deg(sub["phi_deg"].to_numpy(float), float(ra["phi_deg"]))

        keep = (dq <= cfg.dq2) & (dt <= cfg.dt) & (dp <= cfg.dphi)
        if not np.any(keep):
            continue
        #endif

        kept = sub.loc[keep]
        dx = dx[keep]
        dq = dq[keep]
        dt = dt[keep]
        dp = dp[keep]

        score = np.sqrt(
            (dx / cfg.dxb)**2
            + (dq / cfg.dq2)**2
            + (dt / cfg.dt)**2
            + (dp / cfg.dphi)**2
        )

        for j, (ib, rb) in enumerate(kept.iterrows()):
            rows.append({
                "ia": int(ia),
                "ib": int(ib),
                "point_id_a": str(ra["point_id"]),
                "point_id_b": str(rb["point_id"]),
                "dxB": float(dx[j]),
                "dQ2": float(dq[j]),
                "dt_abs": float(dt[j]),
                "dphi_deg": float(dp[j]),
                "match_score": float(score[j]),
            })
        #endfor
    #endfor

    return pd.DataFrame(rows)
#enddef


def greedy_one_to_one_matches(candidates: pd.DataFrame) -> pd.DataFrame:
    """Choose the closest nonconflicting candidate pairs deterministically."""
    if candidates.empty:
        return candidates.copy()
    #endif

    c = candidates.sort_values(
        ["match_score", "dxB", "dQ2", "dt_abs", "dphi_deg", "point_id_a", "point_id_b"]
    )
    used_a = set()
    used_b = set()
    selected = []
    for row in c.itertuples(index=False):
        if row.point_id_a in used_a or row.point_id_b in used_b:
            continue
        #endif
        used_a.add(row.point_id_a)
        used_b.add(row.point_id_b)
        selected.append(row._asdict())
    #endfor
    return pd.DataFrame(selected)
#enddef


def fit_pair_relative_normalization(
        a_value: np.ndarray,
        b_value: np.ndarray,
        sigma_a: np.ndarray,
        sigma_b: np.ndarray,
        transport_unc: np.ndarray,
        norm_a: float,
        norm_b: float) -> Dict[str, float]:
    """
    Fit the *relative* normalization of two matched datasets in log-ratio space.

    Only the relative normalization between two independent experiments is
    identifiable from an A-vs-B comparison.  Fitting independent absolute
    scales for A and B while keeping their point uncertainties fixed creates an
    unphysical common-mode direction in which both datasets can be scaled
    toward zero.  This routine removes that degeneracy.

    For each matched point, define

        z_i = ln(A_i / B_i)

    with approximate point uncertainty

        s_i^2 = (sigma_A/A)^2 + (sigma_B/B)^2
                + (sigma_transport/A)^2.

    A relative multiplicative shift exp(eta) is applied to A, so the residual
    becomes z_i + eta.  The independent quoted normalization uncertainties are
    combined in quadrature in log space and used as the prior width on eta.

    The returned scale_a_to_b = exp(eta) is therefore the factor multiplying
    dataset A (after local kinematic transport) to align it with dataset B.
    """
    A = np.asarray(a_value, dtype=float)
    B = np.asarray(b_value, dtype=float)
    SA = np.asarray(sigma_a, dtype=float)
    SB = np.asarray(sigma_b, dtype=float)
    ST = np.asarray(transport_unc, dtype=float)

    good = (
        np.isfinite(A) & np.isfinite(B) & np.isfinite(SA) & np.isfinite(SB)
        & np.isfinite(ST) & (A > 0.0) & (B > 0.0)
        & (SA >= 0.0) & (SB >= 0.0) & (ST >= 0.0)
    )
    A, B, SA, SB, ST = A[good], B[good], SA[good], SB[good], ST[good]
    if len(A) == 0:
        return {
            "N": 0,
            "eta": np.nan,
            "beta_relative": np.nan,
            "scale_a_to_b": np.nan,
            "relative_norm_frac": np.nan,
            "chi2": np.nan,
            "chi2_per_point": np.nan,
        }
    #endif

    z = np.log(A / B)
    sfrac = np.sqrt((SA / A)**2 + (SB / B)**2 + (ST / A)**2)
    finite = np.isfinite(z) & np.isfinite(sfrac) & (sfrac > 0.0)
    z, sfrac = z[finite], sfrac[finite]
    if len(z) == 0:
        return {
            "N": 0,
            "eta": np.nan,
            "beta_relative": np.nan,
            "scale_a_to_b": np.nan,
            "relative_norm_frac": np.nan,
            "chi2": np.nan,
            "chi2_per_point": np.nan,
        }
    #endif

    # Independent multiplicative normalization uncertainties combine in
    # quadrature.  log1p maps the quoted fractional scale width onto the same
    # additive variable eta used by ln(A/B).
    sigma_eta_a = math.log1p(max(float(norm_a), 0.0))
    sigma_eta_b = math.log1p(max(float(norm_b), 0.0))
    sigma_eta = math.sqrt(sigma_eta_a**2 + sigma_eta_b**2)

    w = 1.0 / (sfrac * sfrac)
    if sigma_eta > 0.0:
        eta = -float(np.sum(w * z)) / float(np.sum(w) + 1.0 / sigma_eta**2)
        penalty = (eta / sigma_eta)**2
        beta_relative = eta / sigma_eta
    else:
        eta = 0.0
        penalty = 0.0
        beta_relative = 0.0
    #endif

    pulls = (z + eta) / sfrac
    chi2 = float(np.sum(pulls**2) + penalty)
    return {
        "N": int(len(z)),
        "eta": float(eta),
        "beta_relative": float(beta_relative),
        "scale_a_to_b": float(math.exp(eta)),
        "relative_norm_frac": float(math.sqrt(norm_a**2 + norm_b**2)),
        "chi2": chi2,
        "chi2_per_point": float(chi2 / len(z)),
    }
#enddef

def build_pairwise_comparisons(
        world: pd.DataFrame,
        cfg: MatchConfig,
        have_gk16: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Match every dataset pair and transport A to B using local model ratios.

    KM15 is nominal.  If GK16 is available, |A_to_B(KM15)-A_to_B(GK16)| is
    added as a separate transport-model uncertainty and included in the pull
    denominator for the summary consistency score.
    """
    all_matches = []
    summaries = []

    for i, key_a in enumerate(DATASET_ORDER):
        for key_b in DATASET_ORDER[i + 1:]:
            a = world.loc[world["dataset"] == key_a].copy().reset_index(drop=True)
            b = world.loc[world["dataset"] == key_b].copy().reset_index(drop=True)
            if a.empty or b.empty:
                continue
            #endif

            candidates = candidate_matches(a, b, cfg)
            selected = greedy_one_to_one_matches(candidates)
            if selected.empty:
                summaries.append({
                    "dataset_a": key_a,
                    "dataset_b": key_b,
                    "N_candidates": len(candidates),
                    "N_matches": 0,
                })
                continue
            #endif

            rows = []
            for m in selected.itertuples(index=False):
                ra = a.iloc[int(m.ia)]
                rb = b.iloc[int(m.ib)]

                c_km15 = float(rb["km15_native"] / ra["km15_native"])
                a_to_b_km15 = float(ra["xs"] * c_km15)
                a_stat = float(ra["stat_abs"] * abs(c_km15))
                a_sys = float(ra["ptp_sys_abs"] * abs(c_km15))
                a_point = float(ra["point_unc_abs"] * abs(c_km15))

                c_gk16 = np.nan
                a_to_b_gk16 = np.nan
                # KM15 is the nominal and sole transport prescription in the
                # current world-data study, so no additional transport-model
                # systematic is assigned here.
                transport_unc = 0.0
                if have_gk16 and np.isfinite(ra["gk16_native"]) and np.isfinite(rb["gk16_native"]):
                    c_gk16 = float(rb["gk16_native"] / ra["gk16_native"])
                    a_to_b_gk16 = float(ra["xs"] * c_gk16)
                    transport_unc = abs(a_to_b_gk16 - a_to_b_km15)
                #endif

                sigma_compare = math.sqrt(
                    a_point**2 + float(rb["point_unc_abs"])**2 + transport_unc**2
                )
                raw_pull = (
                    (a_to_b_km15 - float(rb["xs"])) / sigma_compare
                    if sigma_compare > 0.0 else np.nan
                )

                rows.append({
                    "dataset_a": key_a,
                    "dataset_b": key_b,
                    "point_id_a": str(ra["point_id"]),
                    "point_id_b": str(rb["point_id"]),
                    "match_score": float(m.match_score),
                    "dxB": float(m.dxB),
                    "dQ2": float(m.dQ2),
                    "dt_abs": float(m.dt_abs),
                    "dphi_deg": float(m.dphi_deg),
                    "xB_a": float(ra["xB"]), "xB_b": float(rb["xB"]),
                    "Q2_a": float(ra["Q2"]), "Q2_b": float(rb["Q2"]),
                    "t_abs_a": float(ra["t_abs"]), "t_abs_b": float(rb["t_abs"]),
                    "phi_a": float(ra["phi_deg"]), "phi_b": float(rb["phi_deg"]),
                    "ebeam_a": float(ra["ebeam"]), "ebeam_b": float(rb["ebeam"]),
                    "xs_a": float(ra["xs"]), "xs_b": float(rb["xs"]),
                    "point_unc_a": float(ra["point_unc_abs"]),
                    "point_unc_b": float(rb["point_unc_abs"]),
                    "km15_local_transport_factor": c_km15,
                    "xs_a_to_b_km15": a_to_b_km15,
                    "stat_a_to_b_km15": a_stat,
                    "ptp_sys_a_to_b_km15": a_sys,
                    "point_unc_a_to_b_km15": a_point,
                    "gk16_local_transport_factor": c_gk16,
                    "xs_a_to_b_gk16": a_to_b_gk16,
                    "transport_model_unc_abs": transport_unc,
                    "comparison_unc_abs": sigma_compare,
                    "raw_pull": raw_pull,
                    "norm_frac_a": float(ra["norm_frac"]),
                    "norm_frac_b": float(rb["norm_frac"]),
                })
            #endfor

            pair = pd.DataFrame(rows)
            relfit = fit_pair_relative_normalization(
                pair["xs_a_to_b_km15"].to_numpy(float),
                pair["xs_b"].to_numpy(float),
                pair["point_unc_a_to_b_km15"].to_numpy(float),
                pair["point_unc_b"].to_numpy(float),
                pair["transport_model_unc_abs"].to_numpy(float),
                float(pair["norm_frac_a"].iloc[0]),
                float(pair["norm_frac_b"].iloc[0]),
            )
            relative_scale = float(relfit["scale_a_to_b"])
            pair["relative_norm_beta"] = float(relfit["beta_relative"])
            pair["relative_scale_a_to_b"] = relative_scale
            pair["relative_norm_frac"] = float(relfit["relative_norm_frac"])

            pair["raw_log_pull"] = np.log(pair["xs_a_to_b_km15"] / pair["xs_b"]) / np.sqrt(
                (pair["point_unc_a_to_b_km15"] / pair["xs_a_to_b_km15"])**2
                + (pair["point_unc_b"] / pair["xs_b"])**2
                + (pair["transport_model_unc_abs"] / pair["xs_a_to_b_km15"])**2
            )
            pair["profiled_pull"] = (
                np.log(pair["xs_a_to_b_km15"] / pair["xs_b"]) + math.log(relative_scale)
            ) / np.sqrt(
                (pair["point_unc_a_to_b_km15"] / pair["xs_a_to_b_km15"])**2
                + (pair["point_unc_b"] / pair["xs_b"])**2
                + (pair["transport_model_unc_abs"] / pair["xs_a_to_b_km15"])**2
            )

            summaries.append({
                "dataset_a": key_a,
                "dataset_b": key_b,
                "N_candidates": len(candidates),
                "N_matches": len(pair),
                "relative_norm_beta": float(relfit["beta_relative"]),
                "relative_scale_a_to_b": relative_scale,
                "relative_norm_pct": 100.0 * float(relfit["relative_norm_frac"]),
                "chi2_with_norm_penalty": float(relfit["chi2"]),
                "chi2_per_match": float(relfit["chi2_per_point"]),
                "raw_pull_rms": float(np.sqrt(np.nanmean(pair["raw_pull"]**2))),
                "raw_log_pull_rms": float(np.sqrt(np.nanmean(pair["raw_log_pull"]**2))),
                "profiled_pull_rms": float(np.sqrt(np.nanmean(pair["profiled_pull"]**2))),
                "median_abs_profiled_pull": float(np.nanmedian(np.abs(pair["profiled_pull"]))),
                "median_transport_model_unc_pct": (
                    100.0 * float(np.nanmedian(pair["transport_model_unc_abs"] / np.abs(pair["xs_a_to_b_km15"])))
                ),
            })
            all_matches.append(pair)
        #endfor
    #endfor

    match_table = pd.concat(all_matches, ignore_index=True) if all_matches else pd.DataFrame()
    return match_table, pd.DataFrame(summaries)
#enddef


# =============================================================================
# Plotting
# =============================================================================


# Fixed display windows keep a handful of pathological/very-low-weight points
# from determining the visual scale of summary figures.  The underlying points
# remain in all tables and quantitative calculations; matplotlib simply clips
# markers/error bars outside these display ranges.
NATIVE_RATIO_YLIMS = {
    "KM15": (0.0, 3.0),
    "GK16": (0.0, 3.0),
    "BH": (0.0, 10.0),
}
PAIRWISE_PULL_YLIM = (-5.5, 5.5)


def _annotate_clipped_y(ax, values: pd.Series, ylow: float, yhigh: float) -> None:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(float)
    finite = np.isfinite(arr)
    nlow = int(np.sum(finite & (arr < ylow)))
    nhigh = int(np.sum(finite & (arr > yhigh)))
    nclip = nlow + nhigh
    if nclip <= 0:
        return
    #endif
    ax.text(
        0.012, 0.018,
        f"{nclip} central point(s) outside display range; all retained in analysis",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, alpha=0.75,
    )
#enddef


def plot_native_model_ratios(world: pd.DataFrame, outdir: Path, have_gk16: bool) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    variables = [
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
        ("phi_deg", r"$\phi$ (deg)"),
    ]

    model_specs = [("data_over_km15", "KM15"), ("data_over_bh", "BH")]
    if have_gk16:
        model_specs.insert(1, ("data_over_gk16", "GK16"))
    #endif

    for ratio_col, model_label in model_specs:
        for variable, xlabel in variables:
            fig, ax = plt.subplots(figsize=(8.4, 5.6))
            for key in DATASET_ORDER:
                d = world.loc[world["dataset"] == key]
                if d.empty:
                    continue
                #endif
                relerr = d["point_unc_abs"] / d[
                    "km15_native" if model_label == "KM15" else (
                        "gk16_native" if model_label == "GK16" else "bh_native"
                    )
                ]
                ax.errorbar(
                    d[variable], d[ratio_col], yerr=relerr,
                    fmt="o", ms=2.8, lw=0.7, capsize=0,
                    alpha=0.60, label=DATASET_LABELS[key],
                )
            #endfor
            ax.axhline(1.0, lw=1.1, linestyle="--")
            ylow, yhigh = NATIVE_RATIO_YLIMS[model_label]
            ax.set_ylim(ylow, yhigh)
            _annotate_clipped_y(ax, world[ratio_col], ylow, yhigh)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(rf"$\sigma_{{\rm data}}/\sigma_{{\rm {model_label}}}$")
            ax.set_title(f"Published world data / {model_label} at native kinematics")
            ax.grid(alpha=0.20)
            ax.legend(fontsize=8, ncol=2, frameon=False)
            fig.tight_layout()
            tag = variable.replace("_abs", "")
            fig.savefig(outdir / f"native_data_over_{model_label.lower()}_vs_{tag}.png", dpi=220)
            plt.close(fig)
        #endfor
    #endfor
#enddef


def plot_transport_uncertainty(world: pd.DataFrame, outdir: Path, have_gk16: bool) -> None:
    if not have_gk16:
        return
    #endif

    outdir.mkdir(parents=True, exist_ok=True)
    for variable, xlabel in [
        ("ebeam", r"Native $E_{\rm beam}$ (GeV)"),
        ("xB", r"$x_B$"),
        ("Q2", r"$Q^2$ (GeV$^2$)"),
        ("t_abs", r"$|t|$ (GeV$^2$)"),
    ]:
        fig, ax = plt.subplots(figsize=(8.4, 5.6))
        for key in DATASET_ORDER:
            d = world.loc[world["dataset"] == key]
            ax.scatter(
                d[variable], 100.0 * d["transport_model_unc_frac"],
                s=14, alpha=0.55, label=DATASET_LABELS[key],
            )
        #endfor
        ax.set_xlabel(xlabel)
        ax.set_ylabel("KM15–GK16 transport excursion (%)")
        ax.set_title(rf"Model dependence of transport to $E_{{\rm beam}}={TARGET_EBEAM_GEV:.1f}$ GeV")
        ax.grid(alpha=0.20)
        ax.legend(fontsize=8, ncol=2, frameon=False)
        fig.tight_layout()
        tag = variable.replace("_abs", "")
        fig.savefig(outdir / f"transport_model_uncertainty_vs_{tag}.png", dpi=220)
        plt.close(fig)
    #endfor
#enddef


def plot_pairwise_matrix(summary: pd.DataFrame, outdir: Path) -> None:
    if summary.empty:
        return
    #endif
    outdir.mkdir(parents=True, exist_ok=True)

    n = len(DATASET_ORDER)
    matrix = np.full((n, n), np.nan)
    counts = np.zeros((n, n), dtype=int)
    index = {k: i for i, k in enumerate(DATASET_ORDER)}

    for row in summary.itertuples(index=False):
        i = index[row.dataset_a]
        j = index[row.dataset_b]
        if int(row.N_matches) > 0 and np.isfinite(getattr(row, "chi2_per_match", np.nan)):
            matrix[i, j] = float(row.chi2_per_match)
            matrix[j, i] = float(row.chi2_per_match)
            counts[i, j] = counts[j, i] = int(row.N_matches)
        #endif
    #endfor

    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    im = ax.imshow(matrix, origin="upper", aspect="equal")
    labels = [DATASET_LABELS[k] for k in DATASET_ORDER]
    ax.set_xticks(np.arange(n), labels=labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(n), labels=labels)
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.text(j, i, "—", ha="center", va="center")
            elif np.isfinite(matrix[i, j]):
                ax.text(j, i, f"{matrix[i,j]:.2f}\nN={counts[i,j]}", ha="center", va="center", fontsize=8)
            #endif
        #endfor
    #endfor
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\chi^2/N_{\rm match}$ after normalization shifts")
    ax.set_title("Pairwise world-data compatibility at matched kinematics")
    fig.tight_layout()
    fig.savefig(outdir / "pairwise_compatibility_matrix.png", dpi=220)
    plt.close(fig)
#enddef


def plot_pairwise_pulls(matches: pd.DataFrame, outdir: Path) -> None:
    if matches.empty:
        return
    #endif
    outdir.mkdir(parents=True, exist_ok=True)

    for variable, xlabel, col in [
        ("t", r"$|t|_B$ (GeV$^2$)", "t_abs_b"),
        ("xB", r"$x_{B,B}$", "xB_b"),
        ("Q2", r"$Q_B^2$ (GeV$^2$)", "Q2_b"),
        ("phi", r"$\phi_B$ (deg)", "phi_b"),
    ]:
        fig, ax = plt.subplots(figsize=(8.4, 5.6))
        for (ka, kb), d in matches.groupby(["dataset_a", "dataset_b"], sort=False):
            ax.scatter(
                d[col], d["profiled_pull"],
                s=15, alpha=0.50,
                label=f"{DATASET_LABELS[ka]} → {DATASET_LABELS[kb]}",
            )
        #endfor
        ax.axhline(0.0, lw=1.0)
        ax.axhline(+1.0, lw=0.8, linestyle="--")
        ax.axhline(-1.0, lw=0.8, linestyle="--")
        ax.set_ylim(*PAIRWISE_PULL_YLIM)
        _annotate_clipped_y(
            ax, matches["profiled_pull"], PAIRWISE_PULL_YLIM[0], PAIRWISE_PULL_YLIM[1]
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Matched-data residual / pointwise uncertainty")
        ax.set_title("Pairwise matched world-data residuals")
        ax.grid(alpha=0.20)
        ax.legend(fontsize=6.5, ncol=2, frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / f"pairwise_profiled_pulls_vs_{variable}.png", dpi=220)
        plt.close(fig)
    #endfor
#enddef


# =============================================================================
# Output bookkeeping
# =============================================================================


def save_outputs(
        world: pd.DataFrame,
        dataset_summary: pd.DataFrame,
        model_scores: pd.DataFrame,
        matches: pd.DataFrame,
        pair_summary: pd.DataFrame,
        outdir: Path,
        have_gk16: bool) -> None:
    tables = outdir / "tables"
    figures = outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    world.to_csv(tables / "canonical_world_data_with_models.csv", index=False)
    dataset_summary.to_csv(tables / "dataset_summary.csv", index=False)
    model_scores.to_csv(tables / "native_model_scores.csv", index=False)
    pair_summary.to_csv(tables / "pairwise_overlap_summary.csv", index=False)
    if not matches.empty:
        matches.to_csv(tables / "pairwise_matched_points.csv", index=False)
    #endif

    transport_cols = [
        "point_id", "dataset", "dataset_label",
        "ebeam", "xB", "Q2", "t_abs", "phi_deg",
        "xs", "stat_abs", "ptp_sys_abs", "norm_frac",
        "km15_transport_factor", "xs_10p6_km15",
        "stat_10p6_km15", "ptp_sys_10p6_km15", "point_unc_10p6_km15",
    ]
    if have_gk16:
        transport_cols += [
            "gk16_transport_factor", "xs_10p6_gk16",
            "transport_model_unc_abs", "transport_model_unc_frac",
        ]
    #endif
    world[transport_cols].to_csv(tables / "world_data_transported_to_10p6.csv", index=False)

    plot_native_model_ratios(world, figures / "native_model_ratios", have_gk16)
    plot_transport_uncertainty(world, figures / "transport", have_gk16)
    plot_pairwise_matrix(pair_summary, figures / "pairwise")
    plot_pairwise_pulls(matches, figures / "pairwise")
#enddef


def print_summary(dataset_summary: pd.DataFrame, model_scores: pd.DataFrame, pair_summary: pd.DataFrame, have_gk16: bool) -> None:
    print("\n" + "=" * 80)
    print("DATASET SUMMARY")
    print("=" * 80)
    cols = ["dataset_label", "N", "ebeam_min_GeV", "ebeam_max_GeV", "median_point_unc_pct", "correlated_norm_pct"]
    print(dataset_summary[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\n" + "=" * 80)
    print("NATIVE-KINEMATICS MODEL SCORES")
    print("=" * 80)
    print(model_scores[[
        "dataset_label", "model", "N", "beta", "scale", "chi2_per_point", "median_abs_fractional_residual"
    ]].to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print("\n" + "=" * 80)
    print("PAIRWISE MATCHED-DATA SUMMARY")
    print("=" * 80)
    if pair_summary.empty:
        print("No pairwise matches were found with the configured windows.")
    else:
        print(pair_summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    #endif

    if not have_gk16:
        print("\n[NOTE] World-data transport uses KM15 only by design.")
        print("       No PARTONS/GK16 files are generated by this script.")
    #endif
#enddef


# =============================================================================
# CLI
# =============================================================================


def build_parser() -> argparse.ArgumentParser:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(
        description="World proton DVCS/BH cross-section consistency and model comparison"
    )
    p.add_argument("--outdir", default="../output/world_dvcs_cross_section_comparison")
    p.add_argument("--emff-script", default=str(here / "extract_emff_from_dvcs_bh.py"))
    p.add_argument("--saylor-file", default=str(here / "import" / "saylor_CLAS6.txt"))
    p.add_argument("--georges-file", default=str(here / "import" / "E12-06-114.xlsx"))
    p.add_argument("--lee-file", default=str(here / "import" / "clasdb_E214M1.txt"))
    p.add_argument("--target-ebeam", type=float, default=TARGET_EBEAM_GEV)
    p.add_argument("--workers", type=int, default=1, help="Reserved for future KM15 multiprocessing; current first pass evaluates serially for model safety")
    p.add_argument("--force-km15", action="store_true")

    p.add_argument("--match-dxb", type=float, default=DEFAULT_MATCH_DXB)
    p.add_argument("--match-dq2", type=float, default=DEFAULT_MATCH_DQ2)
    p.add_argument("--match-dt", type=float, default=DEFAULT_MATCH_DT)
    p.add_argument("--match-dphi", type=float, default=DEFAULT_MATCH_DPHI)
    return p
#enddef


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    args.script_dir = Path(__file__).resolve().parent

    outdir = Path(args.outdir).expanduser()
    if not outdir.is_absolute():
        outdir = (args.script_dir / outdir).resolve()
    #endif
    outdir.mkdir(parents=True, exist_ok=True)

    emff_path = resolve_existing(
        Path(args.emff_script),
        [args.script_dir / "extract_emff_from_dvcs_bh.py"],
    )
    emff = load_python_module(emff_path, "emff_world_compare_backend")

    # ---------------------------------------------------------------------
    # 1. Load and canonicalize the six published measurements.
    # ---------------------------------------------------------------------
    world = load_world_data(args, emff)

    # ---------------------------------------------------------------------
    # 2. Native + common-energy KM15/BH calculations.
    # ---------------------------------------------------------------------
    km15_cache = outdir / "cache" / "km15_native_and_target.csv"
    world = evaluate_km15_world(
        world,
        emff,
        cache_path=km15_cache,
        target_ebeam=float(args.target_ebeam),
        force=args.force_km15,
    )

    # ---------------------------------------------------------------------
    # 3. KM15-only transport prescription.
    #
    # The world-data consistency study intentionally uses KM15 alone for
    # model-assisted transport.  PARTONS/GK16 is not run here: the principal
    # external cross-checks are the near-10-GeV Hall-A settings and CLAS12
    # pass-1, while lower-energy measurements are transported with KM15.
    # ---------------------------------------------------------------------
    have_gk16 = False
    print("[TRANSPORT] KM15-only prescription; no PARTONS/GK16 calculation requested")

    # ---------------------------------------------------------------------
    # 4. QA summaries and native-kinematics model comparisons.
    # ---------------------------------------------------------------------
    dataset_summary = make_dataset_summary(world)
    model_scores = make_native_model_scores(world, have_gk16)

    # ---------------------------------------------------------------------
    # 5. Pairwise experiment matching and local model-assisted transport.
    # ---------------------------------------------------------------------
    match_cfg = MatchConfig(
        dxb=float(args.match_dxb),
        dq2=float(args.match_dq2),
        dt=float(args.match_dt),
        dphi=float(args.match_dphi),
    )
    matches, pair_summary = build_pairwise_comparisons(world, match_cfg, have_gk16)

    # ---------------------------------------------------------------------
    # 6. Tables and note-quality first-pass figures.
    # ---------------------------------------------------------------------
    save_outputs(
        world,
        dataset_summary,
        model_scores,
        matches,
        pair_summary,
        outdir,
        have_gk16,
    )

    print_summary(dataset_summary, model_scores, pair_summary, have_gk16)
    print(f"\n[OUTPUT] {outdir}")
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif
