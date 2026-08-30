#!/usr/bin/env python3
"""
Prepare and diagnose KM15 / VGG99 / GK16 BH-purity model selection.

Purpose
-------
This script is deliberately separate from the radius bias/variance study.

Workflow:
  1. Export the exact experimental kinematic points used by the EMFF analysis.
  2. Keep the already-computed KM15/BH decomposition attached to each point.
  3. Generate PARTONS XML scenarios for VGG99 and GK16 at exactly those points.
  4. After PARTONS results have been converted to the simple result CSV schema
     documented below, merge all three models and construct consensus BH-purity
     selections for thresholds from 1% to 10%.
  5. Produce model-agreement diagnostics and a compact selection CSV that the
     EMFF fitter can consume without rerunning any GPD model.

The model calculations are therefore cached once. Threshold scans and later
form-factor fits are cheap.

PARTONS choice
--------------
VGG:
  GPDVGG99 + DVCSCFFStandard(LO) + DVCSProcessVGG99.

GK:
  GPDGK16 + DVCSCFFStandard(LO) + DVCSProcessGV08.

GK16 is used rather than GK11 because PARTONS documents GK16 as the corrected
successor to GK11. GK19 changes the pseudoscalar/transversity sector while H
and E remain those of GK16, so it does not add useful independent information
for this unpolarized DVCS purity test.

IMPORTANT: PARTONS documents GPDGK16 as not thread-safe. Do not evaluate GK16
with threads. Separate processes are acceptable if each process owns an
independent PARTONS instance.

Expected PARTONS result CSV schema
----------------------------------
One CSV per model:
    point_id,ep_xs
where ep_xs is the unpolarized e p -> e p gamma cross section in the SAME
normalization/units as the BH cross section stored in the kinematic table.

The script intentionally validates normalization before accepting results.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd

DEFAULT_MAIN = "extract_emff_from_dvcs_bh.py"
DEFAULT_OUTDIR = "output/emff_from_bh_paper_method/bh_model_selection"
DEFAULT_THRESHOLDS = [0.01 * i for i in range(1, 11)]

# Dataset labels used consistently in the model-selection cache.
DATASET_JO = "jo2015"
DATASET_SAYLOR = "saylor2018"
DATASET_HALLA = "halla_defurne2015"
DATASET_PASS1 = "pass1"


def locate_partons() -> Dict[str, object]:
    """Inspect PATH/environment for a usable PARTONS installation."""
    names = [
        "partons-example",
        "partons",
        "PARTONS",
    ]
    found = {name: shutil.which(name) for name in names}
    env_keys = [
        "PARTONS",
        "PARTONS_ROOT",
        "PARTONS_HOME",
        "LD_LIBRARY_PATH",
        "CMAKE_PREFIX_PATH",
        "PKG_CONFIG_PATH",
    ]
    env = {key: os.environ.get(key, "") for key in env_keys}

    # Lightweight filesystem hints only; avoid expensive recursive scans.
    hints = []
    for base in [
        Path("/cvmfs"),
        Path("/u/scigroup/cvmfs"),
        Path("/usr/local"),
        Path.home(),
    ]:
        if base.exists():
            for candidate in [
                base / "partons",
                base / "PARTONS",
                base / "bin" / "partons-example",
            ]:
                if candidate.exists():
                    hints.append(str(candidate))
                #endif
            #endfor
        #endif
    #endfor

    return {"executables": found, "environment": env, "hints": hints}
#enddef


def print_partons_check() -> None:
    info = locate_partons()
    print("[PARTONS check]")
    for name, path in info["executables"].items():
        print(f"  {name:16s}: {path or 'not found on PATH'}")
    #endfor
    for key, value in info["environment"].items():
        if value:
            print(f"  env {key}: {value}")
        #endif
    #endfor
    if info["hints"]:
        print("  filesystem hints:")
        for item in info["hints"]:
            print(f"    {item}")
        #endfor
    #endif
#enddef


def import_main_analysis(path: Path):
    """Import the production EMFF script without copying its data loaders."""
    import importlib.util

    path = path.resolve()
    spec = importlib.util.spec_from_file_location("emff_main_for_model_selection", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import analysis script: {path}")
    #endif
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
#enddef


def _first_existing(df: pd.DataFrame, names: Iterable[str]) -> str:
    for name in names:
        if name in df.columns:
            return name
        #endif
    #endfor
    raise KeyError(f"None of these columns exist: {list(names)}")
#enddef


def canonicalize_points(
        df: pd.DataFrame,
        dataset: str,
        ebeam_default: float | None = None) -> pd.DataFrame:
    """Convert an analysis dataframe to the common model-evaluation schema."""
    xcol = _first_existing(df, ["xB", "xBavg", "x"])
    qcol = _first_existing(df, ["Q2", "Q2avg"])
    tcol = _first_existing(df, ["t_abs", "t_abs_avg", "-t", "t"])
    pcol = _first_existing(df, ["phi_deg", "phi", "phiavg"])

    if "ebeam" in df.columns:
        ebeam = df["ebeam"].to_numpy(float)
    elif "beam_energy" in df.columns:
        ebeam = df["beam_energy"].to_numpy(float)
    elif ebeam_default is not None:
        ebeam = np.full(len(df), float(ebeam_default))
    else:
        raise KeyError(f"{dataset}: no beam-energy column/default")
    #endif

    t_raw = df[tcol].to_numpy(float)
    t_abs = np.abs(t_raw)

    out = pd.DataFrame({
        "dataset": dataset,
        "source_row": np.arange(len(df), dtype=int),
        "xB": df[xcol].to_numpy(float),
        "Q2": df[qcol].to_numpy(float),
        "t_abs": t_abs,
        "phi_deg": np.mod(df[pcol].to_numpy(float), 360.0),
        "ebeam": ebeam,
    })

    # Carry the exact KM15/BH decomposition when available. This avoids any
    # duplicated KM15 calculation in the model-selection stage.
    for col in [
        "km15_ep", "km15_bh", "km15_dvcs", "km15_int",
        "rbh", "delta_bh", "bh_A", "bh_B", "bh_C",
    ]:
        if col in df.columns:
            out[col] = df[col].to_numpy()
        #endif
    #endfor

    return out
#enddef


def assign_point_ids(points: pd.DataFrame) -> pd.DataFrame:
    out = points.copy()
    out.insert(
        0,
        "point_id",
        [f"{d}:{int(r)}" for d, r in zip(out["dataset"], out["source_row"])],
    )
    if out["point_id"].duplicated().any():
        raise RuntimeError("point_id collision in exported kinematics")
    #endif
    return out
#enddef


def load_points_from_cache(cache: Path) -> pd.DataFrame:
    df = pd.read_csv(cache)
    required = {"point_id", "dataset", "xB", "Q2", "t_abs", "phi_deg", "ebeam"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{cache}: missing columns {sorted(missing)}")
    #endif
    return df
#enddef


def make_partons_module_xml(model: str) -> str:
    """
    Return PARTONS module configuration.

    VGG99 uses its dedicated process implementation.
    GK16 uses the general GV08 process with standard LO CFF convolution.
    """
    if model == "vgg99":
        return """<module type="DVCSObservableModule" name="DVCSCrossSectionUUMinus">
<module type="DVCSProcessModule" name="DVCSProcessVGG99">
<module type="DVCSScalesModule" name="DVCSScalesQ2Multiplier"><param name="lambda" value="1."/></module>
<module type="DVCSXiConverterModule" name="DVCSXiConverterXBToXi"></module>
<module type="DVCSConvolCoeffFunctionModule" name="DVCSCFFStandard">
<param name="qcd_order_type" value="LO"/>
<module type="GPDModule" name="GPDVGG99"></module>
</module>
</module>
</module>"""
    elif model == "gk16":
        return """<module type="DVCSObservableModule" name="DVCSCrossSectionUUMinus">
<module type="DVCSProcessModule" name="DVCSProcessGV08">
<module type="DVCSScalesModule" name="DVCSScalesQ2Multiplier"><param name="lambda" value="1."/></module>
<module type="DVCSXiConverterModule" name="DVCSXiConverterXBToXi"></module>
<module type="DVCSConvolCoeffFunctionModule" name="DVCSCFFStandard">
<param name="qcd_order_type" value="LO"/>
<module type="GPDModule" name="GPDGK16"></module>
</module>
</module>
</module>"""
    else:
        raise ValueError(f"Unknown PARTONS model: {model}")
    #endif
#enddef


def make_single_point_task(row: pd.Series, model: str) -> str:
    """
    Generate one PARTONS task.

    PARTONS documentation defines DVCSObservableKinematic phi in degrees in
    the automation/DAO interface. Our exported phi_deg is the experimental
    Trento azimuth.
    """
    module_xml = make_partons_module_xml(model)
    return f"""<task service="DVCSObservableService" method="computeSingleKinematic" storeInDB="0">
<kinematics type="DVCSObservableKinematic">
<param name="xB" value="{row.xB:.12g}"/>
<param name="t" value="{-abs(row.t_abs):.12g}"/>
<param name="Q2" value="{row.Q2:.12g}"/>
<param name="E" value="{row.ebeam:.12g}"/>
<param name="phi" value="{row.phi_deg:.12g}"/>
</kinematics>
<computation_configuration>
{module_xml}
</computation_configuration>
</task>"""
#enddef


def write_partons_xml_chunks(
        points: pd.DataFrame,
        model: str,
        outdir: Path,
        chunk_size: int) -> pd.DataFrame:
    """
    Write deterministic XML chunks and a manifest.

    Chunking bounds failure recovery and permits process-level parallelism
    without using threads (important for GK16/CLN).
    """
    xml_dir = outdir / "partons_xml" / model
    xml_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    n = len(points)
    for ichunk, start in enumerate(range(0, n, chunk_size)):
        stop = min(start + chunk_size, n)
        part = points.iloc[start:stop]
        xml_path = xml_dir / f"{model}_{ichunk:05d}.xml"
        map_path = xml_dir / f"{model}_{ichunk:05d}_point_map.csv"

        tasks = "\n".join(
            make_single_point_task(row, model)
            for _, row in part.iterrows()
        )
        xml = (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>\n'
            f'<scenario date="2026-08-30" description="{escape(model)} BH-purity evaluation">\n'
            f"{tasks}\n"
            '<task service="DVCSObservableService" method="printResults"></task>\n'
            '</scenario>\n'
        )
        xml_path.write_text(xml)
        part[["point_id"]].to_csv(map_path, index=False)

        manifest_rows.append({
            "model": model,
            "chunk": ichunk,
            "start": start,
            "stop": stop,
            "n_points": len(part),
            "xml": str(xml_path),
            "point_map": str(map_path),
        })
    #endfor

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(xml_dir / f"{model}_manifest.csv", index=False)
    print(f"[PARTONS XML] {model}: {len(manifest)} chunk(s), {n} points -> {xml_dir}")
    return manifest
#enddef


def validate_model_results(
        points: pd.DataFrame,
        result: pd.DataFrame,
        model: str) -> pd.DataFrame:
    required = {"point_id", "ep_xs"}
    missing = required - set(result.columns)
    if missing:
        raise ValueError(f"{model} result missing {sorted(missing)}")
    #endif
    if result["point_id"].duplicated().any():
        raise ValueError(f"{model}: duplicate point_id values")
    #endif

    merged = points[["point_id"]].merge(
        result[["point_id", "ep_xs"]],
        on="point_id",
        how="left",
        validate="one_to_one",
    )
    coverage = np.mean(np.isfinite(merged["ep_xs"].to_numpy(float)))
    print(f"[model results] {model}: finite coverage={coverage:.3%}")
    if coverage < 0.999:
        raise RuntimeError(
            f"{model}: only {coverage:.3%} of points have finite EP cross sections"
        )
    #endif
    return merged.rename(columns={"ep_xs": f"{model}_ep"})
#enddef


def merge_models(
        points: pd.DataFrame,
        vgg_csv: Path,
        gk_csv: Path,
        outdir: Path,
        thresholds: List[float]) -> pd.DataFrame:
    """Merge KM15/VGG99/GK16 and construct all-model consensus selections."""
    out = points.copy()

    if "km15_bh" not in out.columns or "km15_ep" not in out.columns:
        raise RuntimeError(
            "Kinematic cache must contain km15_bh and km15_ep. Export it from "
            "the production analysis after KM15 evaluation."
        )
    #endif

    for model, path in [("vgg99", vgg_csv), ("gk16", gk_csv)]:
        result = validate_model_results(out, pd.read_csv(path), model)
        out = out.merge(result, on="point_id", how="left", validate="one_to_one")
    #endfor

    bh = out["km15_bh"].to_numpy(float)
    for model, ep_col in [
        ("km15", "km15_ep"),
        ("vgg99", "vgg99_ep"),
        ("gk16", "gk16_ep"),
    ]:
        ep = out[ep_col].to_numpy(float)
        out[f"delta_bh_{model}"] = np.abs(1.0 - bh / ep)
    #endfor

    deltas = out[
        ["delta_bh_km15", "delta_bh_vgg99", "delta_bh_gk16"]
    ].to_numpy(float)
    out["delta_bh_consensus_max"] = np.nanmax(deltas, axis=1)
    out["delta_bh_model_spread"] = np.nanmax(deltas, axis=1) - np.nanmin(
        deltas, axis=1
    )

    for threshold in thresholds:
        tag = f"{int(round(100 * threshold)):02d}pct"
        for model in ["km15", "vgg99", "gk16"]:
            out[f"pass_{tag}_{model}"] = out[f"delta_bh_{model}"] <= threshold
        #endfor
        out[f"pass_{tag}_all_models"] = (
            out[f"pass_{tag}_km15"]
            & out[f"pass_{tag}_vgg99"]
            & out[f"pass_{tag}_gk16"]
        )
    #endfor

    outdir.mkdir(parents=True, exist_ok=True)
    out.to_csv(outdir / "bh_model_selection_all_points.csv", index=False)

    summary_rows = []
    for dataset, d in out.groupby("dataset", sort=False):
        for threshold in thresholds:
            tag = f"{int(round(100 * threshold)):02d}pct"
            row = {
                "dataset": dataset,
                "threshold": threshold,
                "n_total": len(d),
            }
            for model in ["km15", "vgg99", "gk16", "all_models"]:
                row[f"n_{model}"] = int(d[f"pass_{tag}_{model}"].sum())
            #endfor
            summary_rows.append(row)
        #endfor
    #endfor
    pd.DataFrame(summary_rows).to_csv(
        outdir / "bh_model_selection_counts.csv", index=False
    )

    nominal = out.loc[out["pass_05pct_all_models"]].copy()
    nominal.to_csv(outdir / "bh_model_selection_consensus_05pct.csv", index=False)

    print(f"[selection] all-point cache -> {outdir / 'bh_model_selection_all_points.csv'}")
    print(f"[selection] nominal 5% all-model consensus: {len(nominal)}/{len(out)} points")
    return out
#enddef


def save_diagnostic_plots(df: pd.DataFrame, outdir: Path) -> None:
    import matplotlib.pyplot as plt

    models = ["km15", "vgg99", "gk16"]
    labels = {"km15": "KM15", "vgg99": "VGG99", "gk16": "GK16"}

    for dataset, d in df.groupby("dataset", sort=False):
        ddir = outdir / "plots" / dataset
        ddir.mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(figsize=(7.2, 6.2))
        for model in ["vgg99", "gk16"]:
            ax.scatter(
                100.0 * d["delta_bh_km15"],
                100.0 * d[f"delta_bh_{model}"],
                s=10,
                alpha=0.6,
                label=labels[model],
            )
        #endfor
        lim = max(
            10.0,
            float(np.nanpercentile(
                100.0 * d[
                    ["delta_bh_km15", "delta_bh_vgg99", "delta_bh_gk16"]
                ].to_numpy(float),
                99.0,
            )),
        )
        ax.plot([0, lim], [0, lim], linestyle="--", linewidth=1.0)
        ax.axvline(5.0, linestyle=":", linewidth=1.0)
        ax.axhline(5.0, linestyle=":", linewidth=1.0)
        ax.set_xlim(0.0, lim)
        ax.set_ylim(0.0, lim)
        ax.set_xlabel(r"KM15 $|1-\sigma_{\rm BH}/\sigma_{\rm EP}|$ (%)")
        ax.set_ylabel("alternative-model deviation (%)")
        ax.set_title(f"{dataset}: BH-purity model comparison")
        ax.legend()
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(ddir / "01_delta_bh_model_comparison.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.2, 5.4))
        thresholds = np.arange(1, 11)
        for model in models:
            counts = [
                int((d[f"delta_bh_{model}"] <= pct / 100.0).sum())
                for pct in thresholds
            ]
            ax.plot(thresholds, counts, marker="o", label=labels[model])
        #endfor
        consensus = [
            int((d["delta_bh_consensus_max"] <= pct / 100.0).sum())
            for pct in thresholds
        ]
        ax.plot(thresholds, consensus, marker="s", label="all-model consensus")
        ax.set_xlabel("BH-purity threshold (%)")
        ax.set_ylabel("selected points")
        ax.set_title(f"{dataset}: selected sample versus purity threshold")
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(ddir / "02_selection_count_vs_threshold.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.2, 5.4))
        disagreement = (
            d[["pass_05pct_km15", "pass_05pct_vgg99", "pass_05pct_gk16"]]
            .astype(int)
            .sum(axis=1)
        )
        labels_x = ["0 models", "1 model", "2 models", "3 models"]
        counts = [(disagreement == i).sum() for i in range(4)]
        ax.bar(np.arange(4), counts)
        ax.set_xticks(np.arange(4))
        ax.set_xticklabels(labels_x)
        ax.set_ylabel("points")
        ax.set_title(f"{dataset}: 5% BH-selection agreement")
        fig.tight_layout()
        fig.savefig(ddir / "03_selection_agreement_05pct.png", dpi=180)
        plt.close(fig)
    #endfor
#enddef


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="KM15/VGG99/GK16 BH-purity model-selection preparation"
    )
    p.add_argument("--check-partons", action="store_true")
    p.add_argument("--main-script", default=DEFAULT_MAIN)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--kinematics-cache", default=None)
    p.add_argument("--vgg-results", default=None)
    p.add_argument("--gk-results", default=None)
    p.add_argument("--make-partons-xml", action="store_true")
    p.add_argument("--chunk-size", type=int, default=250)
    p.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=DEFAULT_THRESHOLDS,
        help="BH-purity thresholds as fractions, default 0.01 ... 0.10",
    )
    return p
#enddef


def main(argv: List[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    outdir = Path(args.outdir)

    if args.check_partons:
        print_partons_check()
    #endif

    cache = (
        Path(args.kinematics_cache)
        if args.kinematics_cache
        else outdir / "bh_model_kinematics.csv"
    )

    if not cache.exists():
        print(
            "\nNo kinematic cache exists yet.\n"
            "The production EMFF script should export its already-evaluated "
            "KM15 points to:\n"
            f"  {cache}\n"
            "This evaluator intentionally does not duplicate the data-loading "
            "and KM15 calculation logic.\n"
        )
        return 2
    #endif

    points = load_points_from_cache(cache)
    print(f"[kinematics] loaded {len(points)} exact experimental points from {cache}")

    if args.make_partons_xml:
        write_partons_xml_chunks(points, "vgg99", outdir, args.chunk_size)
        write_partons_xml_chunks(points, "gk16", outdir, args.chunk_size)
    #endif

    if args.vgg_results and args.gk_results:
        merged = merge_models(
            points,
            Path(args.vgg_results),
            Path(args.gk_results),
            outdir,
            args.thresholds,
        )
        save_diagnostic_plots(merged, outdir)
    #endif

    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif
