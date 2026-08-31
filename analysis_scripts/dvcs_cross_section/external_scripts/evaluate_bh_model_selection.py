#!/usr/bin/env python3
"""
Evaluate KM15 / VGG99 / GK16 BH-purity model selection.

This script consumes the exact kinematic cache exported by
extract_emff_from_dvcs_bh.py and performs the expensive alternative-model
calculations once.  All later BH-threshold scans and EMFF fits can use the
cached CSV products without rerunning PARTONS.

Model definitions
-----------------
KM15:
    Gepard KM15 total cross section and Gepard BH subprocess.  These values are
    already stored in bh_model_kinematics.csv.

GK16:
    PARTONS GPDGK16 + DVCSCFFStandard(LO) + DVCSProcessGV08.

VGG99:
    PARTONS GPDVGG99 + DVCSCFFStandard(LO) + DVCSProcessVGG99.
    GPDVGG99 is configured explicitly with LHAPDF set MSTW2008nlo68cl,
    member 0, matching the PDF set installed in the PARTONS 5.0 farm image.

The purity ratios are kept internally self-consistent:

    delta_KM15 = |1 - BH_Gepard / EP_Gepard,KM15|

    delta_GK16 = |1 - BH_PARTONS,GV08 / EP_PARTONS,GK16|

    delta_VGG99 = |1 - BH_PARTONS,VGG99 / EP_PARTONS,VGG99|

Two PARTONS BH calculations are deliberately retained.  Even though the BH
subprocess should not depend on the GPD, GK16 and VGG99 use different PARTONS
DVCS process modules.  Evaluating BH with the same process module used for each
full cross section removes a needless process-implementation ambiguity.

Azimuth convention
------------------
The imported Gepard Jo/Hall-A datasets use the previously validated mapping

    phi_PARTONS = (180 deg - phi_Gepard) mod 360 deg.

CLAS12/pass1 and CLAS12/pass2 are NOT assumed to share that convention.  By default this script
runs a small BH-only PARTONS/Gepard closure scan over four candidate mappings
(identity, negative, 180-minus, 180-plus), using representative CLAS12 points.
The winning mapping is chosen without measured cross sections and is then used
for the full CLAS12 PARTONS calculation.

Production PARTONS caches are namespaced by the active dataset-specific phi
mappings so a cache generated with an obsolete convention cannot be silently
reused after the convention changes.

PARTONS execution
-----------------
The default farm setup is:

    SIF:      partons_v4.sif
    project:  /scratch/thayward/partons-example
    binary:   ./bin/PARTONS_example

The SIF filename is historical; the executable reports PARTONS 5.0.0.

GK16/CLN is not thread-safe.  This driver parallelizes only by launching
separate Apptainer/PARTONS operating-system processes.  Each PARTONS instance
evaluates its chunk sequentially.

Important output files
----------------------
    bh_model_selection_all_points.csv
    bh_model_selection_counts.csv
    bh_model_selection_pairwise_overlap.csv
    bh_model_selection_km15_05pct.csv
    bh_model_selection_vgg99_05pct.csv
    bh_model_selection_gk16_05pct.csv
    bh_model_selection_consensus_05pct.csv

By default all generated PARTONS/model-prediction artifacts are written under:
    /work/clas12/thayward/CLAS12_exclusive/dvcs/model_predictions/

This includes:
    partons_results/
    partons_logs/
    partons_xml/
    diagnostics and model-selection CSV/PNG products.

The small input kinematic cache remains in the analysis output tree by default:
    output/emff_from_bh_paper_method/bh_model_selection/bh_model_kinematics.csv

Chunk XML, point maps, stdout/stderr, return codes, parsed-result counts, and
warning/error counts are retained for reproducibility and failure recovery.

Performance defaults
--------------------
The default production settings are 4 independent OS processes and 200
sequential points per Apptainer invocation.  This reduces container/module
startup count from 154 invocations per 3846-point calculation at chunk-size 25
to only 20, while retaining enough chunks for four-way load balancing.
Common numerical-library thread counts are constrained to one per PARTONS
process to avoid nested oversubscription.

PARTONS 5.0 behavior verified on the JLab farm:
  * generated XML directories must be explicitly bound into the container;
  * in a multi-point scenario, each computeSingleKinematic task must be followed
    immediately by printResults, otherwise only the final point is printed.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DEFAULT_OUTDIR = "/work/clas12/thayward/CLAS12_exclusive/dvcs/model_predictions"
DEFAULT_KINEMATICS_CACHE = "output/emff_from_bh_paper_method/bh_model_selection/bh_model_kinematics.csv"
DEFAULT_THRESHOLDS = [0.01 * i for i in range(1, 11)]
DEFAULT_PARTONS_SIF = "partons_v4.sif"
DEFAULT_PARTONS_PROJECT = "/scratch/thayward/partons-example"
DEFAULT_PARTONS_EXECUTABLE = "./bin/PARTONS_example"
DEFAULT_PARTONS_WORKERS = 4
DEFAULT_PARTONS_CHUNK_SIZE = 200

# GPDVGG99 requires an LHAPDF forward-PDF set.  The PARTONS 5.0 farm image
# used here contains MSTW2008nlo68cl, and GPDVGG99 exposes the exact parameters
# "setName" and "member".
DEFAULT_VGG99_PDF_SET = "MSTW2008nlo68cl"
DEFAULT_VGG99_PDF_MEMBER = 0

PARTONS_RESULT_RE = re.compile(
    r"Result:\s*([+-]?(?:[\d,]+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*\[([^\]]+)\]"
)
PARTONS_WARNING_TOKEN = "[WARN]"
PARTONS_ERROR_TOKEN = "[ERROR]"
PARTONS_INTEGRATOR_WARNING = "Cannot reach tolerances"


PHI_MAPPING_MODES = (
    "identity",
    "negative",
    "180-minus",
    "180-plus",
)
DEFAULT_IMPORTED_PHI_MAPPING = "180-minus"
DEFAULT_CLAS12_PHI_MAPPING = "auto"
DEFAULT_CLAS12_SCAN_POINTS = 12
DEFAULT_IMPORTED_SCAN_POINTS = 16

# phi_deg in the kinematic cache is not in one universal source convention.
# Gepard-bundled datasets are stored in Gepard's BMK representation and need
# the validated 180-degree shift before they are sent to PARTONS.  Direct
# published tables/spreadsheets and CLAS12 tables retain their published
# Trento-like azimuth and therefore use identity.
#
# Saylor is deliberately listed separately rather than being swept into a
# generic "imported" category.  Its direct supplemental-table convention is
# audited by --validate-phi-conventions before its model comparison is used
# scientifically.
FIXED_DATASET_PHI_MAPPINGS = {
    "jo2015": "180-minus",
    "halla_defurne2015": "180-minus",
    "saylor2018": "identity",
    "halla_georges2022": "identity",
}

FIXED_DATASET_PHI_SOURCES = {
    "jo2015": "Gepard dataset 98 (BMK-stored)",
    "halla_defurne2015": "Gepard datasets 107/108/112/113 (BMK-stored)",
    "saylor2018": "published supplemental table (direct)",
    "halla_georges2022": "published E12-06-114 spreadsheet (Trento)",
    "pass1": "CLAS12 Lee table",
    "pass2": "CLAS12 Hayward table",
}


def transform_phi_to_partons(
        phi_deg: np.ndarray | Sequence[float] | float,
        mode: str):
    """Apply one explicit candidate azimuth transformation."""
    arr = np.asarray(phi_deg, dtype=float)
    if mode == "identity":
        transformed = np.mod(arr, 360.0)
    elif mode == "negative":
        transformed = np.mod(-arr, 360.0)
    elif mode == "180-minus":
        transformed = np.mod(180.0 - arr, 360.0)
    elif mode == "180-plus":
        transformed = np.mod(180.0 + arr, 360.0)
    else:
        raise ValueError(
            f"Unknown phi mapping {mode!r}; choose from {PHI_MAPPING_MODES}"
        )
    #endif

    if np.ndim(phi_deg) == 0:
        return float(transformed)
    #endif
    return transformed
#enddef


def apply_dataset_phi_mappings(
        df: pd.DataFrame,
        imported_mode: str = DEFAULT_IMPORTED_PHI_MAPPING,
        clas12_mode: str = "identity") -> pd.DataFrame:
    """
    Attach phi_partons_deg using the convention of each dataset's source.

    `imported_mode` is retained only as a backwards-compatible override for
    the Gepard-bundled Jo and Defurne datasets.  It is intentionally NOT
    applied to Georges or Saylor, which come from direct publication tables.

    CLAS12 pass1/pass2 share the independently BH-validated CLAS12 mapping.
    """
    out = df.copy()
    out["phi_partons_deg"] = np.nan
    out["phi_mapping_mode"] = ""
    out["phi_source_representation"] = ""

    dataset = out["dataset"].astype(str)

    # Gepard-bundled datasets: the cached phi is the Gepard/BMK representation.
    for key in ["jo2015", "halla_defurne2015"]:
        mask = dataset == key
        if np.any(mask):
            out.loc[mask, "phi_partons_deg"] = transform_phi_to_partons(
                out.loc[mask, "phi_deg"].to_numpy(float),
                imported_mode,
            )
            out.loc[mask, "phi_mapping_mode"] = imported_mode
            out.loc[mask, "phi_source_representation"] = (
                FIXED_DATASET_PHI_SOURCES[key]
            )
        #endif
    #endfor

    # Direct publication inputs must not inherit the Gepard/BMK mapping.
    for key in ["saylor2018", "halla_georges2022"]:
        mask = dataset == key
        if np.any(mask):
            mode = FIXED_DATASET_PHI_MAPPINGS[key]
            out.loc[mask, "phi_partons_deg"] = transform_phi_to_partons(
                out.loc[mask, "phi_deg"].to_numpy(float),
                mode,
            )
            out.loc[mask, "phi_mapping_mode"] = mode
            out.loc[mask, "phi_source_representation"] = (
                FIXED_DATASET_PHI_SOURCES[key]
            )
        #endif
    #endfor

    # Both CLAS12 analyses use the same measured phi convention.
    clas12_mask = dataset.isin(["pass1", "pass2"])
    if np.any(clas12_mask):
        out.loc[clas12_mask, "phi_partons_deg"] = transform_phi_to_partons(
            out.loc[clas12_mask, "phi_deg"].to_numpy(float),
            clas12_mode,
        )
        out.loc[clas12_mask, "phi_mapping_mode"] = clas12_mode
        out.loc[clas12_mask, "phi_source_representation"] = dataset.loc[
            clas12_mask
        ].map(FIXED_DATASET_PHI_SOURCES).to_numpy()
    #endif

    # Fail loudly if a future dataset is added without an explicit convention.
    unresolved = ~np.isfinite(
        pd.to_numeric(out["phi_partons_deg"], errors="coerce").to_numpy(float)
    )
    if np.any(unresolved):
        unknown = sorted(set(dataset.loc[unresolved].astype(str)))
        raise RuntimeError(
            "No explicit PARTONS phi convention is configured for dataset(s): "
            + ", ".join(unknown)
        )
    #endif

    return out
#enddef

def load_points_from_cache(cache: Path) -> pd.DataFrame:
    """Read and validate the exact experimental-point cache."""
    df = pd.read_csv(cache)
    required = {
        "point_id", "dataset", "source_row", "xB", "Q2",
        "t_abs", "phi_deg", "ebeam", "km15_ep", "km15_bh",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{cache}: missing columns {sorted(missing)}")
    #endif

    if df["point_id"].duplicated().any():
        duplicated = df.loc[df["point_id"].duplicated(), "point_id"].tolist()
        raise ValueError(f"{cache}: duplicate point_id values, e.g. {duplicated[:5]}")
    #endif

    for col in ["xB", "Q2", "t_abs", "phi_deg", "ebeam", "km15_ep", "km15_bh"]:
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
        if not np.all(np.isfinite(vals)):
            raise ValueError(f"{cache}: nonfinite values in required column {col}")
        #endif
    #endfor

    df = df.copy()
    # phi_partons_deg is attached only after the CLAS12-specific convention has
    # either been explicitly supplied or validated by the BH-only scan.
    df["delta_bh_km15"] = np.abs(
        1.0 - df["km15_bh"].to_numpy(float) / df["km15_ep"].to_numpy(float)
    )
    return df
#enddef


def locate_partons(sif: Path, project: Path, executable: str) -> None:
    """Print the farm execution configuration and obvious preflight failures."""
    print("[PARTONS preflight]")
    print(f"  apptainer : {shutil.which('apptainer') or 'NOT FOUND'}")
    print(f"  SIF       : {sif}")
    print(f"  project   : {project}")
    print(f"  executable: {executable}")
    print(f"  /scratch  : {'present' if Path('/scratch').exists() else 'NOT FOUND'}")

    failures = []
    if shutil.which("apptainer") is None:
        failures.append("apptainer is not on PATH")
    #endif
    if not sif.exists():
        failures.append(f"SIF does not exist: {sif}")
    #endif
    if not project.exists():
        failures.append(f"PARTONS project copy does not exist: {project}")
    #endif
    if not (project / "bin" / "tmp").exists():
        failures.append(f"writable PARTONS logger directory missing: {project / 'bin' / 'tmp'}")
    #endif

    if failures:
        raise RuntimeError("PARTONS preflight failed:\n  - " + "\n  - ".join(failures))
    #endif
#enddef


def _process_module_xml(process: str, gpd: str) -> str:
    """Return the common PARTONS process-module block."""
    if gpd == "GPDVGG99":
        gpd_xml = f"""<module type="GPDModule" name="{gpd}">
<param name="setName" value="{DEFAULT_VGG99_PDF_SET}" />
<param name="member" value="{DEFAULT_VGG99_PDF_MEMBER}" />
</module>"""
    else:
        gpd_xml = f'<module type="GPDModule" name="{gpd}"></module>'
    #endif

    return f"""<module type="DVCSProcessModule" name="{process}">
<module type="DVCSScalesModule" name="DVCSScalesQ2Multiplier">
<param name="lambda" value="1." />
</module>
<module type="DVCSXiConverterModule" name="DVCSXiConverterXBToXi"></module>
<module type="DVCSConvolCoeffFunctionModule" name="DVCSCFFStandard">
<param name="qcd_order_type" value="LO" />
{gpd_xml}
</module>
</module>"""
#enddef


def partons_observable_xml(calculation: str) -> str:
    """
    Return the PARTONS observable module for one cached quantity.

    Calculations:
      bh_gv08  : BH subprocess with DVCSProcessGV08 / GPDGK16
      gk16     : full negative-charge UU cross section with GV08 / GK16
      bh_vgg99 : BH subprocess with DVCSProcessVGG99 / GPDVGG99
      vgg99    : full negative-charge UU cross section with VGG99 / VGG99
    """
    if calculation == "bh_gv08":
        observable = "DVCSCrossSectionUUBHSubProc"
        process = "DVCSProcessGV08"
        gpd = "GPDGK16"
    elif calculation == "gk16":
        observable = "DVCSCrossSectionUUMinus"
        process = "DVCSProcessGV08"
        gpd = "GPDGK16"
    elif calculation == "bh_vgg99":
        observable = "DVCSCrossSectionUUBHSubProc"
        process = "DVCSProcessVGG99"
        gpd = "GPDVGG99"
    elif calculation == "vgg99":
        observable = "DVCSCrossSectionUUMinus"
        process = "DVCSProcessVGG99"
        gpd = "GPDVGG99"
    else:
        raise ValueError(f"Unknown PARTONS calculation: {calculation}")
    #endif

    return (
        f'<module type="DVCSObservableModule" name="{observable}">\n'
        f'{_process_module_xml(process, gpd)}\n'
        '</module>'
    )
#enddef


def make_single_point_task(row: pd.Series, calculation: str) -> str:
    """Generate one sequential PARTONS compute task for an exact experimental point."""
    module_xml = partons_observable_xml(calculation)
    return f"""<task service="DVCSObservableService" method="computeSingleKinematic" storeInDB="0">
<kinematics type="DVCSObservableKinematic">
<param name="xB" value="{float(row.xB):.15g}" />
<param name="t" value="{-abs(float(row.t_abs)):.15g}" />
<param name="Q2" value="{float(row.Q2):.15g}" />
<param name="E" value="{float(row.ebeam):.15g}" />
<param name="phi" value="{float(row.phi_partons_deg):.15g}" />
</kinematics>
<computation_configuration>
{module_xml}
</computation_configuration>
</task>
<task service="DVCSObservableService" method="printResults"></task>"""
#enddef


def write_partons_xml_chunks(
        points: pd.DataFrame,
        calculation: str,
        outdir: Path,
        chunk_size: int) -> pd.DataFrame:
    """Write deterministic XML chunks and their one-to-one point maps."""
    if chunk_size < 1:
        raise ValueError("--chunk-size must be >= 1")
    #endif

    xml_dir = outdir / "partons_xml" / calculation
    xml_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    n = len(points)
    for ichunk, start in enumerate(range(0, n, chunk_size)):
        stop = min(start + chunk_size, n)
        part = points.iloc[start:stop].copy()

        xml_path = xml_dir / f"{calculation}_{ichunk:05d}.xml"
        map_path = xml_dir / f"{calculation}_{ichunk:05d}_point_map.csv"

        tasks = "\n".join(
            make_single_point_task(row, calculation)
            for _, row in part.iterrows()
        )
        # PARTONS 5.0 keeps only the current observable result in its print
        # buffer.  Therefore each compute task must be followed immediately by
        # printResults; one printResults at the end of a multi-point chunk would
        # print only the final point.
        xml_text = (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>\n'
            f'<scenario date="2026-08-30" description="{calculation} BH-purity evaluation">\n'
            f'{tasks}\n'
            '</scenario>\n'
        )

        xml_path.write_text(xml_text)
        part[
            [
                "point_id", "dataset", "source_row", "xB", "Q2",
                "t_abs", "phi_deg", "phi_partons_deg", "phi_mapping_mode",
                "ebeam",
            ]
        ].to_csv(map_path, index=False)

        manifest_rows.append(
            {
                "calculation": calculation,
                "chunk": ichunk,
                "start": start,
                "stop": stop,
                "n_points": len(part),
                "xml": str(xml_path.resolve()),
                "point_map": str(map_path.resolve()),
            }
        )
    #endfor

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = xml_dir / f"{calculation}_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    print(
        f"[PARTONS XML] {calculation}: {len(manifest)} chunk(s), "
        f"{n} points -> {manifest_path}"
    )
    return manifest
#enddef


def parse_partons_stdout(stdout: str) -> tuple[List[float], List[str], int, int, int]:
    """Parse observable values, units, warnings, integrator warnings, and errors."""
    matches = PARTONS_RESULT_RE.findall(stdout)
    values = [float(value.replace(",", "")) for value, _ in matches]
    units = [unit.strip() for _, unit in matches]
    n_warnings = stdout.count(PARTONS_WARNING_TOKEN)
    n_integrator_warnings = stdout.count(PARTONS_INTEGRATOR_WARNING)
    n_errors = stdout.count(PARTONS_ERROR_TOKEN)
    return values, units, n_warnings, n_integrator_warnings, n_errors
#enddef


def _parse_retained_partons_chunk(job: dict) -> dict | None:
    """
    Reparse and reuse a complete retained PARTONS chunk log.

    Existing stdout/stderr files are treated as the raw cache.  A retained
    chunk is reused only when it contains exactly the expected number of finite
    nb-valued results, contains no PARTONS [ERROR], and closed properly.
    """
    point_map_path = Path(job["point_map"])
    stdout_path = Path(job["stdout_path"])
    stderr_path = Path(job["stderr_path"])

    if not stdout_path.exists():
        return None
    #endif

    stdout = stdout_path.read_text(errors="replace")
    stderr = stderr_path.read_text(errors="replace") if stderr_path.exists() else ""
    combined = stdout + "\n" + stderr

    values, units, n_warnings, n_integrator_warnings, n_errors = parse_partons_stdout(
        combined
    )

    point_map = pd.read_csv(point_map_path)
    n_expected = len(point_map)
    parse_ok = len(values) == n_expected
    finite_ok = parse_ok and np.all(np.isfinite(np.asarray(values, dtype=float)))
    units_ok = parse_ok and all(unit.lower() == "nb" for unit in units)
    closed_ok = "Closed properly" in combined

    if not (n_errors == 0 and parse_ok and finite_ok and units_ok and closed_ok):
        return None
    #endif

    status = {
        "calculation": job["calculation"],
        "chunk": int(job["chunk"]),
        "n_expected": int(n_expected),
        "n_parsed": int(len(values)),
        "returncode": 0,
        "parse_ok": True,
        "finite_ok": True,
        "units_ok": True,
        "n_warnings": int(n_warnings),
        "n_integrator_warnings": int(n_integrator_warnings),
        "n_errors": int(n_errors),
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
        "xml": str(Path(job["xml"])),
        "point_map": str(point_map_path),
        "reused_log": True,
        "elapsed_s": 0.0,
    }

    result_rows = []
    for point_id, value, unit in zip(point_map["point_id"], values, units):
        result_rows.append(
            {
                "point_id": point_id,
                "xs_nb": value,
                "unit": unit,
                "chunk": int(job["chunk"]),
            }
        )
    #endfor

    return {"status": status, "results": result_rows}
#enddef


def _run_partons_chunk(job: dict) -> dict:
    """
    Execute one XML chunk in an independent Apptainer/PARTONS OS process.

    Unless force=True, first attempt to recover a complete prior calculation
    directly from its retained stdout/stderr logs.
    """
    if not job.get("force", False):
        recovered = _parse_retained_partons_chunk(job)
        if recovered is not None:
            return recovered
        #endif
    #endif

    xml_path = Path(job["xml"])
    point_map_path = Path(job["point_map"])
    stdout_path = Path(job["stdout_path"])
    stderr_path = Path(job["stderr_path"])

    # The generated XML lives under the analysis output tree, which is not
    # necessarily visible inside the PARTONS container.  Bind the XML directory
    # at the identical absolute path in addition to /scratch.
    xml_bind = f"{xml_path.parent.resolve()}:{xml_path.parent.resolve()}"
    cmd = [
        "apptainer",
        "exec",
        "--bind",
        "/scratch:/scratch",
        "--bind",
        xml_bind,
        "--pwd",
        job["project"],
        job["sif"],
        job["executable"],
        str(xml_path.resolve()),
    ]

    run_env = os.environ.copy()
    # We parallelize explicitly with independent PARTONS OS processes.  Keep
    # common numerical backends single-threaded inside each process so four
    # workers do not accidentally become 4 x N nested threads.
    run_env.setdefault("OMP_NUM_THREADS", "1")
    run_env.setdefault("OPENBLAS_NUM_THREADS", "1")
    run_env.setdefault("MKL_NUM_THREADS", "1")
    run_env.setdefault("NUMEXPR_NUM_THREADS", "1")

    t0 = time.perf_counter()
    proc = subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=run_env,
    )
    elapsed_s = time.perf_counter() - t0

    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stdout_path.write_text(proc.stdout)
    stderr_path.write_text(proc.stderr)

    values, units, n_warnings, n_integrator_warnings, n_errors = parse_partons_stdout(
        proc.stdout + "\n" + proc.stderr
    )

    point_map = pd.read_csv(point_map_path)
    n_expected = len(point_map)
    parse_ok = len(values) == n_expected
    finite_ok = parse_ok and np.all(np.isfinite(np.asarray(values, dtype=float)))
    units_ok = parse_ok and all(unit.lower() == "nb" for unit in units)

    status = {
        "calculation": job["calculation"],
        "chunk": int(job["chunk"]),
        "n_expected": int(n_expected),
        "n_parsed": int(len(values)),
        "returncode": int(proc.returncode),
        "parse_ok": bool(parse_ok),
        "finite_ok": bool(finite_ok),
        "units_ok": bool(units_ok),
        "n_warnings": int(n_warnings),
        "n_integrator_warnings": int(n_integrator_warnings),
        "n_errors": int(n_errors),
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
        "xml": str(xml_path),
        "point_map": str(point_map_path),
        "reused_log": False,
        "elapsed_s": float(elapsed_s),
    }

    result_rows = []
    if parse_ok:
        for point_id, value, unit in zip(point_map["point_id"], values, units):
            result_rows.append(
                {
                    "point_id": point_id,
                    "xs_nb": value,
                    "unit": unit,
                    "chunk": int(job["chunk"]),
                }
            )
        #endfor
    #endif

    return {"status": status, "results": result_rows}
#enddef


def run_partons_calculation(
        points: pd.DataFrame,
        calculation: str,
        outdir: Path,
        sif: Path,
        project: Path,
        executable: str,
        workers: int,
        chunk_size: int,
        force: bool,
        cache_tag: str = "default") -> pd.DataFrame:
    """
    Run, resume from retained chunk logs, or reuse one complete PARTONS cache.

    cache_tag is part of the cache namespace.  This prevents calculations made
    with one phi convention from being silently reused after the convention is
    changed.
    """
    safe_tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(cache_tag))
    result_dir = outdir / "partons_results" / safe_tag
    result_dir.mkdir(parents=True, exist_ok=True)

    result_path = result_dir / f"{calculation}_results.csv"
    status_path = result_dir / f"{calculation}_chunk_status.csv"

    if result_path.exists() and not force:
        cached = pd.read_csv(result_path)
        required = {"point_id", "xs_nb"}
        if required <= set(cached.columns):
            merged = points[["point_id"]].merge(
                cached[["point_id", "xs_nb"]],
                on="point_id",
                how="left",
                validate="one_to_one",
            )
            coverage = np.mean(np.isfinite(merged["xs_nb"].to_numpy(float)))
            if coverage == 1.0 and len(cached) == len(points):
                print(
                    f"[PARTONS cache] {calculation}: reusing {len(cached)} points "
                    f"from {result_path}"
                )
                return cached
            #endif
        #endif
        print(f"[PARTONS cache] {calculation}: cache incomplete; recalculating")
    #endif

    namespaced_outdir = outdir / "partons_cache_namespaces" / safe_tag
    manifest = write_partons_xml_chunks(
        points=points,
        calculation=calculation,
        outdir=namespaced_outdir,
        chunk_size=chunk_size,
    )

    log_dir = outdir / "partons_logs" / safe_tag / calculation
    log_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    sif_abs = str(sif.resolve())
    project_abs = str(project.resolve())
    for row in manifest.itertuples(index=False):
        jobs.append(
            {
                "calculation": calculation,
                "chunk": int(row.chunk),
                "xml": row.xml,
                "point_map": row.point_map,
                "stdout_path": str(log_dir / f"{calculation}_{int(row.chunk):05d}.stdout.txt"),
                "stderr_path": str(log_dir / f"{calculation}_{int(row.chunk):05d}.stderr.txt"),
                "sif": sif_abs,
                "project": project_abs,
                "executable": executable,
                "force": bool(force),
            }
        )
    #endfor

    statuses = []
    result_rows = []

    print(
        f"[PARTONS run] {calculation}: {len(points)} points, "
        f"{len(jobs)} chunks, {workers} OS process(es)"
    )

    if workers <= 1:
        for ijob, job in enumerate(jobs, start=1):
            payload = _run_partons_chunk(job)
            statuses.append(payload["status"])
            result_rows.extend(payload["results"])
            print(
                f"[PARTONS run] {calculation}: chunk {ijob}/{len(jobs)} "
                f"rc={payload['status']['returncode']} "
                f"parsed={payload['status']['n_parsed']}/{payload['status']['n_expected']} "
                f"warnings={payload['status']['n_warnings']} "
                f"errors={payload['status']['n_errors']} "
                f"time={payload['status'].get('elapsed_s', 0.0):.1f}s "
                f"{'[reused log]' if payload['status'].get('reused_log', False) else '[executed]'}"
            )

            s = payload["status"]
            bad = (
                s["returncode"] != 0
                or s["n_errors"] != 0
                or not s["parse_ok"]
                or not s["finite_ok"]
                or not s["units_ok"]
            )
            if bad:
                pd.DataFrame(statuses).to_csv(status_path, index=False)
                raise RuntimeError(
                    f"{calculation}: chunk {s['chunk']} failed immediately "
                    f"(rc={s['returncode']}, errors={s['n_errors']}, "
                    f"parsed={s['n_parsed']}/{s['n_expected']}). "
                    f"See {s['stdout']} and {s['stderr']}."
                )
            #endif
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            future_map = {
                pool.submit(_run_partons_chunk, job): job for job in jobs
            }
            done = 0
            for future in as_completed(future_map):
                payload = future.result()
                statuses.append(payload["status"])
                result_rows.extend(payload["results"])
                done += 1
                print(
                    f"[PARTONS run] {calculation}: chunk {done}/{len(jobs)} "
                    f"(id={payload['status']['chunk']}) "
                    f"rc={payload['status']['returncode']} "
                    f"parsed={payload['status']['n_parsed']}/{payload['status']['n_expected']} "
                    f"warnings={payload['status']['n_warnings']} "
                    f"errors={payload['status']['n_errors']} "
                    f"time={payload['status'].get('elapsed_s', 0.0):.1f}s "
                    f"{'[reused log]' if payload['status'].get('reused_log', False) else '[executed]'}"
                )
            #endfor
        #endwith
    #endif

    status_df = pd.DataFrame(statuses).sort_values("chunk").reset_index(drop=True)
    status_df.to_csv(status_path, index=False)

    failed = status_df.loc[
        (status_df["returncode"] != 0)
        | (status_df["n_errors"] != 0)
        | (~status_df["parse_ok"])
        | (~status_df["finite_ok"])
        | (~status_df["units_ok"])
    ]
    if len(failed):
        print(f"[PARTONS run] {calculation}: FAILED CHUNKS")
        print(
            failed[
                [
                    "chunk", "returncode", "n_errors", "n_expected", "n_parsed",
                    "parse_ok", "finite_ok", "units_ok", "stdout", "stderr",
                ]
            ].to_string(index=False)
        )
        raise RuntimeError(
            f"{calculation}: {len(failed)} PARTONS chunk(s) failed. "
            f"See {status_path} and the retained logs."
        )
    #endif

    result_df = pd.DataFrame(result_rows)
    if result_df["point_id"].duplicated().any():
        raise RuntimeError(f"{calculation}: duplicate point_id in parsed PARTONS output")
    #endif

    result_df = points[["point_id"]].merge(
        result_df,
        on="point_id",
        how="left",
        validate="one_to_one",
    )
    if len(result_df) != len(points) or not np.all(
        np.isfinite(result_df["xs_nb"].to_numpy(float))
    ):
        raise RuntimeError(f"{calculation}: incomplete/nonfinite final result cache")
    #endif

    result_df.to_csv(result_path, index=False)

    executed = status_df.loc[~status_df["reused_log"].astype(bool)]
    worker_seconds = float(executed["elapsed_s"].sum()) if len(executed) else 0.0
    print(
        f"[PARTONS run] {calculation}: completed {len(result_df)} points; "
        f"warnings={int(status_df['n_warnings'].sum())}, "
        f"integrator warnings={int(status_df['n_integrator_warnings'].sum())}, "
        f"executed_chunks={len(executed)}, worker_seconds={worker_seconds:.1f}"
    )
    print(f"[PARTONS cache] {calculation} -> {result_path}")
    return result_df
#enddef


def _attach_result(
        points: pd.DataFrame,
        result: pd.DataFrame,
        output_column: str) -> pd.DataFrame:
    """Attach one complete PARTONS cache by point_id."""
    result = result[["point_id", "xs_nb"]].rename(columns={"xs_nb": output_column})
    return points.merge(result, on="point_id", how="left", validate="one_to_one")
#enddef



def choose_clas12_phi_scan_points(
        points: pd.DataFrame,
        n_points: int = DEFAULT_CLAS12_SCAN_POINTS) -> pd.DataFrame:
    """
    Choose a small deterministic CLAS12 sample spanning phi and kinematics.

    Points are selected approximately uniformly in sorted phi, with duplicate
    point_ids removed.  This is a convention test, not a physics selection.
    """
    d = points.loc[points["dataset"].astype(str) == "pass1"].copy()
    if len(d) == 0:
        raise RuntimeError(
            "No dataset=='pass1' points are available for the CLAS12 phi scan"
        )
    #endif
    n_points = max(4, min(int(n_points), len(d)))
    d = d.sort_values(["phi_deg", "xB", "Q2", "t_abs"]).reset_index(drop=True)
    indices = np.unique(
        np.rint(np.linspace(0, len(d) - 1, n_points)).astype(int)
    )
    sample = d.iloc[indices].copy().reset_index(drop=True)
    return sample
#enddef


def run_clas12_phi_mapping_scan(
        points: pd.DataFrame,
        outdir: Path,
        sif: Path,
        project: Path,
        executable: str,
        workers: int,
        scan_points: int,
        force: bool = False) -> tuple[str, pd.DataFrame]:
    """
    Determine the CLAS12 phi mapping using BH-only closure against Gepard.

    Four plausible transformations are tested:
      identity   : phi_P = phi_C
      negative   : phi_P = -phi_C
      180-minus  : phi_P = 180 - phi_C
      180-plus   : phi_P = 180 + phi_C

    Both PARTONS BH implementations (GV08 and VGG99) are compared against the
    cached Gepard BH values.  No measured cross sections and no full EP/GPD
    predictions enter the decision.

    The selected mapping minimizes a combined score based on:
      * log-ratio scatter across points (dominant term);
      * median absolute log offset from unity (secondary term).

    For unpolarized BH, phi and -phi are physically indistinguishable because
    the observable is even in phi.  The scan therefore ranks the two BH
    equivalence classes {identity, negative} and {180-minus, 180-plus}, then
    uses identity or 180-minus as the canonical representative.

    A mapping is accepted only when both PARTONS BH processes show a stable
    point-by-point ratio: all finite/positive, >=90% within 0.5--2.0, and
    log-ratio standard deviation <=0.15 for each process.
    """
    scan_root = outdir / "phi_mapping_validation" / "clas12"
    scan_root.mkdir(parents=True, exist_ok=True)

    base = choose_clas12_phi_scan_points(points, n_points=scan_points)
    base[
        [
            "point_id", "dataset", "source_row", "xB", "Q2",
            "t_abs", "phi_deg", "ebeam", "km15_bh",
        ]
    ].to_csv(scan_root / "scan_points.csv", index=False)

    summary_rows = []
    point_rows = []

    for mode in PHI_MAPPING_MODES:
        test = base.copy()
        test["phi_partons_deg"] = transform_phi_to_partons(
            test["phi_deg"].to_numpy(float), mode
        )
        test["phi_mapping_mode"] = mode

        mode_out = scan_root / mode
        mode_out.mkdir(parents=True, exist_ok=True)
        results = {}

        for calculation in ["bh_gv08", "bh_vgg99"]:
            results[calculation] = run_partons_calculation(
                points=test,
                calculation=calculation,
                outdir=mode_out,
                sif=sif,
                project=project,
                executable=executable,
                workers=min(max(1, workers), 2),
                chunk_size=max(1, len(test)),
                force=force,
                cache_tag=f"clas12_phi_{mode}",
            )
        #endfor

        merged = test[
            [
                "point_id", "source_row", "xB", "Q2", "t_abs",
                "phi_deg", "phi_partons_deg", "km15_bh",
            ]
        ].copy()
        merged = merged.merge(
            results["bh_gv08"][["point_id", "xs_nb"]].rename(
                columns={"xs_nb": "partons_bh_gv08"}
            ),
            on="point_id", how="left", validate="one_to_one",
        )
        merged = merged.merge(
            results["bh_vgg99"][["point_id", "xs_nb"]].rename(
                columns={"xs_nb": "partons_bh_vgg99"}
            ),
            on="point_id", how="left", validate="one_to_one",
        )
        merged["phi_mapping_mode"] = mode

        process_metrics = {}
        for process_name, col in [
            ("gv08", "partons_bh_gv08"),
            ("vgg99", "partons_bh_vgg99"),
        ]:
            gepard = merged["km15_bh"].to_numpy(float)
            partons = merged[col].to_numpy(float)
            valid = (
                np.isfinite(gepard)
                & np.isfinite(partons)
                & (gepard > 0.0)
                & (partons > 0.0)
            )
            ratio = np.full(len(merged), np.nan)
            ratio[valid] = partons[valid] / gepard[valid]
            merged[f"{process_name}_over_gepard"] = ratio

            rv = ratio[np.isfinite(ratio) & (ratio > 0.0)]
            if len(rv):
                log_r = np.log(rv)
                median = float(np.median(rv))
                log_std = float(np.std(log_r))
                median_abs_log_offset = float(abs(np.median(log_r)))
                broad_fraction = float(np.mean((rv > 0.5) & (rv < 2.0)))
                rmin = float(np.min(rv))
                rmax = float(np.max(rv))
            else:
                median = np.nan
                log_std = np.inf
                median_abs_log_offset = np.inf
                broad_fraction = 0.0
                rmin = np.nan
                rmax = np.nan
            #endif

            process_metrics[process_name] = {
                "valid_fraction": float(valid.mean()),
                "median": median,
                "log_std": log_std,
                "median_abs_log_offset": median_abs_log_offset,
                "broad_fraction": broad_fraction,
                "min": rmin,
                "max": rmax,
            }
        #endfor

        score = float(sum(
            process_metrics[p]["log_std"]
            + 0.25 * process_metrics[p]["median_abs_log_offset"]
            for p in ["gv08", "vgg99"]
        ))
        accepted = all(
            process_metrics[p]["valid_fraction"] == 1.0
            and process_metrics[p]["broad_fraction"] >= 0.90
            and process_metrics[p]["log_std"] <= 0.15
            for p in ["gv08", "vgg99"]
        )

        summary_rows.append({
            "mode": mode,
            "N": int(len(merged)),
            "score": score,
            "accepted": bool(accepted),
            "gv08_median_ratio": process_metrics["gv08"]["median"],
            "gv08_log_ratio_std": process_metrics["gv08"]["log_std"],
            "gv08_fraction_0p5_to_2": process_metrics["gv08"]["broad_fraction"],
            "gv08_min_ratio": process_metrics["gv08"]["min"],
            "gv08_max_ratio": process_metrics["gv08"]["max"],
            "vgg99_median_ratio": process_metrics["vgg99"]["median"],
            "vgg99_log_ratio_std": process_metrics["vgg99"]["log_std"],
            "vgg99_fraction_0p5_to_2": process_metrics["vgg99"]["broad_fraction"],
            "vgg99_min_ratio": process_metrics["vgg99"]["min"],
            "vgg99_max_ratio": process_metrics["vgg99"]["max"],
        })
        point_rows.append(merged)
    #endfor

    summary = pd.DataFrame(summary_rows).sort_values(
        ["accepted", "score"], ascending=[False, True]
    ).reset_index(drop=True)
    summary.to_csv(scan_root / "phi_mapping_summary.csv", index=False)
    pd.concat(point_rows, ignore_index=True).to_csv(
        scan_root / "phi_mapping_pointwise.csv", index=False
    )

    accepted = summary.loc[summary["accepted"].astype(bool)].copy()
    if len(accepted) == 0:
        print("\n[CLAS12 phi scan] no candidate mapping passed BH closure")
        print(summary.to_string(index=False))
        raise RuntimeError(
            "No CLAS12 phi mapping passed the BH-only PARTONS/Gepard closure. "
            f"Inspect {scan_root / 'phi_mapping_summary.csv'} and "
            f"{scan_root / 'phi_mapping_pointwise.csv'}."
        )
    #endif

    winner = str(accepted.iloc[0]["mode"])

    # Unpolarized BH is even under phi -> -phi, so a BH-only scan cannot
    # distinguish identity from negative.  Likewise, 180-minus and 180-plus
    # form the same BH-equivalent class.  Treat these exact/near-exact
    # degeneracies as physical symmetries rather than an ambiguity.
    equivalence_class = {
        "identity": "zero_shift_class",
        "negative": "zero_shift_class",
        "180-minus": "180_shift_class",
        "180-plus": "180_shift_class",
    }
    canonical_mode = {
        "zero_shift_class": "identity",
        "180_shift_class": "180-minus",
    }

    accepted = accepted.copy()
    accepted["bh_equivalence_class"] = accepted["mode"].map(equivalence_class)

    # Rank equivalence classes by their best score.
    class_best = (
        accepted.groupby("bh_equivalence_class", as_index=False)["score"]
        .min()
        .sort_values("score")
        .reset_index(drop=True)
    )
    winning_class = str(class_best.iloc[0]["bh_equivalence_class"])

    # If two DIFFERENT BH-equivalence classes are comparably good, that is a
    # genuine unresolved ambiguity and should still abort.
    if len(class_best) > 1:
        best = float(class_best.iloc[0]["score"])
        second = float(class_best.iloc[1]["score"])
        if np.isfinite(best) and np.isfinite(second) and second <= 1.05 * best:
            print("\n[CLAS12 phi scan] ambiguous BH-equivalence classes")
            print(accepted.to_string(index=False))
            raise RuntimeError(
                "CLAS12 BH-only phi scan cannot distinguish the zero-shift "
                "class {phi,-phi} from the 180-degree-shift class "
                "{180-phi,180+phi}. Inspect the pointwise scan before "
                "continuing."
            )
        #endif
    #endif

    winner = canonical_mode[winning_class]

    print("\n[CLAS12 phi scan] BH-only PARTONS/Gepard closure")
    print(summary.to_string(index=False))
    same_class_modes = accepted.loc[
        accepted["bh_equivalence_class"] == winning_class, "mode"
    ].tolist()
    if len(same_class_modes) > 1:
        print(
            "[CLAS12 phi scan] BH symmetry equivalence: "
            + ", ".join(same_class_modes)
            + f" are indistinguishable for XUU; using canonical '{winner}'."
        )
    #endif
    print(f"[CLAS12 phi scan] selected mapping: {winner}")
    print(f"[CLAS12 phi scan] diagnostics -> {scan_root}")
    return winner, summary
#enddef


def make_production_cache_tag(
        imported_mode: str,
        clas12_mode: str) -> str:
    """Deterministic cache namespace for the active dataset-specific mappings."""
    return (
        f"gepard-{imported_mode}"
        f"_direct-identity"
        f"_clas12-{clas12_mode}"
        f"_phi-v2"
    )
#enddef


def build_model_selection(
        points: pd.DataFrame,
        bh_gv08: pd.DataFrame,
        gk16: pd.DataFrame,
        bh_vgg99: pd.DataFrame,
        vgg99: pd.DataFrame,
        outdir: Path,
        thresholds: List[float]) -> pd.DataFrame:
    """Construct self-consistent purity quantities and all threshold selections."""
    out = points.copy()
    out = _attach_result(out, bh_gv08, "partons_bh_gv08")
    out = _attach_result(out, gk16, "partons_ep_gk16")
    out = _attach_result(out, bh_vgg99, "partons_bh_vgg99")
    out = _attach_result(out, vgg99, "partons_ep_vgg99")

    for col in [
        "partons_bh_gv08", "partons_ep_gk16",
        "partons_bh_vgg99", "partons_ep_vgg99",
    ]:
        if not np.all(np.isfinite(out[col].to_numpy(float))):
            raise RuntimeError(f"Nonfinite values remain in {col}")
        #endif
    #endfor

    out["delta_bh_gk16"] = np.abs(
        1.0
        - out["partons_bh_gv08"].to_numpy(float)
        / out["partons_ep_gk16"].to_numpy(float)
    )
    out["delta_bh_vgg99"] = np.abs(
        1.0
        - out["partons_bh_vgg99"].to_numpy(float)
        / out["partons_ep_vgg99"].to_numpy(float)
    )

    # Cross-code/process diagnostics only.  These are never used as correction
    # factors and never enter a model's purity definition.
    out["partons_gv08_over_gepard_bh"] = (
        out["partons_bh_gv08"].to_numpy(float)
        / out["km15_bh"].to_numpy(float)
    )
    out["partons_vgg99_over_gepard_bh"] = (
        out["partons_bh_vgg99"].to_numpy(float)
        / out["km15_bh"].to_numpy(float)
    )
    out["partons_vgg99_over_gv08_bh"] = (
        out["partons_bh_vgg99"].to_numpy(float)
        / out["partons_bh_gv08"].to_numpy(float)
    )

    deltas = out[
        ["delta_bh_km15", "delta_bh_vgg99", "delta_bh_gk16"]
    ].to_numpy(float)
    out["delta_bh_model_max"] = np.max(deltas, axis=1)
    out["delta_bh_model_min"] = np.min(deltas, axis=1)
    out["delta_bh_model_spread"] = (
        out["delta_bh_model_max"] - out["delta_bh_model_min"]
    )

    for threshold in thresholds:
        tag = f"{int(round(100.0 * threshold)):02d}pct"
        for model in ["km15", "vgg99", "gk16"]:
            out[f"pass_{tag}_{model}"] = (
                out[f"delta_bh_{model}"].to_numpy(float) <= threshold
            )
        #endfor

        out[f"pass_{tag}_all_models"] = (
            out[f"pass_{tag}_km15"]
            & out[f"pass_{tag}_vgg99"]
            & out[f"pass_{tag}_gk16"]
        )
        out[f"pass_{tag}_any_model"] = (
            out[f"pass_{tag}_km15"]
            | out[f"pass_{tag}_vgg99"]
            | out[f"pass_{tag}_gk16"]
        )
    #endfor

    outdir.mkdir(parents=True, exist_ok=True)
    all_path = outdir / "bh_model_selection_all_points.csv"
    out.to_csv(all_path, index=False)

    count_rows = []
    overlap_rows = []

    group_items = [("ALL", out)] + list(out.groupby("dataset", sort=False))
    for dataset, d in group_items:
        for threshold in thresholds:
            tag = f"{int(round(100.0 * threshold)):02d}pct"

            masks = {
                model: d[f"pass_{tag}_{model}"].to_numpy(bool)
                for model in ["km15", "vgg99", "gk16"]
            }
            all_mask = masks["km15"] & masks["vgg99"] & masks["gk16"]
            any_mask = masks["km15"] | masks["vgg99"] | masks["gk16"]

            count_rows.append(
                {
                    "dataset": dataset,
                    "threshold": threshold,
                    "n_total": len(d),
                    "n_km15": int(masks["km15"].sum()),
                    "n_vgg99": int(masks["vgg99"].sum()),
                    "n_gk16": int(masks["gk16"].sum()),
                    "n_all_models": int(all_mask.sum()),
                    "n_any_model": int(any_mask.sum()),
                    "common_over_union_fraction": (
                        float(all_mask.sum() / any_mask.sum())
                        if any_mask.sum() else np.nan
                    ),
                }
            )

            for a, b in [("km15", "vgg99"), ("km15", "gk16"), ("vgg99", "gk16")]:
                ma = masks[a]
                mb = masks[b]
                inter = ma & mb
                union = ma | mb
                overlap_rows.append(
                    {
                        "dataset": dataset,
                        "threshold": threshold,
                        "model_a": a,
                        "model_b": b,
                        "n_a": int(ma.sum()),
                        "n_b": int(mb.sum()),
                        "n_intersection": int(inter.sum()),
                        "n_union": int(union.sum()),
                        "jaccard": (
                            float(inter.sum() / union.sum())
                            if union.sum() else np.nan
                        ),
                    }
                )
            #endfor
        #endfor
    #endfor

    pd.DataFrame(count_rows).to_csv(
        outdir / "bh_model_selection_counts.csv", index=False
    )
    pd.DataFrame(overlap_rows).to_csv(
        outdir / "bh_model_selection_pairwise_overlap.csv", index=False
    )

    for model in ["km15", "vgg99", "gk16"]:
        selected = out.loc[out[f"pass_05pct_{model}"]].copy()
        selected.to_csv(
            outdir / f"bh_model_selection_{model}_05pct.csv", index=False
        )
        print(
            f"[selection] 5% {model:5s}: {len(selected)}/{len(out)} points -> "
            f"{outdir / f'bh_model_selection_{model}_05pct.csv'}"
        )
    #endfor

    consensus = out.loc[out["pass_05pct_all_models"]].copy()
    consensus.to_csv(
        outdir / "bh_model_selection_consensus_05pct.csv", index=False
    )

    print(f"[selection] all-point cache -> {all_path}")
    print(
        "[selection] note: all-model consensus is a diagnostic/stress-test sample; "
        "it is not the nominal selection prescription."
    )
    return out
#enddef


def _robust_upper(values: np.ndarray, minimum: float = 10.0) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) == 0:
        return minimum
    #endif
    return max(minimum, 1.10 * float(np.nanpercentile(finite, 99.0)))
#enddef


def save_diagnostic_plots(df: pd.DataFrame, outdir: Path) -> None:
    """Write per-dataset model-selection and BH-closure diagnostics."""
    import matplotlib.pyplot as plt

    models = ["km15", "vgg99", "gk16"]
    labels = {"km15": "KM15", "vgg99": "VGG99", "gk16": "GK16"}

    for dataset, d in df.groupby("dataset", sort=False):
        ddir = outdir / "plots" / dataset
        ddir.mkdir(parents=True, exist_ok=True)

        # 01: alternative-model purity deviations against KM15.
        fig, ax = plt.subplots(figsize=(7.2, 6.2))
        for model in ["vgg99", "gk16"]:
            ax.scatter(
                100.0 * d["delta_bh_km15"],
                100.0 * d[f"delta_bh_{model}"],
                s=12,
                alpha=0.65,
                label=labels[model],
            )
        #endfor

        all_delta_pct = 100.0 * d[
            ["delta_bh_km15", "delta_bh_vgg99", "delta_bh_gk16"]
        ].to_numpy(float)
        lim = _robust_upper(all_delta_pct, minimum=10.0)
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

        # 02: selected counts versus threshold.
        fig, ax = plt.subplots(figsize=(8.2, 5.4))
        threshold_pct = np.arange(1, 11)
        for model in models:
            counts = [
                int((d[f"delta_bh_{model}"] <= pct / 100.0).sum())
                for pct in threshold_pct
            ]
            ax.plot(threshold_pct, counts, marker="o", label=labels[model])
        #endfor
        common_counts = [
            int(
                (
                    (d["delta_bh_km15"] <= pct / 100.0)
                    & (d["delta_bh_vgg99"] <= pct / 100.0)
                    & (d["delta_bh_gk16"] <= pct / 100.0)
                ).sum()
            )
            for pct in threshold_pct
        ]
        ax.plot(
            threshold_pct,
            common_counts,
            marker="s",
            label="all-model intersection",
        )
        ax.set_xlabel("BH-purity threshold (%)")
        ax.set_ylabel("selected points")
        ax.set_title(f"{dataset}: selected sample versus purity threshold")
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(ddir / "02_selection_count_vs_threshold.png", dpi=180)
        plt.close(fig)

        # 03: common/union overlap fraction versus threshold.
        fig, ax = plt.subplots(figsize=(8.2, 5.4))
        fractions = []
        for pct in threshold_pct:
            cut = pct / 100.0
            masks = [
                d[f"delta_bh_{model}"].to_numpy(float) <= cut
                for model in models
            ]
            common = masks[0] & masks[1] & masks[2]
            union = masks[0] | masks[1] | masks[2]
            fractions.append(
                common.sum() / union.sum() if union.sum() else np.nan
            )
        #endfor
        ax.plot(threshold_pct, fractions, marker="o")
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel("BH-purity threshold (%)")
        ax.set_ylabel("three-model intersection / union")
        ax.set_title(f"{dataset}: model-selection overlap")
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(ddir / "03_common_overlap_fraction_vs_threshold.png", dpi=180)
        plt.close(fig)

        # 04: number of models accepting each point at nominal 5%.
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        agreement = (
            d[["pass_05pct_km15", "pass_05pct_vgg99", "pass_05pct_gk16"]]
            .astype(int)
            .sum(axis=1)
        )
        x = np.arange(4)
        counts = [(agreement == i).sum() for i in x]
        ax.bar(x, counts)
        ax.set_xticks(x)
        ax.set_xticklabels(["0 models", "1 model", "2 models", "3 models"])
        ax.set_ylabel("points")
        ax.set_title(f"{dataset}: 5% BH-selection agreement")
        fig.tight_layout()
        fig.savefig(ddir / "04_selection_agreement_05pct.png", dpi=180)
        plt.close(fig)

        # 05: full-point PARTONS/Gepard BH closure versus kinematics.
        kin_specs = [
            ("xB", r"$x_B$"),
            ("Q2", r"$Q^2$ (GeV$^2$)"),
            ("t_abs", r"$|t|$ (GeV$^2$)"),
            ("phi_deg", r"$\phi_{\rm Gepard}$ (deg)"),
        ]
        for index, (column, xlabel) in enumerate(kin_specs, start=1):
            fig, ax = plt.subplots(figsize=(7.5, 5.2))
            ax.scatter(
                d[column],
                d["partons_gv08_over_gepard_bh"],
                s=12,
                alpha=0.65,
                label="PARTONS GV08 / Gepard BH",
            )
            ax.scatter(
                d[column],
                d["partons_vgg99_over_gepard_bh"],
                s=12,
                alpha=0.65,
                marker="x",
                label="PARTONS VGG99 / Gepard BH",
            )
            ax.axhline(1.0, linestyle="--", linewidth=1.0)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r"$\sigma_{\rm BH}^{\rm PARTONS}/\sigma_{\rm BH}^{\rm Gepard}$")
            ax.set_title(f"{dataset}: BH implementation closure")
            ax.grid(alpha=0.2)
            ax.legend()
            fig.tight_layout()
            fig.savefig(
                ddir / f"05_{index:02d}_partons_gepard_bh_closure_vs_{column}.png",
                dpi=180,
            )
            plt.close(fig)
        #endfor

        # 06: VGG99-process BH versus GV08-process BH inside PARTONS.
        fig, ax = plt.subplots(figsize=(7.5, 5.2))
        ax.scatter(
            d["t_abs"],
            d["partons_vgg99_over_gv08_bh"],
            s=12,
            alpha=0.65,
        )
        ax.axhline(1.0, linestyle="--", linewidth=1.0)
        ax.set_xlabel(r"$|t|$ (GeV$^2$)")
        ax.set_ylabel(r"$\sigma_{\rm BH}^{\rm VGG99\ process}/\sigma_{\rm BH}^{\rm GV08\ process}$")
        ax.set_title(f"{dataset}: PARTONS BH process-module closure")
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(ddir / "06_partons_bh_process_module_closure.png", dpi=180)
        plt.close(fig)
    #endfor
#enddef



def choose_phi_validation_points(
        points: pd.DataFrame,
        dataset: str,
        n_points: int = DEFAULT_IMPORTED_SCAN_POINTS) -> pd.DataFrame:
    """Choose deterministic points spanning phi for one convention validation."""
    d = points.loc[points["dataset"].astype(str) == str(dataset)].copy()
    if len(d) == 0:
        raise RuntimeError(
            f"No dataset=={dataset!r} points are available for phi validation"
        )
    #endif

    n_points = max(4, min(int(n_points), len(d)))
    d = d.sort_values(["phi_deg", "xB", "Q2", "t_abs"]).reset_index(drop=True)
    indices = np.unique(
        np.rint(np.linspace(0, len(d) - 1, n_points)).astype(int)
    )
    return d.iloc[indices].copy().reset_index(drop=True)
#enddef


def _phi_scan_metrics(merged: pd.DataFrame, column: str) -> Dict[str, float]:
    """Summarize PARTONS/Gepard BH closure for one process and mapping."""
    gepard = merged["km15_bh"].to_numpy(float)
    partons = merged[column].to_numpy(float)
    valid = (
        np.isfinite(gepard)
        & np.isfinite(partons)
        & (gepard > 0.0)
        & (partons > 0.0)
    )
    ratio = np.full(len(merged), np.nan)
    ratio[valid] = partons[valid] / gepard[valid]
    rv = ratio[np.isfinite(ratio) & (ratio > 0.0)]

    if len(rv):
        log_r = np.log(rv)
        return {
            "valid_fraction": float(valid.mean()),
            "median_ratio": float(np.median(rv)),
            "log_ratio_std": float(np.std(log_r)),
            "median_abs_log_offset": float(abs(np.median(log_r))),
            "fraction_0p5_to_2": float(np.mean((rv > 0.5) & (rv < 2.0))),
            "ratio_min": float(np.min(rv)),
            "ratio_max": float(np.max(rv)),
        }
    #endif

    return {
        "valid_fraction": 0.0,
        "median_ratio": np.nan,
        "log_ratio_std": np.inf,
        "median_abs_log_offset": np.inf,
        "fraction_0p5_to_2": 0.0,
        "ratio_min": np.nan,
        "ratio_max": np.nan,
    }
#enddef


def run_imported_phi_convention_validation(
        points: pd.DataFrame,
        outdir: Path,
        sif: Path,
        project: Path,
        executable: str,
        workers: int,
        scan_points: int,
        force: bool = False) -> pd.DataFrame:
    """
    Validate dataset-specific azimuth mappings with BH-only code closure.

    The test intentionally does not use measured electroproduction cross
    sections or VGG/GK full-EP agreement.  Defurne (Gepard/BMK stored),
    Saylor (direct supplemental table), and Georges (direct Trento spreadsheet)
    are each tested under all four plausible phi mappings.

    Because unpolarized BH is even in phi, identity and negative form one
    equivalence class, while 180-minus and 180-plus form the other.  The
    expected class is therefore source-dependent rather than "imported"-data
    dependent.
    """
    dataset_specs = [
        ("halla_defurne2015", "Hall A Defurne 2015"),
        ("saylor2018", "CLAS6 Saylor 2018"),
        ("halla_georges2022", "Hall A Georges 2022"),
    ]
    expected_class = {
        "halla_defurne2015": "180_shift",
        "saylor2018": "zero_shift",
        "halla_georges2022": "zero_shift",
    }
    root = outdir / "phi_convention_validation"
    root.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    point_tables = []

    for dataset, label in dataset_specs:
        available = points.loc[
            points["dataset"].astype(str) == dataset
        ].copy()
        if len(available) == 0:
            print(
                f"[phi validation] {label}: not present in kinematics cache; skipping"
            )
            continue
        #endif

        base = choose_phi_validation_points(
            points, dataset=dataset, n_points=scan_points
        )
        droot = root / dataset
        droot.mkdir(parents=True, exist_ok=True)
        base[
            [
                "point_id", "dataset", "source_row", "xB", "Q2",
                "t_abs", "phi_deg", "ebeam", "km15_bh",
            ]
        ].to_csv(droot / "scan_points.csv", index=False)

        for mode in PHI_MAPPING_MODES:
            test = base.copy()
            test["phi_partons_deg"] = transform_phi_to_partons(
                test["phi_deg"].to_numpy(float), mode
            )
            test["phi_mapping_mode"] = mode

            mode_root = droot / mode
            results = {}
            for calculation in ["bh_gv08", "bh_vgg99"]:
                results[calculation] = run_partons_calculation(
                    points=test,
                    calculation=calculation,
                    outdir=mode_root,
                    sif=sif,
                    project=project,
                    executable=executable,
                    workers=min(max(1, workers), 2),
                    chunk_size=max(1, len(test)),
                    force=force,
                    cache_tag=f"{dataset}_phi_{mode}",
                )
            #endfor

            merged = test[
                [
                    "point_id", "dataset", "source_row", "xB", "Q2",
                    "t_abs", "phi_deg", "phi_partons_deg",
                    "phi_mapping_mode", "ebeam", "km15_bh",
                ]
            ].copy()
            merged = merged.merge(
                results["bh_gv08"][["point_id", "xs_nb"]].rename(
                    columns={"xs_nb": "partons_bh_gv08"}
                ),
                on="point_id", how="left", validate="one_to_one",
            )
            merged = merged.merge(
                results["bh_vgg99"][["point_id", "xs_nb"]].rename(
                    columns={"xs_nb": "partons_bh_vgg99"}
                ),
                on="point_id", how="left", validate="one_to_one",
            )

            for process_name, col in [
                ("gv08", "partons_bh_gv08"),
                ("vgg99", "partons_bh_vgg99"),
            ]:
                metrics = _phi_scan_metrics(merged, col)
                summary_rows.append({
                    "dataset": dataset,
                    "dataset_label": label,
                    "mapping": mode,
                    "process": process_name,
                    "N": int(len(merged)),
                    **metrics,
                })
                ratio = (
                    merged[col].to_numpy(float)
                    / merged["km15_bh"].to_numpy(float)
                )
                merged[f"{process_name}_over_gepard"] = ratio
            #endfor

            point_tables.append(merged)
        #endfor
    #endfor

    if not summary_rows:
        raise RuntimeError(
            "Neither Defurne nor Georges is present in the supplied kinematics cache"
        )
    #endif

    summary = pd.DataFrame(summary_rows)
    points_out = pd.concat(point_tables, ignore_index=True)
    summary.to_csv(root / "phi_mapping_summary.csv", index=False)
    points_out.to_csv(root / "phi_mapping_pointwise.csv", index=False)

    # One publication-audit canvas: each row is a dataset; left is the explicit
    # angle mapping, right is BH code closure for all candidate mappings.
    present_specs = [
        spec for spec in dataset_specs
        if spec[0] in set(points_out["dataset"].astype(str))
    ]
    fig, axes = plt.subplots(
        len(present_specs), 2,
        figsize=(13.5, 4.8 * len(present_specs)),
        squeeze=False,
    )

    for irow, (dataset, label) in enumerate(present_specs):
        ax_map = axes[irow, 0]
        ax_ratio = axes[irow, 1]

        # Mapping visualization using the documented candidate.
        d = points_out.loc[
            (points_out["dataset"].astype(str) == dataset)
            & (points_out["phi_mapping_mode"].astype(str) == "180-minus")
        ].copy()
        d = d.drop_duplicates("point_id").sort_values("phi_deg")
        ax_map.plot(
            d["phi_deg"].to_numpy(float),
            d["phi_partons_deg"].to_numpy(float),
            marker="o",
            linestyle="none",
        )
        ax_map.set_xlabel(r"Stored/published $\phi$ (deg)")
        ax_map.set_ylabel(r"Candidate PARTONS $\phi$ (deg)")
        ax_map.set_title(f"{label}: 180° − φ mapping")
        ax_map.grid(alpha=0.2)

        # GV08 BH closure is sufficient to expose the mapping class; VGG99 BH
        # is retained in the CSV as an independent implementation check.
        for mode in PHI_MAPPING_MODES:
            m = points_out.loc[
                (points_out["dataset"].astype(str) == dataset)
                & (points_out["phi_mapping_mode"].astype(str) == mode)
            ].copy().sort_values("phi_deg")
            ax_ratio.plot(
                m["phi_deg"].to_numpy(float),
                m["gv08_over_gepard"].to_numpy(float),
                marker="o",
                linewidth=1.0,
                label=mode,
            )
        #endfor
        ax_ratio.axhline(1.0, linewidth=0.8)
        ax_ratio.set_xlabel(r"Stored/published $\phi$ (deg)")
        ax_ratio.set_ylabel("PARTONS GV08 BH / Gepard BH")
        ax_ratio.set_yscale("log")
        ax_ratio.set_title(f"{label}: BH-only convention closure")
        ax_ratio.grid(alpha=0.2, which="both")
        ax_ratio.legend(ncol=2)
    #endfor

    fig.suptitle(
        "Dataset-specific azimuth-convention validation using BH-only code closure",
        y=0.995,
    )
    fig.subplots_adjust(
        top=0.94, bottom=0.08, left=0.08, right=0.98,
        hspace=0.34, wspace=0.25,
    )
    fig.savefig(root / "01_halla_phi_convention_validation.png", dpi=190)
    plt.close(fig)

    # Print a compact class-level verdict without choosing a convention from
    # electroproduction chi2.
    equivalence_class = {
        "identity": "zero_shift",
        "negative": "zero_shift",
        "180-minus": "180_shift",
        "180-plus": "180_shift",
    }
    print("\n[phi validation] Dataset-specific BH-only convention summary")
    for dataset, label in present_specs:
        dsum = summary.loc[
            (summary["dataset"].astype(str) == dataset)
            & (summary["process"].astype(str) == "gv08")
        ].copy()
        dsum["class"] = dsum["mapping"].map(equivalence_class)
        class_scores = (
            dsum.assign(
                score=dsum["log_ratio_std"].to_numpy(float)
                + 0.25 * dsum["median_abs_log_offset"].to_numpy(float)
            )
            .groupby("class", as_index=False)["score"].min()
            .sort_values("score")
        )
        best_class = str(class_scores.iloc[0]["class"])
        expected = expected_class[dataset]
        expected_supported = best_class == expected
        canonical_mode = (
            "180-minus" if expected == "180_shift" else "identity"
        )
        canonical = dsum.loc[
            dsum["mapping"].astype(str) == canonical_mode
        ].iloc[0]
        print(
            f"  {label}: source={FIXED_DATASET_PHI_SOURCES[dataset]}; "
            f"best BH class={best_class}; expected={expected}; "
            f"{canonical_mode} GV08 median="
            f"{float(canonical['median_ratio']):.5f}, "
            f"log-std={float(canonical['log_ratio_std']):.5f}; "
            f"expected class supported={expected_supported}"
        )
    #endfor
    print(f"[phi validation] outputs -> {root}")
    return summary
#enddef


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="KM15/VGG99/GK16 BH-purity model-selection evaluator"
    )
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument(
        "--kinematics-cache",
        default=DEFAULT_KINEMATICS_CACHE,
        help=(
            "Exact experimental-point cache exported by extract_emff_from_dvcs_bh.py; "
            f"default: {DEFAULT_KINEMATICS_CACHE}"
        ),
    )
    p.add_argument(
        "--run-partons",
        action="store_true",
        help="Run/reuse PARTONS caches and construct model selections",
    )
    p.add_argument(
        "--make-partons-xml",
        action="store_true",
        help="Only generate all PARTONS XML scenarios and manifests",
    )
    p.add_argument("--partons-sif", default=DEFAULT_PARTONS_SIF)
    p.add_argument("--partons-project", default=DEFAULT_PARTONS_PROJECT)
    p.add_argument("--partons-executable", default=DEFAULT_PARTONS_EXECUTABLE)
    p.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_PARTONS_WORKERS,
        help=(
            "Independent Apptainer/PARTONS OS processes; default "
            f"{DEFAULT_PARTONS_WORKERS}"
        ),
    )
    p.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_PARTONS_CHUNK_SIZE,
        help=(
            "Sequential points per PARTONS process invocation; default "
            f"{DEFAULT_PARTONS_CHUNK_SIZE}"
        ),
    )
    p.add_argument(
        "--force-partons",
        action="store_true",
        help="Ignore complete PARTONS result caches and rerun all calculations",
    )
    p.add_argument(
        "--clas12-phi-mapping",
        choices=("auto",) + PHI_MAPPING_MODES,
        default=DEFAULT_CLAS12_PHI_MAPPING,
        help=(
            "CLAS12/pass1/pass2 phi mapping into PARTONS. Default 'auto' runs a "
            "small BH-only PARTONS/Gepard closure scan over identity, negative, "
            "180-minus, and 180-plus before production."
        ),
    )
    p.add_argument(
        "--imported-phi-mapping",
        choices=PHI_MAPPING_MODES,
        default=DEFAULT_IMPORTED_PHI_MAPPING,
        help=(
            "Backwards-compatible phi mapping override for Gepard-bundled "
            "Jo 2015 and Defurne 2015 only. Georges and Saylor are direct "
            "publication inputs and never inherit this option."
        ),
    )
    p.add_argument(
        "--clas12-phi-scan-points",
        type=int,
        default=DEFAULT_CLAS12_SCAN_POINTS,
        help=(
            "Representative CLAS12 points used by the BH-only automatic phi "
            f"scan; default {DEFAULT_CLAS12_SCAN_POINTS}."
        ),
    )
    p.add_argument(
        "--scan-clas12-phi-only",
        action="store_true",
        help=(
            "Run the CLAS12 BH-only phi-mapping scan, print/write diagnostics, "
            "and exit before the full model calculation."
        ),
    )
    p.add_argument(
        "--validate-phi-conventions",
        action="store_true",
        help=(
            "Run BH-only phi-convention validation for Defurne 2015, "
            "Saylor 2018, and Georges 2022 under all four candidate mappings."
        ),
    )
    p.add_argument(
        "--imported-phi-scan-points",
        type=int,
        default=DEFAULT_IMPORTED_SCAN_POINTS,
        help=(
            "Representative points per non-CLAS12 dataset used by "
            "--validate-phi-conventions; default "
            f"{DEFAULT_IMPORTED_SCAN_POINTS}."
        ),
    )
    p.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=DEFAULT_THRESHOLDS,
        help="BH-purity thresholds as fractions; default 0.01 ... 0.10",
    )
    return p
#enddef


def main(argv: List[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    outdir = Path(args.outdir)
    cache = Path(args.kinematics_cache)

    print(f"[storage] PARTONS/model-prediction output root: {outdir}")
    print(f"[storage] kinematic input cache: {cache}")

    if not cache.exists():
        print(
            "No kinematic cache exists:\n"
            f"  {cache}\n"
            "Run extract_emff_from_dvcs_bh.py first so it exports the exact "
            "KM15-evaluated experimental points."
        )
        return 2
    #endif

    points = load_points_from_cache(cache)
    print(
        f"[kinematics] loaded {len(points)} exact experimental points from {cache}"
    )

    thresholds = sorted(set(float(x) for x in args.thresholds))
    if 0.05 not in thresholds:
        raise ValueError(
            "--thresholds must include 0.05 because nominal 5% products are written"
        )
    #endif

    calculations = ["bh_gv08", "gk16", "bh_vgg99", "vgg99"]

    needs_partons_runtime = (
        args.run_partons
        or args.scan_clas12_phi_only
        or args.validate_phi_conventions
        or args.clas12_phi_mapping == "auto"
    )
    if needs_partons_runtime:
        sif = Path(args.partons_sif).expanduser().resolve()
        project = Path(args.partons_project).expanduser().resolve()
        locate_partons(sif, project, args.partons_executable)
    else:
        sif = Path(args.partons_sif).expanduser().resolve()
        project = Path(args.partons_project).expanduser().resolve()
    #endif

    if args.clas12_phi_mapping == "auto":
        clas12_phi_mapping, _ = run_clas12_phi_mapping_scan(
            points=points,
            outdir=outdir,
            sif=sif,
            project=project,
            executable=args.partons_executable,
            workers=args.workers,
            scan_points=args.clas12_phi_scan_points,
            force=args.force_partons,
        )
    else:
        clas12_phi_mapping = str(args.clas12_phi_mapping)
        print(
            "[CLAS12 phi] using explicitly requested mapping: "
            f"{clas12_phi_mapping}"
        )
    #endif

    if args.validate_phi_conventions:
        run_imported_phi_convention_validation(
            points=points,
            outdir=outdir,
            sif=sif,
            project=project,
            executable=args.partons_executable,
            workers=args.workers,
            scan_points=args.imported_phi_scan_points,
            force=args.force_partons,
        )
    #endif

    if args.scan_clas12_phi_only:
        print("[done] CLAS12 BH-only phi-mapping scan complete")
        return 0
    #endif

    if args.validate_phi_conventions and not args.run_partons:
        print("[done] dataset-specific phi-convention validation complete")
        return 0
    #endif

    points = apply_dataset_phi_mappings(
        points,
        imported_mode=args.imported_phi_mapping,
        clas12_mode=clas12_phi_mapping,
    )
    production_cache_tag = make_production_cache_tag(
        args.imported_phi_mapping, clas12_phi_mapping
    )
    print("[kinematics] PARTONS phi mappings by source:")
    for key in sorted(set(points["dataset"].astype(str))):
        d = points.loc[points["dataset"].astype(str) == key]
        mode = str(d.iloc[0]["phi_mapping_mode"])
        source = str(d.iloc[0]["phi_source_representation"])
        print(f"  {key:22s}: {mode:10s} <- {source}")
    #endfor
    print(f"[PARTONS cache namespace] {production_cache_tag}")

    if args.make_partons_xml:
        xml_preview_outdir = outdir / "xml_preview" / production_cache_tag
        for calculation in calculations:
            write_partons_xml_chunks(
                points,
                calculation,
                xml_preview_outdir,
                args.chunk_size,
            )
        #endfor
        if not args.run_partons:
            return 0
        #endif
    #endif

    if not args.run_partons:
        print(
            "Nothing to run after phi validation. Use --run-partons for the "
            "full cached PARTONS evaluation, --scan-clas12-phi-only for the "
            "CLAS12 BH-only convention test, --validate-phi-conventions for "
            "Defurne/Saylor/Georges BH-only validation, or --make-partons-xml to inspect XML."
        )
        return 0
    #endif

    print(
        f"[PARTONS performance] workers={args.workers}, "
        f"chunk_size={args.chunk_size}; "
        "numerical-library threads constrained to 1 per worker"
    )

    if args.workers < 1:
        raise ValueError("--workers must be >= 1")
    #endif

    results = {}
    for calculation in calculations:
        results[calculation] = run_partons_calculation(
            points=points,
            calculation=calculation,
            outdir=outdir,
            sif=sif,
            project=project,
            executable=args.partons_executable,
            workers=args.workers,
            chunk_size=args.chunk_size,
            force=args.force_partons,
            cache_tag=production_cache_tag,
        )
    #endfor

    mapping_rows = []
    for key in sorted(set(points["dataset"].astype(str))):
        d = points.loc[points["dataset"].astype(str) == key]
        mapping_rows.append({
            "dataset": key,
            "phi_mapping_mode": str(d.iloc[0]["phi_mapping_mode"]),
            "phi_source_representation": str(
                d.iloc[0]["phi_source_representation"]
            ),
            "production_cache_tag": production_cache_tag,
            "mapping_source": (
                "BH-only automatic CLAS12 scan"
                if key in {"pass1", "pass2"}
                and args.clas12_phi_mapping == "auto"
                else "dataset-specific source convention"
            ),
        })
    #endfor
    pd.DataFrame(mapping_rows).to_csv(
        outdir / "partons_phi_mapping_used.csv", index=False
    )

    merged = build_model_selection(
        points=points,
        bh_gv08=results["bh_gv08"],
        gk16=results["gk16"],
        bh_vgg99=results["bh_vgg99"],
        vgg99=results["vgg99"],
        outdir=outdir,
        thresholds=thresholds,
    )
    save_diagnostic_plots(merged, outdir)

    print("[done] BH-purity model-selection evaluation complete")
    return 0
#enddef


if __name__ == "__main__":
    raise SystemExit(main())
#endif
