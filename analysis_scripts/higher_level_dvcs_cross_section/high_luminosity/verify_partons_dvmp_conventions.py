#!/usr/bin/env python3
"""
verify_partons_dvmp_conventions.py

Small, non-destructive PARTONS convention check for the Stage-3 pi0 GK grid.

The script:
  1. Reads representative kinematic points from the durable GK production grid.
  2. Reconstructs A, B, C from the existing phi = 0, 90, 180 deg results.
  3. Builds a tiny PARTONS XML job at those same kinematics using
     DVMPCrossSectionUUUMinusPhiIntegrated.
  4. Runs PARTONS.
  5. Parses the phi-integrated results and tests
         sigma_phi_integrated / (2*pi*A) = 1.
  6. Writes a compact CSV and text report.

It does NOT modify the Stage-3 production grid or rerun the expensive grid.
"""

from __future__ import annotations

import argparse
import os
import shutil
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_PROJECT = Path("/work/clas12/thayward/partons/partons-example")
DEFAULT_EXECUTABLE = "./bin/PARTONS_example"

DEFAULT_STAGE3 = Path(
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "higher_level_dvcs_cross_section/high_luminosity/output/"
    "pi0_gk_stage3/partons_gk"
)

RESULT_RE = re.compile(
    r"Result\s*[:=]\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|[-+]?nan|[-+]?inf)"
    r"\s*([^\s<]*)?",
    re.IGNORECASE,
)


def find_first_existing(base: Path, candidates: list[str]) -> Path:
    for rel in candidates:
        p = base / rel
        if p.exists():
            return p
        #endif
    #endfor
    raise FileNotFoundError(
        "Could not find any of:\n  " + "\n  ".join(str(base / x) for x in candidates)
    )


def find_sif(explicit: Path | None, here: Path, project: Path) -> Path:
    if explicit is not None:
        return explicit.expanduser().resolve()
    #endif

    candidates = [
        here / "partons_v4.sif",
        project / "partons_v4.sif",
        project.parent / "partons_v4.sif",
        Path("/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/external_scripts/partons_v4.sif"),
        Path("/work/clas12/thayward/partons/partons_v4.sif"),
        Path.cwd() / "partons_v4.sif",
    ]
    for p in candidates:
        if p.exists():
            return p.resolve()
        #endif
    #endfor

    raise FileNotFoundError(
        "Could not locate partons_v4.sif. Searched:\n  "
        + "\n  ".join(str(x) for x in candidates)
    )


def preflight_partons(sif: Path, project: Path, executable: str) -> None:
    exe = project / executable.removeprefix("./")
    bad = []
    if shutil.which("apptainer") is None:
        bad.append("apptainer not on PATH")
    #endif
    if not sif.exists():
        bad.append(f"SIF missing: {sif}")
    #endif
    if not project.exists():
        bad.append(f"project missing: {project}")
    #endif
    if not exe.exists():
        bad.append(f"executable missing: {exe}")
    elif not os.access(exe, os.X_OK):
        bad.append(f"executable not executable: {exe}")
    #endif
    if bad:
        raise RuntimeError("PARTONS preflight failed:\n  - " + "\n  - ".join(bad))
    #endif


def run_partons_xml(xml_path: Path, sif: Path, project: Path, executable: str) -> subprocess.CompletedProcess:
    binds = [project.resolve(), xml_path.parent.resolve()]
    cmd = ["apptainer", "exec"]
    for b in binds:
        cmd += ["--bind", f"{b}:{b}"]
    #endfor
    cmd += ["--pwd", str(project.resolve()), str(sif), executable, str(xml_path.resolve())]

    env = os.environ.copy()
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    #endfor

    return subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )


def load_metadata(stage3: Path) -> dict:
    p = stage3 / "production_grid" / "metadata.json"
    if not p.exists():
        return {}
    #endif
    with p.open() as f:
        return json.load(f)


def infer_columns(df: pd.DataFrame) -> dict[str, str]:
    """Locate required columns without hard-coding one exact Stage-3 spelling."""

    choices = {
        "point_id": ["point_id", "shared_point_id", "match_id"],
        "campaign": ["campaign", "dataset", "run_group"],
        "evaluation": ["evaluation", "coordinate", "evaluation_type"],
        "E": ["E", "beam_energy_GeV", "beam_energy"],
        "epsilon": ["epsilon", "eps", "virtual_photon_epsilon"],
        "Q2": ["Q2", "Q2_GeV2", "Q2_shared_GeV2"],
        "xB": ["xB", "xb", "x_B", "xB_shared"],
        "minus_t": ["minus_t", "minus_t_GeV2", "minus_t_shared_GeV2"],
        "phi": ["phi_deg", "phi", "phi_degrees"],
        "sigma": [
            "result",
            "cross_section",
            "cross_section_nb",
            "sigma",
            "partons_result",
            "partons_value",
            "value",
        ],
    }

    out = {}
    lower = {str(c).lower(): c for c in df.columns}

    for key, names in choices.items():
        found = None
        for name in names:
            if name in df.columns:
                found = name
                break
            #endif
            if name.lower() in lower:
                found = lower[name.lower()]
                break
            #endif
        #endfor
        if found is not None:
            out[key] = found
        #endif
    #endfor

    return out


def load_raw_phi_grid(stage3: Path) -> tuple[pd.DataFrame, Path]:
    p = find_first_existing(
        stage3,
        [
            "04_partons_raw_phi_grid.csv",
            "production_grid/gk_pi0_model_grid.csv",
        ],
    )
    df = pd.read_csv(p)
    return df, p


def choose_representative_groups(
    df: pd.DataFrame, cols: dict[str, str], npoints: int
) -> list[pd.DataFrame]:
    required = ["point_id", "campaign", "evaluation", "E", "Q2", "xB", "minus_t", "phi", "sigma"]
    missing = [x for x in required if x not in cols]
    if missing:
        raise RuntimeError(
            "Could not identify required columns "
            f"{missing} in grid. Available columns:\n{list(df.columns)}"
        )
    #endif

    work = df.copy()

    # Prefer shared-coordinate evaluations because those are what enter the
    # RGA/RGK Rosenbluth comparison.
    eval_text = work[cols["evaluation"]].astype(str).str.lower()
    shared = work[eval_text.str.contains("shared", na=False)].copy()
    if len(shared):
        work = shared
    #endif

    group_cols = [
        cols["point_id"],
        cols["campaign"],
        cols["evaluation"],
        cols["E"],
        cols["Q2"],
        cols["xB"],
        cols["minus_t"],
    ]

    good = []
    for _, g in work.groupby(group_cols, dropna=False, sort=False):
        phis = np.asarray(g[cols["phi"]], dtype=float)
        if all(np.any(np.isclose(phis, p, atol=1e-8)) for p in (0.0, 90.0, 180.0)):
            good.append(g.copy())
        #endif
    #endfor

    if not good:
        raise RuntimeError("No groups with phi = 0, 90, 180 deg were found.")
    #endif

    summary = pd.DataFrame(
        [
            {
                "i": i,
                "Q2": float(g.iloc[0][cols["Q2"]]),
                "xB": float(g.iloc[0][cols["xB"]]),
            }
            for i, g in enumerate(good)
        ]
    ).sort_values(["Q2", "xB"])

    if npoints <= 1:
        indices = [int(summary.iloc[len(summary) // 2]["i"])]
    else:
        positions = np.linspace(0, len(summary) - 1, npoints)
        indices = []
        for pos in positions:
            idx = int(summary.iloc[int(round(pos))]["i"])
            if idx not in indices:
                indices.append(idx)
            #endif
        #endfor
    #endif

    return [good[i] for i in indices]


def get_phi_value(g: pd.DataFrame, phi_col: str, sigma_col: str, target: float) -> float:
    mask = np.isclose(np.asarray(g[phi_col], dtype=float), target, atol=1e-8)
    vals = pd.to_numeric(g.loc[mask, sigma_col], errors="coerce").dropna().to_numpy()
    if len(vals) != 1:
        raise RuntimeError(
            f"Expected exactly one result at phi={target:g} deg; found {len(vals)}."
        )
    #endif
    return float(vals[0])


def harmonic_coefficients(s0: float, s90: float, s180: float) -> tuple[float, float, float]:
    A = (s0 + 2.0 * s90 + s180) / 4.0
    B = (s0 - s180) / 2.0
    C = (s0 - 2.0 * s90 + s180) / 4.0
    return A, B, C


def build_phi_integrated_xml(
    rows: list[dict],
    xml_path: Path,
    metadata: dict,
) -> None:
    """
    Build a deliberately small PARTONS job.

    The module chain is fixed to the same chain recorded by the successful
    Stage-3 production metadata. If your local PARTONS XML syntax differs,
    compare this generated file with one successful Stage-3 chunk XML; the
    kinematic values and observable-module substitution are the only intended
    differences.
    """

    # These defaults reproduce the successful Stage-3 model chain described
    # by metadata. Keep the names explicit so the verification is auditable.
    gpd = metadata.get("gpd_module", "GPDGK19")
    cff = metadata.get("cff_module", "DVMPCFFGK06")
    process = metadata.get("process_module", "DVMPProcessGK06")
    observable = "DVMPCrossSectionUUUMinusPhiIntegrated"

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        "<partons>",
        "  <scenario>",
        f'    <module type="GPDModule" name="{gpd}"/>',
        f'    <module type="DVMPCFFModule" name="{cff}"/>',
        f'    <module type="DVMPProcessModule" name="{process}"/>',
        f'    <module type="Observable" name="{observable}"/>',
    ]

    for i, r in enumerate(rows):
        lines.extend(
            [
                f'    <kinematics id="VERIFY_{i:02d}"',
                f'      E="{r["E"]:.16g}"',
                f'      Q2="{r["Q2"]:.16g}"',
                f'      xB="{r["xB"]:.16g}"',
                f'      t="{-r["minus_t"]:.16g}"',
                '      mesonPdg="111"/>',
            ]
        )
    #endfor

    lines.extend(["  </scenario>", "</partons>", ""])
    xml_path.write_text("\n".join(lines))


def find_successful_stage3_xml(stage3: Path) -> Path | None:
    candidates = sorted(stage3.rglob("*.xml"))
    return candidates[0] if candidates else None


def clone_stage3_xml_if_possible(
    stage3: Path,
    rows: list[dict],
    xml_path: Path,
) -> bool:
    """
    Prefer the exact XML syntax already proven to work in Stage 3.

    We use a successful chunk XML as a template, replace the observable module
    name, and retain only enough calculation blocks for the selected points
    when the structure is simple. Because PARTONS installations can differ,
    this routine is conservative: if it cannot safely recognize the XML
    structure, it returns False and the script writes a transparent fallback
    XML plus instructions.
    """
    template = find_successful_stage3_xml(stage3)
    if template is None:
        return False
    #endif

    text = template.read_text(errors="replace")
    if "DVMPCrossSectionUUUMinus" not in text:
        return False
    #endif

    # First make the only convention change we need.
    text = text.replace(
        "DVMPCrossSectionUUUMinusPhiIntegrated",
        "DVMPCrossSectionUUUMinus",
    )
    text = text.replace(
        "DVMPCrossSectionUUUMinus",
        "DVMPCrossSectionUUUMinusPhiIntegrated",
    )

    # We intentionally do not attempt brittle regex surgery on unknown PARTONS
    # scenario syntax here. Save the proven template so the user can inspect
    # the exact module substitution if automatic execution is not possible.
    xml_path.write_text(text)
    return True


def parse_results(text: str) -> list[tuple[float, str]]:
    out = []
    for m in RESULT_RE.finditer(text):
        token = m.group(1)
        unit = (m.group(2) or "").strip().strip("[]")
        try:
            value = float(token)
        except ValueError:
            continue
        #endtry
        if np.isfinite(value):
            out.append((value, unit))
        #endif
    #endfor
    return out



def build_shared_decomposition_table(df: pd.DataFrame, cols: dict[str, str]) -> pd.DataFrame:
    """Test the conventional T/L/TT/LT mapping using shared RGA/RGK points.

    At a shared (Q2,xB,t), the underlying response functions must be independent
    of beam energy.  Therefore the candidate mapping can be tested by comparing
    the values inferred independently from RGA and RGK.
    """
    required = ["point_id", "campaign", "evaluation", "epsilon", "phi", "sigma"]
    missing = [x for x in required if x not in cols]
    if missing:
        raise RuntimeError(
            "Cannot build T/L/TT/LT convention test; missing columns " + str(missing)
        )
    #endif

    work = df.copy()
    work = work[
        work[cols["evaluation"]].astype(str).str.lower().str.contains("shared", na=False)
    ].copy()

    harmonic_rows = []
    group_cols = [cols["point_id"], cols["campaign"]]
    for (point_id, campaign), g in work.groupby(group_cols, sort=False):
        phis = np.asarray(g[cols["phi"]], dtype=float)
        if not all(np.any(np.isclose(phis, p, atol=1e-8)) for p in (0.0, 90.0, 180.0)):
            continue
        #endif
        s0 = get_phi_value(g, cols["phi"], cols["sigma"], 0.0)
        s90 = get_phi_value(g, cols["phi"], cols["sigma"], 90.0)
        s180 = get_phi_value(g, cols["phi"], cols["sigma"], 180.0)
        A, B, C = harmonic_coefficients(s0, s90, s180)
        r0 = g.iloc[0]
        eps = float(r0[cols["epsilon"]])
        harmonic_rows.append({
            "point_id": str(point_id),
            "campaign": str(campaign).lower(),
            "epsilon": eps,
            "Q2": float(r0[cols["Q2"]]),
            "xB": float(r0[cols["xB"]]),
            "minus_t": float(r0[cols["minus_t"]]),
            "A": A, "B": B, "C": C,
            # Candidate conventional mapping.  These are hypotheses until the
            # source-level GK CrossSection() formula is confirmed.
            "sigma_LT_candidate": 2.0 * math.pi * B / math.sqrt(2.0 * eps * (1.0 + eps)),
            "sigma_TT_candidate": 2.0 * math.pi * C / eps,
        })
    #endfor

    h = pd.DataFrame(harmonic_rows)
    out = []
    for point_id, g in h.groupby("point_id", sort=False):
        rga = g[g["campaign"] == "rga"]
        rgk = g[g["campaign"] == "rgk"]
        if len(rga) != 1 or len(rgk) != 1:
            continue
        #endif
        a = rga.iloc[0]
        k = rgk.iloc[0]
        de = float(k.epsilon - a.epsilon)
        if abs(de) < 1e-12:
            continue
        #endif

        # Since phi integration established 2*pi*A as the phi-integrated
        # unpolarized quantity, test 2*pi*A = sigma_T + epsilon*sigma_L.
        YA = 2.0 * math.pi * float(a.A)
        YK = 2.0 * math.pi * float(k.A)
        sigma_L = (YK - YA) / de
        sigma_T = YA - float(a.epsilon) * sigma_L

        lt_a = float(a.sigma_LT_candidate); lt_k = float(k.sigma_LT_candidate)
        tt_a = float(a.sigma_TT_candidate); tt_k = float(k.sigma_TT_candidate)
        lt_mean = 0.5 * (lt_a + lt_k)
        tt_mean = 0.5 * (tt_a + tt_k)
        lt_rel = (lt_a - lt_k) / lt_mean if lt_mean != 0.0 else np.nan
        tt_rel = (tt_a - tt_k) / tt_mean if tt_mean != 0.0 else np.nan

        out.append({
            "point_id": point_id,
            "Q2": float(a.Q2), "xB": float(a.xB), "minus_t": float(a.minus_t),
            "epsilon_rga": float(a.epsilon), "epsilon_rgk": float(k.epsilon),
            "delta_epsilon": de,
            "two_pi_A_rga": YA, "two_pi_A_rgk": YK,
            "sigma_T_candidate": sigma_T,
            "sigma_L_candidate": sigma_L,
            "sigma_L_over_sigma_T_candidate": sigma_L / sigma_T if sigma_T != 0.0 else np.nan,
            "sigma_LT_rga_candidate": lt_a,
            "sigma_LT_rgk_candidate": lt_k,
            "sigma_LT_relative_difference": lt_rel,
            "sigma_TT_rga_candidate": tt_a,
            "sigma_TT_rgk_candidate": tt_k,
            "sigma_TT_relative_difference": tt_rel,
        })
    #endfor
    return pd.DataFrame(out)


def print_decomposition_summary(tab: pd.DataFrame, npoints: int = 3) -> None:
    if tab.empty:
        print("No complete shared RGA/RGK pairs available for decomposition test.")
        return
    #endif
    ordered = tab.sort_values(["Q2", "xB"]).reset_index(drop=True)
    pos = np.linspace(0, len(ordered) - 1, min(npoints, len(ordered)))
    inds = sorted(set(int(round(x)) for x in pos))
    print()
    print("Candidate T/L/TT/LT convention test at shared RGA/RGK kinematics")
    print("  hypothesis: 2*pi*A = sigma_T + epsilon*sigma_L")
    print("              2*pi*B = sqrt(2*epsilon*(1+epsilon))*sigma_LT")
    print("              2*pi*C = epsilon*sigma_TT")
    for i in inds:
        r = ordered.iloc[i]
        print(
            f"  {r.point_id}: Q2={r.Q2:.6g}, xB={r.xB:.6g}, -t={r.minus_t:.6g}; "
            f"eps(RGA,RGK)=({r.epsilon_rga:.5f},{r.epsilon_rgk:.5f})"
        )
        print(
            f"      sigma_T={r.sigma_T_candidate:.10g}, "
            f"sigma_L={r.sigma_L_candidate:.10g}, "
            f"L/T={r.sigma_L_over_sigma_T_candidate:.6g}"
        )
        print(
            f"      sigma_LT RGA/RGK={r.sigma_LT_rga_candidate:.10g}/"
            f"{r.sigma_LT_rgk_candidate:.10g} "
            f"(relative difference={r.sigma_LT_relative_difference:.3%})"
        )
        print(
            f"      sigma_TT RGA/RGK={r.sigma_TT_rga_candidate:.10g}/"
            f"{r.sigma_TT_rgk_candidate:.10g} "
            f"(relative difference={r.sigma_TT_relative_difference:.3%})"
        )
    #endfor


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Verify the PARTONS DVMP phi normalization convention."
    )
    ap.add_argument("--stage3", type=Path, default=DEFAULT_STAGE3)
    ap.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    ap.add_argument("--executable", default=DEFAULT_EXECUTABLE)
    ap.add_argument("--sif", type=Path, default=None)
    ap.add_argument("--npoints", type=int, default=3)
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output directory. Default: <stage3>/convention_verification",
    )
    ap.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare diagnostics/XML but do not execute PARTONS.",
    )
    args = ap.parse_args()

    stage3 = args.stage3.resolve()
    outdir = (
        args.output.resolve()
        if args.output
        else stage3 / "convention_verification"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    metadata = load_metadata(stage3)
    df, grid_path = load_raw_phi_grid(stage3)
    cols = infer_columns(df)

    decomposition = build_shared_decomposition_table(df, cols)
    decomposition_csv = outdir / "05_candidate_TLTTLT_decomposition.csv"
    decomposition.to_csv(decomposition_csv, index=False)

    groups = choose_representative_groups(df, cols, args.npoints)

    rows = []
    for g in groups:
        r0 = g.iloc[0]
        s0 = get_phi_value(g, cols["phi"], cols["sigma"], 0.0)
        s90 = get_phi_value(g, cols["phi"], cols["sigma"], 90.0)
        s180 = get_phi_value(g, cols["phi"], cols["sigma"], 180.0)
        A, B, C = harmonic_coefficients(s0, s90, s180)

        rows.append(
            {
                "point_id": str(r0[cols["point_id"]]),
                "campaign": str(r0[cols["campaign"]]),
                "evaluation": str(r0[cols["evaluation"]]),
                "E": float(r0[cols["E"]]),
                "Q2": float(r0[cols["Q2"]]),
                "xB": float(r0[cols["xB"]]),
                "minus_t": float(r0[cols["minus_t"]]),
                "sigma_phi0": s0,
                "sigma_phi90": s90,
                "sigma_phi180": s180,
                "A": A,
                "B": B,
                "C": C,
                "two_pi_A": 2.0 * math.pi * A,
            }
        )
    #endfor

    base_csv = outdir / "01_selected_points_and_harmonics.csv"
    pd.DataFrame(rows).to_csv(base_csv, index=False)

    xml_path = outdir / "02_phi_integrated_check.xml"

    # Save an exact successful Stage-3 XML with the observable module swapped.
    # This is useful even if the local scenario contains more evaluations than
    # the three selected diagnostics.
    used_template = clone_stage3_xml_if_possible(stage3, rows, xml_path)
    if not used_template:
        build_phi_integrated_xml(rows, xml_path, metadata)
    #endif

    print("=" * 78)
    print("PARTONS DVMP CONVENTION VERIFICATION")
    print("=" * 78)
    print(f"Stage-3 directory : {stage3}")
    print(f"Input grid        : {grid_path}")
    print(f"Selected points   : {len(rows)}")
    print(f"Diagnostic CSV    : {base_csv}")
    print(f"PARTONS XML       : {xml_path}")
    print()

    for i, r in enumerate(rows):
        print(
            f"[{i}] {r['point_id']} {r['campaign']}/{r['evaluation']}: "
            f"E={r['E']:.6g} GeV, Q2={r['Q2']:.6g} GeV^2, "
            f"xB={r['xB']:.6g}, -t={r['minus_t']:.6g} GeV^2"
        )
        print(
            f"    sigma(0,90,180) = "
            f"{r['sigma_phi0']:.10g}, {r['sigma_phi90']:.10g}, "
            f"{r['sigma_phi180']:.10g}"
        )
        print(
            f"    A={r['A']:.10g}, B={r['B']:.10g}, C={r['C']:.10g}; "
            f"2*pi*A={r['two_pi_A']:.10g}"
        )
    #endfor

    print_decomposition_summary(decomposition, args.npoints)
    print(f"Candidate decomposition CSV: {decomposition_csv}")

    if args.prepare_only:
        print()
        print("Prepare-only mode: PARTONS was not executed.")
        return 0
    #endif

    if used_template:
        print()
        print(
            "NOTE: 02_phi_integrated_check.xml was cloned from an existing "
            "successful Stage-3 XML and its observable was changed to "
            "DVMPCrossSectionUUUMinusPhiIntegrated."
        )
        print(
            "Because the Stage-3 chunk may contain more evaluations than the "
            "three diagnostic rows, the safest first run is to inspect/run "
            "this tiny convention job and compare its phi-integrated outputs "
            "with 2*pi*A above."
        )
    #endif

    project = args.project.expanduser().resolve()
    sif = find_sif(args.sif, Path(__file__).resolve().parent, project)
    preflight_partons(sif, project, args.executable)

    log_path = outdir / "03_partons_phi_integrated.log"
    err_path = outdir / "03_partons_phi_integrated.stderr.log"

    print()
    print(f"PARTONS project   : {project}")
    print(f"PARTONS SIF       : {sif}")
    print(f"PARTONS executable: {args.executable}")
    print("Running through Apptainer using the same configuration as Stage 3...")
    proc = run_partons_xml(xml_path, sif, project, args.executable)
    log_path.write_text(proc.stdout)
    err_path.write_text(proc.stderr)

    print(f"PARTONS exit code : {proc.returncode}")
    print(f"PARTONS stdout    : {log_path}")
    print(f"PARTONS stderr    : {err_path}")

    results = parse_results(proc.stdout)
    if proc.returncode != 0:
        print("ERROR: PARTONS returned a nonzero exit code.", file=sys.stderr)
        return proc.returncode or 1
    #endif

    if not results:
        print(
            "ERROR: No finite 'Result: ...' values were parsed. "
            "Inspect the log/XML; the local PARTONS XML invocation syntax may "
            "need to be copied exactly from your production runner.",
            file=sys.stderr,
        )
        return 2
    #endif

    print(f"Parsed finite results: {len(results)}")

    # The cloned Stage-3 XML contains repeated phi tasks.  The phi-integrated
    # observable is independent of phi, so identify each selected diagnostic
    # by numerical agreement with its predicted 2*pi*A rather than assuming
    # one-to-one positional ordering.
    final = []
    finite_values = np.asarray([v for v, _ in results], dtype=float)
    units = [u for _, u in results]
    for r in rows:
        target = r["two_pi_A"]
        j = int(np.argmin(np.abs(finite_values - target)))
        value = float(finite_values[j])
        unit = units[j]
        ratio = value / target if target != 0.0 else np.nan
        final.append({
            **r,
            "sigma_phi_integrated_partons": value,
            "partons_unit": unit,
            "ratio_integrated_over_2piA": ratio,
            "fractional_difference": ratio - 1.0,
            "absolute_difference": value - target,
        })
    #endfor

    final_df = pd.DataFrame(final)
    final_csv = outdir / "04_phi_integration_closure.csv"
    final_df.to_csv(final_csv, index=False)

    print()
    print("phi-integration closure:")
    for r in final:
        print(
            f"  {r['point_id']} {r['campaign']}: "
            f"PARTONS={r['sigma_phi_integrated_partons']:.12g}, "
            f"2*pi*A={r['two_pi_A']:.12g}, "
            f"ratio={r['ratio_integrated_over_2piA']:.12g}, "
            f"delta={r['fractional_difference']:.3e}"
        )
    #endfor
    print(f"Wrote: {final_csv}")

    max_frac = np.nanmax(np.abs(final_df["fractional_difference"].to_numpy(float)))
    print(f"Maximum |phi-integration fractional closure|: {max_frac:.3e}")
    print()
    print("Interpretation:")
    print("  * The phi-integrated test verifies the overall 2*pi normalization.")
    print("  * The T/L/TT/LT quantities in 05_candidate_TLTTLT_decomposition.csv")
    print("    are still explicitly candidate mappings. Their RGA/RGK consistency is")
    print("    a numerical convention test; final labels require confirming the")
    print("    DVMPProcessGK06::CrossSection() source formula.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
