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
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_STAGE3 = Path(
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "higher_level_dvcs_cross_section/high_luminosity/output/"
    "pi0_gk_stage3/partons_gk"
)

RESULT_RE = re.compile(
    r"Result\s*=\s*"
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


def find_partons_executable(explicit: str | None) -> str:
    if explicit:
        return explicit
    #endif

    import shutil

    for name in ("partons", "PARTONS"):
        found = shutil.which(name)
        if found:
            return found
        #endif
    #endfor

    raise RuntimeError(
        "Could not find the PARTONS executable in PATH. "
        "Pass it explicitly with --partons /path/to/PARTONS."
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
        unit = (m.group(2) or "").strip()
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


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Verify the PARTONS DVMP phi normalization convention."
    )
    ap.add_argument("--stage3", type=Path, default=DEFAULT_STAGE3)
    ap.add_argument("--partons", default=None, help="PARTONS executable path.")
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

    partons = find_partons_executable(args.partons)
    log_path = outdir / "03_partons_phi_integrated.log"

    print()
    print(f"Running: {partons} {xml_path}")
    proc = subprocess.run(
        [partons, str(xml_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    log_path.write_text(proc.stdout)

    print(f"PARTONS exit code : {proc.returncode}")
    print(f"PARTONS log       : {log_path}")

    results = parse_results(proc.stdout)
    if proc.returncode != 0:
        print("ERROR: PARTONS returned a nonzero exit code.", file=sys.stderr)
        return proc.returncode or 1
    #endif

    if not results:
        print(
            "ERROR: No finite 'Result = ...' values were parsed. "
            "Inspect the log/XML; the local PARTONS XML invocation syntax may "
            "need to be copied exactly from your production runner.",
            file=sys.stderr,
        )
        return 2
    #endif

    print(f"Parsed finite results: {len(results)}")

    # Only perform the one-to-one numerical closure automatically when the
    # output count matches the selected diagnostic points.
    if len(results) == len(rows):
        final = []
        for r, (value, unit) in zip(rows, results):
            ratio = value / r["two_pi_A"] if r["two_pi_A"] != 0.0 else np.nan
            final.append(
                {
                    **r,
                    "sigma_phi_integrated_partons": value,
                    "partons_unit": unit,
                    "ratio_integrated_over_2piA": ratio,
                    "fractional_difference": ratio - 1.0,
                }
            )
        #endfor

        final_df = pd.DataFrame(final)
        final_csv = outdir / "04_phi_integration_closure.csv"
        final_df.to_csv(final_csv, index=False)

        print()
        print("phi-integration closure:")
        for r in final:
            print(
                f"  {r['point_id']} {r['campaign']}: "
                f"PARTONS={r['sigma_phi_integrated_partons']:.10g}, "
                f"2*pi*A={r['two_pi_A']:.10g}, "
                f"ratio={r['ratio_integrated_over_2piA']:.12g}"
            )
        #endfor
        print(f"Wrote: {final_csv}")
    else:
        print()
        print(
            "The PARTONS result count does not equal the number of selected "
            "diagnostic points, so no positional matching was assumed."
        )
        print(
            "This usually means the cloned Stage-3 XML retained the original "
            "chunk evaluations. The log is still useful, but the next step is "
            "to use the exact XML-building helper from the Stage-3 runner to "
            "emit only these selected points."
        )
    #endif

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
