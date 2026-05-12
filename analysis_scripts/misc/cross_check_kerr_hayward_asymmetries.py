#!/usr/bin/env python3

import os
import re
import sys
import math
import argparse
from array import array
from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional

import ROOT


# -----------------------------------------------------------------------------
# Hard-coded paths
# -----------------------------------------------------------------------------

KERR_DIR = "/u/home/thayward/maggies_fits"

HAYWARD_DIR = (
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "asymmetry_extraction/output/results"
)

OUTPUT_DIR = "output/cross_check"


# -----------------------------------------------------------------------------
# Global configuration
# -----------------------------------------------------------------------------

XB_EDGES = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
XB_CENTERS = [
    0.5 * (XB_EDGES[i] + XB_EDGES[i + 1])
    for i in range(len(XB_EDGES) - 1)
]

N_XB_BINS = len(XB_CENTERS)
N_SECTORS = 6

SECTOR_AXIS_MIN = 0.5
SECTOR_AXIS_MAX = 6.5

ALU_SINPHI_Y_MIN = -0.30
ALU_SINPHI_Y_MAX = 0.30

AUL_SINPHI_Y_MIN = -0.30
AUL_SINPHI_Y_MAX = 0.30

AUL_SIN2PHI_Y_MIN = -0.40
AUL_SIN2PHI_Y_MAX = 0.40

ALL_Y_MIN = -1.00
ALL_Y_MAX = 1.00

ALL_COSPHI_Y_MIN = -0.80
ALL_COSPHI_Y_MAX = 0.80

BAD_EXACT_ABS_VALUE = 1.0
BAD_UNCERTAINTY_THRESHOLD = 0.5
BAD_VALUE_TOLERANCE = 1.0e-12

PERIODS = ["Su22", "Fa22", "Sp23"]

KERR_PERIOD_FILES = {
    "Su22": "pip_fit_summer22.txt",
    "Fa22": "pip_fit_fall22.txt",
    "Sp23": "pip_fit_spring23.txt",
}

HAYWARD_PERIOD_FILE_PATTERNS = {
    "Su22": r"asymmetries_rgc_su22_inb_NH3_epi\+_2_timeStamp_.*\.txt$",
    "Fa22": r"asymmetries_rgc_fa22_inb_NH3_epi\+_2_timeStamp_.*\.txt$",
    "Sp23": r"asymmetries_rgc_sp23_inb_NH3_epi\+_2_timeStamp_.*\.txt$",
}

ROOT_OBJECT_KEEPALIVE = []


# -----------------------------------------------------------------------------
# Asymmetry configuration
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class AsymmetryConfig:
    key: str
    file_name_stem: str
    title: str
    y_title: str
    kerr_value_col: int
    kerr_unc_col: int
    hayward_term: str
    y_min: float
    y_max: float


ASYMMETRIES: List[AsymmetryConfig] = [
    AsymmetryConfig(
        key="ALU_sinphi",
        file_name_stem="BSA_ALU_sinphi",
        title="BSA A_{LU}^{sin#phi}",
        y_title="F_{LU}^{sin#phi}/F_{UU}",
        kerr_value_col=4,
        kerr_unc_col=5,
        hayward_term="ALUsinphi",
        y_min=ALU_SINPHI_Y_MIN,
        y_max=ALU_SINPHI_Y_MAX,
    ),
    AsymmetryConfig(
        key="AUL_sinphi",
        file_name_stem="TSA_AUL_sinphi",
        title="TSA A_{UL}^{sin#phi}",
        y_title="F_{UL}^{sin#phi}/F_{UU}",
        kerr_value_col=8,
        kerr_unc_col=9,
        hayward_term="AULsinphi",
        y_min=AUL_SINPHI_Y_MIN,
        y_max=AUL_SINPHI_Y_MAX,
    ),
    AsymmetryConfig(
        key="AUL_sin2phi",
        file_name_stem="TSA_AUL_sin2phi",
        title="TSA A_{UL}^{sin2#phi}",
        y_title="F_{UL}^{sin2#phi}/F_{UU}",
        kerr_value_col=10,
        kerr_unc_col=11,
        hayward_term="AULsin2phi",
        y_min=AUL_SIN2PHI_Y_MIN,
        y_max=AUL_SIN2PHI_Y_MAX,
    ),
    AsymmetryConfig(
        key="ALL_const",
        file_name_stem="DSA_ALL_const",
        title="DSA A_{LL}^{const}",
        y_title="F_{LL}/F_{UU}",
        kerr_value_col=12,
        kerr_unc_col=13,
        hayward_term="ALL",
        y_min=ALL_Y_MIN,
        y_max=ALL_Y_MAX,
    ),
    AsymmetryConfig(
        key="ALL_cosphi",
        file_name_stem="DSA_ALL_cosphi",
        title="DSA A_{LL}^{cos#phi}",
        y_title="F_{LL}^{cos#phi}/F_{UU}",
        kerr_value_col=14,
        kerr_unc_col=15,
        hayward_term="ALLcosphi",
        y_min=ALL_COSPHI_Y_MIN,
        y_max=ALL_COSPHI_Y_MAX,
    ),
]


# -----------------------------------------------------------------------------
# Data containers
# -----------------------------------------------------------------------------

@dataclass
class FitPoint:
    x: float
    value: float
    uncertainty: float


@dataclass(frozen=True)
class PlotBinConfig:
    label: str
    name_tag: str
    source_xb_bins: List[int]
    x_value: float


# results[period][asym_key][sector][xb_bin_index] = FitPoint
ResultsDict = Dict[str, Dict[str, Dict[int, Dict[int, FitPoint]]]]


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

def fatal(message: str) -> None:
    print("[cross_check] FATAL: {}".format(message), file=sys.stderr)
    sys.exit(1)


def ensure_directory(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def check_file_exists(path: str) -> None:
    if not os.path.isfile(path):
        fatal("Required input file does not exist: {}".format(path))
    #endif


def is_bad_fit_point(point: FitPoint) -> bool:
    if abs(abs(point.value) - BAD_EXACT_ABS_VALUE) < BAD_VALUE_TOLERANCE:
        return True
    #endif

    if point.uncertainty > BAD_UNCERTAINTY_THRESHOLD:
        return True
    #endif

    return False


def sanitize_bad_points_and_sync_uncertainties(
    kerr_results: ResultsDict,
    hayward_results: ResultsDict,
) -> None:
    replaced_hayward = 0
    replaced_kerr = 0
    synced_kerr_uncertainties = 0

    for period in PERIODS:
        for asym in ASYMMETRIES:
            for sector in range(1, N_SECTORS + 1):
                for xb_bin in range(N_XB_BINS):
                    kerr_point = kerr_results[period][asym.key][sector][xb_bin]
                    hayward_point = hayward_results[period][asym.key][sector][xb_bin]

                    original_kerr_value = kerr_point.value
                    original_hayward_value = hayward_point.value
                    original_kerr_uncertainty = kerr_point.uncertainty
                    original_hayward_uncertainty = hayward_point.uncertainty

                    kerr_bad = is_bad_fit_point(kerr_point)
                    hayward_bad = is_bad_fit_point(hayward_point)

                    if hayward_bad and not kerr_bad:
                        hayward_point.value = original_kerr_value
                        hayward_point.uncertainty = original_kerr_uncertainty
                        replaced_hayward += 1
                    elif kerr_bad and not hayward_bad:
                        kerr_point.value = original_hayward_value
                        kerr_point.uncertainty = original_hayward_uncertainty
                        replaced_kerr += 1
                    elif hayward_bad and kerr_bad:
                        print(
                            "[cross_check] WARNING: both points flagged bad; "
                            "setting both to 0 with uncertainty 0 for "
                            "period={} asymmetry={} sector={} xb_bin={}".format(
                                period,
                                asym.key,
                                sector,
                                xb_bin,
                            ),
                            file=sys.stderr,
                        )
                        hayward_point.value = 0.0
                        hayward_point.uncertainty = 0.0
                        kerr_point.value = 0.0
                        kerr_point.uncertainty = 0.0
                        replaced_hayward += 1
                        replaced_kerr += 1
                    #endif

                    if abs(kerr_point.uncertainty - hayward_point.uncertainty) > 0.0:
                        synced_kerr_uncertainties += 1
                    #endif

                    kerr_point.uncertainty = hayward_point.uncertainty
                #endfor
            #endfor
        #endfor
    #endfor

    print(
        "[cross_check] Sanitized bad points: replaced_hayward={} replaced_kerr={} "
        "synced_kerr_uncertainties={}".format(
            replaced_hayward,
            replaced_kerr,
            synced_kerr_uncertainties,
        )
    )


def weighted_average_points(
    points: List[FitPoint],
    period_label: str,
    asym_key: str,
    sector: int,
    xb_bin: int,
    analyzer_name: str,
) -> FitPoint:
    weighted_sum = 0.0
    weight_sum = 0.0

    for point in points:
        if point.uncertainty <= 0.0:
            print(
                "[cross_check] WARNING: skipping zero/negative uncertainty point in RGC average "
                "for analyzer={} period_group={} asymmetry={} sector={} xb_bin={} "
                "value={} uncertainty={}".format(
                    analyzer_name,
                    period_label,
                    asym_key,
                    sector,
                    xb_bin,
                    point.value,
                    point.uncertainty,
                ),
                file=sys.stderr,
            )
            continue
        #endif

        weight = 1.0 / (point.uncertainty * point.uncertainty)
        weighted_sum += weight * point.value
        weight_sum += weight
    #endfor

    if weight_sum <= 0.0:
        fatal(
            "Cannot form RGC weighted average for analyzer={} asymmetry={} sector={} xb_bin={}; "
            "all candidate points had non-positive uncertainties.".format(
                analyzer_name,
                asym_key,
                sector,
                xb_bin,
            )
        )
    #endif

    combined_value = weighted_sum / weight_sum
    combined_uncertainty = math.sqrt(1.0 / weight_sum)

    return FitPoint(
        x=XB_CENTERS[xb_bin],
        value=combined_value,
        uncertainty=combined_uncertainty,
    )


def weighted_average_points_for_plot_bin(
    points: List[FitPoint],
    combined_x: float,
    period_label: str,
    asym_key: str,
    sector: int,
    plot_bin_label: str,
    analyzer_name: str,
) -> FitPoint:
    weighted_sum = 0.0
    weight_sum = 0.0

    for point in points:
        if point.uncertainty <= 0.0:
            print(
                "[cross_check] WARNING: skipping zero/negative uncertainty point in plot-bin average "
                "for analyzer={} period={} asymmetry={} sector={} plot_bin={} "
                "value={} uncertainty={}".format(
                    analyzer_name,
                    period_label,
                    asym_key,
                    sector,
                    plot_bin_label,
                    point.value,
                    point.uncertainty,
                ),
                file=sys.stderr,
            )
            continue
        #endif

        weight = 1.0 / (point.uncertainty * point.uncertainty)
        weighted_sum += weight * point.value
        weight_sum += weight
    #endfor

    if weight_sum <= 0.0:
        fatal(
            "Cannot form weighted plot-bin average for analyzer={} period={} asymmetry={} "
            "sector={} plot_bin={}; all candidate points had non-positive uncertainties.".format(
                analyzer_name,
                period_label,
                asym_key,
                sector,
                plot_bin_label,
            )
        )
    #endif

    combined_value = weighted_sum / weight_sum
    combined_uncertainty = math.sqrt(1.0 / weight_sum)

    return FitPoint(
        x=combined_x,
        value=combined_value,
        uncertainty=combined_uncertainty,
    )


def get_plot_bin_configs(combine_high_xb: bool) -> List[PlotBinConfig]:
    if combine_high_xb:
        return [
            PlotBinConfig(
                label="{:.2f} < x_{{B}} < {:.2f}".format(XB_EDGES[0], XB_EDGES[1]),
                name_tag="xB_0p10_0p20",
                source_xb_bins=[0],
                x_value=0.5 * (XB_EDGES[0] + XB_EDGES[1]),
            ),
            PlotBinConfig(
                label="{:.2f} < x_{{B}} < {:.2f}".format(XB_EDGES[1], XB_EDGES[2]),
                name_tag="xB_0p20_0p30",
                source_xb_bins=[1],
                x_value=0.5 * (XB_EDGES[1] + XB_EDGES[2]),
            ),
            PlotBinConfig(
                label="{:.2f} < x_{{B}} < {:.2f}".format(XB_EDGES[2], XB_EDGES[3]),
                name_tag="xB_0p30_0p40",
                source_xb_bins=[2],
                x_value=0.5 * (XB_EDGES[2] + XB_EDGES[3]),
            ),
            PlotBinConfig(
                label="{:.2f} < x_{{B}} < {:.2f}".format(XB_EDGES[3], XB_EDGES[5]),
                name_tag="xB_0p40_0p60",
                source_xb_bins=[3, 4],
                x_value=0.5 * (XB_EDGES[3] + XB_EDGES[5]),
            ),
        ]
    #endif

    plot_bins = []

    for xb_bin in range(N_XB_BINS):
        plot_bins.append(
            PlotBinConfig(
                label="{:.2f} < x_{{B}} < {:.2f}".format(
                    XB_EDGES[xb_bin],
                    XB_EDGES[xb_bin + 1],
                ),
                name_tag="xB_{:.2f}_{:.2f}".format(
                    XB_EDGES[xb_bin],
                    XB_EDGES[xb_bin + 1],
                ).replace(".", "p"),
                source_xb_bins=[xb_bin],
                x_value=XB_CENTERS[xb_bin],
            )
        )
    #endfor

    return plot_bins


def get_plot_point(
    results: ResultsDict,
    period: str,
    asym: AsymmetryConfig,
    sector: int,
    plot_bin: PlotBinConfig,
    analyzer_name: str,
) -> FitPoint:
    if len(plot_bin.source_xb_bins) == 1:
        xb_bin = plot_bin.source_xb_bins[0]
        return results[period][asym.key][sector][xb_bin]
    #endif

    source_points = []

    for xb_bin in plot_bin.source_xb_bins:
        source_points.append(results[period][asym.key][sector][xb_bin])
    #endfor

    return weighted_average_points_for_plot_bin(
        points=source_points,
        combined_x=plot_bin.x_value,
        period_label=period,
        asym_key=asym.key,
        sector=sector,
        plot_bin_label=plot_bin.label,
        analyzer_name=analyzer_name,
    )


def build_rgc_combined_results(
    input_results: ResultsDict,
    analyzer_name: str,
) -> ResultsDict:
    combined_results: ResultsDict = {
        "RGC": {}
    }

    for asym in ASYMMETRIES:
        combined_results["RGC"][asym.key] = {}

        for sector in range(1, N_SECTORS + 1):
            combined_results["RGC"][asym.key][sector] = {}

            for xb_bin in range(N_XB_BINS):
                period_points = []

                for period in PERIODS:
                    if period not in input_results:
                        fatal(
                            "Cannot build RGC average for analyzer={}; missing period={}.".format(
                                analyzer_name,
                                period,
                            )
                        )
                    #endif

                    if asym.key not in input_results[period]:
                        fatal(
                            "Cannot build RGC average for analyzer={}; missing asymmetry={} in period={}.".format(
                                analyzer_name,
                                asym.key,
                                period,
                            )
                        )
                    #endif

                    if sector not in input_results[period][asym.key]:
                        fatal(
                            "Cannot build RGC average for analyzer={}; missing sector={} for period={} asymmetry={}.".format(
                                analyzer_name,
                                sector,
                                period,
                                asym.key,
                            )
                        )
                    #endif

                    if xb_bin not in input_results[period][asym.key][sector]:
                        fatal(
                            "Cannot build RGC average for analyzer={}; missing xb_bin={} for period={} asymmetry={} sector={}.".format(
                                analyzer_name,
                                xb_bin,
                                period,
                                asym.key,
                                sector,
                            )
                        )
                    #endif

                    period_points.append(input_results[period][asym.key][sector][xb_bin])
                #endfor

                combined_results["RGC"][asym.key][sector][xb_bin] = weighted_average_points(
                    points=period_points,
                    period_label="RGC",
                    asym_key=asym.key,
                    sector=sector,
                    xb_bin=xb_bin,
                    analyzer_name=analyzer_name,
                )
            #endfor
        #endfor
    #endfor

    validate_complete_results_for_periods(
        results=combined_results,
        analyzer_name=analyzer_name,
        periods_to_check=["RGC"],
    )

    return combined_results


def find_unique_file_by_regex(directory: str, pattern: str, period: str) -> str:
    if not os.path.isdir(directory):
        fatal("Hayward input directory does not exist: {}".format(directory))
    #endif

    regex = re.compile(pattern)
    matches = []

    for name in os.listdir(directory):
        if regex.match(name):
            matches.append(os.path.join(directory, name))
        #endif
    #endfor

    matches.sort()

    if len(matches) == 0:
        fatal(
            "Could not find Hayward result file for period {} in directory {} "
            "using regex {}".format(period, directory, pattern)
        )
    #endif

    if len(matches) > 1:
        print(
            "[cross_check] WARNING: multiple Hayward files matched period {}. "
            "Using the lexicographically last one:".format(period),
            file=sys.stderr,
        )

        for path in matches:
            print("[cross_check]   match: {}".format(path), file=sys.stderr)
        #endfor
    #endif

    return matches[-1]


def initialize_results() -> ResultsDict:
    results: ResultsDict = {}

    for period in PERIODS:
        results[period] = {}

        for asym in ASYMMETRIES:
            results[period][asym.key] = {}

            for sector in range(1, N_SECTORS + 1):
                results[period][asym.key][sector] = {}
            #endfor
        #endfor
    #endfor

    return results


def validate_complete_results(results: ResultsDict, analyzer_name: str) -> None:
    validate_complete_results_for_periods(
        results=results,
        analyzer_name=analyzer_name,
        periods_to_check=PERIODS,
    )


def validate_complete_results_for_periods(
    results: ResultsDict,
    analyzer_name: str,
    periods_to_check: List[str],
) -> None:
    missing = []

    for period in periods_to_check:
        for asym in ASYMMETRIES:
            for sector in range(1, N_SECTORS + 1):
                for xb_bin in range(N_XB_BINS):
                    if period not in results:
                        missing.append(
                            "{} period={} missing entirely".format(
                                analyzer_name,
                                period,
                            )
                        )
                        continue
                    #endif

                    if asym.key not in results[period]:
                        missing.append(
                            "{} period={} asymmetry={} missing entirely".format(
                                analyzer_name,
                                period,
                                asym.key,
                            )
                        )
                        continue
                    #endif

                    if sector not in results[period][asym.key]:
                        missing.append(
                            "{} period={} asymmetry={} sector={} missing entirely".format(
                                analyzer_name,
                                period,
                                asym.key,
                                sector,
                            )
                        )
                        continue
                    #endif

                    if xb_bin not in results[period][asym.key][sector]:
                        missing.append(
                            "{} period={} asymmetry={} sector={} xb_bin={}".format(
                                analyzer_name,
                                period,
                                asym.key,
                                sector,
                                xb_bin,
                            )
                        )
                    #endif
                #endfor
            #endfor
        #endfor
    #endfor

    if missing:
        max_print = 100
        print(
            "[cross_check] FATAL: missing {} required fit points for {}.".format(
                len(missing),
                analyzer_name,
            ),
            file=sys.stderr,
        )

        for item in missing[:max_print]:
            print("[cross_check]   missing {}".format(item), file=sys.stderr)
        #endfor

        if len(missing) > max_print:
            print(
                "[cross_check]   ... {} more missing entries".format(
                    len(missing) - max_print
                ),
                file=sys.stderr,
            )
        #endif

        sys.exit(1)
    #endif


# -----------------------------------------------------------------------------
# Kerr parser
# -----------------------------------------------------------------------------

def parse_kerr_file(path: str, period: str, results: ResultsDict) -> None:
    check_file_exists(path)

    with open(path, "r") as infile:
        for line_number, line in enumerate(infile, start=1):
            stripped = line.strip()

            if not stripped:
                continue
            #endif

            if stripped.startswith("#"):
                continue
            #endif

            pieces = stripped.split()

            if len(pieces) < 21:
                fatal(
                    "Kerr file {} line {} has {} columns, expected at least 21. "
                    "Line: {}".format(path, line_number, len(pieces), stripped)
                )
            #endif

            try:
                xb_bin = int(pieces[0])
                kerr_sector_zero_indexed = int(pieces[1])
            except ValueError:
                fatal(
                    "Could not parse Kerr xb_bin/sector on file {} line {}: {}".format(
                        path,
                        line_number,
                        stripped,
                    )
                )
            #endtry

            sector = kerr_sector_zero_indexed + 1

            if xb_bin < 0 or xb_bin >= N_XB_BINS:
                fatal(
                    "Kerr file {} line {} has xb_bin={}, but valid range is 0 to {}.".format(
                        path,
                        line_number,
                        xb_bin,
                        N_XB_BINS - 1,
                    )
                )
            #endif

            if sector < 1 or sector > N_SECTORS:
                fatal(
                    "Kerr file {} line {} has sector_num={}, mapping to sector={}, "
                    "but valid CLAS12 sectors are 1 to 6.".format(
                        path,
                        line_number,
                        kerr_sector_zero_indexed,
                        sector,
                    )
                )
            #endif

            for asym in ASYMMETRIES:
                try:
                    value = float(pieces[asym.kerr_value_col])
                    uncertainty = float(pieces[asym.kerr_unc_col])
                except ValueError:
                    fatal(
                        "Could not parse Kerr value/uncertainty for asymmetry {} "
                        "from file {} line {}.".format(
                            asym.key,
                            path,
                            line_number,
                        )
                    )
                #endtry

                results[period][asym.key][sector][xb_bin] = FitPoint(
                    x=XB_CENTERS[xb_bin],
                    value=value,
                    uncertainty=uncertainty,
                )
            #endfor
        #endfor
    #endwith


def load_kerr_results(kerr_dir: str) -> ResultsDict:
    results = initialize_results()

    for period in PERIODS:
        filename = KERR_PERIOD_FILES[period]
        path = os.path.join(kerr_dir, filename)
        print("[cross_check] Reading Kerr {} from {}".format(period, path))
        parse_kerr_file(path, period, results)
    #endfor

    validate_complete_results(results, "Kerr")

    return results


# -----------------------------------------------------------------------------
# Hayward parser
# -----------------------------------------------------------------------------

def parse_hayward_tuple_list(tuple_text: str) -> List[Tuple[float, float, float]]:
    triples = []

    pattern = re.compile(
        r"\{\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\}"
    )

    for match in pattern.finditer(tuple_text):
        x_value = float(match.group(1))
        asym_value = float(match.group(2))
        uncertainty = float(match.group(3))
        triples.append((x_value, asym_value, uncertainty))
    #endfor

    return triples


def parse_hayward_file(path: str, period: str, results: ResultsDict) -> None:
    check_file_exists(path)

    with open(path, "r") as infile:
        content = infile.read()
    #endwith

    asym_by_hayward_term = {
        asym.hayward_term: asym
        for asym in ASYMMETRIES
    }

    assignment_pattern = re.compile(
        r"enpi_sector([1-6])GEchi2Fits([A-Za-z0-9_]+)\s*=\s*(\{\{.*?\}\})\s*;",
        re.DOTALL,
    )

    found_terms = set()

    for match in assignment_pattern.finditer(content):
        sector = int(match.group(1))
        hayward_term = match.group(2)
        tuple_text = match.group(3)

        if hayward_term not in asym_by_hayward_term:
            continue
        #endif

        asym = asym_by_hayward_term[hayward_term]
        triples = parse_hayward_tuple_list(tuple_text)

        if len(triples) != N_XB_BINS:
            fatal(
                "Hayward file {} period={} sector={} term={} has {} xB points, "
                "expected {}.".format(
                    path,
                    period,
                    sector,
                    hayward_term,
                    len(triples),
                    N_XB_BINS,
                )
            )
        #endif

        for xb_bin, triple in enumerate(triples):
            x_value, value, uncertainty = triple

            results[period][asym.key][sector][xb_bin] = FitPoint(
                x=x_value,
                value=value,
                uncertainty=uncertainty,
            )
        #endfor

        found_terms.add((sector, hayward_term))
    #endfor

    missing_terms = []

    for sector in range(1, N_SECTORS + 1):
        for asym in ASYMMETRIES:
            if (sector, asym.hayward_term) not in found_terms:
                missing_terms.append(
                    "sector={} term={}".format(sector, asym.hayward_term)
                )
            #endif
        #endfor
    #endfor

    if missing_terms:
        print(
            "[cross_check] FATAL: Hayward file {} is missing required terms:".format(
                path
            ),
            file=sys.stderr,
        )

        for item in missing_terms:
            print("[cross_check]   {}".format(item), file=sys.stderr)
        #endfor

        sys.exit(1)
    #endif


def load_hayward_results(hayward_dir: str) -> ResultsDict:
    results = initialize_results()

    for period in PERIODS:
        path = find_unique_file_by_regex(
            hayward_dir,
            HAYWARD_PERIOD_FILE_PATTERNS[period],
            period,
        )
        print("[cross_check] Reading Hayward {} from {}".format(period, path))
        parse_hayward_file(path, period, results)
    #endfor

    validate_complete_results(results, "Hayward")

    return results


# -----------------------------------------------------------------------------
# Plotting and fitting helpers
# -----------------------------------------------------------------------------

def configure_root() -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetTitleFillColor(0)
    ROOT.gStyle.SetPadGridX(False)
    ROOT.gStyle.SetPadGridY(False)
    ROOT.gStyle.SetEndErrorSize(4)


def keep_root_object(obj) -> None:
    ROOT_OBJECT_KEEPALIVE.append(obj)


def make_fit_graph(
    name: str,
    x_values_in: List[float],
    y_values_in: List[float],
    y_errors_in: List[float],
) -> ROOT.TGraphErrors:
    x_values = array("d")
    y_values = array("d")
    ex_values = array("d")
    ey_values = array("d")

    for i in range(len(x_values_in)):
        x_values.append(x_values_in[i])
        y_values.append(y_values_in[i])
        ex_values.append(0.0)

        if y_errors_in[i] > 0.0:
            ey_values.append(y_errors_in[i])
        else:
            ey_values.append(1.0e-9)
        #endif
    #endfor

    graph = ROOT.TGraphErrors(
        len(x_values_in),
        x_values,
        y_values,
        ex_values,
        ey_values,
    )

    graph.SetName(name)
    graph.SetMarkerSize(0.0)
    graph.SetLineWidth(0)
    graph._arrays = (x_values, y_values, ex_values, ey_values)

    return graph


def fit_constant_sector_dependence(
    graph: ROOT.TGraphErrors,
    fit_name: str,
    line_color: int,
    x_min: float,
    x_max: float,
) -> Tuple[ROOT.TF1, float, int, float, float]:
    fit_func = ROOT.TF1(fit_name, "pol0", x_min, x_max)
    fit_func.SetLineColor(line_color)
    fit_func.SetLineStyle(2)
    fit_func.SetLineWidth(1)

    graph.Fit(fit_func, "SQN0")

    chi2 = fit_func.GetChisquare()
    ndf = fit_func.GetNDF()

    if ndf > 0:
        chi2_ndf = chi2 / float(ndf)
        p_value = ROOT.TMath.Prob(chi2, ndf)
    else:
        chi2_ndf = 0.0
        p_value = 1.0
    #endif

    return fit_func, chi2, ndf, chi2_ndf, p_value


def format_fit_quality(chi2_ndf: float, p_value: float) -> str:
    return "c2/ndf={:.2f}, p={:.2f}".format(chi2_ndf, p_value)


def make_single_point_error_graph(
    name: str,
    x_value: float,
    y_value: float,
    y_error: float,
    line_color: int,
) -> ROOT.TGraphErrors:
    x_values = array("d", [x_value])
    y_values = array("d", [y_value])
    ex_values = array("d", [0.0])
    ey_values = array("d", [y_error])

    graph = ROOT.TGraphErrors(
        1,
        x_values,
        y_values,
        ex_values,
        ey_values,
    )

    graph.SetName(name)
    graph.SetLineColor(line_color)
    graph.SetLineWidth(2)
    graph.SetMarkerColor(line_color)
    graph.SetMarkerSize(0.0)
    graph._arrays = (x_values, y_values, ex_values, ey_values)

    return graph


def make_marker(
    x_value: float,
    y_value: float,
    marker_style: int,
    marker_color: int,
    marker_size: float,
) -> ROOT.TMarker:
    marker = ROOT.TMarker(x_value, y_value, marker_style)
    marker.SetMarkerColor(marker_color)
    marker.SetMarkerSize(marker_size)
    return marker


def draw_empty_sixth_pad() -> None:
    pad = ROOT.gPad
    pad.SetFillColor(0)
    pad.SetBorderMode(0)
    pad.SetFrameBorderMode(0)
    pad.Clear()


def draw_xb_bin_title(plot_bin: PlotBinConfig) -> ROOT.TLatex:
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(42)
    latex.SetTextSize(0.055)
    latex.SetTextAlign(22)
    latex.DrawLatex(
        0.50,
        0.91,
        plot_bin.label,
    )
    return latex


def draw_xb_bin_subplot(
    plot_bin_index: int,
    plot_bin: PlotBinConfig,
    period: str,
    asym: AsymmetryConfig,
    kerr_results: ResultsDict,
    hayward_results: Optional[ResultsDict],
    plot_hayward: bool,
) -> None:
    pad = ROOT.gPad
    pad.Clear()
    pad.SetLeftMargin(0.16)
    pad.SetRightMargin(0.04)
    pad.SetBottomMargin(0.16)
    pad.SetTopMargin(0.16)
    pad.SetGridx(False)
    pad.SetGridy(False)

    frame_name = "frame_{}_{}_{}".format(
        period,
        asym.key,
        plot_bin.name_tag,
    )

    frame = ROOT.TH1D(
        frame_name,
        "",
        6,
        SECTOR_AXIS_MIN,
        SECTOR_AXIS_MAX,
    )

    frame.SetDirectory(0)
    frame.SetMinimum(asym.y_min)
    frame.SetMaximum(asym.y_max)
    frame.GetXaxis().SetTitle("Electron sector")
    frame.GetYaxis().SetTitle(asym.y_title)
    frame.GetXaxis().CenterTitle(True)
    frame.GetYaxis().CenterTitle(True)
    frame.GetXaxis().SetTitleSize(0.060)
    frame.GetYaxis().SetTitleSize(0.060)
    frame.GetXaxis().SetLabelSize(0.050)
    frame.GetYaxis().SetLabelSize(0.050)
    frame.GetYaxis().SetTitleOffset(1.18)
    frame.GetXaxis().SetNdivisions(6, False)
    frame.GetYaxis().SetNdivisions(505)

    for sector in range(1, N_SECTORS + 1):
        frame.GetXaxis().SetBinLabel(sector, str(sector))
    #endfor

    frame.Draw("AXIS")
    keep_root_object(frame)

    subplot_title = draw_xb_bin_title(plot_bin)
    keep_root_object(subplot_title)

    zero_line = ROOT.TLine(SECTOR_AXIS_MIN, 0.0, SECTOR_AXIS_MAX, 0.0)
    zero_line.SetLineStyle(2)
    zero_line.SetLineWidth(1)
    zero_line.SetLineColor(ROOT.kGray + 1)
    zero_line.Draw("SAME")
    keep_root_object(zero_line)

    hayward_x_values = []
    hayward_y_values = []
    hayward_y_errors = []

    kerr_x_values = []
    kerr_y_values = []
    kerr_y_errors = []

    first_hayward_marker = None
    first_kerr_marker = None

    for sector in range(1, N_SECTORS + 1):
        kerr_point = get_plot_point(
            results=kerr_results,
            period=period,
            asym=asym,
            sector=sector,
            plot_bin=plot_bin,
            analyzer_name="Kerr",
        )

        if plot_hayward:
            if hayward_results is None:
                fatal("Internal error: plot_hayward is true but hayward_results is None.")
            #endif

            hayward_point = get_plot_point(
                results=hayward_results,
                period=period,
                asym=asym,
                sector=sector,
                plot_bin=plot_bin,
                analyzer_name="Hayward",
            )

            x_kerr = float(sector) - 0.08
            x_hayward = float(sector) + 0.08
        else:
            hayward_point = None
            x_kerr = float(sector)
            x_hayward = float(sector)
        #endif

        kerr_x_values.append(x_kerr)
        kerr_y_values.append(kerr_point.value)
        kerr_y_errors.append(kerr_point.uncertainty)

        kerr_error_graph = make_single_point_error_graph(
            name="err_kerr_{}_{}_{}_sector{}".format(
                period,
                asym.key,
                plot_bin.name_tag,
                sector,
            ),
            x_value=x_kerr,
            y_value=kerr_point.value,
            y_error=kerr_point.uncertainty,
            line_color=ROOT.kBlack,
        )

        kerr_marker = make_marker(
            x_value=x_kerr,
            y_value=kerr_point.value,
            marker_style=20,
            marker_color=ROOT.kBlack,
            marker_size=1.45,
        )

        if plot_hayward:
            hayward_x_values.append(x_hayward)
            hayward_y_values.append(hayward_point.value)
            hayward_y_errors.append(hayward_point.uncertainty)

            hayward_error_graph = make_single_point_error_graph(
                name="err_hayward_{}_{}_{}_sector{}".format(
                    period,
                    asym.key,
                    plot_bin.name_tag,
                    sector,
                ),
                x_value=x_hayward,
                y_value=hayward_point.value,
                y_error=hayward_point.uncertainty,
                line_color=ROOT.kRed + 1,
            )

            hayward_marker = make_marker(
                x_value=x_hayward,
                y_value=hayward_point.value,
                marker_style=21,
                marker_color=ROOT.kRed + 1,
                marker_size=1.45,
            )

            hayward_error_graph.Draw("E1 SAME")
            hayward_marker.Draw("SAME")

            keep_root_object(hayward_error_graph)
            keep_root_object(hayward_marker)

            if first_hayward_marker is None:
                first_hayward_marker = hayward_marker
            #endif
        #endif

        kerr_error_graph.Draw("E1 SAME")
        kerr_marker.Draw("SAME")

        keep_root_object(kerr_error_graph)
        keep_root_object(kerr_marker)

        if first_kerr_marker is None:
            first_kerr_marker = kerr_marker
        #endif
    #endfor

    kerr_fit_graph = make_fit_graph(
        name="fit_graph_kerr_{}_{}_{}".format(
            period,
            asym.key,
            plot_bin.name_tag,
        ),
        x_values_in=kerr_x_values,
        y_values_in=kerr_y_values,
        y_errors_in=kerr_y_errors,
    )

    kerr_fit_func, kerr_chi2, kerr_ndf, kerr_chi2_ndf, kerr_p_value = fit_constant_sector_dependence(
        graph=kerr_fit_graph,
        fit_name="fit_kerr_{}_{}_{}".format(
            period,
            asym.key,
            plot_bin.name_tag,
        ),
        line_color=ROOT.kBlack,
        x_min=SECTOR_AXIS_MIN,
        x_max=SECTOR_AXIS_MAX,
    )

    kerr_fit_func.Draw("SAME")
    keep_root_object(kerr_fit_graph)
    keep_root_object(kerr_fit_func)

    if plot_hayward:
        hayward_fit_graph = make_fit_graph(
            name="fit_graph_hayward_{}_{}_{}".format(
                period,
                asym.key,
                plot_bin.name_tag,
            ),
            x_values_in=hayward_x_values,
            y_values_in=hayward_y_values,
            y_errors_in=hayward_y_errors,
        )

        hayward_fit_func, hayward_chi2, hayward_ndf, hayward_chi2_ndf, hayward_p_value = fit_constant_sector_dependence(
            graph=hayward_fit_graph,
            fit_name="fit_hayward_{}_{}_{}".format(
                period,
                asym.key,
                plot_bin.name_tag,
            ),
            line_color=ROOT.kRed + 1,
            x_min=SECTOR_AXIS_MIN,
            x_max=SECTOR_AXIS_MAX,
        )

        hayward_fit_func.Draw("SAME")
        keep_root_object(hayward_fit_graph)
        keep_root_object(hayward_fit_func)

        legend = ROOT.TLegend(0.53, 0.16, 0.94, 0.31)
        legend.SetBorderSize(1)
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.025)
        legend.AddEntry(
            first_hayward_marker,
            "Hayward, {}".format(format_fit_quality(hayward_chi2_ndf, hayward_p_value)),
            "p",
        )
        legend.AddEntry(
            first_kerr_marker,
            "Kerr, {}".format(format_fit_quality(kerr_chi2_ndf, kerr_p_value)),
            "p",
        )
        legend.Draw("SAME")
        keep_root_object(legend)
    #endif

    pad.Modified()
    pad.Update()


def draw_canvas_title(title_pad: ROOT.TPad, title: str) -> None:
    title_pad.cd()
    title_pad.Clear()

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextAlign(22)
    latex.SetTextFont(42)
    latex.SetTextSize(0.42)
    latex.DrawLatex(0.5, 0.50, title)

    keep_root_object(latex)

    title_pad.Modified()
    title_pad.Update()


def make_comparison_plot(
    output_dir: str,
    period: str,
    asym: AsymmetryConfig,
    kerr_results: ResultsDict,
    hayward_results: Optional[ResultsDict],
    plot_hayward: bool,
    combine_high_xb: bool,
) -> None:
    canvas_name = "c_{}_{}".format(period, asym.key)
    canvas_title = "{}  {}".format(period, asym.title)

    if combine_high_xb:
        canvas = ROOT.TCanvas(canvas_name, canvas_title, 1350, 1000)
    else:
        canvas = ROOT.TCanvas(canvas_name, canvas_title, 1600, 1050)
    #endif

    keep_root_object(canvas)

    title_pad = ROOT.TPad(
        "title_pad_{}_{}".format(period, asym.key),
        "title_pad_{}_{}".format(period, asym.key),
        0.0,
        0.925,
        1.0,
        1.0,
    )
    title_pad.SetFillColor(0)
    title_pad.SetBorderMode(0)
    title_pad.SetFrameBorderMode(0)
    title_pad.Draw()
    keep_root_object(title_pad)

    grid_pad = ROOT.TPad(
        "grid_pad_{}_{}".format(period, asym.key),
        "grid_pad_{}_{}".format(period, asym.key),
        0.0,
        0.0,
        1.0,
        0.925,
    )
    grid_pad.SetFillColor(0)
    grid_pad.SetBorderMode(0)
    grid_pad.SetFrameBorderMode(0)
    grid_pad.Draw()
    keep_root_object(grid_pad)

    plot_bins = get_plot_bin_configs(combine_high_xb)

    grid_pad.cd()

    if combine_high_xb:
        grid_pad.Divide(2, 2, 0.001, 0.001)
    else:
        grid_pad.Divide(3, 2, 0.001, 0.001)
    #endif

    for plot_bin_index, plot_bin in enumerate(plot_bins):
        grid_pad.cd(plot_bin_index + 1)

        draw_xb_bin_subplot(
            plot_bin_index=plot_bin_index,
            plot_bin=plot_bin,
            period=period,
            asym=asym,
            kerr_results=kerr_results,
            hayward_results=hayward_results,
            plot_hayward=plot_hayward,
        )
    #endfor

    if not combine_high_xb:
        grid_pad.cd(6)
        draw_empty_sixth_pad()
    #endif

    if plot_hayward:
        mode_suffix = "sector dependence"
    else:
        mode_suffix = "sector dependence"
    #endif

    if combine_high_xb:
        mode_suffix += ", high-x_{B} combined"
    #endif

    draw_canvas_title(
        title_pad,
        "{}   {}   {}".format(period, asym.title, mode_suffix),
    )

    canvas.cd()
    canvas.Modified()
    canvas.Update()

    ensure_directory(output_dir)

    filename_base = "{}_{}".format(period, asym.file_name_stem)

    if combine_high_xb:
        filename_base += "_xBhigh_combined"
    #endif

    if not plot_hayward:
        filename_base += "_Kerr_only"
    #endif

    png_path = os.path.join(output_dir, filename_base + ".png")

    canvas.SaveAs(png_path)

    print("[cross_check] Wrote {}".format(png_path))


def make_all_plots(
    output_root_dir: str,
    kerr_results: ResultsDict,
    hayward_results: ResultsDict,
    plot_hayward: bool,
    combine_high_xb: bool,
) -> None:
    for period in PERIODS:
        period_dir = os.path.join(output_root_dir, period)
        ensure_directory(period_dir)

        for asym in ASYMMETRIES:
            make_comparison_plot(
                output_dir=period_dir,
                period=period,
                asym=asym,
                kerr_results=kerr_results,
                hayward_results=hayward_results,
                plot_hayward=plot_hayward,
                combine_high_xb=combine_high_xb,
            )
        #endfor
    #endfor

    rgc_kerr_results = build_rgc_combined_results(
        input_results=kerr_results,
        analyzer_name="Kerr",
    )

    rgc_hayward_results = build_rgc_combined_results(
        input_results=hayward_results,
        analyzer_name="Hayward",
    )

    rgc_dir = os.path.join(output_root_dir, "RGC")
    ensure_directory(rgc_dir)

    for asym in ASYMMETRIES:
        make_comparison_plot(
            output_dir=rgc_dir,
            period="RGC",
            asym=asym,
            kerr_results=rgc_kerr_results,
            hayward_results=rgc_hayward_results,
            plot_hayward=plot_hayward,
            combine_high_xb=combine_high_xb,
        )
    #endfor


# -----------------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------------

def write_summary_tables(
    output_root_dir: str,
    kerr_results: ResultsDict,
    hayward_results: ResultsDict,
) -> None:
    summary_path = os.path.join(output_root_dir, "kerr_hayward_cross_check_summary.csv")
    ensure_directory(output_root_dir)

    with open(summary_path, "w") as outfile:
        outfile.write(
            "period,asymmetry,sector,xb_bin,xb_min,xb_max,"
            "kerr_x,kerr_value,kerr_unc,"
            "hayward_x,hayward_value,hayward_unc,"
            "difference,combined_unc,pull\n"
        )

        for period in PERIODS:
            for asym in ASYMMETRIES:
                for sector in range(1, N_SECTORS + 1):
                    for xb_bin in range(N_XB_BINS):
                        k = kerr_results[period][asym.key][sector][xb_bin]
                        h = hayward_results[period][asym.key][sector][xb_bin]

                        difference = h.value - k.value
                        combined_unc = math.sqrt(
                            k.uncertainty * k.uncertainty
                            + h.uncertainty * h.uncertainty
                        )

                        if combined_unc > 0.0:
                            pull = difference / combined_unc
                        else:
                            pull = 0.0
                        #endif

                        outfile.write(
                            "{},{},{},{},{:.8f},{:.8f},"
                            "{:.9f},{:.9f},{:.9f},"
                            "{:.9f},{:.9f},{:.9f},"
                            "{:.9f},{:.9f},{:.9f}\n".format(
                                period,
                                asym.key,
                                sector,
                                xb_bin,
                                XB_EDGES[xb_bin],
                                XB_EDGES[xb_bin + 1],
                                k.x,
                                k.value,
                                k.uncertainty,
                                h.x,
                                h.value,
                                h.uncertainty,
                                difference,
                                combined_unc,
                                pull,
                            )
                        )
                    #endfor
                #endfor
            #endfor
        #endfor
    #endwith

    print("[cross_check] Wrote {}".format(summary_path))


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Kerr and Hayward sector-dependent RGC asymmetry fit results "
            "for BSA, TSA, and DSA terms."
        )
    )

    parser.add_argument(
        "--kerr-dir",
        default=KERR_DIR,
        help="Override Kerr fit text-file directory.",
    )

    parser.add_argument(
        "--hayward-dir",
        default=HAYWARD_DIR,
        help="Override Hayward asymmetry result text-file directory.",
    )

    parser.add_argument(
        "--output-dir",
        default=OUTPUT_DIR,
        help="Override master output directory.",
    )

    parser.add_argument(
        "--kerr-only",
        action="store_true",
        help=(
            "Only plot Kerr points and Kerr constant-fit lines. "
            "Hayward files are still loaded so Kerr uncertainties can be synchronized "
            "to Hayward uncertainties before plotting."
        ),
    )

    parser.add_argument(
        "--combine-high-xb",
        action="store_true",
        help=(
            "Use 2x2 canvases instead of 2x3 canvases by combining the final two "
            "xB bins into one weighted-average bin: 0.40 < x_B < 0.60."
        ),
    )

    parser.add_argument(
        "--no-summary",
        action="store_true",
        help="Disable writing the CSV summary table.",
    )

    return parser


def main() -> int:
    args = build_arg_parser().parse_args()

    configure_root()

    print("[cross_check] Kerr directory: {}".format(args.kerr_dir))
    print("[cross_check] Hayward directory: {}".format(args.hayward_dir))
    print("[cross_check] Output directory: {}".format(args.output_dir))

    if args.kerr_only:
        print("[cross_check] Plotting mode: Kerr only")
    else:
        print("[cross_check] Plotting mode: Kerr and Hayward")
    #endif

    if args.combine_high_xb:
        print("[cross_check] xB plotting mode: combine 0.40-0.50 and 0.50-0.60 into 0.40-0.60")
    else:
        print("[cross_check] xB plotting mode: standard five xB bins")
    #endif

    kerr_results = load_kerr_results(args.kerr_dir)
    hayward_results = load_hayward_results(args.hayward_dir)

    sanitize_bad_points_and_sync_uncertainties(
        kerr_results=kerr_results,
        hayward_results=hayward_results,
    )

    make_all_plots(
        output_root_dir=args.output_dir,
        kerr_results=kerr_results,
        hayward_results=hayward_results,
        plot_hayward=(not args.kerr_only),
        combine_high_xb=args.combine_high_xb,
    )

    if not args.no_summary:
        write_summary_tables(
            output_root_dir=args.output_dir,
            kerr_results=kerr_results,
            hayward_results=hayward_results,
        )
    #endif

    print("[cross_check] Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif