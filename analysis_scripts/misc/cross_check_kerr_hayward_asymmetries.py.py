#!/usr/bin/env python3

import os
import re
import sys
import math
import argparse
from dataclasses import dataclass
from typing import Dict, Tuple, List

import ROOT


# -----------------------------------------------------------------------------
# User-editable defaults
# -----------------------------------------------------------------------------

DEFAULT_KERR_DIR = "/u/home/thayward/maggies_fits"

DEFAULT_HAYWARD_DIR = (
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
    "asymmetry_extraction/output/results"
)

DEFAULT_OUTPUT_DIR = "output/cross_check"

XB_EDGES = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
XB_CENTERS = [
    0.5 * (XB_EDGES[i] + XB_EDGES[i + 1])
    for i in range(len(XB_EDGES) - 1)
]

N_XB_BINS = len(XB_CENTERS)
N_SECTORS = 6

X_AXIS_MIN = 0.10
X_AXIS_MAX = 0.70

SSA_Y_MIN = -0.20
SSA_Y_MAX = 0.20

DSA_Y_MIN = -0.10
DSA_Y_MAX = 1.00

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

PERIODS = ["Su22", "Fa22", "Sp23"]


# -----------------------------------------------------------------------------
# Asymmetry configuration
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class AsymmetryConfig:
    key: str
    directory_name: str
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
        directory_name="BSA_ALU_sinphi",
        title="BSA A_{LU}^{sin#phi}",
        y_title="A_{LU}^{sin#phi}",
        kerr_value_col=4,
        kerr_unc_col=5,
        hayward_term="ALUsinphi",
        y_min=SSA_Y_MIN,
        y_max=SSA_Y_MAX,
    ),
    AsymmetryConfig(
        key="AUL_sinphi",
        directory_name="TSA_AUL_sinphi",
        title="TSA A_{UL}^{sin#phi}",
        y_title="A_{UL}^{sin#phi}",
        kerr_value_col=8,
        kerr_unc_col=9,
        hayward_term="AULsinphi",
        y_min=SSA_Y_MIN,
        y_max=SSA_Y_MAX,
    ),
    AsymmetryConfig(
        key="AUL_sin2phi",
        directory_name="TSA_AUL_sin2phi",
        title="TSA A_{UL}^{sin2#phi}",
        y_title="A_{UL}^{sin2#phi}",
        kerr_value_col=10,
        kerr_unc_col=11,
        hayward_term="AULsin2phi",
        y_min=SSA_Y_MIN,
        y_max=SSA_Y_MAX,
    ),
    AsymmetryConfig(
        key="ALL_const",
        directory_name="DSA_ALL_const",
        title="DSA A_{LL}^{const}",
        y_title="A_{LL}^{const}",
        kerr_value_col=12,
        kerr_unc_col=13,
        hayward_term="ALL",
        y_min=DSA_Y_MIN,
        y_max=DSA_Y_MAX,
    ),
    AsymmetryConfig(
        key="ALL_cosphi",
        directory_name="DSA_ALL_cosphi",
        title="DSA A_{LL}^{cos#phi}",
        y_title="A_{LL}^{cos#phi}",
        kerr_value_col=14,
        kerr_unc_col=15,
        hayward_term="ALLcosphi",
        y_min=DSA_Y_MIN,
        y_max=DSA_Y_MAX,
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


# Structure:
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
    missing = []

    for period in PERIODS:
        for asym in ASYMMETRIES:
            for sector in range(1, N_SECTORS + 1):
                for xb_bin in range(N_XB_BINS):
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
# Plotting helpers
# -----------------------------------------------------------------------------

def configure_root() -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetTitleFillColor(0)
    ROOT.gStyle.SetPadGridX(True)
    ROOT.gStyle.SetPadGridY(True)
    ROOT.gStyle.SetEndErrorSize(4)


def draw_sector_plot(
    sector: int,
    period: str,
    asym: AsymmetryConfig,
    xb_bin: int,
    kerr_point: FitPoint,
    hayward_point: FitPoint,
    show_legend: bool,
) -> None:
    pad = ROOT.gPad
    pad.SetLeftMargin(0.16)
    pad.SetRightMargin(0.04)
    pad.SetBottomMargin(0.16)
    pad.SetTopMargin(0.13)
    pad.SetGridx(True)
    pad.SetGridy(True)

    frame_name = "frame_{}_{}_sector{}_xb{}".format(
        period,
        asym.key,
        sector,
        xb_bin,
    )

    frame_title = "Sector {}".format(sector)

    frame = ROOT.TH1D(
        frame_name,
        frame_title,
        100,
        X_AXIS_MIN,
        X_AXIS_MAX,
    )

    frame.SetMinimum(asym.y_min)
    frame.SetMaximum(asym.y_max)
    frame.GetXaxis().SetTitle("x_{B}")
    frame.GetYaxis().SetTitle(asym.y_title)
    frame.GetXaxis().CenterTitle(True)
    frame.GetYaxis().CenterTitle(True)
    frame.GetXaxis().SetTitleSize(0.060)
    frame.GetYaxis().SetTitleSize(0.060)
    frame.GetXaxis().SetLabelSize(0.050)
    frame.GetYaxis().SetLabelSize(0.050)
    frame.GetYaxis().SetTitleOffset(1.20)
    frame.GetXaxis().SetNdivisions(505)
    frame.GetYaxis().SetNdivisions(505)
    frame.Draw("AXIS")

    zero_line = ROOT.TLine(X_AXIS_MIN, 0.0, X_AXIS_MAX, 0.0)
    zero_line.SetLineStyle(2)
    zero_line.SetLineWidth(1)
    zero_line.Draw("SAME")

    x_kerr = ROOT.std.vector("double")()
    y_kerr = ROOT.std.vector("double")()
    ex_kerr = ROOT.std.vector("double")()
    ey_kerr = ROOT.std.vector("double")()

    x_hayward = ROOT.std.vector("double")()
    y_hayward = ROOT.std.vector("double")()
    ex_hayward = ROOT.std.vector("double")()
    ey_hayward = ROOT.std.vector("double")()

    x_kerr.push_back(kerr_point.x)
    y_kerr.push_back(kerr_point.value)
    ex_kerr.push_back(0.0)
    ey_kerr.push_back(kerr_point.uncertainty)

    x_hayward.push_back(hayward_point.x)
    y_hayward.push_back(hayward_point.value)
    ex_hayward.push_back(0.0)
    ey_hayward.push_back(hayward_point.uncertainty)

    graph_kerr = ROOT.TGraphErrors(1, x_kerr.data(), y_kerr.data(), ex_kerr.data(), ey_kerr.data())
    graph_hayward = ROOT.TGraphErrors(1, x_hayward.data(), y_hayward.data(), ex_hayward.data(), ey_hayward.data())

    graph_kerr.SetName(
        "graph_kerr_{}_{}_sector{}_xb{}".format(
            period,
            asym.key,
            sector,
            xb_bin,
        )
    )
    graph_hayward.SetName(
        "graph_hayward_{}_{}_sector{}_xb{}".format(
            period,
            asym.key,
            sector,
            xb_bin,
        )
    )

    graph_kerr.SetMarkerStyle(20)
    graph_kerr.SetMarkerSize(1.15)
    graph_kerr.SetMarkerColor(ROOT.kBlack)
    graph_kerr.SetLineColor(ROOT.kBlack)
    graph_kerr.SetLineWidth(2)

    graph_hayward.SetMarkerStyle(24)
    graph_hayward.SetMarkerSize(1.25)
    graph_hayward.SetMarkerColor(ROOT.kRed + 1)
    graph_hayward.SetLineColor(ROOT.kRed + 1)
    graph_hayward.SetLineWidth(2)

    graph_kerr.Draw("PE1 SAME")
    graph_hayward.Draw("PE1 SAME")

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.050)
    latex.SetTextFont(42)

    xb_low = XB_EDGES[xb_bin]
    xb_high = XB_EDGES[xb_bin + 1]
    latex.DrawLatex(
        0.20,
        0.80,
        "{:.2f} < x_{{B}} < {:.2f}".format(xb_low, xb_high),
    )

    if show_legend:
        legend = ROOT.TLegend(0.53, 0.68, 0.94, 0.88)
        legend.SetBorderSize(1)
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.045)
        legend.AddEntry(graph_kerr, "Kerr", "pe")
        legend.AddEntry(graph_hayward, "Hayward", "pe")
        legend.Draw()
    #endif

    pad.Update()

    # Keep objects alive by attaching them to the pad.
    pad.GetListOfPrimitives().Add(frame)
    pad.GetListOfPrimitives().Add(zero_line)
    pad.GetListOfPrimitives().Add(graph_kerr)
    pad.GetListOfPrimitives().Add(graph_hayward)


def draw_canvas_title(canvas: ROOT.TCanvas, title: str) -> None:
    canvas.cd()

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextAlign(22)
    latex.SetTextFont(42)
    latex.SetTextSize(0.032)
    latex.DrawLatex(0.5, 0.975, title)

    canvas.Update()


def make_comparison_plot(
    output_dir: str,
    period: str,
    asym: AsymmetryConfig,
    xb_bin: int,
    kerr_results: ResultsDict,
    hayward_results: ResultsDict,
) -> None:
    canvas_name = "c_{}_{}_xbbin{}".format(period, asym.key, xb_bin)
    canvas_title = "{}  {}  xB bin {}".format(period, asym.title, xb_bin)

    canvas = ROOT.TCanvas(canvas_name, canvas_title, 1600, 950)
    canvas.SetTopMargin(0.04)
    canvas.Divide(3, 2, 0.001, 0.001)

    for sector in range(1, N_SECTORS + 1):
        canvas.cd(sector)

        kerr_point = kerr_results[period][asym.key][sector][xb_bin]
        hayward_point = hayward_results[period][asym.key][sector][xb_bin]

        draw_sector_plot(
            sector=sector,
            period=period,
            asym=asym,
            xb_bin=xb_bin,
            kerr_point=kerr_point,
            hayward_point=hayward_point,
            show_legend=(sector == 1),
        )
    #endfor

    draw_canvas_title(
        canvas,
        "{}   {}   {:.2f} < x_{{B}} < {:.2f}".format(
            period,
            asym.title,
            XB_EDGES[xb_bin],
            XB_EDGES[xb_bin + 1],
        ),
    )

    ensure_directory(output_dir)

    filename_base = "{}_{}_xB_{:02d}_{:.2f}_{:.2f}".format(
        period,
        asym.directory_name,
        xb_bin,
        XB_EDGES[xb_bin],
        XB_EDGES[xb_bin + 1],
    )

    filename_base = filename_base.replace(".", "p")

    png_path = os.path.join(output_dir, filename_base + ".png")
    pdf_path = os.path.join(output_dir, filename_base + ".pdf")

    canvas.SaveAs(png_path)
    canvas.SaveAs(pdf_path)

    print("[cross_check] Wrote {}".format(png_path))
    print("[cross_check] Wrote {}".format(pdf_path))


def make_all_plots(
    output_root_dir: str,
    kerr_results: ResultsDict,
    hayward_results: ResultsDict,
) -> None:
    for period in PERIODS:
        period_dir = os.path.join(output_root_dir, period)
        ensure_directory(period_dir)

        for asym in ASYMMETRIES:
            asym_dir = os.path.join(period_dir, asym.directory_name)
            ensure_directory(asym_dir)

            for xb_bin in range(N_XB_BINS):
                make_comparison_plot(
                    output_dir=asym_dir,
                    period=period,
                    asym=asym,
                    xb_bin=xb_bin,
                    kerr_results=kerr_results,
                    hayward_results=hayward_results,
                )
            #endfor
        #endfor
    #endfor


# -----------------------------------------------------------------------------
# Optional text summary
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
                        combined_unc = math.sqrt(k.uncertainty * k.uncertainty + h.uncertainty * h.uncertainty)

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
        default=DEFAULT_KERR_DIR,
        help="Directory containing Kerr fit text files. Default: {}".format(DEFAULT_KERR_DIR),
    )

    parser.add_argument(
        "--hayward-dir",
        default=DEFAULT_HAYWARD_DIR,
        help="Directory containing Hayward asymmetry result text files. Default: {}".format(DEFAULT_HAYWARD_DIR),
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Master output directory. Default: {}".format(DEFAULT_OUTPUT_DIR),
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

    kerr_results = load_kerr_results(args.kerr_dir)
    hayward_results = load_hayward_results(args.hayward_dir)

    make_all_plots(
        output_root_dir=args.output_dir,
        kerr_results=kerr_results,
        hayward_results=hayward_results,
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