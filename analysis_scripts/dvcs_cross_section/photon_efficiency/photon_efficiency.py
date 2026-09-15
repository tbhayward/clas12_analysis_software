#!/usr/bin/env python3
"""
CLAS12 photon-efficiency analysis: Stage 1 selection development.

Purpose
-------
Build and inspect the event selections only. This script does NOT normalize MC
to data and does NOT calculate the photon efficiency.

Stage 1A -- preparation cuts
    y < 0.85
    p_e > 2.0 GeV
    Q2 > 1.0 GeV^2
    W > 2.0 GeV
    -0.30 < Mx2(ep) < 0.40 GeV^2
    -8 < vz_e < 2 cm
    |vz_e - vz_p| < 20 cm

Stage 1B -- common gamma1 cuts, applied after Stage 1A
    0.9 < beta_gamma1 < 1.1
    angle(e,gamma1) > 8 deg
    p_gamma1 >= 0.4 GeV
    tag_pass_fiducial == 1

Stage 1C -- common inferred-probe X selection
    X is the inferred photon stored in probe_raw_*.
    These cuts are common to the denominator and numerator:
    -0.231 < Mx2(ep) < 0.309 GeV^2
    Mx2(e gamma1) > 1.4 GeV^2

    The angle(gamma1,X) requirement is deliberately disabled in this
    diagnostic iteration so its effect can be isolated.

Stage 1D -- reconstructed gamma2 selection, numerator only
    Search neutral_[0..4] for reconstructed PID-22 candidates, excluding the
    tag REC index. Select the PID-22 candidate closest to inferred X.
    p_gamma2 >= 0.4 GeV

    The historical gamma2 beta cut was conditional and its run setting was not
    recorded, so it is NOT imposed here. The skim does not contain a saved
    neutral_pass_fiducial boolean; neutral PCAL coordinates are retained for
    adding the exact gamma2 PCAL fiducial definition once specified.

Diagnostics also retain the nearest neutral candidate regardless of PID so that
PID behavior is visible rather than hidden by the numerator definition.

Inputs
------
Automatically discovers completed skim_version==4 ROOT files beneath:
    /work/clas12/thayward/photon_efficiency/ROOT_trees/<sample>/<period>/

Samples:
    data, clasdis, aaogen, dvcsgen

Unavailable samples are simply skipped.

Outputs
-------
output/1_photon_efficiency_stage1_preparation_<period>.png
output/2_photon_efficiency_stage1_common_photon_<period>.png
output/3_photon_efficiency_stage1_probe_<period>.png
output/5_photon_efficiency_stage1_cutflow_<period>.txt

All plotted distributions are unit-normalized shape comparisons. MC
normalization belongs to Stage 2.
"""

import argparse
import glob
import math
import os
import shutil
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

BASE = "/work/clas12/thayward/photon_efficiency/ROOT_trees"
TREE = "PhotonEfficiency"
EXPECTED_SKIM_VERSION = 4
N_NEUTRAL = 5

SAMPLES = [
    ("data", "Data"),
    ("clasdis", "CLASDIS incl."),
    ("aaogen", "AAOgen"),
    ("dvcsgen", "DVCSgen"),
]

COLORS = {
    "data": ROOT.kBlack,
    "clasdis": ROOT.kBlue + 1,
    "aaogen": ROOT.kRed + 1,
    "dvcsgen": ROOT.kGreen + 2,
}

# Exact requested numerical boundaries.
PROBE_ANGLE_LOW = 3.0
PROBE_ANGLE_HIGH = 4.12 + 3.0 * 1.687       # 9.181 deg
PROBE_MX2_LOW = 0.039 - 3.0 * 0.09          # -0.231 GeV^2
PROBE_MX2_HIGH = 0.039 + 3.0 * 0.09         # +0.309 GeV^2


# ---------------------------------------------------------------------------
# C++ helpers callable from RDataFrame
# ---------------------------------------------------------------------------

ROOT.gInterpreter.Declare(r"""
#include <ROOT/RVec.hxx>
#include <cmath>
#include <limits>

using ROOT::VecOps::RVec;

double pe_angle_deg(double th1, double ph1, double th2, double ph2) {
    const double x1 = std::sin(th1)*std::cos(ph1);
    const double y1 = std::sin(th1)*std::sin(ph1);
    const double z1 = std::cos(th1);
    const double x2 = std::sin(th2)*std::cos(ph2);
    const double y2 = std::sin(th2)*std::sin(ph2);
    const double z2 = std::cos(th2);

    double c = x1*x2 + y1*y2 + z1*z2;
    if (c >  1.0) c =  1.0;
    if (c < -1.0) c = -1.0;
    return std::acos(c) * 180.0/M_PI;
}

double pe_mx2_egamma(double ebeam,
                     double ep, double eth, double eph,
                     double gp, double gth, double gph) {
    const double me = 0.00051099895;
    const double mp = 0.9382720813;

    const double Ee = std::sqrt(ep*ep + me*me);
    const double ex = ep*std::sin(eth)*std::cos(eph);
    const double ey = ep*std::sin(eth)*std::sin(eph);
    const double ez = ep*std::cos(eth);

    const double gx = gp*std::sin(gth)*std::cos(gph);
    const double gy = gp*std::sin(gth)*std::sin(gph);
    const double gz = gp*std::cos(gth);

    const double Em = ebeam + mp - Ee - gp;
    const double px = -ex - gx;
    const double py = -ey - gy;
    const double pz = ebeam - ez - gz;

    return Em*Em - px*px - py*py - pz*pz;
}

double pe_emiss_ep(double ebeam, double ep, double pp) {
    const double me = 0.00051099895;
    const double mp = 0.9382720813;
    const double Ee = std::sqrt(ep*ep + me*me);
    const double Ep = std::sqrt(pp*pp + mp*mp);
    return ebeam + mp - Ee - Ep;
}

double pe_emiss_epg(double ebeam, double ep, double pp, double gp) {
    return pe_emiss_ep(ebeam, ep, pp) - gp;
}

/*
 * Return the retained neutral slot closest to X.
 *
 * pid22_only=false: nearest retained neutral candidate of any PID.
 * pid22_only=true : nearest reconstructed PID-22 candidate.
 *
 * The tag REC index is explicitly excluded.
 */
int pe_best_neutral_slot(const RVec<int>& idx,
                         const RVec<int>& pid,
                         const RVec<double>& da,
                         int tag_idx,
                         bool pid22_only) {
    int best = -1;
    double best_da = std::numeric_limits<double>::infinity();

    const size_t n = std::min(idx.size(), std::min(pid.size(), da.size()));
    for (size_t i=0; i<n; ++i) {
        if (idx[i] < 0) continue;
        if (idx[i] == tag_idx) continue;
        if (pid22_only && pid[i] != 22) continue;
        if (!std::isfinite(da[i]) || da[i] < 0.0) continue;

        if (da[i] < best_da) {
            best_da = da[i];
            best = static_cast<int>(i);
        }
    }
    return best;
}

double pe_vec_value(const RVec<double>& v, int i, double missing=-999.0) {
    return (i >= 0 && static_cast<size_t>(i) < v.size()) ? v[i] : missing;
}

int pe_vec_ivalue(const RVec<int>& v, int i, int missing=-999) {
    return (i >= 0 && static_cast<size_t>(i) < v.size()) ? v[i] : missing;
}

double pe_rad_to_deg(double x) {
    return x * 180.0/M_PI;
}
""")


# ---------------------------------------------------------------------------
# Cuts
# ---------------------------------------------------------------------------

PREP_CUTS = [
    ("y < 0.85", "y < 0.85"),
    ("p_e > 2.0 GeV", "e_p > 2.0"),
    ("Q2 > 1.0 GeV^2", "Q2 > 1.0"),
    ("W > 2.0 GeV", "W > 2.0"),
    ("-0.30 < Mx2(ep) < 0.40 GeV^2",
     "Mx2_ep > -0.30 && Mx2_ep < 0.40"),
    ("-8 < vz_e < 2 cm", "e_vz > -8.0 && e_vz < 2.0"),
    ("|vz_e-vz_p| < 20 cm", "abs(e_vz-p_vz) < 20.0"),
]

COMMON_PHOTON_CUTS = [
    ("0.9 < beta_gamma1 < 1.1",
     "tag_beta > 0.9 && tag_beta < 1.1"),
    ("angle(e,gamma1) > 8 deg",
     "e_gamma1_angle_deg > 8.0"),
    ("p_gamma1 >= 0.4 GeV",
     "tag_corr_p >= 0.4"),
    ("gamma1 fiducial",
     "tag_pass_fiducial == 1"),
]

COMMON_PROBE_CUTS = [
    # angle(gamma1,X) deliberately disabled for this diagnostic iteration.
    ("-0.231 < Mx2(ep) < 0.309 GeV^2",
     f"Mx2_ep > {PROBE_MX2_LOW:.6f} && "
     f"Mx2_ep < {PROBE_MX2_HIGH:.6f}"),
    ("Mx2(e gamma1) > 1.4 GeV^2",
     "Mx2_egamma1 > 1.4"),
]

GAMMA2_CUTS = [
    ("p_gamma2 >= 0.4 GeV",
     "probe_p >= 0.4"),
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--period", default="fa18_inb")
    p.add_argument("--output-dir", default="output")
    p.add_argument(
        "--threads", type=int, default=2,
        help="ROOT implicit-multithreading worker count (default: 2)"
    )
    return p.parse_args()


def usable_file(path):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        return False, "cannot open"

    t = f.Get(TREE)
    if not t:
        f.Close()
        return False, f"missing {TREE}"

    n = int(t.GetEntries())
    if n == 0:
        f.Close()
        return False, "empty tree"

    if not t.GetBranch("skim_version"):
        f.Close()
        return False, "missing skim_version"

    t.GetEntry(0)
    version = int(getattr(t, "skim_version"))
    f.Close()

    if version != EXPECTED_SKIM_VERSION:
        return False, f"skim_version={version}, expected {EXPECTED_SKIM_VERSION}"

    return True, "ok"


def discover(sample, period):
    pattern = os.path.join(
        BASE, sample, period, "*_photon_efficiency.root"
    )
    candidates = sorted(glob.glob(pattern))
    good = []

    print(f"\n[{sample}] {pattern}")
    if not candidates:
        print("  no completed ROOT files found")
        return good

    for path in candidates:
        ok, reason = usable_file(path)
        name = os.path.basename(path)
        if ok:
            good.append(path)
            print(f"  USE   {name}")
        else:
            print(f"  SKIP  {name}: {reason}")

    print(f"  -> {len(good)} usable file(s)")
    return good


def make_chain(paths):
    ch = ROOT.TChain(TREE)
    for path in paths:
        ch.Add(path)
    return ch


def make_dataframe(chain):
    """
    Add derived quantities once. The actual cuts remain explicit below.
    """
    df = ROOT.RDataFrame(chain)

    # RDataFrame histogram actions take column names. Define all expressions
    # that are plotted as explicit columns once here.
    df = df.Define("vz_e_minus_vz_p", "e_vz - p_vz")

    df = df.Define(
        "e_gamma1_angle_deg",
        "pe_angle_deg(e_theta,e_phi,tag_corr_theta,tag_corr_phi)"
    )

    # X is the inferred photon direction from the corrected e'p' system,
    # stored by the producer in probe_raw_*.
    df = df.Define(
        "gamma1_X_angle_deg",
        "pe_angle_deg(tag_corr_theta,tag_corr_phi,probe_raw_theta,probe_raw_phi)"
    )

    df = df.Define(
        "Mx2_egamma1",
        "pe_mx2_egamma(beam_energy,e_p,e_theta,e_phi,"
        "tag_corr_p,tag_corr_theta,tag_corr_phi)"
    )

    df = df.Define(
        "Emiss_ep",
        "pe_emiss_ep(beam_energy,e_p,p_corr_p)"
    )

    df = df.Define(
        "Emiss_epg",
        "pe_emiss_epg(beam_energy,e_p,p_corr_p,tag_corr_p)"
    )

    # Diagnostic nearest neutral candidate, irrespective of PID.
    df = df.Define(
        "nearest_neutral_slot",
        "pe_best_neutral_slot(neutral_idx,neutral_pid,"
        "neutral_delta_alpha,tag_rec_index,false)"
    )
    df = df.Define(
        "nearest_neutral_exists",
        "nearest_neutral_slot >= 0"
    )
    df = df.Define(
        "nearest_neutral_pid",
        "pe_vec_ivalue(neutral_pid,nearest_neutral_slot)"
    )
    df = df.Define(
        "nearest_neutral_angle_deg",
        "pe_rad_to_deg(pe_vec_value(neutral_delta_alpha,"
        "nearest_neutral_slot))"
    )
    df = df.Define(
        "nearest_neutral_p",
        "pe_vec_value(neutral_p,nearest_neutral_slot)"
    )
    df = df.Define(
        "nearest_neutral_beta",
        "pe_vec_value(neutral_beta,nearest_neutral_slot)"
    )
    df = df.Define(
        "nearest_neutral_detector",
        "pe_vec_ivalue(neutral_detector,nearest_neutral_slot)"
    )

    # Actual gamma2 numerator candidate: nearest PID-22 neutral candidate.
    df = df.Define(
        "probe_slot",
        "pe_best_neutral_slot(neutral_idx,neutral_pid,"
        "neutral_delta_alpha,tag_rec_index,true)"
    )
    df = df.Define("probe_exists", "probe_slot >= 0")
    df = df.Define(
        "probe_angle_deg",
        "pe_rad_to_deg(pe_vec_value(neutral_delta_alpha,probe_slot))"
    )
    df = df.Define(
        "probe_p",
        "pe_vec_value(neutral_p,probe_slot)"
    )
    df = df.Define(
        "probe_beta",
        "pe_vec_value(neutral_beta,probe_slot)"
    )
    df = df.Define(
        "probe_detector",
        "pe_vec_ivalue(neutral_detector,probe_slot)"
    )
    df = df.Define(
        "probe_rec_index",
        "pe_vec_ivalue(neutral_idx,probe_slot)"
    )
    df = df.Define(
        "probe_theta",
        "pe_vec_value(neutral_theta,probe_slot)"
    )
    df = df.Define(
        "probe_phi",
        "pe_vec_value(neutral_phi,probe_slot)"
    )
    df = df.Define(
        "e_gamma2_angle_deg",
        "probe_exists ? pe_angle_deg(e_theta,e_phi,probe_theta,probe_phi) : -999.0"
    )

    return df


def apply_cuts(df, cuts, prefix):
    """
    Apply a cut block sequentially and return:
      final dataframe
      [(label, count-result), ...]
    """
    out = df
    results = []
    for i, (label, expression) in enumerate(cuts):
        out = out.Filter(expression, f"{prefix}_{i}_{label}")
        results.append((label, out.Count()))
    return out, results


def make_hist(df, expression, name, title, nbins, xmin, xmax):
    return df.Histo1D((name, title, nbins, xmin, xmax), expression)


def draw_before_after_canvas(sample_dfs,
                             before_dfs,
                             after_dfs,
                             plots,
                             outfile,
                             before_label,
                             after_label):
    """
    Dashed = before; solid = after. Color identifies sample.

    All histogram actions for the canvas are booked first and then evaluated
    with RunGraphs, avoiding a separate TChain scan for each histogram.
    """
    nx = 3
    ny = int(math.ceil(len(plots) / nx))

    cname = "c_" + os.path.basename(outfile).replace(".", "_")
    canvas = ROOT.TCanvas(cname, "", 1800, 540 * ny)
    canvas.Divide(nx, ny, 0.002, 0.002)

    # Book everything before touching GetValue().
    booked = {}
    actions = []
    unique = str(abs(hash(outfile)))

    for ip, plot in enumerate(plots, start=1):
        expr, title, nbins, xmin, xmax, guides = plot
        booked[ip] = {}

        for sample, label in SAMPLES:
            if sample not in sample_dfs:
                continue

            hb_r = make_hist(
                before_dfs[sample], expr,
                f"h_b_{ip}_{sample}_{unique}",
                title, nbins, xmin, xmax
            )
            ha_r = make_hist(
                after_dfs[sample], expr,
                f"h_a_{ip}_{sample}_{unique}",
                title, nbins, xmin, xmax
            )
            booked[ip][sample] = (label, hb_r, ha_r)
            actions.extend([hb_r, ha_r])

    if actions:
        ROOT.RDF.RunGraphs(actions)

    keep = list(actions)

    for ip, plot in enumerate(plots, start=1):
        expr, title, nbins, xmin, xmax, guides = plot

        pad = canvas.cd(ip)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.035)
        pad.SetBottomMargin(0.12)
        pad.SetTopMargin(0.10)

        legend = ROOT.TLegend(0.48, 0.60, 0.94, 0.89)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.022)
        keep.append(legend)

        hist_pairs = []
        ymax = 0.0

        for sample, label in SAMPLES:
            if sample not in booked[ip]:
                continue

            label, hb_r, ha_r = booked[ip][sample]
            hb = hb_r.GetValue()
            ha = ha_r.GetValue()

            hb.SetLineColor(COLORS[sample])
            ha.SetLineColor(COLORS[sample])
            hb.SetLineStyle(2)
            ha.SetLineStyle(1)
            hb.SetLineWidth(2)
            ha.SetLineWidth(3)
            hb.SetStats(0)
            ha.SetStats(0)

            ib = hb.Integral(1, hb.GetNbinsX())
            ia = ha.Integral(1, ha.GetNbinsX())
            if ib > 0:
                hb.Scale(1.0 / ib)
            if ia > 0:
                ha.Scale(1.0 / ia)

            ymax = max(ymax, hb.GetMaximum(), ha.GetMaximum())
            hist_pairs.append((sample, label, hb, ha))
            keep.extend([hb, ha])

            legend.AddEntry(hb, f"{label}: {before_label}", "l")
            legend.AddEntry(ha, f"{label}: {after_label}", "l")

        first = True
        for sample, label, hb, ha in hist_pairs:
            for h in (hb, ha):
                h.SetMaximum(1.35 * ymax if ymax > 0 else 1.0)
                h.SetMinimum(0.0)
                h.GetXaxis().SetTitleSize(0.041)
                h.GetYaxis().SetTitleSize(0.038)
                h.GetXaxis().SetLabelSize(0.033)
                h.GetYaxis().SetLabelSize(0.033)
                h.GetXaxis().SetTitleOffset(1.15)
                h.GetYaxis().SetTitleOffset(1.45)
                h.Draw("HIST" if first else "HIST SAME")
                first = False

        for x in guides:
            line = ROOT.TLine(x, 0.0, x, 1.08 * ymax)
            line.SetLineColor(ROOT.kGray + 2)
            line.SetLineStyle(3)
            line.SetLineWidth(2)
            line.Draw("SAME")
            keep.append(line)

        legend.Draw()

    canvas.SaveAs(outfile)
    keep.append(canvas)
    return keep


def write_cutflow(path, initial_counts, prep_counts,
                  common_counts, common_probe_counts, probe_counts):
    with open(path, "w") as out:
        out.write("CLAS12 photon-efficiency Stage-1 cutflow\n")
        out.write("=======================================\n\n")

        out.write("Stage 1A: preparation\n")
        out.write("Stage 1B: common gamma1 cuts\n")
        out.write("Stage 1C: common inferred-probe X selection\n")
        out.write("Stage 1D: reconstructed PID-22 gamma2 selection\n\n")

        for sample, label in SAMPLES:
            if sample not in initial_counts:
                continue

            n0 = initial_counts[sample]
            out.write(f"[{label}]\n")
            out.write(f"Initial skim hypotheses: {n0}\n")

            previous = n0
            out.write("\n  Stage 1A -- preparation\n")
            for cut_label, count in prep_counts[sample]:
                frac = 100.0 * count / previous if previous else 0.0
                total = 100.0 * count / n0 if n0 else 0.0
                out.write(
                    f"    {cut_label:<46s} "
                    f"{count:12d}  "
                    f"{frac:8.3f}% step  {total:8.3f}% initial\n"
                )
                previous = count

            out.write("\n  Stage 1B -- common gamma1\n")
            for cut_label, count in common_counts[sample]:
                frac = 100.0 * count / previous if previous else 0.0
                total = 100.0 * count / n0 if n0 else 0.0
                out.write(
                    f"    {cut_label:<46s} "
                    f"{count:12d}  "
                    f"{frac:8.3f}% step  {total:8.3f}% initial\n"
                )
                previous = count

            out.write("\n  Stage 1C -- common inferred-probe selection\n")
            for cut_label, count in common_probe_counts[sample]:
                frac = 100.0 * count / previous if previous else 0.0
                total = 100.0 * count / n0 if n0 else 0.0
                out.write(
                    f"    {cut_label:<46s} "
                    f"{count:12d}  "
                    f"{frac:8.3f}% step  {total:8.3f}% initial\n"
                )
                previous = count

            out.write("\n  Stage 1D -- reconstructed gamma2 numerator\n")
            for cut_label, count in probe_counts[sample]:
                frac = 100.0 * count / previous if previous else 0.0
                total = 100.0 * count / n0 if n0 else 0.0
                out.write(
                    f"    {cut_label:<46s} "
                    f"{count:12d}  "
                    f"{frac:8.3f}% step  {total:8.3f}% initial\n"
                )
                previous = count

            out.write("\n\n")


def main():
    args = parse_args()

    if args.threads < 1:
        print("ERROR: --threads must be >= 1", file=sys.stderr)
        return 1

    # ROOT performs the event loops in compiled C++; Python only constructs
    # the computation graph. Keep the default deliberately conservative for
    # shared ifarm nodes.
    ROOT.EnableImplicitMT(args.threads)
    print(f"ROOT implicit multithreading: {args.threads} worker(s)")

    # This analysis owns its output directory: start every run clean so old
    # canvases/cutflows cannot be mistaken for products of the current cuts.
    if os.path.isdir(args.output_dir):
        for name in os.listdir(args.output_dir):
            path = os.path.join(args.output_dir, name)
            if os.path.isdir(path) and not os.path.islink(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
    else:
        os.makedirs(args.output_dir, exist_ok=True)
    print(f"Cleared output directory: {os.path.abspath(args.output_dir)}")

    chains = {}
    dfs = {}

    print(
        f"Looking for completed skim-version-{EXPECTED_SKIM_VERSION} files..."
    )

    for sample, label in SAMPLES:
        files = discover(sample, args.period)
        if not files:
            continue

        chain = make_chain(files)
        chains[sample] = chain
        dfs[sample] = make_dataframe(chain)

        print(
            f"  {label}: {len(files)} file(s), "
            f"{chain.GetEntries():,} PhotonEfficiency rows"
        )

    if not dfs:
        print("ERROR: no usable input files.", file=sys.stderr)
        return 1

    # -----------------------------------------------------------------------
    # Build the three stages.
    # -----------------------------------------------------------------------

    initial_count_handles = {}
    initial_counts = {}
    prep_dfs = {}
    prep_count_handles = {}

    common_dfs = {}
    common_count_handles = {}

    probe_base_dfs = {}
    common_probe_count_handles = {}

    probe_final_dfs = {}
    probe_count_handles = {}

    for sample, df in dfs.items():
        initial_count_handles[sample] = df.Count()

        prep_df, prep_handles = apply_cuts(
            df, PREP_CUTS, "prep"
        )
        prep_dfs[sample] = prep_df
        prep_count_handles[sample] = prep_handles

        common_df, common_handles = apply_cuts(
            prep_df, COMMON_PHOTON_CUTS, "common"
        )
        common_dfs[sample] = common_df
        common_count_handles[sample] = common_handles

        # Common inferred-probe selection: applies to denominator AND numerator.
        common_probe_df, common_probe_handles = apply_cuts(
            common_df, COMMON_PROBE_CUTS, "common_probe"
        )
        probe_base_dfs[sample] = common_probe_df
        common_probe_count_handles[sample] = common_probe_handles

        # Numerator only: require a reconstructed PID-22 gamma2 and apply
        # the explicitly chosen second-photon cuts.
        gamma2_base = common_probe_df.Filter(
            "probe_exists",
            "nearest_PID22_gamma2_exists"
        )
        probe_final, gamma2_handles = apply_cuts(
            gamma2_base, GAMMA2_CUTS, "gamma2"
        )
        probe_final_dfs[sample] = probe_final
        probe_count_handles[sample] = (
            [("nearest retained PID-22 gamma2 exists", gamma2_base.Count())]
            + gamma2_handles
        )

    # Materialize each sample's cutflow in one coordinated event loop rather
    # than triggering one scan per Count(). RunGraphs executes all booked
    # actions that share the same RDataFrame graph together.
    prep_counts = {}
    common_counts = {}
    common_probe_counts = {}
    probe_counts = {}

    for sample in dfs:
        actions = [initial_count_handles[sample]]
        actions += [h for _, h in prep_count_handles[sample]]
        actions += [h for _, h in common_count_handles[sample]]
        actions += [h for _, h in common_probe_count_handles[sample]]
        actions += [h for _, h in probe_count_handles[sample]]

        ROOT.RDF.RunGraphs(actions)

        initial_counts[sample] = int(initial_count_handles[sample].GetValue())
        prep_counts[sample] = [
            (label, int(handle.GetValue()))
            for label, handle in prep_count_handles[sample]
        ]
        common_counts[sample] = [
            (label, int(handle.GetValue()))
            for label, handle in common_count_handles[sample]
        ]
        common_probe_counts[sample] = [
            (label, int(handle.GetValue()))
            for label, handle in common_probe_count_handles[sample]
        ]
        probe_counts[sample] = [
            (label, int(handle.GetValue()))
            for label, handle in probe_count_handles[sample]
        ]

    # -----------------------------------------------------------------------
    # Canvas 1: before vs after preparation.
    # -----------------------------------------------------------------------

    prep_plots = [
        (
            "y",
            "Inelasticity y;y;Unit-normalized entries",
            100, 0.0, 1.0, [0.85]
        ),
        (
            "e_p",
            "Scattered-electron momentum;p_{e} (GeV);Unit-normalized entries",
            120, 0.0, 10.5, [2.0]
        ),
        (
            "Q2",
            "Q^{2};Q^{2} (GeV^{2});Unit-normalized entries",
            120, 0.0, 10.0, [1.0]
        ),
        (
            "W",
            "Invariant mass W;W (GeV);Unit-normalized entries",
            120, 1.0, 4.5, [2.0]
        ),
        (
            "Mx2_ep",
            "M^{2}_{X}(e'p');M^{2}_{X}(e'p') (GeV^{2});Unit-normalized entries",
            140, -1.0, 2.0, [-0.30, 0.40]
        ),
        (
            "Mx2_epg_raw",
            "Missing mass squared e'p'#gamma1;M^{2}_{X}(e'p'#gamma1) (GeV^{2});Unit-normalized entries",
            120, -0.3, 0.3, []
        ),
        (
            "Emiss_epg",
            "Missing energy e'p'#gamma1;E_{miss}(e'p'#gamma1) (GeV);Unit-normalized entries",
            160, -1.0, 10.0, []
        ),
    ]

    prep_file = os.path.join(
        args.output_dir,
        f"1_photon_efficiency_stage1_preparation_{args.period}.png"
    )

    keep_prep = draw_before_after_canvas(
        dfs, dfs, prep_dfs, prep_plots, prep_file,
        "before prep", "after prep"
    )

    # -----------------------------------------------------------------------
    # Canvas 2: preparation -> common gamma1 cuts.
    # -----------------------------------------------------------------------

    common_plots = [
        (
            "e_gamma1_angle_deg",
            "Electron-#gamma1 opening angle;#angle(e,#gamma1) (deg);Unit-normalized entries",
            160, 0.0, 120.0, [8.0]
        ),
        (
            "tag_corr_p",
            "#gamma1 momentum;p_{#gamma1} (GeV);Unit-normalized entries",
            150, 0.0, 10.0, []
        ),
        (
            "Mx2_ep",
            "Missing mass squared e'p';M^{2}_{X}(e'p') (GeV^{2});Unit-normalized entries",
            120, -0.30, 0.40, []
        ),
        (
            "Mx2_epg_raw",
            "Missing mass squared e'p'#gamma1;M^{2}_{X}(e'p'#gamma1) (GeV^{2});Unit-normalized entries",
            140, -0.5, 1.0, []
        ),
        (
            "Mx2_egamma1",
            "Missing mass squared e'#gamma1;M^{2}_{X}(e'#gamma1) (GeV^{2});Unit-normalized entries",
            140, 0.0, 8.0, [1.4]
        ),
    ]

    common_file = os.path.join(
        args.output_dir,
        f"2_photon_efficiency_stage1_common_photon_{args.period}.png"
    )

    keep_common = draw_before_after_canvas(
        dfs, prep_dfs, common_dfs, common_plots, common_file,
        "after prep", "after common #gamma1"
    )

    # -----------------------------------------------------------------------
    # Canvas 3: probe diagnostics and requested probe cuts.
    #
    # "Before" here means after Stage 1B AND after requiring that the relevant
    # diagnostic candidate exists. This avoids sentinel values dominating the
    # distributions.
    # -----------------------------------------------------------------------

    probe_plot_specs = [
        (
            "Mx2_ep",
            "Missing mass squared e'p';M^{2}_{X}(e'p') (GeV^{2});Unit-normalized entries",
            120, -0.5, 0.6,
            [PROBE_MX2_LOW, PROBE_MX2_HIGH]
        ),
        (
            "Mx2_egamma1",
            "Missing mass squared e'#gamma1;M^{2}_{X}(e'#gamma1) (GeV^{2});Unit-normalized entries",
            140, 0.0, 8.0, [1.4]
        ),
        (
            "probe_p",
            "Selected #gamma2 momentum;p_{#gamma2} (GeV);Unit-normalized entries",
            120, 0.0, 8.0, [0.4]
        ),
        (
            "probe_beta",
            "Selected #gamma2 velocity (diagnostic only);#beta_{#gamma2};Unit-normalized entries",
            120, 0.5, 1.5, [0.9, 1.1]
        ),
    ]

    probe_file = os.path.join(
        args.output_dir,
        f"3_photon_efficiency_stage1_probe_{args.period}.png"
    )

    gamma2_plot_base_dfs = {
        sample: probe_base_dfs[sample].Filter(
            "probe_exists", "gamma2_exists_for_probe_canvas"
        )
        for sample in dfs
    }

    keep_probe = draw_before_after_canvas(
        dfs, gamma2_plot_base_dfs, probe_final_dfs,
        probe_plot_specs, probe_file,
        "common probe + #gamma2", "after #gamma2 cuts"
    )

    # -----------------------------------------------------------------------
    # Extra diagnostic: nearest neutral regardless of PID.
    # This is intentionally NOT part of the numerator selection.
    # -----------------------------------------------------------------------

    nearest_before = {}
    nearest_after = {}
    for sample in dfs:
        nearest_before[sample] = probe_base_dfs[sample].Filter(
            "nearest_neutral_exists",
            "nearest_neutral_exists_for_diagnostic"
        )
        nearest_after[sample] = nearest_before[sample].Filter(
            "nearest_neutral_pid == 22",
            "nearest_neutral_is_PID22_for_diagnostic"
        )

    nearest_plots = [
        (
            "nearest_neutral_angle_deg",
            "Nearest neutral match to X;#angle(neutral,X) (deg);Unit-normalized entries",
            120, 0.0, 20.0,
            [PROBE_ANGLE_LOW, PROBE_ANGLE_HIGH]
        ),
        (
            "nearest_neutral_pid",
            "Nearest neutral PID;REC PID;Unit-normalized entries",
            80, -20.0, 2200.0, []
        ),
        (
            "nearest_neutral_p",
            "Nearest neutral momentum;p (GeV);Unit-normalized entries",
            120, 0.0, 8.0, []
        ),
        (
            "nearest_neutral_beta",
            "Nearest neutral velocity;#beta;Unit-normalized entries",
            120, 0.5, 1.5, []
        ),
        (
            "nearest_neutral_detector",
            "Nearest neutral detector;Detector code;Unit-normalized entries",
            4, -0.5, 3.5, []
        ),
        (
            "Mx2_epg_raw",
            "Missing mass squared e'p'#gamma1;M^{2}_{X}(e'p'#gamma1) (GeV^{2});Unit-normalized entries",
            140, -0.5, 1.0, []
        ),
    ]

    nearest_file = os.path.join(
        args.output_dir,
        f"4_photon_efficiency_stage1_nearest_neutral_{args.period}.png"
    )

    keep_nearest = draw_before_after_canvas(
        dfs, nearest_before, nearest_after,
        nearest_plots, nearest_file,
        "nearest neutral", "nearest neutral PID 22"
    )

    # -----------------------------------------------------------------------
    # Cutflow
    # -----------------------------------------------------------------------

    cutflow_file = os.path.join(
        args.output_dir,
        f"5_photon_efficiency_stage1_cutflow_{args.period}.txt"
    )

    write_cutflow(
        cutflow_file,
        initial_counts,
        prep_counts,
        common_counts,
        common_probe_counts,
        probe_counts
    )

    print("\nStage-1 outputs:")
    print(f"  {os.path.abspath(prep_file)}")
    print(f"  {os.path.abspath(common_file)}")
    print(f"  {os.path.abspath(probe_file)}")
    print(f"  {os.path.abspath(nearest_file)}")
    print(f"  {os.path.abspath(cutflow_file)}")

    print("\nNumerical Stage-1 common inferred-probe boundaries:")
    print(
        f"  {PROBE_ANGLE_LOW:.3f} < angle(gamma1,X) "
        f"< {PROBE_ANGLE_HIGH:.3f} deg"
    )
    print(
        f"  {PROBE_MX2_LOW:.3f} < Mx2(ep) "
        f"< {PROBE_MX2_HIGH:.3f} GeV^2"
    )
    print("  Mx2(e gamma1) > 1.400 GeV^2")

    print("\nNo MC normalization or efficiency calculation is performed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
