import ROOT
import os
import math
import sys
import multiprocessing as mp
import concurrent.futures

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# ------------------------------------------------
# runtime event limit
# default = first 5 million events from each tree
# optional command line override:
# python external_radiation_estimation.py 2000000
# python external_radiation_estimation.py all
#
# optional dataset override:
# python external_radiation_estimation.py 2000000 skip_mc
# python external_radiation_estimation.py all skip_mc
# ------------------------------------------------

max_events = 5000000
skip_mc = False

if len(sys.argv) > 1:
    if sys.argv[1].lower() == "all":
        max_events = -1
    else:
        max_events = int(sys.argv[1])
    #endif
#endif

if len(sys.argv) > 2:
    if sys.argv[2].lower() == "skip_mc":
        skip_mc = True
    #endif
#endif

if max_events < 0:
    print("Processing all events in each tree")
else:
    print("Processing up to %d events per tree" % max_events)
#endif

if skip_mc:
    print("Runtime override active: skipping MC dataset")
#endif

# ------------------------------------------------
# input files
# ------------------------------------------------

input_file_data_rga_epi = "/volatile/clas12/thayward/external_radiation_estimate/rga_fa18_inb_epi+.root"
input_file_mc_rga_epi = "/volatile/clas12/thayward/external_radiation_estimate/clasdis_rga_fa18_inb_epi+.root"
input_file_data_rga_epipim = "/volatile/clas12/thayward/external_radiation_estimate/rga_fa18_inb_epi+pi-.root"

input_file_data_rgc_epi = "/volatile/clas12/thayward/external_radiation_estimate/rgc_fa22_inb_epi+.root"
input_file_data_rgc_epipim = "/volatile/clas12/thayward/external_radiation_estimate/rgc_fa22_inb_epi+pi-.root"

tree_name = "PhysicsEvents"
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------
# beam / masses
# ------------------------------------------------

beam_energy = 10.6

m_e = 0.000511
m_p = 0.938272
m_n = 0.939565
m_pi = 0.13957

M_p2 = m_p * m_p
M_n2 = m_n * m_n

# ------------------------------------------------
# binning defaults
# RGA: -6.5 to 1.5
# RGC: -6.0 to 1.5
# ------------------------------------------------

vz_step = 0.5

# ------------------------------------------------
# sector definitions using e_phi in radians
# sector 1: > 5.9 or < 0.7
# sector 2: > 0.7 and < 1.7
# sector 3: > 1.7 and < 2.8
# sector 4: > 2.8 and < 3.9
# sector 5: > 3.9 and < 4.9
# sector 6: > 4.9 and < 5.9
# ------------------------------------------------

N_SECTORS = 6
INTEGRATED_INDEX = 0

sector_colors = [
    ROOT.kRed + 1,
    ROOT.kBlue + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
    ROOT.kCyan + 2
]

sector_labels = [
    "Sector 1",
    "Sector 2",
    "Sector 3",
    "Sector 4",
    "Sector 5",
    "Sector 6"
]

channel_overlay_colors = {
    "epi_plus": ROOT.kRed + 1,
    "epi_plus_pi_minus": ROOT.kBlue + 1
}

channel_overlay_labels = {
    "epi_plus": "e p #rightarrow e #pi^{+} X",
    "epi_plus_pi_minus": "e p #rightarrow e #pi^{+} #pi^{-} X"
}

# ------------------------------------------------
# photon-energy histogram binning
# now in MeV
# ------------------------------------------------

photon_energy_min = 0.0
photon_energy_max = 120.0
photon_energy_bins = 240

# ------------------------------------------------
# combined vz_e histogram binning
# ------------------------------------------------

vz_hist_min = -6.5
vz_hist_max = 1.5
vz_hist_bins = 160

# ------------------------------------------------
# branch configuration
# ------------------------------------------------

channel_configs = {
    "data_rga_epi_plus": {
        "electron": {
            "p": "e_p",
            "theta": "e_theta",
            "phi": "e_phi"
        },
        "hadrons": [
            {
                "name": "pip",
                "p": "p_p",
                "theta": "p_theta",
                "phi": "p_phi",
                "mass": m_pi
            }
        ],
        "observable_type": "Mx2",
        "observable_branch": "Mx2",
        "observable_min": 0.6,
        "observable_max": 1.2,
        "observable_bins": 45,
        "observable_axis_title": "M_{X}^{2} (GeV^{2})",
        "observable_peak_title": "M_{X}^{2}",
        "expected_peak": M_n2,
        "latex_label": "RGA Data: e p #rightarrow e #pi^{+} X",
        "do_energy_correction": True,
        "vz_min": -6.5,
        "vz_max": 1.5,
        "run_period": "RGA",
        "channel_key": "epi_plus"
    },
    "mc_rga_epi_plus": {
        "electron": {
            "p": "e_p",
            "theta": "e_theta",
            "phi": "e_phi"
        },
        "hadrons": [
            {
                "name": "pip",
                "p": "p_p",
                "theta": "p_theta",
                "phi": "p_phi",
                "mass": m_pi
            }
        ],
        "observable_type": "Mx2",
        "observable_branch": "Mx2",
        "observable_min": 0.6,
        "observable_max": 1.2,
        "observable_bins": 45,
        "observable_axis_title": "M_{X}^{2} (GeV^{2})",
        "observable_peak_title": "M_{X}^{2}",
        "expected_peak": M_n2,
        "latex_label": "RGA MC: e p #rightarrow e #pi^{+} X",
        "do_energy_correction": False,
        "vz_min": -6.5,
        "vz_max": 1.5,
        "run_period": "RGA",
        "channel_key": "epi_plus"
    },
    "data_rga_epi_plus_pi_minus": {
        "electron": {
            "p": "e_p",
            "theta": "e_theta",
            "phi": "e_phi"
        },
        "hadrons": [
            {
                "name": "pip",
                "p": "p1_p",
                "theta": "p1_theta",
                "phi": "p1_phi",
                "mass": m_pi
            },
            {
                "name": "pim",
                "p": "p2_p",
                "theta": "p2_theta",
                "phi": "p2_phi",
                "mass": m_pi
            }
        ],
        "observable_type": "Mx2",
        "observable_branch": "Mx2",
        "observable_min": 0.6,
        "observable_max": 1.2,
        "observable_bins": 45,
        "observable_axis_title": "M_{X}^{2} (GeV^{2})",
        "observable_peak_title": "M_{X}^{2}",
        "expected_peak": M_p2,
        "latex_label": "RGA Data: e p #rightarrow e #pi^{+} #pi^{-} X",
        "do_energy_correction": True,
        "vz_min": -6.5,
        "vz_max": 1.5,
        "run_period": "RGA",
        "channel_key": "epi_plus_pi_minus"
    },
    "data_rgc_epi_plus": {
        "electron": {
            "p": "e_p",
            "theta": "e_theta",
            "phi": "e_phi"
        },
        "hadrons": [
            {
                "name": "pip",
                "p": "p_p",
                "theta": "p_theta",
                "phi": "p_phi",
                "mass": m_pi
            }
        ],
        "observable_type": "Mx2",
        "observable_branch": "Mx2",
        "observable_min": 0.6,
        "observable_max": 1.2,
        "observable_bins": 45,
        "observable_axis_title": "M_{X}^{2} (GeV^{2})",
        "observable_peak_title": "M_{X}^{2}",
        "expected_peak": M_n2,
        "latex_label": "RGC Data: e p #rightarrow e #pi^{+} X",
        "do_energy_correction": True,
        "vz_min": -6.5,
        "vz_max": 1.5,
        "run_period": "RGC",
        "channel_key": "epi_plus"
    },
    "data_rgc_epi_plus_pi_minus": {
        "electron": {
            "p": "e_p",
            "theta": "e_theta",
            "phi": "e_phi"
        },
        "hadrons": [
            {
                "name": "pip",
                "p": "p1_p",
                "theta": "p1_theta",
                "phi": "p1_phi",
                "mass": m_pi
            },
            {
                "name": "pim",
                "p": "p2_p",
                "theta": "p2_theta",
                "phi": "p2_phi",
                "mass": m_pi
            }
        ],
        "observable_type": "Mx2",
        "observable_branch": "Mx2",
        "observable_min": 0.6,
        "observable_max": 1.2,
        "observable_bins": 45,
        "observable_axis_title": "M_{X}^{2} (GeV^{2})",
        "observable_peak_title": "M_{X}^{2}",
        "expected_peak": M_p2,
        "latex_label": "RGC Data: e p #rightarrow e #pi^{+} #pi^{-} X",
        "do_energy_correction": True,
        "vz_min": -6.5,
        "vz_max": 1.5,
        "run_period": "RGC",
        "channel_key": "epi_plus_pi_minus"
    }
}

# ------------------------------------------------
# numerical derivative step
# ------------------------------------------------

delta_epsilon = 0.001

# ------------------------------------------------
# helpers
# ------------------------------------------------

def make_safe_tag(text):
    out = text.lower()
    out = out.replace(" ", "_")
    out = out.replace("+", "p")
    out = out.replace("-", "m")
    out = out.replace("/", "_")
    out = out.replace("#", "")
    out = out.replace("{", "")
    out = out.replace("}", "")
    out = out.replace("^", "")
    return out
#enddef

def build_vz_edges(vz_min, vz_max, step):
    edges = []
    current_edge = vz_min
    while current_edge <= vz_max + 1.0e-9:
        edges.append(current_edge)
        current_edge += step
    #endfor
    return edges
#enddef

def print_progress(i_entry, n_entries, dataset_label):
    if n_entries <= 0:
        return
    #endif

    step = max(1, n_entries // 20)

    if i_entry == 0 or (i_entry + 1) % step == 0 or (i_entry + 1) == n_entries:
        frac = 100.0 * float(i_entry + 1) / float(n_entries)
        print("[%s] processed %d / %d entries (%.1f%%)" % (
            dataset_label,
            i_entry + 1,
            n_entries,
            frac
        ))
    #endif
#enddef

def verify_required_branches(tree, cfg, dataset_label):
    branch_names = set()
    branch_list = tree.GetListOfBranches()
    for i in range(branch_list.GetEntries()):
        branch_names.add(branch_list.At(i).GetName())
    #endfor

    required = []

    required.append(cfg["electron"]["p"])
    required.append(cfg["electron"]["theta"])
    required.append(cfg["electron"]["phi"])

    for had in cfg["hadrons"]:
        required.append(had["p"])
        required.append(had["theta"])
        required.append(had["phi"])
    #endfor

    required.append("vz_e")

    if cfg["observable_branch"] is not None:
        required.append(cfg["observable_branch"])
    #endif

    missing = []
    for name in required:
        if name not in branch_names:
            missing.append(name)
        #endif
    #endfor

    if len(missing) > 0:
        raise RuntimeError(
            "Missing required branches for %s: %s" % (
                dataset_label,
                ", ".join(missing)
            )
        )
    #endif
#enddef

def normalize_phi_0_2pi(phi_val):
    two_pi = 2.0 * math.pi
    phi_wrapped = math.fmod(phi_val, two_pi)
    if phi_wrapped < 0.0:
        phi_wrapped += two_pi
    #endif
    return phi_wrapped
#enddef

def get_clas12_sector_from_ephi(e_phi):
    phi = normalize_phi_0_2pi(e_phi)

    if phi > 5.9 or phi < 0.7:
        return 1
    elif phi > 0.7 and phi < 1.7:
        return 2
    elif phi > 1.7 and phi < 2.8:
        return 3
    elif phi > 2.8 and phi < 3.9:
        return 4
    elif phi > 3.9 and phi < 4.9:
        return 5
    elif phi > 4.9 and phi < 5.9:
        return 6
    else:
        return None
    #endif
#enddef

def make_lorentz_vector_from_p_theta_phi(pmag, theta_deg, phi_deg, mass):
    theta = math.radians(theta_deg)
    phi = math.radians(phi_deg)

    px = pmag * math.sin(theta) * math.cos(phi)
    py = pmag * math.sin(theta) * math.sin(phi)
    pz = pmag * math.cos(theta)
    energy = math.sqrt(pmag * pmag + mass * mass)

    lv = ROOT.TLorentzVector()
    lv.SetPxPyPzE(px, py, pz, energy)
    return lv
#enddef

def shifted_electron_lorentz_vector(lv_in, delta_e):
    p_old = lv_in.P()
    e_old = lv_in.E()

    if p_old <= 0.0:
        return None
    #endif

    e_new = e_old + delta_e
    if e_new <= m_e:
        return None
    #endif

    p_new2 = e_new * e_new - m_e * m_e
    if p_new2 <= 0.0:
        return None
    #endif

    p_new = math.sqrt(p_new2)
    scale = p_new / p_old

    lv_out = ROOT.TLorentzVector()
    lv_out.SetPxPyPzE(
        lv_in.Px() * scale,
        lv_in.Py() * scale,
        lv_in.Pz() * scale,
        e_new
    )
    return lv_out
#enddef

def build_event_four_vectors(tree, cfg):
    e_lv = make_lorentz_vector_from_p_theta_phi(
        float(getattr(tree, cfg["electron"]["p"])),
        float(getattr(tree, cfg["electron"]["theta"])),
        float(getattr(tree, cfg["electron"]["phi"])),
        m_e
    )

    hadron_lvs = []
    for had in cfg["hadrons"]:
        had_lv = make_lorentz_vector_from_p_theta_phi(
            float(getattr(tree, had["p"])),
            float(getattr(tree, had["theta"])),
            float(getattr(tree, had["phi"])),
            had["mass"]
        )
        hadron_lvs.append(had_lv)
    #endfor

    return e_lv, hadron_lvs
#enddef

def compute_mx2(e_lv, hadron_lvs):
    beam_lv = ROOT.TLorentzVector()
    beam_lv.SetPxPyPzE(0.0, 0.0, beam_energy, beam_energy)

    target_lv = ROOT.TLorentzVector()
    target_lv.SetPxPyPzE(0.0, 0.0, 0.0, m_p)

    miss_lv = beam_lv + target_lv - e_lv
    for had_lv in hadron_lvs:
        miss_lv = miss_lv - had_lv
    #endfor

    return miss_lv.M2()
#enddef

def compute_observable_value(tree, cfg, e_lv, hadron_lvs):
    if cfg["observable_branch"] is not None:
        return float(getattr(tree, cfg["observable_branch"]))
    #endif

    if cfg["observable_type"] == "Mx2":
        return compute_mx2(e_lv, hadron_lvs)
    else:
        raise RuntimeError("Unknown observable_type: %s" % cfg["observable_type"])
    #endif
#enddef

def compute_observable_from_shifted_electron(cfg, e_shifted, hadron_lvs):
    if cfg["observable_type"] == "Mx2":
        return compute_mx2(e_shifted, hadron_lvs)
    else:
        raise RuntimeError("Unknown observable_type: %s" % cfg["observable_type"])
    #endif
#enddef

def compute_event_dobservable_de(cfg, tree, e_lv, hadron_lvs):
    if cfg["observable_type"] == "Mx2":
        obs_nominal = compute_mx2(e_lv, hadron_lvs)
    else:
        raise RuntimeError("Unknown observable_type: %s" % cfg["observable_type"])
    #endif

    e_shifted = shifted_electron_lorentz_vector(e_lv, delta_epsilon)
    if e_shifted is None:
        return None
    #endif

    obs_shifted = compute_observable_from_shifted_electron(cfg, e_shifted, hadron_lvs)
    slope = (obs_shifted - obs_nominal) / delta_epsilon
    return slope
#enddef

def fit_is_reasonable(fit_func, fit_result, observable_min, observable_max):
    if fit_result is None:
        return False, -999
    #endif

    status = int(fit_result.Status())
    mu_val = fit_func.GetParameter(1)
    sigma_val = abs(fit_func.GetParameter(2))
    mu_err = fit_func.GetParError(1)
    sigma_err = fit_func.GetParError(2)

    fit_ok = False

    if sigma_val > 0.0 and observable_min < mu_val < observable_max:
        fit_ok = True
    #endif

    if not math.isfinite(mu_val) or not math.isfinite(sigma_val):
        fit_ok = False
    #endif

    if mu_err < 0.0 or sigma_err < 0.0:
        fit_ok = False
    #endif

    return fit_ok, status
#enddef

def get_fit_formula(dataset_label):
    label = dataset_label.lower()

    if "rga" in label:
        return "gaus(0)+pol2(3)"
    elif "rgc" in label:
        return "gaus(0)+pol1(3)"
    else:
        raise RuntimeError("Could not determine fit formula for dataset_label = %s" % dataset_label)
    #endif
#enddef

def initialize_fit_function(fit_func, dataset_label, peak_height, expected_peak):
    mu_low = 0.75
    mu_high = 1.00
    sigma_init = 0.03

    fit_func.SetParameter(0, peak_height)
    fit_func.SetParameter(1, expected_peak)
    fit_func.SetParameter(2, sigma_init)

    if "rga" in dataset_label.lower():
        fit_func.SetParameter(3, 0.0)
        fit_func.SetParameter(4, 0.0)
        fit_func.SetParameter(5, 0.0)
    elif "rgc" in dataset_label.lower():
        fit_func.SetParameter(3, 0.0)
        fit_func.SetParameter(4, 0.0)
    else:
        raise RuntimeError("Could not initialize fit parameters for dataset_label = %s" % dataset_label)
    #endif

    fit_func.SetParLimits(0, 0.0, 10.0)
    fit_func.SetParLimits(1, mu_low, mu_high)
    fit_func.SetParLimits(2, 0.003, 0.20)
#enddef

def book_histogram_groups(dataset_label, observable_bins, observable_min, observable_max, n_vz_bins):
    hist_groups = []

    for i_group in range(N_SECTORS + 1):
        this_group = []
        group_tag = "integrated" if i_group == INTEGRATED_INDEX else "sector_%d" % i_group

        for i_vz in range(n_vz_bins):
            h = ROOT.TH1D(
                "h_%s_%s_%d" % (make_safe_tag(dataset_label), group_tag, i_vz),
                "",
                observable_bins,
                observable_min,
                observable_max
            )
            h.Sumw2()
            h.SetDirectory(0)
            this_group.append(h)
        #endfor

        hist_groups.append(this_group)
    #endfor

    return hist_groups
#enddef

def initialize_accumulators(n_vz_bins):
    out = {
        "sum_vz": [0.0] * n_vz_bins,
        "count_vz": [0] * n_vz_bins,
        "sum_slope": [0.0] * n_vz_bins,
        "count_slope": [0] * n_vz_bins
    }
    return out
#enddef

def initialize_all_group_accumulators(n_vz_bins):
    accs = []
    for i_group in range(N_SECTORS + 1):
        accs.append(initialize_accumulators(n_vz_bins))
    #endfor
    return accs
#enddef

def update_group_accumulator(acc, vz_bin, vz_e, slope):
    acc["sum_vz"][vz_bin] += vz_e
    acc["count_vz"][vz_bin] += 1

    if slope is not None and math.isfinite(slope):
        acc["sum_slope"][vz_bin] += slope
        acc["count_slope"][vz_bin] += 1
    #endif
#enddef

def normalize_histogram_list(hist_list):
    for h in hist_list:
        integral = h.Integral()
        if integral > 0.0:
            h.Scale(1.0 / integral)
        #endif
    #endfor
#enddef

def build_graph_from_xy_lists(graph_xy):
    if graph_xy is None:
        return None
    #endif

    x_vals = graph_xy["x"]
    y_vals = graph_xy["y"]

    n = len(x_vals)
    g = ROOT.TGraph(n)
    for i in range(n):
        g.SetPoint(i, x_vals[i], y_vals[i])
    #endfor
    return g
#enddef

def build_graph_errors_from_xyerr_lists(graph_dict):
    x_vals = graph_dict["x"]
    y_vals = graph_dict["y"]
    y_errs = graph_dict["yerr"]

    n = len(x_vals)
    g = ROOT.TGraphErrors(n)
    for i in range(n):
        g.SetPoint(i, x_vals[i], y_vals[i])
        g.SetPointError(i, 0.0, y_errs[i])
    #endfor
    return g
#enddef

def style_integrated_graph(graph):
    graph.SetTitle("")
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.SetLineColor(ROOT.kBlack)
    graph.SetLineWidth(2)
#enddef

def style_sector_graph(graph, color, marker_style):
    graph.SetTitle("")
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(1.0)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetLineWidth(2)
#enddef

def fit_histogram_group(hist_list, acc, vz_edges, n_vz_bins, dataset_label, observable_min, observable_max, expected_peak, do_energy_correction, group_label_for_print):
    x_mu_vals = []
    y_mu_vals = []
    y_mu_errs = []

    x_sigma_vals = []
    y_sigma_vals = []
    y_sigma_errs = []

    fit_functions = []
    fit_statuses = []

    bin_mean_vz = []
    bin_avg_slope = []

    n_fit_success = 0

    print("")
    print("[%s] peak offsets for %s:" % (dataset_label, group_label_for_print))
    print("vz_low   vz_high   <vz_e>      mu           mu_err       offset=mu-ref")
    print("--------------------------------------------------------------")

    for i in range(n_vz_bins):
        h = hist_list[i]

        print("[%s] fitting %s vz bin %2d / %2d : %.1f to %.1f cm" % (
            dataset_label,
            group_label_for_print,
            i + 1,
            n_vz_bins,
            vz_edges[i],
            vz_edges[i + 1]
        ))

        if acc["count_vz"][i] < 25:
            print("[%s]   skipped: only %d entries in this vz bin" % (dataset_label, acc["count_vz"][i]))
            fit_functions.append(None)
            fit_statuses.append(False)
            continue
        #endif

        if h.GetEntries() < 10:
            print("[%s]   skipped: histogram has too few filled entries" % dataset_label)
            fit_functions.append(None)
            fit_statuses.append(False)
            continue
        #endif

        mean_vz = acc["sum_vz"][i] / float(acc["count_vz"][i])

        max_bin = h.GetMaximumBin()
        peak_height = h.GetBinContent(max_bin)

        fit_name = "fit_%s_%s_%d" % (
            make_safe_tag(dataset_label),
            make_safe_tag(group_label_for_print),
            i
        )
        fit_formula = get_fit_formula(dataset_label)
        fit_func = ROOT.TF1(fit_name, fit_formula, observable_min, observable_max)
        fit_func.SetNpx(500)

        initialize_fit_function(fit_func, dataset_label, peak_height, expected_peak)

        fit_result = h.Fit(fit_func, "RQ0S")
        fit_ok, fit_status = fit_is_reasonable(
            fit_func,
            fit_result,
            observable_min,
            observable_max
        )

        print("[%s]   fit formula = %s" % (dataset_label, fit_formula))
        print("[%s]   fit status = %d" % (dataset_label, fit_status))

        if fit_ok and fit_status != 0:
            print("[%s]   fit warning: nonzero Minuit status, but parameters look reasonable; accepting fit" % dataset_label)
        #endif

        fit_functions.append(fit_func)
        fit_statuses.append(fit_ok)

        if fit_ok:
            mu = fit_func.GetParameter(1)
            mu_err = fit_func.GetParError(1)
            sigma = abs(fit_func.GetParameter(2))
            sigma_err = fit_func.GetParError(2)
            offset = mu - expected_peak

            x_mu_vals.append(mean_vz)
            y_mu_vals.append(mu)
            y_mu_errs.append(mu_err)

            x_sigma_vals.append(mean_vz)
            y_sigma_vals.append(sigma)
            y_sigma_errs.append(sigma_err)

            if do_energy_correction and acc["count_slope"][i] > 0:
                avg_slope = acc["sum_slope"][i] / float(acc["count_slope"][i])
                bin_mean_vz.append(mean_vz)
                bin_avg_slope.append(avg_slope)
            #endif

            n_fit_success += 1

            print("[%s]   fit accepted: mu = %.6f, mu_err = %.6f, sigma = %.6f, sigma_err = %.6f" % (
                dataset_label,
                mu,
                mu_err,
                sigma,
                sigma_err
            ))

            print("%6.2f   %7.2f   %8.4f   %10.6f   %10.6f   %+12.6f" % (
                vz_edges[i],
                vz_edges[i + 1],
                mean_vz,
                mu,
                mu_err,
                offset
            ))

            if do_energy_correction and acc["count_slope"][i] > 0:
                print("[%s]   average d(Mx2)/dE in this bin = %.6f" % (
                    dataset_label,
                    avg_slope
                ))
            #endif
        else:
            print("[%s]   fit rejected: parameters not reasonable" % dataset_label)
        #endif
    #endfor

    print("[%s] successful fits for %s: %d / %d" % (
        dataset_label,
        group_label_for_print,
        n_fit_success,
        n_vz_bins
    ))

    graph_mu_dict = {
        "x": x_mu_vals,
        "y": y_mu_vals,
        "yerr": y_mu_errs
    }

    graph_sigma_dict = {
        "x": x_sigma_vals,
        "y": y_sigma_vals,
        "yerr": y_sigma_errs
    }

    cubic_coeffs = None

    if len(x_mu_vals) >= 4:
        graph_mu_tmp = build_graph_errors_from_xyerr_lists(graph_mu_dict)
        style_integrated_graph(graph_mu_tmp)

        cubic_fit = ROOT.TF1(
            "cubic_mu_%s_%s" % (make_safe_tag(dataset_label), make_safe_tag(group_label_for_print)),
            "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
            vz_edges[0],
            vz_edges[-1]
        )
        graph_mu_tmp.Fit(cubic_fit, "Q0")

        c0 = cubic_fit.GetParameter(0)
        c1 = cubic_fit.GetParameter(1)
        c2 = cubic_fit.GetParameter(2)
        c3 = cubic_fit.GetParameter(3)

        cubic_coeffs = [float(c0), float(c1), float(c2), float(c3)]

        print("")
        print("[%s] Cubic fit for mu(vz_e), %s:" % (dataset_label, group_label_for_print))
        print("mu(vz_e) = %.8f + %.8f*vz_e + %.8f*vz_e^2 + %.8f*vz_e^3" % (c0, c1, c2, c3))
    else:
        print("[%s] Not enough fitted points for cubic fit for %s." % (dataset_label, group_label_for_print))
    #endif

    correction_xy = None

    if do_energy_correction and cubic_coeffs is not None and len(bin_mean_vz) > 0:
        print("")
        print("[%s] Computing needed scattered-electron energy correction for %s relative to vz_e = %.2f cm" % (
            dataset_label,
            group_label_for_print,
            vz_edges[0]
        ))

        mu_ref = cubic_coeffs[0] + cubic_coeffs[1] * vz_edges[0] + cubic_coeffs[2] * vz_edges[0] * vz_edges[0] + cubic_coeffs[3] * vz_edges[0] * vz_edges[0] * vz_edges[0]

        x_corr_vals = []
        y_corr_vals = []

        print("vz_value     mu_fit(vz)     mu_ref-mu_fit     <dMx2/dE>     DeltaE_corr")
        print("-----------------------------------------------------------------------")

        for i in range(len(bin_mean_vz)):
            vz_val = bin_mean_vz[i]
            slope_val = bin_avg_slope[i]

            if abs(slope_val) < 1.0e-12:
                continue
            #endif

            mu_val = cubic_coeffs[0] + cubic_coeffs[1] * vz_val + cubic_coeffs[2] * vz_val * vz_val + cubic_coeffs[3] * vz_val * vz_val * vz_val
            delta_mu = mu_ref - mu_val
            delta_e_needed = delta_mu / slope_val

            x_corr_vals.append(vz_val)
            y_corr_vals.append(delta_e_needed)

            print("%7.3f   %12.6f   %13.6f   %12.6f   %12.6f" % (
                vz_val,
                mu_val,
                delta_mu,
                slope_val,
                delta_e_needed
            ))
        #endfor

        correction_xy = {
            "x": x_corr_vals,
            "y": y_corr_vals
        }
    #endif

    return {
        "fit_functions": fit_functions,
        "fit_statuses": fit_statuses,
        "graph_mu_dict": graph_mu_dict,
        "graph_sigma_dict": graph_sigma_dict,
        "cubic_coeffs": cubic_coeffs,
        "correction_xy": correction_xy,
        "y_sigma_vals": y_sigma_vals,
        "y_sigma_errs": y_sigma_errs
    }
#enddef

def draw_integrated_canvas(dataset_label, channel_latex, histograms, fit_results, vz_edges, vz_min, vz_max, expected_peak, observable_axis_title):
    n_vz_bins = len(histograms)

    graph_mu = build_graph_errors_from_xyerr_lists(fit_results["graph_mu_dict"])
    graph_sigma = build_graph_errors_from_xyerr_lists(fit_results["graph_sigma_dict"])
    style_integrated_graph(graph_mu)
    style_integrated_graph(graph_sigma)

    cubic_fit = None
    if fit_results["cubic_coeffs"] is not None:
        cubic_fit = ROOT.TF1(
            "cubic_draw_%s" % make_safe_tag(dataset_label),
            "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
            vz_min,
            vz_max
        )
        cubic_fit.SetParameters(
            fit_results["cubic_coeffs"][0],
            fit_results["cubic_coeffs"][1],
            fit_results["cubic_coeffs"][2],
            fit_results["cubic_coeffs"][3]
        )
        cubic_fit.SetLineColor(ROOT.kRed + 1)
        cubic_fit.SetLineStyle(2)
        cubic_fit.SetLineWidth(2)
    #endif

    canvas = ROOT.TCanvas(
        "canvas_%s" % make_safe_tag(dataset_label),
        dataset_label,
        2200,
        1300
    )

    left_margin_frac = 0.70
    left_pad = ROOT.TPad("left_pad_%s" % make_safe_tag(dataset_label), "", 0.00, 0.00, left_margin_frac, 1.00)
    right_pad = ROOT.TPad("right_pad_%s" % make_safe_tag(dataset_label), "", left_margin_frac, 0.00, 1.00, 1.00)

    left_pad.Draw()
    right_pad.Draw()

    left_pad.cd()

    n_left_cols = 4
    n_left_rows = int(math.ceil(float(n_vz_bins) / float(n_left_cols)))
    left_pad.Divide(n_left_cols, n_left_rows, 0.001, 0.001)

    for i in range(n_vz_bins):
        left_pad.cd(i + 1)

        ROOT.gPad.SetLeftMargin(0.19)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.14)
        ROOT.gPad.SetTopMargin(0.10)

        h = histograms[i]
        h.SetTitle("")
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.55)
        h.SetLineWidth(1)

        local_max = h.GetMaximum()
        y_max = 1.25 * local_max if local_max > 0.0 else 1.0

        h.SetMinimum(0.0)
        h.SetMaximum(y_max)

        h.GetXaxis().SetTitle(observable_axis_title)
        h.GetYaxis().SetTitle("Normalized counts")
        h.GetXaxis().CenterTitle()
        h.GetYaxis().CenterTitle()
        h.GetXaxis().SetTitleSize(0.06)
        h.GetYaxis().SetTitleSize(0.055)
        h.GetXaxis().SetLabelSize(0.05)
        h.GetYaxis().SetLabelSize(0.05)
        h.GetYaxis().SetTitleOffset(1.55)

        h.Draw("E1")

        if fit_results["fit_statuses"][i] and fit_results["fit_functions"][i]:
            fit_results["fit_functions"][i].SetLineColor(ROOT.kRed + 1)
            fit_results["fit_functions"][i].SetLineWidth(2)
            fit_results["fit_functions"][i].Draw("same")
        #endif

        vz_lo = vz_edges[i]
        vz_hi = vz_edges[i + 1]

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.050)
        latex.DrawLatex(0.22, 0.86, "%.1f < v_{z,e} < %.1f (cm)" % (vz_lo, vz_hi))

        if fit_results["fit_statuses"][i]:
            mu = fit_results["fit_functions"][i].GetParameter(1)
            sigma = abs(fit_results["fit_functions"][i].GetParameter(2))
            latex.SetTextSize(0.042)
            latex.DrawLatex(0.22, 0.78, "#mu = %.4f" % mu)
            latex.DrawLatex(0.22, 0.70, "#sigma = %.4f" % sigma)
        #endif
    #endfor

    right_pad.cd()

    right_top = ROOT.TPad("right_top_%s" % make_safe_tag(dataset_label), "", 0.00, 0.50, 1.00, 1.00)
    right_bot = ROOT.TPad("right_bot_%s" % make_safe_tag(dataset_label), "", 0.00, 0.00, 1.00, 0.50)

    right_top.Draw()
    right_bot.Draw()

    right_top.cd()
    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.08)

    frame_right_top = ROOT.TH1D("frame_right_top_%s" % make_safe_tag(dataset_label), "", 100, vz_min, vz_max)
    frame_right_top.SetMinimum(0.8)
    frame_right_top.SetMaximum(1.0)

    frame_right_top.GetXaxis().SetTitle("<v_{z,e}> in bin (cm)")
    frame_right_top.GetYaxis().SetTitle("Fitted peak position #mu of M_{X}^{2} (GeV^{2})")
    frame_right_top.GetXaxis().CenterTitle()
    frame_right_top.GetYaxis().CenterTitle()
    frame_right_top.GetXaxis().SetTitleSize(0.05)
    frame_right_top.GetYaxis().SetTitleSize(0.05)
    frame_right_top.GetXaxis().SetLabelSize(0.04)
    frame_right_top.GetYaxis().SetLabelSize(0.04)
    frame_right_top.GetYaxis().SetTitleOffset(1.35)
    frame_right_top.Draw()

    line_expected = ROOT.TLine(vz_min, expected_peak, vz_max, expected_peak)
    line_expected.SetLineStyle(2)
    line_expected.SetLineWidth(2)
    line_expected.Draw("same")

    graph_mu.Draw("P SAME")

    if cubic_fit:
        cubic_fit.Draw("same")
    #endif

    legend_right_top = ROOT.TLegend(0.16, 0.74, 0.94, 0.90)
    legend_right_top.SetBorderSize(1)
    legend_right_top.SetFillStyle(1001)
    legend_right_top.SetFillColor(ROOT.kWhite)
    legend_right_top.SetTextSize(0.023)
    legend_right_top.AddEntry(graph_mu, "Point = fitted #mu, error bar = fit uncertainty on #mu", "lep")
    legend_right_top.AddEntry(line_expected, "Reference peak position", "l")
    if cubic_fit:
        legend_right_top.AddEntry(cubic_fit, "Cubic fit to #mu(v_{z,e})", "l")
    #endif
    legend_right_top.Draw()

    latex_right_top = ROOT.TLatex()
    latex_right_top.SetNDC()
    latex_right_top.SetTextSize(0.038)
    latex_right_top.DrawLatex(0.16, 0.93, "%s: fitted #mu vs v_{z,e}" % channel_latex)

    right_bot.cd()
    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.08)

    sigma_ymax = 0.10
    for val, err in zip(fit_results["y_sigma_vals"], fit_results["y_sigma_errs"]):
        if val + err > sigma_ymax:
            sigma_ymax = val + err
        #endif
    #endfor
    sigma_ymax *= 1.25

    frame_right_bot = ROOT.TH1D("frame_right_bot_%s" % make_safe_tag(dataset_label), "", 100, vz_min, vz_max)
    frame_right_bot.SetMinimum(0.0)
    frame_right_bot.SetMaximum(sigma_ymax)
    frame_right_bot.GetXaxis().SetTitle("<v_{z,e}> in bin (cm)")
    frame_right_bot.GetYaxis().SetTitle("Fitted #sigma of M_{X}^{2} peak (GeV^{2})")
    frame_right_bot.GetXaxis().CenterTitle()
    frame_right_bot.GetYaxis().CenterTitle()
    frame_right_bot.GetXaxis().SetTitleSize(0.05)
    frame_right_bot.GetYaxis().SetTitleSize(0.05)
    frame_right_bot.GetXaxis().SetLabelSize(0.04)
    frame_right_bot.GetYaxis().SetLabelSize(0.04)
    frame_right_bot.GetYaxis().SetTitleOffset(1.35)
    frame_right_bot.Draw()

    graph_sigma.Draw("P SAME")

    legend_right_bot = ROOT.TLegend(0.17, 0.80, 0.77, 0.91)
    legend_right_bot.SetBorderSize(1)
    legend_right_bot.SetFillStyle(1001)
    legend_right_bot.SetFillColor(ROOT.kWhite)
    legend_right_bot.SetTextSize(0.028)
    legend_right_bot.AddEntry(graph_sigma, "Point = fitted #sigma, error bar = fit uncertainty", "lep")
    legend_right_bot.Draw()

    latex_right_bot = ROOT.TLatex()
    latex_right_bot.SetNDC()
    latex_right_bot.SetTextSize(0.038)
    latex_right_bot.DrawLatex(0.18, 0.95, "%s: fitted #sigma vs v_{z,e}" % channel_latex)

    output_png = os.path.join(output_dir, "external_radiation_estimation_%s.png" % make_safe_tag(dataset_label))
    print("Saving canvas to %s" % output_png)
    canvas.SaveAs(output_png)
    print("Saved: %s" % output_png)

    return canvas
#enddef

def draw_sector_canvas(dataset_label, channel_latex, sector_histograms, sector_fit_results, vz_edges, vz_min, vz_max, expected_peak, observable_axis_title):
    n_vz_bins = len(sector_histograms[0])

    graph_mus = []
    graph_sigmas = []

    for i_sector in range(N_SECTORS):
        g_mu = build_graph_errors_from_xyerr_lists(sector_fit_results[i_sector]["graph_mu_dict"])
        g_sigma = build_graph_errors_from_xyerr_lists(sector_fit_results[i_sector]["graph_sigma_dict"])
        style_sector_graph(g_mu, sector_colors[i_sector], 20 + i_sector)
        style_sector_graph(g_sigma, sector_colors[i_sector], 20 + i_sector)
        graph_mus.append(g_mu)
        graph_sigmas.append(g_sigma)
    #endfor

    canvas = ROOT.TCanvas(
        "canvas_sector_%s" % make_safe_tag(dataset_label),
        "sector_" + dataset_label,
        2200,
        1300
    )

    left_margin_frac = 0.70
    left_pad = ROOT.TPad("left_pad_sector_%s" % make_safe_tag(dataset_label), "", 0.00, 0.00, left_margin_frac, 1.00)
    right_pad = ROOT.TPad("right_pad_sector_%s" % make_safe_tag(dataset_label), "", left_margin_frac, 0.00, 1.00, 1.00)

    left_pad.Draw()
    right_pad.Draw()

    left_pad.cd()

    n_left_cols = 4
    n_left_rows = int(math.ceil(float(n_vz_bins) / float(n_left_cols)))
    left_pad.Divide(n_left_cols, n_left_rows, 0.001, 0.001)

    for i_vz in range(n_vz_bins):
        left_pad.cd(i_vz + 1)

        ROOT.gPad.SetLeftMargin(0.19)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.14)
        ROOT.gPad.SetTopMargin(0.14)

        local_max = 0.0
        for i_sector in range(N_SECTORS):
            h = sector_histograms[i_sector][i_vz]
            if h.GetMaximum() > local_max:
                local_max = h.GetMaximum()
            #endif
        #endfor

        y_max = 1.25 * local_max if local_max > 0.0 else 1.0
        first_drawn = False

        for i_sector in range(N_SECTORS):
            h = sector_histograms[i_sector][i_vz]
            h.SetTitle("")
            h.SetLineColor(sector_colors[i_sector])
            h.SetMarkerColor(sector_colors[i_sector])
            h.SetMarkerStyle(20 + i_sector)
            h.SetMarkerSize(0.40)
            h.SetLineWidth(2)
            h.SetMinimum(0.0)
            h.SetMaximum(y_max)

            h.GetXaxis().SetTitle(observable_axis_title)
            h.GetYaxis().SetTitle("Normalized counts")
            h.GetXaxis().CenterTitle()
            h.GetYaxis().CenterTitle()
            h.GetXaxis().SetTitleSize(0.06)
            h.GetYaxis().SetTitleSize(0.055)
            h.GetXaxis().SetLabelSize(0.05)
            h.GetYaxis().SetLabelSize(0.05)
            h.GetYaxis().SetTitleOffset(1.55)

            if not first_drawn:
                h.Draw("HIST")
                first_drawn = True
            else:
                h.Draw("HIST SAME")
            #endif

            if sector_fit_results[i_sector]["fit_statuses"][i_vz] and sector_fit_results[i_sector]["fit_functions"][i_vz]:
                sector_fit_results[i_sector]["fit_functions"][i_vz].SetLineColor(sector_colors[i_sector])
                sector_fit_results[i_sector]["fit_functions"][i_vz].SetLineWidth(2)
                sector_fit_results[i_sector]["fit_functions"][i_vz].Draw("same")
            #endif
        #endfor

        vz_lo = vz_edges[i_vz]
        vz_hi = vz_edges[i_vz + 1]

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.042)
        latex.DrawLatex(0.20, 0.93, "%.1f < v_{z,e} < %.1f (cm)" % (vz_lo, vz_hi))

        legend = ROOT.TLegend(0.20, 0.64, 0.42, 0.88)
        legend.SetBorderSize(1)
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.024)

        for i_sector in range(N_SECTORS):
            legend.AddEntry(sector_histograms[i_sector][i_vz], sector_labels[i_sector], "l")
        #endfor
        legend.Draw()
    #endfor

    right_pad.cd()

    right_top = ROOT.TPad("right_top_sector_%s" % make_safe_tag(dataset_label), "", 0.00, 0.50, 1.00, 1.00)
    right_bot = ROOT.TPad("right_bot_sector_%s" % make_safe_tag(dataset_label), "", 0.00, 0.00, 1.00, 0.50)

    right_top.Draw()
    right_bot.Draw()

    right_top.cd()
    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.08)

    frame_top = ROOT.TH1D("frame_sector_top_%s" % make_safe_tag(dataset_label), "", 100, vz_min, vz_max)
    frame_top.SetMinimum(0.8)
    frame_top.SetMaximum(1.0)
    frame_top.GetXaxis().SetTitle("<v_{z,e}> in bin (cm)")
    frame_top.GetYaxis().SetTitle("Fitted peak position #mu of M_{X}^{2} (GeV^{2})")
    frame_top.GetXaxis().CenterTitle()
    frame_top.GetYaxis().CenterTitle()
    frame_top.GetXaxis().SetTitleSize(0.05)
    frame_top.GetYaxis().SetTitleSize(0.05)
    frame_top.GetXaxis().SetLabelSize(0.04)
    frame_top.GetYaxis().SetLabelSize(0.04)
    frame_top.GetYaxis().SetTitleOffset(1.35)
    frame_top.Draw()

    line_expected = ROOT.TLine(vz_min, expected_peak, vz_max, expected_peak)
    line_expected.SetLineStyle(2)
    line_expected.SetLineWidth(2)
    line_expected.Draw("same")

    legend_top = ROOT.TLegend(0.15, 0.62, 0.93, 0.90)
    legend_top.SetBorderSize(1)
    legend_top.SetFillStyle(1001)
    legend_top.SetFillColor(ROOT.kWhite)
    legend_top.SetTextSize(0.024)
    legend_top.AddEntry(line_expected, "Reference peak position", "l")

    for i_sector in range(N_SECTORS):
        graph_mu = graph_mus[i_sector]
        graph_mu.Draw("P SAME")
        legend_top.AddEntry(graph_mu, sector_labels[i_sector], "lep")
    #endfor

    legend_top.Draw()

    latex_top = ROOT.TLatex()
    latex_top.SetNDC()
    latex_top.SetTextSize(0.036)
    latex_top.DrawLatex(0.16, 0.93, "%s: sector-by-sector fitted #mu vs v_{z,e}" % channel_latex)

    right_bot.cd()
    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.08)

    sigma_ymax = 0.10
    for i_sector in range(N_SECTORS):
        for val, err in zip(sector_fit_results[i_sector]["y_sigma_vals"], sector_fit_results[i_sector]["y_sigma_errs"]):
            if val + err > sigma_ymax:
                sigma_ymax = val + err
            #endif
        #endfor
    #endfor
    sigma_ymax *= 1.25

    frame_bot = ROOT.TH1D("frame_sector_bot_%s" % make_safe_tag(dataset_label), "", 100, vz_min, vz_max)
    frame_bot.SetMinimum(0.0)
    frame_bot.SetMaximum(sigma_ymax)
    frame_bot.GetXaxis().SetTitle("<v_{z,e}> in bin (cm)")
    frame_bot.GetYaxis().SetTitle("Fitted #sigma of M_{X}^{2} peak (GeV^{2})")
    frame_bot.GetXaxis().CenterTitle()
    frame_bot.GetYaxis().CenterTitle()
    frame_bot.GetXaxis().SetTitleSize(0.05)
    frame_bot.GetYaxis().SetTitleSize(0.05)
    frame_bot.GetXaxis().SetLabelSize(0.04)
    frame_bot.GetYaxis().SetLabelSize(0.04)
    frame_bot.GetYaxis().SetTitleOffset(1.35)
    frame_bot.Draw()

    legend_bot = ROOT.TLegend(0.15, 0.62, 0.70, 0.90)
    legend_bot.SetBorderSize(1)
    legend_bot.SetFillStyle(1001)
    legend_bot.SetFillColor(ROOT.kWhite)
    legend_bot.SetTextSize(0.026)

    for i_sector in range(N_SECTORS):
        graph_sigma = graph_sigmas[i_sector]
        graph_sigma.Draw("P SAME")
        legend_bot.AddEntry(graph_sigma, sector_labels[i_sector], "lep")
    #endfor

    legend_bot.Draw()

    latex_bot = ROOT.TLatex()
    latex_bot.SetNDC()
    latex_bot.SetTextSize(0.036)
    latex_bot.DrawLatex(0.16, 0.93, "%s: sector-by-sector fitted #sigma vs v_{z,e}" % channel_latex)

    output_png = os.path.join(output_dir, "external_radiation_estimation_sector_%s.png" % make_safe_tag(dataset_label))
    print("Saving sector canvas to %s" % output_png)
    canvas.SaveAs(output_png)
    print("Saved: %s" % output_png)

    return canvas
#enddef

def fill_photon_histogram_from_event_list(events_for_photon_hist, cubic_coeffs, vz_ref):
    h = ROOT.TH1D(
        "h_photon_energy_temp",
        "",
        photon_energy_bins,
        photon_energy_min,
        photon_energy_max
    )
    h.Sumw2()
    h.SetDirectory(0)

    if cubic_coeffs is None:
        return h
    #endif

    mu_ref = cubic_coeffs[0] + cubic_coeffs[1] * vz_ref + cubic_coeffs[2] * vz_ref * vz_ref + cubic_coeffs[3] * vz_ref * vz_ref * vz_ref

    for event_tuple in events_for_photon_hist:
        vz_e = event_tuple[0]
        slope = event_tuple[1]

        if slope is None:
            continue
        #endif

        if not math.isfinite(slope):
            continue
        #endif

        if abs(slope) < 1.0e-12:
            continue
        #endif

        mu_val = cubic_coeffs[0] + cubic_coeffs[1] * vz_e + cubic_coeffs[2] * vz_e * vz_e + cubic_coeffs[3] * vz_e * vz_e * vz_e
        delta_mu = mu_ref - mu_val
        delta_e_needed = delta_mu / slope
        delta_e_needed_mev = 1000.0 * delta_e_needed

        h.Fill(delta_e_needed_mev)
    #endfor

    return h
#enddef

def histogram_to_counts_list(hist):
    out = []
    nbins = hist.GetNbinsX()
    for i_bin in range(1, nbins + 1):
        out.append(float(hist.GetBinContent(i_bin)))
    #endfor
    return out
#enddef

def counts_list_to_hist(hist_name, counts, nbins, xmin, xmax):
    h = ROOT.TH1D(
        hist_name,
        "",
        nbins,
        xmin,
        xmax
    )
    h.SetDirectory(0)
    for i_bin in range(len(counts)):
        h.SetBinContent(i_bin + 1, counts[i_bin])
    #endfor
    return h
#enddef

def process_file_worker(task):
    input_file = task["input_file"]
    dataset_label = task["dataset_label"]

    cfg = channel_configs[dataset_label]
    expected_peak = cfg["expected_peak"]
    channel_latex = cfg["latex_label"]

    observable_min = cfg["observable_min"]
    observable_max = cfg["observable_max"]
    observable_bins = cfg["observable_bins"]
    observable_axis_title = cfg["observable_axis_title"]
    observable_peak_title = cfg["observable_peak_title"]

    vz_min = cfg["vz_min"]
    vz_max = cfg["vz_max"]
    vz_edges = build_vz_edges(vz_min, vz_max, vz_step)
    n_vz_bins = len(vz_edges) - 1

    print("")
    print("==============================================================")
    print("Starting dataset: %s" % dataset_label)
    print("Input file: %s" % input_file)
    print("Expected %s peak position = %.6f" % (observable_peak_title, expected_peak))
    print("Using vz_e range %.1f to %.1f cm" % (vz_min, vz_max))
    print("Opening file...")

    f = ROOT.TFile.Open(input_file, "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Could not open file: " + input_file)
    #endif

    print("Reading tree: %s" % tree_name)
    tree = f.Get(tree_name)
    if not tree:
        raise RuntimeError("Could not find tree: " + tree_name + " in " + input_file)
    #endif

    print("Verifying required branches...")
    verify_required_branches(tree, cfg, dataset_label)

    print("Booking %d vz_e histograms for integrated plus 6 sectors..." % n_vz_bins)
    hist_groups = book_histogram_groups(
        dataset_label,
        observable_bins,
        observable_min,
        observable_max,
        n_vz_bins
    )

    group_accs = initialize_all_group_accumulators(n_vz_bins)

    vz_hist = ROOT.TH1D(
        "h_vz_distribution_%s" % make_safe_tag(dataset_label),
        "",
        vz_hist_bins,
        vz_hist_min,
        vz_hist_max
    )
    vz_hist.Sumw2()
    vz_hist.SetDirectory(0)

    n_entries_total = tree.GetEntries()

    if max_events < 0:
        n_entries = n_entries_total
    else:
        n_entries = min(n_entries_total, max_events)
    #endif

    print("Total entries in tree: %d" % n_entries_total)
    print("Processing entries: %d" % n_entries)
    print("Beginning event loop...")

    n_kept_integrated = 0
    n_kept_sector = [0] * N_SECTORS

    events_for_photon_hist = []

    for i_entry in range(n_entries):
        tree.GetEntry(i_entry)
        print_progress(i_entry, n_entries, dataset_label)

        vz_e = float(tree.vz_e)

        if vz_e < vz_min or vz_e >= vz_max:
            continue
        #endif

        e_phi = float(getattr(tree, cfg["electron"]["phi"]))
        sector = get_clas12_sector_from_ephi(e_phi)
        if sector is None:
            continue
        #endif

        need_kinematics = cfg["do_energy_correction"] or cfg["observable_branch"] is None

        e_lv = None
        hadron_lvs = None
        slope = None

        if need_kinematics:
            e_lv, hadron_lvs = build_event_four_vectors(tree, cfg)
        #endif

        if cfg["observable_branch"] is not None:
            observable_val = float(getattr(tree, cfg["observable_branch"]))
        else:
            observable_val = compute_observable_value(tree, cfg, e_lv, hadron_lvs)
        #endif

        if observable_val < observable_min or observable_val > observable_max:
            continue
        #endif

        if cfg["do_energy_correction"]:
            slope = compute_event_dobservable_de(cfg, tree, e_lv, hadron_lvs)
        #endif

        vz_bin = int((vz_e - vz_min) / vz_step)

        if vz_bin < 0 or vz_bin >= n_vz_bins:
            continue
        #endif

        hist_groups[INTEGRATED_INDEX][vz_bin].Fill(observable_val)
        update_group_accumulator(group_accs[INTEGRATED_INDEX], vz_bin, vz_e, slope)
        n_kept_integrated += 1

        sector_index = sector
        hist_groups[sector_index][vz_bin].Fill(observable_val)
        update_group_accumulator(group_accs[sector_index], vz_bin, vz_e, slope)
        n_kept_sector[sector - 1] += 1

        vz_hist.Fill(vz_e)

        if cfg["do_energy_correction"]:
            events_for_photon_hist.append((float(vz_e), None if slope is None else float(slope)))
        #endif
    #endfor

    print("[%s] event loop complete" % dataset_label)
    print("[%s] integrated entries passing plotting cuts: %d" % (dataset_label, n_kept_integrated))
    for i_sector in range(N_SECTORS):
        print("[%s] %s entries passing plotting cuts: %d" % (
            dataset_label,
            sector_labels[i_sector],
            n_kept_sector[i_sector]
        ))
    #endfor

    print("Normalizing histograms...")
    for i_group in range(N_SECTORS + 1):
        normalize_histogram_list(hist_groups[i_group])
    #endfor

    print("Fitting integrated histograms in vz_e bins...")
    integrated_fit_results = fit_histogram_group(
        hist_groups[INTEGRATED_INDEX],
        group_accs[INTEGRATED_INDEX],
        vz_edges,
        n_vz_bins,
        dataset_label,
        observable_min,
        observable_max,
        expected_peak,
        cfg["do_energy_correction"],
        "integrated"
    )

    sector_fit_results = []
    for i_sector in range(N_SECTORS):
        print("Fitting sector histograms for %s..." % sector_labels[i_sector])
        sector_result = fit_histogram_group(
            hist_groups[i_sector + 1],
            group_accs[i_sector + 1],
            vz_edges,
            n_vz_bins,
            dataset_label,
            observable_min,
            observable_max,
            expected_peak,
            cfg["do_energy_correction"],
            sector_labels[i_sector]
        )
        sector_fit_results.append(sector_result)
    #endfor

    print("Building original integrated canvas...")
    draw_integrated_canvas(
        dataset_label,
        channel_latex,
        hist_groups[INTEGRATED_INDEX],
        integrated_fit_results,
        vz_edges,
        vz_min,
        vz_max,
        expected_peak,
        observable_axis_title
    )

    print("Building new sector overlay canvas...")
    draw_sector_canvas(
        dataset_label,
        channel_latex,
        hist_groups[1:],
        sector_fit_results,
        vz_edges,
        vz_min,
        vz_max,
        expected_peak,
        observable_axis_title
    )

    photon_hist_counts = None

    if cfg["do_energy_correction"] and integrated_fit_results["cubic_coeffs"] is not None:
        photon_hist = fill_photon_histogram_from_event_list(
            events_for_photon_hist,
            integrated_fit_results["cubic_coeffs"],
            vz_min
        )
        photon_hist.SetDirectory(0)

        photon_hist_counts = histogram_to_counts_list(photon_hist)

        photon_canvas = ROOT.TCanvas(
            "canvas_photon_energy_%s" % make_safe_tag(dataset_label),
            "photon_energy_" + dataset_label,
            1000,
            800
        )
        photon_canvas.cd()

        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.13)
        ROOT.gPad.SetTopMargin(0.08)

        local_max = photon_hist.GetMaximum()
        y_max = 1.25 * local_max if local_max > 0.0 else 1.0

        photon_hist.SetTitle("")
        photon_hist.SetLineColor(ROOT.kBlue + 1)
        photon_hist.SetLineWidth(2)
        photon_hist.SetMinimum(0.0)
        photon_hist.SetMaximum(y_max)
        photon_hist.GetXaxis().SetTitle("Required e^{-} energy gain (MeV)")
        photon_hist.GetYaxis().SetTitle("Counts")
        photon_hist.GetXaxis().CenterTitle()
        photon_hist.GetYaxis().CenterTitle()
        photon_hist.GetXaxis().SetTitleSize(0.05)
        photon_hist.GetYaxis().SetTitleSize(0.05)
        photon_hist.GetXaxis().SetLabelSize(0.04)
        photon_hist.GetYaxis().SetLabelSize(0.04)
        photon_hist.GetYaxis().SetTitleOffset(1.25)
        photon_hist.Draw("HIST")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.15, 0.93, "%s: event-by-event required e^{-} energy gain" % channel_latex)

        output_png = os.path.join(output_dir, "external_radiation_photon_energy_%s.png" % make_safe_tag(dataset_label))
        photon_canvas.SaveAs(output_png)
        print("Saved: %s" % output_png)
    #endif

    vz_hist_counts = histogram_to_counts_list(vz_hist)

    print("Finished dataset: %s" % dataset_label)
    f.Close()

    return {
        "dataset_label": dataset_label,
        "channel_latex": channel_latex,
        "run_period": cfg["run_period"],
        "channel_key": cfg["channel_key"],
        "do_energy_correction": cfg["do_energy_correction"],
        "integrated_correction_xy": integrated_fit_results["correction_xy"],
        "sector_correction_xys": [sector_fit_results[i]["correction_xy"] for i in range(N_SECTORS)],
        "photon_hist_counts": photon_hist_counts,
        "vz_hist_counts": vz_hist_counts,
        "vz_min": vz_min,
        "vz_max": vz_max
    }
#enddef

def make_combined_energy_correction_plot(results_by_label):
    graphs = []
    labels = []
    colors = []
    styles = []

    x_min = None
    x_max = None

    ordered_labels = [
        "data_rga_epi_plus",
        "data_rga_epi_plus_pi_minus",
        "data_rgc_epi_plus",
        "data_rgc_epi_plus_pi_minus"
    ]

    for dataset_label in ordered_labels:
        if dataset_label not in results_by_label:
            continue
        #endif

        result = results_by_label[dataset_label]
        if result["integrated_correction_xy"] is None:
            continue
        #endif

        graph = build_graph_from_xy_lists(result["integrated_correction_xy"])
        graphs.append(graph)

        if dataset_label == "data_rga_epi_plus":
            labels.append("RGA: e p #rightarrow e #pi^{+} X")
            colors.append(ROOT.kRed + 1)
            styles.append(1)
        elif dataset_label == "data_rga_epi_plus_pi_minus":
            labels.append("RGA: e p #rightarrow e #pi^{+} #pi^{-} X")
            colors.append(ROOT.kBlue + 1)
            styles.append(1)
        elif dataset_label == "data_rgc_epi_plus":
            labels.append("RGC: e p #rightarrow e #pi^{+} X")
            colors.append(ROOT.kRed + 1)
            styles.append(2)
        elif dataset_label == "data_rgc_epi_plus_pi_minus":
            labels.append("RGC: e p #rightarrow e #pi^{+} #pi^{-} X")
            colors.append(ROOT.kBlue + 1)
            styles.append(2)
        #endif

        if x_min is None or result["vz_min"] < x_min:
            x_min = result["vz_min"]
        #endif
        if x_max is None or result["vz_max"] > x_max:
            x_max = result["vz_max"]
        #endif
    #endfor

    if len(graphs) == 0:
        print("Skipping combined energy-correction plot because no correction graphs are available.")
        return
    #endif

    for i in range(len(graphs)):
        graphs[i].SetLineColor(colors[i])
        graphs[i].SetLineStyle(styles[i])
        graphs[i].SetLineWidth(3)
    #endfor

    min_y = 0.0
    max_y = 0.0
    first_point = True

    for graph in graphs:
        n = graph.GetN()
        for i in range(n):
            y_val = float(graph.GetPointY(i))

            if first_point:
                min_y = y_val
                max_y = y_val
                first_point = False
            else:
                if y_val < min_y:
                    min_y = y_val
                #endif
                if y_val > max_y:
                    max_y = y_val
                #endif
            #endif
        #endfor
    #endfor

    if first_point:
        min_y = 0.0
        max_y = 0.01
    #endif

    if min_y > 0.0:
        min_y = 0.0
    #endif

    y_span = max_y - min_y
    if y_span <= 0.0:
        y_span = 0.01
    #endif

    canvas = ROOT.TCanvas("canvas_energy_correction_data_only", "energy correction data only", 1100, 800)
    canvas.cd()

    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.13)
    ROOT.gPad.SetTopMargin(0.08)

    frame = ROOT.TH1D("frame_energy_correction_data_only", "", 100, x_min, x_max)
    frame.SetMinimum(min_y - 0.10 * y_span)
    frame.SetMaximum(max_y + 0.25 * y_span)
    frame.GetXaxis().SetTitle("v_{z,e} (cm)")
    frame.GetYaxis().SetTitle("e^{-} energy correction (GeV)")
    frame.GetXaxis().CenterTitle()
    frame.GetYaxis().CenterTitle()
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.Draw()

    for graph in graphs:
        graph.Draw("L SAME")
    #endfor

    legend = ROOT.TLegend(0.16, 0.66, 0.86, 0.90)
    legend.SetBorderSize(1)
    legend.SetFillStyle(1001)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextSize(0.028)

    for i in range(len(graphs)):
        legend.AddEntry(graphs[i], labels[i], "l")
    #endfor
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.040)
    latex.DrawLatex(0.16, 0.93, "External-radiation energy correction inferred from peak drift vs v_{z,e}")

    output_png = os.path.join(output_dir, "external_radiation_energy_correction_data_only.png")
    canvas.SaveAs(output_png)
    print("Saved: %s" % output_png)
#enddef

def draw_sector_energy_subplot(pad, result, frame_name, title_text):
    pad.cd()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.10)

    graphs = []
    x_min = result["vz_min"]
    x_max = result["vz_max"]

    min_y = None
    max_y = None

    for i_sector in range(N_SECTORS):
        graph_xy = result["sector_correction_xys"][i_sector]
        if graph_xy is None:
            continue
        #endif

        graph = build_graph_from_xy_lists(graph_xy)
        graph.SetLineColor(sector_colors[i_sector])
        graph.SetLineWidth(3)
        graphs.append((i_sector, graph))

        n = graph.GetN()
        for i in range(n):
            y_val = float(graph.GetPointY(i))
            if min_y is None or y_val < min_y:
                min_y = y_val
            #endif
            if max_y is None or y_val > max_y:
                max_y = y_val
            #endif
        #endfor
    #endfor

    if len(graphs) == 0:
        frame = ROOT.TH1D(frame_name, "", 100, x_min, x_max)
        frame.SetMinimum(0.0)
        frame.SetMaximum(0.01)
        frame.GetXaxis().SetTitle("v_{z,e} (cm)")
        frame.GetYaxis().SetTitle("e^{-} energy correction (GeV)")
        frame.GetXaxis().CenterTitle()
        frame.GetYaxis().CenterTitle()
        frame.Draw()
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.050)
        latex.DrawLatex(0.18, 0.90, title_text)
        latex.DrawLatex(0.25, 0.50, "No correction curves available")
        return
    #endif

    if min_y is None:
        min_y = 0.0
    #endif
    if max_y is None:
        max_y = 0.01
    #endif

    if min_y > 0.0:
        min_y = 0.0
    #endif

    y_span = max_y - min_y
    if y_span <= 0.0:
        y_span = 0.01
    #endif

    frame = ROOT.TH1D(frame_name, "", 100, x_min, x_max)
    frame.SetMinimum(min_y - 0.12 * y_span)
    frame.SetMaximum(max_y + 0.22 * y_span)
    frame.GetXaxis().SetTitle("v_{z,e} (cm)")
    frame.GetYaxis().SetTitle("e^{-} energy correction (GeV)")
    frame.GetXaxis().CenterTitle()
    frame.GetYaxis().CenterTitle()
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.Draw()

    for i_sector, graph in graphs:
        graph.Draw("L SAME")
    #endfor

    legend = ROOT.TLegend(0.17, 0.62, 0.87, 0.89)
    legend.SetBorderSize(1)
    legend.SetFillStyle(1001)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextSize(0.030)

    for i_sector, graph in graphs:
        legend.AddEntry(graph, sector_labels[i_sector], "l")
    #endfor
    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.16, 0.93, title_text)
#enddef

def make_sector_energy_correction_multipanel(results_by_label):
    ordered_labels = [
        "data_rga_epi_plus",
        "data_rga_epi_plus_pi_minus",
        "data_rgc_epi_plus",
        "data_rgc_epi_plus_pi_minus"
    ]

    labels = {
        "data_rga_epi_plus": "RGA: e p #rightarrow e #pi^{+} X",
        "data_rga_epi_plus_pi_minus": "RGA: e p #rightarrow e #pi^{+} #pi^{-} X",
        "data_rgc_epi_plus": "RGC: e p #rightarrow e #pi^{+} X",
        "data_rgc_epi_plus_pi_minus": "RGC: e p #rightarrow e #pi^{+} #pi^{-} X"
    }

    canvas = ROOT.TCanvas("canvas_sector_energy_correction_multipanel", "sector energy correction multipanel", 1600, 1200)
    canvas.Divide(2, 2, 0.001, 0.001)

    for i_panel in range(len(ordered_labels)):
        dataset_label = ordered_labels[i_panel]
        canvas.cd(i_panel + 1)

        if dataset_label not in results_by_label:
            ROOT.gPad.SetLeftMargin(0.15)
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetBottomMargin(0.14)
            ROOT.gPad.SetTopMargin(0.10)

            frame = ROOT.TH1D("frame_blank_sector_energy_%d" % i_panel, "", 100, -1.0, 1.0)
            frame.SetMinimum(0.0)
            frame.SetMaximum(1.0)
            frame.Draw()

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.055)
            latex.DrawLatex(0.18, 0.90, labels[dataset_label])
            latex.DrawLatex(0.30, 0.50, "No correction curves available")
            continue
        #endif

        result = results_by_label[dataset_label]

        if not result["do_energy_correction"]:
            ROOT.gPad.SetLeftMargin(0.15)
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetBottomMargin(0.14)
            ROOT.gPad.SetTopMargin(0.10)

            frame = ROOT.TH1D("frame_blank_sector_energy_%d" % i_panel, "", 100, result["vz_min"], result["vz_max"])
            frame.SetMinimum(0.0)
            frame.SetMaximum(0.01)
            frame.GetXaxis().SetTitle("v_{z,e} (cm)")
            frame.GetYaxis().SetTitle("e^{-} energy correction (GeV)")
            frame.Draw()

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.055)
            latex.DrawLatex(0.18, 0.90, labels[dataset_label])
            latex.DrawLatex(0.25, 0.50, "No correction curves available")
            continue
        #endif

        draw_sector_energy_subplot(
            ROOT.gPad,
            result,
            "frame_sector_energy_panel_%d" % i_panel,
            labels[dataset_label]
        )
    #endfor

    output_png = os.path.join(output_dir, "external_radiation_energy_correction_by_sector_multipanel.png")
    canvas.SaveAs(output_png)
    print("Saved: %s" % output_png)
#enddef

def make_photon_energy_multipanel(results_by_label):
    run_periods = []
    for dataset_label, result in results_by_label.items():
        if result["photon_hist_counts"] is None:
            continue
        #endif
        if result["run_period"] not in run_periods:
            run_periods.append(result["run_period"])
        #endif
    #endfor

    if len(run_periods) == 0:
        print("Skipping photon-energy multipanel because no photon histograms are available.")
        return
    #endif

    run_periods = sorted(run_periods)

    n_panels = len(run_periods)
    n_cols = 2 if n_panels > 1 else 1
    n_rows = int(math.ceil(float(n_panels) / float(n_cols)))

    canvas = ROOT.TCanvas("canvas_photon_energy_multipanel", "photon energy multipanel", 900 * n_cols, 700 * n_rows)
    canvas.Divide(n_cols, n_rows, 0.001, 0.001)

    for i_panel in range(n_panels):
        run_period = run_periods[i_panel]
        canvas.cd(i_panel + 1)

        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.13)
        ROOT.gPad.SetTopMargin(0.08)

        hists_to_draw = []
        local_max = 0.0

        for dataset_label, result in results_by_label.items():
            if result["run_period"] != run_period:
                continue
            #endif

            if result["photon_hist_counts"] is None:
                continue
            #endif

            channel_key = result["channel_key"]
            hist_name = "h_photon_mult_%s" % make_safe_tag(dataset_label)
            h = counts_list_to_hist(hist_name, result["photon_hist_counts"], photon_energy_bins, photon_energy_min, photon_energy_max)

            integral = h.Integral()
            if integral > 0.0:
                h.Scale(1.0 / integral)
            #endif

            h.SetLineColor(channel_overlay_colors[channel_key])
            h.SetLineWidth(3)

            if h.GetMaximum() > local_max:
                local_max = h.GetMaximum()
            #endif

            hists_to_draw.append((channel_key, h))
        #endfor

        if len(hists_to_draw) == 0:
            frame = ROOT.TH1D("frame_photon_empty_%d" % i_panel, "", 100, photon_energy_min, photon_energy_max)
            frame.SetMinimum(0.0)
            frame.SetMaximum(1.0)
            frame.GetXaxis().SetTitle("Required e^{-} energy gain (MeV)")
            frame.GetYaxis().SetTitle("Normalized counts")
            frame.Draw()

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.055)
            latex.DrawLatex(0.18, 0.90, "%s" % run_period)
            latex.DrawLatex(0.28, 0.50, "No photon-energy histograms available")
            continue
        #endif

        y_max = 1.25 * local_max if local_max > 0.0 else 1.0

        frame = ROOT.TH1D("frame_photon_%s" % make_safe_tag(run_period), "", 100, photon_energy_min, photon_energy_max)
        frame.SetMinimum(0.0)
        frame.SetMaximum(y_max)
        frame.GetXaxis().SetTitle("Required e^{-} energy gain (MeV)")
        frame.GetYaxis().SetTitle("Normalized counts")
        frame.GetXaxis().CenterTitle()
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().SetTitleOffset(1.25)
        frame.Draw()

        legend = ROOT.TLegend(0.16, 0.72, 0.86, 0.89)
        legend.SetBorderSize(1)
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.030)

        for channel_key, h in hists_to_draw:
            h.Draw("HIST SAME")
            legend.AddEntry(h, channel_overlay_labels[channel_key], "l")
        #endfor

        legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.DrawLatex(0.16, 0.93, "%s: event-by-event required e^{-} energy gain" % run_period)
    #endfor

    output_png = os.path.join(output_dir, "external_radiation_photon_energy_multipanel.png")
    canvas.SaveAs(output_png)
    print("Saved: %s" % output_png)
#enddef

def make_vz_distribution_plot(results_by_label):
    ordered_labels = [
        "data_rga_epi_plus",
        "data_rga_epi_plus_pi_minus",
        "data_rgc_epi_plus",
        "data_rgc_epi_plus_pi_minus",
        "mc_rga_epi_plus"
    ]

    display_labels = {
        "data_rga_epi_plus": "RGA Data: e p #rightarrow e #pi^{+} X",
        "data_rga_epi_plus_pi_minus": "RGA Data: e p #rightarrow e #pi^{+} #pi^{-} X",
        "data_rgc_epi_plus": "RGC Data: e p #rightarrow e #pi^{+} X",
        "data_rgc_epi_plus_pi_minus": "RGC Data: e p #rightarrow e #pi^{+} #pi^{-} X",
        "mc_rga_epi_plus": "RGA MC: e p #rightarrow e #pi^{+} X"
    }

    line_colors = {
        "data_rga_epi_plus": ROOT.kRed + 1,
        "data_rga_epi_plus_pi_minus": ROOT.kBlue + 1,
        "data_rgc_epi_plus": ROOT.kMagenta + 1,
        "data_rgc_epi_plus_pi_minus": ROOT.kGreen + 2,
        "mc_rga_epi_plus": ROOT.kBlack
    }

    line_styles = {
        "data_rga_epi_plus": 1,
        "data_rga_epi_plus_pi_minus": 1,
        "data_rgc_epi_plus": 2,
        "data_rgc_epi_plus_pi_minus": 2,
        "mc_rga_epi_plus": 1
    }

    hist_list = []
    local_max = 0.0

    for dataset_label in ordered_labels:
        if dataset_label not in results_by_label:
            continue
        #endif

        counts = results_by_label[dataset_label]["vz_hist_counts"]
        h = counts_list_to_hist(
            "h_vz_overlay_%s" % make_safe_tag(dataset_label),
            counts,
            vz_hist_bins,
            vz_hist_min,
            vz_hist_max
        )

        integral = h.Integral()
        if integral > 0.0:
            h.Scale(1.0 / integral)
        #endif

        h.SetLineColor(line_colors[dataset_label])
        h.SetLineStyle(line_styles[dataset_label])
        h.SetLineWidth(3)

        if h.GetMaximum() > local_max:
            local_max = h.GetMaximum()
        #endif

        hist_list.append((dataset_label, h))
    #endfor

    if len(hist_list) == 0:
        print("Skipping vz_e distribution plot because no histograms are available.")
        return
    #endif

    y_max = 1.25 * local_max if local_max > 0.0 else 1.0

    canvas = ROOT.TCanvas("canvas_vz_distributions", "vz distributions", 1100, 800)
    canvas.cd()

    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.13)
    ROOT.gPad.SetTopMargin(0.08)

    frame = ROOT.TH1D("frame_vz_distributions", "", 100, vz_hist_min, vz_hist_max)
    frame.SetMinimum(0.0)
    frame.SetMaximum(y_max)
    frame.GetXaxis().SetTitle("v_{z,e} (cm)")
    frame.GetYaxis().SetTitle("Normalized counts")
    frame.GetXaxis().CenterTitle()
    frame.GetYaxis().CenterTitle()
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.Draw()

    legend = ROOT.TLegend(0.16, 0.63, 0.88, 0.89)
    legend.SetBorderSize(1)
    legend.SetFillStyle(1001)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextSize(0.028)

    for dataset_label, h in hist_list:
        h.Draw("HIST SAME")
        legend.AddEntry(h, display_labels[dataset_label], "l")
    #endfor

    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.042)
    latex.DrawLatex(0.16, 0.93, "v_{z,e} distributions for all processed datasets")

    output_png = os.path.join(output_dir, "external_radiation_vz_distributions.png")
    canvas.SaveAs(output_png)
    print("Saved: %s" % output_png)
#enddef

def run_all():
    tasks = [
        {
            "input_file": input_file_data_rga_epi,
            "dataset_label": "data_rga_epi_plus"
        },
        {
            "input_file": input_file_data_rga_epipim,
            "dataset_label": "data_rga_epi_plus_pi_minus"
        },
        {
            "input_file": input_file_data_rgc_epi,
            "dataset_label": "data_rgc_epi_plus"
        },
        {
            "input_file": input_file_data_rgc_epipim,
            "dataset_label": "data_rgc_epi_plus_pi_minus"
        }
    ]

    if not skip_mc:
        tasks.append({
            "input_file": input_file_mc_rga_epi,
            "dataset_label": "mc_rga_epi_plus"
        })
    #endif

    max_workers = min(5, len(tasks))
    print("")
    print("Launching parallel processing with %d worker(s) (capped at 5)" % max_workers)

    results_by_label = {}

    mp_ctx = mp.get_context("spawn")

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers, mp_context=mp_ctx) as executor:
        future_to_label = {}
        for task in tasks:
            future = executor.submit(process_file_worker, task)
            future_to_label[future] = task["dataset_label"]
        #endfor

        for future in concurrent.futures.as_completed(future_to_label):
            dataset_label = future_to_label[future]
            result = future.result()
            results_by_label[dataset_label] = result
            print("Completed dataset: %s" % dataset_label)
        #endfor
    #endwith

    make_combined_energy_correction_plot(results_by_label)
    make_sector_energy_correction_multipanel(results_by_label)
    make_photon_energy_multipanel(results_by_label)
    make_vz_distribution_plot(results_by_label)

    print("")
    print("All processing complete.")
#enddef

if __name__ == "__main__":
    run_all()
#endif