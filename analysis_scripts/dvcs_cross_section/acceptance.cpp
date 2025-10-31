// acceptance.cpp
#include "acceptance.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

// ---------------- style bootstrap ----------------
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);
        const int rf = 42;
        gStyle->SetTitleFont(rf, "XYZ");
        gStyle->SetLabelFont(rf, "XYZ");
        gStyle->SetTextFont(rf);
    }
} _style_guard;

// ---------------- helpers ----------------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

static std::string periodToRunTagKey(const std::string& period) {
    // "DVCS_Fa18_inb" -> "fa18_inb"
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}

static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

static inline int findIndex(const std::pair<double,double>& range,
                            const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < (int)ranges.size(); ++i) if (ranges[i] == range) return i;
    return -1;
}

static inline int phiBinIndex(double phi) {
    double w = std::fmod(phi, TWO_PI);
    if (w < 0) w += TWO_PI;
    const double width = TWO_PI / double(N_PHI_BINS);
    int ip = int(std::floor(w / width));
    if (ip < 0) ip = 0;
    if (ip >= N_PHI_BINS) ip = N_PHI_BINS - 1;
    return ip;
}

static inline int findBin1D(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int k = 0; k < (int)ranges.size(); ++k) if (v >= ranges[k].first && v < ranges[k].second) return k;
    return -1;
}

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
    for (int i = 0; i < N_PHI_BINS; ++i) d[i] = (i + 0.5) * step;
    return d;
}

static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

// ---------------- branches ----------------
struct GenBranch {
    double x = 0.0, Q2 = 0.0, t1 = 0.0, phi2 = 0.0;
    bool has_x = false, has_Q2 = false, has_t1 = false, has_phi = false;
    void bind(TTree* t) {
        auto bindD = [&](const char* n, double* a, bool& f){ if (t && t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };
        bindD("x",   &x,   has_x);
        bindD("Q2",  &Q2,  has_Q2);
        bindD("t1",  &t1,  has_t1);
        bindD("phi2",&phi2,has_phi);
    }
};

// Reconstructed MC needs all variables that exclusivity_cuts used
struct RecBranch {
    double x = 0.0, Q2 = 0.0, t1 = 0.0, phi2 = 0.0;
    double open_angle_ep2 = 0.0, pTmiss = 0.0;
    int detector1 = 0, detector2 = 0;
    // exclusivity variables (bind if present)
    double Delta_phi = 0.0;
    double theta_gamma_gamma = 0.0;
    double xF = 0.0;
    double Emiss2 = 0.0;
    double Mx2 = 0.0;
    double Mx2_1 = 0.0;
    double Mx2_2 = 0.0;

    bool has_x=false, has_Q2=false, has_t1=false, has_phi=false;
    bool hasTopo=false, hasOA=false, hasPT=false;
    bool has_Delta_phi=false, has_theta_gg=false, has_xF=false, has_Emiss2=false, has_Mx2=false, has_Mx2_1=false, has_Mx2_2=false;

    void bind(TTree* t) {
        auto bindD = [&](const char* n, double* a, bool& f){ if (t && t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };
        auto bindI = [&](const char* n, int* a, bool& f){ if (t && t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };

        bindD("x",   &x,   has_x);
        bindD("Q2",  &Q2,  has_Q2);
        bindD("t1",  &t1,  has_t1);
        bindD("phi2",&phi2,has_phi);
        bindD("open_angle_ep2",&open_angle_ep2,hasOA);
        bindD("pTmiss",&pTmiss,hasPT);
        bindI("detector1",&detector1,hasTopo);
        bindI("detector2",&detector2,hasTopo);

        // extras used by exclusivity_cuts
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);
        bindD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gg);
        bindD("xF", &xF, has_xF);
        bindD("Emiss2", &Emiss2, has_Emiss2);
        bindD("Mx2", &Mx2, has_Mx2);
        bindD("Mx2_1", &Mx2_1, has_Mx2_1);
        bindD("Mx2_2", &Mx2_2, has_Mx2_2);
    }

    // Build value map in the same naming as exclusivity_cuts (DVCS channel)
    std::map<std::string,double> valuesMapDVCS() const {
        std::map<std::string,double> m;
        if (has_Delta_phi) m["Delta_phi"] = Delta_phi;
        if (has_theta_gg)  m["theta_gamma_gamma"] = theta_gamma_gamma;
        if (hasPT)         m["pTmiss"] = pTmiss;
        if (has_xF)        m["xF"] = xF;
        if (has_Emiss2)    m["Emiss2"] = Emiss2;
        if (has_Mx2)       m["Mx2"] = Mx2;
        if (has_Mx2_1)     m["Mx2_1"] = Mx2_1;
        if (has_Mx2_2)     m["Mx2_2"] = Mx2_2;
        return m;
    }
};

// ---------------- global DVCS MC-side kinematic cuts (same as exclusivity_cuts) ----------------
struct MCCutsFixed {
    double min_open_angle_deg = 5.0;
    double max_neg_t_GeV2     = 1.0;
    double max_pTmiss_GeV     = 0.20;
};
static inline bool passesGlobalMCCuts(double t1, double oa_deg, double pTmiss, const MCCutsFixed& c) {
    if (oa_deg <= c.min_open_angle_deg) return false;
    if ((-t1) > c.max_neg_t_GeV2)       return false;
    if (pTmiss > c.max_pTmiss_GeV)      return false;
    return true;
}

// ---------------- JSON: load combined_cuts.json (MC stats only) ----------------
struct Stats { double mean=0.0, std=0.0; };
using CutMap  = std::map<std::string, Stats>;          // var -> (mean,std)
using ComboMC = std::map<std::string, CutMap>;         // "DVCS_Fa18_inb_FD_FD" -> CutMap

// crude extractor: find number after a pattern
static bool extract_number_after(const std::string& s, size_t from, const std::string& pat, double& out) {
    size_t p = s.find(pat, from); if (p == std::string::npos) return false;
    p = s.find(':', p); if (p == std::string::npos) return false;
    size_t a = p + 1;
    while (a < s.size() && isspace((unsigned char)s[a])) ++a;
    size_t b = a;
    while (b < s.size() && (isdigit((unsigned char)s[b]) || s[b]=='-' || s[b]=='+' || s[b]=='.' || s[b]=='e' || s[b]=='E')) ++b;
    try { out = std::stod(s.substr(a, b - a)); return true; } catch(...) { return false; }
}

static bool extract_stats_block(const std::string& s, size_t scope_from, const std::string& var, Stats& st) {
    // look for ..."var":{"mean":X,"std":Y}... inside the "mc":{...} object nearest to scope_from
    // 1) find "mc" after scope_from
    size_t pmc = s.find("\"mc\"", scope_from); if (pmc == std::string::npos) return false;
    size_t lbrace = s.find('{', pmc); if (lbrace == std::string::npos) return false;
    // naive bracket matching for mc object
    int depth = 0; size_t i = lbrace;
    for (; i < s.size(); ++i) {
        if (s[i] == '{') ++depth;
        else if (s[i] == '}') { --depth; if (depth == 0) break; }
    }
    if (i >= s.size()) return false;
    size_t mc_end = i;
    // 2) within mc object, find the variable key
    size_t pv = s.find("\"" + var + "\"", lbrace);
    if (pv == std::string::npos || pv > mc_end) return false;
    // 3) extract mean and std near this variable
    double mu=0.0, sg=0.0;
    bool ok1 = extract_number_after(s, pv, "\"mean\"", mu);
    bool ok2 = extract_number_after(s, pv, "\"std\"",  sg);
    if (!(ok1 && ok2)) return false;
    st.mean = mu; st.std = std::abs(sg);
    return true;
}

// Produce possible period variants to tolerate case differences in JSON keys.
// Examples:
//   "DVCS_Sp18_inb" -> {"DVCS_Sp18_inb", "DVCS_Sp18_Inb"}
//   "DVCS_Fa18_out" -> {"DVCS_Fa18_out", "DVCS_Fa18_Out"}
static std::vector<std::string> period_variants(const std::string& period) {
    std::vector<std::string> v;
    v.push_back(period);
    // If ends with _inb or _out, add title-case variant
    if (period.size() >= 4) {
        if (period.rfind("_inb") != std::string::npos) {
            std::string p2 = period;
            p2.replace(p2.size()-4, 4, "_Inb");
            v.push_back(p2);
        } else if (period.rfind("_out") != std::string::npos) {
            std::string p2 = period;
            p2.replace(p2.size()-4, 4, "_Out");
            v.push_back(p2);
        }
    }
    return v;
}

static bool load_combined_mc_cuts(const std::string& json_path, const std::vector<std::string>& periods, ComboMC& out) {
    std::ifstream ifs(json_path);
    if (!ifs) {
        std::cerr << "[acc][ERROR] Cannot open combined cuts JSON: " << json_path << "\n";
        return false;
    }
    const std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    static const std::vector<std::string> topoKeys = {"FD_FD", "CD_FD", "CD_FT"};
    static const std::vector<std::string> vars = {
        "Delta_phi", "theta_gamma_gamma", "pTmiss", "xF",
        "Emiss2", "Mx2", "Mx2_1", "Mx2_2"
    };

    bool any=false;
    for (const auto& period : periods) {
        const auto pvars = period_variants(period); // try both lower and title case suffixes
        for (const auto& period_try : pvars) {
            for (const auto& topo : topoKeys) {
                const std::string comboKey = period_try + "_" + topo; // e.g. "DVCS_Sp18_Inb_FD_FD"
                size_t pkey = s.find("\"" + comboKey + "\"");
                if (pkey == std::string::npos) continue;

                CutMap cmap;
                for (const auto& v : vars) {
                    Stats st;
                    if (extract_stats_block(s, pkey, v, st)) {
                        cmap[v] = st;
                    }
                }
                if (!cmap.empty()) {
                    out[comboKey] = std::move(cmap);
                    any = true;
                }
            }
        }
    }
    return any;
}

static inline bool within3Sigma(double val, const Stats& s) {
    return (val >= s.mean - 3.0*s.std) && (val <= s.mean + 3.0*s.std);
}

static bool passes3SigmaCuts(const CutMap& cuts_mc, const std::map<std::string,double>& values) {
    for (const auto& kv : cuts_mc) {
        const std::string& var = kv.first;
        auto it = values.find(var);
        if (it == values.end()) continue; // variable not present in tree: skip
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// ---------------- topology helpers ----------------
static inline std::string topoKeyFromDet(int d1, int d2) {
    if (d1 == 1 && d2 == 1) return "FD_FD";
    if (d1 == 2 && d2 == 1) return "CD_FD";
    if (d1 == 2 && d2 == 0) return "CD_FT";
    return ""; // unknown
}

static inline bool topoAllowed(int d1, int d2, const std::vector<std::string>& tops) {
    for (const auto& t : tops) {
        if (t == "(FD,FD)" && d1 == 1 && d2 == 1) return true;
        if (t == "(CD,FD)" && d1 == 2 && d2 == 1) return true;
        if (t == "(CD,FT)" && d1 == 2 && d2 == 0) return true;
    }
    return false;
}

// ---------------- accumulation ----------------
struct AccBin {
    double gen = 0.0;
    double rec = 0.0;
};
using AccMap = std::map<std::tuple<int,int,int,int>, AccBin>; // (ix,iQ,it,ip)

static void accumulate_generated(
    TTree* t,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    AccMap& acc)
{
    if (!t) return;
    GenBranch b; b.bind(t);
    if (!(b.has_x && b.has_Q2 && b.has_t1 && b.has_phi)) return;

    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        double xB = b.x, Q2 = b.Q2, tt = std::fabs(b.t1), phi = b.phi2;
        int ix = findBin1D(xB, xB_bins), iQ = findBin1D(Q2, Q2_bins), it = findBin1D(tt, t_bins);
        if (ix < 0 || iQ < 0 || it < 0) continue;
        int ip = phiBinIndex(phi);
        acc[std::make_tuple(ix, iQ, it, ip)].gen += 1.0;
    }
}

static void accumulate_reconstructed(
    TTree* t,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::vector<std::string>& topologies,
    const MCCutsFixed& globalCuts,
    const std::map<std::string, CutMap>& topoCutsForPeriod, // e.g. {"FD_FD":CutMap,...}
    AccMap& acc)
{
    if (!t) return;
    RecBranch b; b.bind(t);
    if (!(b.has_x && b.has_Q2 && b.has_t1 && b.has_phi && b.hasTopo && b.hasOA && b.hasPT)) return;

    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);

        if (!topoAllowed(b.detector1, b.detector2, topologies)) continue;
        if (!passesGlobalMCCuts(b.t1, b.open_angle_ep2, b.pTmiss, globalCuts)) continue;

        // Topology-specific 3σ windows from combined_cuts.json (MC)
        const std::string topoKey = topoKeyFromDet(b.detector1, b.detector2);
        if (topoKey.empty()) continue; // unknown combo

        auto itCuts = topoCutsForPeriod.find(topoKey);
        if (itCuts != topoCutsForPeriod.end()) {
            const CutMap& cuts_mc = itCuts->second;
            if (!passes3SigmaCuts(cuts_mc, b.valuesMapDVCS())) continue;
        } else {
            // If this period had no entry for this topology, veto by default to avoid
            // mixing uncut with cut selections. Comment the next line if you prefer permissive.
            continue;
        }

        double xB = b.x, Q2 = b.Q2, tt = std::fabs(b.t1), phi = b.phi2;
        int ix = findBin1D(xB, xB_bins), iQ = findBin1D(Q2, Q2_bins), it = findBin1D(tt, t_bins);
        if (ix < 0 || iQ < 0 || it < 0) continue;
        int ip = phiBinIndex(phi);
        acc[std::make_tuple(ix, iQ, it, ip)].rec += 1.0;
    }
}

// ---------------- JSON writer ----------------
struct PhiArrays {
    std::vector<double> phi_deg;
    std::vector<double> acc;
    std::vector<double> acc_err;
    std::vector<double> n_gen_phi;
    std::vector<double> n_rec_phi;
    double Ngen_cell = 0.0, Nrec_cell = 0.0;
};
using CellMap = std::map<std::tuple<int,int,int>, PhiArrays>;

static void write_period_json(
    const std::string& out_path,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const CellMap& cells)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "[acc] Cannot open " << out_path << "\n"; return; }
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": " << t_bins.size() << "},\n";
    ofs << "  \"bins\": {\n";
    bool first = true;
    for (const auto& kv : cells) {
        if (!first) ofs << ",\n"; first = false;
        int ix, iQ, it; std::tie(ix, iQ, it) = kv.first;
        const auto& pa = kv.second;
        ofs << "    \"(" << ix << "," << iQ << "," << it << ")\": {";
        ofs << "\"phi\":[";
        for (size_t i = 0; i < pa.phi_deg.size(); ++i) { if (i) ofs << ","; ofs << pa.phi_deg[i]; }
        ofs << "], \"acc\":[";
        for (size_t i = 0; i < pa.acc.size(); ++i) { if (i) ofs << ","; ofs << pa.acc[i]; }
        ofs << "], \"acc_err\":[";
        for (size_t i = 0; i < pa.acc_err.size(); ++i) { if (i) ofs << ","; ofs << pa.acc_err[i]; }
        ofs << "], \"counts_gen_phi\":[";
        for (size_t i = 0; i < pa.n_gen_phi.size(); ++i) { if (i) ofs << ","; ofs << pa.n_gen_phi[i]; }
        ofs << "], \"counts_rec_phi\":[";
        for (size_t i = 0; i < pa.n_rec_phi.size(); ++i) { if (i) ofs << ","; ofs << pa.n_rec_phi[i]; }
        ofs << "], \"total_gen\": " << pa.Ngen_cell << ", \"total_rec\": " << pa.Nrec_cell << "}";
    }
    ofs << "\n  }\n}\n";
}

// ---------------- plotting ----------------
static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list
) {
    std::set<std::pair<double,double>> qs, ts;
    for (const auto& b : scheme) {
        if (std::make_pair(b.xBmin, b.xBmax) == xBrange) {
            qs.emplace(b.Q2min, b.Q2max);
            ts.emplace(b.tmin,  b.tmax);
        }
    }
    Q2_list.assign(qs.begin(), qs.end());
    t_list.assign(ts.begin(), ts.end());
}

static void plot_cells_for_period(
    const std::string& period,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const CellMap& cells,
    const std::string& out_dir_plots)
{
    using std::filesystem::create_directories;
    std::error_code ec;
    create_directories(out_dir_plots, ec);

    static const auto PHI = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        const int W = 280 * ncols + 160;
        const int H = 240 * nrows + 170;

        std::ostringstream cname; cname << "c_acc_" << period << "_xB" << ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // Title
        pTop->cd();
        TLatex head;
        head.SetNDC(); head.SetTextAlign(22);
        head.SetTextFont(42);
        head.SetTextSize(0.36);
        std::ostringstream tit;
        tit << Form("Acceptance  %s   x_{B} #in [%.2g, %.2g]",
                    period.c_str(), xb.first, xb.second);
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        // Panels
        for (int r = 0; r < nrows; ++r) {
            const int it_global = findIndex(t_slice[r], t_bins);
            if (it_global < 0) continue;

            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int iQ_global = findIndex(Q2_slice[ccol], Q2_bins);
                if (iQ_global < 0) continue;

                pGrid->cd(r * ncols + ccol + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                // log-scale y for acceptance
                gPad->SetLogy(1);
                TH1* frame = gPad->DrawFrame(0.0, 1e-4, 360.0, 1.0);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetLabelSize(0.0001);
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("Acceptance");
                ax->CenterTitle(); ay->CenterTitle();
                ax->SetNdivisions(505);
                ax->SetTitleSize(0.060);
                ay->SetTitleSize(0.060);
                ay->SetLabelSize(0.048);
                ax->SetTitleOffset(1.25);
                ay->SetTitleOffset(1.35);

                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                auto itCell = cells.find(std::make_tuple(ix, iQ_global, it_global));
                if (itCell == cells.end()) continue;
                const auto& pa = itCell->second;

                std::vector<double> x, y, ey;
                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    double A = pa.acc[ip];
                    double sA = pa.acc_err[ip];
                    x.push_back(PHI[ip]);
                    y.push_back(A);
                    ey.push_back(std::max(1e-6, sA));
                }
                TGraphErrors* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, ey.data());
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(1.0);
                gr->SetLineWidth(2);
                gr->Draw("P SAME");

                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[ccol].first, Q2_slice[ccol].second,
                         t_slice[r].first,    t_slice[r].second));
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/plot_acceptance_" << period << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ---------------- key resolution ----------------
template <typename MapT>
static TTree* try_keys(const MapT& m, const std::vector<std::string>& keys) {
    for (const auto& k : keys) {
        auto it = m.find(k);
        if (it != m.end()) return it->second;
    }
    return nullptr;
}

} // anon

// =====================================================================
// Public driver
// =====================================================================
void compute_and_plot_acceptance(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& cuts_json_path,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning, 'x');
    const auto Q2_bins = uniqueRanges(binning, 'Q');
    const auto t_bins  = uniqueRanges(binning, 't');

    // Load MC-side 3σ windows from combined_cuts.json for all requested periods
    ComboMC comboCuts;
    const bool haveCuts = load_combined_mc_cuts(cuts_json_path, periods, comboCuts);
    if (!haveCuts) {
        std::ostringstream msg;
        msg << "[acc][FATAL] No MC cut blocks found in combined cuts JSON for requested periods.\n"
            << "  JSON path: " << cuts_json_path << "\n"
            << "  Example expected keys include e.g. \"DVCS_Sp18_inb_FD_FD\" or \"DVCS_Sp18_Inb_FD_FD\".\n"
            << "  Please verify period naming (inb vs Inb, out vs Out) and re-run exclusivity_cuts.\n";
        throw std::runtime_error(msg.str());
    }

    // Output dirs
    const fs::path json_dir  = fs::path(out_root_dir) / "jsons";
    const fs::path plot_root = fs::path(out_root_dir) / "acceptance";
    std::error_code ec;
    fs::create_directories(json_dir, ec);

    // Global DVCS MC thresholds (same as exclusivity_cuts)
    const MCCutsFixed globalCuts;

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);

        // Resolve using the exact period-based keys first (matches load_trees.cpp)
        const std::vector<std::string> gen_keys = {
            period + "_gen",   // e.g., "DVCS_Sp18_inb_gen"
            runTag + "_gen"    // fallback if ever used
        };
        const std::vector<std::string> rec_keys = {
            period + "_rec",       // e.g., "DVCS_Sp18_inb_rec"
            period + "_rec_mc",    // if alt naming exists
            runTag + "_rec",       // fallback
            runTag + "_rec_mc"     // fallback
        };

        TTree* tGen = try_keys(genMcTrees, gen_keys);
        TTree* tRec = try_keys(recMcTrees, rec_keys);

        if (!tGen || !tRec) {
            std::ostringstream tried;
            tried << "[acc][WARN] Missing MC trees for " << period << " (" << runTag << "). Tried gen keys: ";
            for (size_t i = 0; i < gen_keys.size(); ++i) { if (i) tried << ", "; tried << "\"" << gen_keys[i] << "\""; }
            tried << " | rec keys: ";
            for (size_t i = 0; i < rec_keys.size(); ++i) { if (i) tried << ", "; tried << "\"" << rec_keys[i] << "\""; }
            tried << " — skipping.\n";
            std::cerr << tried.str();
            continue;
        }

        // Build topology->CutMap for this specific period from the combined set.
        // We accept either lower-case or title-case period suffix in the JSON.
        std::map<std::string, CutMap> topoCutsForPeriod;
        for (const std::string topoKey : {"FD_FD","CD_FD","CD_FT"}) {
            // Try title-case period key first (matches your JSON), then lower-case.
            const std::string comboKeyTitle = period + "_" + topoKey; // as-passed (could already be *_Inb/*_Out)
            const std::string comboKeyLower = [&](){
                // if period ended with _inb/_out, also try title-case variant
                auto pv = period_variants(period);
                // prefer title-case variant if present
                for (const auto& ptry : pv) {
                    auto it = comboCuts.find(ptry + "_" + topoKey);
                    if (it != comboCuts.end()) return ptry + "_" + topoKey;
                }
                return period + "_" + topoKey;
            }();

            auto it1 = comboCuts.find(comboKeyLower);
            if (it1 != comboCuts.end()) {
                topoCutsForPeriod[topoKey] = it1->second;
                continue;
            }
            auto it2 = comboCuts.find(comboKeyTitle);
            if (it2 != comboCuts.end()) {
                topoCutsForPeriod[topoKey] = it2->second;
            }
        }

        // Accumulate
        AccMap acc;
        accumulate_generated(tGen, xB_bins, Q2_bins, t_bins, acc);
        accumulate_reconstructed(tRec, xB_bins, Q2_bins, t_bins, topologies, globalCuts, topoCutsForPeriod, acc);

        // Build per-cell phi arrays
        const auto PHI_DEG = phiCentersDeg();
        CellMap cells;
        for (int ix = 0; ix < (int)xB_bins.size(); ++ix)
        for (int iQ = 0; iQ < (int)Q2_bins.size(); ++iQ)
        for (int it = 0; it < (int)t_bins.size();  ++it) {
            PhiArrays pa;
            pa.phi_deg = PHI_DEG;
            pa.acc.resize(N_PHI_BINS, 0.0);
            pa.acc_err.resize(N_PHI_BINS, 0.0);
            pa.n_gen_phi.resize(N_PHI_BINS, 0.0);
            pa.n_rec_phi.resize(N_PHI_BINS, 0.0);
            pa.Ngen_cell = 0.0; pa.Nrec_cell = 0.0;

            for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                auto it4 = acc.find(std::make_tuple(ix, iQ, it, ip));
                double g = 0.0, r = 0.0;
                if (it4 != acc.end()) { g = it4->second.gen; r = it4->second.rec; }
                pa.n_gen_phi[ip] = g;
                pa.n_rec_phi[ip] = r;
                pa.Ngen_cell += g; pa.Nrec_cell += r;

                double A = 0.0, sA = 0.0;
                if (g > 0.0) {
                    A  = r / g;
                    sA = std::sqrt(std::max(0.0, A * (1.0 - A) / g)); // binomial
                }
                if (!std::isfinite(A)) A = 0.0;
                A = std::clamp(A, 0.0, 1.2);
                pa.acc[ip]     = A;
                pa.acc_err[ip] = sA;
            }
            cells[std::make_tuple(ix, iQ, it)] = std::move(pa);
        }

        // Write JSON
        const fs::path outJ = fs::path(out_root_dir) / "jsons" / ("acceptance_" + period + ".json");
        write_period_json(outJ.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, cells);
        std::cout << "[acc] Wrote acceptance JSON: " << outJ.string() << "\n";

        // Plots
        const fs::path outPlots = fs::path(out_root_dir) / "acceptance" / periodToRunTagKey(period);
        std::error_code ec2; fs::create_directories(outPlots, ec2);
        plot_cells_for_period(period, binning, xB_bins, Q2_bins, t_bins, cells, outPlots.string());
    }

    std::cout << "[acc] Acceptance computation complete.\n";
}