// total_counts.cpp
// ------------------------------------------------------------
// Compute helicity-separated total counts after exclusivity cuts.
// Strict mode: any missing lookup (trees, branches, dirs, files) causes a fatal exit.
// - Reads: DVCS data trees per period (must have: helicity, x, Q2, t1, phi2)
// - Writes:
//     • <out_root_dir>/jsons/total_counts.json             (master, nested "groups")
//     • <out_root_dir>/jsons/total_counts_<group>.json     (flat per-group file)
// - Produces per-group φ-binned plots under <out_root_dir>/total_counts_plots/<group>/
//
// Per-group JSON:
// {
//   "binning_meta": { "phi_bins": N, "xB_bins": nx, "Q2_bins": nQ, "t_bins": nt },
//   "bins": {
//      "(ix,iQ2,it,ip)": { "helicity": { "+1": Np, "-1": Nm }, "total": Np+Nm }, ...
//   }
// }
//
// Master JSON:
// {
//   "binning_meta": {...},
//   "groups": { "<group>": { "bins": {...} }, ... }
// }
// ------------------------------------------------------------

#include "total_counts.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>

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
#include <string>
#include <tuple>
#include <vector>

// ------------------------------------------------------------
// Helpers and constants
// ------------------------------------------------------------
namespace {

static constexpr int N_PHI_BINS = 12;
static constexpr double TWO_PI  = 2.0 * M_PI;

[[noreturn]] static void fatal(const std::string& msg) {
    std::cerr << "[total_counts][FATAL] " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

static inline double wrapToTwoPi(double phi) {
    if (!std::isfinite(phi)) return std::numeric_limits<double>::quiet_NaN();
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}

static inline int phiToBin(double phi) {
    double w = wrapToTwoPi(phi);
    if (!std::isfinite(w)) return -1;
    const double width = TWO_PI / static_cast<double>(N_PHI_BINS);
    int idx = static_cast<int>(std::floor(w / width));
    if (idx < 0) idx = 0;
    if (idx >= N_PHI_BINS) idx = N_PHI_BINS - 1;
    return idx;
}

static std::vector<std::pair<double,double>>
uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return {s.begin(), s.end()};
}

static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < (int)ranges.size(); ++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

static inline std::string keyStr(int ix,int iQ,int it,int ip) {
    std::ostringstream os; os<<"("<<ix<<","<<iQ<<","<<it<<","<<ip<<")";
    return os.str();
}

struct HelCounts { long long plus=0, minus=0; };

// ------------------------------------------------------------
// JSON writers (strict)
// ------------------------------------------------------------
static void write_total_counts_group_json(
    const std::string& out_path,
    const std::map<std::tuple<int,int,int,int>, HelCounts>& table,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) fatal(std::string("Cannot open for write: ") + out_path);
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size() << "},\n";
    ofs << "  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs << ",\n"; first=false;
        const HelCounts& hc = kv.second;
        int ix,iq,itb,ip; std::tie(ix,iq,itb,ip) = kv.first;
        ofs << "    \"" << keyStr(ix,iq,itb,ip) << "\": {"
            << "\"helicity\":{\"+1\":" << hc.plus
            << ",\"-1\":" << hc.minus << "},"
            << "\"total\":" << (hc.plus + hc.minus)
            << "}";
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote " << out_path << "\n";
}

static void write_total_counts_master_json(
    const std::string& out_path,
    const std::map<std::string, std::map<std::tuple<int,int,int,int>, HelCounts>>& allGroups,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) fatal(std::string("Cannot open for write: ") + out_path);
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size() << "},\n";
    ofs << "  \"groups\": {\n";

    bool firstG = true;
    for (const auto& gkv : allGroups) {
        if (!firstG) ofs << ",\n";
        firstG = false;
        ofs << "    \"" << gkv.first << "\": { \"bins\": {\n";
        bool firstB = true;
        for (const auto& kv : gkv.second) {
            if (!firstB) ofs << ",\n"; firstB=false;
            const HelCounts& hc = kv.second;
            int ix,iq,itb,ip; std::tie(ix,iq,itb,ip)=kv.first;
            ofs << "      \"" << keyStr(ix,iq,itb,ip) << "\": {"
                << "\"helicity\":{\"+1\":" << hc.plus
                << ",\"-1\":" << hc.minus << "},"
                << "\"total\":" << (hc.plus + hc.minus)
                << "}";
        }
        ofs << "\n    }}";
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote master " << out_path << "\n";
}

// ------------------------------------------------------------
// Plotting (strict save)
// ------------------------------------------------------------
static std::vector<double> phiCentersDeg() {
    std::vector<double> v(N_PHI_BINS);
    const double step = 360.0 / static_cast<double>(N_PHI_BINS);
    for (int i = 0; i < N_PHI_BINS; ++i) v[i] = (i + 0.5) * step;
    return v;
}

// ------------------------------------------------------------
// Plotting helper
// ------------------------------------------------------------
static void plot_group_counts(
    const std::string& group,
    const std::map<std::tuple<int,int,int,int>, HelCounts>& table,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::string& out_dir)
{
    namespace fs = std::filesystem;

    static const std::vector<double> PHI = phiCentersDeg();
    std::vector<double> X(N_PHI_BINS), ex(N_PHI_BINS, 0.0);
    for (int i=0;i<N_PHI_BINS;++i) X[i] = PHI[i];

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "[total_counts][FATAL] Cannot create directory: " << out_dir
                  << " (" << ec.message() << ")\n";
        std::exit(EXIT_FAILURE);
    }

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        std::set<std::pair<double,double>> q2set, tset;
        for (const auto& b : binning_scheme)
            if (std::make_pair(b.xBmin,b.xBmax)==xB_bins[ix]) {
                q2set.emplace(b.Q2min,b.Q2max);
                tset.emplace(b.tmin,b.tmax);
            }

        std::vector<std::pair<double,double>> Q2s(q2set.begin(),q2set.end());
        std::vector<std::pair<double,double>> Ts(tset.begin(),tset.end());
        if (Q2s.empty() || Ts.empty()) continue;

        const int nrows = Q2s.size();
        const int ncols = Ts.size();
        const int W = 260 * ncols + 120;
        const int H = 220 * nrows + 140;

        TCanvas* c = new TCanvas(Form("c_counts_%s_xB%d",group.c_str(),ix), "", W,H);
        TPad* pTop = new TPad("pTop","pTop",0.0,0.94,1.0,1.0);
        pTop->SetFillStyle(0); pTop->Draw();
        TPad* pGrid= new TPad("pGrid","pGrid",0.0,0.0,1.0,0.94);
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols,nrows,0.0001,0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextSize(0.55);
        head.SetTextAlign(22);
        head.DrawLatex(0.5,0.55,Form("%s   x_{B} [%.3g,%.3g]",group.c_str(),xB_bins[ix].first,xB_bins[ix].second));

        for (int r=0;r<nrows;++r){
            for (int ccol=0;ccol<ncols;++ccol){
                pGrid->cd(r*ncols+ccol+1);
                gPad->SetGrid(1,1);
                TH1* frame = gPad->DrawFrame(0,0,360,1);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Counts");

                std::vector<double> Yp(N_PHI_BINS,0), Ym(N_PHI_BINS,0);
                for(int ip=0;ip<N_PHI_BINS;++ip){
                    auto it = table.find({ix,r,ccol,ip});
                    if(it == table.end()) continue;
                    Yp[ip] = it->second.plus;
                    Ym[ip] = it->second.minus;
                }

                TGraph* gp = new TGraph(N_PHI_BINS,X.data(),Yp.data());
                gp->SetMarkerStyle(24);
                gp->SetMarkerColor(kRed);
                gp->Draw("LP SAME");

                TGraph* gm = new TGraph(N_PHI_BINS,X.data(),Ym.data());
                gm->SetMarkerStyle(20);
                gm->SetMarkerColor(kBlue);
                gm->Draw("LP SAME");
            }
        }

        // save and verify
        const std::string fpath = out_dir + "/plot_total_counts_" + group + "_xB_" + std::to_string(ix) + ".png";
        c->SaveAs(fpath.c_str());

        // fail-fast verification
        std::error_code fec;
        if (!fs::exists(fpath, fec) || (fs::file_size(fpath, fec) == 0)) {
            std::ostringstream em;
            em << "[total_counts][FATAL] Failed to save plot: " << fpath
               << " (exists=" << fs::exists(fpath)
               << ", size=" << fs::file_size(fpath, fec)
               << ", ec=" << (fec ? fec.message() : "ok") << ")";
            std::cerr << em.str() << std::endl;
            std::exit(EXIT_FAILURE);
        }

        delete c;
    }
}

} // namespace

// ------------------------------------------------------------
// Main compute function (STRICT)
// ------------------------------------------------------------
void compute_total_counts(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& /*topologies*/,     // not used here
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dataTrees,     // keys should match period strings
    const std::string& /*combined_cuts_json*/,          // not applied here
    const std::string& out_json_path,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    // ---- Validate binning axes ----
    const auto xB_bins = uniqueRanges(binning_scheme,'x');
    const auto Q2_bins = uniqueRanges(binning_scheme,'Q');
    const auto t_bins  = uniqueRanges(binning_scheme,'t');
    if (N_PHI_BINS <= 0) fatal("N_PHI_BINS must be > 0");
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty())
        fatal("Empty binning axes (xB/Q2/t). Check binning_scheme.");

    // ---- Prepare output dirs ----
    std::error_code ec;
    const fs::path plots_root = fs::path(out_root_dir) / "total_counts_plots";
    const fs::path jsons_dir  = fs::path(out_root_dir) / "jsons";
    if (!fs::create_directories(plots_root, ec) && ec)
        fatal(std::string("Cannot create plots root: ") + plots_root.string() + " (" + ec.message() + ")");
    ec.clear();
    if (!fs::create_directories(jsons_dir, ec) && ec)
        fatal(std::string("Cannot create jsons dir: ") + jsons_dir.string() + " (" + ec.message() + ")");

    // ---- Per-period processing ----
    std::map<std::string, std::map<std::tuple<int,int,int,int>, HelCounts>> allGroups;

    for (const auto& period : periods) {
        // Tree lookup must succeed (exact match first, then lowercase as a convenience)
        TTree* t = nullptr;
        auto it_exact = dataTrees.find(period);
        auto it_lower = dataTrees.find(toLower(period));
        if (it_exact != dataTrees.end() && it_exact->second) {
            t = it_exact->second;
            std::cout << "[total_counts][INFO] Using tree key \"" << period << "\"\n";
        } else if (it_lower != dataTrees.end() && it_lower->second) {
            t = it_lower->second;
            std::cout << "[total_counts][INFO] Using tree key \"" << toLower(period) << "\" for period \"" << period << "\"\n";
        } else {
            fatal(std::string("Missing DVCS data tree for period: ") + period +
                  " (looked for keys \"" + period + "\" and \"" + toLower(period) + "\")");
        }

        // Required branches (phi2 required; no Delta_phi fallback in this program)
        int    helicity = 0;
        double x = 0.0, Q2 = 0.0, t1 = 0.0, phi2 = std::numeric_limits<double>::quiet_NaN();

        auto mustBindI = [&](const char* name, int* addr){
            if (!t->GetBranch(name)) fatal(std::string("Missing integer branch '") + name + "' in period " + period);
            t->SetBranchAddress(name, addr);
        };
        auto mustBindD = [&](const char* name, double* addr){
            if (!t->GetBranch(name)) fatal(std::string("Missing double branch '") + name + "' in period " + period);
            t->SetBranchAddress(name, addr);
        };

        mustBindI("helicity", &helicity);
        mustBindD("x",       &x);
        mustBindD("Q2",      &Q2);
        mustBindD("t1",      &t1);
        if (!t->GetBranch("phi2"))
            fatal(std::string("Missing required branch 'phi2' in period ") + period + " (this program does not use Delta_phi)");
        t->SetBranchAddress("phi2", &phi2);

        std::map<std::tuple<int,int,int,int>, HelCounts> table;

        const Long64_t nent = t->GetEntries();
        if (nent <= 0) {
            std::cerr << "[total_counts][WARN] Tree has zero entries for period " << period << "\n";
        }

        Long64_t filled_events = 0;
        for (Long64_t i = 0; i < nent; ++i) {
            t->GetEntry(i);
            if (helicity!=+1 && helicity!=-1) continue;
            if (!std::isfinite(x) || !std::isfinite(Q2) || !std::isfinite(t1) || !std::isfinite(phi2)) continue;

            const int ix  = findBin(x,  xB_bins);
            const int iQ  = findBin(Q2, Q2_bins);
            const int itb = findBin(std::fabs(t1), t_bins);
            const int ip  = phiToBin(phi2);

            if (ix<0 || iQ<0 || itb<0 || ip<0) continue;

            auto& hc = table[{ix,iQ,itb,ip}];
            if (helicity==+1) hc.plus++; else hc.minus++;
            ++filled_events;
        }

        std::cout << "[total_counts][INFO] " << period
                  << ": entries=" << nent
                  << ", filled_bins=" << table.size()
                  << ", accepted_events=" << filled_events << "\n";

        allGroups[period] = std::move(table);

        // Plot per group (strict save inside)
        const std::string plot_dir = (plots_root / period).string();
        plot_group_counts(period, allGroups[period], binning_scheme, xB_bins, Q2_bins, t_bins, plot_dir);
    }

    // ---- Write master JSON ----
    write_total_counts_master_json(out_json_path, allGroups, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

    // ---- Write one flat JSON per group ----
    for (const auto& gkv : allGroups) {
        const std::string fname = (jsons_dir / ("total_counts_" + gkv.first + ".json")).string();
        write_total_counts_group_json(fname, gkv.second, N_PHI_BINS, xB_bins, Q2_bins, t_bins);
    }
}