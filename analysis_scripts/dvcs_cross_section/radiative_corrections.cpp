// radiative_corrections.cpp
#include "radiative_corrections.h"

#include "load_binning_scheme.h"

#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TPad.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

// ========== small style bootstrap ==========
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

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

// --------- helpers for bin structures ----------
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
    for (int i=0;i<(int)ranges.size();++i) if (ranges[i]==range) return i;
    return -1;
}

// Pretty period keys for JSON naming
static inline std::string periodToRunTagKey(const std::string& period) {
    // "DVCS_Fa18_inb" -> "fa18_inb"
    auto pos = period.find('_');
    std::string tail = (pos == std::string::npos || pos+1>=period.size())
        ? period : period.substr(pos+1);
    std::transform(tail.begin(), tail.end(), tail.begin(), ::tolower);
    return tail;
}
static inline std::string prettyPeriodKey(const std::string& runTagLower) {
    // "sp18_inb" -> "DVCS_Sp18_inb"
    std::string p = "DVCS_";
    bool upNext = true;
    for (char c : runTagLower) {
        if (upNext) { p.push_back(std::toupper(c)); upNext=false; }
        else if (c=='_') { p.push_back('_'); upNext=true; }
        else p.push_back(c);
    }
    return p;
}

// φ centers in degrees (12 bins)
static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
}
static inline std::vector<double> phiCentersRad() {
    std::vector<double> v(N_PHI_BINS);
    const double step = TWO_PI / double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) v[i] = (i+0.5)*step;
    return v;
}
static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize){
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

// --------- branch binder for GENERATED trees ----------
struct GenBranch {
    double x=0.0, Q2=0.0, t1=0.0, phi2=0.0;
    bool has_x=false, has_Q2=false, has_t1=false, has_phi=false;
    void bind(TTree* t){
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindD("x", &x, has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("t1",&t1,has_t1);
        bindD("phi2",&phi2,has_phi);
    }
};

// --------- accumulation structs ----------
struct CellPhiKey { int ix=0, iQ=0, it=0, ip=0; };
struct CellKey    { int ix=0, iQ=0, it=0; };

struct CountsStore {
    // per-φ counts
    std::map<std::tuple<int,int,int,int>, double> born_phi;
    std::map<std::tuple<int,int,int,int>, double> rad_phi;
    // per-cell totals (across φ)
    std::map<std::tuple<int,int,int>, double> born_tot;
    std::map<std::tuple<int,int,int>, double> rad_tot;
};

static inline int phiBinIndex(double phi){
    double w = std::fmod(phi, TWO_PI); if (w<0) w+=TWO_PI;
    const double width = TWO_PI/double(N_PHI_BINS);
    int ip = int(std::floor(w/width));
    if (ip<0) ip = 0; if (ip>=N_PHI_BINS) ip = N_PHI_BINS-1;
    return ip;
}

static inline int findBin1D(double v, const std::vector<std::pair<double,double>>& ranges){
    for (int k=0;k<(int)ranges.size();++k) if (v>=ranges[k].first && v<ranges[k].second) return k;
    return -1;
}

static void accumulate_generated(
    TTree* t,
    bool isBorn,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    CountsStore& acc)
{
    if (!t) return;
    GenBranch b; b.bind(t);
    const Long64_t n = t->GetEntries();
    for (Long64_t i=0;i<n;++i){
        t->GetEntry(i);
        if (!(b.has_x && b.has_Q2 && b.has_t1 && b.has_phi)) continue;
        const double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi2;

        int ix = findBin1D(xB, xB_bins);
        int iQ = findBin1D(Q2, Q2_bins);
        int it = findBin1D(tt,  t_bins);
        if (ix<0||iQ<0||it<0) continue;

        int ip = phiBinIndex(phi);

        auto key4 = std::make_tuple(ix,iQ,it,ip);
        auto key3 = std::make_tuple(ix,iQ,it);

        if (isBorn){
            acc.born_phi[key4] += 1.0;
            acc.born_tot[key3] += 1.0;
        } else {
            acc.rad_phi[key4]  += 1.0;
            acc.rad_tot[key3]  += 1.0;
        }
    }
}

// --------- RC computation per group ----------
struct PhiArrays {
    std::vector<double> phi_deg;
    std::vector<double> rc;      // Born/Rad
    std::vector<double> rc_err;
    std::vector<double> a_born;  // per-φ counts (Born)
    std::vector<double> b_rad;   // per-φ counts (Rad)
    double A_born = 0.0;         // totals
    double B_rad  = 0.0;
};

using RCPerCell = std::map<std::tuple<int,int,int>, PhiArrays>;

static RCPerCell compute_rc_per_cell(
    const CountsStore& acc,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    RCPerCell out;
    const auto PHI_DEG = phiCentersDeg();

    for (int ix=0; ix<(int)xB_bins.size(); ++ix)
    for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
    for (int it=0; it<(int)t_bins.size();  ++it) {
        auto key3 = std::make_tuple(ix,iQ,it);
        double A = acc.born_tot.count(key3) ? acc.born_tot.at(key3) : 0.0;
        double B = acc.rad_tot.count(key3)  ? acc.rad_tot.at(key3)  : 0.0;

        PhiArrays pa;
        pa.phi_deg = PHI_DEG;
        pa.rc.resize(N_PHI_BINS, 1.0);
        pa.rc_err.resize(N_PHI_BINS, 0.0);
        pa.a_born.resize(N_PHI_BINS, 0.0);
        pa.b_rad.resize(N_PHI_BINS, 0.0);
        pa.A_born = A; pa.B_rad = B;

        for (int ip=0; ip<N_PHI_BINS; ++ip){
            auto key4 = std::make_tuple(ix,iQ,it,ip);
            double a = acc.born_phi.count(key4) ? acc.born_phi.at(key4) : 0.0;
            double b = acc.rad_phi.count(key4)  ? acc.rad_phi.at(key4)  : 0.0;

            pa.a_born[ip]=a; pa.b_rad[ip]=b;

            double RC = 1.0, sRC = 0.0;
            if (A>0.0 && B>0.0 && a>0.0 && b>0.0) {
                RC  = (a*B)/(b*A);
                // conservative error propagation for (a/A)/(b/B)
                sRC = RC * std::sqrt( (1.0/std::max(a,1.0)) + (1.0/std::max(A,1.0))
                                     + (1.0/std::max(b,1.0)) + (1.0/std::max(B,1.0)) );
            } else if (A>0.0 && B>0.0 && (a==0.0 || b==0.0)) {
                // if one φ bin is empty, leave RC=1 with a large error bar (optional)
                RC  = 1.0;
                sRC = 0.0;
            }
            // clip to plotting range
            if (!std::isfinite(RC)) RC = 1.0;
            RC = std::clamp(RC, 0.0, 2.0);

            pa.rc[ip]     = RC;
            pa.rc_err[ip] = sRC;
        }

        out[key3] = std::move(pa);
    }

    return out;
}

// --------- plotting for a beam-energy group ----------
static void plot_group_rc(
    const std::string& energy_label,   // "10.59", "10.60", "10.2"
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const RCPerCell& rcPerCell,
    const std::string& out_dir_plots)  // e.g. output/radiative_correction_plots/10.59
{
    using std::filesystem::create_directories;
    std::error_code ec;
    create_directories(out_dir_plots, ec);

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        // Only the Q² and |t| bins used in this xB slice
        std::set<std::pair<double,double>> qs, ts;
        for (const auto& b : binning_scheme) {
            if (std::make_pair(b.xBmin,b.xBmax) == xb) {
                qs.emplace(b.Q2min,b.Q2max);
                ts.emplace(b.tmin,b.tmax);
            }
        }
        std::vector<std::pair<double,double>> Q2_slice(qs.begin(), qs.end());
        std::vector<std::pair<double,double>> t_slice(ts.begin(), ts.end());
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        // Canvas
        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_rc_E"<<energy_label<<"_xB"<<ix;
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
        tit << Form("Radiative correction  E_{beam}=%s GeV   x_{B} #in [%.2g, %.2g]",
                    energy_label.c_str(), xb.first, xb.second);
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        // Panels
        for (int r = 0; r < nrows; ++r) {
            const int it_global = findIndex(t_slice[r], t_bins);
            if (it_global < 0) continue;

            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int iQ_global = findIndex(Q2_slice[ccol], Q2_bins);
                if (iQ_global < 0) continue;

                pGrid->cd(r*ncols + ccol + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                // frame: x 0..360, y 0..2 with labels
                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 2.0);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetLabelSize(0.0001); // hide default tick labels, overlay custom
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("Born/Rad");
                ax->CenterTitle(); ay->CenterTitle();
                ax->SetNdivisions(505);
                ax->SetTitleSize(0.060);
                ay->SetTitleSize(0.060);
                ay->SetLabelSize(0.048);
                ax->SetTitleOffset(1.25);
                ay->SetTitleOffset(1.35);

                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                auto itCell = rcPerCell.find(std::make_tuple(ix, iQ_global, it_global));
                if (itCell == rcPerCell.end()) continue;
                const auto& pa = itCell->second;

                // make graph with errors
                // convert vectors to C arrays
                int npt = (int)pa.phi_deg.size();
                std::vector<double> ex(npt, 0.0);
                TGraphErrors* gr = new TGraphErrors(npt,
                    const_cast<double*>(pa.phi_deg.data()),
                    const_cast<double*>(pa.rc.data()),
                    ex.data(),
                    const_cast<double*>(pa.rc_err.data()));
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(1.0);
                gr->SetLineWidth(1);
                gr->Draw("P SAME");

                // annotation
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
        fout << out_dir_plots << "/rc_E" << energy_label << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// --------- JSON writers (arrays per φ) ----------
static void write_period_json(
    const std::string& out_path,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const RCPerCell& rcPerCell)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[radcorr] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<N_PHI_BINS<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : rcPerCell){
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
        const auto& pa = kv.second;
        ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {";
        ofs<<"\"phi\":[";
        for (int i=0;i<(int)pa.phi_deg.size();++i){ if(i) ofs<<","; ofs<<pa.phi_deg[i]; }
        ofs<<"], \"rc\":[";
        for (int i=0;i<(int)pa.rc.size();++i){ if(i) ofs<<","; ofs<<pa.rc[i]; }
        ofs<<"], \"rc_err\":[";
        for (int i=0;i<(int)pa.rc_err.size();++i){ if(i) ofs<<","; ofs<<pa.rc_err[i]; }
        ofs<<"], \"counts_born_phi\":[";
        for (int i=0;i<(int)pa.a_born.size();++i){ if(i) ofs<<","; ofs<<pa.a_born[i]; }
        ofs<<"], \"counts_rad_phi\":[";
        for (int i=0;i<(int)pa.b_rad.size();++i){ if(i) ofs<<","; ofs<<pa.b_rad[i]; }
        ofs<<"], \"total_born\": "<<pa.A_born<<", \"total_rad\": "<<pa.B_rad<<"}";
    }
    ofs<<"\n  }\n}\n";
}

static void write_all_groups_json(
    const std::string& out_path,
    const std::map<std::string, RCPerCell>& groupResults) // key: "10.59","10.60","10.2"
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[radcorr] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n  \"groups\": {\n";
    bool firstG=true;
    for (const auto& g : groupResults){
        if (!firstG) ofs<<",\n"; firstG=false;
        ofs<<"    \""<<g.first<<"\": {\"bins\": {";
        bool firstB=true;
        for (const auto& kv : g.second){
            if (!firstB) ofs<<",";
            firstB=false;
            int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
            const auto& pa = kv.second;
            ofs<<"\"("<<ix<<","<<iQ<<","<<it<<")\": {";
            ofs<<"\"phi\":[";
            for (int i=0;i<(int)pa.phi_deg.size();++i){ if(i) ofs<<","; ofs<<pa.phi_deg[i]; }
            ofs<<"], \"rc\":[";
            for (int i=0;i<(int)pa.rc.size();++i){ if(i) ofs<<","; ofs<<pa.rc[i]; }
            ofs<<"], \"rc_err\":[";
            for (int i=0;i<(int)pa.rc_err.size();++i){ if(i) ofs<<","; ofs<<pa.rc_err[i]; }
            ofs<<"], \"counts_born_phi\":[";
            for (int i=0;i<(int)pa.a_born.size();++i){ if(i) ofs<<","; ofs<<pa.a_born[i]; }
            ofs<<"], \"counts_rad_phi\":[";
            for (int i=0;i<(int)pa.b_rad.size();++i){ if(i) ofs<<","; ofs<<pa.b_rad[i]; }
            ofs<<"], \"total_born\": "<<pa.A_born<<", \"total_rad\": "<<pa.B_rad<<"}";
        }
        ofs<<"}}";
    }
    ofs<<"\n  }\n}\n";
}

} // anon

// =====================================================================
// Public driver
// =====================================================================
void compute_radiative_corrections(
    const std::vector<std::string>& periods,                        // e.g. DVCS_Fa18_inb, ...
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& genMcTrees,                // keys: "<tag>_gen"
    const std::map<std::string, TTree*>& radGenMcTrees,             // keys: "<tag>_gen_rad"
    const std::string& out_root_dir)                                // "output"
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // Output dirs
    const fs::path json_dir  = fs::path(out_root_dir)/"jsons";
    const fs::path plot_root = fs::path(out_root_dir)/"radiative_correction_plots";
    std::error_code ec;
    fs::create_directories(json_dir, ec);
    fs::create_directories(plot_root, ec);
    fs::create_directories(plot_root / "10.59", ec);
    fs::create_directories(plot_root / "10.60", ec);
    fs::create_directories(plot_root / "10.2",  ec);

    auto getOrNull = [](const auto& m, const std::string& k)->TTree*{
        auto it=m.find(k); return (it!=m.end()? it->second : nullptr);
    };

    // ---------------- Group definitions (by beam energy) ----------------
    struct Group { std::string label; std::vector<std::string> tags; };
    Group g1059{"10.59", {"sp18_inb","sp18_out"}};
    Group g1060{"10.60", {"fa18_inb","fa18_out"}};
    Group g1020{"10.2",  {"sp19_inb"}};

    std::vector<Group> groups = {g1059, g1060, g1020};
    std::map<std::string, RCPerCell> groupResults; // label -> RCPerCell

    // --------------- Compute RC for each group (pooled statistics) ---------------
    for (const auto& G : groups) {
        CountsStore acc;

        // accumulate generated counts over all runTags in this group
        for (const auto& tag : G.tags) {
            TTree* tBorn = getOrNull(genMcTrees,    tag + "_gen");
            TTree* tRad  = getOrNull(radGenMcTrees, tag + "_gen_rad");
            if (!tBorn || !tRad) {
                std::cerr<<"[radcorr] WARNING: missing gen trees for "<<tag<<" in group "<<G.label<<"\n";
                continue;
            }
            accumulate_generated(tBorn, /*isBorn=*/true,  xB_bins, Q2_bins, t_bins, acc);
            accumulate_generated(tRad,  /*isBorn=*/false, xB_bins, Q2_bins, t_bins, acc);
        }

        // compute per-cell per-φ Born/Rad
        RCPerCell rc = compute_rc_per_cell(acc, xB_bins, Q2_bins, t_bins);
        groupResults[G.label] = rc;

        // plots per group (ONLY energies)
        plot_group_rc(G.label, binning_scheme, xB_bins, Q2_bins, t_bins, rc,
                      (plot_root / G.label).string());

        // write a per-group JSON
        const fs::path outG = json_dir / ("radiative_corrections_group_" + G.label + ".json");
        write_period_json(outG.string(), xB_bins, Q2_bins, t_bins, rc);
        std::cout << "[radcorr] Wrote group JSON: " << outG.string() << "\n";
    }

    // --------------- Also save per-period JSONs (duplicating from group) ---------------
    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);   // "fa18_inb", "fa18_inb_supp", ...
        std::string energyLabel;
        if (runTag.rfind("sp18_",0)==0) energyLabel="10.59";
        else if (runTag.rfind("fa18_",0)==0) energyLabel="10.60";
        else if (runTag.rfind("sp19_",0)==0) energyLabel="10.2";
        else energyLabel="10.60"; // default fallback

        // Use the group's RC (fa18_inb_supp duplicates fa18_inb/fa18_out pool)
        const auto itG = groupResults.find(energyLabel);
        if (itG == groupResults.end()) {
            std::cerr<<"[radcorr] WARNING: no group result found for period "<<period<<" (E="<<energyLabel<<")\n";
            continue;
        }

        const std::string pretty = prettyPeriodKey(runTag);
        const fs::path outP = fs::path(out_root_dir)/"jsons"/("radiative_corrections_"+pretty+".json");
        write_period_json(outP.string(), xB_bins, Q2_bins, t_bins, itG->second);
        std::cout<<"[radcorr] Wrote period JSON (from group "<<energyLabel<<"): "<<outP.string()<<"\n";
    }

    // Optional: write one “all groups” container
    write_all_groups_json((fs::path(out_root_dir)/"jsons"/"radiative_corrections_all_groups.json").string(),
                          groupResults);

    std::cout<<"[radcorr] Radiative corrections (Born/Rad, generated-only) complete.\n";
}