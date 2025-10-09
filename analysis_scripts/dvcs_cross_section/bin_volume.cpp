#include "bin_volume.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TAxis.h>

#include <algorithm>
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
#include <utility>
#include <vector>

namespace {

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

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

// ---------- small helpers ----------
static inline int phiBinIndex(double phi){
    double w = std::fmod(phi, TWO_PI); if (w < 0.0) w += TWO_PI;
    const double width = TWO_PI / double(N_PHI_BINS);
    int ip = int(std::floor(w/width));
    if (ip < 0) ip = 0;
    if (ip >= N_PHI_BINS) ip = N_PHI_BINS - 1;
    return ip;
}
static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
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

// Build unique (min,max) bins from the scheme for a coordinate (x, Q, t).
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
static inline int findBin1D(double v, const std::vector<std::pair<double,double>>& ranges){
    for (int k=0;k<(int)ranges.size();++k) if (v>=ranges[k].first && v<ranges[k].second) return k;
    return -1;
}

// ---------- gen branches ----------
struct GenBranch {
    double x=0.0, Q2=0.0, t1=0.0, phi2=0.0;
    bool has_x=false, has_Q2=false, has_t1=false, has_phi=false;
    void bind(TTree* t){
        auto bindD=[&](const char*n,double*a,bool&f){ if(t && t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindD("x", &x, has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("t1",&t1,has_t1);
        bindD("phi2",&phi2,has_phi);
    }
};

// ---------- accumulation ----------
struct VolBin { double gen=0.0; }; // generator count
using AccMap = std::map<std::tuple<int,int,int,int>, VolBin>; // (ix,iQ,it,ip)

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
    for (Long64_t i=0;i<n;++i){
        t->GetEntry(i);
        double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi2;
        int ix=findBin1D(xB,xB_bins), iQ=findBin1D(Q2,Q2_bins), it=findBin1D(tt,t_bins);
        if (ix<0||iQ<0||it<0) continue;
        int ip = phiBinIndex(phi);
        acc[std::make_tuple(ix,iQ,it,ip)].gen += 1.0;
    }
}

// ---------- JSON writer ----------
struct PhiArrays {
    std::vector<double> phi_deg;
    std::vector<double> vol;     // normalized per cell: sum_phi vol = 1
    std::vector<double> vol_err; // Poisson → σ(N)/ΣN; here σ = sqrt(N)/ΣN
    std::vector<double> n_gen_phi;
    double Ngen_cell=0.0;
};
using CellMap = std::map<std::tuple<int,int,int>, PhiArrays>;

static void write_energy_json(
    const std::string& out_path,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const CellMap& cells)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[binvol] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : cells){
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
        const auto& pa = kv.second;
        ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {";
        auto dumpA=[&](const char* name,const std::vector<double>& v){
            ofs<<"\""<<name<<"\":["; for (size_t i=0;i<v.size();++i){ if(i)ofs<<","; ofs<<v[i]; } ofs<<"],";
        };
        dumpA("phi", pa.phi_deg);
        dumpA("vol", pa.vol);
        dumpA("vol_err", pa.vol_err);
        dumpA("counts_gen_phi", pa.n_gen_phi);
        ofs<<"\"total_gen\": "<<pa.Ngen_cell<<"}";
    }
    ofs<<"\n  }\n}\n";
}

// ---------- slice helpers for plotting ----------
static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list
) {
    std::set<std::pair<double,double>> qs, ts;
    for (const auto& b : scheme) {
        if (std::make_pair(b.xBmin,b.xBmax) == xBrange) {
            qs.emplace(b.Q2min,b.Q2max);
            ts.emplace(b.tmin,b.tmax);
        }
    }
    Q2_list.assign(qs.begin(), qs.end());
    t_list.assign(ts.begin(), ts.end());
}

// ---------- plotting ----------
static void plot_cells_for_energy(
    const std::string& energy_tag,             // "10.59", "10.60", "10.2"
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

    const auto PHI_DEG = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_binvol_"<<energy_tag<<"_xB"<<ix;
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
        tit << Form("Bin Volume (gen)  %s GeV   x_{B} #in [%.2g, %.2g]",
                    energy_tag.c_str(), xb.first, xb.second);
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

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 0.0); // y autoscale
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetLabelSize(0.0001); // custom degree ticks
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("Normalized bin volume");
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

                // Build graph
                std::vector<double> x, y, ey;
                double ymax=0.0;
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    double v  = pa.vol[ip];
                    double ev = std::max(1e-12, pa.vol_err[ip]);
                    x.push_back(PHI_DEG[ip]);
                    y.push_back(v);
                    ey.push_back(ev);
                    ymax = std::max(ymax, v + ev);
                }
                if (ymax <= 0.0) ymax = 1.0;
                frame->GetYaxis()->SetRangeUser(0.0, ymax*1.20);

                TGraphErrors* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, ey.data());
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(1.0);
                gr->SetLineWidth(2);
                gr->Draw("P SAME");

                // annotate Q² and -t
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
        fout << out_dir_plots << "/plot_bin_volume_" << energy_tag << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ---------- driver ----------
} // anon

void compute_and_plot_bin_volume(
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // Energy → list of generator keys
    const std::map<std::string, std::vector<std::string>> energy_to_genkeys = {
        {"10.59", {"sp18_inb_gen", "sp18_out_gen"}},
        {"10.60", {"fa18_inb_gen", "fa18_out_gen"}},
        {"10.2",  {"sp19_inb_gen"}}
    };

    const fs::path json_dir  = fs::path(out_root_dir)/"jsons";
    const fs::path plot_root = fs::path(out_root_dir)/"bin_volume";
    std::error_code ec; fs::create_directories(json_dir, ec);

    auto getOrNull = [](const auto& m, const std::string& k)->TTree*{
        auto it=m.find(k); return (it!=m.end()? it->second : nullptr);
    };

    // Loop over energies
    for (const auto& eg : energy_to_genkeys) {
        const std::string energy = eg.first;

        // Merge-generator accumulation per energy
        AccMap acc;
        for (const auto& tag : eg.second) {
            TTree* tGen = getOrNull(genMcTrees, tag);
            if (!tGen) {
                std::cerr<<"[binvol][WARN] Missing gen tree '"<<tag<<"' for energy "<<energy<<"\n";
                continue;
            }
            accumulate_generated(tGen, xB_bins, Q2_bins, t_bins, acc);
        }

        // Build per-cell normalized φ distributions
        CellMap cells;
        const auto PHI_DEG = phiCentersDeg();

        for (int ix=0; ix<(int)xB_bins.size(); ++ix)
        for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
        for (int it=0; it<(int)t_bins.size();  ++it) {
            PhiArrays pa;
            pa.phi_deg     = PHI_DEG;
            pa.vol.assign(N_PHI_BINS, 0.0);
            pa.vol_err.assign(N_PHI_BINS, 0.0);
            pa.n_gen_phi.assign(N_PHI_BINS, 0.0);
            pa.Ngen_cell = 0.0;

            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                auto it4 = acc.find(std::make_tuple(ix,iQ,it,ip));
                double g = (it4!=acc.end()? it4->second.gen : 0.0);
                pa.n_gen_phi[ip] = g;
                pa.Ngen_cell += g;
            }

            // Normalize to get volume fractions per φ-bin
            if (pa.Ngen_cell > 0.0) {
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    double g = pa.n_gen_phi[ip];
                    double f = g / pa.Ngen_cell;
                    // uncertainty: σ(f) ≈ sqrt(g) / N   (Poisson for g; treating N as fixed large)
                    double ef = (g > 0.0 ? std::sqrt(g) / pa.Ngen_cell : 0.0);
                    pa.vol[ip]     = f;
                    pa.vol_err[ip] = ef;
                }
            }
            cells[std::make_tuple(ix,iQ,it)] = std::move(pa);
        }

        // JSON
        const fs::path outJ = fs::path(out_root_dir)/"jsons"/("bin_volume_"+energy+".json");
        write_energy_json(outJ.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, cells);
        std::cout<<"[binvol] Wrote JSON: "<<outJ.string()<<"\n";

        // Plots (per energy)
        const fs::path outPlots = fs::path(out_root_dir)/"bin_volume"/energy;
        std::error_code ec2; fs::create_directories(outPlots, ec2);
        plot_cells_for_energy(energy, binning_scheme, xB_bins, Q2_bins, t_bins, cells, outPlots.string());
    }

    std::cout<<"[binvol] Bin-volume computation complete.\n";
}