#include "uncorrected_cross_section.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TAxis.h>
#include <TGaxis.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>

#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {
constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

// ------------ style ------------
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetLineWidth(2);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        const int rf = 42;
        gStyle->SetTitleFont(rf, "XYZ");
        gStyle->SetLabelFont(rf, "XYZ");
        gStyle->SetTextFont(rf);
    }
} _style_guard;

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
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

// Build unique (min,max) bins from the scheme for a coordinate
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

// ---------- I/O helpers ----------
static json load_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[json] Failed to open " << path << "\n";
        return json();
    }
    json j; f >> j; return j;
}

// CSV luminosity: run,total,pos,neg,*,*
static double read_total_luminosity(const std::string& filepath) {
    std::ifstream in(filepath);
    if (!in.is_open()) {
        std::cerr << "[luminosity] Cannot open " << filepath << "\n";
        return 0.0;
    }
    auto trim = [](std::string s){
        size_t a = s.find_first_not_of(" \t\r\n");
        size_t b = s.find_last_not_of(" \t\r\n");
        if (a == std::string::npos) return std::string();
        return s.substr(a, b - a + 1);
    };

    double sum_pos=0.0, sum_neg=0.0;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0]=='#') continue;
        std::stringstream ss(line);
        std::string tok; std::vector<std::string> cols;
        while (std::getline(ss, tok, ',')) cols.push_back(trim(tok));
        if (cols.size() < 4) continue;
        try {
            double pos = std::stod(cols[2]);
            double neg = std::stod(cols[3]);
            if (pos > 0) sum_pos += pos;
            if (neg > 0) sum_neg += neg;
        } catch (...) {}
    }
    in.close();
    std::cout << "[luminosity] " << filepath
              << "  total_pos=" << sum_pos
              << "  total_neg=" << sum_neg << "\n";
    return sum_pos + sum_neg;
}

// ---------- slicers ----------
static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list
){
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

static inline std::string periodToRunTagKey(const std::string& period) {
    // "DVCS_Fa18_inb" -> "fa18_inb"
    auto pos = period.find('_');
    std::string tag = period;
    std::transform(tag.begin(), tag.end(), tag.begin(), ::tolower);
    if (pos == std::string::npos || pos+1 >= period.size()) return tag;
    return tag.substr(pos+1);
}

} // anon

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
void compute_uncorrected_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,  // e.g. "output/jsons"
    const std::string& unfolded_counts_dir,  // e.g. "output/jsons"
    const std::string& luminosity_dir,       // e.g. "imports/integrated_luminosity"
    const std::string& output_dir)           // e.g. "output/uncorrected_cross_section"
{
    // ----- bin edges & helpers -----
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    const auto PHI_DEG = phiCentersDeg();

    // ----- luminosity -----
    const double L_fa18_inb = read_total_luminosity(luminosity_dir + "/rga_fa18_inb.txt");
    const double L_fa18_out = read_total_luminosity(luminosity_dir + "/rga_fa18_out.txt");
    const double L_sp19_inb = read_total_luminosity(luminosity_dir + "/rga_sp19_inb.txt");

    const double L_10p6 = L_fa18_inb + L_fa18_out;
    const double L_10p2 = L_sp19_inb;

    const std::map<std::string, double> energyL = {
        {"10.60", L_10p6},
        {"10.2",  L_10p2}
    };

    // ----- energy → periods to combine -----
    const std::map<std::string, std::vector<std::string>> energyPeriods = {
        {"10.60", {"DVCS_Fa18_inb", "DVCS_Fa18_out"}},
        {"10.2",  {"DVCS_Sp19_inb"}}
        // (Sp18* can be added later)
    };

    // ----- load bin-volume JSONs -----
    std::map<std::string, json> jBinVol;
    for (const auto& kv : energyPeriods) {
        const std::string e = kv.first;
        const std::string path = (fs::path(bin_volume_json_dir) / ("bin_volume_" + e + ".json")).string();
        jBinVol[e] = load_json(path);
        if (jBinVol[e].empty()) {
            std::cerr << "[binvol] Missing " << path << "\n";
        } else {
            std::cout << "[binvol] Loaded " << path << "\n";
        }
    }

    // ----- output dirs -----
    fs::create_directories(output_dir);
    fs::create_directories(fs::path(output_dir)/"jsons");
    fs::create_directories(fs::path(output_dir)/"plots");

    // ============================================================
    // Loop energies
    // ============================================================
    for (const auto& eg : energyPeriods) {
        const std::string energy = eg.first;
        const double L_total = energyL.at(energy);
        const json& jvol = jBinVol[energy];

        if (jvol.empty() || L_total <= 0.0) {
            std::cerr << "[xsec] Skipping " << energy << " due to missing bin-volume or luminosity.\n";
            continue;
        }

        // ------------ sum unfolded yields (helicity-resolved) across periods ------------
        // data structure: (ix,iQ,it) -> arrays (len 12) for plus/minus + errors
        struct HelData { std::vector<double> y, ye; HelData(): y(N_PHI_BINS,0.0), ye(N_PHI_BINS,0.0){} };
        using CellKey = std::tuple<int,int,int>;
        std::map<CellKey, HelData> sumPlus, sumMinus;

        for (const std::string& period : eg.second) {
            const std::string path = (fs::path(unfolded_counts_dir) / ("unfolded_"+period+".json")).string();
            json ju = load_json(path);
            if (ju.empty()) {
                std::cerr << "[xsec] Missing unfolded " << path << " (skipping this period)\n";
                continue;
            }

            // iterate bins: key "(ix,iQ,it)"
            if (!ju.contains("bins")) continue;
            for (auto it = ju["bins"].begin(); it != ju["bins"].end(); ++it) {
                const std::string bkey = it.key();
                const json& cell = it.value();

                // parse tuple key
                int ix=0,iQ=0,itb=0;
                if (std::sscanf(bkey.c_str(),"(%d,%d,%d)",&ix,&iQ,&itb) != 3) continue;

                // grab arrays
                if (!cell.contains("yield_plus") || !cell.contains("yield_minus")) continue;
                const auto& yp = cell["yield_plus"];
                const auto& ym = cell["yield_minus"];
                const auto& ep = cell.contains("yield_plus_err")  ? cell["yield_plus_err"]  : json::array();
                const auto& em = cell.contains("yield_minus_err") ? cell["yield_minus_err"] : json::array();
                if (yp.size() != N_PHI_BINS || ym.size() != N_PHI_BINS) continue;

                HelData& hp = sumPlus[CellKey(ix,iQ,itb)];
                HelData& hm = sumMinus[CellKey(ix,iQ,itb)];
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    double vp = yp[ip].get<double>();
                    double vm = ym[ip].get<double>();
                    double epv = (ep.size()==N_PHI_BINS ? std::max(0.0, ep[ip].get<double>()) : std::sqrt(std::max(0.0, vp)));
                    double emv = (em.size()==N_PHI_BINS ? std::max(0.0, em[ip].get<double>()) : std::sqrt(std::max(0.0, vm)));

                    // sum values; add variances
                    hp.y[ip]  += vp;
                    hm.y[ip]  += vm;
                    hp.ye[ip]  = std::sqrt(hp.ye[ip]*hp.ye[ip] + epv*epv);
                    hm.ye[ip]  = std::sqrt(hm.ye[ip]*hm.ye[ip] + emv*emv);
                }
            }
        }

        // ------------ compute cross sections (use absolute per-φ volumes) ------------
        json jxsec;
        jxsec["energy"] = energy;
        jxsec["luminosity_total"] = L_total;
        jxsec["bins"] = json::object();

        if (!jvol.contains("bins")) {
            std::cerr << "[binvol] Malformed bin-volume JSON for energy " << energy << " (no 'bins').\n";
        } else {
            for (auto it = jvol["bins"].begin(); it != jvol["bins"].end(); ++it) {
                const std::string bkey = it.key();
                const json& bv = it.value();

                if (!bv.contains("vol")) continue;
                const auto& volphi = bv["vol"];
                if (volphi.size() != N_PHI_BINS) continue;

                // locate summed yields
                int ix=0,iQ=0,itb=0;
                if (std::sscanf(bkey.c_str(),"(%d,%d,%d)",&ix,&iQ,&itb) != 3) continue;

                auto itP = sumPlus.find(std::make_tuple(ix,iQ,itb));
                auto itM = sumMinus.find(std::make_tuple(ix,iQ,itb));
                if (itP == sumPlus.end() && itM == sumMinus.end()) {
                    // nothing to plot/output for this bin
                    continue;
                }

                // per-phi denominator: L_total * V_phi (absolute volume)
                std::vector<double> denom(N_PHI_BINS, 0.0);
                double Vtot = 0.0; // keep as metadata
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const double Vphi = std::max(0.0, volphi[ip].get<double>());
                    Vtot += Vphi;
                    denom[ip] = L_total * Vphi;
                }

                auto mkOut = [&](const HelData* hd)->json {
                    json o;
                    std::vector<double> xv(N_PHI_BINS), yv(N_PHI_BINS,0.0), ev(N_PHI_BINS,0.0);
                    for (int ip=0; ip<N_PHI_BINS; ++ip) {
                        xv[ip] = PHI_DEG[ip];
                        const double num = (hd ? hd->y[ip]  : 0.0);
                        const double en  = (hd ? hd->ye[ip] : 0.0);
                        const double d   = denom[ip];
                        yv[ip] = (d > 0.0 ? num / d : 0.0);
                        ev[ip] = (d > 0.0 ? en  / d : 0.0);
                    }
                    o["phi"]      = xv;
                    o["xsec"]     = yv;
                    o["xsec_err"] = ev;
                    return o;
                };

                jxsec["bins"][bkey] = {
                    {"helicity_plus",  mkOut(itP==sumPlus.end()? nullptr : &itP->second)},
                    {"helicity_minus", mkOut(itM==sumMinus.end()? nullptr : &itM->second)},
                    {"bin_volume_total", Vtot},
                    {"vol_phi", volphi} // keep what we used
                };
            }
        }

        // ------------ save JSON ------------
        const std::string out_json = (fs::path("output/jsons") / ("uncorrected_xsec_" + energy + ".json")).string();
        std::ofstream ofs(out_json);
        ofs << std::setw(2) << jxsec << "\n";
        ofs.close();
        std::cout << "[xsec] Wrote " << out_json << "\n";

        // ------------ plots (xB slices; overlay + and −) ------------
        const auto Q2_all = uniqueRanges(binning_scheme, 'Q');
        const auto t_all  = uniqueRanges(binning_scheme, 't');

        for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
            const auto xb = xB_bins[ix];
            std::vector<std::pair<double,double>> Q2_slice, t_slice;
            uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
            if (Q2_slice.empty() || t_slice.empty()) continue;

            const int nrows = (int)t_slice.size();
            const int ncols = (int)Q2_slice.size();
            const int W = 280*ncols + 160;
            const int H = 240*nrows + 170;

            std::ostringstream cname; cname<<"c_xsec_"<<energy<<"_xB"<<ix;
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
            tit << Form("Uncorrected Cross Section  %s GeV   x_{B} #in [%.2g, %.2g]",
                        energy.c_str(), xb.first, xb.second);
            head.DrawLatex(0.5, 0.55, tit.str().c_str());

            // Panels
            for (int r=0; r<nrows; ++r) {
                const int it_global = findIndex(t_slice[r], t_all);
                if (it_global < 0) continue;
                for (int cc=0; cc<ncols; ++cc) {
                    const int iQ_global = findIndex(Q2_slice[cc], Q2_all);
                    if (iQ_global < 0) continue;

                    pGrid->cd(r*ncols + cc + 1);
                    gPad->SetGrid(1,1);
                    gPad->SetTopMargin(0.08);
                    gPad->SetBottomMargin(0.18);
                    gPad->SetLeftMargin(0.15);
                    gPad->SetRightMargin(0.10);
                    gPad->SetLogy();

                    TH1* frame = gPad->DrawFrame(0.0, 1e-4, 360.0, 1e3);
                    frame->GetXaxis()->SetLabelSize(0.0001);
                    frame->GetXaxis()->SetTitle("#phi (deg)");
                    frame->GetYaxis()->SetTitle("d#sigma/d#phi (uncorr.)");
                    frame->GetXaxis()->CenterTitle();
                    frame->GetYaxis()->CenterTitle();
                    frame->GetXaxis()->SetNdivisions(505);
                    frame->GetXaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetLabelSize(0.048);
                    frame->GetXaxis()->SetTitleOffset(1.25);
                    frame->GetYaxis()->SetTitleOffset(1.35);

                    drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                    const std::string bkey = "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
                    if (!jxsec["bins"].contains(bkey)) continue;
                    const json& jb = jxsec["bins"][bkey];

                    auto drawHel = [&](const char* node, int mstyle, int color){
                        if (!jb.contains(node)) return (TGraphErrors*)nullptr;
                        const auto& h = jb[node];
                        if (!h.contains("phi") || !h.contains("xsec") || !h.contains("xsec_err")) return (TGraphErrors*)nullptr;
                        const auto& xp = h["phi"];
                        const auto& yp = h["xsec"];
                        const auto& ep = h["xsec_err"];
                        if (xp.size()!=N_PHI_BINS || yp.size()!=N_PHI_BINS || ep.size()!=N_PHI_BINS) return (TGraphErrors*)nullptr;
                        std::vector<double> x(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
                        double ymax=0.0;
                        for (int i=0;i<N_PHI_BINS;++i){
                            x[i]=xp[i].get<double>(); y[i]=yp[i].get<double>(); e[i]=std::max(1e-12, ep[i].get<double>());
                            ymax = std::max(ymax, y[i]+e[i]);
                        }
                        // if (ymax > 0.0) frame->GetYaxis()->SetRangeUser(1e-4, std::max(1.0, ymax*1.5));
                        frame->GetYaxis()->SetRangeUser(1e-4, 1e3);
                        auto* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, e.data());
                        gr->SetMarkerStyle(mstyle);
                        gr->SetMarkerSize(1.0);
                        gr->SetLineWidth(2);
                        gr->SetLineColor(color);
                        gr->SetMarkerColor(color);
                        gr->Draw("P SAME");
                        return gr;
                    };

                    TGraphErrors* gp = drawHel("helicity_plus", 20, kBlue+1);
                    TGraphErrors* gm = drawHel("helicity_minus", 25, kRed+1);

                    TLatex lab;
                    lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                    lab.SetTextFont(42);
                    lab.DrawLatex(0.15, 0.94,
                        Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                             Q2_slice[cc].first, Q2_slice[cc].second,
                             t_slice[r].first,   t_slice[r].second));

                    if (gp || gm){
                        TLegend* leg = new TLegend(0.50, 0.72, 0.90, 0.92);
                        leg->SetBorderSize(1);
                        leg->SetLineColor(kBlack);
                        leg->SetFillColor(kWhite);
                        leg->SetFillStyle(1001);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.040);
                        if (gp) leg->AddEntry(gp, "+ helicity", "lep");
                        if (gm) leg->AddEntry(gm, "- helicity", "lep");
                        leg->Draw();
                    }
                }
            }

            const std::string outP = (fs::path(output_dir)/"plots"/("uncorrected_xsec_"+energy+"_xB_"+std::to_string(ix)+".png")).string();
            c->SaveAs(outP.c_str());
            delete c;
        }
    }

    std::cout << "[xsec] Finished uncorrected cross-section generation.\n";
}


// Add this to uncorrected_cross_section.cpp

void compare_unpolarized_xsec_sp18out_vs_fa18out(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,   // e.g. "output/jsons"
    const std::string& unfolded_counts_dir,   // e.g. "output/jsons"
    const std::string& luminosity_dir,        // e.g. "imports/integrated_luminosity"
    const std::string& output_dir             // e.g. "output/uncorrected_cross_section/compare_sp18out_vs_fa18out"
) {
    namespace fs = std::filesystem;

    // ----- helpers from this file are assumed available:
    // uniqueRanges(...), uniqueQT_for_xB(...), phiCentersDeg(), load_json(...)

    auto read_total_second_col = [&](const std::string& filepath)->double {
        // CSV format: run,total,pos,neg,*,*
        std::ifstream in(filepath);
        if (!in.is_open()) {
            std::cerr << "[luminosity] Cannot open " << filepath << "\n";
            return 0.0;
        }
        auto trim = [](std::string s){
            size_t a = s.find_first_not_of(" \t\r\n");
            size_t b = s.find_last_not_of(" \t\r\n");
            if (a == std::string::npos) return std::string();
            return s.substr(a, b - a + 1);
        };
        double total = 0.0;
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0]=='#') continue;
            std::stringstream ss(line);
            std::string tok; std::vector<std::string> cols;
            while (std::getline(ss, tok, ',')) cols.push_back(trim(tok));
            if (cols.size() < 2) continue;
            try {
                double t = std::stod(cols[1]); // second column
                if (t > 0.0) total += t;
            } catch (...) { /* ignore */ }
        }
        in.close();
        std::cout << "[luminosity] " << filepath << "  total(second_col)=" << total << "\n";
        return total;
    };

    // ----- bin edges -----
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    const auto PHI_DEG = phiCentersDeg();

    // ----- shared bin-volume JSON (10.60 GeV umbrella)
    const std::string binvol_path = (fs::path(bin_volume_json_dir) / "bin_volume_10.60.json").string();
    const json jvol = load_json(binvol_path);
    if (jvol.empty() || !jvol.contains("bins")) {
        std::cerr << "[compare] Missing or malformed bin-volume JSON: " << binvol_path << "\n";
        return;
    }

    // ----- luminosities: use SECOND column totals for both periods -----
    const double L_sp18_out = read_total_second_col(luminosity_dir + "/rga_sp18_out.txt");
    const double L_fa18_out = read_total_second_col(luminosity_dir + "/rga_fa18_out.txt");
    if (L_sp18_out <= 0.0 || L_fa18_out <= 0.0) {
        std::cerr << "[compare] Non-positive luminosities; aborting.\n";
        return;
    }

    // ----- load unfolded yield JSONs for the two periods we compare -----
    const std::string path_sp18 = (fs::path(unfolded_counts_dir) / "unfolded_DVCS_Sp18_out.json").string();
    const std::string path_fa18 = (fs::path(unfolded_counts_dir) / "unfolded_DVCS_Fa18_out.json").string();
    const json ju_sp18 = load_json(path_sp18);
    const json ju_fa18 = load_json(path_fa18);
    if (ju_sp18.empty() || ju_fa18.empty()) {
        std::cerr << "[compare] Missing unfolded JSONs: " << path_sp18 << " or " << path_fa18 << "\n";
        return;
    }

    // ----- output dirs -----
    fs::create_directories(output_dir);

    // ============================================================
    // Build a convenience accessor for cross sections in a bin
    // (unpolarized = (plus + minus) / (L * V_phi))
    // ============================================================
    auto get_xsec_for_bin = [&](const json& jvol_bins, const json& ju, int ix, int iQ, int itb,
                                double L, std::vector<double>& y, std::vector<double>& ye)->bool {
        y.assign(N_PHI_BINS, 0.0);
        ye.assign(N_PHI_BINS, 0.0);
        const std::string bkey = "(" + std::to_string(ix) + "," + std::to_string(iQ) + "," + std::to_string(itb) + ")";
        if (!jvol_bins.contains(bkey)) return false;
        const json& bv = jvol_bins[bkey];
        if (!bv.contains("vol") || bv["vol"].size() != N_PHI_BINS) return false;

        // Unfolded bin
        if (!ju.contains("bins") || !ju["bins"].contains(bkey)) return false;
        const json& cell = ju["bins"][bkey];
        if (!cell.contains("yield_plus") || !cell.contains("yield_minus")) return false;

        const auto& yp = cell["yield_plus"];
        const auto& ym = cell["yield_minus"];
        const auto& ep = cell.contains("yield_plus_err")  ? cell["yield_plus_err"]  : json::array();
        const auto& em = cell.contains("yield_minus_err") ? cell["yield_minus_err"] : json::array();

        for (int ip=0; ip<N_PHI_BINS; ++ip) {
            const double Vphi = std::max(0.0, bv["vol"][ip].get<double>());
            const double denom = L * Vphi;
            if (denom <= 0.0) { y[ip] = 0.0; ye[ip] = 0.0; continue; }

            const double vplus  = yp[ip].get<double>();
            const double vminus = ym[ip].get<double>();
            const double eplus  = (ep.size()==N_PHI_BINS ? std::max(0.0, ep[ip].get<double>()) : std::sqrt(std::max(0.0, vplus)));
            const double eminus = (em.size()==N_PHI_BINS ? std::max(0.0, em[ip].get<double>()) : std::sqrt(std::max(0.0, vminus)));

            const double vtot  = vplus + vminus;
            const double etot2 = eplus*eplus + eminus*eminus;

            y[ip]  = vtot / denom;
            ye[ip] = std::sqrt(std::max(0.0, etot2)) / denom;
        }
        return true;
    };

    // ============================================================
    // Plot per xB-slice, panelled in Q2 x t, with statistics-weighted
    // average percent difference annotation per subplot.
    // ============================================================
    const auto Q2_all = uniqueRanges(binning_scheme, 'Q');
    const auto t_all  = uniqueRanges(binning_scheme, 't');

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname << "c_xsecU_sp18out_vs_fa18out_xB" << ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // Title (smaller font; text is "#sigma_{U}")
        pTop->cd();
        TLatex head;
        head.SetNDC(); head.SetTextAlign(22);
        head.SetTextFont(42);
        head.SetTextSize(0.30); // smaller than before
        std::ostringstream tit;
        tit << Form("#sigma_{U}   x_{B} #in [%.2g, %.2g]   (Sp18 out, Fa18 out)", xb.first, xb.second);
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        for (int r=0; r<nrows; ++r) {
            const int it_global = [&](){
                for (int k=0;k<(int)t_all.size();++k) if (t_all[k]==t_slice[r]) return k;
                return -1;
            }();
            if (it_global < 0) continue;

            for (int cc=0; cc<ncols; ++cc) {
                const int iQ_global = [&](){
                    for (int k=0;k<(int)Q2_all.size();++k) if (Q2_all[k]==Q2_slice[cc]) return k;
                    return -1;
                }();
                if (iQ_global < 0) continue;

                pGrid->cd(r*ncols + cc + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.18);   // a bit more padding on the left
                gPad->SetRightMargin(0.10);
                gPad->SetLogy();

                // Frame: lower bound at 1e-3 per request
                TH1* frame = gPad->DrawFrame(0.0, 1e-3, 360.0, 1e3);
                frame->GetXaxis()->SetLabelSize(0.0001);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("d#sigma_{U}/d#phi (uncorr.)");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.048);
                frame->GetXaxis()->SetTitleOffset(1.25);
                frame->GetYaxis()->SetTitleOffset(1.40); // nudge to avoid clipping

                TGaxis* ax = new TGaxis(gPad->GetUxmin(), gPad->GetUymin(),
                                        gPad->GetUxmax(), gPad->GetUymin(),
                                        0.0, 360.0, 4, "");
                ax->SetLabelFont(42);
                ax->SetLabelSize(0.048);
                ax->SetLabelOffset(0.012);
                ax->SetTitle("");
                ax->SetTickSize(0.02);
                ax->Draw();

                // Build unpolarized cross sections for both periods for this cell
                std::vector<double> y_sp18, e_sp18, y_fa18, e_fa18;
                const bool ok_sp18 = get_xsec_for_bin(jvol["bins"], ju_sp18, ix, iQ_global, it_global,
                                                      L_sp18_out, y_sp18, e_sp18);
                const bool ok_fa18 = get_xsec_for_bin(jvol["bins"], ju_fa18, ix, iQ_global, it_global,
                                                      L_fa18_out, y_fa18, e_fa18);
                if (!(ok_sp18 || ok_fa18)) continue;

                // Graphs
                auto drawSet = [&](const std::vector<double>& y, const std::vector<double>& e,
                                   int mstyle, int color)->TGraphErrors* {
                    std::vector<double> x(N_PHI_BINS), yy(N_PHI_BINS), ee(N_PHI_BINS);
                    for (int i=0;i<N_PHI_BINS;++i) {
                        x[i]  = PHI_DEG[i];
                        yy[i] = (i<(int)y.size()? y[i] : 0.0);
                        ee[i] = (i<(int)e.size()? std::max(1e-12, e[i]) : 0.0);
                    }
                    auto* gr = new TGraphErrors(N_PHI_BINS, x.data(), yy.data(), nullptr, ee.data());
                    gr->SetMarkerStyle(mstyle);
                    gr->SetMarkerSize(1.0);
                    gr->SetLineWidth(2);
                    gr->SetLineColor(color);
                    gr->SetMarkerColor(color);
                    gr->Draw("P SAME");
                    return gr;
                };

                TGraphErrors* g_sp18 = ok_sp18 ? drawSet(y_sp18, e_sp18, 20, kBlue+1) : nullptr; // blue circles
                TGraphErrors* g_fa18 = ok_fa18 ? drawSet(y_fa18, e_fa18, 25, kRed+1)  : nullptr; // red squares

                // Panel subtitle: Q2 and -t ranges
                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[cc].first, Q2_slice[cc].second,
                         t_slice[r].first,   t_slice[r].second));

                // Average percent difference (statistics-weighted)
                // per phi bin: pd_i = |A - B| / ((A + B)/2) * 100
                // weight_i = 1 / ( (eA/A)^2 + (eB/B)^2 ), only if A>0, B>0
                double wsum = 0.0, pd_sum = 0.0;
                int n_used = 0;
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    if (!ok_sp18 || !ok_fa18) break;
                    const double A = y_sp18[ip];
                    const double B = y_fa18[ip];
                    const double eA = e_sp18[ip];
                    const double eB = e_fa18[ip];
                    if (!(A > 0.0) || !(B > 0.0)) continue;
                    if (!std::isfinite(A) || !std::isfinite(B)) continue;

                    const double pd = std::fabs(A - B) / std::max(1e-12, 0.5*(A + B)) * 100.0;

                    const double rA = eA / std::max(1e-12, A);
                    const double rB = eB / std::max(1e-12, B);
                    const double w  = 1.0 / std::max(1e-12, (rA*rA + rB*rB));

                    if (std::isfinite(pd) && std::isfinite(w) && w > 0.0) {
                        pd_sum += w * pd;
                        wsum   += w;
                        ++n_used;
                    }
                }
                if (wsum > 0.0 && n_used > 0) {
                    const double pd_avg = pd_sum / wsum;
                    TLatex stat;
                    stat.SetNDC(); stat.SetTextAlign(13); // left/top-ish
                    stat.SetTextFont(42);
                    stat.SetTextSize(0.034);
                    stat.DrawLatex(0.15, 0.86, Form("#LT#Delta%%#GT_{w} = %.2f%%  (N = %d)", pd_avg, n_used));
                }

                // Legend (top-right, pulled a bit left so it won’t clip)
                if (g_sp18 || g_fa18) {
                    TLegend* leg = new TLegend(0.48, 0.70, 0.93, 0.91);
                    leg->SetBorderSize(1);
                    leg->SetLineColor(kBlack);
                    leg->SetFillColor(kWhite);
                    leg->SetFillStyle(1001);
                    leg->SetTextFont(42);
                    leg->SetTextSize(0.040);
                    if (g_sp18) leg->AddEntry(g_sp18, "Sp18 out (blue)", "lep");
                    if (g_fa18) leg->AddEntry(g_fa18, "Fa18 out (red)",  "lep");
                    leg->Draw();
                }
            }
        }

        const std::string outP = (fs::path(output_dir) / ("xsecU_sp18out_vs_fa18out_xB_" + std::to_string(ix) + ".png")).string();
        c->SaveAs(outP.c_str());
        delete c;
        std::cout << "[compare] Wrote " << outP << "\n";
    }

    std::cout << "[compare] Finished Sp18 out vs Fa18 out unpolarized cross-section comparison.\n";
}