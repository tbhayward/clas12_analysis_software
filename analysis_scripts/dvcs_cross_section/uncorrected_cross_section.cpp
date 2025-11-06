#include "uncorrected_cross_section.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TAxis.h>
#include <TPaveText.h>  
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

// ------------ Target properties for luminosity calculation ------------
constexpr double TARGET_DENSITY = 0.07151;        // g/cm³
constexpr double TARGET_LENGTH = 5.0;             // cm
constexpr double AVOGADRO = 6.022e23;             // mol⁻¹
constexpr double ATOMIC_WEIGHT = 1.00794;         // g/mol
constexpr double ELEMENTARY_CHARGE = 1.60217662e-19; // C

// Cross section conversion factors
constexpr double BARN_TO_CM2 = 1e-24;             // 1 barn = 1e-24 cm²
constexpr double NANOBARN_TO_CM2 = 1e-33;         // 1 nb = 1e-33 cm²
constexpr double CM2_TO_NANOBARN = 1e33;          // 1 cm² = 1e33 nb

// ADD THIS - inverse luminosity conversion
constexpr double CM2LUMI_TO_NANOBARNLUMI = 1e-33; // 1 cm⁻² = 1e-33 nb⁻¹

// Charge unit fix
constexpr double NANOCOULOMB_TO_COULOMB = 1e-9;   // 1 nC = 1e-9 C

// Calculate luminosity from charge (in NANOCOULOMBS)
static inline double charge_to_luminosity(double charge_nC) {
    // Convert nC to C first
    double charge_C = charge_nC * NANOCOULOMB_TO_COULOMB;
    // L = (Q/e) * (N_A * ρ * l) / A_w
    return (charge_C / ELEMENTARY_CHARGE) * (AVOGADRO * TARGET_DENSITY * TARGET_LENGTH) / ATOMIC_WEIGHT;
}

// Convert cross section from cm² to nb
static inline double cm2_to_nanobarn(double xsec_cm2) {
    return xsec_cm2 * CM2_TO_NANOBARN;
}

// Convert cross section from nb to cm²  
static inline double nanobarn_to_cm2(double xsec_nb) {
    return xsec_nb * NANOBARN_TO_CM2;
}

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
struct LuminosityData {
    double total_charge;  // C (from column 2)
    double pos_charge;    // C (from column 3)  
    double neg_charge;    // C (from column 4)
    double total_lumi;    // cm⁻² (calculated)
    double pos_lumi;      // cm⁻² (calculated)
    double neg_lumi;      // cm⁻² (calculated)
};

static LuminosityData read_luminosity_data(const std::string& filepath) {
    LuminosityData data{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::ifstream in(filepath);
    if (!in.is_open()) {
        std::cerr << "[luminosity] Cannot open " << filepath << "\n";
        return data;
    }
    auto trim = [](std::string s){
        size_t a = s.find_first_not_of(" \t\r\n");
        size_t b = s.find_last_not_of(" \t\r\n");
        if (a == std::string::npos) return std::string();
        return s.substr(a, b - a + 1);
    };

    double sum_total = 0.0, sum_pos = 0.0, sum_neg = 0.0;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0]=='#') continue;
        std::stringstream ss(line);
        std::string tok; std::vector<std::string> cols;
        while (std::getline(ss, tok, ',')) cols.push_back(trim(tok));
        if (cols.size() < 4) continue;
        try {
            double total = std::stod(cols[1]);
            double pos = std::stod(cols[2]);
            double neg = std::stod(cols[3]);
            if (total > 0) sum_total += total;
            if (pos > 0) sum_pos += pos;
            if (neg > 0) sum_neg += neg;
        } catch (...) {}
    }
    in.close();
    
    data.total_charge = sum_total;
    data.pos_charge = sum_pos;
    data.neg_charge = sum_neg;
    data.total_lumi = charge_to_luminosity(sum_total);
    data.pos_lumi = charge_to_luminosity(sum_pos);
    data.neg_lumi = charge_to_luminosity(sum_neg);
    
    std::cout << "[luminosity] " << filepath
        << "  total_charge=" << sum_total << " nC"
        << "  pos_charge=" << sum_pos << " nC"  
        << "  neg_charge=" << sum_neg << " nC"
        << "  total_lumi=" << data.total_lumi << " cm⁻²"
        << " (" << data.total_lumi * CM2LUMI_TO_NANOBARNLUMI << " nb⁻¹)"  // CHANGED THIS
        << "  pos_lumi=" << data.pos_lumi << " cm⁻²"
        << " (" << data.pos_lumi * CM2LUMI_TO_NANOBARNLUMI << " nb⁻¹)"    // CHANGED THIS
        << "  neg_lumi=" << data.neg_lumi << " cm⁻²"
        << " (" << data.neg_lumi * CM2LUMI_TO_NANOBARNLUMI << " nb⁻¹)\n"; // CHANGED THIS
    return data;
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
// Main function for polarized cross sections
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
    const auto L_fa18_inb = read_luminosity_data(luminosity_dir + "/rga_fa18_inb.txt");
    const auto L_fa18_out = read_luminosity_data(luminosity_dir + "/rga_fa18_out.txt");
    const auto L_sp19_inb = read_luminosity_data(luminosity_dir + "/rga_sp19_inb.txt");

    // Combine luminosities for 10.6 GeV (Fa18 inb + out)
    const double L_10p6_pos = L_fa18_inb.pos_lumi + L_fa18_out.pos_lumi;
    const double L_10p6_neg = L_fa18_inb.neg_lumi + L_fa18_out.neg_lumi;
    const double L_10p6_total = L_fa18_inb.total_lumi + L_fa18_out.total_lumi;
    
    const double L_10p2_pos = L_sp19_inb.pos_lumi;
    const double L_10p2_neg = L_sp19_inb.neg_lumi;
    const double L_10p2_total = L_sp19_inb.total_lumi;

    // Store luminosity data for each energy
    struct EnergyLuminosity {
        double pos_lumi;
        double neg_lumi;
        double total_lumi;
    };
    
    const std::map<std::string, EnergyLuminosity> energyL = {
        {"10.60", {L_10p6_pos, L_10p6_neg, L_10p6_total}},
        {"10.2",  {L_10p2_pos, L_10p2_neg, L_10p2_total}}
    };

    // ----- energy → periods to combine -----
    const std::map<std::string, std::vector<std::string>> energyPeriods = {
        {"10.60", {"DVCS_Fa18_inb", "DVCS_Fa18_out"}},
        {"10.2",  {"DVCS_Sp19_inb"}}
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
        const EnergyLuminosity& L_data = energyL.at(energy);
        const json& jvol = jBinVol[energy];

        if (jvol.empty() || L_data.pos_lumi <= 0.0 || L_data.neg_lumi <= 0.0) {
            std::cerr << "[xsec] Skipping " << energy << " due to missing bin-volume or luminosity.\n";
            continue;
        }

        // ------------ sum unfolded yields (helicity-resolved) across periods ------------
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

        // ------------ compute cross sections ------------
        json jxsec;
        jxsec["energy"] = energy;
        jxsec["luminosity"] = {
            {"pos_lumi", L_data.pos_lumi},
            {"neg_lumi", L_data.neg_lumi},
            {"total_lumi", L_data.total_lumi}
        };
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
                    continue;
                }

                // Calculate denominators for different cross section types
                std::vector<double> denom_pos(N_PHI_BINS, 0.0);
                std::vector<double> denom_neg(N_PHI_BINS, 0.0);
                std::vector<double> denom_total(N_PHI_BINS, 0.0);
                double Vtot = 0.0;
                
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const double Vphi = std::max(0.0, volphi[ip].get<double>());
                    Vtot += Vphi;
                    denom_pos[ip] = L_data.pos_lumi * Vphi;
                    denom_neg[ip] = L_data.neg_lumi * Vphi;
                    denom_total[ip] = L_data.total_lumi * Vphi;
                }

                auto mkOut = [&](const HelData* hd, const std::vector<double>& denom)->json {
                    json o;
                    std::vector<double> xv(N_PHI_BINS), yv(N_PHI_BINS,0.0), ev(N_PHI_BINS,0.0);
                    for (int ip=0; ip<N_PHI_BINS; ++ip) {
                        xv[ip] = PHI_DEG[ip];
                        const double num = (hd ? hd->y[ip]  : 0.0);
                        const double en  = (hd ? hd->ye[ip] : 0.0);
                        const double d   = denom[ip];
                        // Convert from cm² to nb
                        yv[ip] = (d > 0.0 ? cm2_to_nanobarn(num / d) : 0.0);
                        ev[ip] = (d > 0.0 ? cm2_to_nanobarn(en / d) : 0.0);
                    }
                    o["phi"]      = xv;
                    o["xsec"]     = yv;
                    o["xsec_err"] = ev;
                    return o;
                };

                // Calculate unpolarized cross section (sum of plus and minus)
                HelData unpol_data;
                if (itP != sumPlus.end() && itM != sumMinus.end()) {
                    for (int ip=0; ip<N_PHI_BINS; ++ip) {
                        unpol_data.y[ip] = itP->second.y[ip] + itM->second.y[ip];
                        unpol_data.ye[ip] = std::sqrt(itP->second.ye[ip]*itP->second.ye[ip] + 
                                                     itM->second.ye[ip]*itM->second.ye[ip]);
                    }
                } else if (itP != sumPlus.end()) {
                    unpol_data = itP->second;
                } else if (itM != sumMinus.end()) {
                    unpol_data = itM->second;
                }

                jxsec["bins"][bkey] = {
                    {"helicity_plus",  mkOut(itP==sumPlus.end()? nullptr : &itP->second, denom_pos)},
                    {"helicity_minus", mkOut(itM==sumMinus.end()? nullptr : &itM->second, denom_neg)},
                    {"unpolarized",    mkOut(&unpol_data, denom_total)},
                    {"bin_volume_total", Vtot},
                    {"vol_phi", volphi}
                };
            }
        }

        // ------------ save JSON ------------
        const std::string out_json = (fs::path(output_dir)/"jsons" / ("uncorrected_xsec_" + energy + ".json")).string();
        std::ofstream ofs(out_json);
        ofs << std::setw(2) << jxsec << "\n";
        ofs.close();
        std::cout << "[xsec] Wrote " << out_json << "\n";

        // ------------ plots (xB slices; overlay +, −, and unpolarized) ------------
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
                    frame->GetYaxis()->SetTitle("d#sigma/d#phi (nb/GeV^{4})");
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
                    TGraphErrors* gu = drawHel("unpolarized", 22, kGreen+2);

                    TLatex lab;
                    lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                    lab.SetTextFont(42);
                    lab.DrawLatex(0.15, 0.94,
                        Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                             Q2_slice[cc].first, Q2_slice[cc].second,
                             t_slice[r].first,   t_slice[r].second));

                    if (gp || gm || gu){
                        TLegend* leg = new TLegend(0.50, 0.72, 0.90, 0.92);
                        leg->SetBorderSize(1);
                        leg->SetLineColor(kBlack);
                        leg->SetFillColor(kWhite);
                        leg->SetFillStyle(1001);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.040);
                        if (gp) leg->AddEntry(gp, "+ helicity", "lep");
                        if (gm) leg->AddEntry(gm, "- helicity", "lep");
                        if (gu) leg->AddEntry(gu, "unpolarized", "lep");
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

// ------------------------------------------------------------
// Function for unpolarized cross sections using total charge
// ------------------------------------------------------------
void compute_unpolarized_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,
    const std::string& unfolded_counts_dir,
    const std::string& luminosity_dir,
    const std::string& output_dir)
{
    // ----- bin edges & helpers -----
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    const auto PHI_DEG = phiCentersDeg();

    // ----- luminosity for unpolarized (using total charge) -----
    const auto L_sp18_out = read_luminosity_data(luminosity_dir + "/rga_sp18_out.txt");
    const auto L_fa18_inb = read_luminosity_data(luminosity_dir + "/rga_fa18_inb.txt");
    const auto L_fa18_out = read_luminosity_data(luminosity_dir + "/rga_fa18_out.txt");
    const auto L_sp19_inb = read_luminosity_data(luminosity_dir + "/rga_sp19_inb.txt");

    // Combine for 10.6 GeV (Fa18 inb + out)
    const double L_10p6_total = L_fa18_inb.total_lumi + L_fa18_out.total_lumi;
    const double L_10p2_total = L_sp19_inb.total_lumi;
    const double L_10p59_total = L_sp18_out.total_lumi;

    const std::map<std::string, double> energyLumi = {
        {"10.59", L_10p59_total},
        {"10.60", L_10p6_total},
        {"10.2",  L_10p2_total}
    };

    const std::map<std::string, std::vector<std::string>> energyPeriods = {
        {"10.59", {"DVCS_Sp18_out"}},
        {"10.60", {"DVCS_Fa18_inb", "DVCS_Fa18_out"}},
        {"10.2",  {"DVCS_Sp19_inb"}}
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
    fs::create_directories(fs::path(output_dir)/"jsons_unpol");
    fs::create_directories(fs::path(output_dir)/"plots_unpol");

    // ============================================================
    // Loop energies for unpolarized cross sections
    // ============================================================
    for (const auto& eg : energyPeriods) {
        const std::string energy = eg.first;
        const double L_total = energyLumi.at(energy);
        const json& jvol = jBinVol[energy];

        if (jvol.empty() || L_total <= 0.0) {
            std::cerr << "[xsec_unpol] Skipping " << energy << " due to missing bin-volume or luminosity.\n";
            continue;
        }

        // ------------ sum unfolded yields (total) across periods ------------
        struct YieldData { std::vector<double> y, ye; YieldData(): y(N_PHI_BINS,0.0), ye(N_PHI_BINS,0.0){} };
        using CellKey = std::tuple<int,int,int>;
        std::map<CellKey, YieldData> sumTotal;

        for (const std::string& period : eg.second) {
            const std::string path = (fs::path(unfolded_counts_dir) / ("unfolded_"+period+".json")).string();
            json ju = load_json(path);
            if (ju.empty()) {
                std::cerr << "[xsec_unpol] Missing unfolded " << path << " (skipping this period)\n";
                continue;
            }

            if (!ju.contains("bins")) continue;
            for (auto it = ju["bins"].begin(); it != ju["bins"].end(); ++it) {
                const std::string bkey = it.key();
                const json& cell = it.value();

                int ix=0,iQ=0,itb=0;
                if (std::sscanf(bkey.c_str(),"(%d,%d,%d)",&ix,&iQ,&itb) != 3) continue;

                YieldData& yd = sumTotal[CellKey(ix,iQ,itb)];
                
                // Try to get total yield first, otherwise sum plus and minus
                if (cell.contains("yield") && cell["yield"].size() == N_PHI_BINS) {
                    const auto& yt = cell["yield"];
                    const auto& et = cell.contains("yield_err") ? cell["yield_err"] : json::array();
                    for (int ip=0; ip<N_PHI_BINS; ++ip) {
                        double vt = yt[ip].get<double>();
                        double etv = (et.size()==N_PHI_BINS ? std::max(0.0, et[ip].get<double>()) : std::sqrt(std::max(0.0, vt)));
                        yd.y[ip] += vt;
                        yd.ye[ip] = std::sqrt(yd.ye[ip]*yd.ye[ip] + etv*etv);
                    }
                } else if (cell.contains("yield_plus") && cell.contains("yield_minus")) {
                    const auto& yp = cell["yield_plus"];
                    const auto& ym = cell["yield_minus"];
                    const auto& ep = cell.contains("yield_plus_err") ? cell["yield_plus_err"] : json::array();
                    const auto& em = cell.contains("yield_minus_err") ? cell["yield_minus_err"] : json::array();
                    if (yp.size() == N_PHI_BINS && ym.size() == N_PHI_BINS) {
                        for (int ip=0; ip<N_PHI_BINS; ++ip) {
                            double vp = yp[ip].get<double>();
                            double vm = ym[ip].get<double>();
                            double epv = (ep.size()==N_PHI_BINS ? std::max(0.0, ep[ip].get<double>()) : std::sqrt(std::max(0.0, vp)));
                            double emv = (em.size()==N_PHI_BINS ? std::max(0.0, em[ip].get<double>()) : std::sqrt(std::max(0.0, vm)));
                            yd.y[ip] += (vp + vm);
                            yd.ye[ip] = std::sqrt(yd.ye[ip]*yd.ye[ip] + epv*epv + emv*emv);
                        }
                    }
                }
            }
        }

        // ------------ compute unpolarized cross sections ------------
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

                int ix=0,iQ=0,itb=0;
                if (std::sscanf(bkey.c_str(),"(%d,%d,%d)",&ix,&iQ,&itb) != 3) continue;

                auto itY = sumTotal.find(std::make_tuple(ix,iQ,itb));
                if (itY == sumTotal.end()) continue;

                std::vector<double> denom(N_PHI_BINS, 0.0);
                double Vtot = 0.0;
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const double Vphi = std::max(0.0, volphi[ip].get<double>());
                    Vtot += Vphi;
                    denom[ip] = L_total * Vphi;
                }

                json o;
                std::vector<double> xv(N_PHI_BINS), yv(N_PHI_BINS,0.0), ev(N_PHI_BINS,0.0);
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    xv[ip] = PHI_DEG[ip];
                    const double num = itY->second.y[ip];
                    const double en  = itY->second.ye[ip];
                    const double d   = denom[ip];
                    // Convert from cm² to nb
                    yv[ip] = (d > 0.0 ? cm2_to_nanobarn(num / d) : 0.0);
                    ev[ip] = (d > 0.0 ? cm2_to_nanobarn(en / d) : 0.0);
                }
                o["phi"]      = xv;
                o["xsec"]     = yv;
                o["xsec_err"] = ev;

                jxsec["bins"][bkey] = {
                    {"unpolarized", o},
                    {"bin_volume_total", Vtot},
                    {"vol_phi", volphi}
                };
            }
        }

        // ------------ save JSON ------------
        const std::string out_json = (fs::path(output_dir)/"jsons_unpol" / ("unpol_xsec_" + energy + ".json")).string();
        std::ofstream ofs(out_json);
        ofs << std::setw(2) << jxsec << "\n";
        ofs.close();
        std::cout << "[xsec_unpol] Wrote " << out_json << "\n";

        // ------------ plots ------------
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

            std::ostringstream cname; cname<<"c_xsec_unpol_"<<energy<<"_xB"<<ix;
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
            tit << Form("Unpolarized Cross Section  %s GeV   x_{B} #in [%.2g, %.2g]",
                        energy.c_str(), xb.first, xb.second);
            head.DrawLatex(0.5, 0.55, tit.str().c_str());

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
                    frame->GetYaxis()->SetTitle("d#sigma_{U}/d#phi (nb/GeV^{4})");
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
                    const json& jb = jxsec["bins"][bkey]["unpolarized"];

                    if (!jb.contains("phi") || !jb.contains("xsec") || !jb.contains("xsec_err")) continue;
                    const auto& xp = jb["phi"];
                    const auto& yp = jb["xsec"];
                    const auto& ep = jb["xsec_err"];
                    if (xp.size()!=N_PHI_BINS || yp.size()!=N_PHI_BINS || ep.size()!=N_PHI_BINS) continue;

                    std::vector<double> x(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
                    for (int i=0;i<N_PHI_BINS;++i){
                        x[i]=xp[i].get<double>(); y[i]=yp[i].get<double>(); e[i]=std::max(1e-12, ep[i].get<double>());
                    }
                    frame->GetYaxis()->SetRangeUser(1e-4, 1e3);
                    auto* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, e.data());
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerSize(1.0);
                    gr->SetLineWidth(2);
                    gr->SetLineColor(kBlack);
                    gr->SetMarkerColor(kBlack);
                    gr->Draw("P SAME");

                    TLatex lab;
                    lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                    lab.SetTextFont(42);
                    lab.DrawLatex(0.15, 0.94,
                        Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                             Q2_slice[cc].first, Q2_slice[cc].second,
                             t_slice[r].first,   t_slice[r].second));
                }
            }

            const std::string outP = (fs::path(output_dir)/"plots_unpol"/("unpol_xsec_"+energy+"_xB_"+std::to_string(ix)+".png")).string();
            c->SaveAs(outP.c_str());
            delete c;
        }
    }

    std::cout << "[xsec_unpol] Finished unpolarized cross-section generation.\n";
}

// ------------------------------------------------------------
// Comparison function (updated with proper luminosity calculation)
// ------------------------------------------------------------
void compare_unpolarized_cross_sections_sp18out_vs_fa18out(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,
    const std::string& unfolded_counts_dir,
    const std::string& luminosity_dir,
    const std::string& output_dir)
{
    namespace fs = std::filesystem;

    // Read luminosity using our new function
    auto read_total_luminosity = [](const std::string& filepath)->double {
        LuminosityData data = read_luminosity_data(filepath);
        return data.total_lumi;
    };

    // ---------- bin edges/helpers ----------
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    const auto PHI_DEG = phiCentersDeg();

    // ---------- period and energy keys ----------
    const std::string P_SP18_OUT = "DVCS_Sp18_out";
    const std::string P_FA18_OUT = "DVCS_Fa18_out";
    const std::string E_SP18 = "10.59";
    const std::string E_FA18 = "10.60";

    // ---------- luminosities (using total charge and proper conversion) ----------
    const double L_sp18_out_total = read_total_luminosity(luminosity_dir + "/rga_sp18_out.txt");
    const double L_fa18_out_total = read_total_luminosity(luminosity_dir + "/rga_fa18_out.txt");

    // ---------- load bin-volume JSONs ----------
    const std::string path_vol_sp18 = (fs::path(bin_volume_json_dir) / ("bin_volume_" + E_SP18 + ".json")).string();
    const std::string path_vol_fa18 = (fs::path(bin_volume_json_dir) / ("bin_volume_" + E_FA18 + ".json")).string();
    const json jvol_sp18 = load_json(path_vol_sp18);
    const json jvol_fa18 = load_json(path_vol_fa18);
    if (jvol_sp18.empty()) std::cerr << "[compare] Missing bin-volume JSON " << path_vol_sp18 << "\n";
    if (jvol_fa18.empty()) std::cerr << "[compare] Missing bin-volume JSON " << path_vol_fa18 << "\n";

    // ---------- load unfolded counts per period ----------
    const std::string path_u_sp18 = (fs::path(unfolded_counts_dir) / ("unfolded_" + P_SP18_OUT + ".json")).string();
    const std::string path_u_fa18 = (fs::path(unfolded_counts_dir) / ("unfolded_" + P_FA18_OUT + ".json")).string();
    const json ju_sp18 = load_json(path_u_sp18);
    const json ju_fa18 = load_json(path_u_fa18);
    if (ju_sp18.empty()) std::cerr << "[compare] Missing unfolded " << path_u_sp18 << "\n";
    if (ju_fa18.empty()) std::cerr << "[compare] Missing unfolded " << path_u_fa18 << "\n";

    // ---------- build unpolarized xsec arrays for a period ----------
    struct CellData { std::vector<double> xsec, xerr; bool ok=false; CellData(): xsec(N_PHI_BINS,0.0), xerr(N_PHI_BINS,0.0){} };
    using CellKey = std::tuple<int,int,int>;

    auto build_unpol = [&](const json& jvol, const json& ju, double Ltot)->std::map<CellKey, CellData> {
        std::map<CellKey, CellData> out;
        if (jvol.empty() || ju.empty() || !jvol.contains("bins")) return out;

        auto read_total_yield = [&](const json& cell, std::vector<double>& Y, std::vector<double>& Ye)->bool {
            Y.assign(N_PHI_BINS, 0.0);
            Ye.assign(N_PHI_BINS, 0.0);

            bool have_plus  = cell.contains("yield_plus")  && (int)cell["yield_plus"].size()==N_PHI_BINS;
            bool have_minus = cell.contains("yield_minus") && (int)cell["yield_minus"].size()==N_PHI_BINS;
            bool have_total = cell.contains("yield")       && (int)cell["yield"].size()==N_PHI_BINS;

            if (!(have_plus || have_minus || have_total)) return false;

            if (have_plus || have_minus) {
                const auto& yp = have_plus  ? cell["yield_plus"]  : json::array();
                const auto& ym = have_minus ? cell["yield_minus"] : json::array();
                const auto& ep = (cell.contains("yield_plus_err")  && (int)cell["yield_plus_err"].size()==N_PHI_BINS)  ? cell["yield_plus_err"]  : json::array();
                const auto& em = (cell.contains("yield_minus_err") && (int)cell["yield_minus_err"].size()==N_PHI_BINS) ? cell["yield_minus_err"] : json::array();

                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const double vp = have_plus  ? yp[ip].get<double>() : 0.0;
                    const double vm = have_minus ? ym[ip].get<double>() : 0.0;
                    const double epv = (ep.is_array() ? std::max(0.0, ep[ip].get<double>()) : std::sqrt(std::max(0.0, vp)));
                    const double emv = (em.is_array() ? std::max(0.0, em[ip].get<double>()) : std::sqrt(std::max(0.0, vm)));
                    Y[ip]  = vp + vm;
                    Ye[ip] = std::sqrt(epv*epv + emv*emv);
                }
            } else {
                const auto& y  = cell["yield"];
                const bool has_e = cell.contains("yield_err") && (int)cell["yield_err"].size()==N_PHI_BINS;
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    Y[ip]  = y[ip].get<double>();
                    Ye[ip] = has_e ? std::max(0.0, cell["yield_err"][ip].get<double>())
                                   : std::sqrt(std::max(0.0, Y[ip]));
                }
            }
            return true;
        };

        for (auto it = jvol["bins"].begin(); it != jvol["bins"].end(); ++it) {
            const std::string bkey = it.key();
            const json& bv = it.value();
            if (!bv.contains("vol")) continue;
            if ((int)bv["vol"].size() != N_PHI_BINS) continue;

            int ix=0,iQ=0,itb=0;
            if (std::sscanf(bkey.c_str(),"(%d,%d,%d)",&ix,&iQ,&itb) != 3) continue;

            if (!ju.contains("bins") || !ju["bins"].contains(bkey)) continue;
            const json& cell = ju["bins"][bkey];

            std::vector<double> Y, Ye;
            if (!read_total_yield(cell, Y, Ye)) continue;

            CellData cd;
            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                const double Vphi = std::max(0.0, bv["vol"][ip].get<double>());
                const double denom = Ltot * Vphi;
                if (denom > 0.0) {
                    // Convert from cm^2 to nb
                    cd.xsec[ip] = cm2_to_nanobarn(Y[ip]  / denom);
                    cd.xerr[ip] = cm2_to_nanobarn(Ye[ip] / denom);
                } else {
                    cd.xsec[ip] = 0.0;
                    cd.xerr[ip] = 0.0;
                }
            }
            cd.ok = true;
            out[CellKey(ix,iQ,itb)] = std::move(cd);
        }
        return out;
    };

    const auto map_sp18 = build_unpol(jvol_sp18, ju_sp18, L_sp18_out_total);
    const auto map_fa18 = build_unpol(jvol_fa18, ju_fa18, L_fa18_out_total);

    // ---------- output dir ----------
    fs::create_directories(fs::path(output_dir)/"plots_compare");

    // ---------- helper: fetch a cell ----------
    auto fetchCell = [&](const std::map<CellKey,CellData>& M, int ix, int iQ, int itb)->const CellData* {
        auto it = M.find(std::make_tuple(ix,iQ,itb));
        if (it == M.end() || !it->second.ok) return nullptr;
        return &it->second;
    };

    // ---------- helper: per-subplot phi-averaged ratio and uncertainty ----------
    auto phiAveragedRatio = [](const CellData* A, const CellData* B, double& R, double& eR)->bool {
        if (!A || !B) return false;
        double sumA = 0.0, sumB = 0.0, sumVarA = 0.0, sumVarB = 0.0;
        for (int i = 0; i < N_PHI_BINS; ++i) {
            const double a = A->xsec[i];
            const double b = B->xsec[i];
            const double ea = A->xerr[i];
            const double eb = B->xerr[i];
            if (!(a > 0.0 && b > 0.0)) continue;
            if (!std::isfinite(a) || !std::isfinite(b)) continue;
            sumA += a;
            sumB += b;
            sumVarA += ea * ea;
            sumVarB += eb * eb;
        }
        if (!(sumA > 0.0 && sumB > 0.0)) return false;
        R  = sumA / sumB;
        const double sigmaA = std::sqrt(std::max(0.0, sumVarA));
        const double sigmaB = std::sqrt(std::max(0.0, sumVarB));
        eR = R * std::sqrt( std::pow(sigmaA / std::max(1e-12, sumA), 2.0)
                          + std::pow(sigmaB / std::max(1e-12, sumB), 2.0) );
        return std::isfinite(R) && std::isfinite(eR);
    };

    // ---------- global accumulators for the final single-number ratio ----------
    long long n_points_integrated = 0;
    double sumA = 0.0, sumB = 0.0;

    long long n_points_weighted = 0;
    double wsum = 0.0;
    double wrsum = 0.0;

    // ---------- loop over xB bins (existing canvases) ----------
    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();
        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        // Cross-section canvas
        {
            std::ostringstream cname; cname << "c_xsec_unpol_compare_sp18out_fa18out_xB" << ix;
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
            head.SetTextSize(0.30);
            std::ostringstream tit;
            tit << Form("#sigma_{U}   Sp18 out (blue) vs Fa18 out (red)   x_{B} #in [%.2g, %.2g]",
                        xb.first, xb.second);
            head.DrawLatex(0.5, 0.58, tit.str().c_str());

            for (int r=0; r<nrows; ++r) {
                const int it_global = findIndex(t_slice[r], t_bins);
                if (it_global < 0) continue;

                for (int cc=0; cc<ncols; ++cc) {
                    const int iQ_global = findIndex(Q2_slice[cc], Q2_bins);
                    if (iQ_global < 0) continue;

                    pGrid->cd(r*ncols + cc + 1);
                    gPad->SetGrid(1,1);
                    gPad->SetTopMargin(0.08);
                    gPad->SetBottomMargin(0.18);
                    gPad->SetLeftMargin(0.16);
                    gPad->SetRightMargin(0.10);
                    gPad->SetLogy();

                    TH1* frame = gPad->DrawFrame(0.0, 1e-4, 360.0, 1e3);
                    frame->GetXaxis()->SetLabelSize(0.0001);
                    frame->GetXaxis()->SetTitle("#phi (deg)");
                    frame->GetYaxis()->SetTitle("d#sigma_{U}/d#phi (nb/GeV^{4})");
                    frame->GetXaxis()->CenterTitle();
                    frame->GetYaxis()->CenterTitle();
                    frame->GetXaxis()->SetNdivisions(505);
                    frame->GetXaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetLabelSize(0.048);
                    frame->GetXaxis()->SetTitleOffset(1.25);
                    frame->GetYaxis()->SetTitleOffset(1.42);

                    TGaxis* ax = new TGaxis(gPad->GetUxmin(), gPad->GetUymin(),
                                            gPad->GetUxmax(), gPad->GetUymin(),
                                            0.0, 360.0, 4, "");
                    ax->SetLabelFont(42);
                    ax->SetLabelSize(0.048);
                    ax->SetLabelOffset(0.012);
                    ax->SetTitle("");
                    ax->SetTickSize(0.02);
                    ax->Draw();

                    const CellData* cd_sp18 = fetchCell(map_sp18, ix, iQ_global, it_global);
                    const CellData* cd_fa18 = fetchCell(map_fa18, ix, iQ_global, it_global);
                    if (!cd_sp18 && !cd_fa18) continue;

                    auto drawSet = [&](const CellData* cd, int mstyle, int color)->TGraphErrors* {
                        if (!cd) return (TGraphErrors*)nullptr;
                        std::vector<double> x(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
                        for (int i=0;i<N_PHI_BINS;++i) {
                            x[i] = PHI_DEG[i];
                            y[i] = cd->xsec[i];
                            e[i] = std::max(1e-12, cd->xerr[i]);
                        }
                        auto* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, e.data());
                        gr->SetMarkerStyle(mstyle);
                        gr->SetMarkerSize(1.0);
                        gr->SetLineWidth(2);
                        gr->SetLineColor(color);
                        gr->SetMarkerColor(color);
                        gr->Draw("P SAME");
                        return gr;
                    };

                    TGraphErrors* g_sp18 = drawSet(cd_sp18, 20, kBlue+1);
                    TGraphErrors* g_fa18 = drawSet(cd_fa18, 25, kRed+1);

                    TLatex lab;
                    lab.SetNDC(); lab.SetTextSize(0.046); lab.SetTextAlign(11);
                    lab.SetTextFont(42);
                    lab.DrawLatex(0.15, 0.94,
                        Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                             Q2_slice[cc].first, Q2_slice[cc].second,
                             t_slice[r].first,   t_slice[r].second));

                    if (g_sp18 || g_fa18) {
                        TLegend* leg = new TLegend(0.58, 0.70, 0.90, 0.90);
                        leg->SetBorderSize(1);
                        leg->SetLineColor(kBlack);
                        leg->SetFillColor(kWhite);
                        leg->SetFillStyle(1001);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.042);
                        leg->SetMargin(0.10);
                        if (g_sp18) leg->AddEntry(g_sp18, "Sp18 out (total L)", "lep");
                        if (g_fa18) leg->AddEntry(g_fa18, "Fa18 out (total L)", "lep");
                        leg->Draw();
                    }

                    // Phi-averaged ratio box
                    double Rphi = 0.0, eRphi = 0.0;
                    if (phiAveragedRatio(cd_sp18, cd_fa18, Rphi, eRphi)) {
                        TPaveText* box = new TPaveText(0.16, 0.18, 0.58, 0.30, "NDC");
                        box->SetFillColor(kWhite);
                        box->SetFillStyle(1001);
                        box->SetBorderSize(1);
                        box->SetLineColor(kBlack);
                        box->SetLineWidth(2);
                        box->SetShadowColor(0);
                        box->SetTextFont(42);
                        box->SetTextSize(0.040);
                        box->SetTextAlign(12);
                        box->AddText(Form("Sp18/Fa18 = %.3f #pm %.3f", Rphi, eRphi));
                        box->Draw("same");
                    }
                }
            }

            const std::string outP = (fs::path(output_dir)/"plots_compare"/("uncorr_xsec_compare_sp18out_fa18out_xB_"+std::to_string(ix)+".png")).string();
            c->SaveAs(outP.c_str());
            delete c;
            std::cout << "[compare] Wrote " << outP << "\n";
        }

        // Ratio canvas
        {
            std::ostringstream cname; cname << "c_xsec_ratio_sp18_over_fa18_xB" << ix;
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
            head.SetTextSize(0.30);
            std::ostringstream tit;
            tit << Form("Ratio vs #phi: x_{B} #in [%.2g, %.2g]",
                        xb.first, xb.second);
            head.DrawLatex(0.5, 0.58, tit.str().c_str());

            for (int r=0; r<nrows; ++r) {
                const int it_global = findIndex(t_slice[r], t_bins);
                if (it_global < 0) continue;

                for (int cc=0; cc<ncols; ++cc) {
                    const int iQ_global = findIndex(Q2_slice[cc], Q2_bins);
                    if (iQ_global < 0) continue;

                    pGrid->cd(r*ncols + cc + 1);
                    gPad->SetGrid(1,1);
                    gPad->SetTopMargin(0.08);
                    gPad->SetBottomMargin(0.18);
                    gPad->SetLeftMargin(0.18);
                    gPad->SetRightMargin(0.10);
                    gPad->SetLogy(0);

                    TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 200.0);
                    frame->GetXaxis()->SetLabelSize(0.0001);
                    frame->GetXaxis()->SetTitle("#phi (deg)");
                    frame->GetYaxis()->SetTitle(" #sigma_{U}^{Sp18} / #sigma_{U}^{Fa18} (%)");
                    frame->GetXaxis()->CenterTitle();
                    frame->GetYaxis()->CenterTitle();
                    frame->GetXaxis()->SetNdivisions(505);
                    frame->GetXaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetLabelSize(0.048);
                    frame->GetXaxis()->SetTitleOffset(1.25);
                    frame->GetYaxis()->SetTitleOffset(1.42);

                    TGaxis* ax = new TGaxis(gPad->GetUxmin(), gPad->GetUymin(),
                                            gPad->GetUxmax(), gPad->GetUymin(),
                                            0.0, 360.0, 4, "");
                    ax->SetLabelFont(42);
                    ax->SetLabelSize(0.048);
                    ax->SetLabelOffset(0.012);
                    ax->SetTitle("");
                    ax->SetTickSize(0.02);
                    ax->Draw();

                    // unity reference line
                    TLine* unity = new TLine(0.0, 100.0, 360.0, 100.0);
                    unity->SetLineStyle(2);
                    unity->SetLineWidth(2);
                    unity->SetLineColor(kGray+2);
                    unity->Draw("SAME");

                    const CellData* cd_sp18 = fetchCell(map_sp18, ix, iQ_global, it_global);
                    const CellData* cd_fa18 = fetchCell(map_fa18, ix, iQ_global, it_global);

                    TLatex lab;
                    lab.SetNDC(); lab.SetTextSize(0.046); lab.SetTextAlign(11);
                    lab.SetTextFont(42);
                    lab.DrawLatex(0.15, 0.94,
                        Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                             Q2_slice[cc].first, Q2_slice[cc].second,
                             t_slice[r].first,   t_slice[r].second));

                    if (!cd_sp18 || !cd_fa18) continue;

                    // compute ratio vs phi
                    std::vector<double> xp, yp, ey;
                    xp.reserve(N_PHI_BINS);
                    yp.reserve(N_PHI_BINS);
                    ey.reserve(N_PHI_BINS);

                    for (int ip=0; ip<N_PHI_BINS; ++ip) {
                        const double A  = cd_sp18->xsec[ip];
                        const double B  = cd_fa18->xsec[ip];
                        const double eA = cd_sp18->xerr[ip];
                        const double eB = cd_fa18->xerr[ip];

                        if (!(A > 0.0 && B > 0.0)) continue;
                        if (!std::isfinite(A) || !std::isfinite(B)) continue;

                        const double R  = A / B;
                        const double eR = R * std::sqrt( std::pow(eA/std::max(1e-12,A), 2.0)
                                                       + std::pow(eB/std::max(1e-12,B), 2.0) );

                        xp.push_back(PHI_DEG[ip]);
                        yp.push_back(100.0 * R);
                        ey.push_back(100.0 * eR);

                        // global accumulators
                        sumA += A;
                        sumB += B;
                        ++n_points_integrated;

                        const double varR = eR * eR;
                        if (std::isfinite(varR) && varR > 0.0) {
                            const double w = 1.0 / varR;
                            wsum  += w;
                            wrsum += w * R;
                            ++n_points_weighted;
                        }
                    }

                    if (!xp.empty()) {
                        auto* gr = new TGraphErrors((int)xp.size(), xp.data(), yp.data(), nullptr, ey.data());
                        gr->SetMarkerStyle(21);
                        gr->SetMarkerSize(1.0);
                        gr->SetLineWidth(2);
                        gr->SetLineColor(kBlack);
                        gr->SetMarkerColor(kBlack);
                        gr->Draw("P SAME");
                    }

                    // Phi-averaged ratio box
                    double Rphi = 0.0, eRphi = 0.0;
                    if (phiAveragedRatio(cd_sp18, cd_fa18, Rphi, eRphi)) {
                        TPaveText* box = new TPaveText(0.18, 0.18, 0.58, 0.30, "NDC");
                        box->SetFillColor(kWhite);
                        box->SetFillStyle(1001);
                        box->SetBorderSize(1);
                        box->SetLineColor(kBlack);
                        box->SetLineWidth(2);
                        box->SetShadowColor(0);
                        box->SetTextFont(42);
                        box->SetTextSize(0.040);
                        box->SetTextAlign(12);
                        box->AddText(Form("Sp18/Fa18 = %.3f #pm %.3f", Rphi, eRphi));
                        box->Draw("same");
                    }
                }
            }

            const std::string outP2 = (fs::path(output_dir)/"plots_compare"/("uncorr_xsec_ratio_sp18_over_fa18_xB_"+std::to_string(ix)+".png")).string();
            c->SaveAs(outP2.c_str());
            delete c;
            std::cout << "[compare] Wrote " << outP2 << "\n";
        }
    }

    // ---------- NEW: single "summed low-xB" plot (xB_max <= 0.15), COMMON PHASE SPACE ----------
    {
        // 1) choose xB bins: xmax <= 0.15
        std::vector<int> ix_include;
        for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
            const auto& xb = xB_bins[ix];
            if (xb.second <= 0.15) ix_include.push_back(ix);
        }
        auto ix_allowed = [&](int ix)->bool {
            return std::find(ix_include.begin(), ix_include.end(), ix) != ix_include.end();
        };

        // 2) intersection of cells present in BOTH periods within selected xB
        using CellKey = std::tuple<int,int,int>;
        std::set<CellKey> keys_sp18, keys_fa18, keys_inter;

        for (const auto& kv : map_sp18) {
            int ix,iQ,itb; std::tie(ix,iQ,itb) = kv.first;
            if (ix_allowed(ix)) keys_sp18.insert(kv.first);
        }
        for (const auto& kv : map_fa18) {
            int ix,iQ,itb; std::tie(ix,iQ,itb) = kv.first;
            if (ix_allowed(ix)) keys_fa18.insert(kv.first);
        }
        std::set_intersection(keys_sp18.begin(), keys_sp18.end(),
                              keys_fa18.begin(), keys_fa18.end(),
                              std::inserter(keys_inter, keys_inter.begin()));

        // 3) sum dσU/dφ over the common cells only
        std::vector<double> y_sp18(N_PHI_BINS, 0.0), e2_sp18(N_PHI_BINS, 0.0);
        std::vector<double> y_fa18(N_PHI_BINS, 0.0), e2_fa18(N_PHI_BINS, 0.0);

        for (const auto& key : keys_inter) {
            const auto& cdA = map_sp18.at(key);
            const auto& cdB = map_fa18.at(key);
            for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                y_sp18[ip]  += cdA.xsec[ip];
                e2_sp18[ip] += cdA.xerr[ip] * cdA.xerr[ip];
                y_fa18[ip]  += cdB.xsec[ip];
                e2_fa18[ip] += cdB.xerr[ip] * cdB.xerr[ip];
            }
        }

        bool any = false; for (int ip=0; ip<N_PHI_BINS; ++ip) if (y_sp18[ip] > 0.0 || y_fa18[ip] > 0.0) { any = true; break; }
        if (any) {
            std::vector<double> x(N_PHI_BINS), es_sp18(N_PHI_BINS,0.0), es_fa18(N_PHI_BINS,0.0);
            for (int i=0;i<N_PHI_BINS;++i) {
                x[i] = PHI_DEG[i];
                es_sp18[i] = std::sqrt(std::max(0.0, e2_sp18[i]));
                es_fa18[i] = std::sqrt(std::max(0.0, e2_fa18[i]));
            }

            const int W = 900, H = 650;
            TCanvas* c = new TCanvas("c_sum_lowxB_phi", "c_sum_lowxB_phi", W, H);

            // top title pad
            TPad* pTop  = new TPad("pTop_lowx","pTop_lowx", 0.0, 0.90, 1.0, 1.0);
            pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
            pTop->cd();
            TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.30);
            head.DrawLatex(0.5, 0.55, "#sigma_{U} vs #phi (x_{B} < 0.15, common phase space)");

            // main pad
            c->cd();
            TPad* pMain = new TPad("pMain_lowx","pMain_lowx", 0.0, 0.00, 1.0, 0.90);
            pMain->SetFillStyle(0); pMain->SetBorderSize(0); pMain->Draw();
            pMain->cd();
            gPad->SetGrid(1,1);
            gPad->SetTopMargin(0.08);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.10);
            gPad->SetLogy();

            // choose a sensible log range
            double ymax = 0.0;
            for (int i=0;i<N_PHI_BINS;++i) {
                ymax = std::max(ymax, y_sp18[i] + es_sp18[i]);
                ymax = std::max(ymax, y_fa18[i] + es_fa18[i]);
            }
            if (!(ymax > 0.0)) ymax = 1.0;
            double ymin = 1.0; // lower bound (nb/GeV^4) for readability on log axis
            double ytop = std::pow(10.0, std::ceil(std::log10(ymax*1.5)));

            TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ytop);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("d#sigma_{U}/d#phi (nb/GeV^{4})");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.048);
            frame->GetXaxis()->SetLabelSize(0.048);   // show default labels
            frame->GetXaxis()->SetLabelOffset(0.012); // just below axis
            frame->GetXaxis()->SetTitleOffset(1.25);
            frame->GetYaxis()->SetTitleOffset(1.42);
            // (No TGaxis here -> avoids extra bottom line)

            auto* gr_sp18 = new TGraphErrors(N_PHI_BINS, x.data(), y_sp18.data(), nullptr, es_sp18.data());
            gr_sp18->SetMarkerStyle(20);
            gr_sp18->SetMarkerSize(1.1);
            gr_sp18->SetLineWidth(2);
            gr_sp18->SetLineColor(kBlue+1);
            gr_sp18->SetMarkerColor(kBlue+1);
            gr_sp18->Draw("P SAME");

            auto* gr_fa18 = new TGraphErrors(N_PHI_BINS, x.data(), y_fa18.data(), nullptr, es_fa18.data());
            gr_fa18->SetMarkerStyle(25);
            gr_fa18->SetMarkerSize(1.1);
            gr_fa18->SetLineWidth(2);
            gr_fa18->SetLineColor(kRed+1);
            gr_fa18->SetMarkerColor(kRed+1);
            gr_fa18->Draw("P SAME");

            TLegend* leg = new TLegend(0.58, 0.70, 0.90, 0.90);
            leg->SetBorderSize(1);
            leg->SetLineColor(kBlack);
            leg->SetFillColor(kWhite);
            leg->SetFillStyle(1001);
            leg->SetTextFont(42);
            leg->SetTextSize(0.042);
            leg->AddEntry(gr_sp18, "Sp18 Outb", "lep");
            leg->AddEntry(gr_fa18, "Fa18 Outb", "lep");
            leg->Draw();

            const std::string outSingle = (fs::path(output_dir)/"plots_compare"/"uncorr_xsec_sum_phi_sp18out_vs_fa18out_xB_lt_0p15.png").string();
            c->SaveAs(outSingle.c_str());
            delete c;
            std::cout << "[compare] Wrote " << outSingle << "\n";

            // optional sanity printout for this specific sum
            double sA=0.0, sB=0.0;
            for (int i=0;i<N_PHI_BINS;++i) { sA += y_sp18[i]; sB += y_fa18[i]; }
            if (sB > 0.0) {
                std::cout << "[compare][sum-lowxB] cells=" << keys_inter.size()
                          << " ; sum_sp18=" << sA << " ; sum_fa18=" << sB
                          << " ; ratio=" << (sA/sB) << "\n";
            }
        }
    }

    // ---------- final single-number summary ----------
    double integrated_ratio_pct = std::numeric_limits<double>::quiet_NaN();
    if (sumB > 0.0) {
        integrated_ratio_pct = 100.0 * (sumA / sumB);
    }

    double weighted_mean_ratio_pct = std::numeric_limits<double>::quiet_NaN();
    double weighted_mean_se_pct    = std::numeric_limits<double>::quiet_NaN();
    if (wsum > 0.0) {
        const double Rw  = wrsum / wsum;
        const double sRw = std::sqrt(1.0 / wsum);
        weighted_mean_ratio_pct = 100.0 * Rw;
        weighted_mean_se_pct    = 100.0 * sRw;
    }

    std::cout << "[compare][final] IntegratedRatioPct = " << integrated_ratio_pct
              << " ; points=" << n_points_integrated
              << " ; sumA=" << sumA << " ; sumB=" << sumB << "\n";

    std::cout << "[compare][final] WeightedMeanRatioPct = " << weighted_mean_ratio_pct
              << " +/- " << weighted_mean_se_pct
              << " ; points=" << n_points_weighted
              << " ; wsum=" << wsum << "\n";

    std::cout << "[compare] Finished unpolarized comparison Sp18 out vs Fa18 out.\n";
}