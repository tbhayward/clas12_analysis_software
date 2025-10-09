#include "pi0_corrected_counts.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1.h>
#include <TLatex.h>
#include <TPad.h>
#include <TGaxis.h>

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
#include <vector>

#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {
constexpr int    N_PHI_BINS = 12;

// ---------- style ----------
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

// ---------- bin helpers ----------
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

// ---------- small IO ----------
static json load_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[pi0corr] Failed to open " << path << "\n";
        return json();
    }
    json j; f >> j; return j;
}

// Period like "DVCS_Fa18_inb" -> "fa18_inb"
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    std::string tail = (pos == std::string::npos || pos+1>=period.size())
        ? period : period.substr(pos+1);
    std::transform(tail.begin(), tail.end(), tail.begin(), ::tolower);
    return tail;
}

// ---------- data models ----------
using BinKey4 = std::tuple<int,int,int,int>; // (ix,iQ,it,ip)

struct Hel2 {
    double plus=0.0, minus=0.0;
};
struct Hel2Err {
    double plus=0.0, minus=0.0;
};

// total_counts reader (group -> (ix,iQ,it,ip) -> Hel2 integer counts)
static bool read_total_counts_group(const std::string& total_counts_json,
                                    const std::string& group_key,
                                    std::map<BinKey4, Hel2>& out)
{
    out.clear();
    json j = load_json(total_counts_json);
    if (j.empty() || !j.contains("groups")) return false;
    const auto& G = j["groups"];
    if (!G.contains(group_key)) {
        std::cerr << "[pi0corr] total_counts group '"<<group_key<<"' not found in "<<total_counts_json<<"\n";
        return false;
    }
    const auto& gg = G[group_key];
    for (auto it = gg.begin(); it != gg.end(); ++it) {
        const std::string key = it.key(); // "(ix,iQ,it,ip)"
        int ix=0,iQ=0,itb=0,ip=0;
        if (std::sscanf(key.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&itb,&ip) != 4) continue;
        const auto& obj = it.value();
        if (!obj.contains("helicity")) continue;
        Hel2 h;
        const auto& H = obj["helicity"];
        if (H.contains("+1")) h.plus  = H["+1"].get<double>();
        if (H.contains("-1")) h.minus = H["-1"].get<double>();
        out[BinKey4(ix,iQ,itb,ip)] = h;
    }
    return !out.empty();
}

// contamination reader: period JSON -> (ix,iQ,it,ip) -> c_h and err
struct ContamVal { double val=0.0, err=0.0; };
struct Contam2 { ContamVal plus, minus; };
static bool read_contamination_period(const std::string& contamination_json_path,
                                      std::map<BinKey4, Contam2>& out)
{
    out.clear();
    std::ifstream ifs(contamination_json_path);
    if (!ifs) { std::cerr << "[pi0corr] Cannot open contamination: "<<contamination_json_path<<"\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    // parse "bins" object by hand (structure written by your contam code)
    size_t pos = s.find("\"bins\"");
    if (pos == std::string::npos) return false;
    pos = s.find('{', pos); if (pos==std::string::npos) return false;
    int depth=0; size_t i=pos;
    for (; i<s.size(); ++i){ if(s[i]=='{') depth++; else if(s[i]=='}'){ depth--; if(!depth){ ++i; break; } } }
    std::string binsObj = s.substr(pos, i-pos);

    size_t kpos=0;
    while (true) {
        size_t ks = binsObj.find('"', kpos);
        if (ks == std::string::npos) break;
        size_t ke = binsObj.find('"', ks+1);
        if (ke == std::string::npos) break;
        std::string key = binsObj.substr(ks+1, ke-ks-1);
        int ix,iQ,itb,ip;
        if (std::sscanf(key.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&itb,&ip) != 4) { kpos = ke+1; continue; }

        size_t os = binsObj.find('{', ke);
        if (os == std::string::npos) break;
        int d2=0; size_t j=os;
        for (; j<binsObj.size(); ++j){ if(binsObj[j]=='{') d2++; else if(binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj = binsObj.substr(os, j-os);

        auto findDouble=[&](const std::string& path)->double{
            size_t p=obj.find(path); if(p==std::string::npos) return 0.0;
            p=obj.find(':',p); if(p==std::string::npos) return 0.0;
            size_t a=p+1; while(a<obj.size() && std::isspace((unsigned char)obj[a])) ++a;
            size_t b=a; while(b<obj.size() && (std::isdigit((unsigned char)obj[b])||obj[b]=='-'||obj[b]=='.'||obj[b]=='e'||obj[b]=='E'||obj[b]=='+')) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch(...) { return 0.0; }
        };

        Contam2 c;
        c.plus.val  = findDouble("\"contamination\":{\"+1\":{\"value\"");
        c.plus.err  = findDouble("\"contamination\":{\"+1\":{\"err\"");
        c.minus.val = findDouble("\"contamination\":{\"-1\":{\"value\"");
        c.minus.err = findDouble("\"contamination\":{\"-1\":{\"err\"");
        out[BinKey4(ix,iQ,itb,ip)] = c;

        kpos = j;
    }
    return !out.empty();
}

// ---------- corrected counts + errors per helicity ----------
struct Corr2 { Hel2 val; Hel2Err err; };

static Corr2 correct_one_bin(const Hel2& raw, const Contam2& c)
{
    // N_corr = N_raw * (1 - c),  Var = (1-c)^2 Var(N_raw) + (N_raw)^2 Var(c)
    // Assume Var(N_raw)=N_raw (Poisson) and Var(c)=err^2
    Corr2 out;
    auto one=[&](double N, const ContamVal& k, double& Nc, double& se){
        const double clamped = std::clamp(k.val, 0.0, 1.0);
        const double varN = std::max(0.0, N);        // Poisson
        const double varc = std::max(0.0, k.err*k.err);
        Nc  = N * (1.0 - clamped);
        double var = (1.0 - clamped)*(1.0 - clamped) * varN + N*N * varc;
        se = std::sqrt(std::max(0.0, var));
    };
    one(raw.plus,  c.plus,  out.val.plus,  out.err.plus);
    one(raw.minus, c.minus, out.val.minus, out.err.minus);
    return out;
}

// ---------- JSON writer per period (arrays over φ) ----------
static void write_period_json(const std::string& path,
                              const std::vector<double>& PHI,
                              const std::map<std::tuple<int,int,int>, std::array<Corr2, N_PHI_BINS>>& table,
                              const std::map<std::tuple<int,int,int>, std::array<Hel2,  N_PHI_BINS>>& rawtab)
{
    std::ofstream ofs(path);
    if(!ofs){ std::cerr << "[pi0corr] Cannot write "<<path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table){
        if(!first) ofs<<",\n"; first=false;
        int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
        const auto& corrArr = kv.second;
        const auto& rawArr  = rawtab.at(kv.first);

        ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {";

        auto dump=[&](auto getter, const char* name){
            ofs<<"\""<<name<<"\":[";
            for (int ip=0; ip<N_PHI_BINS; ++ip){ if(ip) ofs<<","; ofs<<getter(ip); }
            ofs<<"],";
        };

        // phi
        ofs<<"\"phi\":[";
        for (int ip=0; ip<N_PHI_BINS; ++ip){ if(ip) ofs<<","; ofs<<PHI[ip]; }
        ofs<<"],";

        dump([&](int ip){ return rawArr[ip].plus;  }, "raw_plus");
        dump([&](int ip){ return std::sqrt(std::max(0.0, rawArr[ip].plus));  }, "raw_plus_err");
        dump([&](int ip){ return rawArr[ip].minus; }, "raw_minus");
        dump([&](int ip){ return std::sqrt(std::max(0.0, rawArr[ip].minus)); }, "raw_minus_err");

        dump([&](int ip){ return corrArr[ip].val.plus;  }, "corr_plus");
        dump([&](int ip){ return corrArr[ip].err.plus;  }, "corr_plus_err");
        dump([&](int ip){ return corrArr[ip].val.minus; }, "corr_minus");
        dump([&](int ip){ return corrArr[ip].err.minus; }, "corr_minus_err");

        ofs.seekp(-1, std::ios_base::cur); // remove trailing comma
        ofs<<"}";
    }
    ofs<<"\n  }\n}\n";
    std::cout<<"[pi0corr] Wrote "<<path<<"\n";
}

// ---------- plotting (RAW vs CORR overlays, ± helicities) ----------
static void plot_period(const std::string& period,
                        const std::string& runTag,
                        const std::vector<Binning>& binning_scheme,
                        const std::vector<std::pair<double,double>>& xB_bins,
                        const std::vector<std::pair<double,double>>& Q2_bins,
                        const std::vector<std::pair<double,double>>& t_bins,
                        const std::map<std::tuple<int,int,int>, std::array<Corr2, N_PHI_BINS>>& table,
                        const std::map<std::tuple<int,int,int>, std::array<Hel2,  N_PHI_BINS>>& rawtab,
                        const std::string& out_dir)
{
    std::error_code ec;
    fs::create_directories(out_dir, ec);

    const auto PHI = phiCentersDeg();

    auto QT_slice = [&](const std::pair<double,double>& xb,
                        std::vector<std::pair<double,double>>& Q2_list,
                        std::vector<std::pair<double,double>>& t_list){
        std::set<std::pair<double,double>> qs, ts;
        for (const auto& b : binning_scheme) {
            if (std::make_pair(b.xBmin,b.xBmax) == xb) {
                qs.emplace(b.Q2min,b.Q2max);
                ts.emplace(b.tmin,b.tmax);
            }
        }
        Q2_list.assign(qs.begin(), qs.end());
        t_list.assign(ts.begin(), ts.end());
    };

    for (int ix=0; ix<(int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        QT_slice(xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();
        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_pi0corr_"<<runTag<<"_xB"<<ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // Title
        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.36);
        std::ostringstream tit;
        tit << Form("#pi^{0}-corrected counts  (%s)   x_{B} #in [%.2g, %.2g]",
                    period.c_str(), xb.first, xb.second);
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        // Panels
        for (int r=0; r<nrows; ++r) {
            const int it_global = findIndex(t_slice[r], t_bins); if (it_global<0) continue;
            for (int cc=0; cc<ncols; ++cc) {
                const int iQ_global = findIndex(Q2_slice[cc], Q2_bins); if (iQ_global<0) continue;

                pGrid->cd(r*ncols + cc + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                // find arrays for this (ix,iQ,it)
                const auto key3 = std::make_tuple(ix, iQ_global, it_global);
                auto itC = table.find(key3);
                auto itR = rawtab.find(key3);
                if (itC==table.end() || itR==rawtab.end()) {
                    TH1* fr = gPad->DrawFrame(0.0, 0.0, 360.0, 1.0);
                    fr->GetXaxis()->SetTitle("#phi (deg)");
                    fr->GetYaxis()->SetTitle("counts");
                    continue;
                }
                const auto& carr = itC->second;
                const auto& rarr = itR->second;

                // build vectors
                std::vector<double> X(N_PHI_BINS), ex(N_PHI_BINS, 0.0);
                std::vector<double> Yrp(N_PHI_BINS,0), Yrpm(N_PHI_BINS,0), Ycp(N_PHI_BINS,0), Ycm(N_PHI_BINS,0);
                std::vector<double> Erp(N_PHI_BINS,0), Erm(N_PHI_BINS,0), Ecp(N_PHI_BINS,0), Ecm(N_PHI_BINS,0);

                double ymax = 0.0;
                for (int ip=0; ip<N_PHI_BINS; ++ip){
                    X[ip]=PHI[ip];
                    Yrp[ip]=rarr[ip].plus;  Erp[ip]=std::sqrt(std::max(0.0, rarr[ip].plus));
                    Yrpm[ip]=rarr[ip].minus; Erm[ip]=std::sqrt(std::max(0.0, rarr[ip].minus));
                    Ycp[ip]=carr[ip].val.plus;  Ecp[ip]=carr[ip].err.plus;
                    Ycm[ip]=carr[ip].val.minus; Ecm[ip]=carr[ip].err.minus;

                    ymax = std::max({ymax, Yrp[ip]+Erp[ip], Yrpm[ip]+Erm[ip], Ycp[ip]+Ecp[ip], Ycm[ip]+Ecm[ip]});
                }
                if (ymax <= 0.0) ymax = 1.0;

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, ymax*1.25);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("counts");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.048);
                frame->GetXaxis()->SetTitleOffset(1.25);
                frame->GetYaxis()->SetTitleOffset(1.35);
                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                auto gr = [&](const std::vector<double>& y, const std::vector<double>& e, int mstyle, int color){
                    auto* g = new TGraphErrors(N_PHI_BINS, X.data(), const_cast<double*>(y.data()), ex.data(), const_cast<double*>(e.data()));
                    g->SetMarkerStyle(mstyle); g->SetMarkerSize(0.9);
                    g->SetLineWidth(2); g->SetLineColor(color); g->SetMarkerColor(color);
                    g->Draw("P SAME");
                    return g;
                };

                // RAW (gray-ish) and CORR (colored)
                TGraphErrors* grp_raw = gr(Yrp, Erp, 24, kGray+2);
                TGraphErrors* grm_raw = gr(Yrpm, Erm, 24, kGray+1);
                TGraphErrors* grp_cor = gr(Ycp, Ecp, 20, kBlue+1);
                TGraphErrors* grm_cor = gr(Ycm, Ecm, 25, kRed+1);

                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11); lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[cc].first, Q2_slice[cc].second,
                         t_slice[r].first,   t_slice[r].second));

                TLegend* leg = new TLegend(0.50, 0.68, 0.90, 0.92);
                leg->SetBorderSize(1); leg->SetLineColor(kBlack); leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
                leg->SetTextFont(42); leg->SetTextSize(0.040);
                leg->AddEntry(grp_cor, "+ helicity (corr)", "lep");
                leg->AddEntry(grm_cor, "- helicity (corr)", "lep");
                leg->AddEntry(grp_raw, "+ helicity (raw)",  "lep");
                leg->AddEntry(grm_raw, "- helicity (raw)",  "lep");
                leg->Draw();
            }
        }

        const std::string outP = (fs::path(out_dir)/("corr_counts_"+period+"_xB_"+std::to_string(ix)+".png")).string();
        c->SaveAs(outP.c_str());
        delete c;
    }
}

} // anon

// =====================================================================
// Public driver
// =====================================================================
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>&     binning_scheme,
    const std::string&              total_counts_json,
    const std::string&              contamination_json_dir,
    const std::string&              out_root_dir)
{
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    const auto PHI     = phiCentersDeg();

    // output dirs
    const fs::path json_dir = fs::path(out_root_dir)/"jsons";
    const fs::path plots_root = fs::path(out_root_dir)/"pi0_corrected_counts";
    std::error_code ec;
    fs::create_directories(json_dir, ec);

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);

        // 1) load RAW counts for this runTag from total_counts.json (group == runTag)
        std::map<BinKey4, Hel2> raw_by_bin4;
        if (!read_total_counts_group(total_counts_json, runTag, raw_by_bin4)) {
            std::cerr << "[pi0corr] Skipping period "<<period<<": no raw counts for group '"<<runTag<<"'.\n";
            continue;
        }

        // 2) load contamination for this period
        const std::string contam_path = (fs::path(contamination_json_dir)/("contamination_"+period+".json")).string();
        std::map<BinKey4, Contam2> c_by_bin4;
        if (!read_contamination_period(contam_path, c_by_bin4)) {
            std::cerr << "[pi0corr] WARNING: no contamination JSON for "<<period<<" at "<<contam_path<<". Using c=0.\n";
        }

        // 3) correct per φ bin; pack arrays per (ix,iQ,it)
        std::map<std::tuple<int,int,int>, std::array<Corr2, N_PHI_BINS>> corr_by_cell;
        std::map<std::tuple<int,int,int>, std::array<Hel2,  N_PHI_BINS>> raw_by_cell;

        for (int ix=0; ix<(int)xB_bins.size(); ++ix)
        for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
        for (int it=0; it<(int)t_bins.size();  ++it) {
            std::array<Corr2, N_PHI_BINS> carr{};
            std::array<Hel2,  N_PHI_BINS> rarr{};
            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                const auto k4 = std::make_tuple(ix,iQ,it,ip);

                Hel2 raw{0.0,0.0};
                auto itR = raw_by_bin4.find(k4);
                if (itR != raw_by_bin4.end()) raw = itR->second;

                Contam2 c; // defaults 0
                auto itC = c_by_bin4.find(k4);
                if (itC != c_by_bin4.end()) c = itC->second;

                rarr[ip]   = raw;
                carr[ip]   = correct_one_bin(raw, c);
            }
            corr_by_cell[std::make_tuple(ix,iQ,it)] = carr;
            raw_by_cell [std::make_tuple(ix,iQ,it)] = rarr;
        }

        // 4) write JSON
        const std::string out_json = (json_dir/("pi0_corrected_counts_"+period+".json")).string();
        write_period_json(out_json, PHI, corr_by_cell, raw_by_cell);

        // 5) plots
        const std::string out_plots_dir = (plots_root/runTag).string();
        plot_period(period, runTag, binning_scheme, xB_bins, Q2_bins, t_bins, corr_by_cell, raw_by_cell, out_plots_dir);
    }

    std::cout << "[pi0corr] Finished π0-corrected counts.\n";
}