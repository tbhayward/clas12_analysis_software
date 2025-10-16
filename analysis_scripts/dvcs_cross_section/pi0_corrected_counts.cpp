// pi0_corrected_counts.cpp
//
// Reads:
//   - total_counts.json (with top-level "groups")
//   - per-period contamination JSONs: output/jsons/contamination/contamination_<period>.json
// Produces:
//   - per-period corrected counts JSONs: output/jsons/pi0_corrected_counts_<period>.json
//   - combined: output/jsons/pi0_corrected_counts_all_periods.json
//   - plots per period under: output/pi0_corrected_plots/<runTag>/plot_pi0corr_<period>_xB_<ix>.png
//
// Corrected per φ, per helicity:
//   N^corr_h = N^raw_h * (1 - c_h)
// where c_h is the π0 contamination fraction (helicity-resolved).
//
// Plots: one canvas per xB slice; Q²×|t| grid; in each panel we draw
// the corrected N(+), corrected N(-) vs φ (deg) as TGraph.

#include "pi0_corrected_counts.h"

#include <TCanvas.h>
#include <TTree.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TStyle.h>
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

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ,it,ip)

struct HelCounts { long long plus=0, minus=0; };
using GroupCounts = std::map<std::string, std::map<BinKey, HelCounts>>;

struct ContamBin { double c_plus=0, c_plus_err=0, c_minus=0, c_minus_err=0; };

// -------------------- style --------------------
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTextFont(42);
    }
} _style_bootstrap;

// -------------------- utils --------------------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}

static inline std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}
static inline int findIndex(const std::pair<double,double>& r,
                            const std::vector<std::pair<double,double>>& R){
    for (int i=0;i<(int)R.size();++i) if (R[i]==r) return i;
    return -1;
}
static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
}

// ------------- JSON I/O -------------
// total_counts reader (matches writer with top-level "groups")
static bool parse_tuple_key4(const std::string& s, BinKey& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey(ix,iQ,it,ip);
    return true;
}

static bool load_total_counts(const std::string& path, GroupCounts& outGroups) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr<<"[pi0_corr][ERROR] Cannot open total_counts JSON: "<<path<<"\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t gpos = s.find("\"groups\"");
    if (gpos==std::string::npos) { std::cerr<<"[pi0_corr][ERROR] 'groups' not found in total_counts.\n"; return false; }
    size_t brace = s.find('{', gpos); if (brace==std::string::npos) return false;
    int d=0; size_t i=brace; for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string groupsObj = s.substr(brace, i-brace);

    size_t kpos=0;
    while (true) {
        size_t q1 = groupsObj.find('"', kpos); if (q1==std::string::npos) break;
        size_t q2 = groupsObj.find('"', q1+1); if (q2==std::string::npos) break;
        std::string gname = groupsObj.substr(q1+1, q2-q1-1);

        size_t objS = groupsObj.find('{', q2); if (objS==std::string::npos) break;
        int d2=0; size_t j=objS;
        for (; j<groupsObj.size(); ++j){ if(groupsObj[j]=='{') d2++; else if(groupsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string binsObj = groupsObj.substr(objS, j-objS);

        std::map<BinKey, HelCounts> gmap;
        size_t bpos=0;
        while (true) {
            size_t bk1 = binsObj.find('"', bpos); if (bk1==std::string::npos) break;
            size_t bk2 = binsObj.find('"', bk1+1); if (bk2==std::string::npos) break;
            std::string key = binsObj.substr(bk1+1, bk2-bk1-1);
            BinKey bk;
            if (!parse_tuple_key4(key, bk)) { bpos=bk2+1; continue; }

            size_t valS = binsObj.find('{', bk2); if (valS==std::string::npos) break;
            int d3=0; size_t jj=valS;
            for (; jj<binsObj.size(); ++jj){ if(binsObj[jj]=='{') d3++; else if(binsObj[jj]=='}'){ d3--; if(!d3){ ++jj; break; } } }
            std::string obj = binsObj.substr(valS, jj-valS);

            auto findLL = [&](const char* pat)->long long {
                size_t p = obj.find(pat); if (p==std::string::npos) return 0;
                p = obj.find(':', p); if (p==std::string::npos) return 0;
                size_t a=p+1; while (a<obj.size() && isspace((unsigned char)obj[a])) ++a;
                size_t b=a; while (b<obj.size() && (isdigit((unsigned char)obj[b])||obj[b]=='-')) ++b;
                try { return std::stoll(obj.substr(a,b-a)); } catch(...) { return 0; }
            };
            HelCounts hc;
            hc.plus  = findLL("\"+1\"");
            hc.minus = findLL("\"-1\"");

            gmap[bk]=hc;
            bpos=jj;
        }

        outGroups[gname]=std::move(gmap);
        kpos=j;
    }
    return !outGroups.empty();
}

namespace {
struct ContamBin { double c_plus=0, c_plus_err=0, c_minus=0, c_minus_err=0; };
using BinKey = std::tuple<int,int,int,int>;
bool parse_tuple_key(const std::string& s, BinKey& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(), "(%d,%d,%d,%d)", &ix,&iQ,&it,&ip) != 4) return false;
    out = BinKey(ix,iQ,it,ip);
    return true;
}

// Helper: remove all whitespace to make pattern-matching robust
static inline std::string squash_ws(std::string z){
    z.erase(std::remove_if(z.begin(), z.end(), [](unsigned char c){
        return std::isspace(c);
    }), z.end());
    return z;
}

static bool load_contam_for_period(const std::string& path,
                                   std::map<BinKey, ContamBin>& out)
{
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[pi0_corr][WARN] Can't open contamination file: " << path << "\n";
        return false;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    // Find the "bins" object
    size_t pos = s.find("\"bins\"");
    if (pos == std::string::npos) {
        std::cerr << "[pi0_corr][WARN] 'bins' not found in " << path << "\n";
        return false;
    }
    size_t br = s.find('{', pos);
    if (br == std::string::npos) return false;
    int d=0; size_t i=br;
    for (; i<s.size(); ++i){
        if (s[i]=='{') d++;
        else if (s[i]=='}'){ d--; if(!d){ ++i; break; } }
    }
    std::string binsObj = s.substr(br, i-br);

    size_t kpos = 0;
    int nbins_loaded = 0;

    while (true) {
        // Key "(ix,iQ,it,ip)"
        size_t q1 = binsObj.find('"', kpos);
        if (q1 == std::string::npos) break;
        size_t q2 = binsObj.find('"', q1+1);
        if (q2 == std::string::npos) break;
        std::string key = binsObj.substr(q1+1, q2-q1-1);

        BinKey bk;
        if (!parse_tuple_key(key, bk)) {
            kpos = q2 + 1;
            continue;
        }

        // Value object { ... } for this bin
        size_t objS = binsObj.find('{', q2);
        if (objS == std::string::npos) break;
        int d2=0; size_t j=objS;
        for (; j<binsObj.size(); ++j){
            if (binsObj[j]=='{') d2++;
            else if (binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } }
        }
        std::string objPretty = binsObj.substr(objS, j-objS);
        std::string obj = squash_ws(objPretty); // strip whitespace

        auto findD = [&](const std::string& pat)->double{
            size_t p = obj.find(pat);
            if (p == std::string::npos) return 0.0;
            p = obj.find(':', p);
            if (p == std::string::npos) return 0.0;
            size_t a = p+1;
            size_t b = a;
            auto isnum = [](char c){
                return std::isdigit((unsigned char)c) || c=='-'||c=='+'||c=='.'||c=='e'||c=='E';
            };
            while (b<obj.size() && isnum(obj[b])) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch (...) { return 0.0; }
        };

        // Read both helicities regardless of spacing
        ContamBin cb;
        cb.c_plus      = findD("\"contamination\":{\"+1\":{\"value\"");
        cb.c_plus_err  = findD("\"contamination\":{\"+1\":{\"err\"");
        cb.c_minus     = findD("\"contamination\":{\"-1\":{\"value\"");
        cb.c_minus_err = findD("\"contamination\":{\"-1\":{\"err\"");

        // If any of these were not found due to slightly different nesting, try a second pattern:
        if (cb.c_plus==0.0 && cb.c_minus==0.0) {
            // Fallback if writer had keys in a different order but same names
            cb.c_plus      = findD("\"+1\":{\"value\"");
            cb.c_plus_err  = findD("\"+1\":{\"err\"");
            cb.c_minus     = findD("\"-1\":{\"value\"");
            cb.c_minus_err = findD("\"-1\":{\"err\"");
        }

        out[bk] = cb;
        ++nbins_loaded;

        kpos = j;
    }

    std::cout << "[pi0_corr] Parsed " << nbins_loaded << " bins from " << path << "\n";
    return nbins_loaded > 0;
}
} // namespace

// ------------- plotting -------------
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

static void plot_corrected_counts_for_period(
    const std::string& period,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<BinKey, CorrBin>& corr,
    const std::string& out_dir_plots)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(out_dir_plots, ec);

    static const auto PHI_DEG = phiCentersDeg();

    for (int ix=0; ix<(int)xB_bins.size(); ++ix){
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_pi0corr_"<<period<<"_xB"<<ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        // title
        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

        // grid
        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // title text
        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.36);
        head.DrawLatex(0.5, 0.55, Form("#pi^{0}-corrected counts  %s   x_{B} #in [%.2g, %.2g]",
                                       period.c_str(), xb.first, xb.second));

        // find global max for Y-scale per canvas (robust)
        double ymax_est = 1.0;
        for (int r=0; r<nrows; ++r){
            const int it_glob = findIndex(t_slice[r], t_bins);
            if (it_glob<0) continue;
            for (int cc=0; cc<ncols; ++cc){
                const int iQ_glob = findIndex(Q2_slice[cc], Q2_bins);
                if (iQ_glob<0) continue;
                for (int ip=0; ip<N_PHI_BINS; ++ip){
                    auto it = corr.find(std::make_tuple(ix,iQ_glob,it_glob,ip));
                    if (it==corr.end()) continue;
                    ymax_est = std::max(ymax_est, std::max(it->second.Np, it->second.Nm));
                }
            }
        }
        const double YMAX = std::max(5.0, ymax_est*1.20);

        for (int r=0; r<nrows; ++r){
            const int it_glob = findIndex(t_slice[r], t_bins);
            if (it_glob<0) continue;
            for (int cc=0; cc<ncols; ++cc){
                const int iQ_glob = findIndex(Q2_slice[cc], Q2_bins);
                if (iQ_glob<0) continue;

                pGrid->cd(r*ncols + cc + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, YMAX);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("N^{corr}");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.048);
                frame->GetXaxis()->SetTitleOffset(1.25);
                frame->GetYaxis()->SetTitleOffset(1.35);

                // series
                std::vector<double> x, yP, yM;
                x.reserve(N_PHI_BINS); yP.reserve(N_PHI_BINS); yM.reserve(N_PHI_BINS);
                for (int ip=0; ip<N_PHI_BINS; ++ip){
                    auto it = corr.find(std::make_tuple(ix,iQ_glob,it_glob,ip));
                    if (it==corr.end()) continue;
                    x.push_back(PHI_DEG[ip]);
                    yP.push_back(std::max(0.0, it->second.Np));
                    yM.push_back(std::max(0.0, it->second.Nm));
                }

                if (!x.empty()){
                    TGraph* gP = new TGraph((int)x.size(), x.data(), yP.data());
                    TGraph* gM = new TGraph((int)x.size(), x.data(), yM.data());
                    gP->SetLineWidth(2); gP->SetMarkerStyle(20); gP->Draw("LP SAME");
                    gM->SetLineWidth(2); gM->SetMarkerStyle(24); gM->Draw("LP SAME");

                    TLegend* leg = new TLegend(0.50, 0.74, 0.90, 0.92);
                    leg->SetTextFont(42); leg->SetTextSize(0.040);
                    leg->AddEntry(gP, "N^{corr}(+1)", "lp");
                    leg->AddEntry(gM, "N^{corr}(-1)", "lp");
                    leg->Draw();
                }

                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11); lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[cc].first, Q2_slice[cc].second,
                         t_slice[r].first,   t_slice[r].second));
            }
        }

        std::ostringstream fout;
        fout<<out_dir_plots<<"/plot_pi0corr_"<<period<<"_xB_"<<ix<<".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

} // namespace

// ======================================================================
// Public driver
// ======================================================================
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json_path,
    const std::string& contamination_dir,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    GroupCounts groups;
    if (!load_total_counts(total_counts_json_path, groups)) {
        std::cerr<<"[pi0_corr][ERROR] Failed to load total_counts json.\n"; return;
    }

    std::map<std::string, std::map<BinKey, CorrBin>> allPeriod;

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);

        auto itG = groups.find(runTag);
        if (itG == groups.end()) {
            std::cerr<<"[pi0_corr][WARN] total_counts has no group '"<<runTag<<"'.\n";
        }

        // load contamination for this period
        std::map<BinKey, ContamBin> contam;
        const fs::path contam_path = fs::path(contamination_dir)/("contamination_"+period+".json");
        if (!load_contam_for_period(contam_path.string(), contam)) {
            std::cerr<<"[pi0_corr][WARN] No contamination file for "<<period<<". Assuming c=0.\n";
        } else {
            std::cout<<"[pi0_corr] Loaded contamination for "<<period<<" ("<<contam.size()<<" bins)\n";
        }

        std::map<BinKey, CorrBin> corr;
        if (itG != groups.end()) {
            const auto& counts = itG->second;
            for (const auto& kv : counts){
                const BinKey bk = kv.first;
                const HelCounts& hc = kv.second;

                ContamBin cb;
                auto itC = contam.find(bk);
                if (itC!=contam.end()) cb = itC->second;

                // apply correction
                const double Np_corr = double(hc.plus ) * (1.0 - std::max(0.0, std::min(0.95, cb.c_plus )));
                const double Nm_corr = double(hc.minus) * (1.0 - std::max(0.0, std::min(0.95, cb.c_minus)));
                CorrBin out; out.Np = Np_corr; out.Nm = Nm_corr; out.Ntot = Np_corr + Nm_corr;
                corr[bk] = out;
            }
        }

        // write per-period JSON
        const fs::path out_json = fs::path(out_root_dir)/"jsons"/("pi0_corrected_counts_"+period+".json");
        std::error_code ec; fs::create_directories(out_json.parent_path(), ec);
        write_period_corrected_json(out_json.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, corr);
        std::cout<<"[pi0_corr] Wrote \""<<out_json.string()<<"\"\n";

        // plots
        const fs::path plots_dir = fs::path(out_root_dir)/"pi0_corrected_plots"/periodToRunTagKey(period);
        plot_corrected_counts_for_period(period, binning_scheme, xB_bins, Q2_bins, t_bins, corr, plots_dir.string());

        allPeriod[period] = std::move(corr);
    }

    // combined file (flat per-period)
    {
        const fs::path out_comb = fs::path(out_root_dir)/"jsons"/"pi0_corrected_counts_all_periods.json";
        std::ofstream ofs(out_comb.string());
        if (ofs){
            ofs<<std::fixed<<std::setprecision(8);
            ofs<<"{\n  \"periods\": {\n";
            bool firstP=true;
            for (const auto& pkv : allPeriod){
                if (!firstP) ofs<<",\n"; firstP=false;
                ofs<<"    \""<<pkv.first<<"\": {\n      \"bins\": {\n";
                bool firstB=true;
                for (const auto& kv : pkv.second){
                    if (!firstB) ofs<<",\n"; firstB=false;
                    int ix,iQ,it,ip; std::tie(ix,iQ,it,ip)=kv.first;
                    const auto& cb = kv.second;
                    ofs<<"        \""<<key4s(ix,iQ,it,ip)<<"\": {"
                       <<"\"helicity\":{\"+1\":"<<cb.Np<<",\"-1\":"<<cb.Nm<<"},\"total\":"<<cb.Ntot<<"}";
                }
                ofs<<"\n      }\n    }";
            }
            ofs<<"\n  }\n}\n";
            std::cout<<"[pi0_corr] Wrote combined \""<<out_comb.string()<<"\"\n";
        }
    }
}