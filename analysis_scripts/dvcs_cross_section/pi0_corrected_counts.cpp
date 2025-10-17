// π0-corrected counts with errors, for periods AND combined groups.
// Uses total_counts.json ("groups": {...}) and
//  - per-period contamination_<period>.json (each has "bins": {...})
//  - combined pi0_contamination_combined.json ("periods": {Group: {bins:{...}}})
//
// Outputs one JSON per group + a master JSON, and error-bar plots
// (xB in canvas title; Q² and −t in subplot headers).

#include "pi0_corrected_counts.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <TString.h>

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
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)

// ────────── styles ──────────
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);
        gStyle->SetTitleFont(42,"XYZ");
        gStyle->SetLabelFont(42,"XYZ");
        gStyle->SetTextFont(42);
    }
} _style_bootstrap;

// ────────── helpers ──────────
static inline std::string toLower(std::string s){
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c){return std::tolower(c);});
    return s;
}
static inline std::string periodToRunTagKey(const std::string& period){
    auto pos = period.find('_');
    if (pos==std::string::npos || pos+1>=period.size()) return toLower(period);
    return toLower(period.substr(pos+1));
}
static inline std::vector<std::pair<double,double>>
uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which=='x') s.emplace(b.xBmin,b.xBmax);
        else if (which=='Q') s.emplace(b.Q2min,b.Q2max);
        else if (which=='t') s.emplace(b.tmin,b.tmax);
    }
    return {s.begin(), s.end()};
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
static inline std::string key4s(int ix,int iQ,int it,int ip){
    std::ostringstream os; os<<"("<<ix<<","<<iQ<<","<<it<<","<<ip<<")"; return os.str();
}

// ────────── structures ──────────
struct HelCounts { long long plus=0, minus=0; };
using GroupCounts = std::map<std::string, std::map<BinKey, HelCounts>>;

struct Contam { double cp=0.0, ep=0.0, cm=0.0, em=0.0; }; // contamination & errors
using ContamTable = std::map<BinKey, Contam>;

struct CorrBin {
    double Np=0.0, ep=0.0;  // corrected +1 and its σ
    double Nm=0.0, em=0.0;  // corrected −1 and its σ
    double Nt=0.0, et=0.0;  // total and its σ (quadrature)
};

// ────────── JSON parsing utilities (tolerant) ──────────
static bool parse_tuple_key(const std::string& s, BinKey& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey(ix,iQ,it,ip);
    return true;
}

// total_counts.json → groups -> (key)->helicity
static bool load_total_counts(const std::string& path, GroupCounts& out) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr<<"[pi0corr][ERROR] open "<<path<<" failed\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t gpos = s.find("\"groups\"");
    if (gpos==std::string::npos) { std::cerr<<"[pi0corr][ERROR] 'groups' not found in total_counts\n"; return false; }
    size_t br = s.find('{', gpos); if (br==std::string::npos) return false;

    int d=0; size_t i=br;
    for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string groupsObj = s.substr(br, i-br);

    size_t k=0;
    while (true) {
        size_t q1 = groupsObj.find('"', k); if (q1==std::string::npos) break;
        size_t q2 = groupsObj.find('"', q1+1); if (q2==std::string::npos) break;
        std::string gname = groupsObj.substr(q1+1, q2-q1-1);

        size_t binsS = groupsObj.find("\"bins\"", q2);
        if (binsS==std::string::npos) { k=q2+1; continue; }
        size_t bbr = groupsObj.find('{', binsS); if (bbr==std::string::npos){ k=q2+1; continue; }

        int d2=0; size_t j=bbr;
        for (; j<groupsObj.size(); ++j){ if(groupsObj[j]=='{') d2++; else if(groupsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string binsObj = groupsObj.substr(bbr, j-bbr);

        std::map<BinKey, HelCounts> table; size_t p=0;
        while (true) {
            size_t k1=binsObj.find('"', p); if (k1==std::string::npos) break;
            size_t k2=binsObj.find('"', k1+1); if (k2==std::string::npos) break;
            BinKey bk; if (!parse_tuple_key(binsObj.substr(k1+1,k2-k1-1), bk)){ p=k2+1; continue; }

            size_t vS = binsObj.find('{', k2); if (vS==std::string::npos) break;
            int d3=0; size_t m=vS;
            for (; m<binsObj.size(); ++m){ if(binsObj[m]=='{') d3++; else if(binsObj[m]=='}'){ d3--; if(!d3){ ++m; break;} } }
            std::string v = binsObj.substr(vS, m-vS);

            auto findLL=[&](const char* pat)->long long{
                size_t pos=v.find(pat); if(pos==std::string::npos) return 0;
                pos=v.find(':',pos); if(pos==std::string::npos) return 0;
                size_t a=pos+1; while(a<v.size() && isspace((unsigned char)v[a])) ++a;
                size_t b=a; while(b<v.size() && (isdigit((unsigned char)v[b])||v[b]=='-')) ++b;
                try { return std::stoll(v.substr(a,b-a)); } catch(...) { return 0; }
            };

            HelCounts hc;
            hc.plus  = findLL("\"+1\"");
            hc.minus = findLL("\"-1\"");
            table[bk] = hc;
            p = m;
        }
        out[gname] = std::move(table);
        k = j;
    }
    return !out.empty();
}

// per-period contamination file
static bool load_contam_period(const std::string& path, ContamTable& out) {
    std::ifstream ifs(path);
    if (!ifs) return false;
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t bpos = s.find("\"bins\"");
    if (bpos==std::string::npos) return false;
    size_t br = s.find('{', bpos); if (br==std::string::npos) return false;

    int d=0; size_t i=br;
    for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string binsObj = s.substr(br, i-br);

    size_t p=0; int nb=0;
    while (true) {
        size_t q1=binsObj.find('"', p); if (q1==std::string::npos) break;
        size_t q2=binsObj.find('"', q1+1); if (q2==std::string::npos) break;
        BinKey bk; if (!parse_tuple_key(binsObj.substr(q1+1,q2-q1-1), bk)) { p=q2+1; continue; }

        size_t vS=binsObj.find('{', q2); if (vS==std::string::npos) break;
        int d2=0; size_t j=vS;
        for (; j<binsObj.size(); ++j){ if(binsObj[j]=='{') d2++; else if(binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj=binsObj.substr(vS, j-vS);

        auto findD=[&](const std::string& pat)->double{
            size_t pos=obj.find(pat); if(pos==std::string::npos) return 0.0;
            pos=obj.find(':',pos); if(pos==std::string::npos) return 0.0;
            size_t a=pos+1, b=a;
            auto isnum=[](char c){return std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E';};
            while(b<obj.size() && isnum(obj[b])) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch(...) { return 0.0; }
        };

        Contam c;
        // prefer nested "contamination" schema; fallback to flat "+1"/"-1"
        c.cp = findD("\"contamination\":{\"+1\":{\"value\"");
        c.ep = findD("\"contamination\":{\"+1\":{\"err\"");
        c.cm = findD("\"contamination\":{\"-1\":{\"value\"");
        c.em = findD("\"contamination\":{\"-1\":{\"err\"");
        if (c.cp==0.0 && c.cm==0.0) {
            c.cp = findD("\"+1\":{\"value\"");
            c.ep = findD("\"+1\":{\"err\"");
            c.cm = findD("\"-1\":{\"value\"");
            c.em = findD("\"-1\":{\"err\"");
        }
        out[bk] = c;
        ++nb; p=j;
    }
    return !out.empty();
}

// combined contamination file: periods -> group -> bins
static bool load_contam_group_from_combined(const std::string& combined_path,
                                            const std::string& group,
                                            ContamTable& out)
{
    std::ifstream ifs(combined_path);
    if (!ifs) return false;
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    // find "periods" then this group, then its "bins"
    size_t P = s.find("\"periods\""); if (P==std::string::npos) return false;
    size_t Pbr = s.find('{', P); if (Pbr==std::string::npos) return false;

    // find group key
    std::string gq = std::string("\"")+group+std::string("\"");
    size_t G = s.find(gq, Pbr); if (G==std::string::npos) return false;

    size_t binsK = s.find("\"bins\"", G); if (binsK==std::string::npos) return false;
    size_t br = s.find('{', binsK); if (br==std::string::npos) return false;

    int d=0; size_t i=br;
    for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string binsObj = s.substr(br, i-br);

    size_t p=0; int nb=0;
    while (true) {
        size_t q1=binsObj.find('"', p); if (q1==std::string::npos) break;
        size_t q2=binsObj.find('"', q1+1); if (q2==std::string::npos) break;
        BinKey bk; if (!parse_tuple_key(binsObj.substr(q1+1,q2-q1-1), bk)) { p=q2+1; continue; }

        size_t vS=binsObj.find('{', q2); if (vS==std::string::npos) break;
        int d2=0; size_t j=vS;
        for (; j<binsObj.size(); ++j){ if(binsObj[j]=='{') d2++; else if(binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj=binsObj.substr(vS, j-vS);

        auto findD=[&](const std::string& pat)->double{
            size_t pos=obj.find(pat); if(pos==std::string::npos) return 0.0;
            pos=obj.find(':',pos); if(pos==std::string::npos) return 0.0;
            size_t a=pos+1, b=a;
            auto isnum=[](char c){return std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E';};
            while(b<obj.size() && isnum(obj[b])) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch(...) { return 0.0; }
        };

        Contam c;
        c.cp = findD("\"contamination\":{\"+1\":{\"value\"");
        c.ep = findD("\"contamination\":{\"+1\":{\"err\"");
        c.cm = findD("\"contamination\":{\"-1\":{\"value\"");
        c.em = findD("\"contamination\":{\"-1\":{\"err\"");
        if (c.cp==0.0 && c.cm==0.0) {
            c.cp = findD("\"+1\":{\"value\"");
            c.ep = findD("\"+1\":{\"err\"");
            c.cm = findD("\"-1\":{\"value\"");
            c.em = findD("\"-1\":{\"err\"");
        }
        out[bk]=c; ++nb; p=j;
    }
    return !out.empty();
}

// ────────── JSON writer ──────────
static void write_group_json(const std::string& out_path,
                             int nPhi,
                             const std::vector<std::pair<double,double>>& xB_bins,
                             const std::vector<std::pair<double,double>>& Q2_bins,
                             const std::vector<std::pair<double,double>>& t_bins,
                             const std::map<BinKey, CorrBin>& table)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[pi0corr][ERROR] cannot write "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi
       <<", \"xB_bins\": "<<xB_bins.size()
       <<", \"Q2_bins\": "<<Q2_bins.size()
       <<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it,ip; std::tie(ix,iQ,it,ip)=kv.first;
        const auto& b = kv.second;
        ofs<<"    \""<<key4s(ix,iQ,it,ip)<<"\": {"
           <<"\"helicity\":{"
           <<"\"+1\":{\"value\":"<<b.Np<<",\"err\":"<<b.ep<<"},"
           <<"\"-1\":{\"value\":"<<b.Nm<<",\"err\":"<<b.em<<"}"
           <<"},"
           <<"\"total\":{\"value\":"<<b.Nt<<",\"err\":"<<b.et<<"}"
           <<"}";
    }
    ofs<<"\n  }\n}\n";
}

// ────────── plotting ──────────
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

static void plot_group(
    const std::string& group,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<BinKey, CorrBin>& corrected,
    const std::map<BinKey, HelCounts>& raw_totals,
    const std::string& out_dir_plots)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(out_dir_plots, ec);

    const auto PHI = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];
        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)Q2_slice.size();
        const int ncols = (int)t_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;
        std::string cname = Form("c_pi0corr_%s_xB%d", group.c_str(), ix);
        TCanvas* c = new TCanvas(cname.c_str(), cname.c_str(), W, H);

        // header pad
        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.92, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
        // grid pad
        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.92);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd(); pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // canvas title with xB
        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.40);
        head.DrawLatex(0.5, 0.55, Form("#pi^{0}-corrected vs raw counts  %s   x_{B} #in [%.3g, %.3g]",
                                       group.c_str(), xb.first, xb.second));

        double ymax = 1.0;
        for (const auto& kv : corrected)
            ymax = std::max(ymax, kv.second.Nt + kv.second.et);
        ymax = std::max(5.0, 1.20*ymax);

        for (int r = 0; r < nrows; ++r) {
            int iQ=findIndex(Q2_slice[r], Q2_bins); if (iQ<0) continue;
            for (int cc=0; cc<ncols; ++cc) {
                int it=findIndex(t_slice[cc], t_bins); if (it<0) continue;

                pGrid->cd(r*ncols + cc + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.18);
                gPad->SetBottomMargin(0.16);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.08);

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Counts");
                frame->GetXaxis()->CenterTitle(); frame->GetYaxis()->CenterTitle();

                std::vector<double> X, Yp_raw, Ym_raw, Yp_corr, Ym_corr, eX, eYp_corr, eYm_corr;
                X.reserve(N_PHI_BINS); eX.assign(N_PHI_BINS, 0.0);

                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    BinKey bk(ix, iQ, it, ip);
                    auto it_raw = raw_totals.find(bk);
                    auto it_corr = corrected.find(bk);
                    if (it_raw == raw_totals.end() && it_corr == corrected.end()) continue;

                    X.push_back(PHI[ip]);
                    double Np_raw = it_raw != raw_totals.end() ? it_raw->second.plus : 0.0;
                    double Nm_raw = it_raw != raw_totals.end() ? it_raw->second.minus : 0.0;
                    double Np_corr = it_corr != corrected.end() ? it_corr->second.Np : 0.0;
                    double Nm_corr = it_corr != corrected.end() ? it_corr->second.Nm : 0.0;
                    double ep_corr = it_corr != corrected.end() ? it_corr->second.ep : 0.0;
                    double em_corr = it_corr != corrected.end() ? it_corr->second.em : 0.0;

                    Yp_raw.push_back(Np_raw);
                    Ym_raw.push_back(Nm_raw);
                    Yp_corr.push_back(Np_corr);
                    Ym_corr.push_back(Nm_corr);
                    eYp_corr.push_back(ep_corr);
                    eYm_corr.push_back(em_corr);
                }

                // --- make graphs ---
                TGraphErrors* gP_corr = new TGraphErrors((int)X.size(), X.data(), Yp_corr.data(), eX.data(), eYp_corr.data());
                TGraphErrors* gM_corr = new TGraphErrors((int)X.size(), X.data(), Ym_corr.data(), eX.data(), eYm_corr.data());
                TGraph* gP_raw = new TGraph((int)X.size(), X.data(), Yp_raw.data());
                TGraph* gM_raw = new TGraph((int)X.size(), X.data(), Ym_raw.data());

                // styles
                gP_corr->SetLineColor(kRed);  gP_corr->SetMarkerColor(kRed);
                gM_corr->SetLineColor(kBlue); gM_corr->SetMarkerColor(kBlue);
                gP_raw->SetLineColor(kRed);   gP_raw->SetMarkerColor(kRed);
                gM_raw->SetLineColor(kBlue);  gM_raw->SetMarkerColor(kBlue);

                gP_corr->SetMarkerStyle(20); // filled circle
                gM_corr->SetMarkerStyle(21);
                gP_raw->SetMarkerStyle(24);  // open square
                gM_raw->SetMarkerStyle(25);

                gP_corr->SetLineWidth(2); gM_corr->SetLineWidth(2);
                gP_raw->SetLineWidth(2);  gM_raw->SetLineWidth(2);

                // draw
                gP_raw->Draw("P SAME");
                gM_raw->Draw("P SAME");
                gP_corr->Draw("P SAME");
                gM_corr->Draw("P SAME");

                if (r==0 && cc==0) {
                    TLegend* leg = new TLegend(0.52, 0.68, 0.94, 0.92);
                    leg->SetBorderSize(1); leg->SetFillStyle(0); leg->SetTextSize(0.037);
                    leg->AddEntry(gP_raw,  "Raw (+1)",   "p");
                    leg->AddEntry(gM_raw,  "Raw (-1)",   "p");
                    leg->AddEntry(gP_corr, "#pi^{0}-corr (+1)", "p");
                    leg->AddEntry(gM_corr, "#pi^{0}-corr (-1)", "p");
                    leg->Draw();
                }

                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.042); lab.SetTextAlign(13);
                lab.DrawLatex(0.14, 0.96, Form("Q^{2} #in [%.3g, %.3g],   -t #in [%.3g, %.3g]",
                    Q2_slice[r].first, Q2_slice[r].second, t_slice[cc].first, t_slice[cc].second));
            }
        }

        std::string fout = (std::filesystem::path(out_dir_plots) /
                           Form("plot_pi0corr_%s_xB_%d.png", group.c_str(), ix)).string();
        c->SaveAs(fout.c_str());
        delete c;
    }
}

// ────────── core per-group computation ──────────
static std::map<BinKey, CorrBin> build_corrected_for_group(
    const std::map<BinKey, HelCounts>& totals,
    const ContamTable& contam)
{
    std::map<BinKey, CorrBin> out;

    for (const auto& kv : totals) {
        const BinKey& bk = kv.first;
        const HelCounts& N = kv.second;

        Contam c{};
        if (auto itc = contam.find(bk); itc != contam.end()) c = itc->second;

        // clamp c into [0, 0.95] for safety
        const double cp = std::clamp(c.cp, 0.0, 0.95);
        const double cm = std::clamp(c.cm, 0.0, 0.95);
        const double ep = std::max(0.0, c.ep);
        const double em = std::max(0.0, c.em);

        // corrected values
        CorrBin b;
        // + helicity
        {
            const double Np = (double)N.plus;
            const double val = Np * (1.0 - cp);
            const double var = (1.0 - cp)*(1.0 - cp) * Np      // Var(N) with Poisson
                             + (Np*Np) * (ep*ep);              // (d/dc)^2 Var(c)
            b.Np = std::max(0.0, val);
            b.ep = std::sqrt(std::max(0.0, var));
        }
        // - helicity
        {
            const double Nm = (double)N.minus;
            const double val = Nm * (1.0 - cm);
            const double var = (1.0 - cm)*(1.0 - cm) * Nm
                             + (Nm*Nm) * (em*em);
            b.Nm = std::max(0.0, val);
            b.em = std::sqrt(std::max(0.0, var));
        }
        b.Nt = b.Np + b.Nm;
        b.et = std::sqrt(b.ep*b.ep + b.em*b.em); // assume independence

        out[bk] = b;
    }
    return out;
}

} // namespace

// ────────────────────────────────────────────────────────────────────
// Public driver
// ────────────────────────────────────────────────────────────────────
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,
    const std::string& contamination_dir,
    const std::string& contamination_combined,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    // bin meta
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // 1) Load all DVCS totals (groups includes: periods + Spring2018/Fall2018/10.6_GeV)
    GroupCounts groups;
    if (!load_total_counts(total_counts_json, groups)) {
        std::cerr<<"[pi0corr][ERROR] could not read total_counts\n";
        return;
    }

    // 2) Build the list of groups we will output:
    //    - all requested periods (by runTag key)
    //    - any combined groups present in total_counts.json (e.g., Spring2018,...)
    std::vector<std::string> group_names;
    group_names.reserve(periods.size()+8);
    for (const auto& p : periods) group_names.push_back(periodToRunTagKey(p)); // runTag keys for per-period totals
    for (const auto& gkv : groups) {
        const std::string& g = gkv.first;
        // If it's already in the list (period key), keep; else (combined like "Spring2018") add it.
        if (std::find(group_names.begin(), group_names.end(), g) == group_names.end())
            group_names.push_back(g);
    }

    // 3) For each group, get contamination:
    //    - per-period groups: read contamination_<DVCS_PeriodName>.json in contamination_dir
    //    - combined groups:   read from contamination_combined.json under "periods"[Group]
    std::map<std::string, std::map<BinKey, CorrBin>> all_out;

    const fs::path out_json_dir = fs::path(out_root_dir)/"jsons";
    const fs::path out_plot_dir = fs::path(out_root_dir)/"pi0_corrected_plots";
    std::error_code ec;
    fs::create_directories(out_json_dir, ec);
    fs::create_directories(out_plot_dir, ec);

    for (const std::string& group : group_names) {
        // find total counts for this group
        auto itG = groups.find(group);
        if (itG == groups.end()) {
            std::cerr<<"[pi0corr][WARN] totals missing for group '"<<group<<"' — skipped\n";
            continue;
        }

        // contamination table for this group
        ContamTable C;

        // detect if group is one of the per-period DVCS_NN style or a combined name
        bool is_period_like = (group.rfind("dvcs_",0)==0) || (group.find('_')!=std::string::npos); // crude
        if (is_period_like) {
            // Need the DVCS Period name, not runTag. Try to re-hydrate:
            // dvcs_<Tag>  → "DVCS_" + TitleCase(Tag)
            // Our contamination file naming is contamination_<DVCS_Period>.json (as produced earlier).
            // Simpler: search for a contamination file that exists for any of the spelled DVCS names
            // but we’ll reconstruct from runTag:
            std::string tag = group;
            // TitleCase after underscores
            if (!tag.empty()) tag[0] = (char)std::toupper((unsigned char)tag[0]);
            for (size_t i=0;i+1<tag.size();++i) if (tag[i]=='_' && i+1<tag.size()) tag[i+1]=(char)std::toupper((unsigned char)tag[i+1]);
            std::string dvcsName = "DVCS_" + tag;

            fs::path cf = fs::path(contamination_dir)/("contamination_"+dvcsName+".json");
            if (!load_contam_period(cf.string(), C)) {
                // also allow raw period name (some trees use exact period strings already)
                fs::path alt = fs::path(contamination_dir)/("contamination_"+group+".json");
                if (!load_contam_period(alt.string(), C)) {
                    std::cerr<<"[pi0corr][WARN] contamination missing for "<<group<<" (looked for "<<cf<<") — assuming c=0\n";
                }
            }
        } else {
            // combined group: read from the combined contamination file
            if (!load_contam_group_from_combined(contamination_combined, group, C)) {
                std::cerr<<"[pi0corr][WARN] combined contamination missing for "<<group<<" — assuming c=0\n";
            }
        }

        // compute corrected + errors
        auto corr = build_corrected_for_group(itG->second, C);
        all_out[group] = corr;

        // write per-group JSON
        const fs::path out_group_json = out_json_dir / ("pi0_corrected_counts_"+group+".json");
        write_group_json(out_group_json.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, corr);
        std::cout<<"[pi0corr] wrote "<<out_group_json.string()<<"\n";

        // plots with error bars
        const fs::path plot_dir = out_plot_dir / group;
        plot_group(group, binning_scheme, xB_bins, Q2_bins, t_bins, corr, itG->second, plot_dir.string());
    }

    // master combined JSON (all groups in one file)
    {
        const fs::path out_master = out_json_dir / "pi0_corrected_counts_all_groups.json";
        std::ofstream ofs(out_master.string());
        if (!ofs) { std::cerr<<"[pi0corr][ERROR] cannot write "<<out_master<<"\n"; return; }
        ofs<<std::fixed<<std::setprecision(8);
        ofs<<"{\n  \"groups\": {\n";
        bool firstG=true;
        for (const auto& gkv : all_out) {
            if (!firstG) ofs<<",\n"; firstG=false;
            ofs<<"    \""<<gkv.first<<"\": {\n      \"bins\": {\n";
            bool firstB=true;
            for (const auto& kv : gkv.second) {
                if (!firstB) ofs<<",\n"; firstB=false;
                int ix,iQ,it,ip; std::tie(ix,iQ,it,ip)=kv.first;
                const auto& b = kv.second;
                ofs<<"        \""<<key4s(ix,iQ,it,ip)<<"\": {"
                   <<"\"helicity\":{"
                   <<"\"+1\":{\"value\":"<<b.Np<<",\"err\":"<<b.ep<<"},"
                   <<"\"-1\":{\"value\":"<<b.Nm<<",\"err\":"<<b.em<<"}"
                   <<"},"
                   <<"\"total\":{\"value\":"<<b.Nt<<",\"err\":"<<b.et<<"}"
                   <<"}";
            }
            ofs<<"\n      }\n    }";
        }
        ofs<<"\n  }\n}\n";
        std::cout<<"[pi0corr] wrote "<<out_master.string()<<"\n";
    }
}