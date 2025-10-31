// unfolding.cpp — uses pi0_corrected_counts_all_groups.json as the canonical source
// Reads corrected helicity counts (value, err) and unfolds by acceptance.
// Acceptance files are expected at: <out_root_dir>/jsons/acceptance_<DVCS_PERIOD>.json
// where <DVCS_PERIOD> is an exact DVCS_* name (e.g. DVCS_Sp18_inb) mapped from the short group key.
//
// Inputs (key ones):
//   - total_counts_json_path  -> path to pi0_corrected_counts_all_groups.json
//   - periods                 -> vector<string> of GROUP KEYS to process exactly as named
//                                (short keys: sp18_inb, sp18_out, fa18_inb_supp, fa18_inb, fa18_out, sp19_inb)
//   - binning_scheme          -> same Binning vector you already use elsewhere
//   - out_root_dir            -> root output dir containing jsons/ and plot dirs
//
// Outputs:
//   - <out_root_dir>/jsons/unfolded_<GROUP>.json
//   - <out_root_dir>/unfolding/<GROUP>/plot_unfolded_<GROUP>_xB_<ix>.png
//
// Notes:
//   - Uses helicity count uncertainties from the pi0-corrected JSON (err fields) instead of Poisson sqrt(N).
//   - Variance propagation: U = N / A, Var(U) ~ (1/A)^2 Var(N) + (N/A^2)^2 Var(A).
//   - Period names are NOT transformed when accessing the corrected master; they must match its group keys exactly.
//   - Acceptance filenames are looked up through a fixed map from short key -> DVCS_* key (exact filenames you already produced).

#include "unfolding.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>

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
        const int rf = 42; // Helvetica
        gStyle->SetTitleFont(rf, "XYZ");
        gStyle->SetLabelFont(rf, "XYZ");
        gStyle->SetTextFont(rf);
    }
} _style_bootstrap;

using BinKey4 = std::tuple<int,int,int,int>; // (ix,iQ,it,ip)

// ---------- tiny helpers ----------
static inline std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
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
static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
}
static inline void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize){
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}
static inline bool parse_tuple_key4(const std::string& s, BinKey4& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey4(ix,iQ,it,ip);
    return true;
}
static inline bool parse_tuple_key3(const std::string& s, std::tuple<int,int,int>& out) {
    int ix,iQ,it;
    if (std::sscanf(s.c_str(),"(%d,%d,%d)",&ix,&iQ,&it)!=3) return false;
    out = std::make_tuple(ix,iQ,it);
    return true;
}

// ---------- acceptance_<DVCS_PERIOD>.json loader ----------
struct AccCell {
    std::vector<double> phi_deg, acc, acc_err;
};
using AccMap3 = std::map<std::tuple<int,int,int>, AccCell>; // (ix,iQ,it)

static bool load_acceptance_json(const std::string& path, AccMap3& out) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr<<"[unf][WARN] Cannot open acceptance JSON: "<<path<<"\n";
        return false;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t bpos = s.find("\"bins\"");
    if (bpos==std::string::npos) return false;
    size_t br = s.find('{', bpos); if (br==std::string::npos) return false;
    int d=0; size_t i=br; for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string binsObj = s.substr(br, i-br);

    auto parseArray = [&](const std::string& obj, const char* key)->std::vector<double>{
        std::vector<double> v;
        size_t p = obj.find(key); if (p==std::string::npos) return v;
        p = obj.find('[', p); if (p==std::string::npos) return v;
        size_t q = obj.find(']', p); if (q==std::string::npos) return v;
        std::string arr = obj.substr(p+1, q-p-1);
        std::stringstream ss(arr);
        while (ss.good()){
            std::string tok; std::getline(ss, tok, ',');
            tok.erase(std::remove_if(tok.begin(), tok.end(), ::isspace), tok.end());
            if (tok.empty()) continue;
            try { v.push_back(std::stod(tok)); } catch(...) {}
        }
        return v;
    };

    size_t kpos=0;
    while (true) {
        size_t q1 = binsObj.find('"', kpos); if (q1==std::string::npos) break;
        size_t q2 = binsObj.find('"', q1+1); if (q2==std::string::npos) break;
        std::string key = binsObj.substr(q1+1, q2-q1-1);
        std::tuple<int,int,int> k3;
        if (!parse_tuple_key3(key, k3)) { kpos = q2+1; continue; }

        size_t objS = binsObj.find('{', q2); if (objS==std::string::npos) break;
        int d2=0; size_t j=objS; for (; j<binsObj.size(); ++j){ if(binsObj[j]=='{') d2++; else if(binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj = binsObj.substr(objS, j-objS);

        AccCell cell;
        cell.phi_deg = parseArray(obj, "\"phi\":[");
        cell.acc     = parseArray(obj, "\"acc\":[");
        cell.acc_err = parseArray(obj, "\"acc_err\":[");
        if (!cell.phi_deg.empty() && cell.acc.size()==cell.phi_deg.size() && cell.acc_err.size()==cell.phi_deg.size())
            out[k3] = std::move(cell);

        kpos = j;
    }
    return !out.empty();
}

// ---------- corrected-counts master loader (pi0_corrected_counts_all_groups.json) ----------
struct HelVals {
    double plus   = 0.0;
    double minus  = 0.0;
    double eplus  = 0.0;  // sigma on plus
    double eminus = 0.0;  // sigma on minus
};
using GroupHelMap = std::map<std::string, std::map<BinKey4, HelVals>>;

static bool load_pi0_corrected_master(const std::string& path, GroupHelMap& outGroups) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr<<"[unf][ERROR] Cannot open pi0_corrected_counts_all_groups JSON: "<<path<<"\n";
        return false;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    // Find "groups" object
    size_t gpos = s.find("\"groups\"");
    if (gpos==std::string::npos) { std::cerr<<"[unf][ERROR] 'groups' not found in corrected master.\n"; return false; }
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
        std::string gObj = groupsObj.substr(objS, j-objS);

        size_t bpos = gObj.find("\"bins\"");
        if (bpos==std::string::npos) { kpos=j; continue; }
        size_t br = gObj.find('{', bpos); if (br==std::string::npos) { kpos=j; continue; }
        int d3=0; size_t m=br; for (; m<gObj.size(); ++m){ if(gObj[m]=='{') d3++; else if(gObj[m]=='}'){ d3--; if(!d3){ ++m; break; } } }
        std::string binsObj = gObj.substr(br, m-br);

        std::map<BinKey4, HelVals> gmap;

        size_t bkey=0;
        while (true) {
            size_t bk1 = binsObj.find('"', bkey); if (bk1==std::string::npos) break;
            size_t bk2 = binsObj.find('"', bk1+1); if (bk2==std::string::npos) break;
            std::string key = binsObj.substr(bk1+1, bk2-bk1-1);
            BinKey4 bk;
            if (!parse_tuple_key4(key, bk)) { bkey=bk2+1; continue; }

            size_t valS = binsObj.find('{', bk2); if (valS==std::string::npos) break;
            int d4=0; size_t jj=valS;
            for (; jj<binsObj.size(); ++jj){ if(binsObj[jj]=='{') d4++; else if(binsObj[jj]=='}'){ d4--; if(!d4){ ++jj; break; } } }
            // NOTE: fix typo: binsObjs -> binsObj
        }
        // The above loop had a small typo; correct and continue parsing:
    }

    // Re-parse bins properly (fixed block)
    {
        // Rewind and parse again cleanly
        size_t kpos2=0;
        outGroups.clear();
        while (true) {
            size_t q1 = groupsObj.find('"', kpos2); if (q1==std::string::npos) break;
            size_t q2 = groupsObj.find('"', q1+1); if (q2==std::string::npos) break;
            std::string gname = groupsObj.substr(q1+1, q2-q1-1);

            size_t objS = groupsObj.find('{', q2); if (objS==std::string::npos) break;
            int d2=0; size_t j=objS;
            for (; j<groupsObj.size(); ++j){ if(groupsObj[j]=='{') d2++; else if(groupsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
            std::string gObj = groupsObj.substr(objS, j-objS);

            size_t bpos2 = gObj.find("\"bins\"");
            if (bpos2==std::string::npos) { kpos2=j; continue; }
            size_t br2 = gObj.find('{', bpos2); if (br2==std::string::npos) { kpos2=j; continue; }
            int d3=0; size_t m2=br2; for (; m2<gObj.size(); ++m2){ if(gObj[m2]=='{') d3++; else if(gObj[m2]=='}'){ d3--; if(!d3){ ++m2; break; } } }
            std::string binsObj = gObj.substr(br2, m2-br2);

            std::map<BinKey4, HelVals> gmap;

            auto extract_block = [&](const std::string& src, const char* label)->std::string {
                size_t p = src.find(label);
                if (p == std::string::npos) return std::string();
                size_t br3 = src.find('{', p);
                if (br3 == std::string::npos) return std::string();
                int dd=0; size_t ii=br3; for (; ii<src.size(); ++ii){ if(src[ii]=='{') dd++; else if(src[ii]=='}'){ dd--; if(!dd){ ++ii; break; } } }
                return src.substr(br3, ii-br3);
            };
            auto find_num = [&](const std::string& src, const char* keyname)->double {
                size_t p = src.find(keyname); if (p==std::string::npos) return 0.0;
                p = src.find(':', p); if (p==std::string::npos) return 0.0;
                size_t a=p+1; while (a<src.size() && std::isspace((unsigned char)src[a])) ++a;
                size_t b=a; while (b<src.size() && (std::isdigit((unsigned char)src[b]) || src[b]=='+' || src[b]=='-' || src[b]=='.' || src[b]=='e' || src[b]=='E')) ++b;
                try { return std::stod(src.substr(a,b-a)); } catch(...) { return 0.0; }
            };

            size_t bkey2=0;
            while (true) {
                size_t bk1 = binsObj.find('"', bkey2); if (bk1==std::string::npos) break;
                size_t bk2 = binsObj.find('"', bk1+1); if (bk2==std::string::npos) break;
                std::string key = binsObj.substr(bk1+1, bk2-bk1-1);
                BinKey4 bk;
                if (!parse_tuple_key4(key, bk)) { bkey2=bk2+1; continue; }

                size_t valS = binsObj.find('{', bk2); if (valS==std::string::npos) break;
                int d4=0; size_t jj=valS;
                for (; jj<binsObj.size(); ++jj){ if(binsObj[jj]=='{') d4++; else if(binsObj[jj]=='}'){ d4--; if(!d4){ ++jj; break; } } }
                std::string obj = binsObj.substr(valS, jj-valS);

                size_t hpos = obj.find("\"helicity\"");
                if (hpos == std::string::npos) { bkey2 = jj; continue; }
                size_t hbr = obj.find('{', hpos); if (hbr == std::string::npos) { bkey2 = jj; continue; }
                int dh=0; size_t hk=hbr; for (; hk<obj.size(); ++hk){ if(obj[hk]=='{') dh++; else if(obj[hk]=='}'){ dh--; if(!dh){ ++hk; break; } } }
                std::string helObj = obj.substr(hbr, hk-hbr);

                std::string plusBlk  = extract_block(helObj, "\"+1\"");
                std::string minusBlk = extract_block(helObj, "\"-1\"");

                HelVals hv;
                hv.plus   = find_num(plusBlk,  "\"value\"");
                hv.eplus  = find_num(plusBlk,  "\"err\"");
                hv.minus  = find_num(minusBlk, "\"value\"");
                hv.eminus = find_num(minusBlk, "\"err\"");

                gmap[bk] = hv;
                bkey2 = jj;
            }

            outGroups[gname] = std::move(gmap);
            kpos2 = j;
        }
    }

    return !outGroups.empty();
}

// ---------- per-cell result ----------
struct UnfoldCell {
    std::vector<double> phi_deg;

    std::vector<double> yield_p;     // +1
    std::vector<double> yield_p_err;

    std::vector<double> yield_m;     // -1
    std::vector<double> yield_m_err;

    std::vector<double> acc, acc_err; // stored for sanity/debug
};

// ---------- JSON writer ----------
static void write_unfolded_json(
    const std::string& out_path,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, UnfoldCell>& cells)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[unf][ERROR] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : cells){
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
        const auto& c = kv.second;
        ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {";
        auto dumpA=[&](const char* name,const std::vector<double>& v){
            ofs<<"\""<<name<<"\":["; for (size_t i=0;i<v.size();++i){ if(i)ofs<<","; ofs<<v[i]; } ofs<<"],";
        };
        dumpA("phi", c.phi_deg);
        dumpA("yield_plus", c.yield_p);
        dumpA("yield_plus_err", c.yield_p_err);
        dumpA("yield_minus", c.yield_m);
        dumpA("yield_minus_err", c.yield_m_err);
        dumpA("acc", c.acc);
        ofs<<"\"acc_err\":["; for (size_t i=0;i<c.acc_err.size();++i){ if(i)ofs<<","; ofs<<c.acc_err[i]; } ofs<<"]";
        ofs<<"}";
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
static void plot_cells_for_group(
    const std::string& group,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, UnfoldCell>& cells,
    const std::string& out_dir_plots)
{
    using std::filesystem::create_directories;
    std::error_code ec;
    create_directories(out_dir_plots, ec);

    static const auto PHI_DEG = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_unf_"<<group<<"_xB"<<ix;
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
        tit << Form("Unfolded Yields  %s   x_{B} in (%.2g, %.2g)",
                    group.c_str(), xb.first, xb.second);
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
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.10);

                // frame (y autoscale from 0)
                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 0.0);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetLabelSize(0.0001);
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("Unfolded yield");
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
                const auto& uc = itCell->second;

                std::vector<double> x, yp, ymp, ym, ymm;
                x.reserve(N_PHI_BINS); yp.reserve(N_PHI_BINS); ymp.reserve(N_PHI_BINS);
                ym.reserve(N_PHI_BINS); ymm.reserve(N_PHI_BINS);

                double ymax=0.0;
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    x.push_back(PHI_DEG[ip]);
                    double vp = uc.yield_p[ip];
                    double vm = uc.yield_m[ip];
                    double ep = std::max(1e-12, uc.yield_p_err[ip]);
                    double em = std::max(1e-12, uc.yield_m_err[ip]);
                    yp.push_back(vp); ymp.push_back(ep);
                    ym.push_back(vm); ymm.push_back(em);
                    ymax = std::max(ymax, std::max(vp+ep, vm+em));
                }

                if (ymax <= 0.0) ymax = 1.0;
                frame->GetYaxis()->SetRangeUser(0.0, ymax*1.20);

                TGraphErrors* grP = new TGraphErrors(N_PHI_BINS, x.data(), yp.data(), nullptr, ymp.data());
                grP->SetMarkerStyle(20);
                grP->SetMarkerSize(1.0);
                grP->SetLineWidth(2);
                grP->SetLineColor(kBlue+1);
                grP->SetMarkerColor(kBlue+1);
                grP->Draw("P SAME");

                TGraphErrors* grM = new TGraphErrors(N_PHI_BINS, x.data(), ym.data(), nullptr, ymm.data());
                grM->SetMarkerStyle(25);
                grM->SetMarkerSize(1.0);
                grM->SetLineWidth(2);
                grM->SetLineColor(kRed+1);
                grM->SetMarkerColor(kRed+1);
                grM->Draw("P SAME");

                // annotate Q2 and -t
                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} in (%.2g, %.2g),   -t in (%.2g, %.2g)",
                         Q2_slice[ccol].first, Q2_slice[ccol].second,
                         t_slice[r].first,    t_slice[r].second));

                // legend
                TLegend* leg = new TLegend(0.50, 0.72, 0.90, 0.92);
                leg->SetBorderSize(1);
                leg->SetLineColor(kBlack);
                leg->SetFillColor(kWhite);
                leg->SetFillStyle(1001);
                leg->SetTextFont(42);
                leg->SetTextSize(0.040);
                leg->AddEntry(grP, "+ helicity", "lep");
                leg->AddEntry(grM, "- helicity", "lep");
                leg->Draw();
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/plot_unfolded_" << group << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ---------- DVCS acceptance filename mapping (short -> exact DVCS_* key) ----------
// These are the ONLY names we will look for on disk, matching what acceptance.cpp wrote.
static const std::map<std::string, std::string> kAccFileKey = {
    {"sp18_inb",      "DVCS_Sp18_inb"},
    {"sp18_out",      "DVCS_Sp18_out"},
    {"fa18_inb_supp", "DVCS_Fa18_inb_supp"},
    {"fa18_inb",      "DVCS_Fa18_inb"},
    {"fa18_out",      "DVCS_Fa18_out"},
    {"sp19_inb",      "DVCS_Sp19_inb"}
};

// ---------- main driver ----------
} // anon

void compute_and_plot_unfolding(
    const std::vector<std::string>& periods,           // exact short group keys to process
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json_path,         // path to pi0_corrected_counts_all_groups.json
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // Load corrected master (value, err per helicity) with short group keys
    GroupHelMap groups;
    if (!load_pi0_corrected_master(total_counts_json_path, groups)) {
        std::cerr<<"[unf][ERROR] Failed to load corrected master json.\n";
        return;
    }

    const fs::path json_dir  = fs::path(out_root_dir)/"jsons";
    const fs::path plot_root = fs::path(out_root_dir)/"unfolding";
    std::error_code ec;
    fs::create_directories(json_dir, ec);

    auto getGroup = [&](const std::string& key)->const std::map<BinKey4,HelVals>*{
        auto it = groups.find(key);
        if (it==groups.end()) return nullptr;
        return &it->second;
    };

    for (const auto& group : periods) {
        const auto* gmap = getGroup(group);
        if (!gmap) {
            std::cerr<<"[unf][WARN] No group '"<<group<<"' in corrected master — skipping\n";
            continue;
        }

        // acceptance JSON for this group (via explicit DVCS_* filename mapping)
        auto itMap = kAccFileKey.find(group);
        if (itMap == kAccFileKey.end()) {
            std::cerr<<"[unf][WARN] No DVCS acceptance mapping for group '"<<group<<"' — skipping.\n";
            continue;
        }
        const std::string dvcs_key = itMap->second; // e.g. "DVCS_Sp18_inb"
        const fs::path acc_path = fs::path(out_root_dir)/"jsons"/("acceptance_"+dvcs_key+".json");

        AccMap3 accCells;
        if (!load_acceptance_json(acc_path.string(), accCells)) {
            std::cerr<<"[unf][WARN] Missing/invalid acceptance for "<<group<<" at "<<acc_path.string()<<" — skipping.\n";
            continue;
        }

        std::map<std::tuple<int,int,int>, UnfoldCell> outCells;

        // Build all cells
        for (int ix=0; ix<(int)xB_bins.size(); ++ix)
        for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
        for (int it=0; it<(int)t_bins.size();  ++it) {
            UnfoldCell uc;
            uc.phi_deg = phiCentersDeg();
            uc.yield_p.assign(N_PHI_BINS, 0.0);
            uc.yield_p_err.assign(N_PHI_BINS, 0.0);
            uc.yield_m.assign(N_PHI_BINS, 0.0);
            uc.yield_m_err.assign(N_PHI_BINS, 0.0);
            uc.acc.assign(N_PHI_BINS, 0.0);
            uc.acc_err.assign(N_PHI_BINS, 0.0);

            // acceptance cell
            auto itAcc = accCells.find(std::make_tuple(ix,iQ,it));
            if (itAcc == accCells.end()) {
                outCells[std::make_tuple(ix,iQ,it)] = std::move(uc);
                continue;
            }
            const auto& ac = itAcc->second;

            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                // acceptance
                double A    = (ip<(int)ac.acc.size()    ? ac.acc[ip]     : 0.0);
                double sA   = (ip<(int)ac.acc_err.size()? ac.acc_err[ip] : 0.0);
                A = std::max(0.0, A);
                const double A_clamp = std::max(A, 1e-12); // protect division

                uc.acc[ip]     = A;
                uc.acc_err[ip] = sA;

                // corrected helicity counts (value, err)
                BinKey4 k4(ix,iQ,it,ip);
                auto itC = gmap->find(k4);

                double Np = 0.0, Nm = 0.0, sNp = 0.0, sNm = 0.0;
                if (itC != gmap->end()) {
                    Np  = std::max(0.0, itC->second.plus);
                    Nm  = std::max(0.0, itC->second.minus);
                    sNp = std::max(0.0, itC->second.eplus);
                    sNm = std::max(0.0, itC->second.eminus);
                }

                // Unfolded (+): U = N/A
                if (A > 0.0) {
                    double U   = Np / A_clamp;
                    double vN  = sNp*sNp;                          // variance of corrected N
                    double vA  = sA*sA;                            // variance of A
                    double varU = (vN/(A_clamp*A_clamp)) + ((Np*Np)/(A_clamp*A_clamp*A_clamp*A_clamp))*vA;
                    uc.yield_p[ip]     = U;
                    uc.yield_p_err[ip] = std::sqrt(std::max(0.0, varU));
                } else {
                    uc.yield_p[ip]     = 0.0;
                    uc.yield_p_err[ip] = 0.0;
                }

                // Unfolded (-): U = N/A
                if (A > 0.0) {
                    double U   = Nm / A_clamp;
                    double vN  = sNm*sNm;
                    double vA  = sA*sA;
                    double varU = (vN/(A_clamp*A_clamp)) + ((Nm*Nm)/(A_clamp*A_clamp*A_clamp*A_clamp))*vA;
                    uc.yield_m[ip]     = U;
                    uc.yield_m_err[ip] = std::sqrt(std::max(0.0, varU));
                } else {
                    uc.yield_m[ip]     = 0.0;
                    uc.yield_m_err[ip] = 0.0;
                }
            }

            outCells[std::make_tuple(ix,iQ,it)] = std::move(uc);
        }

        // JSON
        const fs::path outJ = fs::path(out_root_dir)/"jsons"/("unfolded_"+group+".json");
        write_unfolded_json(outJ.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, outCells);
        std::cout<<"[unf] Wrote unfolded JSON: "<<outJ.string()<<"\n";

        // Plots
        const fs::path outPlots = fs::path(out_root_dir)/"unfolding"/group;
        std::error_code ec2; fs::create_directories(outPlots, ec2);
        plot_cells_for_group(group, binning_scheme, xB_bins, Q2_bins, t_bins, outCells, outPlots.string());
    }

    std::cout<<"[unf] Unfolding complete.\n";
}