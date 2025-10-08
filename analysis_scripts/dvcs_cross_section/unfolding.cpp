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

// ---------- small helpers ----------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
    return s;
}
static std::string periodToRunTagKey(const std::string& period) {
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
    for (int i=0;i<(int)ranges.size();++i) if (ranges[i]==range) return i;
    return -1;
}
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

// ---------- total_counts.json loader (group->(ix,iQ,it,ip)->helicity counts) ----------
struct HelCounts { long long plus=0, minus=0; };
using GroupCounts = std::map<std::string, std::map<BinKey4, HelCounts>>;

static bool parse_tuple_key4(const std::string& s, BinKey4& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey4(ix,iQ,it,ip);
    return true;
}

static bool load_total_counts(const std::string& path, GroupCounts& outGroups) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr<<"[unf][ERROR] Cannot open total_counts JSON: "<<path<<"\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t gpos = s.find("\"groups\"");
    if (gpos==std::string::npos) { std::cerr<<"[unf][ERROR] 'groups' not found in total_counts.\n"; return false; }
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

        std::map<BinKey4, HelCounts> gmap;
        size_t bpos=0;
        while (true) {
            size_t bk1 = binsObj.find('"', bpos); if (bk1==std::string::npos) break;
            size_t bk2 = binsObj.find('"', bk1+1); if (bk2==std::string::npos) break;
            std::string key = binsObj.substr(bk1+1, bk2-bk1-1);
            BinKey4 bk;
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

// ---------- acceptance_<period>.json loader ----------
struct AccCell {
    std::vector<double> phi_deg, acc, acc_err;
};
using AccMap3 = std::map<std::tuple<int,int,int>, AccCell>; // (ix,iQ,it)

static bool parse_tuple_key3(const std::string& s, std::tuple<int,int,int>& out) {
    int ix,iQ,it;
    if (std::sscanf(s.c_str(),"(%d,%d,%d)",&ix,&iQ,&it)!=3) return false;
    out = std::make_tuple(ix,iQ,it);
    return true;
}

static bool load_acceptance_json(const std::string& path, AccMap3& out) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr<<"[unf][WARN] Cannot open acceptance JSON: "<<path<<"\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    // Locate "bins" object
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

// ---------- per-cell result ----------
struct UnfoldCell {
    std::vector<double> phi_deg;

    std::vector<double> yield_p;     // +1
    std::vector<double> yield_p_err;

    std::vector<double> yield_m;     // -1
    std::vector<double> yield_m_err;

    // (optional) store acceptance used for sanity/debug
    std::vector<double> acc, acc_err;
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
static void plot_cells_for_period(
    const std::string& period,
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

        std::ostringstream cname; cname<<"c_unf_"<<period<<"_xB"<<ix;
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
        tit << Form("Unfolded Yields  %s   x_{B} #in [%.2g, %.2g]",
                    period.c_str(), xb.first, xb.second);
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

                // frame (y autoscale from 0)
                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 0.0); // we'll autoscale from data
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

                // Build graphs (+1 blue circle, −1 red square)
                std::vector<double> x, yp, ymp, ym, ymm;
                x.reserve(N_PHI_BINS); yp.reserve(N_PHI_BINS); ymp.reserve(N_PHI_BINS);
                ym.reserve(N_PHI_BINS); ymm.reserve(N_PHI_BINS);

                double ymax=0.0;
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    x.push_back(PHI_DEG[ip]);
                    double vp = uc.yield_p[ip];
                    double vm = uc.yield_m[ip];
                    double ep = std::max(1e-9, uc.yield_p_err[ip]);
                    double em = std::max(1e-9, uc.yield_m_err[ip]);
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

                // annotate Q² and -t
                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
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
        fout << out_dir_plots << "/plot_unfolded_" << period << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ---------- main driver ----------
} // anon

void compute_and_plot_unfolding(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json_path,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // Load total counts
    GroupCounts counts;
    if (!load_total_counts(total_counts_json_path, counts)) {
        std::cerr<<"[unf][ERROR] Failed to load total_counts json.\n"; return;
    }

    const fs::path json_dir  = fs::path(out_root_dir)/"jsons";
    const fs::path plot_root = fs::path(out_root_dir)/"unfolding";
    std::error_code ec;
    fs::create_directories(json_dir, ec);

    auto getCountsGroup = [&](const std::string& runTag)->const std::map<BinKey4,HelCounts>*{
        auto it = counts.find(runTag);
        if (it==counts.end()) return nullptr;
        return &it->second;
    };

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);

        // acceptance JSON for this period
        const fs::path acc_path = fs::path(out_root_dir)/"jsons"/("acceptance_"+period+".json");
        AccMap3 accCells;
        if (!load_acceptance_json(acc_path.string(), accCells)) {
            std::cerr<<"[unf][WARN] Missing/invalid acceptance for "<<period<<" — skipping.\n";
            continue;
        }

        // counts group (individual)
        const auto* gmap = getCountsGroup(runTag);
        if (!gmap) {
            std::cerr<<"[unf][WARN] No counts group '"<<runTag<<"' in total_counts.json — skipping "<<period<<"\n";
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
                // leave zeros; plotter will show empty
                outCells[std::make_tuple(ix,iQ,it)] = std::move(uc);
                continue;
            }
            const auto& ac = itAcc->second;

            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                // acceptance
                double A    = (ip<(int)ac.acc.size()    ? ac.acc[ip]     : 0.0);
                double sA   = (ip<(int)ac.acc_err.size()? ac.acc_err[ip] : 0.0);
                A = std::max(0.0, A);
                const double A_clamp = std::max(A, 1e-9); // protect division

                uc.acc[ip]     = A;
                uc.acc_err[ip] = sA;

                // counts
                BinKey4 k4(ix,iQ,it,ip);
                long long Np = 0, Nm = 0;
                auto itC = gmap->find(k4);
                if (itC != gmap->end()) { Np = itC->second.plus; Nm = itC->second.minus; }

                // unfolded (+)
                if (A > 0.0 && Np >= 0) {
                    double U = double(Np) / A_clamp;
                    // Var(U) ≈ (∂U/∂N)^2 Var(N) + (∂U/∂A)^2 Var(A), with Var(N)=N, Var(A)=sA^2
                    // ∂U/∂N = 1/A, ∂U/∂A = -N/A^2
                    double varU = (double(Np)/(A_clamp*A_clamp)) + (double(Np)*double(Np)/(A_clamp*A_clamp*A_clamp*A_clamp))* (sA*sA);
                    uc.yield_p[ip]     = U;
                    uc.yield_p_err[ip] = std::sqrt(std::max(0.0, varU));
                } else {
                    uc.yield_p[ip]     = 0.0;
                    uc.yield_p_err[ip] = 0.0;
                }

                // unfolded (−)
                if (A > 0.0 && Nm >= 0) {
                    double U = double(Nm) / A_clamp;
                    double varU = (double(Nm)/(A_clamp*A_clamp)) + (double(Nm)*double(Nm)/(A_clamp*A_clamp*A_clamp*A_clamp))* (sA*sA);
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
        const fs::path outJ = fs::path(out_root_dir)/"jsons"/("unfolded_"+period+".json");
        write_unfolded_json(outJ.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, outCells);
        std::cout<<"[unf] Wrote unfolded JSON: "<<outJ.string()<<"\n";

        // Plots
        const fs::path outPlots = fs::path(out_root_dir)/"unfolding"/runTag;
        std::error_code ec2; fs::create_directories(outPlots, ec2);
        plot_cells_for_period(period, binning_scheme, xB_bins, Q2_bins, t_bins, outCells, outPlots.string());
    }

    std::cout<<"[unf] Unfolding complete.\n";
}