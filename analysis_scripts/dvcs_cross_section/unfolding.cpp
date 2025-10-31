// unfolding.cpp - uses pi0_corrected_counts_all_groups.json as the canonical source
// Reads corrected helicity counts (value, err) and unfolds by acceptance.
// Acceptance files are expected at: <out_root_dir>/jsons/acceptance_<PERIOD>.json
//
// Per-period outputs:
//   - <out_root_dir>/jsons/unfolded_<PERIOD>.json
//   - <out_root_dir>/unfolding/<PERIOD>/plot_unfolded_<PERIOD>_xB_<ix>.png
//
// Combined outputs (this file also produces):
//   - <out_root_dir>/jsons/unfolded_Fa18.json
//   - <out_root_dir>/jsons/unfolded_Sp18.json
//   - <out_root_dir>/jsons/unfolded_10p6GeV.json
//   - Plots to corresponding directories under <out_root_dir>/unfolding/.
//
// Notes:
//   - Unfolding: U = N / A with error propagation from both N and A.
//   - Combination: inverse-variance weighting on the unfolded U across requested periods.
//   - PERIOD names are DVCS_* tokens for acceptance files; pi0-corrected master is keyed by runTags.

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

static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

static std::string periodToRunTagKey(const std::string& period) {
    // "DVCS_Fa18_inb" -> "fa18_inb"
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}

// ---------- helpers for bin edges and indices ----------
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
    for (int i = 0; i < (int)ranges.size(); ++i) if (ranges[i] == range) return i;
    return -1;
}

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
    for (int i = 0; i < N_PHI_BINS; ++i) d[i] = (i + 0.5) * step;
    return d;
}

static inline void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

static inline bool parse_tuple_key4(const std::string& s, BinKey4& out) {
    int ix, iQ, it, ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)", &ix, &iQ, &it, &ip) != 4) return false;
    out = BinKey4(ix, iQ, it, ip);
    return true;
}

static inline bool parse_tuple_key3(const std::string& s, std::tuple<int,int,int>& out) {
    int ix, iQ, it;
    if (std::sscanf(s.c_str(),"(%d,%d,%d)", &ix, &iQ, &it) != 3) return false;
    out = std::make_tuple(ix, iQ, it);
    return true;
}

// ---------- acceptance_<PERIOD>.json loader ----------
struct AccCell {
    std::vector<double> phi_deg, acc, acc_err;
};
using AccMap3 = std::map<std::tuple<int,int,int>, AccCell>; // (ix,iQ,it)

static bool load_acceptance_json(const std::string& path, AccMap3& out) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[unf][WARN] Cannot open acceptance JSON: " << path << "\n";
        return false;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t bpos = s.find("\"bins\"");
    if (bpos == std::string::npos) return false;
    size_t br = s.find('{', bpos); if (br == std::string::npos) return false;
    int d = 0; size_t i = br;
    for (; i < s.size(); ++i) {
        if (s[i] == '{') ++d;
        else if (s[i] == '}') { --d; if (!d) { ++i; break; } }
    }
    std::string binsObj = s.substr(br, i - br);

    auto parseArray = [&](const std::string& obj, const char* key)->std::vector<double>{
        std::vector<double> v;
        size_t p = obj.find(key); if (p == std::string::npos) return v;
        p = obj.find('[', p); if (p == std::string::npos) return v;
        size_t q = obj.find(']', p); if (q == std::string::npos) return v;
        std::string arr = obj.substr(p + 1, q - p - 1);
        std::stringstream ss(arr);
        while (ss.good()) {
            std::string tok; std::getline(ss, tok, ',');
            tok.erase(std::remove_if(tok.begin(), tok.end(), ::isspace), tok.end());
            if (tok.empty()) continue;
            try { v.push_back(std::stod(tok)); } catch (...) {}
        }
        return v;
    };

    size_t kpos = 0;
    while (true) {
        size_t q1 = binsObj.find('"', kpos); if (q1 == std::string::npos) break;
        size_t q2 = binsObj.find('"', q1 + 1); if (q2 == std::string::npos) break;
        std::string key = binsObj.substr(q1 + 1, q2 - q1 - 1);
        std::tuple<int,int,int> k3;
        if (!parse_tuple_key3(key, k3)) { kpos = q2 + 1; continue; }

        size_t objS = binsObj.find('{', q2); if (objS == std::string::npos) break;
        int d2 = 0; size_t j = objS;
        for (; j < binsObj.size(); ++j) {
            if (binsObj[j] == '{') ++d2;
            else if (binsObj[j] == '}') { --d2; if (!d2) { ++j; break; } }
        }
        std::string obj = binsObj.substr(objS, j - objS);

        AccCell cell;
        cell.phi_deg = parseArray(obj, "\"phi\":[");
        cell.acc     = parseArray(obj, "\"acc\":[");
        cell.acc_err = parseArray(obj, "\"acc_err\":[");
        if (!cell.phi_deg.empty() && cell.acc.size() == cell.phi_deg.size() && cell.acc_err.size() == cell.phi_deg.size())
            out[k3] = std::move(cell);

        kpos = j;
    }
    return !out.empty();
}

// ---------- corrected-counts master loader (pi0_corrected_counts_all_groups.json) ----------
struct HelVals {
    double plus   = 0.0;
    double minus  = 0.0;
    double eplus  = 0.0;
    double eminus = 0.0;
};
using GroupHelMap = std::map<std::string, std::map<BinKey4, HelVals>>;

static bool extract_braced_block_after_key(const std::string& s, const std::string& quotedKey, std::string& outObj) {
    size_t pkey = s.find(quotedKey);
    if (pkey == std::string::npos) return false;
    size_t pcolon = s.find(':', pkey + quotedKey.size());
    if (pcolon == std::string::npos) return false;
    size_t lbrace = s.find('{', pcolon);
    if (lbrace == std::string::npos) return false;
    int depth = 0;
    for (size_t i = lbrace; i < s.size(); ++i) {
        if (s[i] == '{') ++depth;
        else if (s[i] == '}') {
            --depth;
            if (depth == 0) {
                outObj = s.substr(lbrace, i - lbrace + 1);
                return true;
            }
        }
    }
    return false;
}

static std::string extract_object_member(const std::string& obj, const char* memberKey) {
    size_t p = obj.find(memberKey);
    if (p == std::string::npos) return std::string();
    size_t colon = obj.find(':', p);
    if (colon == std::string::npos) return std::string();
    size_t lbrace = obj.find('{', colon);
    if (lbrace == std::string::npos) return std::string();
    int depth = 0;
    for (size_t i = lbrace; i < obj.size(); ++i) {
        if (obj[i] == '{') ++depth;
        else if (obj[i] == '}') {
            --depth;
            if (depth == 0) {
                return obj.substr(lbrace, i - lbrace + 1);
            }
        }
    }
    return std::string();
}

static double find_numeric_field(const std::string& src, const char* keyname) {
    size_t p = src.find(keyname); if (p == std::string::npos) return 0.0;
    p = src.find(':', p); if (p == std::string::npos) return 0.0;
    size_t a = p + 1; while (a < src.size() && std::isspace((unsigned char)src[a])) ++a;
    size_t b = a; while (b < src.size() && (std::isdigit((unsigned char)src[b]) || src[b]=='+' || src[b]=='-' || src[b]=='.' || src[b]=='e' || src[b]=='E')) ++b;
    try { return std::stod(src.substr(a, b - a)); } catch(...) { return 0.0; }
}

static bool load_pi0_corrected_master(const std::string& path,
                                      const std::vector<std::string>& groupsWanted,
                                      GroupHelMap& outGroups) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[unf][ERROR] Cannot open pi0_corrected_counts_all_groups JSON: " << path << "\n";
        return false;
    }
    const std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    bool any = false;
    for (const auto& gname : groupsWanted) {
        const std::string quoted = "\"" + gname + "\"";
        std::string groupObj;
        if (!extract_braced_block_after_key(s, quoted, groupObj)) {
            std::cerr << "[unf][WARN] Group '" << gname << "' not found in corrected master.\n";
            continue;
        }

        std::string binsObj = extract_object_member(groupObj, "\"bins\"");
        if (binsObj.empty()) {
            std::cerr << "[unf][WARN] Group '" << gname << "' has no 'bins' object — skipping.\n";
            continue;
        }

        std::map<BinKey4, HelVals> gmap;

        size_t cursor = 0;
        while (true) {
            size_t q1 = binsObj.find('"', cursor); if (q1 == std::string::npos) break;
            size_t q2 = binsObj.find('"', q1 + 1);   if (q2 == std::string::npos) break;
            const std::string key = binsObj.substr(q1 + 1, q2 - q1 - 1);

            BinKey4 bk;
            if (!parse_tuple_key4(key, bk)) { cursor = q2 + 1; continue; }

            size_t lbrace = binsObj.find('{', q2);
            if (lbrace == std::string::npos) break;
            int depth = 0; size_t i = lbrace;
            for (; i < binsObj.size(); ++i) {
                if (binsObj[i] == '{') ++depth;
                else if (binsObj[i] == '}') { --depth; if (depth == 0) { ++i; break; } }
            }
            if (i >= binsObj.size()) break;
            const std::string entryObj = binsObj.substr(lbrace, i - lbrace);

            const std::string helObj = extract_object_member(entryObj, "\"helicity\"");
            if (helObj.empty()) { cursor = i; continue; }

            const std::string plusBlk  = extract_object_member(helObj, "\"+1\"");
            const std::string minusBlk = extract_object_member(helObj, "\"-1\"");

            HelVals hv;
            hv.plus   = find_numeric_field(plusBlk,  "\"value\"");
            hv.eplus  = find_numeric_field(plusBlk,  "\"err\"");
            hv.minus  = find_numeric_field(minusBlk, "\"value\"");
            hv.eminus = find_numeric_field(minusBlk, "\"err\"");

            gmap[bk] = hv;
            cursor = i;
        }

        if (!gmap.empty()) {
            outGroups[gname] = std::move(gmap);
            any = true;
        } else {
            std::cerr << "[unf][WARN] Group '" << gname << "' parsed but no bins found.\n";
        }
    }

    if (!any) {
        std::cerr << "[unf][ERROR] No groups parsed from corrected master.\n";
    }
    return any;
}

// ---------- per-cell result ----------
struct UnfoldCell {
    std::vector<double> phi_deg;
    std::vector<double> yield_p;     // +1
    std::vector<double> yield_p_err;
    std::vector<double> yield_m;     // -1
    std::vector<double> yield_m_err;
    std::vector<double> acc, acc_err; // optional for debug
};

// ---------- JSON writer ----------
static void write_unfolded_json(
    const std::string& out_path,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, UnfoldCell>& cells) {
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "[unf][ERROR] Cannot open " << out_path << "\n"; return; }
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": " << t_bins.size() << "},\n";
    ofs << "  \"bins\": {\n";
    bool first = true;
    for (const auto& kv : cells) {
        if (!first) ofs << ",\n"; first = false;
        int ix, iQ, it; std::tie(ix, iQ, it) = kv.first;
        const auto& c = kv.second;
        ofs << "    \"(" << ix << "," << iQ << "," << it << ")\": {";
        auto dumpA=[&](const char* name,const std::vector<double>& v){
            ofs << "\"" << name << "\":[";
            for (size_t i = 0; i < v.size(); ++i) { if (i) ofs << ","; ofs << v[i]; }
            ofs << "],";
        };
        dumpA("phi", c.phi_deg);
        dumpA("yield_plus", c.yield_p);
        dumpA("yield_plus_err", c.yield_p_err);
        dumpA("yield_minus", c.yield_m);
        dumpA("yield_minus_err", c.yield_m_err);
        dumpA("acc", c.acc);
        ofs << "\"acc_err\":[";
        for (size_t i = 0; i < c.acc_err.size(); ++i) { if (i) ofs << ","; ofs << c.acc_err[i]; }
        ofs << "]";
        ofs << "}";
    }
    ofs << "\n  }\n}\n";
}

// ---------- slice helpers for plotting ----------
static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list) {
    std::set<std::pair<double,double>> qs, ts;
    for (const auto& b : scheme) {
        if (std::make_pair(b.xBmin, b.xBmax) == xBrange) {
            qs.emplace(b.Q2min, b.Q2max);
            ts.emplace(b.tmin,  b.tmax);
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
    const std::string& out_dir_plots) {
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

        const int W = 280 * ncols + 160;
        const int H = 240 * nrows + 170;

        std::ostringstream cname; cname << "c_unf_" << group << "_xB" << ix;
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

                pGrid->cd(r * ncols + ccol + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.10);

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

                double ymax = 0.0;
                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    x.push_back(PHI_DEG[ip]);
                    double vp = uc.yield_p[ip];
                    double vm = uc.yield_m[ip];
                    double ep = std::max(1e-12, uc.yield_p_err[ip]);
                    double em = std::max(1e-12, uc.yield_m_err[ip]);
                    yp.push_back(vp); ymp.push_back(ep);
                    ym.push_back(vm); ymm.push_back(em);
                    ymax = std::max(ymax, std::max(vp + ep, vm + em));
                }

                if (ymax <= 0.0) ymax = 1.0;
                frame->GetYaxis()->SetRangeUser(0.0, ymax * 1.20);

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

                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} in (%.2g, %.2g),   -t in (%.2g, %.2g)",
                         Q2_slice[ccol].first, Q2_slice[ccol].second,
                         t_slice[r].first,    t_slice[r].second));

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

// ---------- inverse-variance combiner for unfolded cells ----------
static bool combine_cells_inverse_variance(
    const std::vector<const std::map<std::tuple<int,int,int>, UnfoldCell>*>& parts,
    std::map<std::tuple<int,int,int>, UnfoldCell>& combined) {

    if (parts.empty()) return false;

    // Find union of all keys
    std::set<std::tuple<int,int,int>> keys;
    for (const auto* m : parts) for (const auto& kv : *m) keys.insert(kv.first);

    const int nphi = N_PHI_BINS;
    for (const auto& key : keys) {
        UnfoldCell out;
        out.phi_deg = phiCentersDeg();
        out.yield_p.assign(nphi, 0.0);
        out.yield_p_err.assign(nphi, 0.0);
        out.yield_m.assign(nphi, 0.0);
        out.yield_m_err.assign(nphi, 0.0);
        out.acc.assign(nphi, 0.0);
        out.acc_err.assign(nphi, 0.0);

        for (int ip = 0; ip < nphi; ++ip) {
            double wsum_p = 0.0, wval_p = 0.0;
            double wsum_m = 0.0, wval_m = 0.0;

            double wsum_a = 0.0, wval_a = 0.0; // acceptance "debug" average
            for (const auto* m : parts) {
                auto it = m->find(key);
                if (it == m->end()) continue;
                const auto& c = it->second;

                if (ip < (int)c.yield_p.size()) {
                    double s = std::max(0.0, c.yield_p_err[ip]);
                    if (s > 0.0) { double w = 1.0 / (s*s); wsum_p += w; wval_p += w * c.yield_p[ip]; }
                }
                if (ip < (int)c.yield_m.size()) {
                    double s = std::max(0.0, c.yield_m_err[ip]);
                    if (s > 0.0) { double w = 1.0 / (s*s); wsum_m += w; wval_m += w * c.yield_m[ip]; }
                }
                if (ip < (int)c.acc.size()) {
                    double sa = std::max(0.0, (ip < (int)c.acc_err.size() ? c.acc_err[ip] : 0.0));
                    if (sa > 0.0) { double w = 1.0 / (sa*sa); wsum_a += w; wval_a += w * c.acc[ip]; }
                }
            }

            if (wsum_p > 0.0) { out.yield_p[ip] = wval_p / wsum_p; out.yield_p_err[ip] = std::sqrt(1.0 / wsum_p); }
            if (wsum_m > 0.0) { out.yield_m[ip] = wval_m / wsum_m; out.yield_m_err[ip] = std::sqrt(1.0 / wsum_m); }
            if (wsum_a > 0.0) { out.acc[ip]     = wval_a / wsum_a; out.acc_err[ip]     = std::sqrt(1.0 / wsum_a); }
        }

        combined[key] = std::move(out);
    }
    return !combined.empty();
}

// ---------- main driver ----------
} // anon

void compute_and_plot_unfolding(
    const std::vector<std::string>& periods,           // DVCS_* names for acceptance files
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json_path,         // path to pi0_corrected_counts_all_groups.json
    const std::string& out_root_dir) {
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // Build the runTag list corresponding to DVCS_* periods
    std::vector<std::string> runTagGroups;
    runTagGroups.reserve(periods.size());
    for (const auto& p : periods) runTagGroups.push_back(periodToRunTagKey(p));

    // Load corrected master (value, err per helicity) for the requested runTags.
    GroupHelMap groups;
    if (!load_pi0_corrected_master(total_counts_json_path, runTagGroups, groups)) {
        std::cerr << "[unf][ERROR] Failed to load corrected master json.\n";
        return;
    }

    const fs::path json_dir  = fs::path(out_root_dir) / "jsons";
    std::error_code ec;
    fs::create_directories(json_dir, ec);

    // Store per-period unfolded cells to allow post-facto combinations
    std::map<std::string, std::map<std::tuple<int,int,int>, UnfoldCell>> perPeriodCells;

    auto getGroup = [&](const std::string& runTag)->const std::map<BinKey4,HelVals>*{
        auto it = groups.find(runTag);
        if (it == groups.end()) return nullptr;
        return &it->second;
    };

    // ---------- per-period unfolding ----------
    for (size_t idx = 0; idx < periods.size(); ++idx) {
        const std::string& periodDVCS = periods[idx];          // e.g. DVCS_Fa18_out
        const std::string& runTag     = runTagGroups[idx];     // e.g. fa18_out

        const auto* gmap = getGroup(runTag);
        if (!gmap) {
            std::cerr << "[unf][WARN] No runTag '" << runTag << "' in corrected master — skipping\n";
            continue;
        }

        // acceptance JSON for this DVCS_* period
        const fs::path acc_path = fs::path(out_root_dir) / "jsons" / ("acceptance_" + periodDVCS + ".json");
        AccMap3 accCells;
        if (!load_acceptance_json(acc_path.string(), accCells)) {
            std::cerr << "[unf][WARN] Missing/invalid acceptance for " << periodDVCS << " — skipping.\n";
            continue;
        }

        std::map<std::tuple<int,int,int>, UnfoldCell> outCells;

        for (int ix = 0; ix < (int)xB_bins.size(); ++ix)
        for (int iQ = 0; iQ < (int)Q2_bins.size(); ++iQ)
        for (int it = 0; it < (int)t_bins.size();  ++it) {
            UnfoldCell uc;
            uc.phi_deg = phiCentersDeg();
            uc.yield_p.assign(N_PHI_BINS, 0.0);
            uc.yield_p_err.assign(N_PHI_BINS, 0.0);
            uc.yield_m.assign(N_PHI_BINS, 0.0);
            uc.yield_m_err.assign(N_PHI_BINS, 0.0);
            uc.acc.assign(N_PHI_BINS, 0.0);
            uc.acc_err.assign(N_PHI_BINS, 0.0);

            auto itAcc = accCells.find(std::make_tuple(ix, iQ, it));
            if (itAcc == accCells.end()) {
                outCells[std::make_tuple(ix, iQ, it)] = std::move(uc);
                continue;
            }
            const auto& ac = itAcc->second;

            for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                double A  = (ip < (int)ac.acc.size()     ? ac.acc[ip]     : 0.0);
                double sA = (ip < (int)ac.acc_err.size() ? ac.acc_err[ip] : 0.0);
                A = std::max(0.0, A);
                const double A_clamp = std::max(A, 1e-12);

                uc.acc[ip]     = A;
                uc.acc_err[ip] = sA;

                BinKey4 k4(ix, iQ, it, ip);
                auto itC = gmap->find(k4);

                double Np = 0.0, Nm = 0.0, sNp = 0.0, sNm = 0.0;
                if (itC != gmap->end()) {
                    Np  = std::max(0.0, itC->second.plus);
                    Nm  = std::max(0.0, itC->second.minus);
                    sNp = std::max(0.0, itC->second.eplus);
                    sNm = std::max(0.0, itC->second.eminus);
                }

                if (A > 0.0) {
                    double Up   = Np / A_clamp;
                    double vN   = sNp * sNp;
                    double vA   = sA * sA;
                    double varU = (vN / (A_clamp*A_clamp)) + ((Np*Np) / (A_clamp*A_clamp*A_clamp*A_clamp)) * vA;
                    uc.yield_p[ip]     = Up;
                    uc.yield_p_err[ip] = std::sqrt(std::max(0.0, varU));
                }

                if (A > 0.0) {
                    double Um   = Nm / A_clamp;
                    double vN   = sNm * sNm;
                    double vA   = sA * sA;
                    // FIX: Nm*Nm (typo was "Nm*N m")
                    double varU = (vN / (A_clamp*A_clamp)) + ((Nm*Nm) / (A_clamp*A_clamp*A_clamp*A_clamp)) * vA;
                    uc.yield_m[ip]     = Um;
                    uc.yield_m_err[ip] = std::sqrt(std::max(0.0, varU));
                }
            }

            outCells[std::make_tuple(ix, iQ, it)] = std::move(uc);
        }

        // Persist per-period
        perPeriodCells[periodDVCS] = outCells;

        // JSON and plots
        const fs::path outJ = fs::path(out_root_dir) / "jsons" / ("unfolded_" + periodDVCS + ".json");
        write_unfolded_json(outJ.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, perPeriodCells[periodDVCS]);
        std::cout << "[unf] Wrote unfolded JSON: " << outJ.string() << "\n";

        const fs::path outPlots = fs::path(out_root_dir) / "unfolding" / periodDVCS;
        std::error_code ec2; fs::create_directories(outPlots, ec2);
        plot_cells_for_group(periodDVCS, binning_scheme, xB_bins, Q2_bins, t_bins, perPeriodCells[periodDVCS], outPlots.string());
    }

    // ---------- build combined sets ----------
    auto have = [&](const char* p)->bool {
        return perPeriodCells.find(p) != perPeriodCells.end();
    };

    // Helpers: assemble pointers and write a combined product
    auto combine_and_write = [&](const std::string& label,
                                 const std::vector<std::string>& members) {
        std::vector<const std::map<std::tuple<int,int,int>, UnfoldCell>*> parts;
        parts.reserve(members.size());
        for (const auto& m : members) {
            auto it = perPeriodCells.find(m);
            if (it != perPeriodCells.end()) parts.push_back(&it->second);
        }
        if (parts.empty()) return;

        std::map<std::tuple<int,int,int>, UnfoldCell> combined;
        if (!combine_cells_inverse_variance(parts, combined)) return;

        const fs::path outJ = fs::path(out_root_dir) / "jsons" / ("unfolded_" + label + ".json");
        write_unfolded_json(outJ.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, combined);
        std::cout << "[unf] Wrote unfolded JSON (combined): " << outJ.string() << "\n";

        const fs::path outPlots = fs::path(out_root_dir) / "unfolding" / label;
        std::error_code ec3; fs::create_directories(outPlots, ec3);
        plot_cells_for_group(label, binning_scheme, xB_bins, Q2_bins, t_bins, combined, outPlots.string());
    };

    // Fa18 = inb + out
    if (have("DVCS_Fa18_inb") || have("DVCS_Fa18_out")) {
        combine_and_write("Fa18", {"DVCS_Fa18_inb","DVCS_Fa18_out"});
    }

    // Sp18 = inb + out
    if (have("DVCS_Sp18_inb") || have("DVCS_Sp18_out")) {
        combine_and_write("Sp18", {"DVCS_Sp18_inb","DVCS_Sp18_out"});
    }

    // 10.6 GeV = Sp18 plus Fa18 (ignore Sp19 at 10.2 GeV)
    std::vector<std::string> tenSixMembers;
    if (have("DVCS_Sp18_inb")) tenSixMembers.push_back("DVCS_Sp18_inb");
    if (have("DVCS_Sp18_out")) tenSixMembers.push_back("DVCS_Sp18_out");
    if (have("DVCS_Fa18_inb")) tenSixMembers.push_back("DVCS_Fa18_inb");
    if (have("DVCS_Fa18_out")) tenSixMembers.push_back("DVCS_Fa18_out");
    if (!tenSixMembers.empty()) {
        combine_and_write("10p6GeV", tenSixMembers);
    }

    std::cout << "[unf] Unfolding complete (per-period and combined).\n";
}