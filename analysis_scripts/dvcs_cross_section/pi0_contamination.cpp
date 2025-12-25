#include "pi0_contamination.h"

// ROOT includes (must include concrete histogram/axis headers, not just forward decls)
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TAxis.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using json = nlohmann::json;

struct Triple {
    double v;
    double stat;
    double sys;
};

static void fatal(const std::string &msg) {
    std::cerr << "[pi0_contamination] FATAL: " << msg << "\n";
    std::exit(EXIT_FAILURE);
}

static std::string canonical_period_dir(const std::string &period_display) {
    std::string s = period_display;
    for (char &c : s) {
        if (c == ' ') c = '_';
    }
    return s;
}

static std::string period_display_from_tag_strict(const std::string &tag) {
    if (tag.find("Fa18_inb") != std::string::npos || tag.find("fa18_inb") != std::string::npos) return "Fa18 Inb";
    if (tag.find("Fa18_out") != std::string::npos || tag.find("fa18_out") != std::string::npos) return "Fa18 Out";
    if (tag.find("Sp18_inb") != std::string::npos || tag.find("sp18_inb") != std::string::npos) return "Sp18 Inb";
    if (tag.find("Sp18_out") != std::string::npos || tag.find("sp18_out") != std::string::npos) return "Sp18 Out";
    if (tag.find("Sp19_inb") != std::string::npos || tag.find("sp19_inb") != std::string::npos) return "Sp19 Inb";
    fatal("Cannot infer period display name from tag='" + tag + "'");
    return "";
}

static std::string topo_key_from_det(int detector1, int detector2) {
    if (detector1 == 1 && detector2 == 1) return "FD_FD";
    if (detector1 == 2 && detector2 == 1) return "CD_FD";
    if (detector1 == 2 && detector2 == 3) return "CD_FT";
    return "";
}

static std::string contamination_colname_strict(const std::string &period_display) {
    // Must match existing CSV columns exactly.
    return "contamination ratio, " + period_display;
}

static std::string cutblock_key(const std::string &prefix,
                                const std::string &period_display,
                                const std::string &topo_key) {
    return prefix + "_" + canonical_period_dir(period_display) + "_" + topo_key;
}

static double phi_deg_from_phi2(double phi2_rad) {
    double deg = phi2_rad * 180.0 / M_PI;
    while (deg < 0.0) deg += 360.0;
    while (deg >= 360.0) deg -= 360.0;
    return deg;
}

static bool row_accepts_phi(double phi_deg, double phimin, double phimax) {
    double p = phi_deg;
    while (p < 0.0) p += 360.0;
    while (p >= 360.0) p -= 360.0;

    const double lo = phimin;
    const double hi = phimax;

    if (lo <= hi) {
        return (p >= lo && p < hi);
    }
    return (p >= lo || p < hi);
}

static std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        char ch = line[i];
        if (ch == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (ch == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    out.push_back(cur);
    return out;
}

static std::string join_csv_line(const std::vector<std::string> &cols) {
    std::ostringstream os;
    for (size_t i = 0; i < cols.size(); ++i) {
        const std::string &s = cols[i];
        const bool need_quotes = (s.find(',') != std::string::npos) || (s.find('"') != std::string::npos);
        if (need_quotes) {
            os << '"';
            for (char c : s) {
                if (c == '"') os << "\"\"";
                else os << c;
            }
            os << '"';
        } else {
            os << s;
        }
        if (i + 1 < cols.size()) os << ",";
    }
    return os.str();
}

static std::vector<std::vector<std::string>> read_csv(const std::string &path) {
    std::ifstream in(path);
    if (!in) fatal("Cannot open CSV: " + path);

    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(in, line)) {
        rows.push_back(split_csv_line(line));
    }
    if (rows.empty()) fatal("CSV is empty: " + path);
    return rows;
}

static void write_csv_atomic(const std::string &path, const std::vector<std::vector<std::string>> &rows) {
    const std::string tmp = path + ".tmp";
    {
        std::ofstream out(tmp);
        if (!out) fatal("Cannot write temp CSV: " + tmp);
        for (const auto &r : rows) {
            out << join_csv_line(r) << "\n";
        }
    }
    if (gSystem->Rename(tmp.c_str(), path.c_str()) != 0) {
        fatal("Atomic rename failed for CSV: " + tmp + " -> " + path);
    }
}

static int find_column_exact(const std::vector<std::string> &header, const std::string &name) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) return (int)i;
    }
    return -1;
}

static double to_double_strict(const std::string &s, const std::string &field) {
    char *end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    if (!end || *end != '\0') fatal("Cannot parse double for field '" + field + "': '" + s + "'");
    return v;
}

struct RowBin {
    int row_idx; // index in csv rows vector
    double xBmin, xBmax;
    double Q2min, Q2max;
    double tmin, tmax;       // |t| bins from t_abs_min/t_abs_max
    double phimin, phimax;   // degrees
};

static std::vector<RowBin> load_row_bins_from_csv(const std::vector<std::vector<std::string>> &csv) {
    const auto &hdr = csv[0];

    const int c_xBmin  = find_column_exact(hdr, "xBmin");
    const int c_xBmax  = find_column_exact(hdr, "xBmax");
    const int c_Q2min  = find_column_exact(hdr, "Q2min");
    const int c_Q2max  = find_column_exact(hdr, "Q2max");

    // IMPORTANT: Your schema uses t_abs_min/t_abs_max (not tmin/tmax)
    const int c_tmin   = find_column_exact(hdr, "t_abs_min");
    const int c_tmax   = find_column_exact(hdr, "t_abs_max");

    const int c_phimin = find_column_exact(hdr, "phimin");
    const int c_phimax = find_column_exact(hdr, "phimax");

    std::vector<std::string> missing;
    if (c_xBmin  < 0) missing.push_back("xBmin");
    if (c_xBmax  < 0) missing.push_back("xBmax");
    if (c_Q2min  < 0) missing.push_back("Q2min");
    if (c_Q2max  < 0) missing.push_back("Q2max");
    if (c_tmin   < 0) missing.push_back("t_abs_min");
    if (c_tmax   < 0) missing.push_back("t_abs_max");
    if (c_phimin < 0) missing.push_back("phimin");
    if (c_phimax < 0) missing.push_back("phimax");

    if (!missing.empty()) {
        std::ostringstream os;
        os << "CSV missing required bin columns:";
        for (const auto &m : missing) os << " " << m;
        fatal(os.str());
    }

    std::vector<RowBin> bins;
    bins.reserve(csv.size() - 1);

    for (size_t r = 1; r < csv.size(); ++r) {
        const auto &row = csv[r];
        const int need = std::max({c_xBmin,c_xBmax,c_Q2min,c_Q2max,c_tmin,c_tmax,c_phimin,c_phimax});
        if ((int)row.size() <= need) fatal("CSV row too short at r=" + std::to_string(r));

        RowBin b;
        b.row_idx = (int)r;
        b.xBmin  = to_double_strict(row[c_xBmin],  "xBmin");
        b.xBmax  = to_double_strict(row[c_xBmax],  "xBmax");
        b.Q2min  = to_double_strict(row[c_Q2min],  "Q2min");
        b.Q2max  = to_double_strict(row[c_Q2max],  "Q2max");
        b.tmin   = to_double_strict(row[c_tmin],   "t_abs_min");
        b.tmax   = to_double_strict(row[c_tmax],   "t_abs_max");
        b.phimin = to_double_strict(row[c_phimin], "phimin");
        b.phimax = to_double_strict(row[c_phimax], "phimax");
        bins.push_back(b);
    } //endfor

    return bins;
}

static int find_row_index_linear(const std::vector<RowBin> &bins,
                                 double xB, double Q2, double t_abs, double phi_deg) {
    for (const auto &b : bins) {
        if (!(xB >= b.xBmin && xB < b.xBmax)) continue;
        if (!(Q2 >= b.Q2min && Q2 < b.Q2max)) continue;
        if (!(t_abs >= b.tmin && t_abs < b.tmax)) continue;
        if (!row_accepts_phi(phi_deg, b.phimin, b.phimax)) continue;
        return b.row_idx;
    }
    return -1;
}

static std::string format_triple(const Triple &t) {
    std::ostringstream os;
    os << "("
       << std::fixed << std::setprecision(8) << t.v << ", "
       << std::fixed << std::setprecision(8) << t.stat << ", "
       << std::fixed << std::setprecision(8) << t.sys
       << ")";
    return os.str();
}

// ----------------------
// Cut blocks from JSON
// ----------------------
struct CutVar {
    double mean;
    double sigma;
};

using CutMap = std::unordered_map<std::string, std::unordered_map<std::string, CutVar>>;

static CutMap load_cutmap_strict(const std::string &json_path) {
    std::ifstream in(json_path);
    if (!in) fatal("Cannot open combined cuts JSON: " + json_path);

    json j;
    in >> j;

    if (!j.is_object()) fatal("combined cuts JSON is not an object: " + json_path);

    CutMap out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const json &block = it.value();
        if (!block.is_object()) {
            fatal("Cut block '" + key + "' is not an object");
        }

        std::unordered_map<std::string, CutVar> vars;

        for (auto iv = block.begin(); iv != block.end(); ++iv) {
            const std::string varname = iv.key();
            const json &cfg = iv.value();

            if (!cfg.is_object()) fatal("Cut var '" + varname + "' in block '" + key + "' is not an object");
            if (!cfg.contains("mean") || !cfg.contains("sigma")) {
                fatal("Cut var '" + varname + "' in block '" + key + "' missing mean/sigma");
            }
            const double mean = cfg["mean"].get<double>();
            const double sig  = cfg["sigma"].get<double>();
            vars[varname] = CutVar{mean, sig};
        } //endfor

        out[key] = vars;
    } //endfor

    return out;
}

static void bind_required_branch(TTree *t, const char *name, void *addr) {
    if (!t->GetBranch(name)) {
        fatal(std::string("Missing required branch '") + name + "' in PhysicsEvents");
    }
    t->SetBranchStatus(name, 1);
    t->SetBranchAddress(name, addr);
}

struct BaseVars {
    double x;
    double Q2;
    double t1;
    double phi2;
    double open_angle_ep2;
    double pTmiss;
    int detector1;
    int detector2;
};

static BaseVars setup_base_vars(TTree *t) {
    BaseVars v;
    v.x = 0.0; v.Q2 = 0.0; v.t1 = 0.0; v.phi2 = 0.0;
    v.open_angle_ep2 = 0.0; v.pTmiss = 0.0;
    v.detector1 = -9999; v.detector2 = -9999;

    t->SetBranchStatus("*", 0);

    bind_required_branch(t, "x", &v.x);
    bind_required_branch(t, "Q2", &v.Q2);
    bind_required_branch(t, "t1", &v.t1);
    bind_required_branch(t, "phi2", &v.phi2);
    bind_required_branch(t, "open_angle_ep2", &v.open_angle_ep2);
    bind_required_branch(t, "pTmiss", &v.pTmiss);
    bind_required_branch(t, "detector1", &v.detector1);
    bind_required_branch(t, "detector2", &v.detector2);

    return v;
}

static bool passes_global_cuts(double open_angle_ep2_rad, double t1, double pTmiss) {
    const double open_angle_deg = open_angle_ep2_rad * 180.0 / M_PI;
    const double t_abs = std::fabs(t1);
    if (!(open_angle_deg > 5.0)) return false;
    if (!(t_abs < 1.0)) return false;
    if (!(pTmiss <= 0.20)) return false;
    return true;
}

static bool passes_cutblock_3sigma(const std::unordered_map<std::string, CutVar> &vars,
                                   const std::unordered_map<std::string, double> &values) {
    for (const auto &kv : vars) {
        const std::string &name = kv.first;
        const CutVar &cv = kv.second;

        auto it = values.find(name);
        if (it == values.end()) {
            fatal("Internal error: missing value for cut var '" + name + "'");
        }

        const double x = it->second;
        if (!(std::fabs(x - cv.mean) <= 3.0 * cv.sigma)) return false;
    } //endfor
    return true;
}

struct Counts {
    std::unordered_map<int, double> Ndvcs_exp;
    std::unordered_map<int, double> Npi0_exp;
    std::unordered_map<int, double> Npi0_sim;
    std::unordered_map<int, double> Npi0_mis;
};

static void add_count(std::unordered_map<int,double> &m, int row_idx) {
    m[row_idx] += 1.0;
}

static void accumulate_counts_for_tree(TTree *t,
                                       const std::string &tree_tag,
                                       const std::vector<RowBin> &bins,
                                       const CutMap &cutmap,
                                       const std::string &cut_prefix,
                                       const std::string &period_display,
                                       const std::string &which_counter,
                                       Counts &acc) {
    BaseVars v = setup_base_vars(t);

    const std::vector<std::string> topo_keys = {"FD_FD", "CD_FD", "CD_FT"};
    std::set<std::string> required_cut_vars;

    for (const auto &tk : topo_keys) {
        const std::string key = cutblock_key(cut_prefix, period_display, tk);
        auto it = cutmap.find(key);
        if (it == cutmap.end()) {
            fatal("Missing cut-block key in JSON: '" + key + "'");
        }
        for (const auto &vv : it->second) {
            required_cut_vars.insert(vv.first);
        } //endfor
    } //endfor

    std::unordered_map<std::string, double> cut_values;
    cut_values.reserve(required_cut_vars.size());

    std::unordered_map<std::string, double> cut_storage;
    cut_storage.reserve(required_cut_vars.size());

    for (const auto &name : required_cut_vars) {
        cut_storage[name] = 0.0;
        bind_required_branch(t, name.c_str(), &cut_storage[name]);
    } //endfor

    const Long64_t n = t->GetEntries();
    std::cout << "[pi0_contamination] Counting " << which_counter
              << " from tree tag=" << tree_tag
              << " entries=" << n
              << " using cut-prefix=" << cut_prefix
              << " period=" << period_display
              << "\n";

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);

        const std::string topo_key = topo_key_from_det(v.detector1, v.detector2);
        if (topo_key.empty()) continue;

        if (!passes_global_cuts(v.open_angle_ep2, v.t1, v.pTmiss)) continue;

        cut_values.clear();
        const std::string ck = cutblock_key(cut_prefix, period_display, topo_key);
        const auto &vars_for_topo = cutmap.at(ck);

        for (const auto &kv : vars_for_topo) {
            const std::string &varname = kv.first;
            cut_values[varname] = cut_storage.at(varname);
        } //endfor

        if (!passes_cutblock_3sigma(vars_for_topo, cut_values)) continue;

        const double phi_deg = phi_deg_from_phi2(v.phi2);
        const double t_abs = std::fabs(v.t1);

        const int row_idx = find_row_index_linear(bins, v.x, v.Q2, t_abs, phi_deg);
        if (row_idx < 0) continue;

        if (which_counter == "Ndvcs_exp") {
            add_count(acc.Ndvcs_exp, row_idx);
        } else if (which_counter == "Npi0_exp") {
            add_count(acc.Npi0_exp, row_idx);
        } else if (which_counter == "Npi0_sim") {
            add_count(acc.Npi0_sim, row_idx);
        } else if (which_counter == "Npi0_mis") {
            add_count(acc.Npi0_mis, row_idx);
        } else {
            fatal("Unknown counter name: " + which_counter);
        }
    } //endfor
}

static double get_count(const std::unordered_map<int,double> &m, int row_idx) {
    auto it = m.find(row_idx);
    if (it == m.end()) return 0.0;
    return it->second;
}

static Triple compute_contamination(double Nmis, double Nexp, double Nsim, double Ndvcs) {
    Triple out{0.0, 0.0, 0.0};

    if (Nsim <= 0.0) return out;
    if (Ndvcs <= 0.0) return out;

    const double c = (Nmis * Nexp) / (Nsim * Ndvcs);

    const double sNmis = (Nmis > 0.0) ? std::sqrt(Nmis) : 0.0;
    const double sNexp = (Nexp > 0.0) ? std::sqrt(Nexp) : 0.0;
    const double sNsim = (Nsim > 0.0) ? std::sqrt(Nsim) : 0.0;
    const double sNdvc = (Ndvcs > 0.0) ? std::sqrt(Ndvcs) : 0.0;

    double rel2 = 0.0;
    if (Nmis > 0.0) rel2 += (sNmis / Nmis) * (sNmis / Nmis);
    if (Nexp > 0.0) rel2 += (sNexp / Nexp) * (sNexp / Nexp);
    if (Nsim > 0.0) rel2 += (sNsim / Nsim) * (sNsim / Nsim);
    if (Ndvcs > 0.0) rel2 += (sNdvc / Ndvcs) * (sNdvc / Ndvcs);

    out.v = c;
    out.stat = (c > 0.0) ? c * std::sqrt(rel2) : 0.0;
    out.sys = 0.0;
    return out;
}

static void ensure_dir(const std::string &path) {
    if (gSystem->AccessPathName(path.c_str()) == false) return;
    if (gSystem->mkdir(path.c_str(), true) != 0) {
        fatal("Cannot create directory: " + path);
    }
}

static std::string join_path(const std::string &a, const std::string &b) {
    if (a.empty()) return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

struct XBinKey {
    double xBmin, xBmax;
    bool operator<(const XBinKey &o) const {
        if (xBmin != o.xBmin) return xBmin < o.xBmin;
        return xBmax < o.xBmax;
    }
};

struct CellKey {
    double Q2min, Q2max;
    double tmin, tmax;
    bool operator<(const CellKey &o) const {
        if (Q2min != o.Q2min) return Q2min < o.Q2min;
        if (Q2max != o.Q2max) return Q2max < o.Q2max;
        if (tmin  != o.tmin)  return tmin  < o.tmin;
        return tmax < o.tmax;
    }
};

static void plot_contamination_for_period(const std::string &period_display,
                                          const std::vector<RowBin> &bins,
                                          const std::unordered_map<int, Triple> &c_by_row,
                                          const std::string &output_root_dir) {
    std::map<XBinKey, std::map<CellKey, std::vector<std::pair<double, Triple>>>> table;

    for (const auto &b : bins) {
        auto it = c_by_row.find(b.row_idx);
        if (it == c_by_row.end()) continue;

        XBinKey xb{b.xBmin, b.xBmax};
        CellKey ck{b.Q2min, b.Q2max, b.tmin, b.tmax};

        double phi_center = 0.5 * (b.phimin + b.phimax);
        table[xb][ck].push_back({phi_center, it->second});
    } //endfor

    gStyle->SetOptStat(0);

    const std::string period_dir = canonical_period_dir(period_display);
    const std::string outdir = join_path(join_path(output_root_dir, "pi0_contamination_plots"), period_dir);
    ensure_dir(outdir);

    for (const auto &xb_it : table) {
        const XBinKey &xb = xb_it.first;
        const auto &cells = xb_it.second;

        std::vector<std::pair<double,double>> Q2bins;
        std::vector<std::pair<double,double>> tbins;

        for (const auto &cc : cells) {
            const CellKey &ck = cc.first;
            Q2bins.push_back({ck.Q2min, ck.Q2max});
            tbins.push_back({ck.tmin, ck.tmax});
        } //endfor

        std::sort(Q2bins.begin(), Q2bins.end());
        Q2bins.erase(std::unique(Q2bins.begin(), Q2bins.end()), Q2bins.end());
        std::sort(tbins.begin(), tbins.end());
        tbins.erase(std::unique(tbins.begin(), tbins.end()), tbins.end());

        const int ncols = (int)tbins.size();
        const int nrows = (int)Q2bins.size();
        if (ncols <= 0 || nrows <= 0) continue;

        const int W = 300 * ncols + 160;
        const int H = 260 * nrows + 240;

        const std::string cname = "c_pi0_contam_" + period_dir + "_xB_" + std::to_string((int)std::round(xb.xBmin * 1000.0));
        TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), W, H);

        TPad *ptitle = new TPad("ptitle", "ptitle", 0.0, 0.88, 1.0, 1.0);
        ptitle->SetMargin(0.0, 0.0, 0.0, 0.0);
        ptitle->Draw();
        ptitle->cd();

        TLatex t;
        t.SetNDC(true);
        t.SetTextSize(0.45);
        std::ostringstream ttl;
        ttl << period_display << "   xB=[" << std::fixed << std::setprecision(2) << xb.xBmin
            << "," << xb.xBmax << "]";
        t.DrawLatex(0.05, 0.35, ttl.str().c_str());

        c->cd();
        TPad *pgrid = new TPad("pgrid", "pgrid", 0.0, 0.0, 1.0, 0.88);
        pgrid->SetMargin(0.0, 0.0, 0.0, 0.0);
        pgrid->Draw();
        pgrid->cd();
        pgrid->Divide(ncols, nrows);

        double ymax = 0.0;
        for (const auto &cc : cells) {
            for (const auto &pp : cc.second) {
                ymax = std::max(ymax, pp.second.v + pp.second.stat);
            } //endfor
        } //endfor
        if (ymax <= 0.0) ymax = 1.0;
        ymax *= 1.15;

        int ipad = 0;
        for (int ir = 0; ir < nrows; ++ir) {
            for (int ic = 0; ic < ncols; ++ic) {
                ++ipad;
                pgrid->cd(ipad);

                gPad->SetGrid(1, 1);
                gPad->SetLeftMargin(0.160);
                gPad->SetRightMargin(0.07);
                gPad->SetTopMargin(0.22);
                gPad->SetBottomMargin(0.18);

                const auto &Q2b = Q2bins[ir];
                const auto &tb  = tbins[ic];
                CellKey ck{Q2b.first, Q2b.second, tb.first, tb.second};

                auto itcell = cells.find(ck);
                if (itcell == cells.end()) {
                    TH1F *h = new TH1F("hframe","",1,0.0,360.0);
                    h->SetMinimum(0.0);
                    h->SetMaximum(ymax);
                    h->GetXaxis()->SetTitle("#phi (deg)");
                    h->GetYaxis()->SetTitle("Contamination ratio");
                    h->GetXaxis()->CenterTitle();
                    h->GetYaxis()->CenterTitle();
                    h->GetXaxis()->SetNdivisions(505);
                    h->GetXaxis()->SetLabelSize(0.060);
                    h->GetYaxis()->SetLabelSize(0.060);
                    h->GetXaxis()->SetTitleSize(0.070);
                    h->GetYaxis()->SetTitleSize(0.070);
                    h->Draw("AXIS");
                    continue;
                }

                std::vector<std::pair<double, Triple>> pts = itcell->second;
                std::sort(pts.begin(), pts.end(),
                          [](const auto &a, const auto &b) { return a.first < b.first; });

                const int n = (int)pts.size();
                std::vector<double> x(n), y(n), ex(n), ey(n);
                for (int k = 0; k < n; ++k) {
                    x[k]  = pts[k].first;
                    y[k]  = pts[k].second.v;
                    ex[k] = 0.0;
                    ey[k] = pts[k].second.stat;
                } //endfor

                TGraphErrors *g = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
                g->SetMarkerStyle(20);
                g->SetMarkerSize(0.9);
                g->SetLineWidth(1);

                TH1F *frame = (TH1F*)gPad->DrawFrame(0.0, 0.0, 360.0, ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Contamination ratio");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetLabelSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.060);
                frame->GetXaxis()->SetTitleSize(0.070);
                frame->GetYaxis()->SetTitleSize(0.070);

                g->Draw("PE1");

                TLatex a;
                a.SetNDC(true);
                a.SetTextSize(0.060);
                std::ostringstream lab;
                lab << "Q^{2}=[" << std::fixed << std::setprecision(2) << ck.Q2min << "," << ck.Q2max
                    << "]  |t|=[" << ck.tmin << "," << ck.tmax << "]";
                a.DrawLatex(0.12, 0.83, lab.str().c_str());

                TLegend *leg = new TLegend(0.60, 0.73, 0.93, 0.92);
                leg->SetFillStyle(1001);
                leg->SetBorderSize(1);
                leg->SetTextSize(0.055);
                leg->AddEntry(g, "#pi_{0} contamination ratio", "p");
                leg->Draw();
            } //endfor
        } //endfor

        const int idx = (int)std::round(xb.xBmin * 1000.0);
        const std::string outpng = join_path(outdir, "plot_pi0_contamination_" + period_dir + "_xB_" + std::to_string(idx) + ".png");
        c->SaveAs(outpng.c_str());

        delete c;
    } //endfor
}

} // anonymous namespace

bool compute_pi0_contamination_overall(
    const std::map<std::string, TTree*> &dvcsDataTrees,
    const std::map<std::string, TTree*> &eppi0DataTrees,
    const std::map<std::string, TTree*> &eppi0RecMcTrees,
    const std::map<std::string, TTree*> &eppi0BkgTrees,
    const std::string &combined_cuts_json,
    const std::string &csv_main,
    const std::string &output_root_dir,
    int /*max_workers*/) {

    auto csv = read_csv(csv_main);
    const auto bins = load_row_bins_from_csv(csv);

    const CutMap cutmap = load_cutmap_strict(combined_cuts_json);

    std::set<std::string> periods;
    for (const auto &kv : dvcsDataTrees) {
        periods.insert(period_display_from_tag_strict(kv.first));
    } //endfor

    if (periods.empty()) {
        fatal("No DVCS data trees provided to pi0 contamination stage.");
    }

    auto &hdr = csv[0];
    std::unordered_map<std::string, int> col_for_period;

    for (const auto &p : periods) {
        const std::string col = contamination_colname_strict(p);
        const int idx = find_column_exact(hdr, col);
        if (idx < 0) {
            fatal("CSV is missing required existing column: '" + col + "'");
        }
        col_for_period[p] = idx;
    } //endfor

    for (const auto &period_display : periods) {
        TTree *t_dvcs_exp = nullptr;
        TTree *t_pi0_exp  = nullptr;
        TTree *t_pi0_sim  = nullptr;
        TTree *t_pi0_mis  = nullptr;

        for (const auto &kv : dvcsDataTrees) {
            if (period_display_from_tag_strict(kv.first) == period_display) {
                t_dvcs_exp = kv.second;
                break;
            }
        } //endfor
        for (const auto &kv : eppi0DataTrees) {
            if (period_display_from_tag_strict(kv.first) == period_display) {
                t_pi0_exp = kv.second;
                break;
            }
        } //endfor
        for (const auto &kv : eppi0RecMcTrees) {
            if (period_display_from_tag_strict(kv.first) == period_display) {
                t_pi0_sim = kv.second;
                break;
            }
        } //endfor
        for (const auto &kv : eppi0BkgTrees) {
            if (period_display_from_tag_strict(kv.first) == period_display) {
                t_pi0_mis = kv.second;
                break;
            }
        } //endfor

        if (!t_dvcs_exp) fatal("Missing DVCS data tree for period '" + period_display + "'");
        if (!t_pi0_exp)  fatal("Missing eppi0 DATA tree for period '" + period_display + "'");
        if (!t_pi0_sim)  fatal("Missing eppi0 RECO MC tree for period '" + period_display + "'");
        if (!t_pi0_mis)  fatal("Missing eppi0 BKG-as-epg tree for period '" + period_display + "'");

        for (const std::string &tk : std::vector<std::string>{"FD_FD","CD_FD","CD_FT"}) {
            const std::string k_dvcs  = cutblock_key("DVCS",  period_display, tk);
            const std::string k_eppi0 = cutblock_key("eppi0", period_display, tk);
            if (cutmap.find(k_dvcs) == cutmap.end())  fatal("Missing cut block key: " + k_dvcs);
            if (cutmap.find(k_eppi0) == cutmap.end()) fatal("Missing cut block key: " + k_eppi0);
            std::cout << "[pi0_contamination] Using cut keys: " << k_dvcs << " and " << k_eppi0 << "\n";
        } //endfor

        Counts acc;

        accumulate_counts_for_tree(t_dvcs_exp,
                                   "dvcsData:" + period_display,
                                   bins,
                                   cutmap,
                                   "DVCS",
                                   period_display,
                                   "Ndvcs_exp",
                                   acc);

        accumulate_counts_for_tree(t_pi0_mis,
                                   "pi0BkgAsEpg:" + period_display,
                                   bins,
                                   cutmap,
                                   "DVCS",
                                   period_display,
                                   "Npi0_mis",
                                   acc);

        accumulate_counts_for_tree(t_pi0_exp,
                                   "pi0Data:" + period_display,
                                   bins,
                                   cutmap,
                                   "eppi0",
                                   period_display,
                                   "Npi0_exp",
                                   acc);

        accumulate_counts_for_tree(t_pi0_sim,
                                   "pi0RecMc:" + period_display,
                                   bins,
                                   cutmap,
                                   "eppi0",
                                   period_display,
                                   "Npi0_sim",
                                   acc);

        std::unordered_map<int, Triple> c_by_row;
        c_by_row.reserve(bins.size());

        for (const auto &b : bins) {
            const int r = b.row_idx;

            const double Ndvcs = get_count(acc.Ndvcs_exp, r);
            const double Nmis  = get_count(acc.Npi0_mis,  r);
            const double Nexp  = get_count(acc.Npi0_exp,  r);
            const double Nsim  = get_count(acc.Npi0_sim,  r);

            const Triple tr = compute_contamination(Nmis, Nexp, Nsim, Ndvcs);
            c_by_row[r] = tr;

            const int cidx = col_for_period.at(period_display);
            if (cidx >= (int)csv[r].size()) {
                fatal("CSV row " + std::to_string(r) + " too short for contamination column write");
            }
            csv[r][cidx] = format_triple(tr);
        } //endfor

        plot_contamination_for_period(period_display, bins, c_by_row, output_root_dir);
    } //endfor

    write_csv_atomic(csv_main, csv);

    std::cout << "[pi0_contamination] Done. Updated CSV and wrote plots under "
              << join_path(output_root_dir, "pi0_contamination_plots") << "\n";

    return true;
}