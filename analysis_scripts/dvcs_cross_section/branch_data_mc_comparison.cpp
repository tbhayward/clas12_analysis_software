#include "branch_data_mc_comparison.h"
#include "branch_plot_config.h"
#include "exclusivity_cuts.h"
#include "global_cuts.h"

#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <TROOT.h>
#include <TH1.h>
#include <TString.h>
#include <TDataType.h>
#include <TObjArray.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct PeriodDef {
    std::string pretty;
    std::string period_label;
    std::string data_key;
    std::string mc_key;
};

struct RangeInfo {
    double xmin = 0.0;
    double xmax = 1.0;
    int nbins = 100;
    bool is_integer_like = false;
    bool is_bool_like = false;
    bool from_config = false;
};

struct ChannelHistStore {
    std::vector<std::string> branch_names;
    std::map<std::string, RangeInfo> range_by_branch;
    std::map<std::string, std::vector<TH1D*> > data_hists_by_branch;
    std::map<std::string, std::vector<TH1D*> > mc_hists_by_branch;
};

// -------------------- helpers --------------------

static std::string sanitizeName(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_') out.push_back(c);
        else out.push_back('_');
    }
    return out;
}

static std::string channelToStr(Channel ch) {
    return (ch == Channel::DVCS) ? "dvcs" : "eppi0";
}

static std::string channelPretty(Channel ch) {
    return (ch == Channel::DVCS) ? "DVCS" : "eppi0";
}

static std::string topoToKey(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }
    return "FD_FD";
}

static std::string periodCode(Channel ch, const std::string& runTag) {
    const std::string prefix = (ch == Channel::DVCS) ? "DVCS_" : "eppi0_";

    if (runTag == "fa18_inb") return prefix + "Fa18_Inb";
    if (runTag == "fa18_out") return prefix + "Fa18_Out";
    if (runTag == "sp18_inb") return prefix + "Sp18_Inb";
    if (runTag == "sp18_out") return prefix + "Sp18_Out";
    if (runTag == "sp19_inb") return prefix + "Sp19_Inb";

    std::ostringstream ss;
    ss << "[branch_data_mc_comparison] FATAL: unknown runTag \"" << runTag
       << "\" in periodCode(). Allowed tags are: "
       << "fa18_inb, fa18_out, sp18_inb, sp18_out, sp19_inb.";
    throw std::runtime_error(ss.str());
}

static bool topologyFromDetectors(int detector1, int detector2, Topology& topo_out) {
    if (detector1 == 1 && detector2 == 1) {
        topo_out = Topology::FD_FD;
        return true;
    }
    if (detector1 == 2 && detector2 == 1) {
        topo_out = Topology::CD_FD;
        return true;
    }
    if (detector1 == 2 && detector2 == 0) {
        topo_out = Topology::CD_FT;
        return true;
    }
    return false;
}

static bool within3Sigma(double val, const Stats& s) {
    return (val >= s.mean - 3.0 * s.std) && (val <= s.mean + 3.0 * s.std);
}

static bool passes3SigmaCuts(const std::map<std::string, Stats>& cuts,
                             const std::map<std::string, double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

static std::vector<PeriodDef> getDvcsPeriods() {
    return {
        {"Sp18 Inb", "sp18_inb", "DVCS_Sp18_inb", "DVCS_Sp18_inb_rec"},
        {"Sp18 Out", "sp18_out", "DVCS_Sp18_out", "DVCS_Sp18_out_rec"},
        {"Fa18 Inb", "fa18_inb", "DVCS_Fa18_inb", "DVCS_Fa18_inb_rec"},
        {"Fa18 Out", "fa18_out", "DVCS_Fa18_out", "DVCS_Fa18_out_rec"},
        {"Sp19 Inb", "sp19_inb", "DVCS_Sp19_inb", "DVCS_Sp19_inb_rec"}
    };
}

static std::vector<PeriodDef> getEppi0Periods() {
    return {
        {"Sp18 Inb", "sp18_inb", "DVCS_Sp18_inb_eppi0", "DVCS_Sp18_inb_rec_mc"},
        {"Sp18 Out", "sp18_out", "DVCS_Sp18_out_eppi0", "DVCS_Sp18_out_rec_mc"},
        {"Fa18 Inb", "fa18_inb", "DVCS_Fa18_inb_eppi0", "DVCS_Fa18_inb_rec_mc"},
        {"Fa18 Out", "fa18_out", "DVCS_Fa18_out_eppi0", "DVCS_Fa18_out_rec_mc"},
        {"Sp19 Inb", "sp19_inb", "DVCS_Sp19_inb_eppi0", "DVCS_Sp19_inb_rec_mc"}
    };
}

// -------------------- tiny JSON parser for combined_cuts.json --------------------

class JsonMiniParser {
public:
    explicit JsonMiniParser(const std::string& text) : s(text), pos(0) {}

    std::map<std::string, CutDict> parseCombinedCuts() {
        std::map<std::string, CutDict> out;
        parseObjectBegin();
        skipWs();
        if (peek() == '}') {
            parseObjectEnd();
            return out;
        }

        while (true) {
            std::string key = parseString();
            parseColon();
            CutDict cd = parseCutDict();
            out[key] = cd;
            skipWs();
            if (peek() == ',') {
                ++pos;
                continue;
            }
            break;
        }

        parseObjectEnd();
        return out;
    }

private:
    const std::string& s;
    size_t pos;

    void skipWs() {
        while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos]))) ++pos;
    }

    char peek() {
        skipWs();
        if (pos >= s.size()) {
            throw std::runtime_error("[branch_data_mc_comparison] JSON parse error: unexpected end of input");
        }
        return s[pos];
    }

    void expect(char c) {
        skipWs();
        if (pos >= s.size() || s[pos] != c) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] JSON parse error: expected '" << c
               << "' at position " << pos;
            throw std::runtime_error(ss.str());
        }
        ++pos;
    }

    void parseObjectBegin() { expect('{'); }
    void parseObjectEnd() { expect('}'); }
    void parseColon() { expect(':'); }

    std::string parseString() {
        skipWs();
        expect('"');
        std::string out;
        while (pos < s.size()) {
            char c = s[pos++];
            if (c == '"') return out;
            if (c == '\\') {
                if (pos >= s.size()) {
                    throw std::runtime_error("[branch_data_mc_comparison] JSON parse error: bad escape");
                }
                out.push_back(s[pos++]);
            } else {
                out.push_back(c);
            }
        }
        throw std::runtime_error("[branch_data_mc_comparison] JSON parse error: unterminated string");
    }

    double parseNumber() {
        skipWs();
        size_t start = pos;
        if (pos < s.size() && (s[pos] == '-' || s[pos] == '+')) ++pos;
        while (pos < s.size() &&
               (std::isdigit(static_cast<unsigned char>(s[pos])) ||
                s[pos] == '.' || s[pos] == 'e' || s[pos] == 'E' ||
                s[pos] == '-' || s[pos] == '+')) {
            ++pos;
        }
        std::string token = s.substr(start, pos - start);
        char* endptr = nullptr;
        double val = std::strtod(token.c_str(), &endptr);
        if (endptr == token.c_str()) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] JSON parse error: invalid number token \"" << token << "\"";
            throw std::runtime_error(ss.str());
        }
        return val;
    }

    Stats parseStats() {
        Stats st{0.0, 0.0};
        parseObjectBegin();

        bool have_mean = false;
        bool have_std  = false;

        skipWs();
        if (peek() == '}') {
            parseObjectEnd();
            return st;
        }

        while (true) {
            std::string key = parseString();
            parseColon();
            double val = parseNumber();

            if (key == "mean") {
                st.mean = val;
                have_mean = true;
            } else if (key == "std") {
                st.std = val;
                have_std = true;
            } else {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] JSON parse error: unexpected Stats key \"" << key << "\"";
                throw std::runtime_error(ss.str());
            }

            skipWs();
            if (peek() == ',') {
                ++pos;
                continue;
            }
            break;
        }

        parseObjectEnd();

        if (!have_mean || !have_std) {
            throw std::runtime_error("[branch_data_mc_comparison] JSON parse error: Stats missing mean/std");
        }

        return st;
    }

    std::map<std::string, Stats> parseStatsMap() {
        std::map<std::string, Stats> out;
        parseObjectBegin();

        skipWs();
        if (peek() == '}') {
            parseObjectEnd();
            return out;
        }

        while (true) {
            std::string key = parseString();
            parseColon();
            Stats st = parseStats();
            out[key] = st;

            skipWs();
            if (peek() == ',') {
                ++pos;
                continue;
            }
            break;
        }

        parseObjectEnd();
        return out;
    }

    CutDict parseCutDict() {
        CutDict out;
        parseObjectBegin();

        bool have_data = false;
        bool have_mc   = false;

        skipWs();
        if (peek() == '}') {
            parseObjectEnd();
            return out;
        }

        while (true) {
            std::string key = parseString();
            parseColon();

            if (key == "data") {
                out.data = parseStatsMap();
                have_data = true;
            } else if (key == "mc") {
                out.mc = parseStatsMap();
                have_mc = true;
            } else {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] JSON parse error: unexpected CutDict key \"" << key << "\"";
                throw std::runtime_error(ss.str());
            }

            skipWs();
            if (peek() == ',') {
                ++pos;
                continue;
            }
            break;
        }

        parseObjectEnd();

        if (!have_data || !have_mc) {
            throw std::runtime_error("[branch_data_mc_comparison] JSON parse error: CutDict missing data/mc");
        }

        return out;
    }
};

static std::map<std::string, CutDict> loadCombinedCutsJson(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: cannot open combined cuts JSON: " << path;
        throw std::runtime_error(ss.str());
    }

    std::ostringstream buf;
    buf << ifs.rdbuf();
    JsonMiniParser parser(buf.str());
    return parser.parseCombinedCuts();
}

// -------------------- branch inspection --------------------

static bool isSupportedScalarNumericLeaf(const TLeaf* leaf) {
    if (!leaf) return false;
    if (leaf->GetLenStatic() != 1) return false;
    if (leaf->GetLeafCount() != nullptr) return false;

    const std::string t = leaf->GetTypeName();

    return (
        t == "Bool_t"   ||
        t == "Char_t"   ||
        t == "UChar_t"  ||
        t == "Short_t"  ||
        t == "UShort_t" ||
        t == "Int_t"    ||
        t == "UInt_t"   ||
        t == "Long_t"   ||
        t == "ULong_t"  ||
        t == "Long64_t" ||
        t == "ULong64_t"||
        t == "Float_t"  ||
        t == "Double_t"
    );
}

static bool isIntegerLikeType(const TLeaf* leaf) {
    if (!leaf) return false;
    const std::string t = leaf->GetTypeName();
    return (
        t == "Bool_t"   ||
        t == "Char_t"   ||
        t == "UChar_t"  ||
        t == "Short_t"  ||
        t == "UShort_t" ||
        t == "Int_t"    ||
        t == "UInt_t"   ||
        t == "Long_t"   ||
        t == "ULong_t"  ||
        t == "Long64_t" ||
        t == "ULong64_t"
    );
}

static bool isBoolType(const TLeaf* leaf) {
    if (!leaf) return false;
    return std::string(leaf->GetTypeName()) == "Bool_t";
}

static TLeaf* getScalarNumericLeaf(TTree* tree, const std::string& branch_name) {
    if (!tree) return nullptr;

    TBranch* br = tree->GetBranch(branch_name.c_str());
    if (!br) return nullptr;

    TLeaf* leaf = br->GetLeaf(branch_name.c_str());
    if (!leaf) leaf = tree->GetLeaf(branch_name.c_str());
    if (!leaf) return nullptr;
    if (!isSupportedScalarNumericLeaf(leaf)) return nullptr;
    return leaf;
}

static std::set<std::string> getSupportedScalarNumericBranchNames(TTree* tree) {
    std::set<std::string> out;
    if (!tree) return out;

    TObjArray* branches = tree->GetListOfBranches();
    if (!branches) return out;

    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* br = dynamic_cast<TBranch*>(branches->At(i));
        if (!br) continue;

        const std::string name = br->GetName();
        TLeaf* leaf = br->GetLeaf(name.c_str());
        if (!leaf) leaf = tree->GetLeaf(name.c_str());
        if (!leaf) continue;

        if (isSupportedScalarNumericLeaf(leaf)) {
            out.insert(name);
        }
    }

    return out;
}

// -------------------- event binder for cuts --------------------

struct BranchBinder {
    int runnum = 0;       bool has_runnum = false;
    int detector1 = 0;    bool has_detector1 = false;
    int detector2 = 0;    bool has_detector2 = false;

    double t1 = 0.0;                bool has_t1 = false;
    double open_angle_ep2 = 0.0;    bool has_open_angle_ep2 = false;
    double Emiss2 = 0.0;            bool has_Emiss2 = false;
    double Mx2 = 0.0;               bool has_Mx2 = false;
    double Mx2_1 = 0.0;             bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;             bool has_Mx2_2 = false;
    double pTmiss = 0.0;            bool has_pTmiss = false;
    double xF = 0.0;                bool has_xF = false;
    double Delta_phi = 0.0;         bool has_Delta_phi = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double theta_pi0_pi0 = 0.0;     bool has_theta_pi0_pi0 = false;

    double e_p = 0.0;       bool has_e_p = false;
    double e_theta = 0.0;   bool has_e_theta = false;
    double e_phi = 0.0;     bool has_e_phi = false;

    double p2_p = 0.0;      bool has_p2_p = false;
    double p2_theta = 0.0;  bool has_p2_theta = false;
    double p2_phi = 0.0;    bool has_p2_phi = false;

    void bindCutBranches(TTree* t, Channel ch) {
        if (!t) return;

        auto ena = [&](const char* n) {
            if (t->GetBranch(n)) t->SetBranchStatus(n, 1);
        };

        ena("runnum");
        ena("detector1");
        ena("detector2");
        ena("t1");
        ena("open_angle_ep2");
        ena("Emiss2");
        ena("Mx2");
        ena("Mx2_1");
        ena("Mx2_2");
        ena("pTmiss");
        ena("xF");
        ena("Delta_phi");
        if (ch == Channel::DVCS) ena("theta_gamma_gamma");
        else                     ena("theta_pi0_pi0");

        ena("e_p");
        ena("e_theta");
        ena("e_phi");
        ena("p2_p");
        ena("p2_theta");
        ena("p2_phi");

        auto bI = [&](const char* n, int* a, bool& f) {
            if (t->GetBranch(n)) {
                t->SetBranchAddress(n, a, nullptr, nullptr, kInt_t, false);
                f = true;
            }
        };

        auto bD = [&](const char* n, double* a, bool& f) {
            if (t->GetBranch(n)) {
                t->SetBranchAddress(n, a, nullptr, nullptr, kDouble_t, false);
                f = true;
            }
        };

        bI("runnum",    &runnum,    has_runnum);
        bI("detector1", &detector1, has_detector1);
        bI("detector2", &detector2, has_detector2);

        bD("t1", &t1, has_t1);
        bD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bD("Emiss2", &Emiss2, has_Emiss2);
        bD("Mx2", &Mx2, has_Mx2);
        bD("Mx2_1", &Mx2_1, has_Mx2_1);
        bD("Mx2_2", &Mx2_2, has_Mx2_2);
        bD("pTmiss", &pTmiss, has_pTmiss);
        bD("xF", &xF, has_xF);
        bD("Delta_phi", &Delta_phi, has_Delta_phi);

        if (ch == Channel::DVCS) bD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        else                     bD("theta_pi0_pi0",     &theta_pi0_pi0,     has_theta_pi0_pi0);

        bD("e_p", &e_p, has_e_p);
        bD("e_theta", &e_theta, has_e_theta);
        bD("e_phi", &e_phi, has_e_phi);
        bD("p2_p", &p2_p, has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi", &p2_phi, has_p2_phi);
    }

    std::map<std::string, double> valuesMap(Channel ch) const {
        std::map<std::string, double> m;
        if (has_Delta_phi) m["Delta_phi"] = Delta_phi;
        if (ch == Channel::DVCS) {
            if (has_theta_gamma_gamma) m["theta_gamma_gamma"] = theta_gamma_gamma;
        } else {
            if (has_theta_pi0_pi0) m["theta_pi0_pi0"] = theta_pi0_pi0;
        }
        if (has_pTmiss) m["pTmiss"] = pTmiss;
        if (has_xF) m["xF"] = xF;
        if (has_Emiss2) m["Emiss2"] = Emiss2;
        if (has_Mx2) m["Mx2"] = Mx2;
        if (has_Mx2_1) m["Mx2_1"] = Mx2_1;
        if (has_Mx2_2) m["Mx2_2"] = Mx2_2;
        return m;
    }
};

// -------------------- tree and branch validation --------------------

static void requireAllPeriodTreesPresent(
    const std::string& channel_name,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees)
{
    std::vector<std::string> missing;

    for (const auto& p : periods) {
        auto itd = dataTrees.find(p.data_key);
        if (itd == dataTrees.end() || itd->second == nullptr) {
            missing.push_back("data:" + p.data_key);
        }

        auto itm = recMcTrees.find(p.mc_key);
        if (itm == recMcTrees.end() || itm->second == nullptr) {
            missing.push_back("rec_mc:" + p.mc_key);
        }
    }

    if (!missing.empty()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] Missing required trees for channel "
           << channel_name << ": ";
        for (size_t i = 0; i < missing.size(); ++i) {
            ss << missing[i];
            if (i + 1 < missing.size()) ss << ", ";
        }
        throw std::runtime_error(ss.str());
    }
}

static std::vector<std::string> getCommonBranchesAcrossPeriods(
    const std::string& channel_name,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees)
{
    requireAllPeriodTreesPresent(channel_name, periods, dataTrees, recMcTrees);

    bool first = true;
    std::set<std::string> common;

    for (const auto& p : periods) {
        TTree* dt = dataTrees.at(p.data_key);
        TTree* mt = recMcTrees.at(p.mc_key);

        std::set<std::string> sdata = getSupportedScalarNumericBranchNames(dt);
        std::set<std::string> smc   = getSupportedScalarNumericBranchNames(mt);

        std::set<std::string> both;
        std::set_intersection(
            sdata.begin(), sdata.end(),
            smc.begin(), smc.end(),
            std::inserter(both, both.begin())
        );

        if (first) {
            common = both;
            first = false;
        } else {
            std::set<std::string> next_common;
            std::set_intersection(
                common.begin(), common.end(),
                both.begin(), both.end(),
                std::inserter(next_common, next_common.begin())
            );
            common.swap(next_common);
        }
    }

    std::vector<std::string> out;
    for (const auto& b : common) {
        const BranchHistConfig* cfg = findBranchPlotConfig(b);
        if (cfg && !cfg->enabled) continue;
        out.push_back(b);
    }

    std::sort(out.begin(), out.end());
    return out;
}

// -------------------- histogram range logic --------------------

static RangeInfo determineAutoRangeForBranch(
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& branch_name)
{
    RangeInfo info;

    bool initialized = false;
    double global_min = 0.0;
    double global_max = 0.0;
    bool integer_like = true;
    bool bool_like = true;

    for (const auto& p : periods) {
        TTree* dt = dataTrees.at(p.data_key);
        TTree* mt = recMcTrees.at(p.mc_key);

        TLeaf* dleaf = getScalarNumericLeaf(dt, branch_name);
        TLeaf* mleaf = getScalarNumericLeaf(mt, branch_name);

        if (!dleaf || !mleaf) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] Branch " << branch_name
               << " is not available as scalar numeric in all required trees.";
            throw std::runtime_error(ss.str());
        }

        integer_like = integer_like && isIntegerLikeType(dleaf) && isIntegerLikeType(mleaf);
        bool_like = bool_like && isBoolType(dleaf) && isBoolType(mleaf);

        const double dmin = dt->GetMinimum(branch_name.c_str());
        const double dmax = dt->GetMaximum(branch_name.c_str());
        const double mmin = mt->GetMinimum(branch_name.c_str());
        const double mmax = mt->GetMaximum(branch_name.c_str());

        if (!initialized) {
            global_min = std::min(dmin, mmin);
            global_max = std::max(dmax, mmax);
            initialized = true;
        } else {
            global_min = std::min(global_min, std::min(dmin, mmin));
            global_max = std::max(global_max, std::max(dmax, mmax));
        }
    }

    info.is_integer_like = integer_like;
    info.is_bool_like = bool_like;
    info.from_config = false;

    if (!initialized) {
        info.xmin = 0.0;
        info.xmax = 1.0;
        info.nbins = 100;
        return info;
    }

    if (bool_like) {
        info.xmin = -0.5;
        info.xmax = 1.5;
        info.nbins = 2;
        return info;
    }

    if (integer_like) {
        const long long imin = static_cast<long long>(std::floor(global_min));
        const long long imax = static_cast<long long>(std::ceil(global_max));
        const long long span = imax - imin + 1LL;

        if (span >= 1LL && span <= 200LL) {
            info.xmin = static_cast<double>(imin) - 0.5;
            info.xmax = static_cast<double>(imax) + 0.5;
            info.nbins = static_cast<int>(span);
            return info;
        }
    }

    if (global_min == global_max) {
        double delta = 1.0;
        if (std::abs(global_min) > 0.0) delta = 0.05 * std::abs(global_min);
        info.xmin = global_min - delta;
        info.xmax = global_max + delta;
        info.nbins = 100;
        return info;
    }

    const double width = global_max - global_min;
    const double pad = 0.05 * width;

    info.xmin = global_min - pad;
    info.xmax = global_max + pad;
    info.nbins = 100;

    return info;
}

static RangeInfo determineRangeForBranch(
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& branch_name)
{
    const BranchHistConfig* cfg = findBranchPlotConfig(branch_name);
    if (cfg) {
        if (!cfg->enabled) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] Branch " << branch_name
               << " is explicitly disabled in branch_plot_config.";
            throw std::runtime_error(ss.str());
        }
        if (!(cfg->xhigh > cfg->xlow) || cfg->nbins <= 0) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] Invalid configured histogram range for branch "
               << branch_name << ": nbins=" << cfg->nbins
               << " xlow=" << cfg->xlow
               << " xhigh=" << cfg->xhigh;
            throw std::runtime_error(ss.str());
        }

        RangeInfo info;
        info.xmin = cfg->xlow;
        info.xmax = cfg->xhigh;
        info.nbins = cfg->nbins;
        info.from_config = true;
        return info;
    }

    return determineAutoRangeForBranch(periods, dataTrees, recMcTrees, branch_name);
}

// -------------------- one-pass filling --------------------

static bool passesGlobalCutsForBinder(const BranchBinder& b,
                                      const std::string& period_label)
{
    if (!(b.has_t1 && b.has_open_angle_ep2 && b.has_pTmiss)) return false;
    if (!(b.has_detector1 && b.has_detector2)) return false;

    GlobalCutConfig gcfg = default_global_cuts();

    // We are intentionally combining all topologies into one histogram.
    // Therefore topology must NOT be forced by the global-cuts config here.
    gcfg.enable_topology_filter = false;

    if (gcfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: dvcsgen ycol cut enabled, but missing one or more "
               << "required branches for period_label='" << period_label
               << "': e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi";
            throw std::runtime_error(ss.str());
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  gcfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                              b.detector1, b.detector2,
                              gcfg);
}

static void enablePlotBranches(TTree* tree, const std::vector<std::string>& branch_names) {
    for (const auto& name : branch_names) {
        if (tree->GetBranch(name.c_str())) {
            tree->SetBranchStatus(name.c_str(), 1);
        }
    }
}

static std::map<std::string, TLeaf*> buildPlotLeafMap(
    TTree* tree,
    const std::vector<std::string>& branch_names)
{
    std::map<std::string, TLeaf*> leaf_map;
    for (const auto& name : branch_names) {
        TLeaf* leaf = getScalarNumericLeaf(tree, name);
        if (!leaf) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: expected scalar numeric leaf for branch "
               << name << " but did not find one while building plot leaf map.";
            throw std::runtime_error(ss.str());
        }
        leaf_map[name] = leaf;
    }
    return leaf_map;
}

static void fillHistogramsForTreeSinglePass(
    TTree* tree,
    Channel ch,
    const std::string& period_label,
    const std::vector<std::string>& branch_names,
    const std::map<std::string, TLeaf*>& leaf_map,
    const std::map<std::string, CutDict>& combinedCuts,
    bool use_data_cuts,
    std::map<std::string, TH1D*>& hist_by_branch)
{
    if (!tree) {
        throw std::runtime_error("[branch_data_mc_comparison] Null TTree passed to fillHistogramsForTreeSinglePass.");
    }

    BranchBinder b;
    b.bindCutBranches(tree, ch);

    GlobalCutConfig gcfg = default_global_cuts();

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        if (b.has_runnum && is_excluded_run(b.runnum, gcfg)) continue;
        if (!(b.has_detector1 && b.has_detector2)) continue;

        Topology topo;
        if (!topologyFromDetectors(b.detector1, b.detector2, topo)) continue;

        if (!passesGlobalCutsForBinder(b, period_label)) continue;

        const std::string cut_key = periodCode(ch, period_label) + "_" + topoToKey(topo);
        auto itCuts = combinedCuts.find(cut_key);
        if (itCuts == combinedCuts.end()) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: missing combined cuts entry for key "
               << cut_key;
            throw std::runtime_error(ss.str());
        }

        const std::map<std::string, double> vals = b.valuesMap(ch);
        const std::map<std::string, Stats>& cut_map = use_data_cuts ? itCuts->second.data : itCuts->second.mc;

        if (!passes3SigmaCuts(cut_map, vals)) continue;

        for (const auto& branch_name : branch_names) {
            auto itLeaf = leaf_map.find(branch_name);
            if (itLeaf == leaf_map.end() || itLeaf->second == nullptr) {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] FATAL: missing plot leaf during fill for branch "
                   << branch_name;
                throw std::runtime_error(ss.str());
            }

            auto itHist = hist_by_branch.find(branch_name);
            if (itHist == hist_by_branch.end() || itHist->second == nullptr) {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] FATAL: missing histogram during fill for branch "
                   << branch_name;
                throw std::runtime_error(ss.str());
            }

            itHist->second->Fill(itLeaf->second->GetValue(0));
        }
    }
}

static void normalizeHist(TH1D* h) {
    if (!h) return;
    const double integral = h->Integral();
    if (integral > 0.0) h->Scale(1.0 / integral);
}

static void styleHistogram(TH1D* h, int color, int marker_style) {
    if (!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker_style);
    h->SetMarkerSize(0.75);
    h->SetLineWidth(2);
    h->SetStats(0);
}

static std::string makeCanvasTitle(const std::string& channel_name, const std::string& branch_name) {
    if (channel_name == "eppi0") return "e p -> e p #pi_{0} : " + branch_name;
    return "DVCS : " + branch_name;
}

// -------------------- histogram creation --------------------

static ChannelHistStore createHistStore(
    Channel ch,
    const std::vector<PeriodDef>& periods,
    const std::vector<std::string>& branch_names,
    const std::map<std::string, RangeInfo>& range_by_branch)
{
    ChannelHistStore store;
    store.branch_names = branch_names;
    store.range_by_branch = range_by_branch;

    for (const auto& branch_name : branch_names) {
        const auto itRange = range_by_branch.find(branch_name);
        if (itRange == range_by_branch.end()) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: missing range info for branch "
               << branch_name;
            throw std::runtime_error(ss.str());
        }

        const RangeInfo& rinfo = itRange->second;

        std::vector<TH1D*> data_vec;
        std::vector<TH1D*> mc_vec;
        data_vec.reserve(periods.size());
        mc_vec.reserve(periods.size());

        for (size_t ip = 0; ip < periods.size(); ++ip) {
            const PeriodDef& p = periods[ip];

            TH1D* hd = new TH1D(
                ("h_data_" + sanitizeName(channelToStr(ch)) + "_" +
                 sanitizeName(p.period_label) + "_" + sanitizeName(branch_name)).c_str(),
                "",
                rinfo.nbins,
                rinfo.xmin,
                rinfo.xmax
            );
            hd->SetDirectory(nullptr);

            TH1D* hm = new TH1D(
                ("h_mc_" + sanitizeName(channelToStr(ch)) + "_" +
                 sanitizeName(p.period_label) + "_" + sanitizeName(branch_name)).c_str(),
                "",
                rinfo.nbins,
                rinfo.xmin,
                rinfo.xmax
            );
            hm->SetDirectory(nullptr);

            data_vec.push_back(hd);
            mc_vec.push_back(hm);
        }

        store.data_hists_by_branch[branch_name] = data_vec;
        store.mc_hists_by_branch[branch_name] = mc_vec;
    }

    return store;
}

static void deleteHistStore(ChannelHistStore& store) {
    for (auto& kv : store.data_hists_by_branch) {
        for (TH1D* h : kv.second) delete h;
    }
    for (auto& kv : store.mc_hists_by_branch) {
        for (TH1D* h : kv.second) delete h;
    }
    store.data_hists_by_branch.clear();
    store.mc_hists_by_branch.clear();
}

// -------------------- plotting --------------------

static void saveOneBranchCanvas(
    Channel ch,
    const std::vector<PeriodDef>& periods,
    const ChannelHistStore& store,
    const std::string& branch_name,
    const std::string& channel_out_dir)
{
    auto itRange = store.range_by_branch.find(branch_name);
    if (itRange == store.range_by_branch.end()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing range for branch "
           << branch_name << " أثناء saveOneBranchCanvas";
        throw std::runtime_error(ss.str());
    }
    const RangeInfo& rinfo = itRange->second;

    auto itData = store.data_hists_by_branch.find(branch_name);
    auto itMc   = store.mc_hists_by_branch.find(branch_name);
    if (itData == store.data_hists_by_branch.end() || itMc == store.mc_hists_by_branch.end()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing hist vectors for branch "
           << branch_name;
        throw std::runtime_error(ss.str());
    }

    std::vector<TH1D*> data_hists = itData->second;
    std::vector<TH1D*> mc_hists   = itMc->second;

    double global_ymax = 0.0;

    for (size_t i = 0; i < periods.size(); ++i) {
        normalizeHist(data_hists[i]);
        normalizeHist(mc_hists[i]);

        styleHistogram(data_hists[i], kBlue, 20);
        styleHistogram(mc_hists[i],   kRed,  24);

        global_ymax = std::max(global_ymax, data_hists[i]->GetMaximum());
        global_ymax = std::max(global_ymax, mc_hists[i]->GetMaximum());
    }

    if (global_ymax <= 0.0) global_ymax = 1.0;

    TCanvas c(
        ("c_" + sanitizeName(channelToStr(ch)) + "_" + sanitizeName(branch_name)).c_str(),
        "",
        2100,
        1200
    );
    c.Divide(3, 2, 0.002, 0.002);

    for (size_t i = 0; i < periods.size(); ++i) {
        c.cd(static_cast<int>(i) + 1);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.06);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.10);
        gPad->SetTickx(1);
        gPad->SetTicky(1);

        TH1D* hd = data_hists[i];
        TH1D* hm = mc_hists[i];

        hd->SetTitle(periods[i].pretty.c_str());
        hd->GetXaxis()->SetTitle(branch_name.c_str());
        hd->GetYaxis()->SetTitle("Normalized counts");
        hd->GetXaxis()->SetTitleOffset(1.10);
        hd->GetYaxis()->SetTitleOffset(1.70);
        hd->SetMaximum(1.25 * global_ymax);

        hd->Draw("HIST");
        hm->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.58, 0.74, 0.92, 0.90);
        leg->SetFillStyle(1001);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->AddEntry(hd, "Data", "l");
        leg->AddEntry(hm, "Reconstructed MC", "l");
        leg->Draw();
    }

    c.cd(6);
    gPad->SetLeftMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.05);
    gPad->SetTopMargin(0.05);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAlign(22);
    latex.SetTextFont(42);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.50, 0.72, makeCanvasTitle(channelToStr(ch), branch_name).c_str());
    latex.SetTextSize(0.040);
    latex.DrawLatex(0.50, 0.56, "Global cuts + topology-matched 3#sigma exclusivity cuts");
    latex.DrawLatex(0.50, 0.47, "All topologies combined into one histogram per period");
    latex.DrawLatex(0.50, 0.38, "Each histogram normalized to its own integral");

    std::ostringstream ss;
    ss << "Range: [" << rinfo.xmin << ", " << rinfo.xmax << "], bins = " << rinfo.nbins;
    latex.DrawLatex(0.50, 0.29, ss.str().c_str());

    if (rinfo.from_config) latex.DrawLatex(0.50, 0.20, "Using explicit branch_plot_config range");
    else                   latex.DrawLatex(0.50, 0.20, "Using automatic branch range");

    const std::string out_name =
        channel_out_dir + "/" + sanitizeName(branch_name) + "_data_vs_rec_mc.png";

    c.SaveAs(out_name.c_str());
}

// -------------------- channel driver --------------------

static void runChannelComparisons(
    Channel ch,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::map<std::string, CutDict>& combinedCuts,
    const std::string& outPlotDir)
{
    const std::string channel_name = channelToStr(ch);

    requireAllPeriodTreesPresent(channel_name, periods, dataTrees, recMcTrees);

    const std::vector<std::string> branch_names =
        getCommonBranchesAcrossPeriods(channel_name, periods, dataTrees, recMcTrees);

    if (branch_names.empty()) {
        std::cout << "[branch_data_mc_comparison] No common scalar numeric branches found for "
                  << channel_name << ". Skipping." << std::endl;
        return;
    }

    std::map<std::string, RangeInfo> range_by_branch;
    for (const auto& branch_name : branch_names) {
        range_by_branch[branch_name] = determineRangeForBranch(periods, dataTrees, recMcTrees, branch_name);
    }

    ChannelHistStore store = createHistStore(ch, periods, branch_names, range_by_branch);

    const std::string channel_out_dir = outPlotDir + "/" + channel_name;
    std::filesystem::create_directories(channel_out_dir);

    std::cout << "[branch_data_mc_comparison] " << channel_name
              << ": found " << branch_names.size()
              << " common scalar numeric branches." << std::endl;

    for (size_t ip = 0; ip < periods.size(); ++ip) {
        const PeriodDef& p = periods[ip];

        {
            TTree* tree = dataTrees.at(p.data_key);
            tree->SetBranchStatus("*", 0);
            enablePlotBranches(tree, branch_names);

            BranchBinder binder;
            binder.bindCutBranches(tree, ch);
            std::map<std::string, TLeaf*> leaf_map = buildPlotLeafMap(tree, branch_names);

            std::map<std::string, TH1D*> hist_by_branch;
            for (const auto& branch_name : branch_names) {
                hist_by_branch[branch_name] = store.data_hists_by_branch[branch_name][ip];
            }

            std::cout << "[branch_data_mc_comparison] Filling data tree for "
                      << channel_name << " period " << p.period_label << std::endl;

            fillHistogramsForTreeSinglePass(
                tree,
                ch,
                p.period_label,
                branch_names,
                leaf_map,
                combinedCuts,
                true,
                hist_by_branch
            );

            tree->SetBranchStatus("*", 1);
        }

        {
            TTree* tree = recMcTrees.at(p.mc_key);
            tree->SetBranchStatus("*", 0);
            enablePlotBranches(tree, branch_names);

            BranchBinder binder;
            binder.bindCutBranches(tree, ch);
            std::map<std::string, TLeaf*> leaf_map = buildPlotLeafMap(tree, branch_names);

            std::map<std::string, TH1D*> hist_by_branch;
            for (const auto& branch_name : branch_names) {
                hist_by_branch[branch_name] = store.mc_hists_by_branch[branch_name][ip];
            }

            std::cout << "[branch_data_mc_comparison] Filling rec MC tree for "
                      << channel_name << " period " << p.period_label << std::endl;

            fillHistogramsForTreeSinglePass(
                tree,
                ch,
                p.period_label,
                branch_names,
                leaf_map,
                combinedCuts,
                false,
                hist_by_branch
            );

            tree->SetBranchStatus("*", 1);
        }
    }

    for (size_t i = 0; i < branch_names.size(); ++i) {
        const std::string& branch_name = branch_names[i];
        std::cout << "[branch_data_mc_comparison] " << channel_name
                  << " : (" << (i + 1) << "/" << branch_names.size() << ") "
                  << branch_name << std::endl;

        try {
            saveOneBranchCanvas(
                ch,
                periods,
                store,
                branch_name,
                channel_out_dir
            );
        } catch (const std::exception& e) {
            std::cerr << "[branch_data_mc_comparison] WARNING: skipping branch "
                      << branch_name << " because: " << e.what() << std::endl;
        }
    }

    deleteHistStore(store);
}

} // namespace

void runAllBranchDataMcComparisons(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& outPlotDir
)
{
    TH1::AddDirectory(kFALSE);
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const std::map<std::string, CutDict> combinedCuts = loadCombinedCutsJson(combined_cuts_json);

    const std::string base_out_dir = outPlotDir + "/branch_data_mc_comparisons";
    std::filesystem::create_directories(base_out_dir);

    try {
        runChannelComparisons(
            Channel::DVCS,
            getDvcsPeriods(),
            dvcsDataTrees,
            dvcsRecMcTrees,
            combinedCuts,
            base_out_dir
        );
    } catch (const std::exception& e) {
        std::cerr << "[branch_data_mc_comparison] DVCS skipped: " << e.what() << std::endl;
    }

    try {
        runChannelComparisons(
            Channel::EPPI0,
            getEppi0Periods(),
            eppi0DataTrees,
            eppi0RecMcTrees,
            combinedCuts,
            base_out_dir
        );
    } catch (const std::exception& e) {
        std::cerr << "[branch_data_mc_comparison] eppi0 skipped: " << e.what() << std::endl;
    }

    std::cout << "[branch_data_mc_comparison] Done." << std::endl;
}