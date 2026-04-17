#include "branch_data_mc_comparison.h"

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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct PeriodDef {
    std::string pretty;
    std::string data_key;
    std::string mc_key;
};

struct RangeInfo {
    double xmin = 0.0;
    double xmax = 1.0;
    int nbins = 100;
    bool is_integer_like = false;
    bool is_bool_like = false;
};

static std::string sanitizeName(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_') out.push_back(c);
        else out.push_back('_');
    }
    return out;
}

static std::vector<PeriodDef> getDvcsPeriods() {
    return {
        {"Sp18 Inb", "DVCS_Sp18_inb", "DVCS_Sp18_inb_rec"},
        {"Sp18 Out", "DVCS_Sp18_out", "DVCS_Sp18_out_rec"},
        {"Fa18 Inb", "DVCS_Fa18_inb", "DVCS_Fa18_inb_rec"},
        {"Fa18 Out", "DVCS_Fa18_out", "DVCS_Fa18_out_rec"},
        {"Sp19 Inb", "DVCS_Sp19_inb", "DVCS_Sp19_inb_rec"}
    };
}

static std::vector<PeriodDef> getEppi0Periods() {
    return {
        {"Sp18 Inb", "DVCS_Sp18_inb_eppi0", "DVCS_Sp18_inb_rec_mc"},
        {"Sp18 Out", "DVCS_Sp18_out_eppi0", "DVCS_Sp18_out_rec_mc"},
        {"Fa18 Inb", "DVCS_Fa18_inb_eppi0", "DVCS_Fa18_inb_rec_mc"},
        {"Fa18 Out", "DVCS_Fa18_out_eppi0", "DVCS_Fa18_out_rec_mc"},
        {"Sp19 Inb", "DVCS_Sp19_inb_eppi0", "DVCS_Sp19_inb_rec_mc"}
    };
}

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
    if (!leaf) {
        leaf = tree->GetLeaf(branch_name.c_str());
    }
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

static void requireAllPeriodTreesPresent(
    const std::string& channel_name,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees)
{
    std::vector<std::string> missing;

    for (const auto& p : periods) {
        if (dataTrees.find(p.data_key) == dataTrees.end()) {
            missing.push_back("data:" + p.data_key);
        }
        if (recMcTrees.find(p.mc_key) == recMcTrees.end()) {
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

    std::vector<std::string> out(common.begin(), common.end());
    std::sort(out.begin(), out.end());
    return out;
}

static void normalizeHist(TH1D* h) {
    if (!h) return;
    const double integral = h->Integral();
    if (integral > 0.0) h->Scale(1.0 / integral);
}

static RangeInfo determineRangeForBranch(
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

static TH1D* fillHistogramFromTree(
    TTree* tree,
    const std::string& hist_name,
    const std::string& branch_name,
    const RangeInfo& rinfo)
{
    if (!tree) {
        throw std::runtime_error("[branch_data_mc_comparison] Null TTree passed to fillHistogramFromTree.");
    }

    TLeaf* leaf = getScalarNumericLeaf(tree, branch_name);
    if (!leaf) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] Could not get scalar numeric leaf for branch "
           << branch_name;
        throw std::runtime_error(ss.str());
    }

    TH1D* h = new TH1D(hist_name.c_str(), "", rinfo.nbins, rinfo.xmin, rinfo.xmax);
    h->SetDirectory(nullptr);

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus(branch_name.c_str(), 1);

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        h->Fill(leaf->GetValue(0));
    }

    tree->SetBranchStatus(branch_name.c_str(), 0);

    return h;
}

static std::string makeCanvasTitle(const std::string& channel_name, const std::string& branch_name) {
    if (channel_name == "eppi0") {
        return "e p -> e p #pi_{0} : " + branch_name;
    }
    return "DVCS : " + branch_name;
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

static void saveOneBranchCanvas(
    const std::string& channel_name,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& branch_name,
    const std::string& channel_out_dir)
{
    const RangeInfo rinfo = determineRangeForBranch(periods, dataTrees, recMcTrees, branch_name);

    std::vector<TH1D*> data_hists;
    std::vector<TH1D*> mc_hists;
    data_hists.reserve(periods.size());
    mc_hists.reserve(periods.size());

    double global_ymax = 0.0;

    for (size_t i = 0; i < periods.size(); ++i) {
        const auto& p = periods[i];

        TH1D* hd = fillHistogramFromTree(
            dataTrees.at(p.data_key),
            "h_data_" + sanitizeName(channel_name) + "_" + sanitizeName(branch_name) + "_" + std::to_string(i),
            branch_name,
            rinfo
        );
        TH1D* hm = fillHistogramFromTree(
            recMcTrees.at(p.mc_key),
            "h_mc_" + sanitizeName(channel_name) + "_" + sanitizeName(branch_name) + "_" + std::to_string(i),
            branch_name,
            rinfo
        );

        normalizeHist(hd);
        normalizeHist(hm);

        styleHistogram(hd, kBlue, 20);
        styleHistogram(hm, kRed, 24);

        global_ymax = std::max(global_ymax, hd->GetMaximum());
        global_ymax = std::max(global_ymax, hm->GetMaximum());

        data_hists.push_back(hd);
        mc_hists.push_back(hm);
    }

    if (global_ymax <= 0.0) global_ymax = 1.0;

    TCanvas c(
        ("c_" + sanitizeName(channel_name) + "_" + sanitizeName(branch_name)).c_str(),
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
    latex.DrawLatex(0.50, 0.70, makeCanvasTitle(channel_name, branch_name).c_str());
    latex.SetTextSize(0.040);
    latex.DrawLatex(0.50, 0.54, "Five RGA periods");
    latex.DrawLatex(0.50, 0.45, "Each histogram normalized to its own integral");
    latex.DrawLatex(0.50, 0.36, ("Range: [" + std::to_string(rinfo.xmin) + ", " + std::to_string(rinfo.xmax) + "]").c_str());

    const std::string out_name =
        channel_out_dir + "/" + sanitizeName(branch_name) + "_data_vs_rec_mc.png";

    c.SaveAs(out_name.c_str());

    for (TH1D* h : data_hists) delete h;
    for (TH1D* h : mc_hists) delete h;
}

static void runChannelComparisons(
    const std::string& channel_name,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& outPlotDir)
{
    requireAllPeriodTreesPresent(channel_name, periods, dataTrees, recMcTrees);

    const std::vector<std::string> branches =
        getCommonBranchesAcrossPeriods(channel_name, periods, dataTrees, recMcTrees);

    if (branches.empty()) {
        std::cout << "[branch_data_mc_comparison] No common scalar numeric branches found for "
                  << channel_name << ". Skipping." << std::endl;
        return;
    }

    const std::string channel_out_dir = outPlotDir + "/" + channel_name;
    std::filesystem::create_directories(channel_out_dir);

    std::cout << "[branch_data_mc_comparison] " << channel_name
              << ": found " << branches.size()
              << " common scalar numeric branches." << std::endl;

    for (size_t i = 0; i < branches.size(); ++i) {
        const std::string& branch_name = branches[i];
        std::cout << "[branch_data_mc_comparison] " << channel_name
                  << " : (" << (i + 1) << "/" << branches.size() << ") "
                  << branch_name << std::endl;

        saveOneBranchCanvas(
            channel_name,
            periods,
            dataTrees,
            recMcTrees,
            branch_name,
            channel_out_dir
        );
    }
}

} // namespace

void runAllBranchDataMcComparisons(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& outPlotDir)
{
    TH1::AddDirectory(kFALSE);
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const std::string base_out_dir = outPlotDir + "/branch_data_mc_comparisons";
    std::filesystem::create_directories(base_out_dir);

    try {
        runChannelComparisons(
            "dvcs",
            getDvcsPeriods(),
            dvcsDataTrees,
            dvcsRecMcTrees,
            base_out_dir
        );
    } catch (const std::exception& e) {
        std::cerr << "[branch_data_mc_comparison] DVCS skipped: " << e.what() << std::endl;
    }

    try {
        runChannelComparisons(
            "eppi0",
            getEppi0Periods(),
            eppi0DataTrees,
            eppi0RecMcTrees,
            base_out_dir
        );
    } catch (const std::exception& e) {
        std::cerr << "[branch_data_mc_comparison] eppi0 skipped: " << e.what() << std::endl;
    }

    std::cout << "[branch_data_mc_comparison] Done." << std::endl;
}