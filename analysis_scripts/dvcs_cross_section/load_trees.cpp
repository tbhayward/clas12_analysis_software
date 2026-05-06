// load_trees.cpp

#include "load_trees.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct FileDef {
    std::string tag;
    std::string path;
};

struct ParsedCurrentMcFile {
    bool valid = false;
    std::string kind;
    std::string period_internal;
    std::string period_tag;
    int current_nA = -1;
    std::string beam_energy_token;
    std::string path;
    std::string basename;
};

static constexpr const char* kTreeName = "PhysicsEvents";

static const std::string kDvcsMcDir =
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen";

static std::vector<TFile*>& open_file_registry() {
    static std::vector<TFile*> files;
    return files;
}

static std::vector<std::string> split_string(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::string cur;

    for (char c : s) {
        if (c == delim) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }

    out.push_back(cur);
    return out;
}

static std::string strip_root_extension(const std::string& basename) {
    static const std::string suffix = ".root";

    if (basename.size() < suffix.size()) {
        return basename;
    }

    if (basename.compare(basename.size() - suffix.size(), suffix.size(), suffix) == 0) {
        return basename.substr(0, basename.size() - suffix.size());
    }

    return basename;
}

static bool ends_with(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) {
        return false;
    }

    return s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::string period_tag_from_internal(const std::string& period_internal) {
    if (period_internal == "rga_sp18_inb") {
        return "Sp18_inb";
    }

    if (period_internal == "rga_sp18_out") {
        return "Sp18_out";
    }

    if (period_internal == "rga_fa18_inb") {
        return "Fa18_inb";
    }

    if (period_internal == "rga_fa18_out") {
        return "Fa18_out";
    }

    if (period_internal == "rga_sp19_inb") {
        return "Sp19_inb";
    }

    throw std::runtime_error("[loadTrees] Unknown period_internal in current-study MC filename: " + period_internal);
}

static ParsedCurrentMcFile parse_current_mc_filename(const std::filesystem::path& path) {
    ParsedCurrentMcFile out;

    out.path = path.string();
    out.basename = path.filename().string();

    if (!ends_with(out.basename, ".root")) {
        return out;
    }

    if (!(out.basename.rfind("gen_", 0) == 0 || out.basename.rfind("rec_", 0) == 0)) {
        return out;
    }

    const std::string stem = strip_root_extension(out.basename);
    const std::vector<std::string> tokens = split_string(stem, '_');

    if (tokens.size() != 7) {
        return out;
    }

    const std::string kind = tokens[0];

    if (kind != "gen" && kind != "rec") {
        return out;
    }

    if (tokens[1] != "dvcsgen") {
        return out;
    }

    const std::string period_internal = tokens[2] + "_" + tokens[3] + "_" + tokens[4];
    const std::string current_token = tokens[5];
    const std::string beam_energy_token = tokens[6];

    if (!ends_with(beam_energy_token, "MeV")) {
        throw std::runtime_error("[loadTrees] Could not parse beam-energy token from file: " + out.basename);
    }

    int current_nA = -1;

    if (current_token == "nobkg") {
        current_nA = 0;
    } else {
        if (!ends_with(current_token, "nA")) {
            throw std::runtime_error("[loadTrees] Could not parse current token from file: " + out.basename);
        }

        const std::string current_number = current_token.substr(0, current_token.size() - 2);

        char* endp = nullptr;
        const long parsed = std::strtol(current_number.c_str(), &endp, 10);

        if (endp == current_number.c_str() || *endp != '\0') {
            throw std::runtime_error("[loadTrees] Invalid current token in file: " + out.basename);
        }

        current_nA = (int)parsed;
    }

    out.valid = true;
    out.kind = kind;
    out.period_internal = period_internal;
    out.period_tag = period_tag_from_internal(period_internal);
    out.current_nA = current_nA;
    out.beam_energy_token = beam_energy_token;

    return out;
}

static bool load_single_tree(const FileDef& fd,
                             std::map<std::string, TTree*>& container,
                             const std::string& category_name) {
    if (fd.tag.empty()) {
        std::cerr << "[loadTrees] FATAL: Empty tag in category " << category_name << "\n";
        return false;
    }

    if (fd.path.empty()) {
        std::cerr << "[loadTrees] FATAL: Empty path for tag " << fd.tag
                  << " in category " << category_name << "\n";
        return false;
    }

    if (container.find(fd.tag) != container.end()) {
        std::cerr << "[loadTrees] FATAL: Duplicate tree tag in category "
                  << category_name << ": " << fd.tag << "\n";
        return false;
    }

    TFile* f = TFile::Open(fd.path.c_str(), "READ");

    if (!f || f->IsZombie()) {
        std::cerr << "[loadTrees] FATAL: Cannot open " << fd.path
                  << " for tag " << fd.tag << "\n";
        return false;
    }

    TTree* tree = dynamic_cast<TTree*>(f->Get(kTreeName));

    if (!tree) {
        std::cerr << "[loadTrees] FATAL: " << kTreeName << " not found in "
                  << fd.path << " for tag " << fd.tag << "\n";
        delete f;
        return false;
    }

    open_file_registry().push_back(f);
    container[fd.tag] = tree;

    return true;
}

static bool load_category(const std::string& category_name,
                          const std::vector<FileDef>& files,
                          std::map<std::string, TTree*>& container) {
    std::vector<std::string> loaded_tags;

    for (const auto& fd : files) {
        if (!load_single_tree(fd, container, category_name)) {
            return false;
        }

        loaded_tags.push_back(fd.tag);
    }

    std::cout << "[loadTrees] Loaded " << category_name << ":";

    for (size_t i = 0; i < loaded_tags.size(); ++i) {
        std::cout << " " << loaded_tags[i];

        if (i + 1 < loaded_tags.size()) {
            std::cout << ",";
        }
    }

    std::cout << "\n";

    return true;
}


static bool load_single_tree_optional(const FileDef& fd,
                                      std::map<std::string, TTree*>& container,
                                      const std::string& category_name) {
    if (fd.tag.empty()) {
        std::cerr << "[loadTrees] WARNING: Empty optional tag in category " << category_name << "\n";
        return true;
    }

    if (fd.path.empty()) {
        std::cerr << "[loadTrees] WARNING: Empty optional path for tag " << fd.tag
                  << " in category " << category_name << "\n";
        return true;
    }

    if (container.find(fd.tag) != container.end()) {
        std::cerr << "[loadTrees] FATAL: Duplicate tree tag in optional category "
                  << category_name << ": " << fd.tag << "\n";
        return false;
    }

    if (!std::filesystem::exists(fd.path)) {
        std::cerr << "[loadTrees] WARNING: Optional " << category_name
                  << " file not found; skipping tag " << fd.tag
                  << ": " << fd.path << "\n";
        return true;
    }

    TFile* f = TFile::Open(fd.path.c_str(), "READ");

    if (!f || f->IsZombie()) {
        std::cerr << "[loadTrees] WARNING: Optional " << category_name
                  << " file could not be opened; skipping tag " << fd.tag
                  << ": " << fd.path << "\n";
        if (f) {
            delete f;
        }
        return true;
    }

    TTree* tree = dynamic_cast<TTree*>(f->Get(kTreeName));

    if (!tree) {
        std::cerr << "[loadTrees] WARNING: Optional " << category_name
                  << " file does not contain " << kTreeName
                  << "; skipping tag " << fd.tag
                  << ": " << fd.path << "\n";
        delete f;
        return true;
    }

    open_file_registry().push_back(f);
    container[fd.tag] = tree;

    return true;
}

static bool load_optional_category(const std::string& category_name,
                                   const std::vector<FileDef>& files,
                                   std::map<std::string, TTree*>& container) {
    std::vector<std::string> loaded_tags;
    int skipped = 0;

    for (const auto& fd : files) {
        const size_t before = container.size();

        if (!load_single_tree_optional(fd, container, category_name)) {
            return false;
        }

        if (container.size() > before) {
            loaded_tags.push_back(fd.tag);
        } else {
            ++skipped;
        }
    }

    std::cout << "[loadTrees] Optional " << category_name << ": loaded "
              << loaded_tags.size() << ", skipped " << skipped;

    if (!loaded_tags.empty()) {
        std::cout << " [";
        for (size_t i = 0; i < loaded_tags.size(); ++i) {
            std::cout << loaded_tags[i];
            if (i + 1 < loaded_tags.size()) {
                std::cout << ", ";
            }
        }
        std::cout << "]";
    }

    std::cout << "\n";

    return true;
}

static bool load_current_study_dvcs_mc(std::map<std::string, TTree*>& currentStudyGenMcTrees,
                                       std::map<std::string, TTree*>& currentStudyRecMcTrees) {
    namespace fs = std::filesystem;

    if (!fs::exists(kDvcsMcDir)) {
        std::cerr << "[loadTrees] FATAL: Current-study DVCS MC directory does not exist: "
                  << kDvcsMcDir << "\n";
        return false;
    }

    if (!fs::is_directory(kDvcsMcDir)) {
        std::cerr << "[loadTrees] FATAL: Current-study DVCS MC path is not a directory: "
                  << kDvcsMcDir << "\n";
        return false;
    }

    std::vector<FileDef> gen_files;
    std::vector<FileDef> rec_files;

    std::set<std::string> seen_gen_tags;
    std::set<std::string> seen_rec_tags;

    for (const auto& entry : fs::directory_iterator(kDvcsMcDir)) {
        if (!entry.is_regular_file()) {
            continue;
        }

        ParsedCurrentMcFile parsed = parse_current_mc_filename(entry.path());

        if (!parsed.valid) {
            continue;
        }

        std::ostringstream tag;
        tag << "DVCS_" << parsed.period_tag << "_"
            << parsed.current_nA << "nA_"
            << parsed.kind << "_current";

        FileDef fd;
        fd.tag = tag.str();
        fd.path = parsed.path;

        if (parsed.kind == "gen") {
            if (!seen_gen_tags.insert(fd.tag).second) {
                std::cerr << "[loadTrees] FATAL: Duplicate current-study generated MC tag: "
                          << fd.tag << "\n";
                return false;
            }

            gen_files.push_back(fd);
        } else if (parsed.kind == "rec") {
            if (!seen_rec_tags.insert(fd.tag).second) {
                std::cerr << "[loadTrees] FATAL: Duplicate current-study reconstructed MC tag: "
                          << fd.tag << "\n";
                return false;
            }

            rec_files.push_back(fd);
        }
    }

    auto sort_by_tag = [](const FileDef& a, const FileDef& b) {
        return a.tag < b.tag;
    };

    std::sort(gen_files.begin(), gen_files.end(), sort_by_tag);
    std::sort(rec_files.begin(), rec_files.end(), sort_by_tag);

    if (gen_files.empty()) {
        std::cerr << "[loadTrees] FATAL: No current-study generated DVCS MC files found in "
                  << kDvcsMcDir << "\n";
        return false;
    }

    if (rec_files.empty()) {
        std::cerr << "[loadTrees] FATAL: No current-study reconstructed DVCS MC files found in "
                  << kDvcsMcDir << "\n";
        return false;
    }

    if (!load_category("DVCS current-study generated MC", gen_files, currentStudyGenMcTrees)) {
        return false;
    }

    if (!load_category("DVCS current-study reconstructed MC", rec_files, currentStudyRecMcTrees)) {
        return false;
    }

    return true;
}

} // namespace

bool loadTrees(std::map<std::string, TTree*>& dataTrees,
               std::map<std::string, TTree*>& genMcTrees,
               std::map<std::string, TTree*>& recMcTrees,
               std::map<std::string, TTree*>& eppi0DataTrees,
               std::map<std::string, TTree*>& eppi0GenMcTrees,
               std::map<std::string, TTree*>& eppi0RecMcTrees,
               std::map<std::string, TTree*>& eppi0BkgTrees,
               std::map<std::string, TTree*>& radGenMcTrees,
               std::map<std::string, TTree*>& radRecMcTrees) {
    std::map<std::string, TTree*> unusedCurrentStudyGenMcTrees;
    std::map<std::string, TTree*> unusedCurrentStudyRecMcTrees;

    return loadTrees(dataTrees,
                     genMcTrees,
                     recMcTrees,
                     eppi0DataTrees,
                     eppi0GenMcTrees,
                     eppi0RecMcTrees,
                     eppi0BkgTrees,
                     radGenMcTrees,
                     radRecMcTrees,
                     unusedCurrentStudyGenMcTrees,
                     unusedCurrentStudyRecMcTrees);
}

bool loadTrees(std::map<std::string, TTree*>& dataTrees,
               std::map<std::string, TTree*>& genMcTrees,
               std::map<std::string, TTree*>& recMcTrees,
               std::map<std::string, TTree*>& eppi0DataTrees,
               std::map<std::string, TTree*>& eppi0GenMcTrees,
               std::map<std::string, TTree*>& eppi0RecMcTrees,
               std::map<std::string, TTree*>& eppi0BkgTrees,
               std::map<std::string, TTree*>& radGenMcTrees,
               std::map<std::string, TTree*>& radRecMcTrees,
               std::map<std::string, TTree*>& currentStudyGenMcTrees,
               std::map<std::string, TTree*>& currentStudyRecMcTrees) {
    // ---------------- DVCS data ----------------
    if (!load_category("DVCS data", {
        {"DVCS_Sp18_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root"},
        {"DVCS_Sp18_out",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root"},
        {"DVCS_Fa18_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root"},
        {"DVCS_Fa18_out",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root"},
        {"DVCS_Sp19_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"}
    }, dataTrees)) {
        return false;
    }

    // ---------------- DVCS generated MC for acceptance, no radiative ----------------
    if (!load_category("DVCS generated MC for acceptance", {
        {"DVCS_Sp18_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_50nA_10594MeV.root"},
        {"DVCS_Sp18_out_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_out_45nA_10594MeV.root"},
        {"DVCS_Fa18_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_inb_50nA_10604MeV.root"},
        {"DVCS_Fa18_out_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_50nA_10604MeV.root"},
        {"DVCS_Sp19_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp19_inb_50nA_10200MeV.root"}
    }, genMcTrees)) {
        return false;
    }

    // ---------------- DVCS reconstructed MC for acceptance, no radiative ----------------
    if (!load_category("DVCS reconstructed MC for acceptance", {
        {"DVCS_Sp18_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_50nA_10594MeV.root"},
        {"DVCS_Sp18_out_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_45nA_10594MeV.root"},
        {"DVCS_Fa18_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV.root"},
        {"DVCS_Fa18_out_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV.root"},
        {"DVCS_Sp19_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV.root"}
    }, recMcTrees)) {
        return false;
    }

    // ---------------- DVCS generated MC for radiative studies ----------------
    if (!load_optional_category("DVCS generated MC radiative", {
        {"DVCS_Sp18_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp18_inb_10594MeV.root"},
        {"DVCS_Sp18_out_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp18_out_10594MeV.root"},
        {"DVCS_Fa18_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_inb_10604MeV.root"},
        {"DVCS_Fa18_out_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_out_10604MeV.root"},
        {"DVCS_Sp19_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp19_inb_10200MeV.root"}
    }, radGenMcTrees)) {
        return false;
    }

    // ---------------- DVCS reconstructed MC for radiative studies ----------------
    if (!load_optional_category("DVCS reconstructed MC radiative", {
        {"DVCS_Sp18_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp18_inb_10594MeV.root"},
        {"DVCS_Sp18_out_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp18_out_10594MeV.root"},
        {"DVCS_Fa18_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_fa18_inb_10604MeV.root"},
        {"DVCS_Fa18_out_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_fa18_out_10604MeV.root"},
        {"DVCS_Sp19_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp19_inb_10200MeV.root"}
    }, radRecMcTrees)) {
        return false;
    }

    // ---------------- eppi0 data ----------------
    if (!load_category("eppi0 data", {
        {"EPPI0_Sp18_inb_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root"},
        {"EPPI0_Sp18_out_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root"},
        {"EPPI0_Fa18_inb_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root"},
        {"EPPI0_Fa18_out_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root"},
        {"EPPI0_Sp19_inb_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root"}
    }, eppi0DataTrees)) {
        return false;
    }

    // ---------------- eppi0 generated MC ----------------
    if (!load_category("eppi0 generated MC", {
        {"EPPI0_Sp18_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_inb_50nA_10594MeV.root"},
        {"EPPI0_Sp18_out_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_out_45nA_10594MeV.root"},
        {"EPPI0_Fa18_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_inb_50nA_10604MeV.root"},
        {"EPPI0_Fa18_out_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_out_50nA_10604MeV.root"},
        {"EPPI0_Sp19_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp19_inb_50nA_10200MeV.root"}
    }, eppi0GenMcTrees)) {
        return false;
    }

    // ---------------- eppi0 reconstructed MC ----------------
    if (!load_category("eppi0 reconstructed MC", {
        {"EPPI0_Sp18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_50nA_10594MeV.root"},
        {"EPPI0_Sp18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_45nA_10594MeV.root"},
        {"EPPI0_Fa18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_50nA_10604MeV.root"},
        {"EPPI0_Fa18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_50nA_10604MeV.root"},
        {"EPPI0_Sp19_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_50nA_10200MeV.root"}
    }, eppi0RecMcTrees)) {
        return false;
    }

    // ---------------- eppi0 -> DVCS background MC ----------------
    if (!load_category("eppi0 background MC", {
        {"EPPI0_AS_DVCS_Sp18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_inb_50nA_10594MeV_epgamma.root"},
        {"EPPI0_AS_DVCS_Sp18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_out_45nA_10594MeV_epgamma.root"},
        {"EPPI0_AS_DVCS_Fa18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_50nA_10604MeV_epgamma.root"},
        {"EPPI0_AS_DVCS_Fa18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_50nA_10604MeV_epgamma.root"},
        {"EPPI0_AS_DVCS_Sp19_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_50nA_10200MeV_epgamma.root"}
    }, eppi0BkgTrees)) {
        return false;
    }

    // ---------------- DVCS additional current-study MC ----------------
    if (!load_current_study_dvcs_mc(currentStudyGenMcTrees, currentStudyRecMcTrees)) {
        return false;
    }

    return true;
}