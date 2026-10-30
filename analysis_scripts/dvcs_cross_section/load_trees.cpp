// load_trees.cpp

#include "load_trees.h"

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

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

// REBUILD: Load all part files for a period as one logical PhysicsEvents tree.
// TChain objects remain alive for the duration of the analysis.
static std::vector<TChain*>& rebuilt_chain_registry() {
    static std::vector<TChain*> chains;
    return chains;
}

static bool load_rebuilt_parts(const std::string& category,
                               const std::string& mode,
                               const std::map<std::string,std::string>& tags,
                               std::map<std::string,TTree*>& output) {
    namespace fs = std::filesystem;
    const fs::path root("/volatile/clas12/thayward/dvcs_pass2_full_rebuild");
    for (const auto& item : tags) {
        const fs::path dir = root / category / item.second / mode / "parts";
        if (!fs::is_directory(dir)) {
            std::cerr << "[loadTrees] FATAL: missing rebuilt parts directory: " << dir << "\n";
            return false;
        }
        auto* chain = new TChain("PhysicsEvents");
        std::vector<fs::path> paths;
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".root")
                paths.push_back(entry.path());
        }
        std::sort(paths.begin(), paths.end());
        for (const auto& path : paths) chain->Add(path.string().c_str());
        if (paths.empty() || chain->GetNtrees() == 0) {
            std::cerr << "[loadTrees] FATAL: no usable rebuilt files: " << dir << "\n";
            delete chain;
            return false;
        }
        rebuilt_chain_registry().push_back(chain);
        output[item.first] = chain;
        std::cout << "[loadTrees] REBUILD " << item.first << ": "
                  << paths.size() << " parts, " << chain->GetEntries() << " entries\n";
    }
    return true;
}


// REBUILD paired MC: each logical sample chains MCEvents files at ONE current.
// NOTE: consumers must bind gen_*/rec_* and check validity; no aliases are used.
static bool load_paired_mc(const std::string& category,
                           const std::map<std::string,std::string>& periods,
                           const std::map<std::string,std::string>& currents,
                           std::map<std::string,TTree*>& gen,
                           std::map<std::string,TTree*>& rec,
                           const std::string& tag_prefix,
                           bool pi0=false) {
    namespace fs = std::filesystem;
    const fs::path root("/volatile/clas12/thayward/dvcs_pass2_full_rebuild/mc");
    for (const auto& entry : periods) {
        const std::string& tag = entry.first;
        const std::string& period = entry.second;
        const std::string current = currents.at(period);
        const fs::path dir = root / period / category / (pi0 ? "nominal" : current) / "parts";
        if (!fs::is_directory(dir)) {
            std::cerr << "[loadTrees] FATAL: missing paired MC directory " << dir << "\n";
            return false;
        }
        std::vector<fs::path> files;
        for (const auto& file : fs::directory_iterator(dir)) {
            if (!file.is_regular_file() || file.path().extension() != ".root") continue;
            // AAOgen nominal directories may also contain no-background runs.
            if (pi0 && file.path().filename().string().find("_" + current + "_") == std::string::npos) continue;
            files.push_back(file.path());
        }
        std::sort(files.begin(), files.end());
        if (files.empty()) {
            std::cerr << "[loadTrees] FATAL: no paired MC parts at " << dir << " for " << current << "\n";
            return false;
        }
        auto* generated = new TChain("MCEvents");
        auto* reconstructed = new TChain("MCEvents");
        for (const auto& file : files) {
            if (generated->Add(file.string().c_str()) == 0 ||
                reconstructed->Add(file.string().c_str()) == 0) {
                std::cerr << "[loadTrees] FATAL: failed to chain " << file << "\n";
                delete generated;
                delete reconstructed;
                return false;
            }
        }
        if (!generated->GetBranch("gen_kinematics_valid") ||
            !reconstructed->GetBranch("rec_kinematics_valid")) {
            std::cerr << "[loadTrees] FATAL: paired MC validity branches missing in " << dir << "\n";
            delete generated;
            delete reconstructed;
            return false;
        }
        rebuilt_chain_registry().push_back(generated);
        rebuilt_chain_registry().push_back(reconstructed);
        gen[tag_prefix + tag + "_gen"] = generated;
        rec[tag_prefix + tag + "_rec"] = reconstructed;
        std::cout << "[loadTrees] PAIRED MC " << tag << " " << current << ": "
                  << files.size() << " files, " << generated->GetEntries() << " entries\n";
    }
    return true;
}

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
    // ---------------- DVCS data (REBUILD) ----------------
    // ORIGINAL LEGACY BLOCK PRESERVED BELOW (COMMENTED OUT):
    //     // ---------------- DVCS data ----------------
    //     if (!load_category("DVCS data", {
    //         {"DVCS_Sp18_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root"},
    //         {"DVCS_Sp18_out",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root"},
    //         {"DVCS_Fa18_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root"},
    //         {"DVCS_Fa18_out",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root"},
    //         {"DVCS_Sp19_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"}
    //     }, dataTrees)) {
    //         return false;
    //     }
    if (!load_rebuilt_parts("data", "dvcs_data", {
        {"DVCS_Sp18_inb", "sp18_inb"},
        {"DVCS_Sp18_out", "sp18_out"},
        {"DVCS_Fa18_inb", "fa18_inb"},
        {"DVCS_Fa18_out", "fa18_out"},
        {"DVCS_Sp19_inb", "sp19_inb"}
    }, dataTrees)) return false;

    // REBUILD paired DVCS MC (original generated/reconstructed blocks preserved below).
    if (!load_paired_mc("dvcs_mc", {{"Sp18_inb", "sp18_inb"}, {"Sp18_out", "sp18_out"}, {"Fa18_inb", "fa18_inb"}, {"Fa18_out", "fa18_out"}, {"Sp19_inb", "sp19_inb"}}, {{"sp18_inb", "50nA"}, {"sp18_out", "45nA"}, {"fa18_inb", "50nA"}, {"fa18_out", "50nA"}, {"sp19_inb", "50nA"}}, genMcTrees, recMcTrees, "DVCS_")) return false;
    // ============================================================
    // LEGACY ROOT INPUTS — PRESERVED FOR EASY RESTORATION
    // ============================================================
    //     // ---------------- DVCS generated MC for acceptance, no radiative ----------------
    //     if (!load_category("DVCS generated MC for acceptance", {
    //         {"DVCS_Sp18_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_50nA_10594MeV.root"},
    //         {"DVCS_Sp18_out_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_out_45nA_10594MeV.root"},
    //         {"DVCS_Fa18_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_inb_50nA_10604MeV.root"},
    //         {"DVCS_Fa18_out_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_50nA_10604MeV.root"},
    //         {"DVCS_Sp19_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp19_inb_50nA_10200MeV.root"}
    //     }, genMcTrees)) {
    //         return false;
    //     }
    // 
    //     // ---------------- DVCS reconstructed MC for acceptance, no radiative ----------------
    //     if (!load_category("DVCS reconstructed MC for acceptance", {
    //         {"DVCS_Sp18_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_50nA_10594MeV.root"},
    //         {"DVCS_Sp18_out_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_45nA_10594MeV.root"},
    //         {"DVCS_Fa18_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV.root"},
    //         {"DVCS_Fa18_out_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV.root"},
    //         {"DVCS_Sp19_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV.root"}
    //     }, recMcTrees)) {
    //         return false;
    //     }
    // 

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

    // ---------------- eppi0 data (REBUILD) ----------------
    // ORIGINAL LEGACY BLOCK PRESERVED BELOW (COMMENTED OUT):
    //     // ---------------- eppi0 data ----------------
    //     if (!load_category("eppi0 data", {
    //         {"EPPI0_Sp18_inb_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root"},
    //         {"EPPI0_Sp18_out_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root"},
    //         {"EPPI0_Fa18_inb_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root"},
    //         {"EPPI0_Fa18_out_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root"},
    //         {"EPPI0_Sp19_inb_data",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root"}
    //     }, eppi0DataTrees)) {
    //         return false;
    //     }
    if (!load_rebuilt_parts("data", "pi0_data", {
        {"EPPI0_Sp18_inb_data", "sp18_inb"},
        {"EPPI0_Sp18_out_data", "sp18_out"},
        {"EPPI0_Fa18_inb_data", "fa18_inb"},
        {"EPPI0_Fa18_out_data", "fa18_out"},
        {"EPPI0_Sp19_inb_data", "sp19_inb"}
    }, eppi0DataTrees)) return false;

    // REBUILD paired AAOgen MC (original generated/reconstructed blocks preserved below).
    if (!load_paired_mc("pi0_mc", {{"Sp18_inb", "sp18_inb"}, {"Sp18_out", "sp18_out"}, {"Fa18_inb", "fa18_inb"}, {"Fa18_out", "fa18_out"}, {"Sp19_inb", "sp19_inb"}}, {{"sp18_inb", "50nA"}, {"sp18_out", "45nA"}, {"fa18_inb", "50nA"}, {"fa18_out", "50nA"}, {"sp19_inb", "50nA"}}, eppi0GenMcTrees, eppi0RecMcTrees, "EPPI0_", true)) return false;
    // ============================================================
    // LEGACY ROOT INPUTS — PRESERVED FOR EASY RESTORATION
    // ============================================================
    //     // ---------------- eppi0 generated MC ----------------
    //     if (!load_category("eppi0 generated MC", {
    //         {"EPPI0_Sp18_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_inb_50nA_10594MeV.root"},
    //         {"EPPI0_Sp18_out_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_out_45nA_10594MeV.root"},
    //         {"EPPI0_Fa18_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_inb_50nA_10604MeV.root"},
    //         {"EPPI0_Fa18_out_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_out_50nA_10604MeV.root"},
    //         {"EPPI0_Sp19_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp19_inb_50nA_10200MeV.root"}
    //     }, eppi0GenMcTrees)) {
    //         return false;
    //     }
    // 
    //     // ---------------- eppi0 reconstructed MC ----------------
    //     if (!load_category("eppi0 reconstructed MC", {
    //         {"EPPI0_Sp18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_50nA_10594MeV.root"},
    //         {"EPPI0_Sp18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_45nA_10594MeV.root"},
    //         {"EPPI0_Fa18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_50nA_10604MeV.root"},
    //         {"EPPI0_Fa18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_50nA_10604MeV.root"},
    //         {"EPPI0_Sp19_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_50nA_10200MeV.root"}
    //     }, eppi0RecMcTrees)) {
    //         return false;
    //     }
    // 

    // ---------------- eppi0 -> DVCS background MC (REBUILD) ----------------
    // ORIGINAL LEGACY BLOCK PRESERVED BELOW (COMMENTED OUT):
    //     // ---------------- eppi0 -> DVCS background MC ----------------
    //     if (!load_category("eppi0 background MC", {
    //         {"EPPI0_AS_DVCS_Sp18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_inb_50nA_10594MeV_epgamma.root"},
    //         {"EPPI0_AS_DVCS_Sp18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_out_45nA_10594MeV_epgamma.root"},
    //         {"EPPI0_AS_DVCS_Fa18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_50nA_10604MeV_epgamma.root"},
    //         {"EPPI0_AS_DVCS_Fa18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_50nA_10604MeV_epgamma.root"},
    //         {"EPPI0_AS_DVCS_Sp19_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_50nA_10200MeV_epgamma.root"}
    //     }, eppi0BkgTrees)) {
    //         return false;
    //     }
    if (!load_rebuilt_parts("leakage", "pi0_as_dvcs", {
        {"EPPI0_AS_DVCS_Sp18_inb_rec", "sp18_inb"},
        {"EPPI0_AS_DVCS_Sp18_out_rec", "sp18_out"},
        {"EPPI0_AS_DVCS_Fa18_inb_rec", "fa18_inb"},
        {"EPPI0_AS_DVCS_Fa18_out_rec", "fa18_out"},
        {"EPPI0_AS_DVCS_Sp19_inb_rec", "sp19_inb"}
    }, eppi0BkgTrees)) return false;

    // REBUILD: current-study MC comes from the same paired MCEvents format.
    // Keep nominal samples separate from other current configurations.
    const std::map<std::string, std::string> period_tags = {
        {"Sp18_inb", "sp18_inb"}, {"Sp18_out", "sp18_out"},
        {"Fa18_inb", "fa18_inb"}, {"Fa18_out", "fa18_out"},
        {"Sp19_inb", "sp19_inb"}
    };
    const std::map<std::string, std::vector<std::string>> current_configs = {
        {"sp18_inb", {"35nA", "50nA"}},
        {"sp18_out", {"45nA"}},
        {"fa18_inb", {"45nA", "50nA", "55nA"}},
        {"fa18_out", {"40nA", "50nA"}},
        {"sp19_inb", {"50nA"}}
    };
    for (const auto& period : period_tags) {
        for (const auto& current : current_configs.at(period.second)) {
            std::map<std::string, TTree*> gen_tmp, rec_tmp;
            if (!load_paired_mc("dvcs_mc", {{period.first, period.second}},
                                {{period.second, current}}, gen_tmp, rec_tmp, "DVCS_"))
                return false;
            const std::string old_gen = "DVCS_" + period.first + "_gen";
            const std::string old_rec = "DVCS_" + period.first + "_rec";
            const std::string stem = "DVCS_" + period.first + "_" + current + "_";
            currentStudyGenMcTrees[stem + "gen_current"] = gen_tmp.at(old_gen);
            currentStudyRecMcTrees[stem + "rec_current"] = rec_tmp.at(old_rec);
        }
    }

    return true;
}