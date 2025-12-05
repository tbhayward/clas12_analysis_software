// load_trees.cpp  — canonical-key strict version
#include "load_trees.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>

struct FileDef { std::string tag; std::string path; };

// Helper to load a set of files under one category
static void loadCategory(const std::vector<FileDef>& files,
                         std::map<std::string, TTree*>& container) {
    std::vector<std::string> loadedTags;
    for (const auto& fd : files) {
        TFile* f = TFile::Open(fd.path.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "[Error] Cannot open " << fd.path << "\n";
            continue;
        }
        TTree* tree = dynamic_cast<TTree*>(f->Get("PhysicsEvents"));
        if (!tree) {
            std::cerr << "[Error] PhysicsEvents not found in " << fd.path << "\n";
            continue;
        }
        container[fd.tag] = tree;
        loadedTags.push_back(fd.tag);
    }
    if (!loadedTags.empty()) {
        std::cout << "[Loaded]";
        for (size_t i = 0; i < loadedTags.size(); ++i) {
            std::cout << " " << loadedTags[i];
            if (i + 1 < loadedTags.size()) std::cout << ",";
        }
        std::cout << "\n";
    }
}

void loadTrees(std::map<std::string, TTree*>& dataTrees,
               std::map<std::string, TTree*>& genMcTrees,
               std::map<std::string, TTree*>& recMcTrees,
               std::map<std::string, TTree*>& eppi0DataTrees,
               std::map<std::string, TTree*>& eppi0GenMcTrees,
               std::map<std::string, TTree*>& eppi0RecMcTrees,
               std::map<std::string, TTree*>& eppi0BkgTrees,
               std::map<std::string, TTree*>& radGenMcTrees,
               std::map<std::string, TTree*>& radRecMcTrees)
{
    // ---------------- DVCS data ----------------
    loadCategory({
        {"DVCS_Sp18_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root"},
        {"DVCS_Sp18_out",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root"},
        {"DVCS_Fa18_inb_supp", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_supplemental_epgamma.root"},
        {"DVCS_Fa18_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root"},
        {"DVCS_Fa18_out",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root"},
        {"DVCS_Sp19_inb",      "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"}
    }, dataTrees);

    // ---------------- DVCS generated MC (no radiative) ----------------
    loadCategory({
        {"DVCS_Sp18_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root"},
        {"DVCS_Sp18_out_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_out_10594MeV.root"},
        {"DVCS_Fa18_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_inb_10604MeV.root"},
        {"DVCS_Fa18_out_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_10604MeV.root"},
        {"DVCS_Sp19_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp19_inb_10200MeV.root"}
    }, genMcTrees);

    // ---------------- DVCS reconstructed MC (no radiative) ----------------
    loadCategory({
        {"DVCS_Sp18_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_10594MeV.root"},
        {"DVCS_Sp18_out_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_10594MeV.root"},
        {"DVCS_Fa18_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_10604MeV.root"},
        {"DVCS_Fa18_out_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_10604MeV.root"},
        {"DVCS_Sp19_inb_rec",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_10200MeV.root"}
    }, recMcTrees);

    // ---------------- DVCS generated MC (radiative) ----------------
    loadCategory({
        {"DVCS_Sp18_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp18_inb_10594MeV.root"},
        {"DVCS_Sp18_out_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp18_out_10594MeV.root"},
        {"DVCS_Fa18_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_inb_10604MeV.root"},
        {"DVCS_Fa18_out_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_out_10604MeV.root"},
        {"DVCS_Sp19_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp19_inb_10200MeV.root"}
    }, radGenMcTrees);

    // ---------------- DVCS reconstructed MC (radiative) ----------------
    loadCategory({
        {"DVCS_Sp18_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp18_inb_10594MeV.root"},
        {"DVCS_Sp18_out_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp18_out_10594MeV.root"},
        {"DVCS_Fa18_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_fa18_inb_10604MeV.root"},
        {"DVCS_Fa18_out_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_fa18_out_10604MeV.root"},
        {"DVCS_Sp19_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp19_inb_10200MeV.root"}
    }, radRecMcTrees);

    // ---------------- eπ0 data ----------------
    loadCategory({
        {"DVCS_Sp18_inb_eppi0",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root"},
        {"DVCS_Sp18_out_eppi0",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root"},
        {"DVCS_Fa18_inb_supp_eppi0", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_supplemental_eppi0.root"},
        {"DVCS_Fa18_inb_eppi0",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root"},
        {"DVCS_Fa18_out_eppi0",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root"},
        {"DVCS_Sp19_inb_eppi0",      "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root"}
    }, eppi0DataTrees);

    // ---------------- eπ0 generated MC ----------------
    loadCategory({
        {"DVCS_Sp18_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_inb_10594MeV.root"},
        {"DVCS_Sp18_out_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_out_10594MeV.root"},
        {"DVCS_Fa18_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_inb_10604MeV.root"},
        {"DVCS_Fa18_out_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_out_10604MeV.root"},
        {"DVCS_Sp19_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp19_inb_10200MeV.root"}
    }, eppi0GenMcTrees);

    // ---------------- eπ0 reconstructed MC ----------------
    loadCategory({
        {"DVCS_Sp18_inb_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_10594MeV.root"},
        {"DVCS_Sp18_out_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_10594MeV.root"},
        {"DVCS_Fa18_inb_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_10604MeV.root"},
        {"DVCS_Fa18_out_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_10604MeV.root"},
        {"DVCS_Sp19_inb_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_10200MeV.root"}
    }, eppi0RecMcTrees);

    // ---------------- eπ0 → DVCS background MC ----------------
    loadCategory({
        {"DVCS_Sp18_inb_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_inb_epgamma.root"},
        {"DVCS_Sp18_out_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_out_epgamma.root"},
        {"DVCS_Fa18_inb_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_epgamma.root"},
        {"DVCS_Fa18_out_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_epgamma.root"},
        {"DVCS_Sp19_inb_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_epgamma.root"}
    }, eppi0BkgTrees);
}


// // load_trees.cpp  — canonical-key strict version
// #include "load_trees.h"
// #include <TFile.h>
// #include <TTree.h>
// #include <iostream>
// #include <vector>
// #include <map>
// #include <string>

// struct FileDef { std::string tag; std::string path; };

// // Helper to load a set of files under one category
// static void loadCategory(const std::vector<FileDef>& files,
//                          std::map<std::string, TTree*>& container) {
//     std::vector<std::string> loadedTags;
//     for (const auto& fd : files) {
//         TFile* f = TFile::Open(fd.path.c_str());
//         if (!f || f->IsZombie()) {
//             std::cerr << "[Error] Cannot open " << fd.path << "\n";
//             continue;
//         }
//         TTree* tree = dynamic_cast<TTree*>(f->Get("PhysicsEvents"));
//         if (!tree) {
//             std::cerr << "[Error] PhysicsEvents not found in " << fd.path << "\n";
//             continue;
//         }
//         container[fd.tag] = tree;
//         loadedTags.push_back(fd.tag);
//     }
//     if (!loadedTags.empty()) {
//         std::cout << "[Loaded]";
//         for (size_t i = 0; i < loadedTags.size(); ++i) {
//             std::cout << " " << loadedTags[i];
//             if (i + 1 < loadedTags.size()) std::cout << ",";
//         }
//         std::cout << "\n";
//     }
// }

// void loadTrees(std::map<std::string, TTree*>& dataTrees,
//                std::map<std::string, TTree*>& genMcTrees,
//                std::map<std::string, TTree*>& recMcTrees,
//                std::map<std::string, TTree*>& eppi0DataTrees,
//                std::map<std::string, TTree*>& eppi0GenMcTrees,
//                std::map<std::string, TTree*>& eppi0RecMcTrees,
//                std::map<std::string, TTree*>& eppi0BkgTrees,
//                std::map<std::string, TTree*>& radGenMcTrees,
//                std::map<std::string, TTree*>& radRecMcTrees)
// {
//     // ---------------- DVCS data ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb",      "/volatile/clas12/thayward/dvcs_event_builder_test/dvcs_rga_sp18_inb.root"},
//         {"DVCS_Sp18_out",      "/volatile/clas12/thayward/dvcs_event_builder_test/dvcs_rga_sp18_out.root"},
//         {"DVCS_Fa18_inb_supp", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_supplemental_epgamma.root"},
//         {"DVCS_Fa18_inb",      "/volatile/clas12/thayward/dvcs_event_builder_test/dvcs_rga_fa18_inb.root"},
//         {"DVCS_Fa18_out",      "/volatile/clas12/thayward/dvcs_event_builder_test/dvcs_rga_fa18_out.root"},
//         {"DVCS_Sp19_inb",      "/volatile/clas12/thayward/dvcs_event_builder_test/dvcs_rga_sp19_inb.root"}
//     }, dataTrees);

//     // ---------------- DVCS generated MC (no radiative) ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root"},
//         {"DVCS_Sp18_out_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_out_10594MeV.root"},
//         {"DVCS_Fa18_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_inb_10604MeV.root"},
//         {"DVCS_Fa18_out_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_10604MeV.root"},
//         {"DVCS_Sp19_inb_gen",  "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp19_inb_10200MeV.root"}
//     }, genMcTrees);

//     // ---------------- DVCS reconstructed MC (no radiative) ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb_rec",  "/volatile/clas12/thayward/dvcs_event_builder_test/rec_dvcsgen_sp18_inb.root"},
//         {"DVCS_Sp18_out_rec",  "/volatile/clas12/thayward/dvcs_event_builder_test/rec_dvcsgen_sp18_out.root"},
//         {"DVCS_Fa18_inb_rec",  "/volatile/clas12/thayward/dvcs_event_builder_test/rec_dvcsgen_fa18_inb.root"},
//         {"DVCS_Fa18_out_rec",  "/volatile/clas12/thayward/dvcs_event_builder_test/rec_dvcsgen_fa18_out.root"},
//         {"DVCS_Sp19_inb_rec",  "/volatile/clas12/thayward/dvcs_event_builder_test/rec_dvcsgen_sp19_inb.root"}
//     }, recMcTrees);

//     // ---------------- DVCS generated MC (radiative) ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp18_inb_10594MeV.root"},
//         {"DVCS_Sp18_out_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp18_out_10594MeV.root"},
//         {"DVCS_Fa18_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_inb_10604MeV.root"},
//         {"DVCS_Fa18_out_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_out_10604MeV.root"},
//         {"DVCS_Sp19_inb_gen_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp19_inb_10200MeV.root"}
//     }, radGenMcTrees);

//     // ---------------- DVCS reconstructed MC (radiative) ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp18_inb_10594MeV.root"},
//         {"DVCS_Sp18_out_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp18_out_10594MeV.root"},
//         {"DVCS_Fa18_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_fa18_inb_10604MeV.root"},
//         {"DVCS_Fa18_out_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_fa18_out_10604MeV.root"},
//         {"DVCS_Sp19_inb_rec_rad", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp19_inb_10200MeV.root"}
//     }, radRecMcTrees);

//     // ---------------- eπ0 data ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb_eppi0",      "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_rga_sp18_inb.root"},
//         {"DVCS_Sp18_out_eppi0",      "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_rga_sp18_out.root"},
//         {"DVCS_Fa18_inb_supp_eppi0", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_supplemental_eppi0.root"},
//         {"DVCS_Fa18_inb_eppi0",      "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_rga_fa18_inb.root"},
//         {"DVCS_Fa18_out_eppi0",      "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_rga_fa18_out.root"},
//         {"DVCS_Sp19_inb_eppi0",      "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_rga_sp19_inb.root"}
//     }, eppi0DataTrees);

//     // ---------------- eπ0 generated MC ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_inb_10594MeV.root"},
//         {"DVCS_Sp18_out_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_out_10594MeV.root"},
//         {"DVCS_Fa18_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_inb_10604MeV.root"},
//         {"DVCS_Fa18_out_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_out_10604MeV.root"},
//         {"DVCS_Sp19_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp19_inb_10200MeV.root"}
//     }, eppi0GenMcTrees);

//     // ---------------- eπ0 reconstructed MC ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb_rec_mc", "/volatile/clas12/thayward/dvcs_event_builder_test/rec_aaogen_sp18_inb.root"},
//         {"DVCS_Sp18_out_rec_mc", "/volatile/clas12/thayward/dvcs_event_builder_test/rec_aaogen_sp18_out.root"},
//         {"DVCS_Fa18_inb_rec_mc", "/volatile/clas12/thayward/dvcs_event_builder_test/rec_aaogen_fa18_inb.root"},
//         {"DVCS_Fa18_out_rec_mc", "/volatile/clas12/thayward/dvcs_event_builder_test/rec_aaogen_fa18_out.root"},
//         {"DVCS_Sp19_inb_rec_mc", "/volatile/clas12/thayward/dvcs_event_builder_test/rec_aaogen_sp19_inb.root"}
//     }, eppi0RecMcTrees);

//     // ---------------- eπ0 → DVCS background MC ----------------
//     loadCategory({
//         {"DVCS_Sp18_inb_bkg", "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_bkg_sp18_inb.root"},
//         {"DVCS_Sp18_out_bkg", "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_bkg_sp18_out.root"},
//         {"DVCS_Fa18_inb_bkg", "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_bkg_fa18_inb.root"},
//         {"DVCS_Fa18_out_bkg", "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_bkg_fa18_out.root"},
//         {"DVCS_Sp19_inb_bkg", "/volatile/clas12/thayward/dvcs_event_builder_test/eppi0_bkg_sp19_inb.root"}
//     }, eppi0BkgTrees);
// }