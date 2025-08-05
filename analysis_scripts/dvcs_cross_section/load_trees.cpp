#include "load_trees.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>

struct FileDef { std::string tag; std::string path; };

void loadCategory(const std::vector<FileDef>& files,
                  std::map<std::string, TTree*>& container) {
    for (const auto& fd : files) {
        std::cout << "Opening: " << fd.path << std::endl;
        TFile* f = TFile::Open(fd.path.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "  [Error] Cannot open " << fd.path << std::endl;
            continue;
        }
        TTree* tree = dynamic_cast<TTree*>(f->Get("PhysicsEvents"));
        if (!tree) {
            std::cerr << "  [Error] PhysicsEvents not found in " << fd.path << std::endl;
            continue;
        }
        container[fd.tag] = tree;
        std::cout << "  [Loaded] " << fd.tag << std::endl;
    }
}

void loadTrees(std::map<std::string, TTree*>& dataTrees,
               std::map<std::string, TTree*>& genMcTrees,
               std::map<std::string, TTree*>& recMcTrees,
               std::map<std::string, TTree*>& eppi0DataTrees,
               std::map<std::string, TTree*>& eppi0GenMcTrees,
               std::map<std::string, TTree*>& eppi0RecMcTrees) {
    // DVCS data files
    std::vector<FileDef> dvcsData = {
        {"sp18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root"},
        {"sp18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root"},
        {"fa18_inb_supp", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_supplemental_epgamma.root"},
        {"fa18_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root"},
        {"fa18_out", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root"},
        {"sp19_inb", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"}
    };
    loadCategory(dvcsData, dataTrees);

    // DVCS generated MC
    std::vector<FileDef> dvcsGen = {
        {"sp18_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_inb_10594MeV.root"},
        {"sp18_out_gen", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp18_out_10594MeV.root"},
        {"fa18_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_inb_10604MeV.root"},
        {"fa18_out_gen", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_out_10604MeV.root"},
        {"sp19_inb_gen", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp19_inb_10200MeV.root"}
    };
    loadCategory(dvcsGen, genMcTrees);

    // DVCS reconstructed MC
    std::vector<FileDef> dvcsRec = {
        {"sp18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_10594MeV.root"},
        {"sp18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_10594MeV.root"},
        {"fa18_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_10604MeV.root"},
        {"fa18_out_rec", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_10604MeV.root"},
        {"sp19_inb_rec", "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_10200MeV.root"}
    };
    loadCategory(dvcsRec, recMcTrees);

    // eppi0 data files
    std::vector<FileDef> eppi0Data = {
        {"sp18_inb_eppi0", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root"},
        {"sp18_out_eppi0", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root"},
        {"fa18_inb_supp_eppi0", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_supplemental_eppi0.root"},
        {"fa18_inb_eppi0", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root"},
        {"fa18_out_eppi0", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root"},
        {"sp19_inb_eppi0", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root"}
    };
    loadCategory(eppi0Data, eppi0DataTrees);

    // eppi0 generated MC
    std::vector<FileDef> eppi0Gen = {
        {"sp18_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_inb_10594MeV.root"},
        {"sp18_out_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp18_out_10594MeV.root"},
        {"fa18_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_inb_10604MeV.root"},
        {"fa18_out_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_fa18_out_10604MeV.root"},
        {"sp19_inb_gen_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_sp19_inb_10200MeV.root"}
    };
    loadCategory(eppi0Gen, eppi0GenMcTrees);

    // eppi0 reconstructed MC
    std::vector<FileDef> eppi0Rec = {
        {"sp18_inb_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_10594MeV.root"},
        {"sp18_out_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_10594MeV.root"},
        {"fa18_inb_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_10604MeV.root"},
        {"fa18_out_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_10604MeV.root"},
        {"sp19_inb_rec_mc", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_10200MeV.root"}
    };
    loadCategory(eppi0Rec, eppi0RecMcTrees);

    // eppi0 reconstructed MC identified as DVCS background
    std::vector<FileDef> eppi0Bkg = {
        {"sp18_inb_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_inb_epgamma.root"},
        {"sp18_out_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_out_epgamma.root"},
        {"fa18_inb_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_epgamma.root"},
        {"fa18_out_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_epgamma.root"},
        {"sp19_inb_bkg", "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_epgamma.root"}
    };
    loadCategory(eppi0Bkg, eppi0RecMcTrees);
}