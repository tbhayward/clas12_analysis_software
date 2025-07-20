// plot_normalized_yields.C
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <cmath>

int main() {
    std::cout << "[Start] C++ per-run normalization\n";
    // make sure output directory exists
    gSystem->mkdir("output", true);

    // --- hard-coded ROOT file paths ---
    std::map<std::string, std::map<std::string, std::string>> filePaths = {
        {"RGC_Su22", {
            {"NH3","/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_NH3_eX.root"},
            {"C",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_C_eX.root"},
            {"CH2","/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_CH2_eX.root"},
            {"He", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_He_eX.root"},
            {"ET", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_ET_eX.root"}
        }},
        {"RGC_Fa22", {
            {"NH3","/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_NH3_eX.root"},
            {"C",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_C_eX.root"},
            {"CH2","/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_CH2_eX.root"},
            {"He", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_He_eX.root"},
            {"ET", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_ET_eX.root"}
        }},
        {"RGC_Sp23", {
            {"NH3","/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_NH3_eX.root"},
            {"C",  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_C_eX.root"},
            {"CH2","/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_CH2_eX.root"},
            {"He", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_He_eX.root"},
            {"ET", "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_ET_eX.root"}
        }}
    };

    // --- read run‐by‐run charge map from CSV ---
    const char* runinfo = "/u/home/thayward/clas12_analysis_software/analysis_scripts/"
                          "asymmetry_extraction/imports/clas12_run_info.csv";
    std::ifstream csv(runinfo);
    if(!csv) {
        std::cerr << "[Error] could not open " << runinfo << "\n";
        return 1;
    }
    std::map<int,double> chargeMap;
    std::string line;
    while(std::getline(csv, line)) {
        if(line.empty() || line[0]=='#') continue;
        std::stringstream ss(line);
        int run; double charge;
        char comma;
        ss >> run >> comma >> charge;
        chargeMap[run] = charge;
    }
    std::cout << "[Loaded] " << chargeMap.size() << " run charges\n";

    // --- prepare output text file ---
    std::ofstream log("output/per_run_integrals.txt");
    log << "Per-run integrals & stats\n"
        << "==========================\n\n";

    // --- histogram settings ---
    const int    nBins = 100;
    const double xmin  = 0.0;
    const double xmax  = 1.0;
    const double wBin  = (xmax - xmin)/nBins;

    // loop periods & targets
    for(auto const& [period, targets] : filePaths) {
      for(auto const& [target, path] : targets) {
        std::cout << "[Process] " << period << " / " << target 
                  << "  ("<< path <<")\n";

        // open file & tree
        TFile f(path.c_str());
        if(!f.IsOpen()) {
          std::cerr<<" [Error] could not open "<<path<<"\n"; 
          continue;
        }
        TTree* tree = (TTree*)f.Get("PhysicsEvents");
        if(!tree) {
          std::cerr<<" [Error] no PhysicsEvents in "<<path<<"\n"; 
          continue;
        }

        // first pass: collect unique runs
        std::set<int> runs;
        int runnum; float x;
        tree->SetBranchAddress("runnum",&runnum);
        tree->SetBranchAddress("x",&x);
        Long64_t nEntries = tree->GetEntries();
        std::cout<<"  [Scan1] "<<nEntries<<" entries for run list\n";
        for(Long64_t i=0;i<nEntries;++i){
          tree->GetEntry(i);
          runs.insert(runnum);
        }
        std::cout<<"  [Found] "<<runs.size()<<" unique runs\n";

        // prepare histograms: map run->bin counts
        std::map<int,std::vector<long long>> hist;
        for(int r : runs){
          hist[r] = std::vector<long long>(nBins,0);
        }

        // second pass: fill counts
        std::cout<<"  [Scan2] filling histograms\n";
        for(Long64_t i=0;i<nEntries;++i){
          tree->GetEntry(i);
          if(x < xmin || x>=xmax) continue;
          int bin = int((x - xmin)/wBin);
          hist[runnum][bin] += 1;
        }

        // compute integrals & collect them
        std::vector<double> integrals;
        log<<"=== Period "<<period<<", Target "<<target<<" ===\n";
        for(auto const& [r, binsVec] : hist){
          auto it = chargeMap.find(r);
          if(it==chargeMap.end()) {
            std::cerr<<"   [Warn] no charge for run "<<r<<"\n";
            continue;
          }
          double totalCounts = 0;
          for(auto c : binsVec) totalCounts += c;
          double integ = totalCounts / it->second;
          integrals.push_back(integ);
          log<<"Run="<<r<<", Integral="<<integ<<"\n";
        }
        // compute mean & stddev
        double mean=0, stddev=0;
        if(!integrals.empty()){
          for(auto v: integrals) mean += v;
          mean /= integrals.size();
          for(auto v: integrals) stddev += (v-mean)*(v-mean);
          stddev = std::sqrt(stddev/integrals.size());
        }
        log<<"Mean="<<mean<<", Std="<<stddev<<"\n\n";

        f.Close();
      }
    }

    log.close();
    std::cout << "[Done] written output/per_run_integrals.txt\n";
    return 0;
}