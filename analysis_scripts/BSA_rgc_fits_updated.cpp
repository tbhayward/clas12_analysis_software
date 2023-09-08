#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TH1D.h>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <TSystem.h>
#include <iomanip> // Include this header for std::fixed and std::setprecision

size_t currentFits = 0;
size_t currentBin = 0;
int n = 1;
float total_charge_carbon;
float cpp, cpm, cmp, cmm;

int main(int argc, char *argv[]) {

  // Check for correct number of command line arguments
    if (argc != 5) {
        cout << "Usage: " << argv[0];
        cout << " <data_root_file> <mc_root_file> ";
        cout << " <output_asymmetry_file> <output_kinematic_file>" << endl;
        return 1;
    }

  // Clear the contents of the output_file
  std::ofstream ofs(output_file, std::ios::trunc);
  ofs.close();

  // Clear the contents of the kinematic_file
  std::ofstream ofs2(kinematic_file, std::ios::trunc);
  ofs2.close();

  // load bins from external csv file
  load_bins_from_csv("bins_single_hadron.csv");
  cout<< endl <<"-- Loaded information from bins.csv. " << endl;

  cout<< "Found " << allBins.size() << " sets of bins: " << endl;
  for (size_t i = 0; i < binNames.size(); ++i) {
    cout << binNames[i];
    if (i == binNames.size() - 1) { cout << "."; }
    else { cout << ", "; }
  }
  std::cout << std::endl;

  cout<< "Found " << allBins[currentFits].size() << " bin indices for: " << endl;
  for (size_t i = 0; i < allBins[currentFits].size(); ++i) {
    cout << allBins[currentFits][i];
    if (i == allBins[currentFits].size() - 1) { cout << "."; }
    else { cout << ", "; }
  }
  std::cout << std::endl;

  cout << "Found " << variable_names.size() << " variables: " << endl;
  for (size_t i = 0; i < variable_names.size(); ++i) {
    cout << i << ":" << variable_names[i] << std::flush;
    if (i == variable_names.size() - 1) {
      // cout << ". ";
    } else {
      cout << ", ";
    }
  }
  cout << endl;

  // load run infrom from external csv file
  load_run_info_from_csv("run_info_rgc.csv");
  cout<< endl << endl <<"-- Loaded information from run_info_rgc.csv" << endl;

  cout << "Found " << run_info_list.size() << " runs." << endl;
  cpp = 0; // total accumulated charge of positive beam - positive target
  cpm = 0; // total accumulated charge of positive beam - negative target
  cmp = 0; // total accumulated charge of negative beam - positive target
  cmm = 0; // total accumulated charge of negative beam - negative target
  total_charge_carbon = 0; // total accumulated charge of carbon target
  for (const auto& run_info : run_info_list) {
      if (run_info.target_polarization > 0) {
        cpp += run_info.positive_charge;
        cmp += run_info.negative_charge;
      } else if (run_info.target_polarization < 0) {
        cpm += run_info.positive_charge;
        cmm += run_info.negative_charge;
      } else if (run_info.target_polarization == 0) {
        total_charge_carbon += run_info.total_charge;
      }
  }

  cout << "Total pos-pos (beam-target) charge: " << cpp << " (nc). ";
  cout << "Total pos-neg charge: " << cpm << " (nc). ";
  cout << "Total neg-pos charge: " << cmp << " (nc). ";
  cout << "Total neg-neg charge: " << cmm << " (nc). ";
  cout << "Total unpolarized (carbon) charge: " << total_charge_carbon << " (nc)." << endl;

  // Load data and mc root files
  TFile* gData = new TFile(argv[1], "READ");
  TFile* gMC = new TFile(argv[2], "READ");

  if (!gData->IsOpen() || !gMC->IsOpen()) {
        cout << "Error opening ROOT files (is the location correct?). Exiting." << endl;
        return 2;
    }

}