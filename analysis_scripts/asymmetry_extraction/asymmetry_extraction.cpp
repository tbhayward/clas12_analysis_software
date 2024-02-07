// Standard C++ Library Headers
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
// ROOT Library Headers
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TSystem.h>
#include <TText.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include "TH2.h"
#include "TMinuit.h"
#include "TStyle.h"
#include <TBranch.h>
#include <cstdio>
// tbhayward libraries
#include "common_vars.h"  // Include the common header
#include "load_bins_from_csv.h"
#include "load_run_info_from_csv.h"
#include "dilution_factor.h"
#include "asymmetry_fits.h"
#include "BaseKinematicCuts.h"
#include "KinematicCuts.h"
#include "InclusiveKinematicCuts.h"
#include "SingleHadronKinematicCuts.h"
#include "B2BDihadronKinematicCuts.h"
#include "DihadronKinematicCuts.h"
#include "formatLabelName.h"
#include "readChi2Fits.h"
#include "histConfigs.h"
#include "charge_accumulation.h"
#include "plot_data.h"
#include "modifyTree.h"
#include "fitting_process.h"

// Using namespace declaration
using namespace std;

TTreeReader dataReader;  // Declare as global variable
TTreeReader mcReader;  // Declare as global variable

BaseKinematicCuts* kinematicCuts = nullptr;
BaseKinematicCuts* mckinematicCuts = nullptr;

int currentFits = 0;
int currentBin = 0;
int n = 1;
double total_charge_carbon = 0;
double cmm = 0; 
double cpm = 0; 
double cmp = 0; 
double cpp = 0; 
int channel = 1;
std::string mlmPrefix = "xF";

int main(int argc, char *argv[]) {
  // Start the timer
  auto start_time = std::chrono::high_resolution_clock::now();

  // initialize ROOT application for graphics
  TApplication theApp("App", nullptr, nullptr);
  // Set ROOT to batch mode
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // Check for correct number of command line arguments
  if (argc < 3) {
      std::cout << "Usage: " << argv[0];
      std::cout << " <data_root_file> <mc_root_file> [channel]" << std::endl;
      return 1;
  }

  cout << endl << endl;
  // Set default channel to 1
  channel = 1;  // Default value
  if (argc >= 4) {
      try {
          channel = std::stoi(argv[3]);
          if (channel < 0 || channel > 3) {
              cout << "Invalid channel specified. Defaulting to single hadron." << endl;
              channel = 1;
          }
      } catch (const std::exception& e) {
          cout << "Error parsing channel. Defaulting to channel single hadron." << endl;
      }
  } else {
      cout << "No channel specified. Defaulting to channel single hadron." << endl;
  }
  cout << "Using channel: " << channel;

  // Allocate kinematicCuts and mckinematicCuts based on the channel
  switch (channel) {
      case 0:
          kinematicCuts = new InclusiveKinematicCuts(dataReader);
          mckinematicCuts = new InclusiveKinematicCuts(mcReader);
          break;
      case 1:
          kinematicCuts = new SingleHadronKinematicCuts(dataReader);
          mckinematicCuts = new SingleHadronKinematicCuts(mcReader);
          break;
      case 2:
          kinematicCuts = new B2BDihadronKinematicCuts(dataReader);
          mckinematicCuts = new B2BDihadronKinematicCuts(mcReader);
          break;
      case 3:
          kinematicCuts = new DihadronKinematicCuts(dataReader);
          mckinematicCuts = new DihadronKinematicCuts(mcReader);
          break;
  }

  cout << endl << endl;
  std::string inputFileName(argv[2]);
  std::size_t found = inputFileName.find_last_of("/\\");
  std::string directoryPath = inputFileName.substr(0, found);
  std::string outputFileName = directoryPath + "/temp_mc.root";

  modifyTree(argv[2], outputFileName.c_str());

  // Load data and mc root files
  TFile* data_file = new TFile(argv[1], "READ");
  TFile* mc_file = new TFile(outputFileName.c_str(), "READ");
  if (!data_file->IsOpen() || !mc_file->IsOpen()) {
    cout << "Error opening ROOT files (is the location correct?). Exiting." << endl;
    return 2;
  } else {
    cout << "-- ROOT files opened successfully." << endl;
  }
  
  TTree* data = (TTree*)data_file->Get("PhysicsEvents");
  TTree* mc = (TTree*)mc_file->Get("PhysicsEvents");

  if (!data || !mc) {
    cout << "-- Error getting trees from ROOT files." << endl;
    return 3;
  } else {
    cout << "-- Trees successfully extracted from ROOT files." << endl << endl;
  }

  dataReader.SetTree(data);  // Initialize the global variable
  mcReader.SetTree(mc);  // Initialize the global variable

  // Generate output file names based on the input data file name and current time
  std::string dataRootFileName = argv[1];
  std::string baseName = dataRootFileName.substr(dataRootFileName.find_last_of("/\\") + 1);
  baseName = baseName.substr(0, baseName.find_last_of(".")); // Remove the file extension

  // Check if the file extension is .root
  if (dataRootFileName.substr(dataRootFileName.find_last_of(".") + 1) != "root") {
      std::cerr << "The input file must be a .root file." << std::endl;
      return 1;
  }

  // Get the current time in UTC
  auto time_now = std::chrono::system_clock::now();
  std::time_t time_now_t = std::chrono::system_clock::to_time_t(time_now);
  // Convert it to UTC (GMT)
  std::tm gmt = *std::gmtime(&time_now_t);
  // Apply the Eastern Time Zone offset
  // EST: UTC-5, EDT: UTC-4. Assuming EST here
  gmt.tm_hour -= 5;
  // Correct for out-of-range hours
  std::mktime(&gmt);
  char buffer[80];
  std::strftime(buffer, sizeof(buffer), "%m_%d_%H%M%S", &gmt);

  std::string timeStamp(buffer);
  std::string output_file = "output/results/asymmetries_" + baseName + 
    "_timeStamp_" + timeStamp + ".txt";
  std::string kinematic_file = "output/results/kinematics_" + baseName + 
    "_timeStamp_" + timeStamp + ".txt";

  // Clear the contents of the output files
  std::ofstream ofs(output_file, std::ios::trunc);
  ofs.close();
  std::ofstream ofs2(kinematic_file, std::ios::trunc);
  ofs2.close();

  // load bins from external csv file
  load_bins_from_csv("imports/bins_single_hadron.csv");
  cout<< endl <<"-- Loaded information from bins.csv. " << endl;

  cout<< "Found " << allBins.size() << " sets of bins: " << endl;
  for (size_t i = 0; i < binNames.size(); ++i) {
    cout << binNames[i];
    if (i == binNames.size() - 1) { cout << "."; }
    else { cout << ", "; }
  }
  cout << endl;

  cout<< "Found " << allBins[currentFits].size() << " bin indices for: " << endl;
  for (size_t i = 0; i < allBins[currentFits].size(); ++i) {
    cout << allBins[currentFits][i];
    if (i == allBins[currentFits].size() - 1) { cout << "."; }
    else { cout << ", "; }
  }
  cout << endl;

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

  // load run info from external csv file
  load_run_info_from_csv("imports/clas12_run_info.csv");
  cout<< endl << endl <<"-- Loaded information from run_info_rgc.csv" << endl;
  cout << "Found " << run_info_list.size() << " runs." << endl;
  // charge_acuumulation determines the total charges from the runs in supplied dataFile
  // by comparing to master list of CLAS12 runs 
  charge_accumulation(dataReader, run_info_list);
  cout << "Total pos-pos (beam-target) charge: " << cpp << " (nC). ";
  cout << "Total pos-neg charge: " << cpm << " (nC). ";
  cout << "Total neg-pos charge: " << cmp << " (nC). ";
  cout << "Total neg-neg charge: " << cmm << " (nC). ";
  cout << "Total unpolarized (carbon) charge: " << total_charge_carbon << " (nC)."<< endl << endl;

  createIntegratedKinematicPlots();
  // createIntegratedKinematicPlotsForBinsAndFits();
  // createCorrelationPlotsforrunnum();
  // createCorrelationPlots();
  currentFits=0;
  dataReader.Restart(); mcReader.Restart();

  for (size_t i = 0; i < allBins.size(); ++i) {
    cout << "-- Beginning kinematic fits." << endl;
    for (int asymmetry = 0; asymmetry < 3; ++asymmetry){
      if (asymmetry > 0 && cpp == 1) {
        cout << "Skipping TSA and DSA for unpolarized target data." << endl;
        continue;
      }
      switch (asymmetry) {
        case 0: cout << "    Beginning chi2 BSA." << endl; break;
        case 1: cout << "    Beginning chi2 TSA." << endl; break;
        case 2: cout << "    Beginning chi2 DSA." << endl; break;
      }
      switch (channel) {
        case 0: calculate_inclusive(output_file.c_str(), kinematic_file.c_str(), 
        binNames[i], asymmetry); break;
        case 1: performChi2Fits_single_hadron(output_file.c_str(), kinematic_file.c_str(), 
        binNames[i], asymmetry); break;
        case 2: performChi2Fits_b2b_dihadron(output_file.c_str(), kinematic_file.c_str(), 
        binNames[i], asymmetry); break;
      }
    }
    cout << endl << "     Completed " << binNames[i] << " chi2 fits." << endl;
    // switch (channel) {
    //   case 0: cout << "No MLM fit for inclusive." << endl; break;
    //   case 1: performMLMFits_single_hadron(output_file.c_str(), 
    //     kinematic_file.c_str(), binNames[i]); break;
    //   case 2: performMLMFits_b2b_dihadron(output_file.c_str(), 
    //     kinematic_file.c_str(), binNames[i]); break;
    //   case 3: cout << "No dihadron MLM fit (yet)." << endl; break;
    // }
    cout << endl << "     Completed " << binNames[i] << " MLM fits." << endl;
    cout << endl << endl;
    currentFits++;
  }

  mc_file->Close();
  delete mc_file;
  delete kinematicCuts;
  delete mckinematicCuts;
  kinematicCuts = nullptr;
  mckinematicCuts = nullptr;

  // Stop the timer
  auto end_time = std::chrono::high_resolution_clock::now();
  // Calculate the elapsed time in seconds and microseconds
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - 
    start_time).count();
  double seconds = duration / 1e6;
  // Convert to hours, minutes, and seconds
  int hours = static_cast<int>(seconds) / 3600;
  int remaining_time = static_cast<int>(seconds) % 3600;
  int mins = remaining_time / 60;
  int remaining_seconds = remaining_time % 60;
  // Print the elapsed time
  cout << "Time elapsed: ";
  cout << hours << " hours, " << mins << " mins, " << remaining_seconds << " seconds." << endl;
  return 0;
}