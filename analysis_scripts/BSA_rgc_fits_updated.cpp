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

void load_bins_from_csv(const std::string& filename) {
  // Open the input file with the given filename
  std::ifstream file(filename);

  // Declare a string to store each line read from the file
  std::string line;

  // Declare a boolean flag to check if we have reached the bin declarations in the file
  bool reached_bins = false;

  // Loop through each line in the file until there are no more lines left to read
  while (std::getline(file, line)) {
    // If the line is empty or starts with a '#' (comment), skip to the next line
    if (line.empty() || line[0] == '#') { continue; }

    // If we have not reached the bin declarations yet
    if (!reached_bins) {
      // If the line contains a '-', it marks the start of bin declarations; 
      // set the flag to true and continue to the next line
      if (line.find("-") != std::string::npos) {
        reached_bins = true;
        continue;
      }
      // Use a stringstream to split the line by commas and store variable names
      std::stringstream ss_var(line);
      std::string var_name;
      while (std::getline(ss_var, var_name, ',')) {
        variable_names.push_back(var_name);
      }
    } else {
      // If we have reached the bin declarations, use a stringstream to split the line by commas
      std::stringstream ss(line);

      // Declare strings to store the bin name and property
      std::string bin_name, property;

      // Read the bin name and property from the stringstream
      std::getline(ss, bin_name, ',');
      binNames.push_back(bin_name);
      std::getline(ss, property, ',');
      propertyNames.push_back(property);

      // Declare a vector to store the bin values
      std::vector<float> bin_values;

      // Declare a string to store each value read from the stringstream
      std::string value;

      // Read the values from the stringstream and store them in the bin_values vector
      while (std::getline(ss, value, ',')) {
        bin_values.push_back(std::stof(value));
      }

      // Add the bin_values vector to the bins_map and allBins containers
      bins_map[bin_name] = bin_values;
      allBins.push_back(bin_values);
    }
  }

  // Loop through each variable name to remove newline and carriage return characters
  for (size_t i = 0; i < variable_names.size(); ++i) {
    // Remove newline characters from the variable name
    variable_names[i].erase(std::remove(variable_names[i].begin(),
      variable_names[i].end(), '\n'), variable_names[i].end());
    // Remove carriage return characters from the variable name
    variable_names[i].erase(std::remove(variable_names[i].begin(),
      variable_names[i].end(), '\r'), variable_names[i].end());
  }
}

int main(int argc, char *argv[]) {

  // Check for correct number of command line arguments
    if (argc != 5) {
        std::cout << "Usage: " << argv[0];
        std::cout << " <data_root_file> <mc_root_file> ";
        std::cout << " <output_asymmetry_file> <output_kinematic_file>" << std::endl;
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
  std::cout<< std::endl <<"-- Loaded information from bins.csv. " << std::endl;

  std::cout<< "Found " << allBins.size() << " sets of bins: " << std::endl;
  for (size_t i = 0; i < binNames.size(); ++i) {
    std::cout << binNames[i];
    if (i == binNames.size() - 1) { std::cout << "."; }
    else { std::cout << ", "; }
  }
  std::cout << std::endl;

  std::cout<< "Found " << allBins[currentFits].size() << " bin indices for: " << std::endl;
  for (size_t i = 0; i < allBins[currentFits].size(); ++i) {
    std::cout << allBins[currentFits][i];
    if (i == allBins[currentFits].size() - 1) { std::cout << "."; }
    else { std::cout << ", "; }
  }
  std::cout << std::endl;

  std::cout << "Found " << variable_names.size() << " variables: " << std::endl;
  for (size_t i = 0; i < variable_names.size(); ++i) {
    std::cout << i << ":" << variable_names[i] << std::flush;
    if (i == variable_names.size() - 1) {
      // std::cout << ". ";
    } else {
      std::cout << ", ";
    }
  }
  std::cout << std::endl;

  // load run infrom from external csv file
  load_run_info_from_csv("run_info_rgc.csv");
  std::cout<< std::endl << std::endl <<"-- Loaded information from run_info_rgc.csv" << std::endl;

  std::cout << "Found " << run_info_list.size() << " runs." << std::endl;
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

  std::cout << "Total pos-pos (beam-target) charge: " << cpp << " (nc). ";
  std::cout << "Total pos-neg charge: " << cpm << " (nc). ";
  std::cout << "Total neg-pos charge: " << cmp << " (nc). ";
  std::cout << "Total neg-neg charge: " << cmm << " (nc). ";
  std::cout << "Total unpolarized (carbon) charge: " << total_charge_carbon << " (nc)."<<std::endl;

  // Load data and mc root files
  TFile* gData = new TFile(argv[1], "READ");
  TFile* gMC = new TFile(argv[2], "READ");

  if (!gData->IsOpen() || !gMC->IsOpen()) {
        std::cout << "Error opening ROOT files (is the location correct?). Exiting." << std::endl;
        return 2;
    }

}