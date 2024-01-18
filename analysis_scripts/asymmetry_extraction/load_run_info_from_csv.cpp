// load_run_info_from_csv.cpp

// Include the header file
#include "load_run_info_from_csv.h"

// Standard C++ Library Headers
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

// tbhayward includes
#include "load_run_info_from_csv.h"

// Declare a vector to store the run information
std::vector<RunInfo> run_info_list;

void load_run_info_from_csv(const std::string& filename) {
  // Open the input file with the given filename
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  // Declare a string to store each line read from the file
  std::string line;

  // Loop through each line in the file until there are no more lines left to read
  while (std::getline(file, line)) {
    std::cout << "Reading line: " << line << std::endl;
    // If the line is empty or starts with a '#' (comment), skip to the next line
    if (line.empty() || line[0] == '#') { continue; }

    // Use a stringstream to split the line by commas
    std::stringstream ss(line);

    // Declare a struct to store the run information
    RunInfo run_info;

    // Declare a string to store each piece of information read from the stringstream
    std::string info;

    // Read the run number from the stringstream and convert it to an integer
    std::getline(ss, info, ',');
    // run_info.runnum = std::stoi(info);
    try {
        run_info.runnum = std::stoi(info);
        // ... similarly for std::stof calls ...
    } catch (const std::exception& e) {
        std::cerr << "Exception caught while parsing line: " << line << "\n";
        std::cerr << "Exception message: " << e.what() << "\n";
        continue; // Skip to next line or handle error appropriately
    }

    // Read the total charge from the stringstream and convert it to a double
    std::getline(ss, info, ',');
    run_info.total_charge = std::stof(info);

    // Read the positive charge from the stringstream and convert it to a double
    std::getline(ss, info, ',');
    run_info.positive_charge = std::stof(info);

    // Read the negative charge from the stringstream and convert it to a double
    std::getline(ss, info, ',');
    run_info.negative_charge = std::stof(info);

    // Read the target polarization from the stringstream and convert it to a double
    std::getline(ss, info, ',');
    run_info.target_polarization = std::stof(info);

    // Read the target polarization from the stringstream and convert it to a double
    std::getline(ss, info, ',');
    run_info.target_polarization_uncertainty = std::stof(info);

    // Add the struct to the run_info_list vector
    run_info_list.push_back(run_info);
  }
}