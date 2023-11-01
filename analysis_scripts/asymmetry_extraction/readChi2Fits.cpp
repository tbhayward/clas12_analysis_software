#include "readChi2Fits.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::map<std::string, std::vector<std::vector<double>>> readChi2Fits(const std::string& filepath){
  std::ifstream infile(filepath);
  std::string line;
  std::map<std::string, std::vector<std::vector<double>>> chi2Fits;

  while (std::getline(infile, line)) {
    std::size_t start = line.find("{{") + 2;
    std::size_t end = line.find("}}");
    std::string sub = line.substr(start, end - start);
    std::stringstream ss(sub);
    double mean, value, error;
    char comma;
    ss >> mean >> comma >> value >> comma >> error;
    
    std::size_t eq_pos = line.find("=");
    std::string key = line.substr(0, eq_pos - 1);
    
    std::cout << "Key: " << key << " Mean: " << mean;
    std::cout << " Value: " << value << " Error: " << error << std::endl;
    
    chi2Fits[key].push_back({mean, value, error});  // Store in the map
    
    std::cout << "Inserted into map. Current size: " << chi2Fits.size() << std::endl;
  }

  return chi2Fits;
}