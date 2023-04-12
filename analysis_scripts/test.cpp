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
size_t currentFits = 0;
int n = 1;

std::vector<std::string> variable_names;
std::map<std::string, std::vector<float>> bins_map;
std::vector<std::vector<float>> allBins;
std::vector<std::string> binNames;
std::vector<int> variable_indices;
std::vector<std::string> propertyNames;

string trim_newline(const string &str) {
  if (!str.empty() && str.back() == '\n') {
    return str.substr(0, str.size() - 1);
  }
  return str;
}

void load_bins_from_csv(const std::string& filename) {
  std::ifstream file(filename);
  std::string line;
  bool reached_bins = false; // Flag to check if we have reached the bin declarations

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') { continue; } // Ignore comment lines

    if (!reached_bins) {
      if (line.find("-") != std::string::npos) { // set flag to true and continue to next line
        reached_bins = true;
        continue;
      }
      std::stringstream ss_var(line);
      std::string var_name;
      while (std::getline(ss_var, var_name, ',')) {
        variable_names.push_back(trim_newline(var_name));
      }
    } else {
      std::stringstream ss(line);
      std::string bin_name, property;
      std::getline(ss, bin_name, ',');
      binNames.push_back(trim_newline(bin_name));

      // Retrieve the index of the variable to be used for this bin
      std::string index_str;
      std::getline(ss, index_str, ',');
      int variable_index = std::stoi(index_str);

      property = variable_names[variable_index];
      propertyNames.push_back(property);

      std::vector<float> bin_values;
      std::string value;
      while (std::getline(ss, value, ',')) {
        bin_values.push_back(std::stof(value));
      }
      bins_map[bin_name] = bin_values;
      allBins.push_back(bin_values);
    }
  }
}

struct eventData {
  std::map<std::string, float> data;
};

std::vector<eventData> gData;
size_t currentBin = 0;

eventData parseLine(const std::string& line, const std::vector<std::string>& variable_names) {
  std::istringstream iss(line);
  eventData data;

  float value;
  std::string value_str;
  // size_t var_name_index = 0;
  // for (const auto& var_name : variable_names) {
  //   if (var_name_index == variable_names.size() - 1) {
  //     std::getline(iss, value_str); // Read the remaining value without specifying a delimiter
  //   } else {
  //     std::getline(iss, value_str, ' '); // Use space as the delimiter
  //   }
  //   value_str.erase(std::remove(value_str.begin(), value_str.end(), '\n'), value_str.end());
  //   value = std::stof(value_str);
  //   data.data[var_name] = value;

  //   std::cout << "Var: " << var_name << ", Value_str: " << value_str << ", Value: " << value << std::endl;
  //   var_name_index++;
  // }
  size_t var_name_index = 0;
  for (const auto& var_name : variable_names) {
    cout << var_name << endl;
    if (!(iss >> value)) {
      break;
    }
    data.data[var_name] = value;

    // std::cout << "Var: " << var_name << ", Value: " << value << std::endl;
    var_name_index++;
  }



  // Print the final values of status, runnum, and evnum
  cout << data.data["status"] << " " << data.data["runnum"] << " " << data.data["Delta_phi"] << endl;

  data.data["pol"] = 0.86;
  // Calculate b2b_factor
  const float M = 0.938272088; // proton mass
  float gamma = (2 * M * data.data["x"]) / sqrt(data.data["Q2"]);
  float epsilon = (1 - data.data["y"] - (0.25) * gamma * gamma * data.data["y"] * data.data["y"]) /
    (1 - data.data["y"] + (0.50) * data.data["y"] * data.data["y"] + (0.25) * gamma * gamma *
    data.data["y"] * data.data["y"]);
  float depolarization_factor = sqrt(1 - epsilon * epsilon);
  data.data["b2b_factor"] = (depolarization_factor * data.data["PTPT"]) / (M * M);

  return data;
}

std::vector<eventData> readData(const std::string& filename,
  const std::vector<std::string>& variable_names) {
  std::ifstream infile(filename);
  std::string line;
  std::vector<eventData> data;
  while (std::getline(infile, line)) {
    data.push_back(parseLine(line, variable_names));
  }
  return data;
}

double getEventProperty(const eventData& event, int currentFits) {
  std::string property = propertyNames[currentFits];
  int variable_index = variable_indices[currentFits]; // Use the variable index from the bin
  std::string variable_name = variable_names[variable_index]; // Get the var name using the index
  return event.data.at(variable_name); // Access the property value using the map's indexing
}

// Apply kinematic cuts to the data
bool applyKinematicCuts(const eventData& data, int currentFits) {
  //  if (currentFits <= 4) { return data.data.at("status") <= 1e2; }
  return true;
}

TH1D* createHistogramForBin(const std::vector<eventData>& data, const char* histName,
  int binIndex) {

  double varMin = allBins[currentFits][binIndex];
  double varMax = allBins[currentFits][binIndex + 1];

  TH1D* histPos = new TH1D(Form("%s_pos", histName), "", 24, 0, 2 * TMath::Pi());
  TH1D* histNeg = new TH1D(Form("%s_neg", histName), "", 24, 0, 2 * TMath::Pi());

  double sumPol = 0;
  int numEvents = 0;

  for (const eventData& event : data) {

    double currentVariable = getEventProperty(event, currentFits);
    if (applyKinematicCuts(event, currentFits) && currentVariable >= varMin &&
      currentVariable < varMax) {
      if (event.data.at("helicity") > 0) {
        histPos->Fill(event.data.at("Delta_phi"));
      } else {
        histNeg->Fill(event.data.at("Delta_phi"));
      }
      sumPol += event.data.at("pol");
      numEvents++;
    }
  }

  double meanPol = sumPol / numEvents;

  int numBins = histPos->GetNbinsX();
  TH1D* histAsymmetry = new TH1D(Form("%s_asymmetry", histName), "",
    numBins, 0, 2 * TMath::Pi());

  for (int iBin = 1; iBin <= numBins; ++iBin) {
    double Np = histPos->GetBinContent(iBin);
    double Nm = histNeg->GetBinContent(iBin);
    double asymmetry = (1 / meanPol) * (Np - Nm) / (Np + Nm);
    double error = (2 / meanPol) * std::sqrt(Np * Nm / TMath::Power(Np + Nm, 3));

    histAsymmetry->SetBinContent(iBin, asymmetry);
    histAsymmetry->SetBinError(iBin, error);
  }

  delete histPos;
  delete histNeg;
  return histAsymmetry;
}

// Function to fit
double funcToFit(double* x, double* par) {
  double A = par[0];
  double B = par[1];
  double Delta_phi = x[0];
  return A * sin(Delta_phi) + B * sin(2 * Delta_phi);
}

void performChi2Fits(const char *filename, const char* output_file, const std::string& prefix) {
  gData = readData(filename, variable_names);

  TF1* fitFunction = new TF1("fitFunction", funcToFit, 0, 2 * TMath::Pi(), 2);

  size_t numBins = allBins[currentFits].size() - 1;

  for (size_t i = 0; i < numBins; ++i) {
      cout << "Beginning chi2 fit for " << binNames[currentFits]
        << " bin " << i << ". ";
      char histName[32];
      snprintf(histName, sizeof(histName), "hist_%zu", i);

      TH1D* hist = createHistogramForBin(gData, histName, i);
      hist->Fit(fitFunction, "Q");

      double sumVariable = 0;
      double sumb2b = 0;
      double numEvents = 0;
      for (const eventData& event : gData) {
        double currentVariable = getEventProperty(event, currentFits);
        if (applyKinematicCuts(event, currentFits) && currentVariable >= allBins[currentFits][i] &&
          currentVariable < allBins[currentFits][i + 1]) {
            sumVariable += currentVariable;
            sumb2b += event.data.at("b2b_factor");
            numEvents += 1;
        }
      }
      cout << "Found " << numEvents << " events." << endl;
  }
}

void test(const char* data_file, const char* output_file) {

  // Clear the contents of the output_file
  std::ofstream ofs(output_file, std::ios::trunc);
  ofs.close();

  // load bins from external csv file
  load_bins_from_csv("bins.csv");
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
  cout << endl << endl << endl;


  // cout << endl << endl;
  // for (size_t i = 0; i < allBins.size(); ++i) {
  for (size_t i = 0; i < 1; ++i) {
    cout << "-- Beginning kinematic fits." << endl;
    performChi2Fits(data_file, output_file, binNames[i]);
    cout << endl << "     Completed " << binNames[i] << " chi2 fits." << endl;
    // performMLMFits(data_file, output_file, binNames[i]);
    // cout << endl << "     Completed " << binNames[i] << " MLM fits." << endl;
    cout << endl << endl;
    currentFits++;
  }
  cout << endl << endl;
}