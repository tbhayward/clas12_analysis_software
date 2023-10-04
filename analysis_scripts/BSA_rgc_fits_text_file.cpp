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
#include <chrono>
#include <iomanip> // Include this header for std::fixed and std::setprecision

size_t currentFits = 0;
size_t currentBin = 0;
int n = 1;

std::map<std::string, std::vector<double>> bins_map;
std::vector<std::vector<double>> allBins;
std::vector<std::string> binNames;
std::vector<std::string> propertyNames;
std::vector<std::string> variable_names;
double total_charge_carbon;
double cpp, cpm, cmp, cmm;

string trim_newline(const string &str) {
  if (!str.empty() && str.back() == '\n') {
    return str.substr(0, str.size() - 1);
  }
  return str;
}

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
      std::vector<double> bin_values;

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

struct RunInfo {
  int runnum;
  double total_charge;
  double positive_charge;
  double negative_charge;
  double target_polarization;
  double target_polarization_uncertainty;
};

// Declare a vector to store the run information
std::vector<RunInfo> run_info_list;

void load_run_info_from_csv(const std::string& filename) {
  // Open the input file with the given filename
  std::ifstream file(filename);

  // Declare a string to store each line read from the file
  std::string line;

  // Loop through each line in the file until there are no more lines left to read
  while (std::getline(file, line)) {
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
    run_info.runnum = std::stoi(info);

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

// function to get the polarization value
double getPol(int runnum) {
  double pol = 0.86; 
    if (runnum == 11 ) { pol = 0.86; } // runnum == 11 indicates Monte Carlo in CLAS12
    else if (runnum >= 5032 && runnum < 5333) { pol = 0.8592; } 
    else if (runnum >= 5333 && runnum <= 5666) { pol = 0.8922; }
    else if (runnum >= 6616 && runnum <= 6783) { pol = 0.8453; }
    else if (runnum >= 6142 && runnum <= 6149) { pol = 0.81132; }
    else if (runnum >= 6150 && runnum <= 6188) { pol = 0.82137; }
    else if (runnum >= 6189 && runnum <= 6260) { pol = 0.83598; }
    else if (runnum >= 6261 && runnum <= 6339) { pol = 0.80770; }
    else if (runnum >= 6340 && runnum <= 6342) { pol = 0.85536; }
    else if (runnum >= 6344 && runnum <= 6399) { pol = 0.87038; }
    else if (runnum >= 6420 && runnum <= 6476) { pol = 0.88214; }
    else if (runnum >= 6479 && runnum <= 6532) { pol = 0.86580; }
    else if (runnum >= 6533 && runnum <= 6603) { pol = 0.87887; }
    else if (runnum >= 11013 && runnum <= 11309) { pol = 0.84983; }
    else if (runnum >= 11323 && runnum <= 11334) { pol = 0.87135; }
    else if (runnum >= 11335 && runnum <= 11387) { pol = 0.85048; }
    else if (runnum >= 11389 && runnum <= 11571) { pol = 0.84262; }
    else if (runnum >= 16000) { pol = 0.83534; } // RGC +/- 0.01440
  return pol;
}

struct eventData {
  std::unordered_map<std::string, double> data;
};

std::vector<eventData> gData;
std::vector<eventData> gMC;

eventData parseLine(const std::string& line, const std::vector<std::string>& variable_names) {
  // Create a stringstream from the input line
  std::istringstream iss(line);

  // Initialize an eventData object to store the parsed data
  eventData data;

  // Declare value for storing the extracted double and an index for variable names
  double value;
  std::string value_str;
  size_t var_name_index = 0;

  // Iterate through the provided variable names
  for (const auto& var_name : variable_names) {
    // Attempt to extract a double value from the stringstream
    if (!(iss >> value)) {
      break;
    }
    // Insert the extracted value into the data map with the corresponding variable name
    data.data.emplace(var_name, value);

    // Increment the variable name index
    var_name_index++;
  }

  // Calculate the polarization value and store it in the data map
  data.data["pol"] = getPol(data.data["runnum"]);

  // Get the target polarization value from the run_info_list and store it in the data map
  int runnum = static_cast<int>(data.data["runnum"]);
  if (data.data["runnum"] == 11) { data.data["target_pol"] = 0; } // MC
  else { for (const auto& run_info : run_info_list) {
      if (run_info.runnum == runnum) {
        data.data["target_pol"] = run_info.target_polarization;
        break;
      }
    }
  }

  // Add in tmin and t
  double me = 0.000511; // electron mass in GeV
  double mp = 0.938272; // proton mass in GeV
  // double E = 10.55; // beam energy
  double E = mp;
  // double Ep = data.data["e_p"]; // scattered electron momentum 
  // double theta = data.data["e_theta"]; // scattered electron angle
  double Ep = data.data["p_p"]; // scattered proton momentum 
  double theta = data.data["p_theta"]; // scattered proton angle 
  data.data["tmin"] = -pow((mp*data.data["x"]),2)/(1-data.data["x"]);
  // data.data["t"] = 2*me*(E - Ep) - 2*sqrt(me*me + E*E)*sqrt(me*me + Ep*Ep) +
  //         2*sqrt(me*me + E*E)*sqrt(me*me + Ep*Ep)*cos(theta);
  data.data["t"] = 2*mp*(E - Ep) - 2*sqrt(mp*mp + E*E)*sqrt(mp*mp + Ep*Ep) +
          2*sqrt(mp*mp + E*E)*sqrt(mp*mp + Ep*Ep)*cos(theta);

  // Return the populated eventData object
  return data;
}


std::vector<eventData> readData(const std::string& filename,
  const std::vector<std::string>& variable_names) {
  // Open the input file using the provided filename
  std::ifstream infile(filename);

  // Count the number of lines in the file
  size_t numberOfLines = std::count(std::istreambuf_iterator<char>(infile),
    std::istreambuf_iterator<char>(), '\n');

  // Reset the file stream to the beginning
  infile.clear();
  infile.seekg(0, std::ios::beg);

  // Declare a string to store each line and a vector to store the parsed eventData objects
  std::string line;
  std::vector<eventData> data;

  // Read the input file line by line
  while (std::getline(infile, line)) {
    // Parse the line using the parseLine function and the provided variable names
    eventData parsedData = parseLine(line, variable_names);

    // Add the parsed eventData object to the data vector
    data.push_back(parsedData);
  }

  // Close the input file
  infile.close();

  // Return the vector of eventData objects
  return data;
}

double getEventProperty(const eventData& event, int currentFits) {

  std::string property = propertyNames[currentFits];
  // Access the property value using the map's indexing
  return event.data.at(property);
}

// Apply kinematic cuts to the data
bool applyKinematicCuts(const eventData& data, int currentFits, bool isMC) {

  bool goodEvent = 0;
  std::string property = binNames[currentFits];
  //
  //
  //
  // epX
  if (property == "xF") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.4 &&
      data.data.at("y")<0.75;
  }
  if (property == "Mx") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("y")<0.75;
  }
  if (property == "Q2bin") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.4 &&
      data.data.at("y")<0.75 && data.data.at("x")>0.2 && data.data.at("x")<0.3 &&
      data.data.at("pT")>0.25 && data.data.at("pT")<0.35 && data.data.at("xF")<0;
  }
  if (property == "PTTFR" || property ==  "xTFR" || property == "zetaTFR" || 
    property == "Q2TFR" || property ==  "x") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.4 &&
      data.data.at("y")<0.75 && data.data.at("xF")<0;
  }
  if (property == "PTCFR" || property == "xCFR" || property == "zetaCFR" ||
    property == "Q2TFR") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.4 &&
      data.data.at("y")<0.75 && data.data.at("xF")>0;
  } 
  //
  //
  //
  // epi+X
  if (property == "xFpip") { 
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.5 &&
      data.data.at("y")<0.75;
  }
  if (property == "PTTFRpip" || property ==  "xTFRpip" || property == "zTFRpip" || 
    property == "Q2TFRpip" || property ==  "xpip") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.5 &&
      data.data.at("y")<0.75 && data.data.at("xF")<0;
  }
  if (property == "PTCFRpip" || property == "xCFRpip" || property == "zCFRpip" ||
    property == "Q2TFRpip") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.5 &&
      data.data.at("y")<0.75 && data.data.at("xF")>0;
  }
  //
  //
  //
  // epi-X
  if (property == "xFpim") { 
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.5 &&
      data.data.at("y")<0.75;
  }
  if (property == "PTTFRpim" || property ==  "xTFRpim" || property == "zTFRpim" || 
    property == "Q2TFRpim" || property ==  "xpim") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.5 &&
      data.data.at("y")<0.75 && data.data.at("xF")<0;
  }
  if (property == "PTCFRpim" || property == "xCFRpim" || property == "zCFRpim" ||
    property == "Q2TFRpim") {
    goodEvent = data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.5 &&
      data.data.at("y")<0.75 && data.data.at("xF")>0;
  } 
  if (isMC) { return goodEvent; }
  else {return goodEvent && data.data.at("target_pol") != 0; } // if data, skip Pt = 0 (carbon)

  return goodEvent;  
}

double dilution_factor(double currentVariable, const std::string& prefix) {

  // epX
  if (prefix == "xF") { 
    return 0.186121-0.0263337*currentVariable-0.175587*std::pow(currentVariable,2)+
      0.0522814*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "xFpip") { 
    return 0.122453+0.189509*currentVariable-0.133621*std::pow(currentVariable,2)-
      0.0401427*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "xFpim") { 
    return 0.128348+0.195055*currentVariable-0.242617*std::pow(currentVariable,2)-
      0.204807*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "Mx") { 
    return 0.0847657+0.0762168*currentVariable-0.0128988*std::pow(currentVariable,2)+
      0.00274429*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "Q2TFR") {
    return 0.0884319+0.0414953*currentVariable-0.00584857*std::pow(currentVariable,2)+
      0.000500127*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "Q2bin") {
    return -0.341032+0.762811*currentVariable-0.399944*std::pow(currentVariable,2)+
      0.0686534*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "xTFR") {
    return 0.111702+0.0858432*currentVariable+0.880331*std::pow(currentVariable,2)-
      0.990298*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "xTFRpip") {
    return 0.117706-0.194421*currentVariable+0.977489*std::pow(currentVariable,2)-
      0.926193*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "xTFRpim") {
    return 0.0787795-0.263136*currentVariable+1.378*std::pow(currentVariable,2)-
      1.65335*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "PTTFR") {
    return 0.184491-0.161007*currentVariable+0.298733*std::pow(currentVariable,2)-
      0.187826*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "PTTFRpip") {
    return 0.176079-0.328598*currentVariable+0.475598*std::pow(currentVariable,2)-
      0.167004*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "PTTFRpim") {
    return 0.0275594+0.276354*currentVariable-0.471179*std::pow(currentVariable,2)-
      0.21497*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "zetaTFR") {
    return 1.52544-7.07359*currentVariable+12.5954*std::pow(currentVariable,2)-
      7.72548*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "zTFRpip") {
    return 0.0565765+0.882732*currentVariable-3.33409*std::pow(currentVariable,2)+
      5.51154*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "zTFRpim") {
    return -0.0253779+1.62183*currentVariable-6.76455*std::pow(currentVariable,2)+
      8.56005*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "Q2TFR") {
    return 0.093586+0.0370678*currentVariable-0.00373394*std::pow(currentVariable,2)+
      0.000215739*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "xCFR") {
    return 0.089331+0.429008*currentVariable-0.617364*std::pow(currentVariable,2)+
      0.7584*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "xCFRpip") {
    return 0.119971+0.416041*currentVariable-0.922544*std::pow(currentVariable,2)+
      1.01908*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "xCFRpim") {
    return 0.121553-0.12187*currentVariable+0.923064*std::pow(currentVariable,2)-
      0.949773*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "PTCFR") {
    return 0.151263+0.170759*currentVariable-0.439815*std::pow(currentVariable,2)+
      0.278509*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "PTCFRpip") {
    return 0.184542-0.0499585*currentVariable+0.163844*std::pow(currentVariable,2)+
      0.157106*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "PTCFRpim") {
    return 0.147254-0.134125*currentVariable+0.407317*std::pow(currentVariable,2)-
      0.339619*std::pow(currentVariable,3);
  }
  // epX
  if (prefix == "zetaCFR") {
    return 1.32783-6.22826*currentVariable+11.2985*std::pow(currentVariable,2)-
      7.01171*std::pow(currentVariable,3);
  }
  // epi+X
  if (prefix == "zCFRpip") {
    return 1.32783-6.22826*currentVariable+11.2985*std::pow(currentVariable,2)-
      7.01171*std::pow(currentVariable,3);
  }
  // epi-X
  if (prefix == "zCFRpim") {
    return 0.0997319+0.28069*currentVariable-0.547782*std::pow(currentVariable,2)+
      0.244802*std::pow(currentVariable,3);
  }
  return 0.14;
}

// Negative log-likelihood function
void negLogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    // npar: number of parameters
    // gin: an array of derivatives (if needed)
    // f: the value of the function
    // par: an array of the parameter values
    // iflag: a flag (see TMinuit documentation for details)

    // Extract parameters from the input parameter array
    double ALU_sinphi = par[0];
    double AUL_sinphi = par[1];
    double AUL_sin2phi = par[2];
    double ALL = par[3];
    double ALL_cosphi = par[4];
    double AUU_cosphi = par[5];
    double AUU_cos2phi = par[6];

    // Initialize variables for counting events (N), positive helicity sum (sum_P), 
    // and negative helicity sum (sum_N)
    double N = 0;
    double NUU = 0; // normalization integral
    double sum_PP = 0; // positive beam -- positive target
    double sum_PM = 0; // positive beam -- negative target
    double sum_MP = 0; // negative beam -- positive target
    double sum_MM = 0; // negative beam -- negative target

    // Iterate through the global event data (gData)
    for (const eventData &event : gData) {
        // Get the value of the current variable of interest for the event
        double currentVariable = getEventProperty(event, currentFits);

        // Apply kinematic cuts and check if the current variable is within the specified bin range
        if (applyKinematicCuts(event, currentFits, 0) && 
          currentVariable >= allBins[currentFits][currentBin] && 
          currentVariable < allBins[currentFits][currentBin + 1]) {

          // Increment the event count
          N += 1;

          // Extract phi and polarization (pol) from the event data
          double phi = event.data.at("phi");
          double Pb = event.data.at("pol");
          double Pt = std::abs(event.data.at("target_pol"));
          double DepA = event.data.at("DepA");
          double DepB = event.data.at("DepB");
          double DepC = event.data.at("DepC");
          double DepV = event.data.at("DepV");
          double DepW = event.data.at("DepW");

          double Df = dilution_factor(currentVariable, binNames[currentFits]); // dilution factor
          // Check if the helicities are positive or negative and update the corresponding sum
          if (event.data.at("helicity") > 0 && event.data.at("target_pol") > 0) {
            sum_PP += log(1 
              + (DepV/DepA)*AUU_cosphi*cos(phi) + (DepB/DepA)*AUU_cos2phi*cos(2*phi) // UU 
              + Pb*((DepW/DepA)*ALU_sinphi*sin(phi)) // BSA
              + Df*Pt*((DepV/DepA)*AUL_sinphi*sin(phi) + (DepB/DepA)*AUL_sin2phi*sin(2*phi)) // TSA
              + Df*Pb*Pt*((DepC/DepA)*ALL + (DepW/DepA)*ALL_cosphi*cos(phi)) ); // DSA
          } else if (event.data.at("helicity") > 0 && event.data.at("target_pol") < 0 ) {
            sum_PM += log(1 
              + (DepV/DepA)*AUU_cosphi*cos(phi) + (DepB/DepA)*AUU_cos2phi*cos(2*phi) // UU
              + Pb*((DepW/DepA)*ALU_sinphi*sin(phi)) // BSA
              - Df*Pt*((DepV/DepA)*AUL_sinphi*sin(phi) + (DepB/DepA)*AUL_sin2phi*sin(2*phi)) // TSA
              - Df*Pb*Pt*((DepC/DepA)*ALL + (DepW/DepA)*ALL_cosphi*cos(phi)) ); // DSA
          } else if (event.data.at("helicity") < 0 && event.data.at("target_pol") > 0 ) {
            sum_MP += log(1
              + (DepV/DepA)*AUU_cosphi*cos(phi) + (DepB/DepA)*AUU_cos2phi*cos(2*phi) // UU
              - Pb*((DepW/DepA)*ALU_sinphi*sin(phi)) // BSA
              + Df*Pt*((DepV/DepA)*AUL_sinphi*sin(phi) + (DepB/DepA)*AUL_sin2phi*sin(2*phi)) // TSA
              - Df*Pb*Pt*((DepC/DepA)*ALL + (DepW/DepA)*ALL_cosphi*cos(phi)) ); // DSA
          } else if (event.data.at("helicity") < 0 && event.data.at("target_pol") < 0 ) {
            sum_MM += log(1 
              + (DepV/DepA)*AUU_cosphi*cos(phi) + (DepB/DepA)*AUU_cos2phi*cos(2*phi) // UU
              - Pb*((DepW/DepA)*ALU_sinphi*sin(phi)) // BSA
              - Df*Pt*((DepV/DepA)*AUL_sinphi*sin(phi) + (DepB/DepA)*AUL_sin2phi*sin(2*phi)) // TSA
              + Df*Pb*Pt*((DepC/DepA)*ALL + (DepW/DepA)*ALL_cosphi*cos(phi)) ); // DSA
          }
        }
    }

    // Iterate through the global event mc (gMC)
    for (const eventData &event : gMC) {
        // Get the value of the current variable of interest for the event
        double currentVariable = getEventProperty(event, currentFits);

        // Apply kinematic cuts and check if the current variable is within the specified bin range
        if (applyKinematicCuts(event, currentFits, 1) && 
          currentVariable >= allBins[currentFits][currentBin] && 
          currentVariable < allBins[currentFits][currentBin + 1]) {

          // Extract Delta_phi and polarization (pol) from the event data
          double phi = event.data.at("phi");
          double DepA = event.data.at("DepA");
          double DepB = event.data.at("DepB");
          double DepC = event.data.at("DepC");
          double DepV = event.data.at("DepV");
          double DepW = event.data.at("DepW");

          NUU+= 1 + (DepV/DepA)*AUU_cosphi*cos(phi) + (DepB/DepA)*AUU_cos2phi*cos(2*phi); // UU
        }
    }

    // determine min pos or neg beam helicity accumulated charge to scale down higher one
    double minBeamCharge = std::min({(cpp+cpm),(cmp+cmm)}); 
    // determine min pos or neg target helicity accumulated charge to scale down higher one
    double minTargetCharge = std::min({(cpp+cmp),(cpm+cmm)}); 

    double nll = N * log(NUU) - 
      minBeamCharge*minTargetCharge/((cpp+cpm)*(cpp+cmp))*sum_PP -
      minBeamCharge*minTargetCharge/((cpp+cpm)*(cpm+cmm))*sum_PM - 
      minBeamCharge*minTargetCharge/((cmp+cmm)*(cpp+cmp))*sum_MP - 
      minBeamCharge*minTargetCharge/((cmp+cmm)*(cpm+cmm))*sum_MM;
    cout << "On MLM fit " << binNames[currentFits] << " " << currentFits << ", " << nll << endl;
    cout << "AUU_cosphi = " << AUU_cosphi << ", AUU_cos2phi = " << AUU_cos2phi;
    cout << ", ALU_sinphi = " << ALU_sinphi;
    cout << ", AUL_sinphi = " << AUL_sinphi << ", AUL_sin2phi = " << AUL_sin2phi;
    cout << ", ALL = " << ALL << ", ALL_cosphi = " << ALL_cosphi << "." << endl;
    // Calculate the negative log-likelihood value and store it in the output variable f
    f = nll;
}

void performMLMFits(const char *filename, const char* output_file, const char* kinematic_file,
  const std::string& prefix) {
  // Read the event data from the input file and store it in the global variable gData

  // Determine the number of bins
  size_t numBins = allBins[currentFits].size() - 1;

  // Initialize TMinuit
  double arglist[10]; arglist[0] = 1;
  int ierflg = 0;
  TMinuit minuit(7); // parameter numbers
  minuit.SetPrintLevel(-1);
  minuit.SetErrorDef(0.5); // error definition for MLE, 1 for chi2
  // This is due to the fact that −logL = chi2/2. 
  // The default value of ErrorDef=1 corresponds to one standard deviation for chi2 function.
  minuit.SetFCN(negLogLikelihood);

  // Declare string streams for storing the MLM fit results
  std::ostringstream mlmFitsAStream; std::ostringstream mlmFitsBStream; 
  std::ostringstream mlmFitsCStream; std::ostringstream mlmFitsDStream; 
  std::ostringstream mlmFitsEStream; std::ostringstream mlmFitsFStream;
  std::ostringstream mlmFitsGStream; 

  // Initialize the string streams with the output variable names
  mlmFitsAStream << prefix << "MLMFitsALUsinphi = {";
  mlmFitsBStream << prefix << "MLMFitsAULsinphi = {";
  mlmFitsCStream << prefix << "MLMFitsAULsin2phi = {";
  mlmFitsDStream << prefix << "MLMFitsALL = {";
  mlmFitsEStream << prefix << "MLMFitsALLcosphi = {";
  mlmFitsFStream << prefix << "MLMFitsAUUcosphi = {";
  mlmFitsGStream << prefix << "MLMFitsAUUcos2phi = {";

  // Initialize string streams to store the mean variables for each bin and asymmetries
  std::ostringstream asymmetryStream;
  asymmetryStream << "\\begin{table}[h]" << std::endl;
  asymmetryStream << "\\centering" << std::endl;
  asymmetryStream << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|} \\hline" << std::endl;
  asymmetryStream << "Bin & $<" << prefix << ">$ & $F_{UU}^{\\cos(\\phi)}/F_{UU}$ & ";
  asymmetryStream << "$F_{UU}^{\\cos(2\\phi)}/F_{UU}$ ";
  asymmetryStream << "& $F_{LU}^{\\sin(\\phi)}/F_{UU}$ & $F_{UL}^{\\sin(\\phi)}/F_{UU}$ & ";
  asymmetryStream << "$F_{UL}^{\\sin(2\\phi)}/F_{UU}$ & $F_{LL}/F_{UU}$ &";
  asymmetryStream << "$F_{LL}^{\\cos(\\phi)}/F_{UU}$ \\\\ \\hline" << std::endl;

  // Iterate through each bin
  for (size_t i = 0; i < numBins; ++i) {
    cout << endl << "Beginning MLM fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    currentBin = i;

    // Define the parameters with initial values and limits
    minuit.DefineParameter(0, "ALU_sinphi", -0.015, 0.01, -1, 1);
    minuit.DefineParameter(1, "AUL_sinphi", -0.020, 0.01, -1, 1);
    minuit.DefineParameter(2, "AUL_sin2phi", -0.010, 0.01, -1, 1);
    minuit.DefineParameter(3, "ALL", 0.40, 0.01, -1, 1);
    minuit.DefineParameter(4, "ALL_cosphi", 0.01, 0.01, -1, 1);
    minuit.DefineParameter(5, "AUU_cosphi", -0.1, 0.01, -1, 1);
    minuit.DefineParameter(6, "AUU_cos2phi", 0.10, 0.01, -1, 1);

    // Minimize the negative log-likelihood function
    minuit.Migrad(); cout << endl;
    // minuit.mnexcm("MINImize",arglist,1,ierflg);
    // minuit.mnexcm("MINOS",arglist,1,ierflg); // running MINOS and HESSE recommended
    // minuit.mnexcm("HESSE",arglist,0,ierflg);

    // Extract the fitted parameter values and errors
    double ALU_sinphi, ALU_sinphi_error;
    minuit.GetParameter(0, ALU_sinphi, ALU_sinphi_error);
    double AUL_sinphi, AUL_sinphi_error;
    minuit.GetParameter(1, AUL_sinphi, AUL_sinphi_error);
    double AUL_sin2phi, AUL_sin2phi_error;
    minuit.GetParameter(2, AUL_sin2phi, AUL_sin2phi_error);
    double ALL, ALL_error;
    minuit.GetParameter(3, ALL, ALL_error);
    double ALL_cosphi, ALL_cosphi_error;
    minuit.GetParameter(4, ALL_cosphi, ALL_cosphi_error);
    double AUU_cosphi, AUU_cosphi_error;
    minuit.GetParameter(5, AUU_cosphi, AUU_cosphi_error);
    double AUU_cos2phi, AUU_cos2phi_error;
    minuit.GetParameter(6, AUU_cos2phi, AUU_cos2phi_error);

    // Calculate the mean values of the current variable 
    double sumVariable = 0;
    double numEvents = 0;
    for (const eventData &event : gData) {
      double currentVariable = getEventProperty(event, currentFits);
        if (applyKinematicCuts(event, currentFits, 0) && currentVariable >= 
          allBins[currentFits][i] && currentVariable < allBins[currentFits][i + 1]) {
            sumVariable += currentVariable;
            numEvents += 1;
        }
    }
    double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;

    // output to text file
    mlmFitsAStream << "{" << meanVariable << ", " << ALU_sinphi << ", " << ALU_sinphi_error << "}";
    mlmFitsBStream << "{" << meanVariable << ", " << AUL_sinphi << ", " << AUL_sinphi_error << "}";
    mlmFitsCStream << "{" << meanVariable << ", " << AUL_sin2phi << ", "<<AUL_sin2phi_error << "}";
    mlmFitsDStream << "{" << meanVariable << ", " << ALL << ", " << ALL_error << "}";
    mlmFitsEStream << "{" << meanVariable << ", " << ALL_cosphi << ", "<<ALL_cosphi_error << "}";
    mlmFitsFStream << "{" << meanVariable << ", " << AUU_cosphi << ", "<<AUU_cosphi_error << "}";
    mlmFitsGStream << "{" << meanVariable << ", " << AUU_cos2phi << ", "<<AUU_cos2phi_error << "}";

    if (i < numBins - 1) {
        mlmFitsAStream << ", "; mlmFitsBStream << ", "; mlmFitsCStream << ", ";
        mlmFitsDStream << ", "; mlmFitsEStream << ", "; mlmFitsFStream << ", "; 
        mlmFitsGStream << ", ";
    }

    // outputs of asymmetries for LaTeX tables
    // Set fixed-point notation and one digit past the decimal
    asymmetryStream << std::fixed << std::setprecision(2); 
    asymmetryStream << (i+1) << " & " << meanVariable << " & ";
    // AUU cosphi
    asymmetryStream << "$" << 100*AUU_cosphi << "_{" << TMath::Abs(100*0.5*AUU_cosphi) << "}^{";
    asymmetryStream << 100*AUU_cosphi_error << "}$ &";
    // AUU cos2phi
    asymmetryStream << "$" << 100*AUU_cos2phi << "_{" << TMath::Abs(100*0.5*AUU_cos2phi) << "}^{";
    asymmetryStream << 100*AUU_cos2phi_error << "}$ &";
    // ALU sinphi
    asymmetryStream << "$" << 100*ALU_sinphi << "_{" << TMath::Abs(100*0.068*ALU_sinphi) << "}^{";
    asymmetryStream << 100*ALU_sinphi_error << "}$ &";
    // AUL sinphi
    asymmetryStream << "$" << 100*AUL_sinphi << "_{" << TMath::Abs(100*0.092*AUL_sinphi) << "}^{";
    asymmetryStream << 100*AUL_sinphi_error << "}$ &";
    // AUL sin2phi
    asymmetryStream << "$" << 100*AUL_sin2phi << "_{" << TMath::Abs(100*0.092*AUL_sin2phi) << "}^{";
    asymmetryStream << 100*AUL_sin2phi_error << "}$ &";
    // ALL 
    asymmetryStream << "$" << 100*ALL << "_{" << TMath::Abs(100*0.097*ALL) << "}^{";
    asymmetryStream << 100*ALL_error << "}$ &";
    // ALL cosphi
    asymmetryStream << "$" << 100*ALL_cosphi << "_{" << TMath::Abs(100*0.097*ALL_cosphi) << "}^{";
    asymmetryStream << 100*ALL_cosphi << "}$";
    asymmetryStream << std::string(" \\\\ \\hline ");
  }

  mlmFitsAStream << "};"; mlmFitsBStream << "};"; mlmFitsCStream << "};";
  mlmFitsDStream << "};"; mlmFitsEStream << "};"; mlmFitsFStream << "};"; 
  mlmFitsGStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << mlmFitsAStream.str() << std::endl;
  outputFile << mlmFitsBStream.str() << std::endl;
  outputFile << mlmFitsCStream.str() << std::endl;
  outputFile << mlmFitsDStream.str() << std::endl;
  outputFile << mlmFitsEStream.str() << std::endl;
  outputFile << mlmFitsFStream.str() << std::endl;
  outputFile << mlmFitsGStream.str() << std::endl;

  outputFile.close();


  // Finally, close the table
  asymmetryStream << "\\end{tabular}" << std::endl;
  asymmetryStream << "\\caption{The mean kinematic value and the final ";
  asymmetryStream << "extracted structure function ratios for " << prefix;
  asymmetryStream << ". Asymmetries are given as ";
  asymmetryStream << "$100{A}_{\\pm100\\Delta\\text{sys}}^";
  asymmetryStream << "{\\pm100\\Delta\\text{stat}}$.}" << std::endl;
  asymmetryStream << "\\label{table:kinematics_" << prefix << "}" << std::endl;
  asymmetryStream << "\\end{table}" << std::endl;
  asymmetryStream << endl << endl << endl;
  std::ofstream kinematicFile(kinematic_file, std::ios_base::app);
  // Write the string stream content to the file
  kinematicFile << asymmetryStream.str() << std::endl; 
  kinematicFile.close();
}

double asymmetry_value_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, double Ptp, double Ptm, 
  int asymmetry_index) {
  double Df = dilution_factor(currentVariable, prefix); // dilution factor
  // return the asymmetry value 
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      return (1 / meanPol) * (Ptm*(Npp-Nmp)+Ptp*(Npm-Nmm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    case 1: // target-spin asymmetry
      return (1 / Df) * ((Npp+Nmp)-(Npm+Nmm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    case 2: // double-spin asymmetry
      return (1 / (Df*meanPol)) * ((Npp-Nmp)+(Nmm-Npm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    default:
      cout << "Invalid asymmetry_index!" << endl;
      return 0;
  }
}

double asymmetry_error_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, double Ptp, double Ptm, 
  int asymmetry_index) {
  double Df = dilution_factor(currentVariable, prefix); // dilution factor
  // return the asymmetry error 
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      return (2 / meanPol) * std::sqrt(
        ((cmm*cpm*cpp*Nmp*std::pow(Ptm,2)*std::pow(Npp*Ptm+Npm*Ptp,2))+
        (cmp*cpm*cpp*Nmm*std::pow(Ptp,2)*std::pow(Npp*Ptm+Npm*Ptp,2))+
        (cmm*cmp*std::pow(Nmp*Ptm+Nmm*Ptp,2)*(cpm*Npp*std::pow(Ptm,2)+cpp*Npm*std::pow(Ptp,2))))/
        (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
    case 1: // target-spin asymmetry
      return (1 / Df) * std::sqrt(
        (((cmp*cpm*cpp*Nmm*std::pow(Nmp+Npp,2)+cmm*cmp*cpp*Npm*std::pow(Nmp+Npp,2)+
        cmm*cpm*std::pow(Nmm+Npm,2)*(cpp*Nmp+cmp*Npp))*std::pow(Ptm+Ptp,2))) /
        (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
    case 2: // double-spin asymmetry
      return (1 / (Df*meanPol)) * std::sqrt(
        (cmp*cpm*cpp*Nmm*std::pow((Nmp+Npp)*Ptm+(Nmp+2*Npm-Npp)*Ptp,2) + 
        cmm*cmp*cpp*Npm*std::pow(Nmp*(Ptm-Ptp)+2*Nmm*Ptp+Npp*(Ptm+Ptp),2) +
        cmm*cpm*(cmp*Npp*std::pow((-Nmm+2*Nmp+Npm)*Ptm+(Nmm+Npm)*Ptp,2) +
        cpp*Nmp*std::pow((Nmm-Npm+2*Npp)*Ptm+(Nmm+Npm)*Ptp,2))) / 
        (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
    default:
      cout << "Invalid asymmetry_index!" << endl;
      return 0;
  }
}

TH1D* createHistogramForBin(const std::vector<eventData>& data, const char* histName,
  int binIndex, const std::string& prefix, int asymmetry_index) {
  // Determine the variable range for the specified bin
  double varMin = allBins[currentFits][binIndex];
  double varMax = allBins[currentFits][binIndex + 1];

  // Create positive and negative helicity histograms
  TH1D* histPosPos = new TH1D(Form("%s_pospos", histName), "", 12, 0, 2 * TMath::Pi());
  TH1D* histPosNeg = new TH1D(Form("%s_posneg", histName), "", 12, 0, 2 * TMath::Pi());
  TH1D* histNegPos = new TH1D(Form("%s_negpos", histName), "", 12, 0, 2 * TMath::Pi());
  TH1D* histNegNeg = new TH1D(Form("%s_negneg", histName), "", 12, 0, 2 * TMath::Pi());

  // Initialize variables to store the sums and event counts
  double sumVariable = 0;
  double numEvents = 0;
  // Variables to calculate the mean polarization
  double sumPol = 0; // sum of the beam polarization
  double test = 0;
  double sumTargetPosPol = 0; // sum of the target positive polarization
  double sumTargetNegPol = 0; // sum of the target negative polarization
  int numEventsPosTarget = 0;
  int numEventsNegTarget = 0;

  // Fill the positive and negative helicity histograms
  for (const eventData& event : data) {

    double currentVariable = getEventProperty(event, currentFits);
    if (applyKinematicCuts(event, currentFits, 0) && currentVariable >= varMin && 
      currentVariable < varMax) {
      sumVariable+=currentVariable;

      if (event.data.at("helicity") > 0 && event.data.at("target_pol") > 0) {
        histPosPos->Fill(event.data.at("phi"));
      } else if (event.data.at("helicity") > 0 && event.data.at("target_pol") < 0) {
        histPosNeg->Fill(event.data.at("phi"));
      } else if (event.data.at("helicity") < 0 && event.data.at("target_pol") > 0) {
        histNegPos->Fill(event.data.at("phi"));
      } else if (event.data.at("helicity") < 0 && event.data.at("target_pol") < 0) {
        histNegNeg->Fill(event.data.at("phi"));
      }

      // Accumulate polarization and event count for mean polarization calculation
      sumPol += event.data.at("pol");
      test += 0.83534;
      if (event.data.at("target_pol") > 0) {
        sumTargetPosPol+=event.data.at("target_pol");
        numEventsPosTarget++;
      } else if (event.data.at("target_pol") < 0) {
        sumTargetNegPol+=event.data.at("target_pol");
        numEventsNegTarget++;
      }
      numEvents++;
    }
  }

  cout << endl << test << " " << sumPol;
  // Calculate the mean polarization
  double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
  double meanPol = sumPol / numEvents; // mean beam polarization for data 
  double Ptp = sumTargetPosPol / numEventsPosTarget;// mean positive target polarization for data
  double Ptm = - sumTargetNegPol / numEventsNegTarget;// mean negative target polarization for data
  // the negative sign here is correct; RGC lists the polarizations with signs to tell which is 
  // which but the polarization really should just be "percent of polarized nucleii"
  cout << endl;
  cout << numEvents << " " << sumVariable << " " << meanVariable << endl;
  cout << numEvents << " " << sumPol << " " << meanPol << endl;
  cout << numEvents << " " << " " << numEventsPosTarget << " " << sumTargetPosPol << " " << Ptp << endl;
  cout << numEvents << " " << " " << numEventsNegTarget << " " << sumTargetNegPol << " " << Ptm << endl;


  // Create the asymmetry histogram
  int numBins = histPosPos->GetNbinsX();
  TH1D* histAsymmetry = new TH1D(Form("%s_asymmetry", histName), "", 
    numBins, 0, 2 * TMath::Pi());

  // Calculate the asymmetry and its error for each bin, and fill the asymmetry histogram
  for (int iBin = 1; iBin <= numBins; ++iBin) {
    double Npp = histPosPos->GetBinContent(iBin)/cpp;
    double Npm = histPosNeg->GetBinContent(iBin)/cpm;
    double Nmp = histNegPos->GetBinContent(iBin)/cmp;
    double Nmm = histNegNeg->GetBinContent(iBin)/cmm;

    // Calculate the asymmetry and error for the current bin
    double asymmetry = asymmetry_value_calculation(meanVariable, prefix, Npp, Npm, Nmp, Nmm, 
      meanPol, Ptp, Ptm, asymmetry_index);
    double error = asymmetry_error_calculation(meanVariable, prefix, Npp, Npm, Nmp, Nmm, meanPol, 
      Ptp, Ptm, asymmetry_index);

    // Fill the asymmetry histogram with the calculated values
    histAsymmetry->SetBinContent(iBin, asymmetry);
    histAsymmetry->SetBinError(iBin, error);
  }

  // Delete the temporary positive and negative helicity histograms
  delete histPosPos;
  delete histPosNeg;
  delete histNegPos;
  delete histNegNeg;

  // Return the final asymmetry histogram
  return histAsymmetry;
}

// Function to fit the beam-spin asymmetry histogram
double BSA_funcToFit(double* x, double* par) {
  // Retrieve the parameters 
  double ALU_offset = par[0];
  double ALU_sinphi = par[1];
  // double AUU_cosphi = par[2];
  // double AUU_cos2phi = par[3];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALU_offset + ALU_sinphi*sin(phi);
  // return (ALU_offset + ALU_sinphi*sin(phi)) / (1 + AUU_cosphi*cos(phi) + AUU_cos2phi*cos(2*phi));
}

// Function to fit the target-spin asymmetry histogram
double TSA_funcToFit(double* x, double* par) {
  // Retrieve the parameters A
  double AUL_offset = par[0];
  double AUL_sinphi = par[1];
  double AUL_sin2phi = par[2];
  // double AUU_cosphi = par[3];
  // double AUU_cos2phi = par[4];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return AUL_offset + AUL_sinphi*sin(phi)+AUL_sin2phi*sin(2*phi);
  // return (AUL_offset + AUL_sinphi*sin(phi)+AUL_sin2phi*sin(2*phi)) /
  //   (1 + AUU_cosphi*cos(phi) + AUU_cos2phi*cos(2*phi));
}

// Function to fit the double-spin asymmetry histogram
double DSA_funcToFit(double* x, double* par) {
  // Retrieve the parameters A
  double ALL = par[0];
  double ALL_cosphi = par[1];
  // double AUU_cosphi = par[2];
  // double AUU_cos2phi = par[3];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALL+ALL_cosphi*cos(phi);
  // return (ALL+ALL_cosphi*cos(phi)) / (1 + AUU_cosphi*cos(phi) + AUU_cos2phi*cos(2*phi));
}

void plotHistogramAndFit(TH1D* histogram, TF1* fitFunction, int binIndex, int asymmetryIndex, 
  const std::string& prefix) {
  // Define the label for the y-axis
  std::string yAxisLabel, fileNameSuffix;
  switch (asymmetryIndex) {
      case 0: yAxisLabel = "A_{LU}"; fileNameSuffix = "ALU"; break;
      case 1: yAxisLabel = "A_{UL}"; fileNameSuffix = "AUL"; break;
      case 2: yAxisLabel = "A_{LL}"; fileNameSuffix = "ALL"; break;
      default: std::cerr << "Invalid asymmetry index!" << std::endl; return;
  }

  // Create a canvas to draw on
  TCanvas* canvas = new TCanvas("canvas", "", 800, 600);

  // Adjust the canvas margins to ensure axis labels are not cut off
  canvas->SetLeftMargin(0.16); canvas->SetBottomMargin(0.16);

  // Create a TGraphErrors manually from the histogram
  TGraphErrors *graph = new TGraphErrors();
  
  // Add points to the TGraphErrors
  for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
    double x = histogram->GetBinCenter(i);
    double y = histogram->GetBinContent(i);
    double ex = 0;  // we don't want horizontal error bars
    double ey = histogram->GetBinError(i);
    graph->SetPoint(i - 1, x, y);
    graph->SetPointError(i - 1, ex, ey);
  }
  
  // Set the point color to black
  graph->SetMarkerColor(kBlack);
  graph->SetMarkerStyle(kFullCircle);

  // Set the fit function's line color to red
  fitFunction->SetLineColor(kRed);

  // Set the labels of the x and y axis
  graph->GetXaxis()->SetTitle("#phi");
  graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

  // Set the range of the x-axis to be from 0 to 2pi
  graph->GetXaxis()->SetRangeUser(0, 2*TMath::Pi());

  // Draw the graph using the AP option to draw axis and points
  graph->Draw("AP");

  // Set the range of the fit function to match the range of the x-axis
  fitFunction->SetRange(0, 2*TMath::Pi());
  // Draw the fit function on top of the graph
  fitFunction->Draw("same");

  // Center the labels and increase the font size
  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->CenterTitle();
  graph->GetXaxis()->SetTitleSize(0.05);
  graph->GetYaxis()->SetTitleSize(0.05);

  // Create a new TPaveStats object which will serve as our custom statistics box.
  // Adjusted the box position and size to ensure it doesn't overlap with the axes labels
  TPaveStats *statBox = new TPaveStats(0.16, 0.73, 0.45, 0.9, "brNDC");
  // changed coordinates for top left position
  statBox->SetFillColor(0);
  statBox->SetTextSize(0.0225);
  statBox->SetTextAlign(12);
  statBox->SetTextColor(1);
  statBox->SetShadowColor(0); // remove shadow
  
  // Iterate over each parameter in the fit function.
  for (int i = 0; i < fitFunction->GetNpar(); ++i) {
    TText *text=statBox->AddText(Form("Param %d: %.4f +/- %.4f",i,fitFunction->GetParameter(i), 
      fitFunction->GetParError(i))); 
    text->SetTextColor(1);
  }

  TText *text = statBox->AddText(Form("#chi^{2}/Ndf: %.4f", fitFunction->GetChisquare() / 
    fitFunction->GetNDF()));
  text->SetTextColor(1);
  statBox->Draw();

    // Create the filename for the PNG
  std::string filename = "output/" + prefix + "_" + std::to_string(binIndex) + "_" + 
    fileNameSuffix + ".png";
  
  // Create a title string for the graph by removing the "output/" and ".png" portions 
  // of the filename
  std::string title = filename.substr(7, filename.size()-7-4);  
  // start from the 7th index (after "output/") and take (filename.size()-7-4) characters

  // Set the title to the title string
  graph->SetTitle(title.c_str());


  // Save the canvas as a PNG
  canvas->SaveAs(filename.c_str());

  // Clean up
  delete canvas;
  delete graph;
}


void performChi2Fits(const char *filename, const char* output_file, const char* kinematic_file,
  const std::string& prefix, int asymmetry_index) {

  // Initialize string streams to store the results for each bin
  std::ostringstream chi2FitsAStream, chi2FitsBStream, chi2FitsCStream;
  // std::ostringstream chi2FitsDStream, chi2FitsEStream;

  // Initialize string streams to store the mean variables for each bin
  std::ostringstream meanVariablesStream;
  meanVariablesStream << "\\begin{table}[h]" << std::endl;
  meanVariablesStream << "\\centering" << std::endl;
  meanVariablesStream << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|} \\hline" << std::endl;
  meanVariablesStream << "Bin & $<Q^2>$ & $<W>$ ";
  meanVariablesStream << "& $<x_B>$ & $<y>$ & $<z>$ & $<\\zeta>$ & $<P_T>$ ";
  meanVariablesStream << "& $<x_F>$ & $<t>$ & ";
  meanVariablesStream << "$<t_{\\text{min}}>$\\\\ \\hline" << std::endl; 


  // Create a new TF1 object called fitFunction representing the function to fit
  // and create string stream prefix depending on current asymmetry we're fitting
  TF1* fitFunction;
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      fitFunction = new TF1("fitFunction", BSA_funcToFit, 0, 2 * TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALUoffset = {";
      chi2FitsBStream << prefix << "chi2FitsALUsinphi = {";
      // chi2FitsCStream << prefix << "chi2FitsALUAUUcosphi = {";
      // chi2FitsDStream << prefix << "chi2FitsALUAUUcos2phi = {";
      break;
    case 1: // target-spin asymmetry
      fitFunction = new TF1("fitFunction", TSA_funcToFit, 0, 2 * TMath::Pi(), 3);
      chi2FitsAStream << prefix << "chi2FitsAULoffset = {";
      chi2FitsBStream << prefix << "chi2FitsAULsinphi = {";
      chi2FitsCStream << prefix << "chi2FitsAULsin2phi = {";
      // chi2FitsDStream << prefix << "chi2FitsAULAUUcosphi = {";
      // chi2FitsEStream << prefix << "chi2FitsAULAUUcos2phi = {";
      break;
    case 2: // double-spin asymmetry
      fitFunction = new TF1("fitFunction", DSA_funcToFit, 0, 2 * TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALL = {";
      chi2FitsBStream << prefix << "chi2FitsALLcosphi = {";
      // chi2FitsCStream << prefix << "chi2FitsALLAUUcosphi = {";
      // chi2FitsDStream << prefix << "chi2FitsALLAUUcos2phi = {";
      break;
    default:
      cout << "Invalid asymmetry_index! Using default function form of BSA." << endl;
      fitFunction = new TF1("fitFunction", BSA_funcToFit, 0, 2 * TMath::Pi(), 2);
  }

  // Determine the number of bins
  size_t numBins = allBins[currentFits].size() - 1;

  // Loop over each bin
  for (size_t i = 0; i < numBins; ++i) {
    cout << "Beginning chi2 fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    char histName[32];
    snprintf(histName, sizeof(histName), "hist_%zu", i);

    // Create a histogram for the current bin
    TH1D* hist = createHistogramForBin(gData, histName, i, prefix, asymmetry_index);
    // Fit the histogram using the fitFunction and get the fit result
    hist->Fit(fitFunction, "QS");
    plotHistogramAndFit(hist, fitFunction, i, asymmetry_index, prefix);

    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double numEvents = 0;
    // Variables to calculate the mean depolarization factor
    double sumDepA = 0; double sumDepB = 0; double sumDepC = 0; double sumDepV = 0; double sumDepW = 0;

    // Variables to calculate the mean kinematics in each bin
    double sumQ2 = 0; double sumW = 0; double sumx = 0; double sumy = 0;
    double sumz = 0; double sumzeta = 0; double sumpT = 0; double sumxF = 0;
    double sumt = 0; double sumtmin = 0;

    // Loop over all events and calculate the sums and event counts
    int counter = 0;
    for (const eventData& event : gData) {
      // if (counter > 100000) {
      //   break;
      // }
      double currentVariable = getEventProperty(event, currentFits);
      if (applyKinematicCuts(event, currentFits, 0) && currentVariable>=allBins[currentFits][i] && 
        currentVariable < allBins[currentFits][i + 1]) {
          sumVariable += currentVariable;

          // sum the depolarization values
          sumDepA += event.data.at("DepA");
          sumDepB += event.data.at("DepB");
          sumDepC += event.data.at("DepC");
          sumDepV += event.data.at("DepV");
          sumDepW += event.data.at("DepW");

          // sum the kinematic variable values
          sumQ2 += event.data.at("Q2");
          sumW += event.data.at("W");
          sumx += event.data.at("x");
          sumy += event.data.at("y");
          sumz += event.data.at("z");
          sumzeta += event.data.at("zeta");
          sumpT += event.data.at("pT");
          sumxF += event.data.at("xF");
          sumt += event.data.at("t");
          sumtmin += event.data.at("tmin");

          numEvents += 1;
          counter++;
      }
    }
    cout << "Found " << numEvents << " events in this bin." << endl;

    // Calculate the mean values for the variable and depolarization factors
    double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
    double meanDepA = numEvents > 0 ? sumDepA / numEvents : 0.0;
    double meanDepB = numEvents > 0 ? sumDepB / numEvents : 0.0;
    double meanDepC = numEvents > 0 ? sumDepC / numEvents : 0.0;
    double meanDepV = numEvents > 0 ? sumDepV / numEvents : 0.0;
    double meanDepW = numEvents > 0 ? sumDepW / numEvents : 0.0;

    // Calculate the mean values for the kinematic variables
    double meanQ2 = numEvents > 0 ? sumQ2 / numEvents : 0.0;
    double meanW = numEvents > 0 ? sumW / numEvents : 0.0;
    double meanx = numEvents > 0 ? sumx / numEvents : 0.0;
    double meany = numEvents > 0 ? sumy / numEvents : 0.0;
    double meanz = numEvents > 0 ? sumz / numEvents : 0.0;
    double meanzeta = numEvents > 0 ? sumzeta / numEvents : 0.0;
    double meanpT = numEvents > 0 ? sumpT / numEvents : 0.0;
    double meanxF = numEvents > 0 ? sumxF / numEvents : 0.0;
    double meant = numEvents > 0 ? sumt / numEvents : 0.0;
    double meantmin = numEvents > 0 ? sumtmin / numEvents : 0.0;

    switch (asymmetry_index) {
      case 0: {// beam-spin asymmetry
        // Get the fitted parameters and their errors
        double ALU_offset = fitFunction->GetParameter(0);
        double ALU_offset_error = fitFunction->GetParError(0);
        double ALU_sinphi = fitFunction->GetParameter(1); 
        double ALU_sinphi_error = fitFunction->GetParError(1);
        // double AUU_cosphi = fitFunction->GetParameter(2); 
        // double AUU_cosphi_error = fitFunction->GetParError(2);
        // double AUU_cos2phi = fitFunction->GetParameter(3); 
        // double AUU_cos2phi_error = fitFunction->GetParError(3);
        ALU_sinphi = (meanDepA/meanDepW)*ALU_sinphi;
        ALU_sinphi_error = (meanDepA/meanDepW)*ALU_sinphi_error;
        // AUU_cosphi = (meanDepA/meanDepV)*AUU_cosphi;
        // AUU_cosphi_error = (meanDepA/meanDepV)*AUU_cosphi_error;
        // AUU_cos2phi = (meanDepB/meanDepV)*AUU_cos2phi;
        // AUU_cos2phi_error = (meanDepB/meanDepV)*AUU_cos2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALU_offset << ", " << ALU_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALU_sinphi << ", " << ALU_sinphi_error <<"}";
        // chi2FitsCStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", " << AUU_cosphi_error <<"}";
        // chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", " << AUU_cos2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; 
            // chi2FitsCStream << ", "; chi2FitsDStream << ", ";
        }
        break;
      }
      case 1: {// target-spin asymmetry
        // Get the fitted parameters and their errors
        double AUL_offset = fitFunction->GetParameter(0);
        double AUL_offset_error = fitFunction->GetParError(0);
        double AUL_sinphi = fitFunction->GetParameter(1);
        double AUL_sinphi_error = fitFunction->GetParError(1);
        double AUL_sin2phi = fitFunction->GetParameter(2);
        double AUL_sin2phi_error = fitFunction->GetParError(2);
        // double AUU_cosphi = fitFunction->GetParameter(3); 
        // double AUU_cosphi_error = fitFunction->GetParError(3);
        // double AUU_cos2phi = fitFunction->GetParameter(4); 
        // double AUU_cos2phi_error = fitFunction->GetParError(4);
        AUL_sinphi = (meanDepA/meanDepV)*AUL_sinphi;
        AUL_sinphi_error = (meanDepA/meanDepV)*AUL_sinphi_error;
        AUL_sin2phi = (meanDepA/meanDepB)*AUL_sin2phi;
        AUL_sin2phi_error = (meanDepA/meanDepB)*AUL_sin2phi_error;
        // AUU_cosphi = (meanDepA/meanDepV)*AUU_cosphi;
        // AUU_cosphi_error = (meanDepA/meanDepV)*AUU_cosphi_error;
        // AUU_cos2phi = (meanDepB/meanDepV)*AUU_cos2phi;
        // AUU_cos2phi_error = (meanDepB/meanDepV)*AUU_cos2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< AUL_offset << ", " << AUL_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< AUL_sinphi << ", " << AUL_sinphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUL_sin2phi << ", " << AUL_sin2phi_error <<"}";
        // chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", " << AUU_cosphi_error <<"}";
        // chi2FitsEStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", " << AUU_cos2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; chi2FitsCStream << ", ";
            // chi2FitsDStream << ", "; chi2FitsEStream << ", ";
        }
        break;
      }
      case 2: {// double-spin asymmetry
        // Get the fitted parameters and their errors
        double ALL = fitFunction->GetParameter(0);
        double ALL_error = fitFunction->GetParError(0);
        double ALL_cosphi = fitFunction->GetParameter(1);
        double ALL_cosphi_error = fitFunction->GetParError(1);
        // double AUU_cosphi = fitFunction->GetParameter(2); 
        // double AUU_cosphi_error = fitFunction->GetParError(2);
        // double AUU_cos2phi = fitFunction->GetParameter(3); 
        // double AUU_cos2phi_error = fitFunction->GetParError(3);
        ALL = (meanDepA/meanDepC)*ALL;
        ALL_error = (meanDepA/meanDepC)*ALL_error;
        ALL_cosphi = (meanDepA/meanDepW)*ALL_cosphi;
        ALL_cosphi_error = (meanDepA/meanDepW)*ALL_cosphi_error;
        // AUU_cosphi = (meanDepA/meanDepV)*AUU_cosphi;
        // AUU_cosphi_error = (meanDepA/meanDepV)*AUU_cosphi_error;
        // AUU_cos2phi = (meanDepB/meanDepV)*AUU_cos2phi;
        // AUU_cos2phi_error = (meanDepB/meanDepV)*AUU_cos2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALL << ", " << ALL_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALL_cosphi << ", " << ALL_cosphi_error <<"}";
        // chi2FitsCStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", " << AUU_cosphi_error <<"}";
        // chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", " << AUU_cos2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", ";
            // chi2FitsCStream << ", "; chi2FitsDStream << ", ";
        }
        break;
      }
    }

    delete hist;

    // outputs of mean kinematic variables
    meanVariablesStream << std::fixed << std::setprecision(3); // Set precision to 3 digits 
    meanVariablesStream << (i+1) << "~&~" << meanQ2 << "~&~" << meanW << "~&~" << meanx << "~&~";
    meanVariablesStream << meany << "~&~" << meanz << "~&~" << meanzeta << "~&~";
    meanVariablesStream << meanpT << "~&~" << meanxF << "~&~" << meant << "~&~" << meantmin; 
    meanVariablesStream << std::string(" \\\\ \\hline ");

  }

  chi2FitsAStream << "};";  chi2FitsBStream << "};";  chi2FitsCStream << "};"; 
  // chi2FitsDStream << "};";  chi2FitsEStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << chi2FitsAStream.str() << std::endl;
  outputFile << chi2FitsBStream.str() << std::endl;
  if (asymmetry_index==1) { outputFile << chi2FitsCStream.str() << std::endl; }

  // outputFile << chi2FitsCStream.str() << std::endl;
  // outputFile << chi2FitsDStream.str() << std::endl;
  // if (asymmetry_index==1) { outputFile << chi2FitsEStream.str() << std::endl; }
  outputFile.close();

  meanVariablesStream << "\\end{tabular}\n";
  meanVariablesStream << "\\caption{The mean kinematic variables in each of the bins ";
  meanVariablesStream << "for the extracted $" << prefix << "$ asymmetries.}\n";
  meanVariablesStream << "\\label{table:kinematics_" << prefix << "}\n";
  meanVariablesStream << "\\end{table}. Values given in GeV or GeV$^2$ where appropriate.\n";
  meanVariablesStream << endl << endl << endl;
  if (asymmetry_index == 0) {
    std::ofstream kinematicFile(kinematic_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicFile << meanVariablesStream.str() << std::endl; 
    kinematicFile.close();
  }

}

void BSA_rgc_fits_text_file(const char* data_file, const char* mc_file, const char* output_file, 
  const char* kinematic_file) {
  // Start the timer
  auto start_time = std::chrono::high_resolution_clock::now();

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

  // Read data from the input file and store it in the global variable gData
  gData = readData(data_file, variable_names);
  // Read mc from the input file and store it in the global variable gData
  gMC = readData(mc_file, variable_names);

  cout << endl << endl;
  for (size_t i = 0; i < allBins.size(); ++i) {
    cout << "-- Beginning kinematic fits." << endl;
    for (int asymmetry = 0; asymmetry < 1; ++asymmetry){
      switch (asymmetry) {
        case 0: cout << "    chi2 BSA." << endl; break;
        case 1: cout << "    chi2 TSA." << endl; break;
        case 2: cout << "    chi2 DSA." << endl; break;
      }
      performChi2Fits(data_file, output_file, kinematic_file, binNames[i], asymmetry);
    }
    cout << endl << "     Completed " << binNames[i] << " chi2 fits." << endl;
    // performMLMFits(data_file, output_file, kinematic_file, binNames[i]);
    // cout << endl << "     Completed " << binNames[i] << " MLM fits." << endl;
    cout << endl << endl;
    currentFits++;
  }
  cout << endl << endl;

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
}