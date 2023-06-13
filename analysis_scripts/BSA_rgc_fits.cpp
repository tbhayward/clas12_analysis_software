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
size_t currentBin = 0;
int n = 1;

std::map<std::string, std::vector<float>> bins_map;
std::vector<std::vector<float>> allBins;
std::vector<std::string> binNames;
std::vector<std::string> propertyNames;
std::vector<std::string> variable_names;
float total_charge_carbon;
float cpp, cpm, cmp, cmm;

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

struct RunInfo {
  int runnum;
  float total_charge;
  float positive_charge;
  float negative_charge;
  float target_polarization;
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

    // Read the total charge from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.total_charge = std::stof(info);

    // Read the positive charge from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.positive_charge = std::stof(info);

    // Read the negative charge from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.negative_charge = std::stof(info);

    // Read the target polarization from the stringstream and convert it to a float
    std::getline(ss, info, ',');
    run_info.target_polarization = std::stof(info);

    // Add the struct to the run_info_list vector
    run_info_list.push_back(run_info);
  }
}

// function to get the polarization value
float getPol(int runnum) {
  float pol = 0.86; 
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
    else if (runnum >= 16000) { pol = 0.81; } // RGC
  return pol;
}

struct eventData {
  std::unordered_map<std::string, float> data;
};

std::vector<eventData> gData;

eventData parseLine(const std::string& line, const std::vector<std::string>& variable_names) {
  // Create a stringstream from the input line
  std::istringstream iss(line);

  // Initialize an eventData object to store the parsed data
  eventData data;

  // Declare value for storing the extracted float and an index for variable names
  float value;
  std::string value_str;
  size_t var_name_index = 0;

  // Iterate through the provided variable names
  for (const auto& var_name : variable_names) {
    // Attempt to extract a float value from the stringstream
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
  for (const auto& run_info : run_info_list) {
    if (run_info.runnum == runnum) {
      data.data["target_pol"] = run_info.target_polarization;
      break;
    }
  }

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
bool applyKinematicCuts(const eventData& data, int currentFits) {

    return data.data.at("Q2")>1 && data.data.at("W")>2 && data.data.at("Mx")>1.5 &&
      data.data.at("y")<0.8 && data.data.at("xF")<0 && data.data.at("target_pol") != 0;
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

    // Initialize variables for counting events (N), positive helicity sum (sum_P), 
    // and negative helicity sum (sum_N)
    double N = 0;
    double sum_PP = 0; // positive beam -- positive target
    double sum_PM = 0; // positive beam -- negative target
    double sum_MP = 0; // negative beam -- positive target
    double sum_MM = 0; // negative beam -- negative target

    float Df = 0.18; // dilution factor, placeholder from MC studies from proposal

    // Iterate through the global event data (gData)
    for (const eventData &event : gData) {
        // Get the value of the current variable of interest for the event
        double currentVariable = getEventProperty(event, currentFits);

        // Apply kinematic cuts and check if the current variable is within the specified bin range
        if (applyKinematicCuts(event, currentFits) && 
          currentVariable >= allBins[currentFits][currentBin] && 
          currentVariable < allBins[currentFits][currentBin + 1]) {

          // Increment the event count
          N += 1;

          // Extract Delta_phi and polarization (pol) from the event data
          double phi = event.data.at("phi");
          double Pb = event.data.at("pol");
          double Pt = std::abs(event.data.at("target_pol"));

          // Check if the helicities is positive or negative and update the corresponding sum
          if (event.data.at("helicity") > 0 && event.data.at("target_pol") > 0) {
            sum_PP += log(1 + 0*Pb*(ALU_sinphi*sin(phi)) // BSA
              + Df*Pt*(AUL_sinphi*sin(phi) + AUL_sin2phi*sin(2*phi)) // TSA
              + 0*Df*Pb*Pt*(ALL + ALL_cosphi*cos(phi)) ); // DSA
          } else if (event.data.at("helicity") > 0 && event.data.at("target_pol") < 0 ) {
            sum_PM += log(1 + 0*Pb*(ALU_sinphi*sin(phi)) // BSA
              - Df*Pt*(AUL_sinphi*sin(phi) + AUL_sin2phi*sin(2*phi)) // TSA
              - 0*Df*Pb*Pt*(ALL + ALL_cosphi*cos(phi)) ); // DSA
          } else if (event.data.at("helicity") < 0 && event.data.at("target_pol") > 0 ) {
            sum_MP += log(1 - 0*Pb*(ALU_sinphi*sin(phi)) // BSA
              + Df*Pt*(AUL_sinphi*sin(phi) + AUL_sin2phi*sin(2*phi)) // TSA
              - 0*Df*Pb*Pt*(ALL + ALL_cosphi*cos(phi)) ); // DSA
          } else if (event.data.at("helicity") < 0 && event.data.at("target_pol") < 0 ) {
            sum_MM += log(1 - 0*Pb*(ALU_sinphi*sin(phi)) // BSA
              - Df*Pt*(AUL_sinphi*sin(phi) + AUL_sin2phi*sin(2*phi)) // TSA
              + 0*Df*Pb*Pt*(ALL + ALL_cosphi*cos(phi)) ); // DSA
          }
        }
    }

    // determine min pos or neg beam helicity accumulated charge to scale down higher one
    float minBeamCharge = std::min({(cpp+cpm),(cmp+cmm)}); 
    // determine min pos or neg target helicity accumulated charge to scale down higher one
    float minTargetCharge = std::min({(cpp+cmp),(cpm+cmm)}); 
    
    // Calculate the negative log-likelihood value and store it in the output variable f
    f = N * log(N) - 
      minBeamCharge*minTargetCharge/((cpp+cpm)*(cpp+cmp))*sum_PP -
      minBeamCharge*minTargetCharge/((cpp+cpm)*(cmp+cmm))*sum_PM - 
      minBeamCharge*minTargetCharge/((cmp+cmm)*(cpp+cmp))*sum_MP - 
      minBeamCharge*minTargetCharge/((cmp+cmm)*(cmp+cmm))*sum_MM;
}

void performMLMFits(const char *filename, const char* output_file, const std::string& prefix) {
  // Read the event data from the input file and store it in the global variable gData
  gData = readData(filename, variable_names);

  // Determine the number of bins
  size_t numBins = allBins[currentFits].size() - 1;

  // Initialize TMinuit with 3 parameters 
  double arglist[10]; arglist[0] = 1;
  int ierflg = 0;
  TMinuit minuit(5);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(negLogLikelihood);

  // Declare string streams for storing the MLM fit results
  std::ostringstream mlmFitsAStream;
  std::ostringstream mlmFitsBStream; std::ostringstream mlmFitsCStream;
  std::ostringstream mlmFitsDStream; std::ostringstream mlmFitsEStream;

  // Initialize the string streams with the output variable names
  mlmFitsAStream << prefix << "MLMFits_ALU_sinphi = {";
  mlmFitsBStream << prefix << "MLMFits_AUL_sinphi = {";
  mlmFitsCStream << prefix << "MLMFits_AUL_sin2phi = {";
  mlmFitsDStream << prefix << "MLMFits_ALL = {";
  mlmFitsEStream << prefix << "MLMFits_ALL_cosphi = {";

  // Iterate through each bin
  for (size_t i = 0; i < numBins; ++i) {
    cout << endl << "Beginning MLM fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    currentBin = i;

    // Define the parameters with initial values and limits
    minuit.DefineParameter(0, "ALU_sinphi", -0.02, 0.01, -1, 1);
    minuit.DefineParameter(1, "AUL_sinphi", -0.02, 0.01, -1, 1);
    minuit.DefineParameter(2, "AUL_sin2phi", -0.01, 0.01, -1, 1);
    minuit.DefineParameter(3, "ALL", 0.3, 0.01, -1, 1);
    minuit.DefineParameter(4, "ALL_cosphi", 0.01, 0.01, -1, 1);

    // Minimize the negative log-likelihood function
    minuit.Migrad();

    // Calculate the mean values of the current variable and the back-to-back factor (b2b_factor)
    double sumVariable = 0;
    double numEvents = 0;
    for (const eventData &event : gData) {
      double currentVariable = getEventProperty(event, currentFits);
        if (applyKinematicCuts(event, currentFits) && currentVariable >= 
          allBins[currentFits][i] && currentVariable < allBins[currentFits][i + 1]) {
            sumVariable += currentVariable;
            numEvents += 1;
        }
    }
    double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;


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

    // output to text file
    mlmFitsAStream << "{" << meanVariable << ", " << ALU_sinphi << ", " << ALU_sinphi_error << "}";
    mlmFitsBStream << "{" << meanVariable << ", " << AUL_sinphi << ", " << AUL_sinphi_error << "}";
    mlmFitsCStream << "{" << meanVariable << ", " << AUL_sin2phi << ", "<<AUL_sin2phi_error << "}";
    mlmFitsDStream << "{" << meanVariable << ", " << ALL << ", " << ALL_error << "}";
    mlmFitsEStream << "{" << meanVariable << ", " << ALL_cosphi << ", "<<ALL_cosphi_error << "}";

    if (i < numBins - 1) {
        mlmFitsAStream << ", "; mlmFitsBStream << ", "; mlmFitsCStream << ", ";
        mlmFitsDStream << ", "; mlmFitsEStream << ", "; 
    }
  }

  mlmFitsAStream << "};"; mlmFitsBStream << "};"; mlmFitsCStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << mlmFitsAStream.str() << std::endl;
  outputFile << mlmFitsBStream.str() << std::endl;
  outputFile << mlmFitsCStream.str() << std::endl;
  outputFile << mlmFitsDStream.str() << std::endl;
  outputFile << mlmFitsEStream.str() << std::endl;

  outputFile.close();
}

float asymmetry_value_calculation(float Npp, float Npm, float Nmp, float Nmm, float meanPol, 
  float Ptp, float Ptm, int asymmetry_index) {
  float Df = 0.18; // dilution factor, placeholder from MC studies from proposal
  // return the asymmetry value 
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      return (1 / meanPol) * (Ptm*(Npp-Nmp)+Ptp*(Npm-Nmm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    case 1: // target-spin asymmetry
      cout << endl << Npp << " " << Npm << " " << Nmp << " " << Nmm << " " << Ptp << " " << Ptm << endl;
      cout << (1 / Df) * ((Npp+Nmp)-(Npm+Nmm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm)) << endl << endl;
      return (1 / Df) * ((Npp+Nmp)-(Npm+Nmm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    case 2: // double-spin asymmetry
      return (1 / (Df*meanPol)) * ((Npp-Nmp)+(Nmm-Npm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    default:
      cout << "Invalid asymmetry_index!" << endl;
      return 0;
  }
}

float asymmetry_error_calculation(float Npp, float Npm, float Nmp, float Nmm, float meanPol, 
  float Ptp, float Ptm, int asymmetry_index) {
  float Df = 0.18; // dilution factor, placeholder from MC studies from proposal
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
  int binIndex, int asymmetry_index) {

  // Determine the variable range for the specified bin
  double varMin = allBins[currentFits][binIndex];
  double varMax = allBins[currentFits][binIndex + 1];

  // Create positive and negative helicity histograms
  TH1D* histPosPos = new TH1D(Form("%s_pospos", histName), "", 12, 0, 2 * TMath::Pi());
  TH1D* histPosNeg = new TH1D(Form("%s_posneg", histName), "", 12, 0, 2 * TMath::Pi());
  TH1D* histNegPos = new TH1D(Form("%s_negpos", histName), "", 12, 0, 2 * TMath::Pi());
  TH1D* histNegNeg = new TH1D(Form("%s_negneg", histName), "", 12, 0, 2 * TMath::Pi());

  // Variables to calculate the mean polarization
  double sumPol = 0; // sum of the beam polarization
  double sumTargetPosPol = 0; // sum of the target positive polarization
  double sumTargetNegPol = 0; // sum of the target negative polarization
  int numEvents = 0;
  int numEventsPosTarget = 0;
  int numEventsNegTarget = 0;

  // Fill the positive and negative helicity histograms
  for (const eventData& event : data) {
    double currentVariable = getEventProperty(event, currentFits);
    if (applyKinematicCuts(event, currentFits) && currentVariable >= varMin && 
      currentVariable < varMax) {
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
  // scale the histograms by the accumulated faraday cup charge
  // histPosPos->Scale(1.0 / cpp);
  // histPosNeg->Scale(1.0 / cpm);
  // histNegPos->Scale(1.0 / cmp);
  // histNegNeg->Scale(1.0 / cmm);

  // Calculate the mean polarization
  float meanPol = sumPol / numEvents; // mean beam polarization for data 
  float Ptp = sumTargetPosPol / numEventsPosTarget;// mean positive target polarization for data
  float Ptm = - sumTargetNegPol / numEventsNegTarget;// mean negative target polarization for data
  // the negative sign here is correct; RGC lists the polarizations with signs to tell which is 
  // which but the polarization really should just be "percent of polarized nucleii"

  // Create the asymmetry histogram
  int numBins = histPosPos->GetNbinsX();
  TH1D* histAsymmetry = new TH1D(Form("%s_asymmetry", histName), "", 
    numBins, 0, 2 * TMath::Pi());

  // Calculate the asymmetry and its error for each bin, and fill the asymmetry histogram
  for (int iBin = 1; iBin <= numBins; ++iBin) {
    float Npp = histPosPos->GetBinContent(iBin)/cpp;
    float Npm = histPosNeg->GetBinContent(iBin)/cpm;
    float Nmp = histNegPos->GetBinContent(iBin)/cmp;
    float Nmm = histNegNeg->GetBinContent(iBin)/cmm;

    // Calculate the asymmetry and error for the current bin
    float asymmetry = asymmetry_value_calculation(Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, 
      asymmetry_index);
    float error = asymmetry_error_calculation(Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, 
      asymmetry_index);

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
  // Retrieve the parameters A
  double ALU_sinphi = par[0];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALU_sinphi*sin(phi);
}

// Function to fit the target-spin asymmetry histogram
double TSA_funcToFit(double* x, double* par) {
  // Retrieve the parameters A
  double AUL_sinphi = par[0];
  double AUL_sin2phi = par[1];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return AUL_sinphi*sin(phi)+AUL_sin2phi*sin(2*phi);
}

// Function to fit the double-spin asymmetry histogram
double DSA_funcToFit(double* x, double* par) {
  // Retrieve the parameters A
  double ALL = par[0];
  double ALL_cosphi = par[1];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALL+ALL_cosphi*cos(phi);
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
  canvas->SetLeftMargin(0.12); canvas->SetBottomMargin(0.12);

  // Set the histogram's line and point color to black
  histogram->SetLineColor(kBlack);
  histogram->SetMarkerColor(kBlack);
  histogram->SetMarkerStyle(kFullCircle);  

  // Set the fit function's line color to red
  fitFunction->SetLineColor(kRed);

  // Draw the histogram using the E option to draw just the points with error bars
  histogram->Draw("E1");

  // Draw the fit function on top of the histogram
  fitFunction->Draw("same");

  // Set the labels of the x and y axis
  histogram->GetXaxis()->SetTitle("#phi");
  histogram->GetYaxis()->SetTitle(yAxisLabel.c_str());

  // Center the labels and increase the font size
  histogram->GetXaxis()->CenterTitle();
  histogram->GetYaxis()->CenterTitle();
  histogram->GetXaxis()->SetTitleSize(0.05);
  histogram->GetYaxis()->SetTitleSize(0.05);

  // Customize the stat box
  histogram->SetStats(0);  // Turn off automatic stats
  TPaveStats *statBox = new TPaveStats(0.6, 0.8, 0.9, 0.9, "brNDC");
  statBox->SetFillColor(0);
  statBox->SetTextSize(0.035);
  statBox->SetTextAlign(12);
  statBox->SetTextColor(1);
  TText *text = statBox->AddText(Form("Entries = %.0f", histogram->GetEntries()));
  text->SetTextColor(1);
  text = statBox->AddText(Form("Chi^2/Ndf = %.4f", fitFunction->GetChisquare() / 
    fitFunction->GetNDF()));
  text->SetTextColor(1);
  statBox->Draw();

  // Create the filename for the PNG
  std::string filename = "output/" + prefix + "_" + std::to_string(binIndex) + "_" + 
    fileNameSuffix + ".png";

  // Save the canvas as a PNG
  canvas->SaveAs(filename.c_str());

  // Clean up
  delete canvas;
}



void performChi2Fits(const char *filename, const char* output_file, const std::string& prefix, 
  int asymmetry_index) {
  // Read data from the input file and store it in the global variable gData
  gData = readData(filename, variable_names);

  // Initialize string streams to store the results for each bin
  std::ostringstream chi2FitsAStream, chi2FitsBStream;

  // Create a new TF1 object called fitFunction representing the function to fit
  // and create string stream prefix depending on current asymmetry we're fitting
  TF1* fitFunction;
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      fitFunction = new TF1("fitFunction", BSA_funcToFit, 0, 2 * TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALUsinphi = {";
      break;
    case 1: // target-spin asymmetry
      fitFunction = new TF1("fitFunction", TSA_funcToFit, 0, 2 * TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsAULsinphi = {";
      chi2FitsBStream << prefix << "chi2FitsAULsin2phi = {";
      break;
    case 2: // double-spin asymmetry
      fitFunction = new TF1("fitFunction", DSA_funcToFit, 0, 2 * TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALL = {";
      chi2FitsBStream << prefix << "chi2FitsALLcosphi = {";
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
    TH1D* hist = createHistogramForBin(gData, histName, i, asymmetry_index);
    // Fit the histogram using the fitFunction and get the fit result
    hist->Fit(fitFunction, "QS");
    plotHistogramAndFit(hist, fitFunction, i, asymmetry_index, prefix);

    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double numEvents = 0;
    // Loop over all events and calculate the sums and event counts
    for (const eventData& event : gData) {
      double currentVariable = getEventProperty(event, currentFits);
      if (applyKinematicCuts(event, currentFits) && currentVariable >= allBins[currentFits][i] && 
        currentVariable < allBins[currentFits][i + 1]) {
          sumVariable += currentVariable;
          numEvents += 1;
      }
    }
    cout << "Found " << numEvents << " events in this bin." << endl;

    // Calculate the mean values for the variable and b2b factor
    double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;

    switch (asymmetry_index) {
      case 0: {// beam-spin asymmetry
        // Get the fitted parameters and their errors
        double ALU_sinphi = fitFunction->GetParameter(0);
        double ALU_sinphi_error = fitFunction->GetParError(0);
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALU_sinphi << ", " << ALU_sinphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", ";
        }
        break;
      }
      case 1: {// target-spin asymmetry
        // Get the fitted parameters and their errors
        double AUL_sinphi = fitFunction->GetParameter(0);
        double AUL_sinphi_error = fitFunction->GetParError(0);
        double AUL_sin2phi = fitFunction->GetParameter(1);
        double AUL_sin2phi_error = fitFunction->GetParError(1);
        chi2FitsAStream<<"{"<<meanVariable<<", "<< AUL_sinphi << ", " << AUL_sinphi_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< AUL_sin2phi << ", " << AUL_sin2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", ";
            chi2FitsBStream << ", ";
        }
        break;
      }
      case 2: {// double-spin asymmetry
        // Get the fitted parameters and their errors
        double ALL = fitFunction->GetParameter(0);
        double ALL_error = fitFunction->GetParError(0);
        double ALL_cosphi = fitFunction->GetParameter(1);
        double ALL_cosphi_error = fitFunction->GetParError(1);
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALL << ", " << ALL_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALL_cosphi << ", " << ALL_cosphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", ";
            chi2FitsBStream << ", ";
        }
        break;
      }
    }

    delete hist;
  }

  chi2FitsAStream << "};"; 
  chi2FitsBStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << chi2FitsAStream.str() << std::endl;
  if (asymmetry_index>0) { outputFile << chi2FitsBStream.str() << std::endl; }

  outputFile.close();
}

void BSA_rgc_fits(const char* data_file, const char* output_file) {

  // Clear the contents of the output_file
  std::ofstream ofs(output_file, std::ios::trunc);
  ofs.close();

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

  cout << endl << endl;
  for (size_t i = 0; i < allBins.size(); ++i) {
    cout << "-- Beginning kinematic fits." << endl;
    for (int asymmetry = 1; asymmetry < 2; ++asymmetry){
      switch (asymmetry) {
        case 0: cout << "    chi2 BSA." << endl; break;
        case 1: cout << "    chi2 TSA." << endl; break;
        case 2: cout << "    chi2 DSA." << endl; break;
      }
      performChi2Fits(data_file, output_file, binNames[i], asymmetry);
    }
    cout << endl << "     Completed " << binNames[i] << " chi2 fits." << endl;
    // performMLMFits(data_file, output_file, binNames[i]);
    // cout << endl << "     Completed " << binNames[i] << " MLM fits." << endl;
    cout << endl << endl;
    currentFits++;
  }
  cout << endl << endl;
}