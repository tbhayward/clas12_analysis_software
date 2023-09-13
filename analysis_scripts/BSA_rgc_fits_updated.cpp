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
#include <iomanip> 
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TF1.h> 
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TText.h>
#include <TApplication.h>
using namespace std;


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

// Apply kinematic cuts to the data
bool applyKinematicCuts(TTree* data, int entry, int currentFits, bool isMC) {

  bool goodEvent = 0;
  std::string property = binNames[currentFits];

  double target_pol, Q2, W, Mx, x, y, pT, xF;
  data->SetBranchAddress("target_pol", &target_pol);
  data->SetBranchAddress("Q2", &Q2);
  data->SetBranchAddress("W", &W);
  data->SetBranchAddress("Mx", &Mx);
  data->SetBranchAddress("x", &x);
  data->SetBranchAddress("y", &y);
  data->SetBranchAddress("pT", &pT);
  data->SetBranchAddress("xF", &xF);

  data->GetEntry(entry);
  //
  //
  //
  // epX
  if (property == "xF") {
    goodEvent = Q2>1 && W>2 && Mx>1.4 && y <0.75;
  }
  if (property == "Q2bin") {
    goodEvent = Q2>1 && W>2 && Mx>1.4 && y<0.75 && x>0.2 && x<0.3 && pT>0.25 && pT<0.35 && 
      xF<0;
  }
  if (property == "PTTFR" || property ==  "xTFR" || property == "zetaTFR" || 
    property == "Q2TFR" || property ==  "x") {
    goodEvent = Q2>1 && W>2 && Mx>1.4 && y<0.75 && xF<0;
  }
  if (property == "PTCFR" || property == "xCFR" || property == "zetaCFR" ||
    property == "Q2TFR") {
    goodEvent = Q2>1 && W>2 && Mx>1.4 && y<0.75 && xF>0;
  } 
  //
  //
  //
  // epi+X
  if (property == "xFpip") { 
    goodEvent = Q2>1 && W>2 && Mx>1.5 && y<0.75;
  }
  if (property == "PTTFRpip" || property ==  "xTFRpip" || property == "zTFRpip" || 
    property == "Q2TFRpip" || property ==  "xpip") {
    goodEvent = Q2>1 && W>2 && Mx>1.5 && y<0.75 && xF<0;
  }
  if (property == "PTCFRpip" || property == "xCFRpip" || property == "zCFRpip" ||
    property == "Q2TFRpip") {
    goodEvent = Q2>1 && W>2 && Mx>1.5 && y<0.75 && xF>0;
  }
  //
  //
  //
  // epi-X
  if (property == "xFpim") { 
    goodEvent = Q2>1 && W>2 && Mx>1.5 && y<0.75;
  }
  if (property == "PTTFRpim" || property ==  "xTFRpim" || property == "zTFRpim" || 
    property == "Q2TFRpim" || property ==  "xpim") {
    goodEvent = Q2>1 && W>2 && Mx>1.5 && y<0.75 && xF<0;
  }
  if (property == "PTCFRpim" || property == "xCFRpim" || property == "zCFRpim" ||
    property == "Q2TFRpim") {
    goodEvent = Q2>1 && W>2 && Mx>1.5 && y<0.75 && xF>0;
  } 
  if (isMC) { return goodEvent; }
  else {return goodEvent && target_pol != 0; } // if data, skip Pt = 0 (carbon)

  return goodEvent;  
}

float dilution_factor(float currentVariable, const std::string& prefix) {

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

TH1D* createHistogramForBin(TTree* data, const char* histName, int binIndex, 
  const std::string& prefix, int asymmetry_index) {

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
  double sumTargetPosPol = 0; // sum of the target positive polarization
  double sumTargetNegPol = 0; // sum of the target negative polarization
  int numEventsPosTarget = 0;
  int numEventsNegTarget = 0;

  int helicity; data->SetBranchAddress("helicity", &helicity); // beam helicity 
  double beam_pol; data->SetBranchAddress("beam_pol", &beam_pol); // beam polarization
  double target_pol; data->SetBranchAddress("target_pol", &target_pol); // target polarization
  double phi; data->SetBranchAddress("phi", &phi); // trento phi
  double xF; data->SetBranchAddress("xF", &xF); // xF
  TTree* data_copied = data;
  double currentVariable; 
  data_copied->SetBranchAddress(propertyNames[currentFits].c_str(), &currentVariable); 

  cout << endl;
  for (int entry = 0; entry < 10; ++entry) {  // Just read the first 10 entries for debugging
    data->GetEntry(entry);
    data_copied->GetEntry(entry);
    cout << "Entry " << entry << " : " << propertyNames[currentFits].c_str();
    cout << " " << currentVariable << " " << xF << " " << helicity << endl;
  }

  cout << "End first loop" << endl;

  // for (int entry = 0; entry < data->GetEntries(); ++entry) {
  for (int entry = 0; entry < 10; ++entry) {
    data->GetEntry(entry);
    data_copied->GetEntry(entry);
    cout << "Entry " << entry << " : " << propertyNames[currentFits].c_str();
    cout << " " << currentVariable << " " << xF << " " << helicity << endl;
    if (applyKinematicCuts(data, entry, currentFits, 0) && currentVariable >= varMin && 
      currentVariable < varMax) {
    //   sumVariable+=currentVariable;

    //   if (helicity > 0 && target_pol > 0) { histPosPos->Fill(phi); } 
    //   else if (helicity > 0 && target_pol < 0) { histPosNeg->Fill(phi); } 
    //   else if (helicity < 0 && target_pol > 0) { histNegPos->Fill(phi); } 
    //   else if (helicity < 0 && target_pol < 0) { histNegNeg->Fill(phi); }

    //   // Accumulate polarization and event count for mean polarization calculation
    //   sumPol += beam_pol;
    //   if (target_pol > 0) {
    //     sumTargetPosPol+=target_pol;
    //     numEventsPosTarget++;
    //   } else if (target_pol < 0) {
    //     sumTargetNegPol+=target_pol;
    //     numEventsNegTarget++;
    //   }
    //   numEvents++;
    }
  }

  // Calculate the mean polarization
  double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
  double meanPol = sumPol / numEvents; // mean beam polarization for data 
  double Ptp = sumTargetPosPol / numEventsPosTarget;// mean positive target polarization for data
  double Ptm = - sumTargetNegPol / numEventsNegTarget;// mean negative target polarization for data
  // the negative sign here is correct; RGC lists the polarizations with signs to tell which is 
  // which but the polarization really should just be "percent of polarized nucleii"

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
    float x = histogram->GetBinCenter(i);
    float y = histogram->GetBinContent(i);
    float ex = 0;  // we don't want horizontal error bars
    float ey = histogram->GetBinError(i);
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
  // canvas->SaveAs(filename.c_str());

  // Clean up
  delete canvas;
  delete graph;
}

void performChi2Fits(TTree* data, const char* output_file, const char* kinematic_file,
  const std::string& prefix, int asymmetry_index) {

  // Initialize string streams to store the results for each bin
  std::ostringstream chi2FitsAStream, chi2FitsBStream, chi2FitsCStream;
  // std::ostringstream chi2FitsDStream, chi2FitsEStream;

  // Initialize string streams to store the mean variables for each bin
  std::ostringstream meanVariablesStream;
  meanVariablesStream << "\\begin{table}[h]" << endl;
  meanVariablesStream << "\\centering" << endl;
  meanVariablesStream << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|} \\hline" << endl;
  meanVariablesStream << "Bin & $<Q^2>$ & $<W>$ ";
  meanVariablesStream << "& $<x_B>$ & $<y>$ & $<z>$ & $<\\zeta>$ & $<P_T>$ ";
  meanVariablesStream << "& $<x_F>$ & $<t>$ & ";
  meanVariablesStream << "$<t_{\\text{min}}>$\\\\ \\hline" << endl; 

  // Create a new TF1 object called fitFunction representing the function to fit
  // and create string stream prefix depending on current asymmetry we're fitting
  TF1* fitFunction;
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      fitFunction = new TF1("fitFunction", BSA_funcToFit, 0, 2*TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALUoffset = {";
      chi2FitsBStream << prefix << "chi2FitsALUsinphi = {";
      // chi2FitsCStream << prefix << "chi2FitsALUAUUcosphi = {";
      // chi2FitsDStream << prefix << "chi2FitsALUAUUcos2phi = {";
      break;
    case 1: // target-spin asymmetry
      fitFunction = new TF1("fitFunction", TSA_funcToFit, 0, 2*TMath::Pi(), 3);
      chi2FitsAStream << prefix << "chi2FitsAULoffset = {";
      chi2FitsBStream << prefix << "chi2FitsAULsinphi = {";
      chi2FitsCStream << prefix << "chi2FitsAULsin2phi = {";
      // chi2FitsDStream << prefix << "chi2FitsAULAUUcosphi = {";
      // chi2FitsEStream << prefix << "chi2FitsAULAUUcos2phi = {";
      break;
    case 2: // double-spin asymmetry
      fitFunction = new TF1("fitFunction", DSA_funcToFit, 0, 2*TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALL = {";
      chi2FitsBStream << prefix << "chi2FitsALLcosphi = {";
      // chi2FitsCStream << prefix << "chi2FitsALLAUUcosphi = {";
      // chi2FitsDStream << prefix << "chi2FitsALLAUUcos2phi = {";
      break;
    default:
      cout << "Invalid asymmetry_index! Using default function form of BSA." << endl;
      fitFunction = new TF1("fitFunction", BSA_funcToFit, 0, 2*TMath::Pi(), 2);
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
    TH1D* hist = createHistogramForBin(data, histName, i, prefix, asymmetry_index);
    // Fit the histogram using the fitFunction and get the fit result
    hist->Fit(fitFunction, "QS");
    plotHistogramAndFit(hist, fitFunction, i, asymmetry_index, prefix);

    // Initialize variables to store the sums and event counts
    float sumVariable = 0;
    float numEvents = 0;
    // Variables to calculate the mean depolarization factor
    float sumDepA = 0; float sumDepB = 0; float sumDepC = 0; float sumDepV = 0; float sumDepW = 0;

    // Variables to calculate the mean kinematics in each bin
    float sumQ2 = 0; float sumW = 0; float sumx = 0; float sumy = 0;
    float sumz = 0; float sumzeta = 0; float sumpT = 0; float sumxF = 0;
    float sumt = 0; float sumtmin = 0;



    delete hist;
  }


}

int main(int argc, char *argv[]) {
  // initialize ROOT application for graphics
  TApplication theApp("App", nullptr, nullptr);
  // Set ROOT to batch mode
  gROOT->SetBatch(kTRUE);

  // Check for correct number of command line arguments
  if (argc != 5) {
      cout << "Usage: " << argv[0];
      cout << " <data_root_file> <mc_root_file> ";
      cout << " <output_asymmetry_file> <output_kinematic_file>" << endl;
      return 1;
  }

  const char* output_file = argv[3];
  // Clear the contents of the output_file
  std::ofstream ofs(output_file, std::ios::trunc);
  ofs.close();
  const char* kinematic_file = argv[4];
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
  cout << "Total unpolarized (carbon) charge: " << total_charge_carbon << " (nc)."<< endl << endl;

  // Load data and mc root files
  TFile* data_file = new TFile(argv[1], "READ");
  TFile* mc_file = new TFile(argv[2], "READ");
  if (!data_file->IsOpen() || !mc_file->IsOpen()) {
    cout << "Error opening ROOT files (is the location correct?). Exiting." << endl;
    return 2;
  } else {
    cout << "ROOT files opened successfully." << endl;
  }
  
  TTree* data = (TTree*)data_file->Get("PhysicsEvents");
  TTree* mc = (TTree*)mc_file->Get("PhysicsEvents");

  if (!data || !mc) {
    cout << "Error getting trees from ROOT files." << endl;
    return 3;
  } else {
    cout << "Trees successfully extracted from ROOT files." << endl << endl;
  }

  for (size_t i = 0; i < allBins.size(); ++i) {
    cout << "-- Beginning kinematic fits." << endl;
    for (int asymmetry = 0; asymmetry < 1; ++asymmetry){
      switch (asymmetry) {
        case 0: cout << "    Beginning chi2 BSA." << endl; break;
        case 1: cout << "    Beginning chi2 TSA." << endl; break;
        case 2: cout << "    Beginning chi2 DSA." << endl; break;
      }
      performChi2Fits(data, output_file, kinematic_file, binNames[i], asymmetry);
    }
  }

  cout << endl; return 0;
}