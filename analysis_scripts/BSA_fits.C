#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TH1D.h>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
size_t currentFits = 0;
int n = 1;

std::map<std::string, std::vector<float>> bins_map;
std::vector<std::vector<float>> allBins;
std::vector<std::string> binNames;
std::vector<std::string> variable_names;

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
        variable_names.push_back(var_name);
      }
    } else {
      std::stringstream ss(line);
      std::string bin_name, property;
      std::getline(ss, bin_name, ',');
      std::getline(ss, property, ',');
      binNames.push_back(bin_name);

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

// function to get the polarization value
float getPol(int runnum) {
  float pol; 
    if (runnum == 11 ) { pol = 0.86; } // MC
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
  return pol;
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

  for (const auto& var_name : variable_names) {
    iss >> value;
    data.data[var_name] = value;
  }

  // Debug print statement to ensure "Delta_phi" is present
  if (data.data.find("helicity") != data.data.end()) {
      std::cout << "Delta_phi found: " << data.data["Delta_phi"] << std::endl;
  } else {
      std::cout << "Delta_phi not found" << std::endl;
  }

  int runnum = static_cast<int>(data.data["runnum"]);
  data.data["pol"] = getPol(runnum);

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
  std::string property = binNames[currentFits];

  // Access the property value using the map's indexing
  return event.data.at(property);
}


// Apply kinematic cuts to the data
bool applyKinematicCuts(const eventData& data, int currentFits) {
    // if (data.helicity == 0 || 
    //   data.y > 0.75 || 
    //   data.p2_p < 1.20 || data.p2_p > 4.0 ||
    //   data.z2 < 0.2 || 
    //   data.xF1 > 0.0 || data.xF2 < 0.0 || 
    //   data.Mx < 0.95 || data.Mx1 < 1.35 || data.Mx2 < 1.80) {
    //     return false;
    // }
    // if ((data.runnum >= 5032 && data.runnum <= 5419) || 
    //   (data.runnum >= 6616 && data.runnum <= 6783) || 
    //   (data.runnum >= 6156 && data.runnum <= 6603) || 
    //   (data.runnum >= 11284 && data.runnum <= 11300) || 
    //     (data.runnum >= 11323 && data.runnum <= 11571)) { // inbending data
    //     if (data.vz_e < -8 || data.vz_e > 3 || data.vz_p2 < -10 || 
    //       data.vz_p2 > 2.5 || data.vz_p1 < -10 || data.vz_p1 > 2.5) {
    //         return false;
    //     }
    // } else if ((data.runnum >= 5422 && data.runnum <= 5666) || 
    //   (data.runnum >= 11093 &&  data.runnum <= 11283)) { 
    //     // outbending data
    //     if (data.vz_e < -10 || data.vz_e > 2.5 || data.vz_p2 < -8 || 
    //       data.vz_p2 > 3 || data.vz_p1 < -8 || data.vz_p1 > 3) {
    //         return false;
    //     }
    // }
     if (currentFits <= 4) { return data.data.at("status") <= 1e2;
    } else if (currentFits == 5) { return data.data.at("status") == 1e0;
    } else if (currentFits == 6) { return data.data.at("status") == 1e1;
    } else if (currentFits == 7) { return data.data.at("status") == 1e2;
    } else if (currentFits == 8) { return data.data.at("status") == 1e0;
    } else if (currentFits == 9) { return data.data.at("status") == 1e1;
    } else if (currentFits == 10) { return data.data.at("status") == 1e2;
    } else if (currentFits == 11) { return data.data.at("status") <= 1e2 || 
      data.data.at("status") == 1e3;
    } else if (currentFits == 12) { return data.data.at("status") <= 1e2 || 
      data.data.at("status") == 1e4;
    } else if (currentFits == 13) { return data.data.at("status") <= 1e2 || 
      data.data.at("status") == 1e5;
    } else { return true;
    }
  }


// Negative log-likelihood function
void negLogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    // npar: number of parameters
    // gin: an array of derivatives (if needed)
    // f: the value of the function
    // par: an array of the parameter values
    // iflag: a flag (see TMinuit documentation for details)

    double A = par[0]; double B = par[1];
    double N = 0; double sum_N = 0; double sum_P = 0;

    for (const eventData &event : gData) {
        double currentVariable = getEventProperty(event, currentFits);
        if (applyKinematicCuts(event, currentFits) && 
          currentVariable >= allBins[currentFits][currentBin] && 
          currentVariable < allBins[currentFits][currentBin + 1]) {
          N += 1;
          double Delta_phi = event.data.at("Delta_phi");
          double pol = event.data.at("pol");
          if (event.data.at("helicity") > 0) {
            sum_P += log(1 + pol * (A * sin(Delta_phi) + B * sin(2 * Delta_phi)));
          } else if (event.data.at("helicity") < 0) {
            sum_N += log(1 - pol * (A * sin(Delta_phi) + B * sin(2 * Delta_phi)));
          }
        }
    }
    f = N * log(N) - sum_P - sum_N;
}

void performMLMFits(const char *filename, const char* output_file, const std::string& prefix) {
  gData = readData(filename, variable_names);

  size_t numBins = allBins[currentFits].size() - 1;

  double arglist[10]; arglist[0] = 1; 
  int ierflg = 0;
  TMinuit minuit(2); // 2 is the number of parameters (A and B)
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(negLogLikelihood);

  std::ostringstream mlmFitsAStream; std::ostringstream mlmFitsAScaledStream;
  std::ostringstream mlmFitsBStream; std::ostringstream mlmFitsScaledBStream;

  mlmFitsAStream << prefix << "MLMFitsA = {"; mlmFitsAScaledStream << prefix << "MLMFitsScaledA = {";
  mlmFitsBStream << prefix << "MLMFitsB = {"; mlmFitsScaledBStream << prefix << "MLMFitsScaledB = {";

  for (size_t i = 0; i < numBins; ++i) {
    cout << "Beginning MLM fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    currentBin = i;

    // Define the parameters
    minuit.DefineParameter(0, "A", -0.02, 0.1, -1, 1);
    minuit.DefineParameter(1, "B", 0.00, 0.1, -1, 1);

    // Minimize the function
    minuit.Migrad();
    // minuit.mnexcm("MINOS",arglist,1,ierflg); // running MINOS and HESSE recommended
    // minuit.mnexcm("HESSE",arglist,0,ierflg);

    double sumVariable = 0;
    double sumb2b = 0;
    double numEvents = 0;
    for (const eventData &event : gData) {
      double currentVariable = getEventProperty(event, currentFits);
        if (applyKinematicCuts(event, currentFits) && currentVariable >= 
          allBins[currentFits][i] && currentVariable < allBins[currentFits][i + 1]) {
            sumVariable += currentVariable;
            sumb2b += event.data.at("b2b_factor");
            numEvents += 1;
        }
    }
    cout << "Found " << numEvents << " events." << endl;
    double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
    double meanb2b = numEvents > 0 ? sumb2b / numEvents : 0.0;

    double A, A_error, B, B_error;
    minuit.GetParameter(0, A, A_error);
    minuit.GetParameter(1, B, B_error);

    double scaled_A = A / meanb2b; double scaled_A_error = A_error / meanb2b;
    double scaled_B = B / meanb2b; double scaled_B_error = B_error / meanb2b;

    mlmFitsAStream << "{" << meanVariable << ", " << A << ", " << A_error << "}";
    mlmFitsAScaledStream << "{" << meanVariable << ", " << scaled_A << ", " << 
      scaled_A_error << "}";

    mlmFitsBStream << "{" << meanVariable << ", " << B << ", " << B_error << "}";
    mlmFitsScaledBStream << "{" << meanVariable << ", " << scaled_B << ", " << 
      scaled_B_error << "}";

    if (i < numBins - 1) {
        mlmFitsAStream << ", "; mlmFitsAScaledStream << ", ";
        mlmFitsBStream << ", "; mlmFitsScaledBStream << ", ";
    }
  }

  mlmFitsAStream << "};"; mlmFitsAScaledStream << "};"; 
  mlmFitsBStream << "};"; mlmFitsScaledBStream << "};";

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << mlmFitsAStream.str() << std::endl;
  outputFile << mlmFitsAScaledStream.str() << std::endl;
  outputFile << mlmFitsBStream.str() << std::endl;
  outputFile << mlmFitsScaledBStream.str() << std::endl << std::endl;

  outputFile.close();
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
    cout << "Hey we made it this far" << endl;
    cout << event.data.at("Delta_phi") << endl;
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

  std::ostringstream chi2FitsAStream;
  std::ostringstream chi2FitsAScaledStream;
  std::ostringstream chi2FitsBStream;
  std::ostringstream chi2FitsBScaledStream;

  chi2FitsAStream << prefix << "chi2FitsA = {"; 
  chi2FitsAScaledStream << prefix << "chi2FitsScaledA = {";
  chi2FitsBStream << prefix << "chi2FitsB = {"; 
  chi2FitsBScaledStream << prefix << "chi2FitsScaledB = {";

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
      
      double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
      double meanb2b = numEvents > 0 ? sumb2b / numEvents : 0.0;

      double A = fitFunction->GetParameter(0); double A_error = fitFunction->GetParError(0);
      double B = fitFunction->GetParameter(1); double B_error = fitFunction->GetParError(1);

      double scaled_A = A / meanb2b; double scaled_A_error = A_error / meanb2b;
      double scaled_B = B / meanb2b; double scaled_B_error = B_error / meanb2b;

      chi2FitsAStream << "{" << meanVariable << ", " << A << ", " << A_error << "}";
      chi2FitsAScaledStream << "{" << meanVariable << ", " << scaled_A << ", " << 
        scaled_A_error << "}";
      chi2FitsBStream << "{" << meanVariable << ", " << B << ", " << B_error << "}";
      chi2FitsBScaledStream << "{" << meanVariable << ", " << scaled_B << ", " << 
        scaled_B_error << "}";

      if (i < numBins - 1) {
          chi2FitsAStream << ", ";
          chi2FitsAScaledStream << ", ";
          chi2FitsBStream << ", ";
          chi2FitsBScaledStream << ", ";
      }
      delete hist;
    }

    chi2FitsAStream << "};"; chi2FitsAScaledStream << "};"; 
    chi2FitsBStream << "};"; chi2FitsBScaledStream << "};";

    std::ofstream outputFile(output_file, std::ios_base::app);
    outputFile << chi2FitsAStream.str() << std::endl;
    outputFile << chi2FitsAScaledStream.str() << std::endl;
    outputFile << chi2FitsBStream.str() << std::endl;
    outputFile << chi2FitsBScaledStream.str() << std::endl;

    outputFile.close();
  }

void BSA_fits(const char* data_file, const char* output_file) {

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

  cout<< "Found " << variable_names.size() << " variables: " << endl;
  for (size_t i = 0; i < variable_names.size(); ++i) {
    cout << variable_names[i];
    if (i == variable_names.size() - 1) { cout << "."; } 
    else { cout << ", "; }
  }
  std::cout << std::endl;

  cout << endl << endl;
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