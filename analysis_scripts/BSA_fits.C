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


// function to get the polarization value
float getPol(int runnum) {
  float pol; 
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

  // Calculate the back-to-back factor and store it in the data map
  const float M = 0.938272088; // proton mass
  float gamma = (2 * M * data.data["x"]) / sqrt(data.data["Q2"]);
  float epsilon = (1 - data.data["y"] - (0.25) * gamma * gamma * data.data["y"] * data.data["y"]) /
    (1 - data.data["y"] + (0.50) * data.data["y"] * data.data["y"] + (0.25) * gamma * gamma *
    data.data["y"] * data.data["y"]);
  float depolarization_factor = sqrt(1 - epsilon * epsilon);
  data.data["b2b_factor"] = (depolarization_factor * data.data["PTPT"]) / (M * M);

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
    if (currentFits <= 4) {
        return data.data.at("status") <= 1e2;
    } else if (currentFits == 5) {
        return data.data.at("status") == 1e0;
    } else if (currentFits == 6) {
        return data.data.at("status") == 1e1;
    } else if (currentFits == 7) {
        return data.data.at("status") == 1e2;
    } else if (currentFits == 8) {
        return data.data.at("status") == 1e0;
    } else if (currentFits == 9) {
        return data.data.at("status") == 1e1;
    } else if (currentFits == 10) {
        return data.data.at("status") == 1e2;
    } else if (currentFits == 11) {
        return data.data.at("status") <= 1e2 || data.data.at("status") == 1e3;
    } else if (currentFits == 12) {
        return data.data.at("status") <= 1e2 || data.data.at("status") == 1e4;
    } else if (currentFits == 13) {
        return data.data.at("status") <= 1e2 || data.data.at("status") == 1e5;
    }
    return true;
}

// Negative log-likelihood function
void negLogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    // npar: number of parameters
    // gin: an array of derivatives (if needed)
    // f: the value of the function
    // par: an array of the parameter values
    // iflag: a flag (see TMinuit documentation for details)

    // Extract parameters A and B from the input parameter array
    double A = par[0];
    double B = par[1];

    // Initialize variables for counting events (N), positive helicity sum (sum_P), 
    // and negative helicity sum (sum_N)
    double N = 0;
    double sum_N = 0;
    double sum_P = 0;

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
          double Delta_phi = event.data.at("Delta_phi");
          double pol = event.data.at("pol");

          // Check if the helicity is positive or negative and update the corresponding sum
          if (event.data.at("helicity") > 0) {
            sum_P += log(1 + pol * (A * sin(Delta_phi) + B * sin(2 * Delta_phi)));
          } else if (event.data.at("helicity") < 0) {
            sum_N += log(1 - pol * (A * sin(Delta_phi) + B * sin(2 * Delta_phi)));
          }
        }
    }

    // Calculate the negative log-likelihood value and store it in the output variable f
    f = N * log(N) - sum_P - sum_N;
}


void performMLMFits(const char *filename, const char* output_file, const std::string& prefix) {
  // Read the event data from the input file and store it in the global variable gData
  gData = readData(filename, variable_names);

  // Determine the number of bins
  size_t numBins = allBins[currentFits].size() - 1;

  // Initialize TMinuit with 2 parameters (A and B)
  double arglist[10]; arglist[0] = 1;
  int ierflg = 0;
  TMinuit minuit(2);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(negLogLikelihood);

  // Declare string streams for storing the MLM fit results
  std::ostringstream mlmFitsAStream;
  std::ostringstream mlmFitsAScaledStream;
  std::ostringstream mlmFitsBStream;
  std::ostringstream mlmFitsScaledBStream;

  // Initialize the string streams with the output variable names
  mlmFitsAStream << prefix << "MLMFitsA = {";
  mlmFitsAScaledStream << prefix << "MLMFitsScaledA = {";
  mlmFitsBStream << prefix << "MLMFitsB = {";
  mlmFitsScaledBStream << prefix << "MLMFitsScaledB = {";

  // Iterate through each bin
  for (size_t i = 0; i < numBins; ++i) {
    cout << endl << "Beginning MLM fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    currentBin = i;

    // Define the parameters A and B with initial values and limits
    minuit.DefineParameter(0, "A", -0.02, 0.1, -1, 1);
    minuit.DefineParameter(1, "B", 0.00, 0.1, -1, 1);

    // Minimize the negative log-likelihood function
    minuit.Migrad();

    // Calculate the mean values of the current variable and the back-to-back factor (b2b_factor)
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
    double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
    double meanb2b = numEvents > 0 ? sumb2b / numEvents : 0.0;

    // Extract the fitted parameter values and errors
    double A, A_error, B, B_error;
    minuit.GetParameter(0, A, A_error);
    minuit.GetParameter(1, B, B_error);

    // Calculate the scaled parameter values and errors
    double scaled_A = A / meanb2b;
    double scaled_A_error = A_error / meanb2b;
    double scaled_B = B / meanb2b;
    double scaled_B_error = B_error;

    // output to text file
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

  // Determine the variable range for the specified bin
  double varMin = allBins[currentFits][binIndex];
  double varMax = allBins[currentFits][binIndex + 1];

  // Create positive and negative helicity histograms
  TH1D* histPos = new TH1D(Form("%s_pos", histName), "", 24, 0, 2 * TMath::Pi());
  TH1D* histNeg = new TH1D(Form("%s_neg", histName), "", 24, 0, 2 * TMath::Pi());

  // Variables to calculate the mean polarization
  double sumPol = 0;
  int numEvents = 0;

  // Fill the positive and negative helicity histograms
  for (const eventData& event : data) {
    double currentVariable = getEventProperty(event, currentFits);
    if (applyKinematicCuts(event, currentFits) && currentVariable >= varMin && 
      currentVariable < varMax) {
      if (event.data.at("helicity") > 0) {
        histPos->Fill(event.data.at("Delta_phi"));
      } else {
        histNeg->Fill(event.data.at("Delta_phi"));
      }
      // Accumulate polarization and event count for mean polarization calculation
      sumPol += event.data.at("pol");
      numEvents++;
    }
  }

  // Calculate the mean polarization
  double meanPol = sumPol / numEvents;

  // Create the asymmetry histogram
  int numBins = histPos->GetNbinsX();
  TH1D* histAsymmetry = new TH1D(Form("%s_asymmetry", histName), "", 
    numBins, 0, 2 * TMath::Pi());

  // Calculate the asymmetry and its error for each bin, and fill the asymmetry histogram
  for (int iBin = 1; iBin <= numBins; ++iBin) {
    double Np = histPos->GetBinContent(iBin);
    double Nm = histNeg->GetBinContent(iBin);

    // Calculate the asymmetry and error for the current bin
    double asymmetry = (1 / meanPol) * (Np - Nm) / (Np + Nm);
    double error = (2 / meanPol) * std::sqrt(Np * Nm / TMath::Power(Np + Nm, 3));

    // Fill the asymmetry histogram with the calculated values
    histAsymmetry->SetBinContent(iBin, asymmetry);
    histAsymmetry->SetBinError(iBin, error);
  }

  // Delete the temporary positive and negative helicity histograms
  delete histPos;
  delete histNeg;

  // Return the final asymmetry histogram
  return histAsymmetry;
}

// Function to fit the asymmetry histogram
double funcToFit(double* x, double* par) {
  // Retrieve the parameters A and B
  double A = par[0];
  double B = par[1];
  
  // Retrieve the Delta_phi variable from the input x array
  double Delta_phi = x[0];

  // Calculate and return the value of the function for the given Delta_phi and parameters A, B
  return A * sin(Delta_phi) + B * sin(2 * Delta_phi);
}

void performChi2Fits(const char *filename, const char* output_file, const std::string& prefix) {
  // Read data from the input file and store it in the global variable gData
  gData = readData(filename, variable_names);

  // Create a new TF1 object called fitFunction representing the function to fit
  TF1* fitFunction = new TF1("fitFunction", funcToFit, 0, 2 * TMath::Pi(), 2);

  // Initialize string streams to store the results for each bin
  std::ostringstream chi2FitsAStream;
  std::ostringstream chi2FitsAScaledStream;
  std::ostringstream chi2FitsBStream;
  std::ostringstream chi2FitsBScaledStream;

  // Add prefix to each string stream
  chi2FitsAStream << prefix << "chi2FitsA = {";
  chi2FitsAScaledStream << prefix << "chi2FitsScaledA = {";
  chi2FitsBStream << prefix << "chi2FitsB = {";
  chi2FitsBScaledStream << prefix << "chi2FitsScaledB = {";

  // Determine the number of bins
  size_t numBins = allBins[currentFits].size() - 1;

  // Loop over each bin
  for (size_t i = 0; i < numBins; ++i) {
    cout << "Beginning chi2 fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    char histName[32];
    snprintf(histName, sizeof(histName), "hist_%zu", i);

    // Create a histogram for the current bin
    TH1D* hist = createHistogramForBin(gData, histName, i);
    // Fit the histogram using the fitFunction
    hist->Fit(fitFunction, "Q");

    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double sumb2b = 0;
    double numEvents = 0;

    // Loop over all events and calculate the sums and event counts
    for (const eventData& event : gData) {
      double currentVariable = getEventProperty(event, currentFits);
      if (applyKinematicCuts(event, currentFits) && currentVariable >= allBins[currentFits][i] && 
        currentVariable < allBins[currentFits][i + 1]) {
          sumVariable += currentVariable;
          sumb2b += event.data.at("b2b_factor");
          numEvents += 1;
      }
    }
    cout << "Found " << numEvents << " events in this bin." << endl;

    // Calculate the mean values for the variable and b2b factor
    double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
    double meanb2b = numEvents > 0 ? sumb2b / numEvents : 0.0;

    // Get the fitted parameters and their errors
    double A = fitFunction->GetParameter(0);
    double A_error = fitFunction->GetParError(0);
    double B = fitFunction->GetParameter(1);
    double B_error = fitFunction->GetParError(1);

    // Calculate the scaled parameters and their errors
    double scaled_A = A / meanb2b;
    double scaled_A_error = A_error / meanb2b;
    double scaled_B = B / meanb2b;
    double scaled_B_error = B_error / meanb2b;

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
  for (size_t i = 0; i < allBins.size(); ++i) {
  // for (size_t i = 0; i < 1; ++i) {
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