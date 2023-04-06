#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TH1D.h>
#include <algorithm>

std::vector<float> xBins = {0.05, 0.12, 0.20, 0.28, 0.36, 0.44, 0.60};
std::vector<float> zetaBins = {0.30, 0.400, 0.475, 0.550, 0.625, 0.70, 0.80};
std::vector<float> PT1Bins = {0, 0.18, 0.36, 0.54, 0.72, 0.90, 1.30};
std::vector<float> PT2Bins = {0, 0.18, 0.36, 0.54, 0.72, 0.90, 1.30};
std::vector<float> PTPTBins = {0, 0.09, 0.18, 0.27, 0.36, 0.45, 0.6};
std::vector<float> zeta00Bins = {0.40, 0.49, 0.58, 0.67, 0.80};
std::vector<float> zeta20Bins = {0.40, 0.49, 0.58, 0.67, 0.80};
std::vector<float> zeta32Bins = {0.30, 0.37, 0.44, 0.51, 0.60};
std::vector<float> Q200Bins = {1.00, 1.50, 2.00, 2.50, 3.00};
std::vector<float> Q220Bins = {1.00, 2.00, 3.00, 4.00, 5.00};
std::vector<float> Q232Bins = {2.00, 3.50, 5.00, 6.50, 8.00};
std::vector<float> z1Bins = {0.10, 0.20, 0.28, 0.36, 0.44, 0.52, 0.7};
std::vector<float> xF1Bins = {-0.10, 0.00, 0.08, 0.16, 0.26, 0.36, 0.50};
std::vector<float> xF2Bins = {-0.82, -0.60, -0.38, -0.16, 0.06, 0.28, 0.50};

std::vector<std::vector<float>> allBins = {xBins, zetaBins, PT1Bins, PT2Bins, PTPTBins,
  zeta00Bins, zeta20Bins, zeta32Bins, Q200Bins, Q220Bins, Q232Bins, z1Bins, xF1Bins, xF2Bins};
size_t currentFits = 0;
std::vector<std::string> binNames = {"x", "zeta", "PT1", "PT2", "PTPT", "zeta00", "zeta20", 
  "zeta32", "Q200", "Q220", "Q232", "z1", "xF1", "xF2"};

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
  int status, runnum, evnum, helicity;
  float e_p, e_theta, e_phi, vz_e;
  float p2_p, p2_theta, p2_phi, vz_p2;
  float p1_p, p1_theta, p1_phi, vz_p1;
  float Q2, W, x, y, z2, z1;
  float Mx, Mx2, Mx1;
  float zeta, Mh;
  float PT2, PT1, PTPT;
  float xF2, xF1, eta2, eta1, Delta_eta;
  float phi2, phi1, Delta_phi;
  float pol;
  float b2b_factor;
};

std::vector<eventData> gData;
size_t currentBin = 0;

eventData parseLine(const std::string& line) {
  std::istringstream iss(line);
  eventData data;
  iss >> data.status >> data.runnum >> data.evnum >> data.helicity
    >> data.e_p >> data.e_theta >> data.e_phi >> data.vz_e
    >> data.p2_p >> data.p2_theta >> data.p2_phi >> data.vz_p2
    >> data.p1_p >> data.p1_theta >> data.p1_phi >> data.vz_p1
    >> data.Q2 >> data.W >> data.x >> data.y >> data.z2 >> data.z1
    >> data.Mx >> data.Mx2 >> data.Mx1
    >> data.zeta >> data.Mh
    >> data.PT2 >> data.PT1 >> data.PTPT
    >> data.xF2 >> data.xF1 >> data.eta2 >> data.eta1 >> data.Delta_eta
    >> data.phi2 >> data.phi1 >> data.Delta_phi;
    data.pol = getPol(data.runnum);

    // Calculate b2b_factor
    const float M = 0.938272088; // proton mass
    float gamma = (2 * M * data.x) / sqrt(data.Q2);
    float epsilon = (1-event.y-(0.25)*gamma*gamma*event.y*event.y)/(1-y+(0.50)*y*y+(0.25)*gamma*gamma*event.y*event.y)
    float depolarization_factor = sqrt(1-epsilon*epsilon);
    data.b2b_factor = (depolarization_factor*data.PTPT)/(M*M);

  return data;
}

std::vector<eventData> readData(const std::string& filename) {
  std::ifstream infile(filename);
  std::string line;
  std::vector<eventData> data;

  while (std::getline(infile, line)) {
    data.push_back(parseLine(line));
  }

  return data;
}

double getEventProperty(const eventData& event, int currentFits) {
  switch (currentFits) {
    case 0: return event.x;
    case 1: return event.zeta;
    case 2: return event.PT2;
    case 3: return event.PT1;
    case 4: return event.PTPT;
    case 5: return event.zeta;
    case 6: return event.zeta;
    case 7: return event.zeta;
    case 8: return event.Q2;
    case 9: return event.Q2;
    case 10: return event.Q2;
    case 11: return event.z2;
    case 12: return event.xF2;
    case 13: return event.xF1;
    default: return 0.0;
  }
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
    return (currentFits <= 4) ? (data.status <= 1e2) : true; // x, zeta, PT1, PT2, PTPT
    return (currentFits == 5) ? (data.status == 1e0) : true; // 1st zeta-x bin
    return (currentFits == 6) ? (data.status == 1e1) : true; // 2nd zeta-x bin
    return (currentFits == 7) ? (data.status == 1e2) : true; // 3rd zeta-x bin
    return (currentFits == 8) ? (data.status == 1e0) : true; // 1st Q2-x bin
    return (currentFits == 9) ? (data.status == 1e1) : true; // 2nd Q2-x bin
    return (currentFits == 10) ? (data.status == 1e2) : true; // 3rd Q2-x bin
    return (currentFits == 11) ? (data.status <= 1e2 || data.status == 1e3) : true; // z1
    return (currentFits == 12) ? (data.status <= 1e2 || data.status == 1e4) : true; // xF1
    return (currentFits == 13) ? (data.status <= 1e2 || data.status == 1e5) : true; // xF2
}

// Negative log-likelihood function
void negLogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    // npar: number of parameters
    // gin: an array of derivatives (if needed)
    // f: the value of the function
    // par: an array of the parameter values
    // iflag: a flag (see TMinuit documentation for details)

    double A = par[0];
    double B = par[1];

    double N = 0;
    double sum_N = 0;
    double sum_P = 0;

    for (const eventData &event : gData) {
        double currentVariable = getEventProperty(event, currentFits);
        if (applyKinematicCuts(event, currentFits) && 
          currentVariable >= allBins[currentFits][currentBin] && 
          currentVariable < allBins[currentFits][currentBin + 1]) {
          N += 1;
          double Delta_phi = event.Delta_phi;
          double pol = event.pol;
          if (event.helicity > 0) {
            sum_P += log(1 + pol * (A * sin(Delta_phi) + B * sin(2 * Delta_phi)));
          } else if (event.helicity < 0) {
            sum_N += log(1 - pol * (A * sin(Delta_phi) + B * sin(2 * Delta_phi)));
          }
        }
    }

    f = N * log(N) - sum_P - sum_N;
}

void performMLMFits(const char *filename, const char* output_file, const std::string& prefix) {
  gData = readData(filename);

  size_t numBins = allBins[currentFits].size() - 1;

  double arglist[10]; arglist[0] = 1; 
  int ierflg = 0;
  TMinuit minuit(2); // 2 is the number of parameters (A and B)
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(negLogLikelihood);

  std::ostringstream mlmFitsAStream;
  std::ostringstream mlmFitsAScaledStream;
  std::ostringstream mlmFitsBStream;
  std::ostringstream mlmFitsScaledBStream;

  mlmFitsAStream << prefix << "MLMFitsA = {";
  mlmFitsAScaledStream << prefix << "MLMFitsScaledA = {";
  mlmFitsBStream << prefix << "MLMFitsB = {";
  mlmFitsScaledBStream << prefix << "MLMFitsScaledB = {";

  for (size_t i = 0; i < numBins; ++i) {
    cout << endl << "Beginning MLM fit for " << binNames[currentFits]
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
            sumb2b += event.b2b_factor;
            numEvents += 1;
        }
    }
    double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
    double meanb2b = numEvents > 0 ? sumb2b / numEvents : 0.0;

    double A, A_error, B, B_error;
    minuit.GetParameter(0, A, A_error);
    minuit.GetParameter(1, B, B_error);

    double scaled_A = A / meanb2b;
    double scaled_A_error = A_error / meanb2b;
    double scaled_B = B / meanb2b;
    double scaled_B_error = B_error / meanb2b;

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
    double currentVariable = getEventProperty(event, currentFits);
    if (applyKinematicCuts(event, currentFits) && currentVariable >= varMin && 
      currentVariable < varMax) {
      if (event.helicity > 0) {
        histPos->Fill(event.Delta_phi);
      } else {
        histNeg->Fill(event.Delta_phi);
      }
      sumPol += event.pol;
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
  gData = readData(filename);

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
      cout << endl << "Beginning chi2 fit for " << binNames[currentFits]
        << " bin " << i << ". ";
      char histName[32];
      snprintf(histName, sizeof(histName), "hist_%zu", i);

      TH1D* hist = createHistogramForBin(gData, histName, i);
      hist->Fit(fitFunction, "Q");

      double sumVariable = 0;
      double sumb2b = 0;
      double numEvents = 0;
      double sumPTPT = 0;
      double sumDepolarization_factor = 0;
      for (const eventData& event : gData) {
        double currentVariable = getEventProperty(event, currentFits);
        if (applyKinematicCuts(event, currentFits) && currentVariable >= allBins[currentFits][i] && 
          currentVariable < allBins[currentFits][i + 1]) {
            sumVariable += currentVariable;
            sumb2b += event.b2b_factor;
            numEvents += 1;
            sumPTPT += event.PTPT;
            float gamma = (2 * 0.938272088 * event.x) / sqrt(event.Q2);
            float epsilon = (1-event.y-(0.25)*gamma*gamma*event.y*event.y)/(1-y+(0.50)*y*y+(0.25)*gamma*gamma*event.y*event.y);
            float depolarization_factor = sqrt(1-epsilon*epsilon);
            sumDepolarization_factor += depolarization_factor;
        }
      }
      
      double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
      double meanb2b = numEvents > 0 ? sumb2b / numEvents : 0.0;
      double meanPTPT = numEvents > 0 ? sumPTPT / numEvents : 0.0;
      double meanDepolarization_factor = numEvents > 0 ? sumDepolarization_factor / numEvents : 0.0;

      cout << numEvents << " " << meanPTPT << " " << meanDepolarization_factor << " " << meanb2b << endl;

      double A = fitFunction->GetParameter(0);
      double A_error = fitFunction->GetParError(0);
      double B = fitFunction->GetParameter(1);
      double B_error = fitFunction->GetParError(1);

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

  // for (size_t i = 0; i < allBins.size(); ++i) {
  for (size_t i = 0; i < 1; ++i) {
    performChi2Fits(data_file, output_file, binNames[i]);
    cout << endl << "     Completed " << binNames[i] << " chi2 fits." << endl;
    // performMLMFits(data_file, output_file, binNames[i]);
    // cout << endl << "     Completed " << binNames[i] << " MLM fits." << endl;
    cout << endl << endl << endl;
    currentFits++;
  }
}


