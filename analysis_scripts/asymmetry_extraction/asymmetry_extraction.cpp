// Standard C++ Library Headers
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
// ROOT Library Headers
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TSystem.h>
#include <TText.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include "TH2.h"
#include "TMinuit.h"
#include "TStyle.h"
// tbhayward libraries
#include "common_vars.h"  // Include the common header
#include "load_bins_from_csv.h"
#include "load_run_info_from_csv.h"
#include "dilution_factor.h"
#include "asymmetry_fits.h"
#include "KinematicCuts.h"
#include "formatLabelName.h"
#include "readChi2Fits.h"
#include "histConfigs.h"
#include "charge_accumulation.h"

// Using namespace declaration
using namespace std;

TTreeReader dataReader;  // Declare as global variable
TTreeReader mcReader;  // Declare as global variable

int currentFits = 0;
size_t currentBin = 0;
int n = 1;
double total_charge_carbon = 0;
double cmm = 0; 
double cpm = 0; 
double cmp = 0; 
double cpp = 0; 
std::string mlmPrefix = "xF";

int num_data_elec = 626802;
int num_mc_elec = 300366;


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

  TTreeReaderValue<int> runnum(dataReader, "runnum");
  TTreeReaderValue<int> evnum(dataReader, "evnum");
  TTreeReaderValue<double> xF(dataReader, "xF");
  TTreeReaderValue<double> Mx(dataReader, "Mx");
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> phi(dataReader, "phi");
  // TTreeReaderValue<double> phi(dataReader, "phi23");
  TTreeReaderValue<double> DepA(dataReader, "DepA");
  TTreeReaderValue<double> DepB(dataReader, "DepB");
  TTreeReaderValue<double> DepC(dataReader, "DepC");
  TTreeReaderValue<double> DepV(dataReader, "DepV");
  TTreeReaderValue<double> DepW(dataReader, "DepW");
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  KinematicCuts data_kinematicCuts(dataReader);  // Create an instance of the KinematicCuts class
  while (dataReader.Next()) {
    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = data_kinematicCuts.applyCuts(currentFits, false);
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= allBins[currentFits][currentBin] && 
          *currentVariable < allBins[currentFits][currentBin + 1] && passedKinematicCuts) {

      // Increment the event count
      N += 1;

      double Df = dilution_factor(*currentVariable, mlmPrefix); // dilution factor
      double Pb = *beam_pol;
      double Pt = std::abs(*target_pol);

      if (*helicity > 0 && *target_pol > 0) { 
        sum_PP = sum_PP + log(1 
          + (*DepV / *DepA)*AUU_cosphi*cos(*phi) + (*DepB / *DepA)*AUU_cos2phi*cos(2 * *phi) // UU 
          + Pb*((*DepW / *DepA)*ALU_sinphi*sin(*phi)) // BSA
          + Df*Pt*((*DepV / *DepA)*AUL_sinphi*sin(*phi)+ // TSA
            (*DepB / *DepA)*AUL_sin2phi*sin(2 * *phi))//TSA
          + Df*Pb*Pt*((*DepC / *DepA)*ALL + (*DepW / *DepA)*ALL_cosphi*cos(*phi)) ); // DSA
      } else if (*helicity > 0 && *target_pol < 0) { 
        sum_PM = sum_PM + log(1 
          + (*DepV / *DepA)*AUU_cosphi*cos(*phi) + (*DepB / *DepA)*AUU_cos2phi*cos(2 * *phi) // UU
          + Pb*((*DepW / *DepA)*ALU_sinphi*sin(*phi)) // BSA
          - Df*Pt*((*DepV / *DepA)*AUL_sinphi*sin(*phi)+ // TSA
            (*DepB / *DepA)*AUL_sin2phi*sin(2 * *phi)) // TSA
          - Df*Pb*Pt*((*DepC / *DepA)*ALL + (*DepW / *DepA)*ALL_cosphi*cos(*phi)) ); // DSA
      } else if (*helicity < 0 && *target_pol > 0) { 
        sum_MP = sum_MP + log(1 
          + (*DepV / *DepA)*AUU_cosphi*cos(*phi) + (*DepB / *DepA)*AUU_cos2phi*cos(2 * *phi) // UU 
          - Pb*((*DepW / *DepA)*ALU_sinphi*sin(*phi)) // BSA
          + Df*Pt*((*DepV / *DepA)*AUL_sinphi*sin(*phi)+ // TSA
            (*DepB / *DepA)*AUL_sin2phi*sin(2 * *phi))//TSA
          - Df*Pb*Pt*((*DepC / *DepA)*ALL + (*DepW / *DepA)*ALL_cosphi*cos(*phi)) ); // DSA
      } else if (*helicity < 0 && *target_pol < 0) { 
        sum_MM = sum_MM + log(1 
          + (*DepV / *DepA)*AUU_cosphi*cos(*phi) + (*DepB / *DepA)*AUU_cos2phi*cos(2 * *phi) // UU 
          - Pb*((*DepW / *DepA)*ALU_sinphi*sin(*phi)) // BSA
          - Df*Pt*((*DepV / *DepA)*AUL_sinphi*sin(*phi)+ // TSA
            (*DepB / *DepA)*AUL_sin2phi*sin(2 * *phi))//TSA
          + Df*Pb*Pt*((*DepC / *DepA)*ALL + (*DepW / *DepA)*ALL_cosphi*cos(*phi)) ); // DSA
      }
    }
  }
  dataReader.Restart();  // Reset the TTreeReader at the end of the function
  
  TTreeReaderValue<double> mc_phi(mcReader, "phi");
  // TTreeReaderValue<double> mc_phi(mcReader, "phi23");
  TTreeReaderValue<double> mc_DepA(mcReader, "DepA");
  TTreeReaderValue<double> mc_DepB(mcReader, "DepB");
  TTreeReaderValue<double> mc_DepC(mcReader, "DepC");
  TTreeReaderValue<double> mc_DepV(mcReader, "DepV");
  TTreeReaderValue<double> mc_DepW(mcReader, "DepW");
  TTreeReaderValue<double> mc_currentVariable(mcReader, propertyNames[currentFits].c_str());

  KinematicCuts mc_kinematicCuts(mcReader);  // Create an instance of the KinematicCuts class
  while (mcReader.Next()) {
    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = mc_kinematicCuts.applyCuts(currentFits, true);
    // Check if the currentVariable is within the desired range
    if (*mc_currentVariable >= allBins[currentFits][currentBin] && 
          *mc_currentVariable < allBins[currentFits][currentBin + 1] && passedKinematicCuts) {
      NUU+=1+(*mc_DepV / *mc_DepA)*AUU_cosphi*cos(*mc_phi)+
        (*mc_DepB / *mc_DepA)*AUU_cos2phi*cos(2 * *mc_phi); // UU
    }
  }
  mcReader.Restart();  // Reset the TTreeReader at the end of the function

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

void performMLMFits(const char* output_file, const char* kinematic_file,
  const std::string& prefix) {
  // Read the event data from the input file and store it in the global variable gData
  mlmPrefix = prefix;

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

    // Read the chi2 fits into a map
    std::map<std::string, std::vector<std::vector<double>>> chi2Fits = 
        readChi2Fits(std::string(output_file));

    // Construct the key based on the prefix and the fit name
    // For now, let's assume fitName is a string that contains the fit name like "ALUsinphi"
    std::string fitName = "ALUsinphi";// Replace this with the logic to determine the fit name
    std::string key = std::string(prefix) + "chi2Fits" + fitName; 

    std::vector<double> chi2Result = chi2Fits[key][currentFits];
    // Define the parameters with initial values and limits
    minuit.DefineParameter(0, "ALU_sinphi", 
      chi2Fits[std::string(prefix)+"chi2FitsALUsinphi"][currentFits][1], 0.01, -1, 1);
    minuit.DefineParameter(1, "AUL_sinphi", 
      chi2Fits[std::string(prefix)+"chi2FitsAULsinphi"][currentFits][1], 0.01, -1, 1);
    minuit.DefineParameter(2, "AUL_sin2phi", 
      chi2Fits[std::string(prefix)+"chi2FitsAULsin2phi"][currentFits][1], 0.01, -1, 1);
    minuit.DefineParameter(3, "ALL", 
      chi2Fits[std::string(prefix)+"chi2FitsALL"][currentFits][1], 0.01, -1, 1);
    minuit.DefineParameter(4, "ALL_cosphi", 
      chi2Fits[std::string(prefix)+"chi2FitsALLcosphi"][currentFits][1], 0.01, -1, 1);
    minuit.DefineParameter(5, "AUU_cosphi", -0.1, 0.01, -1, 1);
    minuit.DefineParameter(6, "AUU_cos2phi", 0.10, 0.01, -1, 1);

    // Minimize the negative log-likelihood function
    minuit.Migrad(); cout << endl;

    // Extract the fitted parameter values and errors
    double ALU_sinphi, ALU_sinphi_error; minuit.GetParameter(0, ALU_sinphi, ALU_sinphi_error);
    double AUL_sinphi, AUL_sinphi_error; minuit.GetParameter(1, AUL_sinphi, AUL_sinphi_error);
    double AUL_sin2phi, AUL_sin2phi_error; minuit.GetParameter(2, AUL_sin2phi, AUL_sin2phi_error);
    double ALL, ALL_error; minuit.GetParameter(3, ALL, ALL_error);
    double ALL_cosphi, ALL_cosphi_error; minuit.GetParameter(4, ALL_cosphi, ALL_cosphi_error);
    double AUU_cosphi, AUU_cosphi_error; minuit.GetParameter(5, AUU_cosphi, AUU_cosphi_error);
    double AUU_cos2phi, AUU_cos2phi_error; minuit.GetParameter(6, AUU_cos2phi, AUU_cos2phi_error);

    // Calculate the mean values of the current variable 
    double sumVariable = 0;
    double numEvents = 0;
    KinematicCuts data_kinematicCuts(dataReader);  // Create an instance of the KinematicCuts class
    TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = data_kinematicCuts.applyCuts(currentFits, false);
      // bool passedKinematicCuts = true;
      // Check if the currentVariable is within the desired range
      if (*currentVariable >= allBins[currentFits][i] && 
        *currentVariable < allBins[currentFits][i + 1] && passedKinematicCuts) {
        sumVariable += *currentVariable;
        numEvents += 1;
      }
    }
    dataReader.Restart();  // Reset the TTreeReader at the end of the function
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

  // Create the legend
  // TLegend *leg = new TLegend(0.16171, 0.7, 0.4, 0.9);  // Adjusted to the upper-left corner
  TLegend *leg = new TLegend(0.19, 0.675, 0.45, 0.875);  // Adjusted to the upper-left corner
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->SetTextSize(0.025);  // Reduced text size
  // leg->SetTextAlign(12);  // Left-align text

  // Add fit parameters as legend entries based on the value of 'asymmetry'.
  const char* paramName;
  for (int i = 0; i < fitFunction->GetNpar(); ++i) {
      if (i == 0 && (asymmetryIndex == 0 || asymmetryIndex == 1)) {
        paramName = "offset";
      } else if (i == 0 && asymmetryIndex == 2) {
        paramName = "#it{A}_{LL}";
      } else if (asymmetryIndex == 0) {
        if (i == 1) paramName = "#it{A}_{LU}^{sin#phi}";
      } else if (asymmetryIndex == 1) {
        if (i == 1) paramName = "#it{A}_{UL}^{sin#phi}";
        if (i == 2) paramName = "#it{A}_{UL}^{sin2#phi}";
      } else if (asymmetryIndex == 2) {
        if (i == 1) paramName = "#it{A}_{LL}^{cos#phi}";
      }
      leg->AddEntry((TObject*)0, Form("%s: %.4f #pm %.4f", paramName, 
        fitFunction->GetParameter(i), fitFunction->GetParError(i)), "");
  }

  // Add the chi-squared per degree of freedom to the legend
  leg->AddEntry((TObject*)0, Form("#chi^{2}/Ndf: %.4f", 
    fitFunction->GetChisquare() / fitFunction->GetNDF()), "");

  // Draw the legend
  leg->Draw("same");

  // Create the filename for the PNG
  string filename = "output/" + prefix + "_" + 
    fileNameSuffix + "_" + std::to_string(binIndex) + ".png";
  
  // Determine the variable range for the specified bin
  double varMin = allBins[currentFits][binIndex];
  double varMax = allBins[currentFits][binIndex + 1];
  // Create a title string for the graph 
  string formattedLabelName = formatLabelName(prefix);
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(3) << varMin << " #leq ";
  oss << formattedLabelName << " < " << varMax;
  std::string title = oss.str();

  // Set the title to the title string
  graph->SetTitle(title.c_str());

  // Save the canvas as a PNG
  canvas->SaveAs(filename.c_str());

  // Clean up
  delete canvas;
  delete graph;
}

TH1D* createHistogramForBin(const char* histName, int binIndex, 
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

  TTreeReaderValue<int> runnum(dataReader, "runnum");
  TTreeReaderValue<int> evnum(dataReader, "evnum");
  TTreeReaderValue<double> xF(dataReader, "xF");
  TTreeReaderValue<double> Mx(dataReader, "Mx");
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> phi(dataReader, "phi");
  // TTreeReaderValue<double> phi(dataReader, "phi23");
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  KinematicCuts kinematicCuts(dataReader);  // Create an instance of the KinematicCuts class
  // Counter to limit the number of processed entries
  while (dataReader.Next()) {

    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts.applyCuts(currentFits, false);
    // bool passedKinematicCuts = true;
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
      sumVariable += *currentVariable;

      if (*helicity > 0 && *target_pol < 0) { histPosNeg->Fill(*phi); } 
      else if (*helicity < 0 && *target_pol > 0) {  histNegPos->Fill(*phi); }

      if (*helicity > 0 && (*target_pol > 0 || *runnum < 11571) ) { histPosPos->Fill(*phi); } 
      else if (*helicity < 0 && (*target_pol < 0 || *runnum < 11571) ) {  histNegNeg->Fill(*phi); } 
      // this structure allows the same script to run for both polarized and unpolarized targets
      // if it is an RGC run with a polarized target (runnum > 11571) then we assign all four
      // combinations, if it is an earlier experiment then we only assign PosPos and NegNeg
      // and set the Ptp and Ptm below to 1, this allows for a regular BSA calculation


      // Accumulate polarization and event count for mean polarization calculation
      sumPol += *beam_pol;
      if (*target_pol > 0) {
        sumTargetPosPol += *target_pol;
        numEventsPosTarget++;
      } else if (*target_pol < 0) {
        sumTargetNegPol += *target_pol;
        numEventsNegTarget++;
      }
      numEvents++; // Increment the numEvents
    }
  }
  dataReader.Restart();  // Reset the TTreeReader at the end of the function

  // Calculate the mean polarization
  double meanVariable = numEvents > 0 ? sumVariable / numEvents : 0.0;
  double meanPol = sumPol / numEvents; // mean beam polarization for data 
  double Ptp = numEventsPosTarget > 0 ? sumTargetPosPol / numEventsPosTarget : 1;
  double Ptm = numEventsNegTarget > 0 ? -sumTargetNegPol / numEventsNegTarget : 1;
  // the negative sign here is correct; RGC lists the polarizations with signs to tell which is 
  // which but the polarization really should just be "percent of polarized nucleii"
  // the sign gives the helicity

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

void performChi2Fits(const char* output_file, const char* kinematic_file,
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
    TH1D* hist = createHistogramForBin(histName, i, prefix, asymmetry_index);
    // Fit the histogram using the fitFunction and get the fit result
    hist->Fit(fitFunction, "QS");
    plotHistogramAndFit(hist, fitFunction, i, asymmetry_index, prefix);

    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double numEvents = 0;
    // Variables to calculate the mean depolarization factor
    double sumDepA = 0; double sumDepB = 0; 
    double sumDepC = 0; double sumDepV = 0; double sumDepW = 0;

    // Variables to calculate the mean kinematics in each bin
    double sumQ2 = 0; double sumW = 0; double sumx = 0; double sumy = 0;
    double sumz = 0; double sumzeta = 0; double sumpT = 0; double sumxF = 0;
    double sumt = 0; double sumtmin = 0;

    // Declare reader locations
    TTreeReaderValue<int> runnum(dataReader, "runnum");
    TTreeReaderValue<int> evnum(dataReader, "evnum");
    TTreeReaderValue<double> Q2(dataReader, "Q2");
    TTreeReaderValue<double> W(dataReader, "W");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> y(dataReader, "y");
    TTreeReaderValue<double> z(dataReader, "z");
    TTreeReaderValue<double> zeta(dataReader, "zeta");
    TTreeReaderValue<double> pT(dataReader, "pT");
    TTreeReaderValue<double> xF(dataReader, "xF");
    TTreeReaderValue<double> Mx(dataReader, "Mx");
    TTreeReaderValue<double> t(dataReader, "t");
    // TTreeReaderValue<double> t(dataReader, "t1");
    TTreeReaderValue<double> tmin(dataReader, "tmin");
    TTreeReaderValue<double> DepA(dataReader, "DepA");
    TTreeReaderValue<double> DepB(dataReader, "DepB");
    TTreeReaderValue<double> DepC(dataReader, "DepC");
    TTreeReaderValue<double> DepV(dataReader, "DepV");
    TTreeReaderValue<double> DepW(dataReader, "DepW");
    TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

    // Determine the variable range for the specified bin
    double varMin = allBins[currentFits][i];
    double varMax = allBins[currentFits][i + 1];
    int counter = 0;
    KinematicCuts kinematicCuts(dataReader);  // Create an instance of the KinematicCuts class
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = kinematicCuts.applyCuts(currentFits, false);
      // Check if the currentVariable is within the desired range
      if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
        // sum the kinematic variable values
        sumVariable += *currentVariable;
        sumQ2 += *Q2;
        sumW += *W;
        sumx += *x;
        sumy += *y;
        sumz += *z;
        sumzeta += *zeta;
        sumpT += *pT;
        sumxF += *xF;
        sumt += *t;
        sumtmin += *tmin;

        // sum the depolarization values
        sumDepA += *DepA;
        sumDepB += *DepB;
        sumDepC += *DepC;
        sumDepV += *DepV;
        sumDepW += *DepW;

        numEvents += 1; 
        counter++;
      }

    }
    dataReader.Restart();  // Reset the TTreeReader at the end of the function
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
        // chi2FitsCStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", "<<AUU_cosphi_error <<"}";
        // chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", "<<AUU_cos2phi_error <<"}";
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
        // chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", "<<AUU_cosphi_error <<"}";
        // chi2FitsEStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", "<<AUU_cos2phi_error <<"}";
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
        // chi2FitsCStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", "<<AUU_cosphi_error <<"}";
        // chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", "<<AUU_cos2phi_error <<"}";
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

template<typename T>
void FillHistogram(TTreeReader& reader, const std::string& branchName, TH1D* hist, 
  KinematicCuts& kinematicCuts, int fitIndex) {
    TTreeReaderValue<T> val(reader, branchName.c_str());
    while (reader.Next()) {
        if (kinematicCuts.applyCuts(fitIndex, false)) {
            hist->Fill(*val);
        }
    }
}


void createIntegratedKinematicPlots() {
    const std::string outputDir = "output/integrated_plots/";
    const std::vector<std::string> branchesToSkip = {"helicity", "beam_pol", "target_pol", "DepA", "DepB", "DepC", "DepV", "DepW", "evnum"};

    TObjArray* branches = dataReader.GetTree()->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: Unable to retrieve branch list from data TTree." << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetTextSize(0.05); // Increase the text size globally
    bool restart = true;
    for (Int_t i = 0; i < branches->GetEntries(); ++i) {

        TBranch* branch = (TBranch*)branches->At(i);
        std::string branchName = branch->GetName();
        if (branchName == "e_p" && restart) {
          // stupid hack to get it to do the runnum plot instead of it being blank 
          // due to reader restarts
          i = 0; 
          restart = false;
        }
        branch = (TBranch*)branches->At(i);
        branchName = branch->GetName();

        if (std::find(branchesToSkip.begin(), branchesToSkip.end(), branchName) != 
          branchesToSkip.end()) {
            continue; // Skip this branch
        }

        // TTreeReaderValue<Double_t> dataVal(dataReader, branchName.c_str());
        // TTreeReaderValue<Double_t> mcVal(mcReader, branchName.c_str());

        HistConfig config = {100, 0, 1}; // Default configuration
        if (histConfigs.find(branchName) != histConfigs.end()) {
            config = histConfigs[branchName];
        }

        TH1D* dataHist = new TH1D((branchName + "_data").c_str(), "", config.nBins, config.xMin, config.xMax);
        TH1D* mcHist = new TH1D((branchName + "_mc").c_str(), "", config.nBins, config.xMin, config.xMax);

        // Set x-axis title
        dataHist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
        mcHist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());

        dataHist->GetXaxis()->SetTitleSize(0.05); // Increase the x-axis title font size
        dataHist->GetYaxis()->SetTitleSize(0.05); // Increase the y-axis title font size
        mcHist->GetXaxis()->SetTitleSize(0.05);
        mcHist->GetYaxis()->SetTitleSize(0.05);

        // Set y-axis title and center it
        dataHist->GetYaxis()->SetTitle("Normalized Counts");
        dataHist->GetYaxis()->CenterTitle();
        mcHist->GetYaxis()->SetTitle("Normalized Counts");
        mcHist->GetYaxis()->CenterTitle();

        // Set y-axis title offset to make room for centering
        dataHist->GetYaxis()->SetTitleOffset(1.6);
        mcHist->GetYaxis()->SetTitleOffset(1.6);

        // // Loop over dataReader to fill the histogram
        // KinematicCuts kinematicCuts(dataReader);
        // while (dataReader.Next()) {
        //     bool passedKinematicCuts = kinematicCuts.applyCuts(0, false);
        //     if (*dataVal >= config.xMin && *dataVal < config.xMax && passedKinematicCuts) {
        //         dataHist->Fill(*dataVal);
        //     }
        // }

        // // Loop over mcReader to fill the histogram
        // KinematicCuts mc_kinematicCuts(mcReader);
        // while (mcReader.Next()) {
        //     bool passedKinematicCuts = mc_kinematicCuts.applyCuts(0, true);
        //     if (*mcVal >= config.xMin && *mcVal < config.xMax && passedKinematicCuts) {
        //         mcHist->Fill(*mcVal);
        //     }
        // }

        KinematicCuts dataKinematicCuts(dataReader);
        KinematicCuts mcKinematicCuts(mcReader);

        if (branchName == "runnum") {
          // Declare TTreeReaderValue for integers for dataReader
          TTreeReaderValue<int> dataVal(dataReader, branchName.c_str());

          // Fill histogram for dataReader
          FillHistogram<int>(dataReader, branchName, dataHist, dataKinematicCuts, 0);

          // Check if the "runnum" branch exists in mcReader
          if (mcReader.GetTree()->GetBranch(branchName.c_str())) {
              // "runnum" branch exists, declare TTreeReaderValue for mcReader
              TTreeReaderValue<int> mcVal(mcReader, branchName.c_str());

              // Fill histogram for mcReader
              FillHistogram<int>(mcReader, branchName, mcHist, mcKinematicCuts, 0);
          } else {
              // "runnum" branch does not exist, use default value
              int defaultRunNum = 11;
              mcHist->Fill(defaultRunNum);
          }
        } else {
          // Declare TTreeReaderValue for doubles
          TTreeReaderValue<double> dataVal(dataReader, branchName.c_str());
          TTreeReaderValue<double> mcVal(mcReader, branchName.c_str());

          // Fill histograms for double values
          FillHistogram<double>(dataReader, branchName, dataHist, dataKinematicCuts, 0);
          FillHistogram<double>(mcReader, branchName, mcHist, mcKinematicCuts, 0);
        }

        // Normalize the histograms
        dataHist->Scale(1.0 / dataHist->Integral());
        mcHist->Scale(1.0 / mcHist->Integral());

        // // Normalize the histograms
        // dataHist->Scale(1.0 / num_data_elec);
        // mcHist->Scale(1.0 / num_mc_elec);

        // Find the maximum value for y-axis
        double maxY = 1.2*std::max(dataHist->GetMaximum(), mcHist->GetMaximum());
        dataHist->SetMaximum(1.1 * maxY);
        mcHist->SetMaximum(1.1 * maxY);

        // Create a canvas for drawing the histograms
        TCanvas* c = new TCanvas((branchName + "_canvas").c_str(), branchName.c_str(), 800, 600);
        // Adjust the margins to avoid cutting off labels
        c->SetLeftMargin(0.15);
        c->SetBottomMargin(0.15);

        // Create a legend and adjust its font size
        TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9);
        leg->SetTextSize(0.04); // Increase the legend text size

        // Add entries to the legend with scientific notation for the number of entries
        dataHist->SetEntries(dataHist->GetEntries());
        mcHist->SetEntries(mcHist->GetEntries());
        leg->AddEntry(dataHist, (std::string("Data (") + std::to_string((int)dataHist->GetEntries()) + " entries)").c_str(), "l");
        leg->AddEntry(mcHist, (std::string("MC (") + std::to_string((int)mcHist->GetEntries()) + " entries)").c_str(), "l");

        // Set line colors for histograms
        dataHist->SetLineColor(kBlack);
        mcHist->SetLineColor(kRed);

        // Draw histograms on the canvas
        dataHist->Draw("HIST");
        mcHist->Draw("HISTSAME");
        leg->Draw();

        // Save the canvas to a file
        c->SaveAs((outputDir + branchName + ".png").c_str());

        // Clean up the created objects to avoid memory leaks
        delete dataHist;
        delete mcHist;
        delete c;
        delete leg;

        // Restart the TTreeReaders for the next branch
        dataReader.Restart();
        mcReader.Restart();
    }
}


void createIntegratedKinematicPlotsForBinsAndFits() {
    const std::string outputDir = "output/binned_plots/";
    const std::vector<std::string> branchesToSkip = {
        "helicity", "beam_pol", "target_pol", "DepA", "DepB", "DepC", "DepV", "DepW", "evnum"
    };

    TObjArray* branches = dataReader.GetTree()->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: Unable to retrieve branch list from data TTree." << std::endl;
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetTextSize(0.05); // Increase the text size globally

    // Loop over all sets of supplied kinematic variables
    for (size_t fitIndex = 0; fitIndex < allBins.size(); ++fitIndex) {
        std::string currentVariable = binNames[fitIndex]; // Assuming binNames is a vector<string> with the same size as allBins
        std::string branchVariable = propertyNames[fitIndex]; // Corresponding branch name for the current variable

        // Loop over all possible bins within the current set
        for (size_t binIndex = 0; binIndex < allBins[fitIndex].size() - 1; ++binIndex) {
            double binLowerEdge = allBins[fitIndex][binIndex];
            double binUpperEdge = allBins[fitIndex][binIndex + 1];

            // Format the bin edges with three decimal places
            std::ostringstream lowerEdgeStream, upperEdgeStream;
            lowerEdgeStream << std::fixed << std::setprecision(3) << binLowerEdge;
            upperEdgeStream << std::fixed << std::setprecision(3) << binUpperEdge;

            std::string binIndexLabel = "bin_" + std::to_string(binIndex + 1);

            // Now we iterate over all branches, except those we wish to skip
            for (Int_t i = 0; i < branches->GetEntries(); ++i) {
                TBranch* branch = (TBranch*)branches->At(i);
                std::string branchName = branch->GetName();

                if (std::find(branchesToSkip.begin(), branchesToSkip.end(), branchName) != branchesToSkip.end()) {
                    continue; // Skip this branch
                }

                // Determine histogram configuration, default or specific
                HistConfig config = {100, 0, 1}; // Default configuration
                if (histConfigs.find(branchName) != histConfigs.end()) {
                    config = histConfigs[branchName];
                } else {
                    std::cerr << "Warning: No specific histogram configuration found for " << branchName << ". Using default configuration." << std::endl;
                }

                // // Set up the data and MC values to be read from the trees
                // TTreeReaderValue<Double_t> dataVal(dataReader, branchName.c_str());
                // TTreeReaderValue<Double_t> mcVal(mcReader, branchName.c_str());
                // TTreeReaderValue<Double_t> binVariable(dataReader, branchVariable.c_str());
                // TTreeReaderValue<Double_t> mcBinVariable(mcReader, branchVariable.c_str());

                // Create histogram title with formatted bin edges
                std::string formattedVariableName = formatLabelName(branchVariable);
                std::string plotTitle = lowerEdgeStream.str() + " < " + formattedVariableName + " < " + upperEdgeStream.str();

                // Create histograms with titles reflecting the bin edges and plotted variable
                std::string histName = currentVariable + "_" + branchName + "_" + binIndexLabel;
                TH1D* dataHist = new TH1D((histName + "_data").c_str(), plotTitle.c_str(), config.nBins, config.xMin, config.xMax);
                TH1D* mcHist = new TH1D((histName + "_mc").c_str(), plotTitle.c_str(), config.nBins, config.xMin, config.xMax);

                // Set histogram titles and styles
                dataHist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
                mcHist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());

                dataHist->GetYaxis()->SetTitle("Normalized Counts");
                mcHist->GetYaxis()->SetTitle("Normalized Counts");

                dataHist->GetYaxis()->CenterTitle();
                mcHist->GetYaxis()->CenterTitle();

                dataHist->GetYaxis()->SetTitleOffset(1.6);
                mcHist->GetYaxis()->SetTitleOffset(1.6);

                dataHist->SetLineColor(kBlack);
                mcHist->SetLineColor(kRed);

                KinematicCuts kinematicCuts(dataReader);
                KinematicCuts mcKinematicCuts(mcReader);

                if (branchName == "runnum") {
                  TTreeReaderValue<int> dataVal(dataReader, branchName.c_str());
                  TTreeReaderValue<double> binVariable(dataReader, branchVariable.c_str());

                  if (mcReader.GetTree()->GetBranch(branchName.c_str())) {
                      TTreeReaderValue<int> mcVal(mcReader, branchName.c_str());
                      TTreeReaderValue<double> mcBinVariable(mcReader, branchVariable.c_str());
                      FillHistogram<int>(dataReader, branchName, dataHist, kinematicCuts, fitIndex);
                      FillHistogram<int>(mcReader, branchName, mcHist, mcKinematicCuts, fitIndex);
                  } else {
                      int defaultRunNum = 11;
                      FillHistogram<int>(dataReader, branchName, dataHist, kinematicCuts, fitIndex);
                      mcHist->Fill(defaultRunNum);
                  }
                } else {
                  TTreeReaderValue<double> dataVal(dataReader, branchName.c_str());
                  TTreeReaderValue<double> binVariable(dataReader, branchVariable.c_str());
                  TTreeReaderValue<double> mcVal(mcReader, branchName.c_str());
                  TTreeReaderValue<double> mcBinVariable(mcReader, branchVariable.c_str());
                  FillHistogram<double>(dataReader, branchName, dataHist, kinematicCuts, fitIndex);
                  FillHistogram<double>(mcReader, branchName, mcHist, mcKinematicCuts, fitIndex);
                }

                // Normalize the histograms
                if (dataHist->Integral() != 0) {
                    dataHist->Scale(1.0 / dataHist->Integral());
                }
                if (mcHist->Integral() != 0) {
                    mcHist->Scale(1.0 / mcHist->Integral());
                }

                // // Normalize the histograms
                // if (dataHist->Integral() != 0) {
                //     dataHist->Scale(1.0 / num_data_elec);
                // }
                // if (mcHist->Integral() != 0) {
                //     mcHist->Scale(1.0 / num_mc_elec);
                // }

                // Find the maximum y-value between both histograms to set the y-axis range
                double maxY = std::max(dataHist->GetMaximum(), mcHist->GetMaximum());
                dataHist->SetMaximum(1.2 * maxY);
                mcHist->SetMaximum(1.2 * maxY);

                // Create a canvas for drawing the histograms
                TCanvas* c = new TCanvas((histName + "_canvas").c_str(), branchName.c_str(), 800, 600);
                c->SetLeftMargin(0.15);
                c->SetBottomMargin(0.15);

                // Draw histograms on the canvas
                dataHist->Draw("HIST");
                mcHist->Draw("HIST SAME");

                // Create a legend for the histograms
                TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9);
                leg->AddEntry(dataHist, ("Data (" + std::to_string(static_cast<int>(dataHist->GetEntries())) + " entries)").c_str(), "l");
                leg->AddEntry(mcHist, ("MC (" + std::to_string(static_cast<int>(mcHist->GetEntries())) + " entries)").c_str(), "l");
                leg->Draw();

                // Save the canvas to a file
                std::string outputFileName = outputDir + histName + ".png";
                c->SaveAs(outputFileName.c_str());

                // Clean up
                delete dataHist;
                delete mcHist;
                delete c;
                delete leg;

                // Restart the TTreeReaders for the next branch
                dataReader.Restart();
                mcReader.Restart();
            }
        }
        // Increment the currentFits to process the next set of kinematic variables
        currentFits++;
    }
}

void createCorrelationPlots() {
    const std::string outputDir = "output/correlation_plots/";
    const std::vector<std::string> branchesToSkip = {"helicity", "beam_pol", "target_pol", "DepA", "DepB", "DepC", "DepV", "DepW", "evnum"};

    // Retrieve the list of branches
    TObjArray* branches = dataReader.GetTree()->GetListOfBranches();
    if (!branches) {
        std::cerr << "Error: Unable to retrieve branch list from data TTree." << std::endl;
        return;
    }

    std::vector<std::string> branchNames;
    for (Int_t i = 0; i < branches->GetEntries(); ++i) {
        std::string name = branches->At(i)->GetName();
        if (std::find(branchesToSkip.begin(), branchesToSkip.end(), name) == branchesToSkip.end()) {
            branchNames.push_back(name);
        }
    }

    KinematicCuts kinematicCuts(dataReader);

    for (size_t i = 0; i < branchNames.size(); ++i) {
        for (size_t j = i + 1; j < branchNames.size(); ++j) {
            std::string branchX = branchNames[i];
            std::string branchY = branchNames[j];

            // Use the double version for runnum
            if (branchX == "runnum") branchX = "runnum_double";
            if (branchY == "runnum") branchY = "runnum_double";

            HistConfig configX = histConfigs.count(branchX) ? histConfigs[branchX] : HistConfig{100, 0, 1};
            HistConfig configY = histConfigs.count(branchY) ? histConfigs[branchY] : HistConfig{100, 0, 1};

            std::string histName = branchX + "_vs_" + branchY;
            TH2D* hist = new TH2D(histName.c_str(), "", configX.nBins, configX.xMin, configX.xMax, configY.nBins, configY.xMin, configY.xMax);

            hist->GetXaxis()->SetTitle(formatLabelName(branchX).c_str());
            hist->GetYaxis()->SetTitle(formatLabelName(branchY).c_str());
            hist->GetXaxis()->CenterTitle();
            hist->GetYaxis()->CenterTitle();
            hist->GetXaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->SetTitleSize(0.05);
            hist->GetXaxis()->SetTitleOffset(1.3);
            hist->GetYaxis()->SetTitleOffset(1.6);

            dataReader.Restart();
            while (dataReader.Next()) {
                if (kinematicCuts.applyCuts(0, false)) {
                    double xValue = branchX == "runnum" ? static_cast<double>(*runnumVal) : *TTreeReaderValue<double>(dataReader, branchX.c_str());
                    double yValue = branchY == "runnum" ? static_cast<double>(*runnumVal) : *TTreeReaderValue<double>(dataReader, branchY.c_str());
                    hist->Fill(xValue, yValue);
                }
            }

            TCanvas* c = new TCanvas(histName.c_str(), histName.c_str(), 800, 600);
            c->SetLeftMargin(0.15);
            c->SetBottomMargin(0.15);
            c->SetRightMargin(0.15);
            hist->Draw("COLZ");
            c->SaveAs((outputDir + histName + ".png").c_str());
            delete hist;
            delete c;
        }
    }
}



int main(int argc, char *argv[]) {
  // Start the timer
  auto start_time = std::chrono::high_resolution_clock::now();

  // initialize ROOT application for graphics
  TApplication theApp("App", nullptr, nullptr);
  // Set ROOT to batch mode
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // Check for correct number of command line arguments
  if (argc < 3) {
      std::cout << "Usage: " << argv[0];
      std::cout << " <data_root_file> <mc_root_file>" << std::endl;
      return 1;
  }

  // Load data and mc root files
  TFile* data_file = new TFile(argv[1], "READ");
  TFile* mc_file = new TFile(argv[2], "READ");
  if (!data_file->IsOpen() || !mc_file->IsOpen()) {
    cout << "Error opening ROOT files (is the location correct?). Exiting." << endl;
    return 2;
  } else {
    cout << "-- ROOT files opened successfully." << endl;
  }
  
  TTree* data = (TTree*)data_file->Get("PhysicsEvents");
  TTree* mc = (TTree*)mc_file->Get("PhysicsEvents");

  if (!data || !mc) {
    cout << "-- Error getting trees from ROOT files." << endl;
    return 3;
  } else {
    cout << "-- Trees successfully extracted from ROOT files." << endl << endl;
  }

  dataReader.SetTree(data);  // Initialize the global variable
  mcReader.SetTree(mc);  // Initialize the global variable

  // Generate output file names based on the input data file name and current time
  std::string dataRootFileName = argv[1];
  std::string baseName = dataRootFileName.substr(dataRootFileName.find_last_of("/\\") + 1);
  baseName = baseName.substr(0, baseName.find_last_of(".")); // Remove the file extension

  // Check if the file extension is .root
  if (dataRootFileName.substr(dataRootFileName.find_last_of(".") + 1) != "root") {
      std::cerr << "The input file must be a .root file." << std::endl;
      return 1;
  }

  // Get the current time in UTC
  auto time_now = std::chrono::system_clock::now();
  std::time_t time_now_t = std::chrono::system_clock::to_time_t(time_now);
  // Convert it to UTC (GMT)
  std::tm gmt = *std::gmtime(&time_now_t);
  // Apply the Eastern Time Zone offset
  // EST: UTC-5, EDT: UTC-4. Assuming EST here
  gmt.tm_hour -= 5;
  // Correct for out-of-range hours
  std::mktime(&gmt);
  char buffer[80];
  std::strftime(buffer, sizeof(buffer), "%m_%d_%H%M%S", &gmt);

  std::string timeStamp(buffer);
  std::string output_file = "output/results/asymmetries_" + baseName + 
    "_timeStamp_" + timeStamp + ".txt";
  std::string kinematic_file = "output/results/kinematics_" + baseName + 
    "_timeStamp_" + timeStamp + ".txt";

  // Clear the contents of the output files
  std::ofstream ofs(output_file, std::ios::trunc);
  ofs.close();
  std::ofstream ofs2(kinematic_file, std::ios::trunc);
  ofs2.close();

  // load bins from external csv file
  load_bins_from_csv("imports/bins_single_hadron.csv");
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

  // load run info from external csv file
  load_run_info_from_csv("imports/clas12_run_info.csv");
  cout<< endl << endl <<"-- Loaded information from run_info_rgc.csv" << endl;
  cout << "Found " << run_info_list.size() << " runs." << endl;
  // cpp = 0; // total accumulated charge of positive beam - positive target
  // cpm = 0; // total accumulated charge of positive beam - negative target
  // cmp = 0; // total accumulated charge of negative beam - positive target
  // cmm = 0; // total accumulated charge of negative beam - negative target
  // total_charge_carbon = 0; // total accumulated charge of carbon target
  // charge_acuumulation determines the total charges from the runs in supplied dataFile
  // by comparing to master list of CLAS12 runs 
  charge_accumulation(dataReader, run_info_list);
  cout << "Total pos-pos (beam-target) charge: " << cpp << " (nC). ";
  cout << "Total pos-neg charge: " << cpm << " (nC). ";
  cout << "Total neg-pos charge: " << cmp << " (nC). ";
  cout << "Total neg-neg charge: " << cmm << " (nC). ";
  cout << "Total unpolarized (carbon) charge: " << total_charge_carbon << " (nC)."<< endl << endl;

  createIntegratedKinematicPlots();
  // createIntegratedKinematicPlotsForBinsAndFits();
  createCorrelationPlots();
  currentFits=0;
  dataReader.Restart(); mcReader.Restart();

  for (size_t i = 0; i < allBins.size(); ++i) {
    cout << "-- Beginning kinematic fits." << endl;
    for (int asymmetry = 0; asymmetry < 1; ++asymmetry){
      switch (asymmetry) {
        case 0: cout << "    Beginning chi2 BSA." << endl; break;
        case 1: cout << "    Beginning chi2 TSA." << endl; break;
        case 2: cout << "    Beginning chi2 DSA." << endl; break;
      }
      performChi2Fits(output_file.c_str(), kinematic_file.c_str(), binNames[i], asymmetry);
    }
    cout << endl << "     Completed " << binNames[i] << " chi2 fits." << endl;
    // performMLMFits(output_file.c_str(), kinematic_file.c_str(), binNames[i]);
    cout << endl << "     Completed " << binNames[i] << " MLM fits." << endl;
    cout << endl << endl;
    currentFits++;
  }

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
  return 0;
}


// git pull; 
// g++ -o BSA_rgc_fits_root_tree BSA_rgc_fits_root_tree.cpp 
// -L/site/12gev_phys/2.4/Linux_CentOS7.9.2009-gcc9.2.0/root/6.20.04/lib 
// `root-config --cflags --libs` -lMinuit;
// ./BSA_rgc_fits_root_tree /work/clas12/thayward/CLAS12_SIDIS/RGC/p/rgc_8.7.0_epX_Mx-1.4.root 
// /work/clas12/thayward/CLAS12_SIDIS/RGC/p/rgc_nh3_mc.root /u/home/thayward/test_asymmetries.txt
// /u/home/thayward/test_kinematics.txt