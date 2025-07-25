// tbhayward libraries
#include "common_vars.h"  // Include the common header
#include "load_bins_from_csv.h"
#include "load_run_info_from_csv.h"
#include "dilution_factor.h"
#include "asymmetry_fits.h"
#include "BaseKinematicCuts.h"
#include "KinematicCuts.h"
#include "InclusiveKinematicCuts.h"
#include "SingleHadronKinematicCuts.h"
#include "B2BDihadronKinematicCuts.h"
#include "DihadronKinematicCuts.h"
#include "dvcsKinematicCuts.h"
#include "eppi0KinematicCuts.h"
#include "formatLabelName.h"
#include "readChi2Fits.h"
#include "histConfigs.h"
#include "charge_accumulation.h"
#include "plot_data.h"
#include "modifyTree.h"
#include "fitting_process.h" // Include your header file
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
#include <random>
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
#include <TBranch.h>
#include <cstdio>
#include <cmath>

// Using namespace declaration
using namespace std;

/******************** INCLUSIVE DIS CASE ********************/

std::tuple<int, int, int, int, double, double, double> getInclusiveCounts(int binIndex) {

  // Determine the variable range for the specified bin
  double varMin = allBins[currentFits][binIndex];
  double varMax = allBins[currentFits][binIndex + 1];

  // Create positive and negative helicity histograms
  int npp = 0; 
  int npm = 0;
  int nmp = 0;
  int nmm = 0;

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
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  // Counter to limit the number of processed entries
  while (dataReader.Next()) {

    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // bool passedKinematicCuts = true;
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
      sumVariable += *currentVariable;

      if (*helicity > 0 && *target_pol < 0) { npm++; } 
      else if (*helicity < 0 && *target_pol > 0) {  nmp++; }

      if (*helicity > 0 && (*target_pol > 0 || *runnum < 11571) ) { npp++; } 
      else if (*helicity < 0 && (*target_pol < 0 || *runnum < 11571) ) {  nmm++; } 
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
  double meanPol = sumPol / numEvents; // mean beam polarization for data 
  double Ptp = numEventsPosTarget > 0 ? sumTargetPosPol / numEventsPosTarget : 1;
  double Ptm = numEventsNegTarget > 0 ? -sumTargetNegPol / numEventsNegTarget : 1;
  // the negative sign here is correct; RGC lists the polarizations with signs to tell which is 
  // which but the polarization really should just be "percent of polarized nucleii"
  // the sign gives the helicity

  return std::make_tuple(npp, npm, nmp, nmm, meanPol, Ptp, Ptm);
}

void calculate_inclusive(const char* output_file, const char* kinematic_file,
  const std::string& prefix, int asymmetry_index) {

  // Initialize string streams to store the results for each bin
  std::ostringstream chi2FitsAStream; 
  chi2FitsAStream << std::fixed << std::setprecision(9);

  // Initialize string streams to store the mean variables for each bin
  std::ostringstream meanVariablesStream;
  meanVariablesStream << "\\begin{table}[h]" << endl;
  meanVariablesStream << "\\centering" << endl;
  meanVariablesStream << "\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline" << endl;
  meanVariablesStream << "Bin & $<Q^2>$ & $<W>$ & $<x_B>$ & $<y>$ & $<t>$ &";
  meanVariablesStream << "$<t_{\\text{min}}>$\\\\ \\hline" << endl; 

  // and create string stream prefix depending on current asymmetry we're fitting
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      chi2FitsAStream << prefix << "chi2FitsALUoffset = {";
      break;
    case 1: // target-spin asymmetry
      chi2FitsAStream << prefix << "chi2FitsAULoffset = {";
      break;
    case 2: // double-spin asymmetry
      chi2FitsAStream << prefix << "chi2FitsALL = {";
      break;
    default:
      cout << "Invalid asymmetry_index! Using default function form of BSA." << endl;
  }

  // Determine the number of bins
  size_t numBins = allBins[currentFits].size() - 1;

  // Loop over each bin
  for (size_t i = 0; i < numBins; ++i) {
    cout << "Beginning chi2 fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    char histName[32];
    snprintf(histName, sizeof(histName), "hist_%zu", i);

    // Get counts for the current bin
    auto [npp, npm, nmp, nmm, meanPol, Ptp, Ptm] = getInclusiveCounts(i);
    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double numEvents = 0;
    // Variables to calculate the mean depolarization factor
    double sumDepA = 0; 
    double sumDepC = 0; double sumDepV = 0; double sumDepW = 0;

    // Variables to calculate the mean kinematics in each bin
    double sumQ2 = 0; double sumW = 0; double sumx = 0; double sumy = 0;
    double sumt = 0; double sumtmin = 0;

    // Declare reader locations
    TTreeReaderValue<int> runnum(dataReader, "runnum");
    TTreeReaderValue<int> evnum(dataReader, "evnum");
    TTreeReaderValue<double> Q2(dataReader, "Q2");
    TTreeReaderValue<double> W(dataReader, "W");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> y(dataReader, "y");
    TTreeReaderValue<double> t(dataReader, "t");
    TTreeReaderValue<double> tmin(dataReader, "tmin");
    TTreeReaderValue<double> DepA(dataReader, "DepA");
    TTreeReaderValue<double> DepC(dataReader, "DepC");
    TTreeReaderValue<double> DepV(dataReader, "DepV");
    TTreeReaderValue<double> DepW(dataReader, "DepW");
    TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

    // Determine the variable range for the specified bin
    double varMin = allBins[currentFits][i];
    double varMax = allBins[currentFits][i + 1];
    int counter = 0;
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
      // Check if the currentVariable is within the desired range
      if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
        // sum the kinematic variable values
        sumVariable += *currentVariable;
        sumQ2 += *Q2;
        sumW += *W;
        sumx += *x;
        sumy += *y;
        sumt += *t;
        sumtmin += *tmin;

        // sum the depolarization values
        sumDepA += *DepA;
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
    double meanDepC = numEvents > 0 ? sumDepC / numEvents : 0.0;
    double meanDepV = numEvents > 0 ? sumDepV / numEvents : 0.0;
    double meanDepW = numEvents > 0 ? sumDepW / numEvents : 0.0;

    // Calculate the mean values for the kinematic variables
    double meanQ2 = numEvents > 0 ? sumQ2 / numEvents : 0.0;
    double meanW = numEvents > 0 ? sumW / numEvents : 0.0;
    double meanx = numEvents > 0 ? sumx / numEvents : 0.0;
    double meany = numEvents > 0 ? sumy / numEvents : 0.0;
    double meant = numEvents > 0 ? sumt / numEvents : 0.0;
    double meantmin = numEvents > 0 ? sumtmin / numEvents : 0.0;


    switch (asymmetry_index) {
      case 0: {// beam-spin asymmetry
        // Get the fitted parameters and their errors
        double ALU_offset = asymmetry_value_calculation(meanVariable, 
          prefix, npp, npm, nmp, nmm, meanPol, Ptp, Ptm, asymmetry_index);
        double ALU_offset_error = asymmetry_error_calculation(meanVariable,
          prefix, npp, npm, nmp, nmm, meanPol, Ptp, Ptm, asymmetry_index);
        ALU_offset = (meanDepA/meanDepW)*ALU_offset;
        ALU_offset_error = (meanDepA/meanDepW)*ALU_offset_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALU_offset << ", " << ALU_offset_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; 
        }
        break;
      }
      case 1: {// target-spin asymmetry
        // Get the fitted parameters and their errors
        double AUL_offset = asymmetry_value_calculation(meanVariable, 
          prefix, npp, npm, nmp, nmm, meanPol, Ptp, Ptm, asymmetry_index);
        double AUL_offset_error = asymmetry_error_calculation(meanVariable,
          prefix, npp, npm, nmp, nmm, meanPol, Ptp, Ptm, asymmetry_index);
        AUL_offset = (meanDepA/meanDepV)*AUL_offset;
        AUL_offset_error = (meanDepA/meanDepV)*AUL_offset_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< AUL_offset << ", " << AUL_offset_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; 
        }
        break;
      }
      case 2: {// double-spin asymmetry
        // Get the fitted parameters and their errors
        double ALL = asymmetry_value_calculation(meanVariable, prefix, 
          npp, npm, nmp, nmm, meanPol, Ptp, Ptm, asymmetry_index);
        double ALL_error = asymmetry_error_calculation(meanVariable, prefix, 
          npp, npm, nmp, nmm, meanPol, Ptp, Ptm, asymmetry_index);
        ALL = -(meanDepA/meanDepC)*ALL;
        ALL_error = (meanDepA/meanDepC)*ALL_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALL << ", " << ALL_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; 
        }
        break;
      }
    }

    // outputs of mean kinematic variables
    meanVariablesStream << std::fixed << std::setprecision(3); // Set precision to 3 digits 
    meanVariablesStream << (i+1) << "~&~" << meanQ2 << "~&~" << meanW << "~&~" << meanx << "~&~";
    meanVariablesStream << meany << "~&~";
    meanVariablesStream << meant << "~&~" << meantmin; 
    meanVariablesStream << std::string(" \\\\ \\hline ");
  }

  chi2FitsAStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << chi2FitsAStream.str() << std::endl;

  outputFile.close();

  meanVariablesStream << "\\end{tabular}\n";
  meanVariablesStream << "\\caption{The mean kinematic variables in each of the bins ";
  meanVariablesStream << "for the extracted $" << prefix << "$ asymmetries.";
  meanVariablesStream << " Values given in GeV or GeV$^2$ where appropriate.}\n";
  meanVariablesStream << "\\label{table:kinematics_" << prefix << "}\n";
  meanVariablesStream << "\\end{table}\n";
  meanVariablesStream << endl << endl << endl;
  if (asymmetry_index == 0) {
    std::ofstream kinematicFile(kinematic_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicFile << meanVariablesStream.str() << std::endl; 
    kinematicFile.close();
  }
}

/******************** SINGLE HADRON CASE ********************/

// Negative log-likelihood function
void negLogLikelihood_single_hadron(Int_t &npar, Double_t *gin, Double_t &f, 
  Double_t *par, Int_t iflag) {
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
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> phi(dataReader, "phi");
  TTreeReaderValue<double> Q2(dataReader, "Q2");
  TTreeReaderValue<double> x(dataReader, "x");
  TTreeReaderValue<double> z(dataReader, "z");
  TTreeReaderValue<double> pT(dataReader, "pT");
  TTreeReaderValue<double> DepA(dataReader, "DepA");
  TTreeReaderValue<double> DepB(dataReader, "DepB");
  TTreeReaderValue<double> DepC(dataReader, "DepC");
  TTreeReaderValue<double> DepV(dataReader, "DepV");
  TTreeReaderValue<double> DepW(dataReader, "DepW");
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  // Initial definitions (move outside the loop)
  double dilution_factor = dilutionFactors[currentBin].first;
  double sigmaDf = dilutionFactors[currentBin].second;
  double sigmaPb = 0.015;
  double sigmaPtp = 0.025;
  double sigmaPtm = 0.025;

  // Random number generation setup (outside the loop)
  std::random_device rd;
  std::mt19937 gen(rd());

  // Normal distributions (outside the loop)
  std::normal_distribution<> distDf(0.0, sigmaDf);
  std::normal_distribution<> distPb(0.0, sigmaPb);
  std::normal_distribution<> distStandard(0.0, 1.0); // Standard normal distribution

  while (dataReader.Next()) {
    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= allBins[currentFits][currentBin] && 
          *currentVariable < allBins[currentFits][currentBin + 1] && passedKinematicCuts) {

      // Increment the event count
      N += 1;

      // // Get per-event values
      double Df = dilution_factor + distDf(gen);
      double Pb = *beam_pol;
      double Pt = std::abs(*target_pol);
      // Adjust Pb with its uncertainty
      Pb += distPb(gen);
      // Select sigma for Pt based on the sign of *target_pol
      double sigmaPt = (*target_pol >= 0) ? sigmaPtp : sigmaPtm;
      // Adjust Pt with its uncertainty
      Pt += sigmaPt * distStandard(gen);
      // Restore the sign of Pt
      double signPt = (*target_pol >= 0) ? 1.0 : -1.0;
      Pt = signPt * Pt;

      // Proceed with your calculations
      if (*helicity > 0 && *target_pol >= 0) { 
        sum_PP += log(1 
          + (*DepV / *DepA)*AUU_cosphi*cos(*phi) + (*DepB / *DepA)*AUU_cos2phi*cos(2 * *phi) // UU 
          + Pb*((*DepW / *DepA)*ALU_sinphi*sin(*phi)) // BSA
          + Df*Pt*((*DepV / *DepA)*AUL_sinphi*sin(*phi)+ // TSA
            (*DepB / *DepA)*AUL_sin2phi*sin(2 * *phi))//TSA
          + Df*Pb*Pt*((*DepC / *DepA)*ALL + (*DepW / *DepA)*ALL_cosphi*cos(*phi)) ); // DSA
      } else if (*helicity > 0 && *target_pol < 0) { 
        sum_PM += log(1 
          + (*DepV / *DepA)*AUU_cosphi*cos(*phi) + (*DepB / *DepA)*AUU_cos2phi*cos(2 * *phi) // UU
          + Pb*((*DepW / *DepA)*ALU_sinphi*sin(*phi)) // BSA
          - Df*Pt*((*DepV / *DepA)*AUL_sinphi*sin(*phi)+ // TSA
            (*DepB / *DepA)*AUL_sin2phi*sin(2 * *phi)) // TSA
          - Df*Pb*Pt*((*DepC / *DepA)*ALL + (*DepW / *DepA)*ALL_cosphi*cos(*phi)) ); // DSA
      } else if (*helicity < 0 && *target_pol >= 0) { 
        sum_MP += log(1 
          + (*DepV / *DepA)*AUU_cosphi*cos(*phi) + (*DepB / *DepA)*AUU_cos2phi*cos(2 * *phi) // UU 
          - Pb*((*DepW / *DepA)*ALU_sinphi*sin(*phi)) // BSA
          + Df*Pt*((*DepV / *DepA)*AUL_sinphi*sin(*phi)+ // TSA
            (*DepB / *DepA)*AUL_sin2phi*sin(2 * *phi))//TSA
          - Df*Pb*Pt*((*DepC / *DepA)*ALL + (*DepW / *DepA)*ALL_cosphi*cos(*phi)) ); // DSA
      } else if (*helicity < 0 && *target_pol < 0) { 
        sum_MM += log(1 
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
  TTreeReaderValue<double> mc_DepA(mcReader, "DepA");
  TTreeReaderValue<double> mc_DepB(mcReader, "DepB");
  TTreeReaderValue<double> mc_DepC(mcReader, "DepC");
  TTreeReaderValue<double> mc_DepV(mcReader, "DepV");
  TTreeReaderValue<double> mc_DepW(mcReader, "DepW");
  TTreeReaderValue<double> mc_currentVariable(mcReader, propertyNames[currentFits].c_str());

  while (mcReader.Next()) {
    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = mckinematicCuts->applyCuts(currentFits, true);
    // Check if the currentVariable is within the desired range
    if (*mc_currentVariable >= allBins[currentFits][currentBin] && 
          *mc_currentVariable < allBins[currentFits][currentBin + 1] && passedKinematicCuts) {
      NUU += 1 + (*mc_DepV / *mc_DepA)*AUU_cosphi*cos(*mc_phi) +
        (*mc_DepB / *mc_DepA)*AUU_cos2phi*cos(2 * *mc_phi); // UU
    }
  }
  mcReader.Restart();  // Reset the TTreeReader at the end of the function

  // Determine min positive or negative beam helicity accumulated charge to scale down higher one
  double minBeamCharge = std::min({(cpp+cpm),(cmp+cmm)}); 
  // Determine min positive or negative target helicity accumulated charge to scale down higher one
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

void performMLMFits_single_hadron(const char* output_file, const char* kinematic_file,
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
  minuit.SetFCN(negLogLikelihood_single_hadron);

  // Declare string streams for storing the MLM fit results
  std::ostringstream mlmFitsAStream; std::ostringstream mlmFitsBStream; 
  std::ostringstream mlmFitsCStream; std::ostringstream mlmFitsDStream; 
  std::ostringstream mlmFitsEStream; std::ostringstream mlmFitsFStream;
  std::ostringstream mlmFitsGStream; 

  mlmFitsAStream << std::fixed << std::setprecision(9);
  mlmFitsBStream << std::fixed << std::setprecision(9);
  mlmFitsCStream << std::fixed << std::setprecision(9);
  mlmFitsDStream << std::fixed << std::setprecision(9);
  mlmFitsEStream << std::fixed << std::setprecision(9);
  mlmFitsFStream << std::fixed << std::setprecision(9);
  mlmFitsGStream << std::fixed << std::setprecision(9);

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
    cout << endl << "Beginning MLM fit for " << binNames[currentFits]+1
      << " bin " << i << ". ";
    currentBin = i;

    // // Read the chi2 fits into a map
    // std::map<std::string, std::vector<std::vector<double>>> chi2Fits = 
    //     readChi2Fits(std::string(output_file));

    // Construct the key based on the prefix and the fit name
    // For now, let's assume fitName is a string that contains the fit name like "ALUsinphi"
    // std::string fitName = "ALUsinphi";// Replace this with the logic to determine the fit name
    // std::string key = std::string(prefix) + "chi2Fits" + fitName; 

    // std::vector<double> chi2Result = chi2Fits[key][currentFits];
    // Define the parameters with initial values and limits
    minuit.DefineParameter(0, "ALU_sinphi", -0.02, 0.01, -1, 1);
    minuit.DefineParameter(1, "AUL_sinphi", -0.02, 0.01, -1, 1);
    minuit.DefineParameter(2, "AUL_sin2phi", -0.02, 0.01, -1, 1);
    minuit.DefineParameter(3, "ALL", 0.30, 0.01, -1, 1);
    minuit.DefineParameter(4, "ALL_cosphi", 0.01, 0.01, -1, 1);
    minuit.DefineParameter(5, "AUU_cosphi", 0.00, 0.00, -1, 1);
    minuit.DefineParameter(6, "AUU_cos2phi", 0.00, 0.00, -1, 1);

    // After defining parameters
    minuit.Migrad(); cout << endl; // First attempt to find the minimum

    // If you decide to use MINImize, replace Migrad with the following lines:
    // arglist[0] = 500; // Max calls
    // arglist[1] = 1.;  // Tolerance
    // minuit.mnexcm("MINImize", arglist, 2, ierflg);


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
    TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
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

void plotHistogramAndFit_single_hadron(TH1D* histogram, TF1* fitFunction, int binIndex, 
  int asymmetryIndex, const std::string& prefix) {
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
  string filename = "output/individual_chi2_fits/" + prefix + "_" + 
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
  if (canvas->GetListOfPrimitives()->GetSize() > 0) {
      // There's something in the canvas, save it
      canvas->SaveAs(filename.c_str());
  } else {
      std::cout << "Canvas is empty, not saving to file." << std::endl;
  }

  // Clean up
  delete canvas;
  delete graph;
}

TH1D* createHistogramForBin_single_hadron(const char* histName, int binIndex, 
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
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> x(dataReader, "x");
  TTreeReaderValue<double> Q2(dataReader, "Q2");
  TTreeReaderValue<double> z(dataReader, "z");
  TTreeReaderValue<double> pT(dataReader, "pT");
  TTreeReaderValue<double> DepA(dataReader, "DepA");
  TTreeReaderValue<double> DepB(dataReader, "DepB");
  TTreeReaderValue<double> phi(dataReader, "phi");
  // TTreeReaderValue<double> phi(dataReader, "phi2");
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());
  // TTreeReaderValue<int> currentVariable(dataReader, propertyNames[currentFits].c_str());

  while (dataReader.Next()) {
    
    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // bool passedKinematicCuts = true;
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
      sumVariable += *currentVariable;
      if (*helicity > 0 && *target_pol < 0) { histPosNeg->Fill(*phi); } 
      else if (*helicity < 0 && *target_pol > 0) {  histNegPos->Fill(*phi); }

      if (*helicity > 0 && (*target_pol >= 0) ) { histPosPos->Fill(*phi); } 
      else if (*helicity < 0 && (*target_pol <= 0) ) {  histNegNeg->Fill(*phi); } 
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
    double asymmetry = asymmetry_value_calculation(meanVariable, prefix, 
      Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, asymmetry_index);
    double error = asymmetry_error_calculation(meanVariable, prefix, 
      Npp*cpp, Npm*cpm, Nmp*cmp, Nmm*cmm, meanPol, Ptp, Ptm, asymmetry_index);

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

void performChi2Fits_single_hadron(const char* output_file, const char* kinematic_file,
  const char* kinematicPlot_file, const std::string& prefix, int asymmetry_index) {

  // Initialize string streams to store the results for each bin
  std::ostringstream chi2FitsAStream, chi2FitsBStream, chi2FitsCStream;
  chi2FitsAStream << std::fixed << std::setprecision(9);
  chi2FitsBStream << std::fixed << std::setprecision(9);
  chi2FitsCStream << std::fixed << std::setprecision(9);

  // Initialize string stream to store the kinematics in each bin for use in LaTeX 
  std::ostringstream meanVariablesStream;
  meanVariablesStream << "\\begin{table}[h]" << endl;
  meanVariablesStream << "\\centering" << endl;
  meanVariablesStream << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|} \\hline" << endl;
  meanVariablesStream << "Bin & $<Q^2>$ & $<W>$ ";
  meanVariablesStream << "& $<x_B>$ & $<y>$ & $<z>$ & $<\\xi>$ & $<P_T>$ ";
  meanVariablesStream << "& $<x_F>$ & $<t>$ & ";
  meanVariablesStream << "$<t_{\\text{min}}>$\\\\ \\hline" << endl; 

  // Initalize string stream to store the kinematics in each bin for use in plotting 
  std::ostringstream meanVariablesPlotStream;
  meanVariablesPlotStream << prefix << "Kinematics = {";

  // Create a new TF1 object called fitFunction representing the function to fit
  // and create string stream prefix depending on current asymmetry we're fitting
  TF1* fitFunction;
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      fitFunction = new TF1("fitFunction", BSA_single_hadron, 0, 2*TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALUoffset = {";
      chi2FitsBStream << prefix << "chi2FitsALUsinphi = {";
      break;
    case 1: // target-spin asymmetry
      fitFunction = new TF1("fitFunction", TSA_single_hadron, 0, 2*TMath::Pi(), 3);
      chi2FitsAStream << prefix << "chi2FitsAULoffset = {";
      chi2FitsBStream << prefix << "chi2FitsAULsinphi = {";
      chi2FitsCStream << prefix << "chi2FitsAULsin2phi = {";
      break;
    case 2: // double-spin asymmetry
      fitFunction = new TF1("fitFunction", DSA_single_hadron, 0, 2*TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALL = {";
      chi2FitsBStream << prefix << "chi2FitsALLcosphi = {";
      break;
    default:
      cout << "Invalid asymmetry_index! Using default function form of BSA." << endl;
      fitFunction = new TF1("fitFunction", BSA_single_hadron, 0, 2*TMath::Pi(), 2);
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
    TH1D* hist = createHistogramForBin_single_hadron(histName, i, prefix, asymmetry_index);
    // Fit the histogram using the fitFunction and get the fit result
    hist->Fit(fitFunction, "QS");
    plotHistogramAndFit_single_hadron(hist, fitFunction, i, asymmetry_index, prefix);

    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double numEvents = 0;
    // Variables to calculate the mean depolarization factor
    double sumDepA = 0; double sumDepB = 0; 
    double sumDepC = 0; double sumDepV = 0; double sumDepW = 0;

    // Variables to calculate the mean kinematics in each bin
    double sumQ2 = 0; double sumW = 0; double sumx = 0; double sumy = 0;
    double sumz = 0; double sumxi = 0; double sumpT = 0; double sumxF = 0;
    double sumt = 0; double sumtmin = 0;

    // Declare reader locations
    TTreeReaderValue<int> runnum(dataReader, "runnum");
    TTreeReaderValue<int> evnum(dataReader, "evnum");
    TTreeReaderValue<double> Q2(dataReader, "Q2");
    TTreeReaderValue<double> W(dataReader, "W");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> y(dataReader, "y");
    TTreeReaderValue<double> z(dataReader, "z");
    TTreeReaderValue<double> xi(dataReader, "xi");
    TTreeReaderValue<double> pT(dataReader, "pT");
    TTreeReaderValue<double> xF(dataReader, "xF");
    TTreeReaderValue<double> t(dataReader, "t");
    TTreeReaderValue<double> tmin(dataReader, "tmin");
    TTreeReaderValue<double> DepA(dataReader, "DepA");
    TTreeReaderValue<double> DepB(dataReader, "DepB");
    TTreeReaderValue<double> DepC(dataReader, "DepC");
    TTreeReaderValue<double> DepV(dataReader, "DepV");
    TTreeReaderValue<double> DepW(dataReader, "DepW");
    TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());
    // Declare a pointer for currentVariable


    // Determine the variable range for the specified bin
    double varMin = allBins[currentFits][i];
    double varMax = allBins[currentFits][i + 1];
    int counter = 0;
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
      // Check if the currentVariable is within the desired range
      if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
        // sum the kinematic variable values
        sumVariable += *currentVariable;
        sumQ2 += *Q2;
        sumW += *W;
        sumx += *x;
        sumy += *y;
        sumz += *z;
        sumxi += *xi;
        sumpT += *pT;
        sumxF += *xF;
        sumt += *t;
        sumtmin += *tmin;
        // double epsilonNum = *DepB; double epsilonDen = *DepA; 
        // sumtmin += epsilonNum/epsilonDen;

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
    double meanxi = numEvents > 0 ? sumxi / numEvents : 0.0;
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
        ALU_sinphi = (meanDepA/meanDepW)*ALU_sinphi;
        ALU_sinphi_error = (meanDepA/meanDepW)*ALU_sinphi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALU_offset << ", " << ALU_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALU_sinphi << ", " << ALU_sinphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; 
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
        AUL_sinphi = (meanDepA/meanDepV)*AUL_sinphi;
        AUL_sinphi_error = (meanDepA/meanDepV)*AUL_sinphi_error;
        AUL_sin2phi = (meanDepA/meanDepB)*AUL_sin2phi;
        AUL_sin2phi_error = (meanDepA/meanDepB)*AUL_sin2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< AUL_offset << ", " << AUL_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< AUL_sinphi << ", " << AUL_sinphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUL_sin2phi << ", " << AUL_sin2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; chi2FitsCStream << ", ";
        }
        break;
      }
      case 2: {// double-spin asymmetry
        // Get the fitted parameters and their errors
        double ALL = fitFunction->GetParameter(0);
        double ALL_error = fitFunction->GetParError(0);
        double ALL_cosphi = fitFunction->GetParameter(1);
        double ALL_cosphi_error = fitFunction->GetParError(1);
        ALL = (meanDepA/meanDepC)*ALL;
        ALL_error = (meanDepA/meanDepC)*ALL_error;
        ALL_cosphi = (meanDepA/meanDepW)*ALL_cosphi;
        ALL_cosphi_error = (meanDepA/meanDepW)*ALL_cosphi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALL << ", " << ALL_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALL_cosphi << ", " << ALL_cosphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", ";
        }
        break;
      }
    }

    delete hist;

    // outputs of mean kinematic variables for LaTeX
    meanVariablesStream << std::fixed << std::setprecision(3); // Set precision to 3 digits 
    meanVariablesStream << (i+1) << "~&~" << meanQ2 << "~&~" << meanW << "~&~" << meanx << "~&~";
    meanVariablesStream << meany << "~&~" << meanz << "~&~" << meanxi << "~&~";
    meanVariablesStream << meanpT << "~&~" << meanxF << "~&~" << meant << "~&~" << meantmin; 
    meanVariablesStream << std::string(" \\\\ \\hline ");

    // outputs of mean kinematic variables for plotting
    meanVariablesPlotStream << "{" << meanQ2 << ", " << meanW << ", " << meanx << ", ";
    meanVariablesPlotStream << meany << ", " << meanz << ", " << meanxi << ", ";
    meanVariablesPlotStream << meanpT << ", " << meanxF << ", " << meant << "}";
    if (i < numBins - 1) {
        meanVariablesPlotStream << ", "; 
    }
  }

  chi2FitsAStream << "};";  chi2FitsBStream << "};";  chi2FitsCStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << chi2FitsAStream.str() << std::endl;
  outputFile << chi2FitsBStream.str() << std::endl;
  if (asymmetry_index==1) { outputFile << chi2FitsCStream.str() << std::endl; }

  outputFile.close();

  meanVariablesStream << "\\end{tabular}\n";
  meanVariablesStream << "\\caption{The mean kinematic variables in each of the bins ";
  meanVariablesStream << "for the extracted $" << prefix << "$ asymmetries.";
  meanVariablesStream << " Values given in GeV or GeV$^2$ where appropriate.}\n";
  meanVariablesStream << "\\end{table}\n";
  meanVariablesStream << endl << endl << endl;
  if (asymmetry_index == 0) {
    std::ofstream kinematicFile(kinematic_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicFile << meanVariablesStream.str() << std::endl; 
    kinematicFile.close();

    meanVariablesPlotStream << "};";
    std::ofstream kinematicPlot_File(kinematicPlot_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicPlot_File << meanVariablesPlotStream.str() << std::endl;
    kinematicPlot_File.close();
  }
}


/******************** B2B DIHADRON (dSIDIS) CASE ********************/

// Negative log-likelihood function
void negLogLikelihood_b2b_dihadron(Int_t &npar, Double_t *gin, Double_t &f, 
  Double_t *par, Int_t iflag) {
  // npar: number of parameters
  // gin: an array of derivatives (if needed)
  // f: the value of the function
  // par: an array of the parameter values
  // iflag: a flag (see TMinuit documentation for details)

  // Extract parameters from the input parameter array
  // LU
  double ALU_sinphi1 = par[0];
  double ALU_sinphi2 = par[1];
  double ALU_sinDeltaphi = par[2];
  double ALU_sin2Deltaphi = par[3];
  // UL 
  double AUL_sinphi1 = par[4];
  double AUL_sinphi2 = par[5];
  double AUL_sin2phi1 = par[6];
  double AUL_sin2phi2 = par[7];
  double AUL_sinDeltaphi = par[8];
  double AUL_sin2Deltaphi = par[9];
  double AUL_sinSumphi = par[10];
  // LL
  double ALL = par[11];
  double ALL_cosphi1 = par[12]; 
  double ALL_cosphi2 = par[13]; 
  // UU
  double AUU_cosphi1 = par[14];
  double AUU_cosphi2 = par[15]; 
  double AUU_cos2phi1 = par[16];
  double AUU_cos2phi2 = par[17];
  double AUU_cosSumphi = par[18]; 

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
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> Q2(dataReader, "Q2");
  TTreeReaderValue<double> x(dataReader, "x");
  TTreeReaderValue<double> z(dataReader, "z");
  TTreeReaderValue<double> pT(dataReader, "pT");
  TTreeReaderValue<double> phi1(dataReader, "phi1");
  TTreeReaderValue<double> phi2(dataReader, "phi2");
  TTreeReaderValue<double> DepA(dataReader, "DepA");
  TTreeReaderValue<double> DepB(dataReader, "DepB");
  TTreeReaderValue<double> DepC(dataReader, "DepC");
  TTreeReaderValue<double> DepV(dataReader, "DepV");
  TTreeReaderValue<double> DepW(dataReader, "DepW");
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  while (dataReader.Next()) {
    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= allBins[currentFits][currentBin] && 
          *currentVariable < allBins[currentFits][currentBin + 1] && passedKinematicCuts) {

      // Increment the event count
      N += 1;

      // Initial definitions
      double Df = dilutionFactors[currentBin].first;
      double sigmaDf = dilutionFactors[currentBin].second;
      double Pb = *beam_pol;
      double sigmaPb = 0.015;
      double Pt = std::abs(*target_pol);
      double sigmaPtp = 0.025;
      double sigmaPtm = 0.025;

      // Random number generation setup
      std::random_device rd;
      std::mt19937 gen(rd());

      // Normal distributions
      std::normal_distribution<> distDf(0.0, sigmaDf);
      std::normal_distribution<> distPb(0.0, sigmaPb);

      // Select sigma for Pt based on the sign of *target_pol
      double sigmaPt = (*target_pol >= 0) ? sigmaPtp : sigmaPtm;
      std::normal_distribution<> distPt(0.0, sigmaPt);

      // Adjust the values
      Df += distDf(gen);
      Pb += distPb(gen);
      Pt += distPt(gen);

      // Restore the sign of Pt
      double signPt = (*target_pol >= 0) ? 1.0 : -1.0;
      Pt = signPt * Pt;

      if (*helicity > 0 && *target_pol >= 0) { 
        sum_PP = sum_PP + log(1 +
          // UU
          (*DepV / *DepA)*AUU_cosphi1*cos(*phi1) +
          (*DepV / *DepA)*AUU_cosphi2*cos(*phi2) +
          (*DepB / *DepA)*AUU_cos2phi1*cos(2**phi1) +
          (*DepB / *DepA)*AUU_cos2phi2*cos(2**phi2) +
          (*DepB / *DepA)*AUU_cosSumphi*cos(*phi1 + *phi2) 
          // LU
          + Pb*(
            (*DepW / *DepA)*ALU_sinphi1*sin(*phi1) +
            (*DepW / *DepA)*ALU_sinphi2*sin(*phi2) +
            (*DepC / *DepA)*ALU_sinDeltaphi*sin(*phi1 - *phi2) +
            (*DepC / *DepA)*ALU_sin2Deltaphi*sin(2**phi1 - 2**phi2)
          )
          // UL
          + Df*Pt*(
            (*DepV / *DepA)*AUL_sinphi1*sin(*phi1) +
            (*DepV / *DepA)*AUL_sinphi2*sin(*phi2) +
            (*DepB / *DepA)*AUL_sin2phi1*sin(2**phi1) +
            (*DepB / *DepA)*AUL_sin2phi2*sin(2**phi2) +
            AUL_sinDeltaphi*sin(*phi1 - *phi2) +  // yes, there is no depolarization factor
            AUL_sin2Deltaphi*sin(2**phi1 - 2**phi2) +  // yes, there is no depolarization factor
            (*DepB / *DepA)*AUL_sinSumphi*sin(*phi1 - *phi2) 
          )
          // LL
          + Df*Pb*Pt*(
            (*DepC / *DepA)*ALL + 
            (*DepW / *DepA)*ALL_cosphi1*cos(*phi1) + 
            (*DepW / *DepA)*ALL_cosphi2*cos(*phi2) 
          )
        );
      }
      if (*helicity < 0 && *target_pol >= 0) { 
        sum_MP = sum_MP + log(1 +
          // UU
          (*DepV / *DepA)*AUU_cosphi1*cos(*phi1) +
          (*DepV / *DepA)*AUU_cosphi2*cos(*phi2) +
          (*DepB / *DepA)*AUU_cos2phi1*cos(2**phi1) +
          (*DepB / *DepA)*AUU_cos2phi2*cos(2**phi2) +
          (*DepB / *DepA)*AUU_cosSumphi*cos(*phi1 + *phi2) 
          // LU
          - Pb*(
            (*DepW / *DepA)*ALU_sinphi1*sin(*phi1) +
            (*DepW / *DepA)*ALU_sinphi2*sin(*phi2) +
            (*DepC / *DepA)*ALU_sinDeltaphi*sin(*phi1 - *phi2) +
            (*DepC / *DepA)*ALU_sinDeltaphi*sin(2**phi1 - 2**phi2)
          )
          // UL
          + Df*Pt*(
            (*DepV / *DepA)*AUL_sinphi1*sin(*phi1) +
            (*DepV / *DepA)*AUL_sinphi2*sin(*phi2) +
            (*DepB / *DepA)*AUL_sin2phi1*sin(2**phi1) +
            (*DepB / *DepA)*AUL_sin2phi2*sin(2**phi2) +
            AUL_sinDeltaphi*sin(*phi1 - *phi2) +  // yes, there is no depolarization factor
            AUL_sin2Deltaphi*sin(2**phi1 - 2**phi2) +  // yes, there is no depolarization factor
            (*DepB / *DepA)*AUL_sinSumphi*sin(*phi1 - *phi2) 
          )
          // LL
          - Df*Pb*Pt*(
            (*DepC / *DepA)*ALL + 
            (*DepW / *DepA)*ALL_cosphi1*cos(*phi1) + 
            (*DepW / *DepA)*ALL_cosphi2*cos(*phi2) 
          )
        );
      }
      if (*helicity > 0 && *target_pol < 0) { 
        sum_PP = sum_PM + log(1 +
          // UU
          (*DepV / *DepA)*AUU_cosphi1*cos(*phi1) +
          (*DepV / *DepA)*AUU_cosphi2*cos(*phi2) +
          (*DepB / *DepA)*AUU_cos2phi1*cos(2**phi1) +
          (*DepB / *DepA)*AUU_cos2phi2*cos(2**phi2) +
          (*DepB / *DepA)*AUU_cosSumphi*cos(*phi1 + *phi2) 
          // LU
          + Pb*(
            (*DepW / *DepA)*ALU_sinphi1*sin(*phi1) +
            (*DepW / *DepA)*ALU_sinphi2*sin(*phi2) +
            (*DepC / *DepA)*ALU_sinDeltaphi*sin(*phi1 - *phi2) +
            (*DepC / *DepA)*ALU_sinDeltaphi*sin(2**phi1 - 2**phi2)
          )
          // UL
          - Df*Pt*(
            (*DepV / *DepA)*AUL_sinphi1*sin(*phi1) +
            (*DepV / *DepA)*AUL_sinphi2*sin(*phi2) +
            (*DepB / *DepA)*AUL_sin2phi1*sin(2**phi1) +
            (*DepB / *DepA)*AUL_sin2phi2*sin(2**phi2) +
            AUL_sinDeltaphi*sin(*phi1 - *phi2) +  // yes, there is no depolarization factor
            AUL_sin2Deltaphi*sin(2**phi1 - 2**phi2) +  // yes, there is no depolarization factor
            (*DepB / *DepA)*AUL_sinSumphi*sin(*phi1 - *phi2) 
          )
          // LL
          - Df*Pb*Pt*(
            (*DepC / *DepA)*ALL + 
            (*DepW / *DepA)*ALL_cosphi1*cos(*phi1) + 
            (*DepW / *DepA)*ALL_cosphi2*cos(*phi2) 
          )
        );
      }
      if (*helicity < 0 && *target_pol < 0) { 
        sum_MM = sum_MM + log(1 +
          // UU
          (*DepV / *DepA)*AUU_cosphi1*cos(*phi1) +
          (*DepV / *DepA)*AUU_cosphi2*cos(*phi2) +
          (*DepB / *DepA)*AUU_cos2phi1*cos(2**phi1) +
          (*DepB / *DepA)*AUU_cos2phi2*cos(2**phi2) +
          (*DepB / *DepA)*AUU_cosSumphi*cos(*phi1 + *phi2) 
          // LU
          - Pb*(
            (*DepW / *DepA)*ALU_sinphi1*sin(*phi1) +
            (*DepW / *DepA)*ALU_sinphi2*sin(*phi2) +
            (*DepC / *DepA)*ALU_sinDeltaphi*sin(*phi1 - *phi2) +
            (*DepC / *DepA)*ALU_sinDeltaphi*sin(2**phi1 - 2**phi2)
          )
          // UL
          - Df*Pt*(
            (*DepV / *DepA)*AUL_sinphi1*sin(*phi1) +
            (*DepV / *DepA)*AUL_sinphi2*sin(*phi2) +
            (*DepB / *DepA)*AUL_sin2phi1*sin(2**phi1) +
            (*DepB / *DepA)*AUL_sin2phi2*sin(2**phi2) +
            AUL_sinDeltaphi*sin(*phi1 - *phi2) +  // yes, there is no depolarization factor
            AUL_sin2Deltaphi*sin(2**phi1 - 2**phi2) +  // yes, there is no depolarization factor
            (*DepB / *DepA)*AUL_sinSumphi*sin(*phi1 - *phi2) 
          )
          // LL
          + Df*Pb*Pt*(
            (*DepC / *DepA)*ALL + 
            (*DepW / *DepA)*ALL_cosphi1*cos(*phi1) + 
            (*DepW / *DepA)*ALL_cosphi2*cos(*phi2) 
          )
        );
      }
    }
  }
  dataReader.Restart();  // Reset the TTreeReader at the end of the function
  
  TTreeReaderValue<double> mc_phi1(mcReader, "phi1");
  TTreeReaderValue<double> mc_phi2(mcReader, "phi2");
  TTreeReaderValue<double> mc_DepA(mcReader, "DepA");
  TTreeReaderValue<double> mc_DepB(mcReader, "DepB");
  TTreeReaderValue<double> mc_DepC(mcReader, "DepC");
  TTreeReaderValue<double> mc_DepV(mcReader, "DepV");
  TTreeReaderValue<double> mc_DepW(mcReader, "DepW");
  TTreeReaderValue<double> mc_currentVariable(mcReader, propertyNames[currentFits].c_str());

  while (mcReader.Next()) {
    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = mckinematicCuts->applyCuts(currentFits, true);
    // Check if the currentVariable is within the desired range
    if (*mc_currentVariable >= allBins[currentFits][currentBin] && 
          *mc_currentVariable < allBins[currentFits][currentBin + 1] && passedKinematicCuts) {
      NUU+=1+
        (*DepV / *DepA)*AUU_cosphi1*cos(*mc_phi1) +
        (*DepV / *DepA)*AUU_cosphi2*cos(*mc_phi2) +
        (*DepB / *DepA)*AUU_cos2phi1*cos(2**mc_phi1) +
        (*DepB / *DepA)*AUU_cos2phi2*cos(2**mc_phi2) +
        (*DepB / *DepA)*AUU_cosSumphi*cos(*mc_phi1 + *mc_phi2);
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
  cout << sum_PP << " " << sum_PM << " " << sum_MP << " " << sum_MM << endl;
  cout << "On MLM fit " << binNames[currentFits] << " " << currentFits << ", " << nll << endl;
  cout << "ALU_sinphi1 = " << ALU_sinphi1;
  cout << ", ALU_sinphi2 = " << ALU_sinphi2;
  cout << ", ALU_sinDeltaphi = " << ALU_sinDeltaphi;
  cout << ", AUL_sinDeltaphi = " << AUL_sinDeltaphi;
  cout << ", AUL_sinSumphi = " << AUL_sinSumphi;
  cout << ", ALL = " << ALL << "." << endl;
  // Calculate the negative log-likelihood value and store it in the output variable f
  f = nll;
}

void performMLMFits_b2b_dihadron(const char* output_file, const char* kinematic_file,
  const std::string& prefix) {
  // Read the event data from the input file and store it in the global variable gData
  mlmPrefix = prefix;

  // Determine the number of bins
  size_t numBins = allBins[currentFits].size() - 1;

  // Initialize TMinuit
  double arglist[10]; arglist[0] = 1;
  int ierflg = 0;
  TMinuit minuit(14); // parameter numbers
  minuit.SetPrintLevel(-1);
  minuit.SetErrorDef(0.5); // error definition for MLE, 1 for chi2
  // This is due to the fact that −logL = chi2/2. 
  // The default value of ErrorDef=1 corresponds to one standard deviation for chi2 function.
  minuit.SetFCN(negLogLikelihood_b2b_dihadron);

  // Initialize string streams for results 
  std::ostringstream mlmFitsStreams[19]; // For maximum 19 parameters (TSA w/ UU case)
  for (auto& stream : mlmFitsStreams) {
    stream << std::fixed << std::setprecision(9);
  }

  mlmFitsStreams[0] << std::fixed << std::setprecision(9);
  mlmFitsStreams[1] << std::fixed << std::setprecision(9);
  mlmFitsStreams[2] << std::fixed << std::setprecision(9);
  mlmFitsStreams[3] << std::fixed << std::setprecision(9);
  mlmFitsStreams[4] << std::fixed << std::setprecision(9);
  mlmFitsStreams[5] << std::fixed << std::setprecision(9);
  mlmFitsStreams[6] << std::fixed << std::setprecision(9);
  mlmFitsStreams[7] << std::fixed << std::setprecision(9);
  mlmFitsStreams[8] << std::fixed << std::setprecision(9);
  mlmFitsStreams[9] << std::fixed << std::setprecision(9);
  mlmFitsStreams[10] << std::fixed << std::setprecision(9);
  mlmFitsStreams[11] << std::fixed << std::setprecision(9);
  mlmFitsStreams[12] << std::fixed << std::setprecision(9);
  mlmFitsStreams[13] << std::fixed << std::setprecision(9);
  mlmFitsStreams[14] << std::fixed << std::setprecision(9);
  mlmFitsStreams[15] << std::fixed << std::setprecision(9);
  mlmFitsStreams[16] << std::fixed << std::setprecision(9);
  mlmFitsStreams[17] << std::fixed << std::setprecision(9);
  mlmFitsStreams[18] << std::fixed << std::setprecision(9);

  // Initialize the string streams with the output variable names
  mlmFitsStreams[0] << prefix << "MLMFitsALUsinphi1 = {";
  mlmFitsStreams[1] << prefix << "MLMFitsALUsinphi2 = {";
  mlmFitsStreams[2] << prefix << "MLMFitsALUsinDeltaphi = {";
  mlmFitsStreams[3] << prefix << "MLMFitsALUsin2Deltaphi = {";
  mlmFitsStreams[4] << prefix << "MLMFitsAULsinphi1 = {";
  mlmFitsStreams[5] << prefix << "MLMFitsAULsinphi2 = {";
  mlmFitsStreams[6] << prefix << "MLMFitsAULsin2phi1 = {";
  mlmFitsStreams[7] << prefix << "MLMFitsAULsin2phi2 = {";
  mlmFitsStreams[8] << prefix << "MLMFitsAULsinDeltaphi = {";
  mlmFitsStreams[9] << prefix << "MLMFitsAULsin2Deltaphi = {";
  mlmFitsStreams[10] << prefix << "MLMFitsAULsinSumphi = {";
  mlmFitsStreams[11] << prefix << "MLMFitsALL = {";
  mlmFitsStreams[12] << prefix << "MLMFitsALLcosphi1 = {";
  mlmFitsStreams[13] << prefix << "MLMFitsALLcosphi2 = {";
  mlmFitsStreams[14] << prefix << "MLMFitsAUUcosphi1 = {";
  mlmFitsStreams[15] << prefix << "MLMFitsAUUcosphi2 = {";
  mlmFitsStreams[16] << prefix << "MLMFitsAUUcos2phi1 = {";
  mlmFitsStreams[17] << prefix << "MLMFitsAUUcos2phi2 = {";
  mlmFitsStreams[18] << prefix << "MLMFitsAUUcosSumphi = {";

  // Initialize string streams to store the mean variables for each bin and asymmetries
  std::ostringstream ALUStream;
  ALUStream << "\\begin{table}[h]" << std::endl;
  ALUStream << "\\centering" << std::endl;
  ALUStream << "\\begin{tabular}{|c|c|c|c|c|c|} \\hline" << std::endl;
  ALUStream << "Bin & $<" << prefix << ">$ & ";
  ALUStream << "$F_{LU}^{\\sin(\\phi_1)}/F_{UU}$ & ";
  ALUStream << "$F_{LU}^{\\sin(\\phi_2)}/F_{UU}$ & ";
  ALUStream << "$F_{LU}^{\\sin(\\phi_1 - \\phi_2)}/F_{UU}$ & ";
  ALUStream << "$F_{LU}^{\\sin(2\\phi_1 - 2\\phi_2)}/F_{UU}$ ";
  ALUStream << "\\\\ \\hline" << std::endl;

  std::ostringstream AULStream;
  AULStream << "\\begin{table}[h]" << std::endl;
  AULStream << "\\centering" << std::endl;
  AULStream << "\\begin{tabular}{|c|c|c|c|c|c|} \\hline" << std::endl;
  AULStream << "Bin & $<" << prefix << ">$ & ";
  AULStream << "$F_{UL}^{\\sin(\\phi_1)}/F_{UU}$ & ";
  AULStream << "$F_{UL}^{\\sin(\\phi_2)}/F_{UU}$ & ";
  AULStream << "$F_{UL}^{\\sin(2\\phi_1)}/F_{UU}$ & ";
  AULStream << "$F_{UL}^{\\sin(2\\phi_2)}/F_{UU}$ & ";
  AULStream << "$F_{UL}^{\\sin(\\phi_1 - \\phi_2)}/F_{UU}$ & ";
  AULStream << "$F_{UL}^{\\sin(2\\phi_1 - 2\\phi_2)}/F_{UU}$ ";
  AULStream << "$F_{UL}^{\\sin(2\\phi_1 + 2\\phi_2)}/F_{UU}$ ";
  AULStream << "\\\\ \\hline" << std::endl;

  std::ostringstream ALLStream;
  ALLStream << "\\begin{table}[h]" << std::endl;
  ALLStream << "\\centering" << std::endl;
  ALLStream << "\\begin{tabular}{|c|c|c|c|c|c|} \\hline" << std::endl;
  ALLStream << "Bin & $<" << prefix << ">$ & ";
  ALLStream << "$F_{LL}/F_{UU}$ & ";
  ALLStream << "$F_{LL}^{\\cos(\\phi_1)}/F_{UU}$ & ";
  ALLStream << "$F_{LL}^{\\cos(\\phi_2)}/F_{UU}$ & ";
  ALUStream << "\\\\ \\hline" << std::endl;

  std::ostringstream AUUStream;
  AUUStream << "\\begin{table}[h]" << std::endl;
  AUUStream << "\\centering" << std::endl;
  AUUStream << "\\begin{tabular}{|c|c|c|c|c|c|} \\hline" << std::endl;
  AUUStream << "Bin & $<" << prefix << ">$ & ";
  AUUStream << "$F_{UU}^{\\cos(\\phi_1)}/F_{UU}$ & ";
  AUUStream << "$F_{UU}^{\\cos(\\phi_2)}/F_{UU}$ & ";
  AUUStream << "$F_{UU}^{\\cos(2\\phi_1)}/F_{UU}$ & ";
  AUUStream << "$F_{UU}^{\\cos(2\\phi_2)}/F_{UU}$ & ";
  AUUStream << "$F_{UU}^{\\cos(\\phi_1 + \\phi_2)}/F_{UU}$ & ";
  AUUStream << "\\\\ \\hline" << std::endl;

  // Iterate through each bin
  for (size_t i = 0; i < numBins; ++i) {
    cout << endl << "Beginning MLM fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    currentBin = i;

    // Define the parameters with initial values and limits
    minuit.DefineParameter(0, "ALU_sinphi1", 0, 0.001, -1, 1);
    minuit.DefineParameter(1, "ALU_sinphi2", 0, 0.001, -1, 1);
    minuit.DefineParameter(2, "ALU_sinDeltaphi", 0, 0.001, -1, 1);
    minuit.DefineParameter(3, "ALU_sin2Deltaphi", 0, 0.001, -1, 1);
    minuit.DefineParameter(4, "AUL_sinphi1", 0, 0.001, -1, 1);
    minuit.DefineParameter(5, "AUL_sinphi2", 0, 0.001, -1, 1);
    minuit.DefineParameter(6, "AUL_sin2phi1", 0.0, 0.001, -1, 1);
    minuit.DefineParameter(7, "AUL_sin2phi2", 0, 0.001, -1, 1);
    minuit.DefineParameter(8, "AUL_sinDeltaphi", 0, 0.001, -1, 1);
    minuit.DefineParameter(9, "AUL_sin2Deltaphi", 0, 0.001, -1, 1);
    minuit.DefineParameter(10, "AUL_sinSumphi", 0, 0.001, -1, 1);
    minuit.DefineParameter(11, "ALL", 0, 0.001, -1, 1);
    minuit.DefineParameter(12, "ALL_cosphi1", 0, 0.001, -1, 1);
    minuit.DefineParameter(13, "ALL_cosphi2", 0.0, 0.001, -1, 1);
    minuit.DefineParameter(14, "AUU_cosphi1", 0, 0.001, -1, 1);
    minuit.DefineParameter(15, "AUU_cosphi2", 0, 0.001, -1, 1);
    minuit.DefineParameter(16, "AUU_cos2phi1", 0, 0.001, -1, 1);
    minuit.DefineParameter(17, "AUU_cos2phi2", 0, 0.001, -1, 1);
    minuit.DefineParameter(18, "AUU_cosSumphi", 0, 0.001, -1, 1);

    // After defining parameters
    minuit.Migrad(); cout << endl; // First attempt to find the minimum

    // If you decide to use MINImize, replace Migrad with the following lines:
    arglist[0] = 500; // Max calls
    arglist[1] = 1.;  // Tolerance
    minuit.mnexcm("MINImize", arglist, 2, ierflg);

    // Extract the fitted parameter values and errors
    double ALU_sinphi1, ALU_sinphi1_error; minuit.GetParameter(0, ALU_sinphi1, ALU_sinphi1_error);
    double ALU_sinphi2, ALU_sinphi2_error; minuit.GetParameter(1, ALU_sinphi2, ALU_sinphi2_error);
    double ALU_sinDeltaphi, ALU_sinDeltaphi_error; minuit.GetParameter(2, ALU_sinDeltaphi, ALU_sinDeltaphi_error);
    double ALU_sin2Deltaphi, ALU_sin2Deltaphi_error; minuit.GetParameter(3, ALU_sin2Deltaphi, ALU_sin2Deltaphi_error);
    double AUL_sinphi1, AUL_sinphi1_error; minuit.GetParameter(4, AUL_sinphi1, AUL_sinphi1_error);
    double AUL_sinphi2, AUL_sinphi2_error; minuit.GetParameter(5, AUL_sinphi2, AUL_sinphi2_error);
    double AUL_sin2phi1, AUL_sin2phi1_error; minuit.GetParameter(6, AUL_sin2phi1, AUL_sin2phi1_error);
    double AUL_sin2phi2, AUL_sin2phi2_error; minuit.GetParameter(7, AUL_sin2phi2, AUL_sin2phi2_error);
    double AUL_sinDeltaphi, AUL_sinDeltaphi_error; minuit.GetParameter(8, AUL_sinDeltaphi, AUL_sinDeltaphi_error);
    double AUL_sin2Deltaphi, AUL_sin2Deltaphi_error; minuit.GetParameter(9, AUL_sin2Deltaphi, AUL_sin2Deltaphi_error);
    double AUL_sinSumphi, AUL_sinSumphi_error; minuit.GetParameter(10, AUL_sinSumphi, AUL_sinSumphi_error);
    double ALL, ALL_error; minuit.GetParameter(11, ALL, ALL_error);
    double ALL_cosphi1, ALL_cosphi1_error; minuit.GetParameter(12, ALL_cosphi1, ALL_cosphi1_error);
    double ALL_cosphi2, ALL_cosphi2_error; minuit.GetParameter(13, ALL_cosphi2, ALL_cosphi2_error);
    double AUU_cosphi1, AUU_cosphi1_error; minuit.GetParameter(14, AUU_cosphi1, AUU_cosphi1_error);
    double AUU_cosphi2, AUU_cosphi2_error; minuit.GetParameter(15, AUU_cosphi2, AUU_cosphi2_error);
    double AUU_cos2phi1, AUU_cos2phi1_error; minuit.GetParameter(16, AUU_cos2phi1, AUU_cos2phi1_error);
    double AUU_cos2phi2, AUU_cos2phi2_error; minuit.GetParameter(17, AUU_cos2phi2, AUU_cos2phi2_error);
    double AUU_cosSumphi, AUU_cosSumphi_error; minuit.GetParameter(18, AUU_cosSumphi, AUU_cosSumphi_error);

    // Calculate the mean values of the current variable 
    double sumVariable = 0;
    double numEvents = 0;
    TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
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
    mlmFitsStreams[0] << "{" << meanVariable << ", " << ALU_sinphi1 << ", " << ALU_sinphi1_error << "}";
    mlmFitsStreams[1] << "{" << meanVariable << ", " << ALU_sinphi2 << ", " << ALU_sinphi2_error << "}";
    mlmFitsStreams[2] << "{" << meanVariable << ", " << ALU_sinDeltaphi << ", " << ALU_sinDeltaphi_error << "}";
    mlmFitsStreams[3] << "{" << meanVariable << ", " << ALU_sin2Deltaphi << ", " << ALU_sin2Deltaphi_error << "}";
    mlmFitsStreams[4] << "{" << meanVariable << ", " << AUL_sinphi1 << ", " << AUL_sinphi1_error << "}";
    mlmFitsStreams[5] << "{" << meanVariable << ", " << AUL_sinphi2 << ", " << AUL_sinphi2_error << "}";
    mlmFitsStreams[6] << "{" << meanVariable << ", " << AUL_sin2phi1 << ", " << AUL_sin2phi1_error << "}";
    mlmFitsStreams[7] << "{" << meanVariable << ", " << AUL_sin2phi2 << ", " << AUL_sin2phi2_error << "}";
    mlmFitsStreams[8] << "{" << meanVariable << ", " << AUL_sinDeltaphi << ", " << AUL_sinDeltaphi_error << "}";
    mlmFitsStreams[9] << "{" << meanVariable << ", " << AUL_sin2Deltaphi << ", " << AUL_sin2Deltaphi_error << "}";
    mlmFitsStreams[10] << "{" << meanVariable << ", " << AUL_sinSumphi << ", " << AUL_sinSumphi_error << "}";
    mlmFitsStreams[11] << "{" << meanVariable << ", " << ALL << ", " << ALL_error << "}";
    mlmFitsStreams[12] << "{" << meanVariable << ", " << ALL_cosphi1 << ", " << ALL_cosphi1_error << "}";
    mlmFitsStreams[13] << "{" << meanVariable << ", " << ALL_cosphi2 << ", " << ALL_cosphi2_error << "}";
    mlmFitsStreams[14] << "{" << meanVariable << ", " << AUU_cosphi1 << ", " << AUU_cosphi1_error << "}";
    mlmFitsStreams[15] << "{" << meanVariable << ", " << AUU_cosphi2 << ", " << AUU_cosphi2_error << "}";
    mlmFitsStreams[16] << "{" << meanVariable << ", " << AUU_cos2phi1 << ", " << AUU_cos2phi1_error << "}";
    mlmFitsStreams[17] << "{" << meanVariable << ", " << AUU_cos2phi2 << ", " << AUU_cos2phi2_error << "}";
    mlmFitsStreams[18] << "{" << meanVariable << ", " << AUU_cosSumphi << ", " << AUU_cosSumphi_error << "}";

   
    if (i < numBins - 1) {
        mlmFitsStreams[0] << ", "; mlmFitsStreams[1] << ", "; mlmFitsStreams[2] << ", "; 
        mlmFitsStreams[3] << ", "; mlmFitsStreams[4] << ", "; mlmFitsStreams[5] << ", "; 
        mlmFitsStreams[6] << ", "; mlmFitsStreams[7] << ", "; mlmFitsStreams[8] << ", "; 
        mlmFitsStreams[9] << ", "; mlmFitsStreams[10] << ", "; mlmFitsStreams[11] << ", "; 
        mlmFitsStreams[12] << ", "; mlmFitsStreams[13] << ", "; mlmFitsStreams[14] << ", "; 
        mlmFitsStreams[15] << ", "; mlmFitsStreams[16] << ", "; mlmFitsStreams[17] << ", "; 
        mlmFitsStreams[18] << ", "; 
    }

    // outputs of asymmetries for LaTeX tables
    // Set fixed-point notation and one digit past the decimal
    ALUStream << std::fixed << std::setprecision(2); 
    ALUStream << (i+1) << " & " << meanVariable << " & ";
    // ALU sinphi1
    ALUStream << "$" << 100*ALU_sinphi1 << "_{" << TMath::Abs(100*0.068*ALU_sinphi1) << "}^{";
    ALUStream << 100*ALU_sinphi1_error << "}$ &";
    // ALU sinphi2
    ALUStream << "$" << 100*ALU_sinphi2 << "_{" << TMath::Abs(100*0.068*ALU_sinphi2) << "}^{";
    ALUStream << 100*ALU_sinphi2_error << "}$ &";
    // ALU sinDeltaphi
    ALUStream << "$" << 100*ALU_sinDeltaphi << "_{" << TMath::Abs(100*0.068*ALU_sinDeltaphi) << "}^{";
    ALUStream << 100*ALU_sinDeltaphi << "}$ &";
    // ALU sin2Deltaphi
    ALUStream << "$" << 100*ALU_sin2Deltaphi << "_{" << TMath::Abs(100*0.068*ALU_sin2Deltaphi) << "}^{";
    ALUStream << 100*ALU_sin2Deltaphi << "}$ &";
    //
    ALUStream << std::string(" \\\\ \\hline ");
    //

    // Set fixed-point notation and one digit past the decimal
    AULStream << std::fixed << std::setprecision(2); 
    AULStream << (i+1) << " & " << meanVariable << " & ";
    // AUL sinphi1
    AULStream << "$" << 100*AUL_sinphi1 << "_{" << TMath::Abs(100*0.092*AUL_sinphi1) << "}^{";
    AULStream << 100*AUL_sinphi1_error << "}$ &";
    //
    // AUL sinphi2
    AULStream << "$" << 100*AUL_sinphi2 << "_{" << TMath::Abs(100*0.092*AUL_sinphi2) << "}^{";
    AULStream << 100*AUL_sinphi2_error << "}$ &";
    //
    // AUL sin2phi1
    AULStream << "$" << 100*AUL_sin2phi1 << "_{" << TMath::Abs(100*0.092*AUL_sin2phi1) << "}^{";
    AULStream << 100*AUL_sin2phi1_error << "}$ &";
    //
    // AUL sin2phi2
    AULStream << "$" << 100*AUL_sin2phi2 << "_{" << TMath::Abs(100*0.092*AUL_sin2phi2) << "}^{";
    AULStream << 100*AUL_sin2phi2_error << "}$ &";
    // AUL sinDeltaphi
    AULStream << "$" << 100*AUL_sinDeltaphi << "_{" << TMath::Abs(100*0.092*AUL_sinDeltaphi) << "}^{";
    AULStream << 100*AUL_sinDeltaphi_error << "}$ &";
    //
    // AUL sin2Deltaphi
    AULStream << "$" << 100*AUL_sin2Deltaphi << "_{" << TMath::Abs(100*0.092*AUL_sin2phi2) << "}^{";
    AULStream << 100*AUL_sin2phi2_error << "}$ &";
    //
    AULStream << std::string(" \\\\ \\hline ");
    //


    // Set fixed-point notation and one digit past the decimal
    ALLStream << std::fixed << std::setprecision(2); 
    ALLStream << (i+1) << " & " << meanVariable << " & ";
    // ALL 
    ALLStream << "$" << 100*ALL << "_{" << TMath::Abs(100*0.092*ALL) << "}^{";
    ALLStream << 100*ALL_error << "}$ &";
    // ALL cosphi1
    ALLStream << "$" << 100*ALL_cosphi1 << "_{" << TMath::Abs(100*0.092*ALL_cosphi1) << "}^{";
    ALLStream << 100*ALL_cosphi1_error << "}$ &";
    // ALL cosphi2
    ALLStream << "$" << 100*ALL_cosphi2 << "_{" << TMath::Abs(100*0.092*ALL_cosphi2) << "}^{";
    ALLStream << 100*ALL_cosphi2_error << "}$ &";
    //
    ALUStream << std::string(" \\\\ \\hline ");
    //



    // Set fixed-point notation and one digit past the decimal
    AUUStream << std::fixed << std::setprecision(2); 
    AUUStream << (i+1) << " & " << meanVariable << " & ";
    // AUU_cosphi1 
    AUUStream << "$" << 100*AUU_cosphi1 << "_{" << TMath::Abs(100*0.5*AUU_cosphi1) << "}^{";
    AUUStream << 100*AUU_cosphi1_error << "}$ &";
    //
    // AUU_cosphi2 
    AUUStream << "$" << 100*AUU_cosphi2 << "_{" << TMath::Abs(100*0.5*AUU_cosphi2) << "}^{";
    AUUStream << 100*AUU_cosphi2_error << "}$ &";
    //
    // AUU_cos2phi1 
    AUUStream << "$" << 100*AUU_cos2phi1 << "_{" << TMath::Abs(100*0.5*AUU_cos2phi1) << "}^{";
    AUUStream << 100*AUU_cos2phi1_error << "}$ &";
    //
    // AUU_cos2phi2 
    AUUStream << "$" << 100*AUU_cos2phi2 << "_{" << TMath::Abs(100*0.5*AUU_cos2phi2) << "}^{";
    AUUStream << 100*AUU_cos2phi2_error << "}$ &";
    //
    // AUU_cosSumphi 
    AUUStream << "$" << 100*AUU_cosSumphi << "_{" << TMath::Abs(100*0.5*AUU_cosSumphi) << "}^{";
    AUUStream << 100*AUU_cosSumphi_error << "}$ &";
    //
    AUUStream << std::string(" \\\\ \\hline ");
    //
  }

  for (int i = 0; i < 19; ++i) {
    mlmFitsStreams[i] << "}; ";
  }
  std::ofstream outputFile(output_file, std::ios_base::app);
  for (int i = 0; i < 19; ++i) {
    outputFile << mlmFitsStreams[i].str() << std::endl;
  }

  outputFile.close();

  // Finally, close the table
  ALUStream << "\\end{tabular}" << std::endl;
  ALUStream << "\\caption{The mean kinematic value and the final ";
  ALUStream << "extracted structure function ratios for the beam-spin asymmetries as ";
  ALUStream << "a function of "  << prefix << ". Structure function ratios are given as ";
  ALUStream << "$100{A}_{\\pm100\\Delta\\text{sys}}^";
  ALUStream << "{\\pm100\\Delta\\text{stat}}$.}" << std::endl;
  ALUStream << "\\label{table:kinematics_" << prefix << "}" << std::endl;
  ALUStream << "\\end{table}" << std::endl;
  ALUStream << endl << endl << endl;
  //
  AULStream << "\\end{tabular}" << std::endl;
  AULStream << "\\caption{The mean kinematic value and the final ";
  AULStream << "extracted structure function ratios for the target-spin asymmetries as ";
  AULStream << "a function of "  << prefix << ". Structure function ratios are given as ";
  AULStream << "$100{A}_{\\pm100\\Delta\\text{sys}}^";
  AULStream << "{\\pm100\\Delta\\text{stat}}$.}" << std::endl;
  AULStream << "\\label{table:kinematics_" << prefix << "}" << std::endl;
  AULStream << "\\end{table}" << std::endl;
  AULStream << endl << endl << endl;
  //
  ALLStream << "\\end{tabular}" << std::endl;
  ALLStream << "\\caption{The mean kinematic value and the final ";
  ALLStream << "extracted structure function ratios for the double-spin asymmetries as ";
  ALLStream << "a function of "  << prefix << ". Structure function ratios are given as ";
  ALLStream << "$100{A}_{\\pm100\\Delta\\text{sys}}^";
  ALLStream << "{\\pm100\\Delta\\text{stat}}$.}" << std::endl;
  ALLStream << "\\label{table:kinematics_" << prefix << "}" << std::endl;
  ALLStream << "\\end{table}" << std::endl;
  ALLStream << endl << endl << endl;
  //
  AUUStream << "\\end{tabular}" << std::endl;
  AUUStream << "\\caption{The mean kinematic value and the final ";
  AUUStream << "extracted structure function ratios for the unpolarized modulations as ";
  AUUStream << "a function of "  << prefix << ". Structure function ratios are given as ";
  AUUStream << "$100{A}_{\\pm100\\Delta\\text{sys}}^";
  AUUStream << "{\\pm100\\Delta\\text{stat}}$.}" << std::endl;
  AUUStream << "\\label{table:kinematics_" << prefix << "}" << std::endl;
  AUUStream << "\\end{table}" << std::endl;
  AUUStream << endl << endl << endl;
  //
  std::ofstream kinematicFile(kinematic_file, std::ios_base::app);
  // Write the string stream content to the file
  kinematicFile << ALUStream.str() << std::endl; 
  kinematicFile << AULStream.str() << std::endl;
  kinematicFile << ALLStream.str() << std::endl;
  kinematicFile << AUUStream.str() << std::endl;
  kinematicFile.close();
}

TH2D* createHistogramForBin_b2b_dihadron(const char* histName, int binIndex, 
  const std::string& prefix, int asymmetry_index) {

  // Determine the variable range for the specified bin
  double varMin = allBins[currentFits][binIndex];
  double varMax = allBins[currentFits][binIndex + 1];

  // Create positive and negative helicity histograms
  int nPhiBins = 12; // Number of bins for phi variables
  TH2D* histPosPos = new TH2D(Form("%s_pospos", histName), "", 
    nPhiBins, 0, 2 * TMath::Pi(), nPhiBins, 0, 2 * TMath::Pi());
  TH2D* histPosNeg = new TH2D(Form("%s_posneg", histName), "", 
    nPhiBins, 0, 2 * TMath::Pi(), nPhiBins, 0, 2 * TMath::Pi());
  TH2D* histNegPos = new TH2D(Form("%s_negpos", histName), "", 
    nPhiBins, 0, 2 * TMath::Pi(), nPhiBins, 0, 2 * TMath::Pi());
  TH2D* histNegNeg = new TH2D(Form("%s_negneg", histName), "", 
    nPhiBins, 0, 2 * TMath::Pi(), nPhiBins, 0, 2 * TMath::Pi());

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
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> Q2(dataReader, "Q2");
  TTreeReaderValue<double> x(dataReader, "x");
  TTreeReaderValue<double> z(dataReader, "z");
  TTreeReaderValue<double> pT(dataReader, "pT");
  TTreeReaderValue<double> phi1(dataReader, "phi1");
  TTreeReaderValue<double> phi2(dataReader, "phi2");
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  // Counter to limit the number of processed entries
  while (dataReader.Next()) {

    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // bool passedKinematicCuts = true;
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
      sumVariable += *currentVariable;

      if (*helicity > 0 && *target_pol < 0) { histPosNeg->Fill(*phi1, *phi2); } 
      else if (*helicity < 0 && *target_pol > 0) {  histNegPos->Fill(*phi1, *phi2); }

      if (*helicity>0 && (*target_pol>0 || *runnum<11571) ) { histPosPos->Fill(*phi1,*phi2);} 
      else if (*helicity<0 && (*target_pol<0 || *runnum<11571) ) {histNegNeg->Fill(*phi1,*phi2);} 
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

  // Calculate and return the asymmetry histogram
  TH2D* histAsymmetry = new TH2D(Form("%s_asymmetry", histName), "", 
    nPhiBins, 0, 2 * TMath::Pi(), nPhiBins, 0, 2 * TMath::Pi());
  for (int iBinX = 1; iBinX <= nPhiBins; ++iBinX) {
    for (int iBinY = 1; iBinY <= nPhiBins; ++iBinY) {
      // Calculate asymmetry and error for each bin
      double Npp = histPosPos->GetBinContent(iBinX, iBinY) / cpp;
      double Npm = histPosNeg->GetBinContent(iBinX, iBinY) / cpp;
      double Nmp = histNegPos->GetBinContent(iBinX, iBinY) / cpp;
      double Nmm = histNegNeg->GetBinContent(iBinX, iBinY) / cpp;
      double asymmetry = asymmetry_value_calculation(meanVariable, prefix, 
        Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, asymmetry_index);
      double error = asymmetry_error_calculation(meanVariable, prefix, 
        Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, asymmetry_index);

      histAsymmetry->SetBinContent(iBinX, iBinY, asymmetry);
      histAsymmetry->SetBinError(iBinX, iBinY, error);
    }
  }

  // Delete the temporary positive and negative helicity histograms
  delete histPosPos;
  delete histPosNeg;
  delete histNegPos;
  delete histNegNeg;

  // Return the final asymmetry histogram
  return histAsymmetry;
}

void performChi2Fits_b2b_dihadron(const char* output_file, const char* kinematic_file,
  const std::string& prefix, int asymmetry_index) {

  // Initialize string streams for results and mean variables
  std::ostringstream chi2FitsStreams[8]; // For maximum 8 parameters (TSA case)
  for (auto& stream : chi2FitsStreams) {
    stream << std::fixed << std::setprecision(9);
  }

  // Initialize string streams to store the mean variables for each bin
  std::ostringstream meanVariablesStream;
  meanVariablesStream << "\\begin{table}[h]" << endl;
  meanVariablesStream << "\\centering" << endl;
  meanVariablesStream << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|} \\hline" << endl;
  meanVariablesStream << "Bin & $<Q^2>$ & $<W>$ ";
  meanVariablesStream << "& $<x_B>$ & $<y>$ & $<z1>$ & $<\\xi_2>$ & $<P_{1T}>$ ";
  meanVariablesStream << "& $<P_{2T}>$ & $<x_{F1}>$ & ";
  meanVariablesStream << "$<x_{F2}>$\\\\ \\hline" << endl; 

  // Create a new TF1 object called fitFunction representing the function to fit
  // and create string stream prefix depending on current asymmetry we're fitting
  TF2* fitFunction;
  // switch (asymmetry_index) {
  //   case 0: // beam-spin asymmetry
  //     fitFunction = new TF2("fitFunction", BSA_b2b_dihadron,0,2*TMath::Pi(),0,2*TMath::Pi(),5);
  //     chi2FitsStreams[0] << prefix << "chi2FitsALUoffset = {";
  //     chi2FitsStreams[1] << prefix << "chi2FitsALUsinphi1 = {";
  //     chi2FitsStreams[2] << prefix << "chi2FitsALUsinphi2 = {";
  //     chi2FitsStreams[3] << prefix << "chi2FitsALUsinDeltaphi = {";
  //     chi2FitsStreams[4] << prefix << "chi2FitsALUsin2Deltaphi = {";
  //     break;
  //   case 1: // target-spin asymmetry
  //     fitFunction = new TF2("fitFunction", TSA_b2b_dihadron,0,2*TMath::Pi(),0,2*TMath::Pi(),8);
  //     chi2FitsStreams[0] << prefix << "chi2FitsAULoffset = {";
  //     chi2FitsStreams[1] << prefix << "chi2FitsAULsinphi1 = {";
  //     chi2FitsStreams[2] << prefix << "chi2FitsAULsinphi2 = {";
  //     chi2FitsStreams[3] << prefix << "chi2FitsAULsin2phi1 = {";
  //     chi2FitsStreams[4] << prefix << "chi2FitsAULsin2phi2 = {";
  //     chi2FitsStreams[5] << prefix << "chi2FitsAULsinDeltaphi = {";
  //     chi2FitsStreams[6] << prefix << "chi2FitsAULsin2Deltaphi = {";
  //     chi2FitsStreams[7] << prefix << "chi2FitsAULsinSumphi = {";
  //     break;
  //   case 2: // double-spin asymmetry
  //     fitFunction = new TF2("fitFunction", DSA_b2b_dihadron,0,2*TMath::Pi(),0,2*TMath::Pi(),3);
  //     chi2FitsStreams[0] << prefix << "chi2FitsALL = {";
  //     chi2FitsStreams[1] << prefix << "chi2FitsALLcosphi1 = {";
  //     chi2FitsStreams[2] << prefix << "chi2FitsALLcosphi2 = {";
  //     break;
  // }
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      fitFunction = new TF2("fitFunction", BSA_b2b_dihadron,0,2*TMath::Pi(),0,2*TMath::Pi(),3);
      chi2FitsStreams[0] << prefix << "chi2FitsALUoffset = {";
      chi2FitsStreams[1] << prefix << "chi2FitsALUsinDeltaphi = {";
      chi2FitsStreams[2] << prefix << "chi2FitsALUsin2Deltaphi = {";
      break;
    case 1: // target-spin asymmetry
      fitFunction = new TF2("fitFunction", TSA_b2b_dihadron,0,2*TMath::Pi(),0,2*TMath::Pi(),4);
      chi2FitsStreams[0] << prefix << "chi2FitsAULoffset = {";
      chi2FitsStreams[1] << prefix << "chi2FitsAULsinDeltaphi = {";
      chi2FitsStreams[2] << prefix << "chi2FitsAULsin2Deltaphi = {";
      chi2FitsStreams[3] << prefix << "chi2FitsAULsinSumphi = {";
      break;
    case 2: // double-spin asymmetry
      fitFunction = new TF2("fitFunction", DSA_b2b_dihadron,0,2*TMath::Pi(),0,2*TMath::Pi(),1);
      chi2FitsStreams[0] << prefix << "chi2FitsALL = {";
      break;
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
    TH2D* hist = createHistogramForBin_b2b_dihadron(histName, i, prefix, asymmetry_index);
    // Fit the histogram using the fitFunction and get the fit result
    hist->Fit(fitFunction, "QS");
    // not plotting function here for 2D dihadron cases

    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double numEvents = 0;
    // Variables to calculate the mean depolarization factor
    double sumDepA = 0; double sumDepB = 0; 
    double sumDepC = 0; double sumDepV = 0; double sumDepW = 0;

    // Variables to calculate the mean kinematics in each bin
    double sumQ2 = 0; double sumW = 0; double sumx = 0; double sumy = 0;
    double sumz1 = 0; double sumxi2 = 0; double sumpT1 = 0; double sumpT2 = 0; 
    double sumxF1 = 0; double sumxF2 = 0;

    // Declare reader locations
    TTreeReaderValue<int> runnum(dataReader, "runnum");
    TTreeReaderValue<int> evnum(dataReader, "evnum");
    TTreeReaderValue<double> Q2(dataReader, "Q2");
    TTreeReaderValue<double> W(dataReader, "W");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> y(dataReader, "y");
    TTreeReaderValue<double> z1(dataReader, "z1");
    TTreeReaderValue<double> xi2(dataReader, "xi2");
    TTreeReaderValue<double> pT1(dataReader, "pT1");
    TTreeReaderValue<double> pT2(dataReader, "pT2");
    TTreeReaderValue<double> xF1(dataReader, "xF1");
    TTreeReaderValue<double> xF2(dataReader, "xF2");
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
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
      // Check if the currentVariable is within the desired range
      if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
        // sum the kinematic variable values
        sumVariable += *currentVariable;
        sumQ2 += *Q2;
        sumW += *W;
        sumx += *x;
        sumy += *y;
        sumz1 += *z1;
        sumxi2 += *xi2;
        sumpT1 += *pT1;
        sumpT2 += *pT2;
        sumxF1 += *xF1;
        sumxF2 += *xF2;

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
    double meanz1 = numEvents > 0 ? sumz1 / numEvents : 0.0;
    double meanxi2 = numEvents > 0 ? sumxi2 / numEvents : 0.0;
    double meanpT1 = numEvents > 0 ? sumpT1 / numEvents : 0.0;
    double meanpT2 = numEvents > 0 ? sumpT2 / numEvents : 0.0;
    double meanxF1 = numEvents > 0 ? sumxF1 / numEvents : 0.0;
    double meanxF2 = numEvents > 0 ? sumxF2 / numEvents : 0.0;

    switch (asymmetry_index) {
      case 0: {// beam-spin asymmetry
        // Get the fitted parameters and their errors
        double ALU_offset = fitFunction->GetParameter(0);
        double ALU_offset_error = fitFunction->GetParError(0);
        // double ALU_sinphi1 = fitFunction->GetParameter(1); 
        // double ALU_sinphi1_error = fitFunction->GetParError(1);
        // double ALU_sinphi2 = fitFunction->GetParameter(2); 
        // double ALU_sinphi2_error = fitFunction->GetParError(2);
        double ALU_sinDeltaphi = fitFunction->GetParameter(1); 
        double ALU_sinDeltaphi_error = fitFunction->GetParError(1);
        double ALU_sin2Deltaphi = fitFunction->GetParameter(2); 
        double ALU_sin2Deltaphi_error = fitFunction->GetParError(2);
        //
        // ALU_sinphi1 = (meanDepA/meanDepW)*ALU_sinphi1;
        // ALU_sinphi1_error = (meanDepA/meanDepW)*ALU_sinphi1_error;
        // //
        // ALU_sinphi2 = (meanDepA/meanDepW)*ALU_sinphi2;
        // ALU_sinphi2_error = (meanDepA/meanDepW)*ALU_sinphi2_error;
        //
        ALU_sinDeltaphi = (meanDepA/meanDepC)*ALU_sinDeltaphi;
        ALU_sinDeltaphi_error = (meanDepA/meanDepC)*ALU_sinDeltaphi_error;
        //
        ALU_sin2Deltaphi = (meanDepA/meanDepC)*ALU_sin2Deltaphi;
        ALU_sin2Deltaphi_error = (meanDepA/meanDepC)*ALU_sin2Deltaphi_error;
        //
        chi2FitsStreams[0]<<"{"<<meanVariable<<", "<< ALU_offset << ", " << ALU_offset_error <<"}";
        // chi2FitsStreams[1]<<"{"<<meanVariable<<", "<< ALU_sinphi1 << ", " << ALU_sinphi1_error <<"}";
        // chi2FitsStreams[2]<<"{"<<meanVariable<<", "<< ALU_sinphi2 << ", " << ALU_sinphi2_error <<"}";
        chi2FitsStreams[1]<<"{"<<meanVariable<<", "<< ALU_sinDeltaphi << ", " << ALU_sinDeltaphi_error <<"}";
        chi2FitsStreams[2]<<"{"<<meanVariable<<", "<< ALU_sin2Deltaphi << ", " << ALU_sin2Deltaphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsStreams[0] << ", "; chi2FitsStreams[1] << ", "; 
            chi2FitsStreams[2] << ", "; 
            // chi2FitsStreams[3] << ", "; 
            // chi2FitsStreams[4] << ", ";
        }
        break;
      }
      case 1: {// target-spin asymmetry
        // Get the fitted parameters and their errors
        double AUL_offset = fitFunction->GetParameter(0);
        double AUL_offset_error = fitFunction->GetParError(0);
        // double AUL_sinphi1 = fitFunction->GetParameter(1);
        // double AUL_sinphi1_error = fitFunction->GetParError(1);
        // double AUL_sinphi2 = fitFunction->GetParameter(2);
        // double AUL_sinphi2_error = fitFunction->GetParError(2);
        // double AUL_sin2phi1 = fitFunction->GetParameter(3);
        // double AUL_sin2phi1_error = fitFunction->GetParError(3);
        // double AUL_sin2phi2 = fitFunction->GetParameter(4);
        // double AUL_sin2phi2_error = fitFunction->GetParError(4);
        double AUL_sinDeltaphi = fitFunction->GetParameter(1);
        double AUL_sinDeltaphi_error = fitFunction->GetParError(1);
        double AUL_sin2Deltaphi = fitFunction->GetParameter(2);
        double AUL_sin2Deltaphi_error = fitFunction->GetParError(2);
        double AUL_sinSumphi = fitFunction->GetParameter(3);
        double AUL_sinSumphi_error = fitFunction->GetParError(3);
        //
        // AUL_sinphi1 = (meanDepA/meanDepV)*AUL_sinphi1;
        // AUL_sinphi1_error = (meanDepA/meanDepV)*AUL_sinphi1_error;
        // //
        // AUL_sinphi2 = (meanDepA/meanDepV)*AUL_sinphi2;
        // AUL_sinphi2_error = (meanDepA/meanDepV)*AUL_sinphi2_error;
        // //
        // AUL_sin2phi1 = (meanDepA/meanDepB)*AUL_sin2phi1;
        // AUL_sin2phi1_error = (meanDepA/meanDepB)*AUL_sin2phi1_error;
        // //
        // AUL_sin2phi2 = (meanDepA/meanDepB)*AUL_sin2phi2;
        // AUL_sin2phi2_error = (meanDepA/meanDepB)*AUL_sin2phi2_error;
        //
        // No depolarization factor for the sin(Deltaphi) and sin(2Deltaphi) asymmetries
        // in the TSA
        //
        AUL_sinSumphi = (meanDepA/meanDepB)*AUL_sinSumphi;
        AUL_sinSumphi_error = (meanDepA/meanDepB)*AUL_sinSumphi_error;
        //
        chi2FitsStreams[0]<<"{"<<meanVariable<<", "<< AUL_offset << ", " << AUL_offset_error <<"}";
        // chi2FitsStreams[1]<<"{"<<meanVariable<<", "<< AUL_sinphi1 << ", " << AUL_sinphi1_error <<"}";
        // chi2FitsStreams[2]<<"{"<<meanVariable<<", "<< AUL_sinphi2 << ", " << AUL_sinphi2_error <<"}";
        // chi2FitsStreams[3]<<"{"<<meanVariable<<", "<< AUL_sin2phi1 << ", " << AUL_sin2phi1_error <<"}";
        // chi2FitsStreams[4]<<"{"<<meanVariable<<", "<< AUL_sin2phi2 << ", " << AUL_sin2phi2_error <<"}";
        chi2FitsStreams[1]<<"{"<<meanVariable<<", "<< AUL_sinDeltaphi << ", " << AUL_sinDeltaphi_error <<"}";
        chi2FitsStreams[2]<<"{"<<meanVariable<<", "<< AUL_sin2Deltaphi << ", " << AUL_sin2Deltaphi_error <<"}";
        chi2FitsStreams[3]<<"{"<<meanVariable<<", "<< AUL_sinSumphi << ", " << AUL_sinSumphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsStreams[0] << ", "; chi2FitsStreams[1] << ", "; chi2FitsStreams[2] << ", ";
            chi2FitsStreams[3] << ", "; 
            // chi2FitsStreams[4] << ", "; chi2FitsStreams[5] << ", ";
            // chi2FitsStreams[6] << ", "; chi2FitsStreams[7] << ", "; 
        }
        break;
      }
      case 2: {// double-spin asymmetry
        // Get the fitted parameters and their errors
        double ALL = fitFunction->GetParameter(0);
        double ALL_error = fitFunction->GetParError(0);
        // double ALL_cosphi1 = fitFunction->GetParameter(1);
        // double ALL_cosphi1_error = fitFunction->GetParError(1);
        // double ALL_cosphi2 = fitFunction->GetParameter(2);
        // double ALL_cosphi2_error = fitFunction->GetParError(2);
        //
        ALL = (meanDepA/meanDepC)*ALL;
        ALL_error = (meanDepA/meanDepC)*ALL_error;
        //
        // ALL_cosphi1 = (meanDepA/meanDepW)*ALL_cosphi1;
        // ALL_cosphi1_error = (meanDepA/meanDepW)*ALL_cosphi1_error;
        // //
        // ALL_cosphi2 = (meanDepA/meanDepW)*ALL_cosphi2;
        // ALL_cosphi2_error = (meanDepA/meanDepW)*ALL_cosphi2_error;
        //
        chi2FitsStreams[0]<<"{"<<meanVariable<<", "<< ALL << ", " << ALL_error <<"}";
        // chi2FitsStreams[1]<<"{"<<meanVariable<<", "<< ALL_cosphi1 << ", " << ALL_cosphi1_error <<"}";
        // chi2FitsStreams[2]<<"{"<<meanVariable<<", "<< ALL_cosphi2 << ", " << ALL_cosphi2_error <<"}";
        if (i < numBins - 1) {
            chi2FitsStreams[0] << ", "; 
            // chi2FitsStreams[1] << ", "; chi2FitsStreams[2] << ", ";
            // chi2FitsStreams[3] << ", "; chi2FitsStreams[4] << ", "; chi2FitsStreams[5] << ", ";
            // chi2FitsStreams[6] << ", "; chi2FitsStreams[7] << ", "; 
        }
        break;
      }
    }

    delete hist;

    // outputs of mean kinematic variables
    meanVariablesStream << std::fixed << std::setprecision(3); // Set precision to 3 digits 
    meanVariablesStream << (i+1) << "~&~" << meanQ2 << "~&~" << meanW << "~&~" << meanx << "~&~";
    meanVariablesStream << meany << "~&~" << meanz1 << "~&~" << meanxi2 << "~&~";
    meanVariablesStream << meanpT1 << "~&~" << meanpT2 << "~&~" << meanxF1 << "~&~" << meanxF2; 
    meanVariablesStream << std::string(" \\\\ \\hline ");
  }

  chi2FitsStreams[0] << "};";  chi2FitsStreams[1] << "};";  chi2FitsStreams[2] << "};"; 
  chi2FitsStreams[3] << "};";  chi2FitsStreams[4] << "};";  chi2FitsStreams[5] << "};"; 
  chi2FitsStreams[6] << "};";  chi2FitsStreams[7] << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << chi2FitsStreams[0].str() << std::endl;
  outputFile << chi2FitsStreams[1].str() << std::endl;
  outputFile << chi2FitsStreams[2].str() << std::endl;
  if (asymmetry_index<2) { 
    outputFile << chi2FitsStreams[3].str() << std::endl;
    outputFile << chi2FitsStreams[4].str() << std::endl;
  }
  if (asymmetry_index==1) {
    outputFile << chi2FitsStreams[5].str() << std::endl;
    outputFile << chi2FitsStreams[6].str() << std::endl;
    outputFile << chi2FitsStreams[7].str() << std::endl;
  }

  outputFile.close();

  meanVariablesStream << "\\end{tabular}\n";
  meanVariablesStream << "\\caption{The mean kinematic variables in each of the bins ";
  meanVariablesStream << "for the extracted $" << prefix << "$ asymmetries.";
  meanVariablesStream << " Values given in GeV or GeV$^2$ where appropriate.}\n";
  meanVariablesStream << "\\end{table}\n";
  meanVariablesStream << endl << endl << endl;
  if (asymmetry_index == 0) {
    std::ofstream kinematicFile(kinematic_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicFile << meanVariablesStream.str() << std::endl; 
    kinematicFile.close();
  }
}

/******************** DVCS CASE ********************/

void plotHistogramAndFit_dvcs(TH1D* histogram, TF1* fitFunction, int binIndex, 
  int asymmetryIndex, const std::string& prefix) {
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
  TLegend *leg = new TLegend(1-0.19, 0.675, 1-0.45, 0.875);  // Adjusted to the upper-right corner
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
        if (i == 2) paramName = "#it{A}_{UU}^{cos#phi}";
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
  string filename = "output/individual_chi2_fits/" + prefix + "_" + 
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
  if (canvas->GetListOfPrimitives()->GetSize() > 0) {
      // There's something in the canvas, save it
      canvas->SaveAs(filename.c_str());
  } else {
      std::cout << "Canvas is empty, not saving to file." << std::endl;
  }

  // Clean up
  delete canvas;
  delete graph;
}

TH1D* createHistogramForBin_dvcs(const char* histName, int binIndex, 
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
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> Q2(dataReader, "Q2");
  TTreeReaderValue<double> x(dataReader, "x");
  TTreeReaderValue<double> z(dataReader, "z");
  TTreeReaderValue<double> pT(dataReader, "pT");
  TTreeReaderValue<double> phi(dataReader, "phi2"); 
  // this is phi2 because we're using processing_dihadron to identify proton and photon 
  // (which isn't really a hadron of course)
  // so phi2 is the dvcs photon angle
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  // Counter to limit the number of processed entries
  while (dataReader.Next()) {

    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // bool passedKinematicCuts = true;
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
      sumVariable += *currentVariable;

      if (*helicity > 0 && *target_pol < 0) { histPosNeg->Fill(*phi); } 
      else if (*helicity < 0 && *target_pol > 0) {  histNegPos->Fill(*phi); }

      if (*helicity > 0 && (*target_pol >= 0) ) { histPosPos->Fill(*phi); } 
      else if (*helicity < 0 && (*target_pol <= 0) ) {  histNegNeg->Fill(*phi); } 
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
    double asymmetry = asymmetry_value_calculation(meanVariable, prefix, 
      Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, asymmetry_index);
    double error = asymmetry_error_calculation(meanVariable, prefix, 
      Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, asymmetry_index);

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

void performChi2Fits_dvcs(const char* output_file, const char* kinematic_file,
  const char* kinematicPlot_file, const std::string& prefix, int asymmetry_index) {

  // Initialize string streams to store the results for each bin
  std::ostringstream chi2FitsAStream, chi2FitsBStream, chi2FitsCStream;
  chi2FitsAStream << std::fixed << std::setprecision(9);
  chi2FitsBStream << std::fixed << std::setprecision(9);
  chi2FitsCStream << std::fixed << std::setprecision(9);

  // Initialize string stream to store the kinematics in each bin for use in LaTeX 
  std::ostringstream meanVariablesStream;
  meanVariablesStream << "\\begin{table}[h]" << endl;
  meanVariablesStream << "\\centering" << endl;
  meanVariablesStream << "\\begin{tabular}{|c|c|c|c|c|} \\hline" << endl;
  meanVariablesStream << "Bin & $<Q^2>$ & $<W>$ ";
  meanVariablesStream << "& $<x_B>$ & $<y>$ & ";
  meanVariablesStream << "& $<t>$ \\hline" << endl; 

  // Initalize string stream to store the kinematics in each bin for use in plotting 
  std::ostringstream meanVariablesPlotStream;
  meanVariablesPlotStream << prefix << "Kinematics = {";

  // Create a new TF1 object called fitFunction representing the function to fit
  // and create string stream prefix depending on current asymmetry we're fitting
  TF1* fitFunction;
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      fitFunction = new TF1("fitFunction", BSA_dvcs, 0, 2*TMath::Pi(), 3);
      chi2FitsAStream << prefix << "chi2FitsALUoffset = {";
      chi2FitsBStream << prefix << "chi2FitsALUsinphi = {";
      chi2FitsCStream << prefix << "chi2FitsAUUcosphi = {";
      break;
    case 1: // target-spin asymmetry
      fitFunction = new TF1("fitFunction", TSA_single_hadron, 0, 2*TMath::Pi(), 3);
      chi2FitsAStream << prefix << "chi2FitsAULoffset = {";
      chi2FitsBStream << prefix << "chi2FitsAULsinphi = {";
      chi2FitsCStream << prefix << "chi2FitsAULsin2phi = {";
      break;
    case 2: // double-spin asymmetry
      fitFunction = new TF1("fitFunction", DSA_single_hadron, 0, 2*TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALL = {";
      chi2FitsBStream << prefix << "chi2FitsALLcosphi = {";
      break;
    default:
      cout << "Invalid asymmetry_index! Using default function form of BSA." << endl;
      fitFunction = new TF1("fitFunction", BSA_single_hadron, 0, 2*TMath::Pi(), 2);
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
    TH1D* hist = createHistogramForBin_dvcs(histName, i, prefix, asymmetry_index);
    // Fit the histogram using the fitFunction and get the fit result
    hist->Fit(fitFunction, "QS");
    plotHistogramAndFit_dvcs(hist, fitFunction, i, asymmetry_index, prefix);

    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double numEvents = 0;
    // Variables to calculate the mean depolarization factor
    double sumDepA = 0; double sumDepB = 0; 
    double sumDepC = 0; double sumDepV = 0; double sumDepW = 0;

    // Variables to calculate the mean kinematics in each bin
    double sumQ2 = 0; double sumW = 0; double sumx = 0; double sumy = 0;
    double sumt = 0; 

    // Declare reader locations
    TTreeReaderValue<int> runnum(dataReader, "runnum");
    TTreeReaderValue<int> evnum(dataReader, "evnum");
    TTreeReaderValue<double> Q2(dataReader, "Q2");
    TTreeReaderValue<double> W(dataReader, "W");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> y(dataReader, "y");
    TTreeReaderValue<double> t(dataReader, "t");
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
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
      // Check if the currentVariable is within the desired range
      if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
        // sum the kinematic variable values
        sumVariable += *currentVariable;
        sumQ2 += *Q2;
        sumW += *W;
        sumx += *x;
        sumy += *y;
        sumt += *t;

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
    double meant = numEvents > 0 ? sumt / numEvents : 0.0;

    switch (asymmetry_index) {
      case 0: {// beam-spin asymmetry
        // Get the fitted parameters and their errors
        double ALU_offset = fitFunction->GetParameter(0);
        double ALU_offset_error = fitFunction->GetParError(0);
        double ALU_sinphi = fitFunction->GetParameter(1); 
        double ALU_sinphi_error = fitFunction->GetParError(1);
        double AUU_cosphi = fitFunction->GetParameter(2); 
        double AUU_cosphi_error = fitFunction->GetParError(2);
        // ALU_sinphi = (meanDepA/meanDepW)*ALU_sinphi;
        // ALU_sinphi_error = (meanDepA/meanDepW)*ALU_sinphi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALU_offset << ", " << ALU_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALU_sinphi << ", " << ALU_sinphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", " << AUU_cosphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; chi2FitsCStream << ", "; 
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
        // AUL_sinphi = (meanDepA/meanDepV)*AUL_sinphi;
        // AUL_sinphi_error = (meanDepA/meanDepV)*AUL_sinphi_error;
        // AUL_sin2phi = (meanDepA/meanDepB)*AUL_sin2phi;
        // AUL_sin2phi_error = (meanDepA/meanDepB)*AUL_sin2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< AUL_offset << ", " << AUL_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< AUL_sinphi << ", " << AUL_sinphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUL_sin2phi << ", " << AUL_sin2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; chi2FitsCStream << ", ";
        }
        break;
      }
      case 2: {// double-spin asymmetry
        // Get the fitted parameters and their errors
        double ALL = fitFunction->GetParameter(0);
        double ALL_error = fitFunction->GetParError(0);
        double ALL_cosphi = fitFunction->GetParameter(1);
        double ALL_cosphi_error = fitFunction->GetParError(1);
        // ALL = (meanDepA/meanDepC)*ALL;
        // ALL_error = (meanDepA/meanDepC)*ALL_error;
        // ALL_cosphi = (meanDepA/meanDepW)*ALL_cosphi;
        // ALL_cosphi_error = (meanDepA/meanDepW)*ALL_cosphi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALL << ", " << ALL_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALL_cosphi << ", " << ALL_cosphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", ";
        }
        break;
      }
    }

    delete hist;

    // outputs of mean kinematic variables for LaTeX
    meanVariablesStream << std::fixed << std::setprecision(3); // Set precision to 3 digits 
    meanVariablesStream << (i+1) << "~&~" << meanQ2 << "~&~" << meanW << "~&~" << meanx << "~&~";
    meanVariablesStream << meany << "~&~" << meant; 
    meanVariablesStream << std::string(" \\\\ \\hline ");

    // outputs of mean kinematic variables for plotting
    meanVariablesPlotStream << "{" << meanQ2 << ", " << meanW << ", " << meanx << ", ";
    meanVariablesPlotStream << meany << ", "; meanVariablesPlotStream << meant << "}";
    if (i < numBins - 1) {
        meanVariablesPlotStream << ", "; 
    }
  }

  chi2FitsAStream << "};";  chi2FitsBStream << "};";  chi2FitsCStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << chi2FitsAStream.str() << std::endl;
  outputFile << chi2FitsBStream.str() << std::endl;
  if (asymmetry_index==0 || asymmetry_index==1) { outputFile << chi2FitsCStream.str() << std::endl; }

  outputFile.close();

  meanVariablesStream << "\\end{tabular}\n";
  meanVariablesStream << "\\caption{The mean kinematic variables in each of the bins ";
  meanVariablesStream << "for the extracted $" << prefix << "$ asymmetries.";
  meanVariablesStream << " Values given in GeV or GeV$^2$ where appropriate.}\n";
  meanVariablesStream << "\\end{table}\n";
  meanVariablesStream << endl << endl << endl;
  if (asymmetry_index == 0) {
    std::ofstream kinematicFile(kinematic_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicFile << meanVariablesStream.str() << std::endl; 
    kinematicFile.close();

    meanVariablesPlotStream << "};";
    std::ofstream kinematicPlot_File(kinematicPlot_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicPlot_File << meanVariablesPlotStream.str() << std::endl;
    kinematicPlot_File.close();
  }
}

/******************** exclusive eppi0 CASE ********************/

void plotHistogramAndFit_eppi0(TH1D* histogram, TF1* fitFunction, int binIndex, 
  int asymmetryIndex, const std::string& prefix) {
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
  TLegend *leg = new TLegend(1-0.19, 0.675, 1-0.45, 0.875);  // Adjusted to the upper-right corner
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
        if (i == 2) paramName = "#it{A}_{UU}^{cos#phi}";
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
  string filename = "output/individual_chi2_fits/" + prefix + "_" + 
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
  if (canvas->GetListOfPrimitives()->GetSize() > 0) {
      // There's something in the canvas, save it
      canvas->SaveAs(filename.c_str());
  } else {
      std::cout << "Canvas is empty, not saving to file." << std::endl;
  }

  // Clean up
  delete canvas;
  delete graph;
}

TH1D* createHistogramForBin_eppi0(const char* histName, int binIndex, 
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
  TTreeReaderValue<int> helicity(dataReader, "helicity");
  TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> Q2(dataReader, "Q2");
  TTreeReaderValue<double> x(dataReader, "x");
  TTreeReaderValue<double> z(dataReader, "z");
  TTreeReaderValue<double> pT(dataReader, "pT");
  TTreeReaderValue<double> phi(dataReader, "phi2"); 
  // TTreeReaderValue<double> phi(dataReader, "gamma_phi1"); 
  // TTreeReaderValue<double> phi(dataReader, "gamma_phi2"); 
  // this is phi2 because we're using processing_dihadron to identify proton and photon/eppi0 
  // (which isn't really a hadron of course)
  // so phi2 is the dvcs/eppi0 photon angle
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  // Counter to limit the number of processed entries
  while (dataReader.Next()) {

    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // bool passedKinematicCuts = true;
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
      sumVariable += *currentVariable;

      if (*helicity > 0 && *target_pol < 0) { histPosNeg->Fill(*phi); } 
      else if (*helicity < 0 && *target_pol > 0) {  histNegPos->Fill(*phi); }

      if (*helicity > 0 && (*target_pol >= 0) ) { histPosPos->Fill(*phi); } 
      else if (*helicity < 0 && (*target_pol <= 0) ) {  histNegNeg->Fill(*phi); } 
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
    double asymmetry = asymmetry_value_calculation(meanVariable, prefix, 
      Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, asymmetry_index);
    double error = asymmetry_error_calculation(meanVariable, prefix, 
      Npp, Npm, Nmp, Nmm, meanPol, Ptp, Ptm, asymmetry_index);

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

void performChi2Fits_eppi0(const char* output_file, const char* kinematic_file,
  const char* kinematicPlot_file, const std::string& prefix, int asymmetry_index) {

  // Initialize string streams to store the results for each bin
  std::ostringstream chi2FitsAStream, chi2FitsBStream, chi2FitsCStream;
  chi2FitsAStream << std::fixed << std::setprecision(9);
  chi2FitsBStream << std::fixed << std::setprecision(9);
  chi2FitsCStream << std::fixed << std::setprecision(9);

  // Initialize string stream to store the kinematics in each bin for use in LaTeX 
  std::ostringstream meanVariablesStream;
  meanVariablesStream << "\\begin{table}[h]" << endl;
  meanVariablesStream << "\\centering" << endl;
  meanVariablesStream << "\\begin{tabular}{|c|c|c|c|c|} \\hline" << endl;
  meanVariablesStream << "Bin & $<Q^2>$ & $<W>$ ";
  meanVariablesStream << "& $<x_B>$ & $<y>$ & ";
  meanVariablesStream << "& $<t>$ \\hline" << endl; 

  // Initalize string stream to store the kinematics in each bin for use in plotting 
  std::ostringstream meanVariablesPlotStream;
  meanVariablesPlotStream << prefix << "Kinematics = {";

  // Create a new TF1 object called fitFunction representing the function to fit
  // and create string stream prefix depending on current asymmetry we're fitting
  TF1* fitFunction;
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      fitFunction = new TF1("fitFunction", BSA_dvcs, 0, 2*TMath::Pi(), 3);
      chi2FitsAStream << prefix << "chi2FitsALUoffset = {";
      chi2FitsBStream << prefix << "chi2FitsALUsinphi = {";
      chi2FitsCStream << prefix << "chi2FitsAUUcosphi = {";
      break;
    case 1: // target-spin asymmetry
      fitFunction = new TF1("fitFunction", TSA_single_hadron, 0, 2*TMath::Pi(), 3);
      chi2FitsAStream << prefix << "chi2FitsAULoffset = {";
      chi2FitsBStream << prefix << "chi2FitsAULsinphi = {";
      chi2FitsCStream << prefix << "chi2FitsAULsin2phi = {";
      break;
    case 2: // double-spin asymmetry
      fitFunction = new TF1("fitFunction", DSA_single_hadron, 0, 2*TMath::Pi(), 2);
      chi2FitsAStream << prefix << "chi2FitsALL = {";
      chi2FitsBStream << prefix << "chi2FitsALLcosphi = {";
      break;
    default:
      cout << "Invalid asymmetry_index! Using default function form of BSA." << endl;
      fitFunction = new TF1("fitFunction", BSA_single_hadron, 0, 2*TMath::Pi(), 2);
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
    TH1D* hist = createHistogramForBin_eppi0(histName, i, prefix, asymmetry_index);
    // Fit the histogram using the fitFunction and get the fit result
    hist->Fit(fitFunction, "QS");
    plotHistogramAndFit_eppi0(hist, fitFunction, i, asymmetry_index, prefix);

    // Initialize variables to store the sums and event counts
    double sumVariable = 0;
    double numEvents = 0;
    // Variables to calculate the mean depolarization factor
    double sumDepA = 0; double sumDepB = 0; 
    double sumDepC = 0; double sumDepV = 0; double sumDepW = 0;

    // Variables to calculate the mean kinematics in each bin
    double sumQ2 = 0; double sumW = 0; double sumx = 0; double sumy = 0;
    double sumt = 0; 

    // Declare reader locations
    TTreeReaderValue<int> runnum(dataReader, "runnum");
    TTreeReaderValue<int> evnum(dataReader, "evnum");
    TTreeReaderValue<double> Q2(dataReader, "Q2");
    TTreeReaderValue<double> W(dataReader, "W");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> y(dataReader, "y");
    TTreeReaderValue<double> t(dataReader, "t");
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
    while (dataReader.Next()) {
      // Apply kinematic cuts (this function will need to be adapted)
      bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
      // Check if the currentVariable is within the desired range
      if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
        // sum the kinematic variable values
        sumVariable += *currentVariable;
        sumQ2 += *Q2;
        sumW += *W;
        sumx += *x;
        sumy += *y;
        sumt += *t;

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
    double meant = numEvents > 0 ? sumt / numEvents : 0.0;

    switch (asymmetry_index) {
      case 0: {// beam-spin asymmetry
        // Get the fitted parameters and their errors
        double ALU_offset = fitFunction->GetParameter(0);
        double ALU_offset_error = fitFunction->GetParError(0);
        double ALU_sinphi = fitFunction->GetParameter(1); 
        double ALU_sinphi_error = fitFunction->GetParError(1);
        double AUU_cosphi = fitFunction->GetParameter(2); 
        double AUU_cosphi_error = fitFunction->GetParError(2);
        // ALU_sinphi = (meanDepA/meanDepW)*ALU_sinphi;
        // ALU_sinphi_error = (meanDepA/meanDepW)*ALU_sinphi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALU_offset << ", " << ALU_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALU_sinphi << ", " << ALU_sinphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", " << AUU_cosphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; chi2FitsCStream << ", "; 
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
        // AUL_sinphi = (meanDepA/meanDepV)*AUL_sinphi;
        // AUL_sinphi_error = (meanDepA/meanDepV)*AUL_sinphi_error;
        // AUL_sin2phi = (meanDepA/meanDepB)*AUL_sin2phi;
        // AUL_sin2phi_error = (meanDepA/meanDepB)*AUL_sin2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< AUL_offset << ", " << AUL_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< AUL_sinphi << ", " << AUL_sinphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUL_sin2phi << ", " << AUL_sin2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; chi2FitsCStream << ", ";
        }
        break;
      }
      case 2: {// double-spin asymmetry
        // Get the fitted parameters and their errors
        double ALL = fitFunction->GetParameter(0);
        double ALL_error = fitFunction->GetParError(0);
        double ALL_cosphi = fitFunction->GetParameter(1);
        double ALL_cosphi_error = fitFunction->GetParError(1);
        // ALL = (meanDepA/meanDepC)*ALL;
        // ALL_error = (meanDepA/meanDepC)*ALL_error;
        // ALL_cosphi = (meanDepA/meanDepW)*ALL_cosphi;
        // ALL_cosphi_error = (meanDepA/meanDepW)*ALL_cosphi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALL << ", " << ALL_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALL_cosphi << ", " << ALL_cosphi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", ";
        }
        break;
      }
    }

    delete hist;

    // outputs of mean kinematic variables for LaTeX
    meanVariablesStream << std::fixed << std::setprecision(3); // Set precision to 3 digits 
    meanVariablesStream << (i+1) << "~&~" << meanQ2 << "~&~" << meanW << "~&~" << meanx << "~&~";
    meanVariablesStream << meany << "~&~" << meant; 
    meanVariablesStream << std::string(" \\\\ \\hline ");

    // outputs of mean kinematic variables for plotting
    meanVariablesPlotStream << "{" << meanQ2 << ", " << meanW << ", " << meanx << ", ";
    meanVariablesPlotStream << meany << ", "; meanVariablesPlotStream << meant << "}";
    if (i < numBins - 1) {
        meanVariablesPlotStream << ", "; 
    }
  }

  chi2FitsAStream << "};";  chi2FitsBStream << "};";  chi2FitsCStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << chi2FitsAStream.str() << std::endl;
  outputFile << chi2FitsBStream.str() << std::endl;
  if (asymmetry_index==0 || asymmetry_index==1) { outputFile << chi2FitsCStream.str() << std::endl; }

  outputFile.close();

  meanVariablesStream << "\\end{tabular}\n";
  meanVariablesStream << "\\caption{The mean kinematic variables in each of the bins ";
  meanVariablesStream << "for the extracted $" << prefix << "$ asymmetries.";
  meanVariablesStream << " Values given in GeV or GeV$^2$ where appropriate.}\n";
  meanVariablesStream << "\\end{table}\n";
  meanVariablesStream << endl << endl << endl;
  if (asymmetry_index == 0) {
    std::ofstream kinematicFile(kinematic_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicFile << meanVariablesStream.str() << std::endl; 
    kinematicFile.close();

    meanVariablesPlotStream << "};";
    std::ofstream kinematicPlot_File(kinematicPlot_file, std::ios_base::app);
    // Write the string stream content to the file
    kinematicPlot_File << meanVariablesPlotStream.str() << std::endl;
    kinematicPlot_File.close();
  }
}
