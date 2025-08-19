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
        ALL = (meanDepA/meanDepC)*ALL;
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
void negLogLikelihood_single_hadron(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
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
  double Df = dilution_factor;
  double sigmaDf = dilutionFactors[currentBin].second;
  // // Su22
  // double sigmaPb = 0.0086;
  // double sigmaPtp = 0.0368;
  // double sigmaPtm = 0.03542;
  // // Fa22
  // double sigmaPb = 0.0045;
  // double sigmaPtp = 0.0272;
  // double sigmaPtm = 0.0404;
  // // Sp23
  // double sigmaPb = 0.0061;
  // double sigmaPtp = 0.0404;
  // double sigmaPtm = 0.0376;

  // // Random number generation setup (outside the loop)
  // std::random_device rd;
  // std::mt19937 gen(rd());

  // // Normal distributions (outside the loop)
  // std::normal_distribution<> distDf(0.0, sigmaDf);
  // std::normal_distribution<> distPb(0.0, sigmaPb);
  // std::normal_distribution<> distStandard(0.0, 1.0); // Standard normal distribution

  while (dataReader.Next()) {
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= allBins[currentFits][currentBin] && 
          *currentVariable < allBins[currentFits][currentBin + 1] && passedKinematicCuts) {

      // Increment the event count
      N += 1;

      // // // Get per-event values
      // double Df = dilution_factor + distDf(gen);
      // double Pb = *beam_pol;
      // double Pt = std::abs(*target_pol);
      // // Adjust Pb with its uncertainty
      // Pb += distPb(gen);
      // // Select sigma for Pt based on the sign of *target_pol
      // double sigmaPt = (*target_pol >= 0) ? sigmaPtp : sigmaPtm;
      // // Adjust Pt with its uncertainty
      // Pt += sigmaPt * distStandard(gen);
      // Restore the sign of Pt
      // double signPt = (*target_pol >= 0) ? 1.0 : -1.0;
      // Pt = signPt * Pt;

      // Proceed with your calculations
      double Pb = *beam_pol;                 // moved inside the loop
      double Pt = std::abs(*target_pol);     // moved inside the loop

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
    cout << endl << "Beginning MLM fit for " << binNames[currentFits]
      << " bin " << i << ". ";
    currentBin = i;

    // std::vector<double> chi2Result = chi2Fits[key][currentFits];
    // Define the parameters with initial values and limits
    minuit.DefineParameter(0, "ALU_sinphi", 0.00, 0.01, -1, 1);
    minuit.DefineParameter(1, "AUL_sinphi", 0.00, 0.00, -1, 1);
    minuit.DefineParameter(2, "AUL_sin2phi", 0.00, 0.00, -1, 1);
    minuit.DefineParameter(3, "ALL", 0.00, 0.00, -1, 1);
    minuit.DefineParameter(4, "ALL_cosphi", 0.00, 0.00, -1, 1);
    minuit.DefineParameter(5, "AUU_cosphi", 0.01, 0.01, -1, 1);
    minuit.DefineParameter(6, "AUU_cos2phi", 0.01, 0.01, -1, 1);

    // After defining parameters
    minuit.Migrad(); cout << endl; // First attempt to find the minimum

    // If you decide to use MINImize, replace Migrad with the following lines:
    arglist[0] = 500; // Max calls
    arglist[1] = 1.;  // Tolerance
    minuit.mnexcm("MINImize", arglist, 2, ierflg);


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
    asymmetryStream << "$" << 100*ALU_sinphi << "_{" << TMath::Abs(100*0.022*ALU_sinphi) << "}^{";
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
        if (i == 2) paramName = "#it{A}_{UU}^{cos#phi}";
        if (i == 3) paramName = "#it{A}_{UU}^{cos2#phi}";
      } else if (asymmetryIndex == 1) {
        if (i == 1) paramName = "#it{A}_{UL}^{sin#phi}";
        if (i == 2) paramName = "#it{A}_{UL}^{sin2#phi}";
        if (i == 3) paramName = "#it{A}_{UU}^{cos#phi}";
        if (i == 4) paramName = "#it{A}_{UU}^{cos2#phi}";
      } else if (asymmetryIndex == 2) {
        if (i == 1) paramName = "#it{A}_{LL}^{cos#phi}";
        if (i == 2) paramName = "#it{A}_{UU}^{cos#phi}";
        if (i == 3) paramName = "#it{A}_{UU}^{cos2#phi}";
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
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  while (dataReader.Next()) {
    
    // Apply kinematic cuts (this function will need to be adapted)
    bool passedKinematicCuts = kinematicCuts->applyCuts(currentFits, false);
    // bool passedKinematicCuts = true;
    // Check if the currentVariable is within the desired range
    if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
      sumVariable += *currentVariable;
      if (*helicity > 0 && *target_pol < 0) { histPosNeg->Fill(*phi); } 
      else if (*helicity < 0 && *target_pol > 0) {  histNegPos->Fill(*phi);}

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
    std::cout << Npp << " " << Npm << " " << Nmp << " " << Nmm << std::endl;

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
  std::ostringstream chi2FitsDStream, chi2FitsEStream;
  chi2FitsAStream << std::fixed << std::setprecision(9);
  chi2FitsBStream << std::fixed << std::setprecision(9);
  chi2FitsCStream << std::fixed << std::setprecision(9);
  chi2FitsDStream << std::fixed << std::setprecision(9);
  chi2FitsEStream << std::fixed << std::setprecision(9);

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
      chi2FitsCStream << prefix << "chi2FitsAUUcosphi = {";
      chi2FitsDStream << prefix << "chi2FitsAUUcos2phi = {";
      break;
    case 1: // target-spin asymmetry
      fitFunction = new TF1("fitFunction", TSA_single_hadron, 0, 2*TMath::Pi(), 5);
      chi2FitsAStream << prefix << "chi2FitsAULoffset = {";
      chi2FitsBStream << prefix << "chi2FitsAULsinphi = {";
      chi2FitsCStream << prefix << "chi2FitsAULsin2phi = {";
      chi2FitsDStream << prefix << "chi2FitsAUUcosphi = {";
      chi2FitsEStream << prefix << "chi2FitsAUUcos2phi = {";
      break;
    case 2: // double-spin asymmetry
      fitFunction = new TF1("fitFunction", DSA_single_hadron, 0, 2*TMath::Pi(), 4);
      chi2FitsAStream << prefix << "chi2FitsALL = {";
      chi2FitsBStream << prefix << "chi2FitsALLcosphi = {";
      chi2FitsCStream << prefix << "chi2FitsAUUcosphi = {";
      chi2FitsDStream << prefix << "chi2FitsAUUcos2phi = {";
      break;
    default:
      cout << "Invalid asymmetry_index! Using default function form of BSA." << endl;
      fitFunction = new TF1("fitFunction", BSA_single_hadron, 0, 2*TMath::Pi(), 4);
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
        double AUU_cosphi = fitFunction->GetParameter(2);
        double AUU_cosphi_error = fitFunction->GetParError(2);
        double AUU_cos2phi = fitFunction->GetParameter(3);
        double AUU_cos2phi_error = fitFunction->GetParError(3);
        ALU_sinphi = (meanDepA/meanDepW)*ALU_sinphi;
        ALU_sinphi_error = (meanDepA/meanDepW)*ALU_sinphi_error;
        AUU_cosphi = (meanDepA/meanDepV)*AUU_cosphi;
        AUU_cosphi_error = (meanDepA/meanDepV)*AUU_cosphi_error;
        AUU_cos2phi = (meanDepA/meanDepB)*AUU_cos2phi;
        AUU_cos2phi_error = (meanDepA/meanDepB)*AUU_cos2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALU_offset << ", " << ALU_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALU_sinphi << ", " << ALU_sinphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", " << AUU_cosphi_error <<"}";
        chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", " << AUU_cos2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; 
            chi2FitsCStream << ", "; chi2FitsDStream << ", "; 
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
        double AUU_cosphi = fitFunction->GetParameter(3);
        double AUU_cosphi_error = fitFunction->GetParError(3);
        double AUU_cos2phi = fitFunction->GetParameter(4);
        double AUU_cos2phi_error = fitFunction->GetParError(4);
        AUL_sinphi = (meanDepA/meanDepV)*AUL_sinphi;
        AUL_sinphi_error = (meanDepA/meanDepV)*AUL_sinphi_error;
        AUL_sin2phi = (meanDepA/meanDepB)*AUL_sin2phi;
        AUL_sin2phi_error = (meanDepA/meanDepB)*AUL_sin2phi_error;
        AUU_cosphi = (meanDepA/meanDepV)*AUU_cosphi;
        AUU_cosphi_error = (meanDepA/meanDepV)*AUU_cosphi_error;
        AUU_cos2phi = (meanDepA/meanDepB)*AUU_cos2phi;
        AUU_cos2phi_error = (meanDepA/meanDepB)*AUU_cos2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< AUL_offset << ", " << AUL_offset_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< AUL_sinphi << ", " << AUL_sinphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUL_sin2phi << ", " << AUL_sin2phi_error <<"}";
        chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", " << AUU_cosphi_error <<"}";
        chi2FitsEStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", " << AUU_cos2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", "; chi2FitsCStream << ", ";
            chi2FitsDStream << ", "; chi2FitsEStream << ", ";
        }
        break;
      }
      case 2: {// double-spin asymmetry
        // Get the fitted parameters and their errors
        double ALL = fitFunction->GetParameter(0);
        double ALL_error = fitFunction->GetParError(0);
        double ALL_cosphi = fitFunction->GetParameter(1);
        double ALL_cosphi_error = fitFunction->GetParError(1);
        double AUU_cosphi = fitFunction->GetParameter(2);
        double AUU_cosphi_error = fitFunction->GetParError(2);
        double AUU_cos2phi = fitFunction->GetParameter(3);
        double AUU_cos2phi_error = fitFunction->GetParError(3);
        ALL = (meanDepA/meanDepC)*ALL;
        ALL_error = (meanDepA/meanDepC)*ALL_error;
        ALL_cosphi = (meanDepA/meanDepW)*ALL_cosphi;
        ALL_cosphi_error = (meanDepA/meanDepW)*ALL_cosphi_error;
        AUU_cosphi = (meanDepA/meanDepV)*AUU_cosphi;
        AUU_cosphi_error = (meanDepA/meanDepV)*AUU_cosphi_error;
        AUU_cos2phi = (meanDepA/meanDepB)*AUU_cos2phi;
        AUU_cos2phi_error = (meanDepA/meanDepB)*AUU_cos2phi_error;
        chi2FitsAStream<<"{"<<meanVariable<<", "<< ALL << ", " << ALL_error <<"}";
        chi2FitsBStream<<"{"<<meanVariable<<", "<< ALL_cosphi << ", " << ALL_cosphi_error <<"}";
        chi2FitsCStream<<"{"<<meanVariable<<", "<< AUU_cosphi << ", " << AUU_cosphi_error <<"}";
        chi2FitsDStream<<"{"<<meanVariable<<", "<< AUU_cos2phi << ", " << AUU_cos2phi_error <<"}";
        if (i < numBins - 1) {
            chi2FitsAStream << ", "; chi2FitsBStream << ", ";
            chi2FitsCStream << ", "; chi2FitsDStream << ", ";
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

  chi2FitsAStream << "};";  chi2FitsBStream << "};";  
  chi2FitsCStream << "};"; chi2FitsDStream << "};"; 

  std::ofstream outputFile(output_file, std::ios_base::app);
  outputFile << chi2FitsAStream.str() << std::endl;
  outputFile << chi2FitsBStream.str() << std::endl;
  outputFile << chi2FitsCStream.str() << std::endl;
  outputFile << chi2FitsDStream.str() << std::endl;
  if (asymmetry_index==1) { outputFile << chi2FitsEStream.str() << std::endl; }

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

      // // Random number generation setup
      // std::random_device rd;
      // std::mt19937 gen(rd());

      // // Normal distributions
      // std::normal_distribution<> distDf(0.0, sigmaDf);
      // std::normal_distribution<> distPb(0.0, sigmaPb);

      // // Select sigma for Pt based on the sign of *target_pol
      // double sigmaPt = (*target_pol >= 0) ? sigmaPtp : sigmaPtm;
      // std::normal_distribution<> distPt(0.0, sigmaPt);

      // // Adjust the values
      // Df += distDf(gen);
      // Pb += distPb(gen);
      // Pt += distPt(gen);

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
          - Df*Pb*Pt*(
            (*DepC / *DepA)*ALL + 
            (*DepW / *DepA)*ALL_cosphi1*cos(*phi1) + 
            (*DepW / *DepA)*ALL_cosphi2*cos(*phi2) 
          )
        );
      }
      if (*helicity > 0 && *target_pol < 0) { 
        sum_PM = sum_PM + log(1 +
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
            (*DepC / *DepA)*ALU_sin2Deltaphi*sin(2**phi1 - 2**phi2)
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
        (*mc_DepV / *mc_DepA)*AUU_cosphi1*cos(*mc_phi1) +
        (*mc_DepV / *mc_DepA)*AUU_cosphi2*cos(*mc_phi2) +
        (*mc_DepB / *mc_DepA)*AUU_cos2phi1*cos(2**mc_phi1) +
        (*mc_DepB / *mc_DepA)*AUU_cos2phi2*cos(2**mc_phi2) +
        (*mc_DepB / *mc_DepA)*AUU_cosSumphi*cos(*mc_phi1 + *mc_phi2);
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
      double Npm = histPosNeg->GetBinContent(iBinX, iBinY) / cpm;
      double Nmp = histNegPos->GetBinContent(iBinX, iBinY) / cmp;
      double Nmm = histNegNeg->GetBinContent(iBinX, iBinY) / cmm;
      double asymmetry = asymmetry_value_calculation(meanVariable, prefix, 
        Npp*cpp, Npm*cpm, Nmp*cmp, Nmm*cmm, meanPol, Ptp, Ptm, asymmetry_index);
      double error = asymmetry_error_calculation(meanVariable, prefix, 
        Npp*cpp, Npm*cpm, Nmp*cmp, Nmm*cmm, meanPol, Ptp, Ptm, asymmetry_index);

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
      Npp*cpp, Npm*cpm, Nmp*cmp, Nmm*cmm, meanPol, Ptp, Ptm, asymmetry_index);
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
  meanVariablesStream << "& $<x_B>$ & $<y>$ & $<t>$";
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


/******************** General Exclusive CASE ********************/

/******************** GENERAL EXCLUSIVE (simultaneous ALU/AUL/ALL) ********************/

// ===================== GeneralExclusive (begin) =====================

// Context passed implicitly to the FCN
struct GEContext {
  TH1D* hLU  = nullptr;  // ALU(phi) histogram (data)
  TH1D* hUL  = nullptr;  // AUL(phi) histogram (data)
  TH1D* hLL  = nullptr;  // ALL(phi) histogram (data)
  // Mean depolarization ratios used in the model:
  // rVA = <DepV>/<DepA>, rBA = <DepB>/<DepA>, rWA = <DepW>/<DepA>, rCA = <DepC>/<DepA>
  double rVA = 1.0, rBA = 1.0, rWA = 1.0, rCA = 1.0;
};

// Global ctx used by the Minuit FCN
static GEContext g_ge_ctx;

// FCN: chi2 across ALU, AUL, ALL histograms with shared UU denominator (cos and cos2)
// Parameters (all are structure-function ratios unless noted):
//  p[0] = ALU_offset
//  p[1] = AUL_offset
//  p[2] = F_LU^{sin(phi)} / F_UU
//  p[3] = F_UL^{sin(phi)} / F_UU
//  p[4] = F_UL^{sin(2phi)} / F_UU
//  p[5] = F_LL / F_UU
//  p[6] = F_LL^{cos(phi)} / F_UU
//  p[7] = F_UU^{cos(phi)} / F_UU         (shared UU modulation)
//  p[8] = F_UU^{cos(2phi)} / F_UU        (shared UU modulation)
//  p[9] = F_UU^{sin(phi)} / F_{UU}      (shared acceptance modulation)
//  p[10] = F_UU^{sin(2phi)} / F_{UU}      (shared acceptance modulation)
static void chi2Fcn_GeneralExclusive(Int_t& /*npar*/, Double_t* /*gin*/, Double_t& f,
                                     Double_t* par, Int_t /*iflag*/) {
  const double a0 = par[0], a1 = par[1];
  const double aLU = par[2], aUL1 = par[3], aUL2 = par[4];
  const double aLL = par[5], aLLc = par[6];
  const double aUUc = par[7], aUUc2 = par[8];
  const double aUUs = par[9], aUUs2 = par[10];

  auto denom = [&](double phi) {
    return 1.0 + g_ge_ctx.rVA * aUUc  * std::cos(phi)
               + g_ge_ctx.rBA * aUUc2 * std::cos(2.0*phi)
               + aUUs * std::sin(phi)
               + aUUs2 * std::sin(2.0*phi);
  };

  auto modelALU = [&](double phi) {
    return a0 + (g_ge_ctx.rWA * aLU * std::sin(phi)) / denom(phi);
  };
  auto modelAUL = [&](double phi) {
    return a1 + (g_ge_ctx.rVA * aUL1 * std::sin(phi)
              +  g_ge_ctx.rBA * aUL2 * std::sin(2.0*phi)) / denom(phi);
  };
  auto modelALL = [&](double phi) {
    return (g_ge_ctx.rCA * aLL + g_ge_ctx.rWA * aLLc * std::cos(phi)) / denom(phi);
  };

  auto chi2_from_hist = [&](TH1D* h, auto model) -> double {
    if (!h) return 0.0;
    double c2 = 0.0;
    const int nb = h->GetNbinsX();
    for (int i=1;i<=nb;++i){
      const double y  = h->GetBinContent(i);
      const double ey = h->GetBinError(i);
      if (ey<=0) continue;
      const double x  = h->GetBinCenter(i);
      const double yhat = model(x);
      const double pull = (y - yhat) / ey;
      c2 += pull*pull;
    }
    return c2;
  };

  double chi2 = 0.0;
  chi2 += chi2_from_hist(g_ge_ctx.hLU, modelALU);
  chi2 += chi2_from_hist(g_ge_ctx.hUL, modelAUL);
  chi2 += chi2_from_hist(g_ge_ctx.hLL, modelALL);

  f = chi2;
}


// Build three asymmetry histograms (BSA/TSA/DSA) for one kinematic bin
// Returns { hALU(phi), hAUL(phi), hALL(phi) }
static std::tuple<TH1D*, TH1D*, TH1D*>
createHistogramForBin_GeneralExclusive(const char* histBaseName, int binIndex, const std::string& prefix) {

  // Variable range for this kinematic bin
  const double varMin = allBins[currentFits][binIndex];
  const double varMax = allBins[currentFits][binIndex + 1];

  // 4-yield histograms (per helicity/target state) to build asymmetries per phi-bin
  const int nPhiBins = 12;
  TH1D* ppp = new TH1D(Form("%s_pp", histBaseName), "", nPhiBins, 0, 2*TMath::Pi()); // +beam +target
  TH1D* ppm = new TH1D(Form("%s_pm", histBaseName), "", nPhiBins, 0, 2*TMath::Pi()); // +beam -target
  TH1D* pmp = new TH1D(Form("%s_mp", histBaseName), "", nPhiBins, 0, 2*TMath::Pi()); // -beam +target
  TH1D* pmm = new TH1D(Form("%s_mm", histBaseName), "", nPhiBins, 0, 2*TMath::Pi()); // -beam -target

  // to compute mean polarizations for this bin (used by asymmetry_value/error_calculation)
  double sumPol = 0.0;
  double sumTargetPosPol = 0.0, sumTargetNegPol = 0.0;
  int    nPosT = 0, nNegT = 0;
  double sumVar = 0.0, nEvt = 0.0;

  // Readers
  TTreeReaderValue<int>    runnum(dataReader, "runnum");
  TTreeReaderValue<int>    helicity(dataReader, "helicity");
  TTreeReaderValue<double> target_pol(dataReader, "target_pol");
  TTreeReaderValue<double> beam_pol  (dataReader, "beam_pol");
  TTreeReaderValue<double> phi       (dataReader, "phi");
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  while (dataReader.Next()) {
    if (!kinematicCuts->applyCuts(currentFits, false)) continue;
    if (*currentVariable < varMin || *currentVariable >= varMax) continue;

    sumVar += *currentVariable; nEvt += 1.0;
    sumPol += *beam_pol;
    if (*target_pol > 0) { sumTargetPosPol += *target_pol; ++nPosT; }
    else if (*target_pol < 0) { sumTargetNegPol += *target_pol; ++nNegT; }

    if (*helicity > 0 && *target_pol < 0) { ppm->Fill(*phi); }
    else if (*helicity < 0 && *target_pol > 0) { pmp->Fill(*phi); }

    if (*helicity > 0 && (*target_pol >= 0)) { ppp->Fill(*phi); }
    else if (*helicity < 0 && (*target_pol <= 0)) { pmm->Fill(*phi); }
  }
  dataReader.Restart();

  const double meanVar = (nEvt>0)? sumVar/nEvt : 0.0;
  const double meanPb  = (nEvt>0)? sumPol/nEvt : 0.0;
  const double Ptp     = (nPosT>0)? (sumTargetPosPol/nPosT) : 1.0;
  const double Ptm     = (nNegT>0)? -(sumTargetNegPol/nNegT) : 1.0;

  // Build final asymmetry histograms
  TH1D* hALU = new TH1D(Form("%s_ALU", histBaseName), "", nPhiBins, 0, 2*TMath::Pi());
  TH1D* hAUL = new TH1D(Form("%s_AUL", histBaseName), "", nPhiBins, 0, 2*TMath::Pi());
  TH1D* hALL = new TH1D(Form("%s_ALL", histBaseName), "", nPhiBins, 0, 2*TMath::Pi());

  for (int ib = 1; ib <= nPhiBins; ++ib) {
    // Normalize to accumulated charges (guard 0)
    const double Npp = ppp->GetBinContent(ib) / std::max(cpp, 1.0);
    const double Npm = ppm->GetBinContent(ib) / std::max(cpm, 1.0);
    const double Nmp = pmp->GetBinContent(ib) / std::max(cmp, 1.0);
    const double Nmm = pmm->GetBinContent(ib) / std::max(cmm, 1.0);

    // BSA
    {
      const int asym = 0;
      const double val = asymmetry_value_calculation(meanVar, prefix, Npp, Npm, Nmp, Nmm, meanPb, Ptp, Ptm, asym);
      const double err = asymmetry_error_calculation(meanVar, prefix, Npp*cpp, Npm*cpm, Nmp*cmp, Nmm*cmm, meanPb, Ptp, Ptm, asym);
      hALU->SetBinContent(ib, val); hALU->SetBinError(ib, err);
    }
    // TSA
    {
      const int asym = 1;
      const double val = asymmetry_value_calculation(meanVar, prefix, Npp, Npm, Nmp, Nmm, meanPb, Ptp, Ptm, asym);
      const double err = asymmetry_error_calculation(meanVar, prefix, Npp*cpp, Npm*cpm, Nmp*cmp, Nmm*cmm, meanPb, Ptp, Ptm, asym);
      hAUL->SetBinContent(ib, val); hAUL->SetBinError(ib, err);
    }
    // DSA
    {
      const int asym = 2;
      const double val = asymmetry_value_calculation(meanVar, prefix, Npp, Npm, Nmp, Nmm, meanPb, Ptp, Ptm, asym);
      const double err = asymmetry_error_calculation(meanVar, prefix, Npp*cpp, Npm*cpm, Nmp*cmp, Nmm*cmm, meanPb, Ptp, Ptm, asym);
      hALL->SetBinContent(ib, val); hALL->SetBinError(ib, err);
    }
  }

  delete ppp; delete ppm; delete pmp; delete pmm;
  return {hALU, hAUL, hALL};
}


// Plot the three histograms with model curves that include the shared UU denominator.
// Saves a single 3-panel PNG.  Filename includes run/timestamp suffix derived from output_file.
// Plot the three histograms with model curves that include the shared UU denominator.
// Saves a single 3-panel PNG.  Filename includes run/timestamp suffix derived from output_file.
static void plotHistogramAndFit_GeneralExclusive(
  TH1D* hALU, TH1D* hAUL, TH1D* hALL,
  const double par[11], int binIndex, const std::string& prefix,
  const std::string& runSuffix) {

  // Unpack all 11 fit parameters
  const double a0   = par[0],  a1   = par[1];
  const double aLU  = par[2],  aUL1 = par[3],  aUL2 = par[4];
  const double aLL  = par[5],  aLLc = par[6];
  const double aUUc = par[7],  aUUc2= par[8];
  const double aUUs = par[9],  aUUs2= par[10];

  auto denom = [&](double phi) {
    return 1.0
      + g_ge_ctx.rVA * aUUc  * std::cos(phi)
      + g_ge_ctx.rBA * aUUc2 * std::cos(2.0*phi)
      + aUUs  * std::sin(phi)
      + aUUs2 * std::sin(2.0*phi);
  };

  auto yALU = [&](double phi){
    return a0 + (g_ge_ctx.rWA * aLU * std::sin(phi)) / denom(phi);
  };
  auto yAUL = [&](double phi){
    return a1 + (g_ge_ctx.rVA * aUL1 * std::sin(phi)
               +  g_ge_ctx.rBA * aUL2 * std::sin(2.0*phi)) / denom(phi);
  };
  auto yALL = [&](double phi){
    return (g_ge_ctx.rCA * aLL + g_ge_ctx.rWA * aLLc * std::cos(phi)) / denom(phi);
  };

  TCanvas* c = new TCanvas(Form("cGE_%d",binIndex), "", 1600, 560);
  c->Divide(3,1);

  auto drawOne = [&](int pad, TH1D* h, auto ymodel, const char* ytitle){
    c->cd(pad);
    gPad->SetLeftMargin(0.20);
    gPad->SetRightMargin(0.06);
    gPad->SetBottomMargin(0.16);

    TGraphErrors* gr = new TGraphErrors();
    const int nb = h->GetNbinsX();
    for (int i=1;i<=nb;++i){
      gr->SetPoint(i-1, h->GetBinCenter(i), h->GetBinContent(i));
      gr->SetPointError(i-1, 0.0, h->GetBinError(i));
    }
    gr->SetMarkerStyle(kFullCircle);
    gr->SetMarkerColor(kBlack);
    gr->GetXaxis()->SetTitle("#phi");
    gr->GetYaxis()->SetTitle(ytitle);
    gr->GetXaxis()->SetLimits(0, 2*TMath::Pi());
    gr->GetYaxis()->SetTitleOffset(1.6);
    gr->Draw("AP");

    const int np = 360;
    TGraph* gm = new TGraph(np);
    for (int j=0; j<np; ++j){
      double phi = (2.0*TMath::Pi()) * (j/(double)(np-1));
      gm->SetPoint(j, phi, ymodel(phi));
    }
    gm->SetLineColor(kRed);
    gm->Draw("L same");

    TLegend* L = new TLegend(0.10, 0.64, 0.90, 0.90);
    L->SetBorderSize(1); L->SetFillColor(0); L->SetTextSize(0.030);

    const std::string yt(ytitle);
    if (yt == "A_{LU}") {
      L->AddEntry((TObject*)0, Form("offset: %.6f", a0), "");
      L->AddEntry((TObject*)0, Form("A_{LU}^{sin#phi}: %.6f", aLU), "");
    } else if (yt == "A_{UL}") {
      L->AddEntry((TObject*)0, Form("offset: %.6f", a1), "");
      L->AddEntry((TObject*)0, Form("A_{UL}^{sin#phi}: %.6f", aUL1), "");
      L->AddEntry((TObject*)0, Form("A_{UL}^{sin2#phi}: %.6f", aUL2), "");
    } else { // A_LL
      L->AddEntry((TObject*)0, Form("A_{LL}: %.6f", aLL), "");
      L->AddEntry((TObject*)0, Form("A_{LL}^{cos#phi}: %.6f", aLLc), "");
    }
    // Shared UU terms (these two were the culprits—must use aUUs / aUUs2, not aUUc2)
    L->AddEntry((TObject*)0, Form("A_{UU}^{cos#phi}: %.6f",  aUUc),  "");
    L->AddEntry((TObject*)0, Form("A_{UU}^{cos2#phi}: %.6f", aUUc2), "");
    L->AddEntry((TObject*)0, Form("A_{UU}^{sin#phi}: %.6f",  aUUs),  "");
    L->AddEntry((TObject*)0, Form("A_{UU}^{sin2#phi}: %.6f", aUUs2), "");
    L->Draw("same");
  };

  drawOne(1, hALU, yALU, "A_{LU}");
  drawOne(2, hAUL, yAUL, "A_{UL}");
  drawOne(3, hALL, yALL, "A_{LL}");

  const double vminB = allBins[currentFits][binIndex];
  const double vmaxB = allBins[currentFits][binIndex+1];
  std::ostringstream ttl; ttl<<std::fixed<<std::setprecision(3)
    << vminB << " \\leq " << formatLabelName(prefix) << " < " << vmaxB;
  c->SetTitle(ttl.str().c_str());

  std::string fname = "output/individual_chi2_fits/" + prefix +
                      "_GE_bin_" + std::to_string(binIndex) + "_" + runSuffix + ".png";
  c->SaveAs(fname.c_str());
  delete c;
}


// Simultaneous chi2 fits across ALU/AUL/ALL per bin; writes:
//  - per-parameter arrays to output_file
//  - mean-kinematics LaTeX table to kinematic_file
//  - compact kinematics list to kinematicPlot_file
//  - combined, labeled covariance and correlation files to output/results/ with run+timestamp in the name
void performChi2Fits_GeneralExclusive(const char* output_file,
                                      const char* kinematic_file,
                                      const char* kinematicPlot_file,
                                      const string& prefix) {
  // Prepare output streams for parameter arrays
  std::ostringstream sALUoff, sAULoff, sALU, sAUL, sAUL2, sALL, sALLc;
  std::ostringstream sAUUc, sAUUc2, sAUUs, sAUUs2;
  for (auto* s : {&sALUoff,&sAULoff,&sALU,&sAUL,&sAUL2,&sALL,&sALLc,&sAUUc,&sAUUc2,&sAUUs,&sAUUs2})
    (*s) << std::fixed << std::setprecision(9);

  sALUoff << prefix << "GEchi2FitsALUoffset = {";
  sAULoff << prefix << "GEchi2FitsAULoffset = {";
  sALU    << prefix << "GEchi2FitsALUsinphi = {";
  sAUL    << prefix << "GEchi2FitsAULsinphi = {";
  sAUL2   << prefix << "GEchi2FitsAULsin2phi = {";
  sALL    << prefix << "GEchi2FitsALL = {";
  sALLc   << prefix << "GEchi2FitsALLcosphi = {";
  sAUUc   << prefix << "GEchi2FitsAUUcosphi = {";
  sAUUc2  << prefix << "GEchi2FitsAUUcos2phi = {";
  sAUUs  << prefix << "GEchi2FitsAUUsinphi = {";
  sAUUs2  << prefix << "GEchi2FitsAUUsin2phi = {";

  // Kinematic LaTeX and list
  std::ostringstream kinLatex;
  kinLatex << "\\begin{table}[h]\n\\centering\n"
           << "\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline\n"
           << "Bin & $\\langle Q^2\\rangle$ & $\\langle W\\rangle$ & $\\langle x_B\\rangle$ & "
           << "$\\langle y\\rangle$ & $\\langle t\\rangle$ & $\\langle t_{\\min}\\rangle$ \\\\ \\hline\n";

  std::ostringstream kinList;
  kinList << prefix << "GEKinematics = {";

  // Results directory for matrices
  gSystem->mkdir("output/results", kTRUE);

  // ---- Build combined file names using the same run+timestamp suffix as other outputs ----
  auto deriveSuffixFromOut = [](const char* outPath)->std::string{
    std::string s = outPath ? outPath : "";
    size_t slash = s.find_last_of("/\\");
    std::string base = (slash==std::string::npos)? s : s.substr(slash+1); // filename
    size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0,dot);              // drop extension
    const std::string prefixes[] = {"asymmetries_", "kinematics_", "kinematicPlots_"};
    for (const auto& pfx : prefixes) {
      if (base.rfind(pfx, 0) == 0) { base = base.substr(pfx.size()); break; }
    }
    return base; // "<root>_timeStamp_<ts>"
  };

  const std::string suffix = deriveSuffixFromOut(output_file);
  const std::string covPath  = "output/results/GE_" + prefix + "_cov_"  + suffix + ".txt";
  const std::string corrPath = "output/results/GE_" + prefix + "_corr_" + suffix + ".txt";

  // Parameter names (order used everywhere)
  const int npar = 11;
  const char* names[npar] = {
    "ALU_offset","AUL_offset",
    "F_LU_sin/F_UU","F_UL_sin/F_UU","F_UL_sin2/F_UU",
    "F_LL/F_UU","F_LL_cos/F_UU",
    "F_UU_cos/F_UU","F_UU_cos2/F_UU",
    "F_UU_sin/F_UU","F_UU_sin2/F_UU"
  };

  // Open the combined matrix files (truncate once, then append per bin)
  {
    std::ofstream of(covPath, std::ios::out | std::ios::trunc);
    of << std::setprecision(9);
    of << "# Covariance matrices for GeneralExclusive simultaneous fit\n";
    of << "# Prefix (kinematic variable): " << prefix << "\n";
    of << "# Run+timestamp key: " << suffix << "\n";
    of << "# Parameters (order): ";
    for (int ip=0; ip<npar; ++ip) of << names[ip] << (ip<npar-1 ? ", " : "");
    of << "\n\n";
  }
  {
    std::ofstream of(corrPath, std::ios::out | std::ios::trunc);
    of << std::setprecision(9);
    of << "# Correlation matrices for GeneralExclusive simultaneous fit\n";
    of << "# Prefix (kinematic variable): " << prefix << "\n";
    of << "# Run+timestamp key: " << suffix << "\n";
    of << "# rho_{ij} = cov_{ij} / (sigma_i sigma_j)\n";
    of << "# Parameters (order): ";
    for (int ip=0; ip<npar; ++ip) of << names[ip] << (ip<npar-1 ? ", " : "");
    of << "\n\n";
  }

  const size_t numBins = allBins[currentFits].size() - 1;

  // Loop over kinematic bins
  for (size_t i = 0; i < numBins; ++i) {
    std::cout << "Beginning simultaneous chi2 GE fit for " << binNames[currentFits]
              << " bin " << i << ". " << std::endl;

    // Build asymmetry histograms
    char hname[64]; snprintf(hname, sizeof(hname), "GE_%zu", i);
    TH1D *hALU, *hAUL, *hALL;
    std::tie(hALU, hAUL, hALL) = createHistogramForBin_GeneralExclusive(hname, (int)i, prefix);

    // ---- Compute mean depolarization & mean kinematics (and count events) ----
    double sumQ2=0, sumW=0, sumx=0, sumy=0, sumt=0, sumtmin=0, nEvt=0;
    double sumDepA=0, sumDepB=0, sumDepC=0, sumDepV=0, sumDepW=0, sumVar=0;

    {
      TTreeReaderValue<double> Q2(dataReader,"Q2");
      TTreeReaderValue<double> W (dataReader,"W");
      TTreeReaderValue<double> x (dataReader,"x");
      TTreeReaderValue<double> y (dataReader,"y");
      TTreeReaderValue<double> t (dataReader,"t");
      TTreeReaderValue<double> tmin(dataReader,"tmin");
      TTreeReaderValue<double> DepA(dataReader,"DepA");
      TTreeReaderValue<double> DepB(dataReader,"DepB");
      TTreeReaderValue<double> DepC(dataReader,"DepC");
      TTreeReaderValue<double> DepV(dataReader,"DepV");
      TTreeReaderValue<double> DepW(dataReader,"DepW");
      TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

      const double vmin = allBins[currentFits][i];
      const double vmax = allBins[currentFits][i+1];
      while (dataReader.Next()){
        if (!kinematicCuts->applyCuts(currentFits,false)) continue;
        if (*currentVariable < vmin || *currentVariable >= vmax) continue;
        sumQ2 += *Q2; sumW += *W; sumx += *x; sumy += *y; sumt += *t; sumtmin += *tmin;
        sumDepA += *DepA; sumDepB += *DepB; sumDepC += *DepC; sumDepV += *DepV; sumDepW += *DepW;
        sumVar += *currentVariable; nEvt += 1.0;
      }
      dataReader.Restart();
    }

    // Print number of events in this bin (like other channels)
    std::cout << "Found " << nEvt << " events in this bin." << std::endl;

    const double meanVar  = (nEvt>0)? sumVar/nEvt : 0.0;
    const double meanQ2   = (nEvt>0)? sumQ2/nEvt  : 0.0;
    const double meanW    = (nEvt>0)? sumW/nEvt   : 0.0;
    const double meanx    = (nEvt>0)? sumx/nEvt   : 0.0;
    const double meany    = (nEvt>0)? sumy/nEvt   : 0.0;
    const double meant    = (nEvt>0)? sumt/nEvt   : 0.0;
    const double meantmin = (nEvt>0)? sumtmin/nEvt: 0.0;

    const double depA = (nEvt>0)? (sumDepA/nEvt) : 1.0;
    const double depB = (nEvt>0)? (sumDepB/nEvt) : 1.0;
    const double depC = (nEvt>0)? (sumDepC/nEvt) : 1.0;
    const double depV = (nEvt>0)? (sumDepV/nEvt) : 1.0;
    const double depW = (nEvt>0)? (sumDepW/nEvt) : 1.0;

    // Pass dep ratios to FCN (so fitted amplitudes are directly S.F. ratios)
    g_ge_ctx.hLU = hALU;
    g_ge_ctx.hUL = hAUL;
    g_ge_ctx.hLL = hALL;
    g_ge_ctx.rVA = (depA!=0.0)? (depV/depA) : 1.0;
    g_ge_ctx.rBA = (depA!=0.0)? (depB/depA) : 1.0;
    g_ge_ctx.rWA = (depA!=0.0)? (depW/depA) : 1.0;
    g_ge_ctx.rCA = (depA!=0.0)? (depC/depA) : 1.0;

    // ───────────────────────────────────────────────────────────────────
    // Fit with Minuit (robust sequence + proper error evaluation)
    // ───────────────────────────────────────────────────────────────────

    TMinuit minuit(11);

    // Print level: -1 = quiet, 0 = terse, 1 = normal, 2 = verbose.
    // Keep it quiet for production; bump to 1 while debugging.
    minuit.SetPrintLevel(-1);

    // ERROR DEF = 1.0 means your FCN returns χ² (so Δχ² = 1 ⇒ “1σ”).
    // If your FCN were −log L, you would use 0.5 instead.
    minuit.SetErrorDef(1.0);

    // Your χ² functor:
    minuit.SetFCN(chi2Fcn_GeneralExclusive);

    // -------------------------------------------------------------------
    // Parameter definitions
    // name, initial value, step size, lower bound, upper bound
    //
    // NOTE on step sizes:
    //   - A realistic step size helps MIGRAD find the curvature and speeds up HESSE.
    //   - If a parameter should be *fixed*, you must call FixParameter(i).
    //     Setting step=0 does NOT fix it in TMinuit.
    // -------------------------------------------------------------------
    minuit.DefineParameter(0,  "ALU_offset",      0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(1,  "AUL_offset",      0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(2,  "F_LU_sin/F_UU",   0.00,  0.001, -1.0,  1.0);
    minuit.DefineParameter(3,  "F_UL_sin/F_UU",   0.00,  0.001, -1.0,  1.0);
    minuit.DefineParameter(4,  "F_UL_sin2/F_UU",  0.00,  0.001, -1.0,  1.0);
    minuit.DefineParameter(5,  "F_LL/F_UU",       0.00,  0.001, -1.0,  1.0);
    minuit.DefineParameter(6,  "F_LL_cos/F_UU",   0.00,  0.001, -1.0,  1.0);
    minuit.DefineParameter(7,  "F_UU_cos/F_UU",   0.00,  0.001, -0.5,  0.5);
    minuit.DefineParameter(8,  "F_UU_cos2/F_UU",  0.00,  0.001, -0.5,  0.5);

    // These two are meant to be zero (no sine harmonics in UU); FIX them:
    minuit.DefineParameter(9,  "F_UU_sin/F_UU",   0.00,  0.001, -1.0,  1.0);
    minuit.DefineParameter(10, "F_UU_sin2/F_UU",  0.00,  0.001, -1.0,  1.0);
    minuit.FixParameter(9);
    minuit.FixParameter(10);

    // -------------------------------------------------------------------
    // Robust minimization control
    //  - Strategy 2: more thorough (slower) but improves reliability.
    //  - Tolerance ~ 0.1 is a good default; it’s the EDM target in ERRDEF units.
    //  - Use MINImize (MIGRAD (+ SIMPLEX fallback) in one command).
    // -------------------------------------------------------------------
    double arglist[10];
    int ierflg = 0;

    // Set strategy = 2
    arglist[0] = 2;                      // 0=fast, 1=default, 2=robust
    minuit.mnexcm("SET STR", arglist, 1, ierflg);

    // Set tolerance (target EDM)
    arglist[0] = 0.1;                    // tighten (e.g. 0.05 or 0.01) if needed
    minuit.mnexcm("SET TOL", arglist, 1, ierflg);

    // Run MINImize with generous call budget
    arglist[0] = 5000;                   // max function calls
    arglist[1] = 0.1;                    // same meaning as SET TOL (can omit if already set)
    minuit.mnexcm("MINImize", arglist, 2, ierflg);

    // If MINImize struggled, try SIMPLEX to move closer, then MIGRAD again
    if (ierflg != 0) {
      arglist[0] = 2000;
      minuit.mnexcm("SIMPLEX", arglist, 1, ierflg);

      arglist[0] = 5000; arglist[1] = 0.1;
      minuit.mnexcm("MIGRAD",  arglist, 2, ierflg);
    }

    // -------------------------------------------------------------------
    // Error evaluation
    //  - HESSE computes the covariance matrix at the minimum (parabolic errors).
    //  - MINOS (optional) gives *asymmetric* errors; helpful near bounds
    //    or when uncertainties are non-parabolic. It’s slower—enable as needed.
    // -------------------------------------------------------------------
    minuit.mnexcm("HESSE", nullptr, 0, ierflg);

    // ----- Optional MINOS (asymmetric) errors: physics-only, near-bounds gate -----
    bool doMINOS               = true;   // master switch
    bool minosOnlyPhysics      = true;   // only amplitudes (indices 2..8)
    bool minosNearBoundsOnly   = true;  // run MINOS only if value ~ bound

    if (doMINOS) {
      double arglist[2]; int ierflg = 0;

      auto run_minos_for = [&](int ip) {
        // TMinuit mnexcm("MINOS", ...) expects 1-based parameter number!
        arglist[0] = ip + 1;
        minuit.mnexcm("MINOS", arglist, 1, ierflg);
      };

      // Helper: should we run MINOS for this parameter?
      auto should_run_minos = [&](int ip)->bool {
        if (!minosNearBoundsOnly) return true;

        // Query current value, parabolic error, and bounds
        TString pname; Double_t val=0, err=0, lo=0, up=0; Int_t iv=0;
        minuit.mnpout(ip, pname, val, err, lo, up, iv);
        if (iv == 0) return false;             // not variable (fixed/const)

        bool hasLimits = (lo < up);            // MINUIT uses lo==up when no limits
        if (!hasLimits) return false;          // near-bounds test only makes sense with limits

        // "Near bound" heuristic: within max(2*σ_parab, 10% of range)
        double range  = up - lo;
        double margin = std::max(2.0*err, 0.10*range);
        return ((val - lo) < margin) || ((up - val) < margin);
      };

      if (minosOnlyPhysics) {
        const int physIdx[] = {2,3,4,5,6,7,8}; // F_LU_sin/F_UU ... F_UU_cos2/F_UU
        for (int ip : physIdx) {
          // skip if fixed (shouldn't be), or not near bounds (if gated)
          TString pname; Double_t v,e,lo,up; Int_t iv;
          minuit.mnpout(ip, pname, v, e, lo, up, iv);
          if (iv == 0) continue;               // fixed or not variable
          if (!should_run_minos(ip)) continue;
          run_minos_for(ip);
        }
      } else {
        // All free parameters
        for (int ip = 0; ip < 11; ++ip) {
          TString pname; Double_t v,e,lo,up; Int_t iv;
          minuit.mnpout(ip, pname, v, e, lo, up, iv);
          if (iv == 0) continue;               // fixed or not variable
          if (!should_run_minos(ip)) continue;
          run_minos_for(ip);
        }
      }
    }

    // -------------------------------------------------------------------
    // Convergence diagnostics (ALWAYS check these before trusting errors)
    //  - fmin: best-fit FCN value (χ² here).
    //  - edm: estimated distance to minimum (want ≪ 1).
    //  - istat: 3=good, 2=covariance made pos-def, 1=forced pos-def, 0=not calculated.
    // -------------------------------------------------------------------
    double fmin, edm, errdef;
    int npari, nparx, istat;
    minuit.mnstat(fmin, edm, errdef, npari, nparx, istat);

    // Optionally warn if convergence is marginal
    if (edm > 1e-3*errdef || istat < 2) {
      // Consider: better start values, looser bounds, higher strategy,
      // or tighter tolerance; also inspect correlations/parameterization.
      if (minuit.GetPrintLevel() <= 0) {
        std::cerr << "[WARN] Minuit convergence is marginal: "
                  << "EDM=" << edm << ", istat=" << istat << "\n";
      }
    }

    // -------------------------------------------------------------------
    // Retrieve fit results (values + symmetric/asymmetric errors)
    //  - GetParameter(i, val, err) returns the parabolic (HESSE) error.
    //  - mnerrs(i, eplus, eminus, eparab, gcc) returns MINOS asym. errors
    //    if MINOS ran; eparab is the parabolic error; gcc is global corr. coeff.
    // -------------------------------------------------------------------
    for (int i = 0; i < 11; ++i) {
      double val, err;
      minuit.GetParameter(i, val, err);

      double eplus=0, eminus=0, eparab=0, gcc=0;
      minuit.mnerrs(i, eplus, eminus, eparab, gcc);  // OK even if MINOS skipped

      // Store or print as you like; example:
      // printf("%2d %-18s  val=% .6f  err(parab)=% .6f  MINOS(+/−)=(% .6f,% .6f)  GCC=% .3f\n",
      //        i, minuit.GetParName(i).Data(), val, err, eplus, eminus, gcc);
    }

    // Extract all 11 parameters and errors
    double pval[11], perr[11];
    for (int ip=0; ip<11; ++ip) minuit.GetParameter(ip, pval[ip], perr[ip]);

    // Plot summary (three panels) with model including shared UU denominator
    plotHistogramAndFit_GeneralExclusive(hALU, hAUL, hALL, pval, (int)i, prefix, suffix);

    // Append to parameter arrays (use mean of binning variable as x)
    sALUoff << "{" << meanVar << ", " << pval[0]  << ", " << perr[0]  << "}";
    sAULoff << "{" << meanVar << ", " << pval[1]  << ", " << perr[1]  << "}";
    sALU    << "{" << meanVar << ", " << pval[2]  << ", " << perr[2]  << "}";
    sAUL    << "{" << meanVar << ", " << pval[3]  << ", " << perr[3]  << "}";
    sAUL2   << "{" << meanVar << ", " << pval[4]  << ", " << perr[4]  << "}";
    sALL    << "{" << meanVar << ", " << pval[5]  << ", " << perr[5]  << "}";
    sALLc   << "{" << meanVar << ", " << pval[6]  << ", " << perr[6]  << "}";
    sAUUc   << "{" << meanVar << ", " << pval[7]  << ", " << perr[7]  << "}";
    sAUUc2  << "{" << meanVar << ", " << pval[8]  << ", " << perr[8]  << "}";
    sAUUs   << "{" << meanVar << ", " << pval[9]  << ", " << perr[9]  << "}";
    sAUUs2  << "{" << meanVar << ", " << pval[10] << ", " << perr[10] << "}";

    if (i < numBins - 1) {
      sALUoff << ", "; sAULoff << ", "; sALU << ", "; sAUL << ", "; sAUL2 << ", ";
      sALL << ", "; sALLc << ", "; sAUUc << ", "; sAUUc2 << ", "; sAUUs << ", ";
      sAUUs2 << ", ";
    }

    // Kinematics LaTeX row
    kinLatex << std::fixed << std::setprecision(3)
             << (i+1) << " ~&~ " << meanQ2 << " ~&~ " << meanW << " ~&~ " << meanx
             << " ~&~ " << meany << " ~&~ " << meant << " ~&~ " << meantmin
             << " \\\\ \\hline ";

    // Kinematics list row
    kinList << "{" << meanQ2 << ", " << meanW << ", " << meanx << ", "
            << meany << ", " << meant << ", " << meantmin << "}";
    if (i < numBins - 1) kinList << ", ";

    // ---- Save labeled matrices by APPENDING to the combined files ----
    std::vector<double> cov(npar*npar, 0.0);
    minuit.mnemat(cov.data(), npar);

    std::vector<double> errv(npar);
    for (int ip=0; ip<npar; ++ip) errv[ip] = std::sqrt(std::max(cov[ip*npar+ip], 0.0));

    const double vminB = allBins[currentFits][i];
    const double vmaxB = allBins[currentFits][i+1];

    // Covariance block
    {
      std::ofstream of(covPath, std::ios::out | std::ios::app);
      of << std::setprecision(9);
      of << "## Bin " << i << "  Range: [" << vminB << ", " << vmaxB << ")  Events: " << nEvt << "\n";
      of << std::left << std::setw(22) << "#";
      for (int c=0;c<npar;++c) of << std::setw(22) << names[c];
      of << "\n";
      for (int r=0; r<npar; ++r) {
        of << std::left << std::setw(22) << names[r];
        for (int c=0; c<npar; ++c) {
          of << std::setw(22) << cov[r*npar + c];
        }
        of << "\n";
      }
      of << "\n";
    }

    // Correlation block
    {
      std::ofstream of(corrPath, std::ios::out | std::ios::app);
      of << std::setprecision(9);
      of << "## Bin " << i << "  Range: [" << vminB << ", " << vmaxB << ")  Events: " << nEvt << "\n";
      of << std::left << std::setw(22) << "#";
      for (int c=0;c<npar;++c) of << std::setw(22) << names[c];
      of << "\n";
      for (int r=0; r<npar; ++r) {
        of << std::left << std::setw(22) << names[r];
        for (int c=0; c<npar; ++c) {
          double denom = (errv[r]>0 && errv[c]>0) ? (errv[r]*errv[c]) : 1.0;
          double rho = cov[r*npar + c] / denom;
          of << std::setw(22) << rho;
        }
        of << "\n";
      }
      of << "\n";
    }

    // cleanup hists
    delete hALU; delete hAUL; delete hALL;
  }

  // Close arrays and write to files
  sALUoff << "};"; sAULoff << "};"; sALU << "};"; sAUL << "};";
  sAUL2  << "};"; sALL    << "};"; sALLc<< "};"; sAUUc << "};"; sAUUc2 << "};";
  sAUUs << "};"; sAUUs2 << "};";

  {
    std::ofstream out(output_file, std::ios::app);
    out << sALUoff.str() << "\n";
    out << sAULoff.str() << "\n";
    out << sALU.str()    << "\n";
    out << sAUL.str()    << "\n";
    out << sAUL2.str()   << "\n";
    out << sALL.str()    << "\n";
    out << sALLc.str()   << "\n";
    out << sAUUc.str()   << "\n";
    out << sAUUc2.str()  << "\n";
    out << sAUUs.str()  << "\n";
    out << sAUUs2.str()  << "\n";
  }

  // Finish LaTeX/table and kinematics list
  kinLatex << "\\end{tabular}\n"
           << "\\caption{Mean kinematics per bin for the simultaneous BSA/TSA/DSA "
           << "(GeneralExclusive) fit vs $" << prefix << "$.}\n"
           << "\\label{table:GE_kinematics_" << prefix << "}\n"
           << "\\end{table}\n\n\n";
  {
    std::ofstream kf(kinematic_file, std::ios::app);
    kf << kinLatex.str() << std::endl;
  }

  kinList << "};";
  {
    std::ofstream kp(kinematicPlot_file, std::ios::app);
    kp << kinList.str() << "\n";
  }
}

// ===================== GeneralExclusive (end) =====================