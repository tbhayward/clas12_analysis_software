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
#include "TSpline.h"
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

// ------------------------------------------------------------------
// Global toggles for the MLM fits (defaults are conservative)
// ------------------------------------------------------------------
struct MLMFitOptions {
  bool includeMCinLL        = false;  // turn MC normalization on/off in the NLL
  bool useChargeScaling     = true;  // apply quadrant charge reweighting
  bool runMINOS             = true;  // enable MINOS asymmetric errors
  bool minosOnlyPhysics     = true;  // restrict MINOS to physics amplitudes
  bool minosNearBoundsOnly  = true;  // run MINOS only when param near its bounds
  int  printLevel           = -1;    // TMinuit print level (-1 quiet)
  int  strategy             = 2;     // 0 fast, 1 default, 2 robust
  double tol                = 0.1;   // target EDM (ERRDEF units)
  int  maxCalls             = 5000;  // function call budget for MINImize
};

// Make this visible to both functions
static MLMFitOptions g_mlm_opts;

// Negative log-likelihood function
// Safe log1p helper: guard against tiny/negative arguments to log
inline double SAFE_LOG1P(double x) {
  // x = model - 1  (we compute log(1 + x))
  const double y = 1.0 + x;
  const double eps = 1e-12;
  return std::log(y > eps ? y : eps);
}

// Negative log-likelihood function (single-hadron, unbinned)
void negLogLikelihood_single_hadron(Int_t &npar, Double_t * /*gin*/, Double_t &f,
                                    Double_t *par, Int_t /*iflag*/) {
  // Parameters:
  //  0: ALU_sinphi
  //  1: AUL_sinphi
  //  2: AUL_sin2phi
  //  3: ALL
  //  4: ALL_cosphi
  //  5: AUU_cosphi
  //  6: AUU_cos2phi

  const double ALU_sinphi  = par[0];
  const double AUL_sinphi  = par[1];
  const double AUL_sin2phi = par[2];
  const double ALL         = par[3];
  const double ALL_cosphi  = par[4];
  const double AUU_cosphi  = par[5];
  const double AUU_cos2phi = par[6];

  // Dilution factor for this bin
  const double Df      = dilutionFactors[currentBin].first;
  // const double sigmaDf = dilutionFactors[currentBin].second; // available if you later add nuisance profiling

  // Readers (data)
  TTreeReaderValue<int>    helicity     (dataReader, "helicity");
  TTreeReaderValue<double> beam_pol     (dataReader, "beam_pol");
  TTreeReaderValue<double> target_pol   (dataReader, "target_pol");
  TTreeReaderValue<double> phi          (dataReader, "phi");
  TTreeReaderValue<double> DepA         (dataReader, "DepA");
  TTreeReaderValue<double> DepB         (dataReader, "DepB");
  TTreeReaderValue<double> DepC         (dataReader, "DepC");
  TTreeReaderValue<double> DepV         (dataReader, "DepV");
  TTreeReaderValue<double> DepW         (dataReader, "DepW");
  TTreeReaderValue<double> currentVar   (dataReader, propertyNames[currentFits].c_str());

  // Accumulators
  double N = 0.0;          // number of selected data events
  double sum_PP = 0.0;     // (beam +, target +)
  double sum_PM = 0.0;     // (beam +, target -)
  double sum_MP = 0.0;     // (beam -, target +)
  double sum_MM = 0.0;     // (beam -, target -)

  // Select events and accumulate log-likelihood parts
  while (dataReader.Next()) {
    if (!kinematicCuts->applyCuts(currentFits, false)) continue;
    if (*currentVar < allBins[currentFits][currentBin] ||
        *currentVar >= allBins[currentFits][currentBin + 1]) continue;

    ++N;

    const double Pb = *beam_pol;
    const double Pt_abs = std::abs(*target_pol);
    const int    sB = (*helicity > 0) ? +1 : -1;      // beam helicity sign (+/-)
    const int    sT = (*target_pol >= 0) ? +1 : -1;   // target polarization sign (+/-)
    const double ph = *phi;

    // Depolarization ratios (compute once per event)
    const double rVA = (*DepV)/(*DepA);
    const double rBA = (*DepB)/(*DepA);
    const double rWA = (*DepW)/(*DepA);
    const double rCA = (*DepC)/(*DepA);

    // Model increment excluding the "1": i.e. terms added to 1 in the log(1 + …)
    // UU terms always add, BSA flips with beam helicity, TSA with target sign,
    // DSA with the product.
    const double deltaUU = rVA*AUU_cosphi*std::cos(ph) + rBA*AUU_cos2phi*std::cos(2.0*ph);
    const double deltaB  = sB * Pb * (rWA * ALU_sinphi * std::sin(ph));
    const double deltaT  = sT * Df * Pt_abs * (rVA*AUL_sinphi*std::sin(ph) + rBA*AUL_sin2phi*std::sin(2.0*ph));
    const double deltaD  = (sB*sT) * Df * Pb * Pt_abs * (rCA*ALL + rWA*ALL_cosphi*std::cos(ph));

    const double d = deltaUU + deltaB + deltaT + deltaD;

    // Route to quadrant sum (so we can apply your charge scaling per quadrant later)
    if (sB > 0 && sT > 0)      sum_PP += SAFE_LOG1P(d);
    else if (sB > 0 && sT < 0) sum_PM += SAFE_LOG1P(d);
    else if (sB < 0 && sT > 0) sum_MP += SAFE_LOG1P(d);
    else                       sum_MM += SAFE_LOG1P(d);
  }
  dataReader.Restart();

  // MC normalization integral over UU terms only (as in your original)
  double NUU = 1.0;  // default neutral normalization if MC is disabled or empty
  if (g_mlm_opts.includeMCinLL) {
    double accum = 0.0;
    Long64_t mcount = 0;

    TTreeReaderValue<double> mc_phi   (mcReader, "phi");
    TTreeReaderValue<double> mc_DepA  (mcReader, "DepA");
    TTreeReaderValue<double> mc_DepB  (mcReader, "DepB");
    TTreeReaderValue<double> mc_DepV  (mcReader, "DepV");
    TTreeReaderValue<double> mc_currentVar(mcReader, propertyNames[currentFits].c_str());

    while (mcReader.Next()) {
      if (!mckinematicCuts->applyCuts(currentFits, true)) continue;
      if (*mc_currentVar < allBins[currentFits][currentBin] ||
          *mc_currentVar >= allBins[currentFits][currentBin + 1]) continue;

      const double rVA = (*mc_DepV)/(*mc_DepA);
      const double rBA = (*mc_DepB)/(*mc_DepA);
      const double ph  = *mc_phi;
      // Only UU modulation in normalization
      const double normTerm = 1.0 + rVA*AUU_cosphi*std::cos(ph) + rBA*AUU_cos2phi*std::cos(2.0*ph);

      accum += (normTerm > 0.0 ? normTerm : 0.0);  // avoid negative if fluctuates
      ++mcount;
    }
    mcReader.Restart();

    if (mcount > 0) {
      // Using the mean value over MC as normalization factor
      NUU = accum / static_cast<double>(mcount);
      if (!(NUU > 0.0)) NUU = 1.0;  // safety
    } else {
      // No MC entries in this bin; fall back gracefully
      NUU = 1.0;
    }
  }

  // Quadrant charge scaling (optional)
  double wPP = 1.0, wPM = 1.0, wMP = 1.0, wMM = 1.0;
  if (g_mlm_opts.useChargeScaling) {
    // cpp,cpm,cmp,cmm are assumed global (as in your original code)
    const double bp = (cpp + cpm);  // total + beam charge
    const double bm = (cmp + cmm);  // total - beam charge
    const double tp = (cpp + cmp);  // total + target charge
    const double tm = (cpm + cmm);  // total - target charge

    const double minBeamCharge   = std::min(bp, bm);
    const double minTargetCharge = std::min(tp, tm);

    auto safe_div = [](double a, double b) { return (b != 0.0) ? (a / b) : 0.0; };

    wPP = safe_div(minBeamCharge*minTargetCharge, bp*tp);
    wPM = safe_div(minBeamCharge*minTargetCharge, bp*tm);
    wMP = safe_div(minBeamCharge*minTargetCharge, bm*tp);
    wMM = safe_div(minBeamCharge*minTargetCharge, bm*tm);
  }

  // Build the negative log-likelihood:
  //   NLL = N * log(NUU) - (weighted) sum log-likelihood contributions by quadrant
  const double nll = N * std::log(NUU > 0.0 ? NUU : 1.0)
                   - (wPP*sum_PP + wPM*sum_PM + wMP*sum_MP + wMM*sum_MM);

  // (Optional) quick monitor
  // std::cout << "MLM NLL (bin " << currentBin << "): " << nll << "  N=" << N << "  NUU=" << NUU << std::endl;

  f = nll;  // return to Minuit
}

void performMLMFits_single_hadron(const char* output_file,
                                  const char* kinematic_file,
                                  const std::string& prefix) {
  mlmPrefix = prefix;  // your existing global usage

  const size_t numBins = allBins[currentFits].size() - 1;

  // Streams for outputs
  std::ostringstream mlmA, mlmB, mlmC, mlmD, mlmE, mlmF, mlmG;
  for (auto* s : {&mlmA,&mlmB,&mlmC,&mlmD,&mlmE,&mlmF,&mlmG})
    (*s) << std::fixed << std::setprecision(9);

  mlmA << prefix << "MLMFitsALUsinphi = {";
  mlmB << prefix << "MLMFitsAULsinphi = {";
  mlmC << prefix << "MLMFitsAULsin2phi = {";
  mlmD << prefix << "MLMFitsALL = {";
  mlmE << prefix << "MLMFitsALLcosphi = {";
  mlmF << prefix << "MLMFitsAUUcosphi = {";
  mlmG << prefix << "MLMFitsAUUcos2phi = {";

  // LaTeX table
  std::ostringstream asymTab;
  asymTab << "\\begin{table}[h]\n\\centering\n"
          << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|} \\hline\n"
          << "Bin & $<" << prefix << ">$ & $F_{UU}^{\\cos\\phi}/F_{UU}$ & "
          << "$F_{UU}^{\\cos2\\phi}/F_{UU}$ & $F_{LU}^{\\sin\\phi}/F_{UU}$ & "
          << "$F_{UL}^{\\sin\\phi}/F_{UU}$ & $F_{UL}^{\\sin2\\phi}/F_{UU}$ & "
          << "$F_{LL}/F_{UU}$ & $F_{LL}^{\\cos\\phi}/F_{UU}$ \\\\ \\hline\n";

  // Loop bins
  for (size_t i = 0; i < numBins; ++i) {
    currentBin = static_cast<int>(i);
    std::cout << "\nBeginning MLM fit for " << binNames[currentFits]
              << " bin " << i << "… " << std::endl;

    // Create a fresh Minuit per bin (simplifies state handling)
    TMinuit minuit(7);
    minuit.SetPrintLevel(g_mlm_opts.printLevel);
    minuit.SetErrorDef(0.5); // -log L uses ERRDEF=0.5 for 1σ

    minuit.SetFCN(negLogLikelihood_single_hadron);

    // Parameter definitions: (name, init, step, low, high)
    // Use realistic step sizes (not zero); bounds as in χ² fit for UU cosines.
    minuit.DefineParameter(0, "ALU_sinphi",  0.00, 0.005, -1.0,  1.0);
    minuit.DefineParameter(1, "AUL_sinphi",  0.00, 0.005, -1.0,  1.0);
    minuit.DefineParameter(2, "AUL_sin2phi", 0.00, 0.005, -1.0,  1.0);
    minuit.DefineParameter(3, "ALL",         0.00, 0.005, -1.0,  1.0);
    minuit.DefineParameter(4, "ALL_cosphi",  0.00, 0.005, -1.0,  1.0);
    minuit.DefineParameter(5, "AUU_cosphi",  0.00, 0.003, -0.5,  0.5);
    minuit.DefineParameter(6, "AUU_cos2phi", 0.00, 0.003, -0.5,  0.5);

    // Robust control
    {
      double arglist[10]; int ierflg = 0;
      // Strategy
      arglist[0] = g_mlm_opts.strategy;
      minuit.mnexcm("SET STR", arglist, 1, ierflg);

      // MINImize (MIGRAD+SIMPLEX fallback)
      arglist[0] = g_mlm_opts.maxCalls;
      arglist[1] = g_mlm_opts.tol;       // target EDM (in ERRDEF units)
      minuit.mnexcm("MINImize", arglist, 2, ierflg);

      if (ierflg != 0) {
        // Fallback: SIMPLEX then MIGRAD
        arglist[0] = g_mlm_opts.maxCalls/3.0;
        minuit.mnexcm("SIMPLEX", arglist, 1, ierflg);
        arglist[0] = g_mlm_opts.maxCalls;
        arglist[1] = g_mlm_opts.tol;
        minuit.mnexcm("MIGRAD",  arglist, 2, ierflg);
      }

      // HESSE for covariances
      minuit.mnexcm("HESSE", nullptr, 0, ierflg);

      // Optional: MINOS (asymmetric), physics-only and/or near-bounds
      if (g_mlm_opts.runMINOS) {
        auto near_bounds = [&](int ip)->bool{
          if (!g_mlm_opts.minosNearBoundsOnly) return true;
          TString pname; Double_t v,e,lo,up; Int_t iv;
          minuit.mnpout(ip, pname, v, e, lo, up, iv);
          if (iv == 0) return false;
          const bool hasLim = (lo < up);
          if (!hasLim) return false;
          const double range = up - lo;
          const double margin = std::max(2.0*double(e), 0.10*range);
          return ((v - lo) < margin) || ((up - v) < margin);
        };

        double arg1[2]; int ifl=0;
        if (g_mlm_opts.minosOnlyPhysics) {
          const int physIdx[] = {0,1,2,3,4,5,6};
          for (int ip : physIdx) {
            TString pname; Double_t v,e,lo,up; Int_t iv;
            minuit.mnpout(ip, pname, v, e, lo, up, iv);
            if (iv == 0) continue;
            if (!near_bounds(ip)) continue;
            arg1[0] = ip + 1; // 1-based index
            minuit.mnexcm("MINOS", arg1, 1, ifl);
          }
        } else {
          for (int ip=0; ip<7; ++ip) {
            TString pname; Double_t v,e,lo,up; Int_t iv;
            minuit.mnpout(ip, pname, v, e, lo, up, iv);
            if (iv == 0) continue;
            if (!near_bounds(ip)) continue;
            arg1[0] = ip + 1;
            minuit.mnexcm("MINOS", arg1, 1, ifl);
          }
        }
      }

      // Diagnostics
      double fmin, edm, errdef; int npari, nparx, istat;
      minuit.mnstat(fmin, edm, errdef, npari, nparx, istat);
      if (edm > 1e-3*errdef || istat < 2) {
        if (minuit.SetPrintLevel() <= 0) {
          std::cerr << "[WARN][MLM] Convergence marginal in bin " << i
                    << " : EDM=" << edm << "  istat=" << istat << std::endl;
        }
      }
    }

    // Extract results
    double ALU_sinphi,  e_ALU_sinphi;  minuit.GetParameter(0, ALU_sinphi,  e_ALU_sinphi);
    double AUL_sinphi,  e_AUL_sinphi;  minuit.GetParameter(1, AUL_sinphi,  e_AUL_sinphi);
    double AUL_sin2phi, e_AUL_sin2phi; minuit.GetParameter(2, AUL_sin2phi, e_AUL_sin2phi);
    double ALL,         e_ALL;         minuit.GetParameter(3, ALL,         e_ALL);
    double ALL_cosphi,  e_ALL_cosphi;  minuit.GetParameter(4, ALL_cosphi,  e_ALL_cosphi);
    double AUU_cosphi,  e_AUU_cosphi;  minuit.GetParameter(5, AUU_cosphi,  e_AUU_cosphi);
    double AUU_cos2phi, e_AUU_cos2phi; minuit.GetParameter(6, AUU_cos2phi, e_AUU_cos2phi);

    // Compute mean <variable> for the bin
    double sumV = 0.0; double nV = 0.0;
    {
      TTreeReaderValue<double> v(dataReader, propertyNames[currentFits].c_str());
      while (dataReader.Next()) {
        if (!kinematicCuts->applyCuts(currentFits, false)) continue;
        if (*v < allBins[currentFits][i] || *v >= allBins[currentFits][i+1]) continue;
        sumV += *v; nV += 1.0;
      }
      dataReader.Restart();
    }
    const double meanV = (nV > 0.0) ? (sumV / nV) : 0.0;

    // Append to text arrays
    mlmA << "{" << meanV << ", " << ALU_sinphi  << ", " << e_ALU_sinphi  << "}";
    mlmB << "{" << meanV << ", " << AUL_sinphi  << ", " << e_AUL_sinphi  << "}";
    mlmC << "{" << meanV << ", " << AUL_sin2phi << ", " << e_AUL_sin2phi << "}";
    mlmD << "{" << meanV << ", " << ALL         << ", " << e_ALL         << "}";
    mlmE << "{" << meanV << ", " << ALL_cosphi  << ", " << e_ALL_cosphi  << "}";
    mlmF << "{" << meanV << ", " << AUU_cosphi  << ", " << e_AUU_cosphi  << "}";
    mlmG << "{" << meanV << ", " << AUU_cos2phi << ", " << e_AUU_cos2phi << "}";

    if (i < numBins - 1) {
      mlmA << ", "; mlmB << ", "; mlmC << ", "; mlmD << ", ";
      mlmE << ", "; mlmF << ", "; mlmG << ", ";
    }

    // Table row (percent formatting with simple sys placeholders, as in your code)
    asymTab << std::fixed << std::setprecision(2);
    asymTab << (i+1) << " & " << meanV << " & "
            << "$" << 100*AUU_cosphi  << "_{" << std::abs(100*0.5*AUU_cosphi)  << "}^{" << 100*e_AUU_cosphi  << "}$ & "
            << "$" << 100*AUU_cos2phi << "_{" << std::abs(100*0.5*AUU_cos2phi) << "}^{" << 100*e_AUU_cos2phi << "}$ & "
            << "$" << 100*ALU_sinphi  << "_{" << std::abs(100*0.022*ALU_sinphi) << "}^{" << 100*e_ALU_sinphi << "}$ & "
            << "$" << 100*AUL_sinphi  << "_{" << std::abs(100*0.092*AUL_sinphi) << "}^{" << 100*e_AUL_sinphi << "}$ & "
            << "$" << 100*AUL_sin2phi << "_{" << std::abs(100*0.092*AUL_sin2phi)<< "}^{" << 100*e_AUL_sin2phi << "}$ & "
            << "$" << 100*ALL         << "_{" << std::abs(100*0.097*ALL)        << "}^{" << 100*e_ALL         << "}$ & "
            << "$" << 100*ALL_cosphi  << "_{" << std::abs(100*0.097*ALL_cosphi) << "}^{" << 100*e_ALL_cosphi  << "}$ \\\\ \\hline\n";
  }

  // Close arrays
  mlmA << "};"; mlmB << "};"; mlmC << "};"; mlmD << "};";
  mlmE << "};"; mlmF << "};"; mlmG << "};";

  // Write arrays
  {
    std::ofstream out(output_file, std::ios::app);
    out << mlmA.str() << "\n"
        << mlmB.str() << "\n"
        << mlmC.str() << "\n"
        << mlmD.str() << "\n"
        << mlmE.str() << "\n"
        << mlmF.str() << "\n"
        << mlmG.str() << "\n";
  }

  // Close table and write
  asymTab << "\\end{tabular}\n"
          << "\\caption{Mean kinematics and structure function ratios for " << prefix
          << ". Asymmetries shown as $100 A_{\\pm100\\Delta\\mathrm{sys}}^{\\pm100\\Delta\\mathrm{stat}}$.}\n"
          << "\\label{table:kinematics_" << prefix << "}\n"
          << "\\end{table}\n\n\n";
  {
    std::ofstream kf(kinematic_file, std::ios::app);
    kf << asymTab.str() << std::endl;
  }
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

// ─────────────────────────────────────────────────────────────────────
// Context and globals
// ─────────────────────────────────────────────────────────────────────
struct GEContext {
  TH1D* hLU  = nullptr;
  TH1D* hUL  = nullptr;
  TH1D* hLL  = nullptr;

  // Mean depolarization ratios <DepX>/<DepA> used in the model
  double rVA = 1.0, rBA = 1.0, rWA = 1.0, rCA = 1.0;

  // Per-φ-bin means of sin(theta_gamma) and centered/normalized version
  std::vector<double> sTG_phi_mean;     // size = nPhiBins
  std::vector<double> sTG_phi_centered; // (m_i - mbar)/sigma_w

  // Weighted mean and std (for diagnostics/guarding T-leakage params)
  double sTG_wmean = 0.0;
  double sTG_wstd  = 0.0;

  // Smoothed interpolation of centered ⟨sinθγ⟩(φ)
  TGraph*   sTGc_graph   = nullptr;   // owns the (x,y) points
  TSpline3* sTGc_spline  = nullptr;   // cubic spline over those points

  int    nPhiBins = 12;
  double phiMin   = 0.0;
  double phiMax   = 2.0*TMath::Pi();
};
static GEContext g_ge_ctx;

// Global switches to enable/disable leakage fits (fixed to 0 if disabled)
static bool   g_fit_enable_TUL  = true;   // TSA leakage amplitude A_T-UL
static bool   g_fit_enable_TLL  = true;   // DSA leakage amplitude A_T-LL
static double g_fit_fixed_AT_UL = 0.0;    // value used if fixed
static double g_fit_fixed_AT_LL = 0.0;    // value used if fixed

// ─────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────

// Wrap phi into [phiMin, phiMax)
static inline double GE_wrap_phi(double phi, double phiMin, double phiMax) {
  const double L = phiMax - phiMin; // should be 2π
  double p = phi;
  while (p <  phiMin) p += L;
  while (p >= phiMax) p -= L;
  return p;
}

// Fallback: centered value by φ-bin (piecewise-constant)
static inline double GE_sTG_centered_bin(double phi, TH1D* hRef) {
  if (!hRef || g_ge_ctx.sTG_phi_centered.empty()) return 0.0;
  const int ib  = hRef->GetXaxis()->FindBin(phi); // 1..n
  const int idx = std::max(1, std::min(ib, (int)g_ge_ctx.sTG_phi_centered.size())) - 1;
  return g_ge_ctx.sTG_phi_centered[idx];
}

// Preferred: smoothed centered ⟨sinθγ⟩ via spline; falls back to binned
static inline double GE_sTG_centered_interp(double phi, TH1D* hRef) {
  const double pw = GE_wrap_phi(phi, g_ge_ctx.phiMin, g_ge_ctx.phiMax);
  if (g_ge_ctx.sTGc_spline) return g_ge_ctx.sTGc_spline->Eval(pw);
  return GE_sTG_centered_bin(pw, hRef);
}

// ─────────────────────────────────────────────────────────────────────
// FCN: global χ² across ALU, AUL, ALL with shared UU denominator (cos, cos2)
// Parameters (structure-function ratios unless noted):
//  p[0]  = ALU_offset
//  p[1]  = AUL_offset
//  p[2]  = F_LU^{sinφ} / F_UU
//  p[3]  = F_UL^{sinφ} / F_UU
//  p[4]  = F_UL^{sin2φ} / F_UU
//  p[5]  = F_LL / F_UU
//  p[6]  = F_LL^{cosφ} / F_UU
//  p[7]  = F_UU^{cosφ} / F_UU         (shared UU modulation)
//  p[8]  = F_UU^{cos2φ} / F_UU        (shared UU modulation)
//  p[9]  = A_T-UL                     (T-leakage in TSA; multiplies m^c(φ)·sinφ; no depol)
//  p[10] = A_T-LL                     (T-leakage in DSA; multiplies m^c(φ)·sinφ; no depol)
// ─────────────────────────────────────────────────────────────────────
// FCN: χ² with (i) denominator floor barrier and
//      (ii) amplitude box on rVA*aUUc and rVA*aUL1 so their
//           physical asymmetries stay within [-A_MAX_AMP, +A_MAX_AMP].
static void chi2Fcn_GeneralExclusive(Int_t& /*npar*/, Double_t* /*gin*/, Double_t& f,
                                     Double_t* par, Int_t /*iflag*/)
{
  // ─── knobs you can adjust ─────────────────────────────────────────
  const double EPS_FLOOR   = 1e-2;  // demand D(φ) ≥ EPS_FLOOR on [0,2π)
  const double LAMBDA_DEN  = 1e6;   // penalty strength for violating D floor
  const double EPS_EVAL    = 1e-6;  // tiny clamp during evaluation

  const double A_MAX_AMP   = 0.999; // |A| ≤ this (applied to selected amplitudes)
  const bool   LIMIT_A_UUCOS = true;  // enforce on A_UU^cos = rVA*aUUc
  const bool   LIMIT_A_ULSIN = true;  // enforce on A_UL^sin = rVA*aUL1
  const bool   LIMIT_A_ULSIN2= false; // optional
  const bool   LIMIT_A_LU    = false; // optional
  const bool   LIMIT_A_LLCOS = false; // optional
  const bool   LIMIT_A_LL0   = false; // optional
  const double LAMBDA_AMP  = 1e5;   // penalty weight for amplitude box
  // ──────────────────────────────────────────────────────────────────

  // parameters (ratios unless noted)
  const double a0    = par[0],  a1    = par[1];
  const double aLU   = par[2],  aUL1  = par[3],  aUL2 = par[4];
  const double aLL   = par[5],  aLLc  = par[6];
  const double aUUc  = par[7],  aUUc2 = par[8];
  const double aTUL  = par[9];      // A_T-UL (no depol)
  const double aTLL  = par[10];     // A_T-LL (no depol)

  // effective UU coefficients in the denominator (include depol ratios)
  const double B_eff = g_ge_ctx.rVA * aUUc;
  const double C_eff = g_ge_ctx.rBA * aUUc2;

  // denominator and models
  auto denom = [&](double phi) {
    const double d = 1.0 + B_eff*std::cos(phi) + C_eff*std::cos(2.0*phi);
    return (d < EPS_EVAL ? EPS_EVAL : d);
  };
  auto modelALU = [&](double phi, TH1D* /*hRef*/) {
    return a0 + (g_ge_ctx.rWA * aLU * std::sin(phi)) / denom(phi);
  };
  auto modelAUL = [&](double phi, TH1D* hRef) {
    const double sTGc = GE_sTG_centered_interp(phi, hRef);  // centered/σ ⟨sinθγ⟩
    const double num = g_ge_ctx.rVA * aUL1 * std::sin(phi)
                     + g_ge_ctx.rBA * aUL2 * std::sin(2.0*phi)
                     + aTUL * sTGc * std::sin(phi);          // leakage (no depol)
    return a1 + num / denom(phi);
  };
  auto modelALL = [&](double phi, TH1D* hRef) {
    const double sTGc = GE_sTG_centered_interp(phi, hRef);  // reuse same m^c(φ)
    const double num = (g_ge_ctx.rCA * aLL)
                     + (g_ge_ctx.rWA * aLLc * std::cos(phi))
                     + (aTLL * sTGc * std::sin(phi));        // DSA leakage (no depol)
    return num / denom(phi);
  };

  // χ² from one histogram
  auto chi2_from_hist = [&](TH1D* h, auto model) -> double {
    if (!h) return 0.0;
    double c2 = 0.0;
    const int nb = h->GetNbinsX();
    for (int i=1;i<=nb;++i){
      const double y  = h->GetBinContent(i);
      const double ey = h->GetBinError(i);
      if (ey<=0 || !std::isfinite(y) || !std::isfinite(ey)) continue;
      const double phi  = h->GetBinCenter(i);
      const double yhat = model(phi, h);
      const double pull = (y - yhat) / ey;
      c2 += pull*pull;
    }
    return c2;
  };

  // base χ²
  double chi2 = 0.0;
  chi2 += chi2_from_hist(g_ge_ctx.hLU, modelALU);
  chi2 += chi2_from_hist(g_ge_ctx.hUL, modelAUL);
  chi2 += chi2_from_hist(g_ge_ctx.hLL, modelALL);

  // (1) denominator soft barrier over full domain
  auto Dmin_over_domain = [&](double B, double C) -> double {
    auto quad = [&](double x){ return 2.0*C*x*x + B*x + (1.0 - C); }; // x = cosφ ∈ [-1,1]
    double m = std::min(quad(-1.0), quad(+1.0));                      // 1 + C ∓ B
    if (C > 0.0) {
      const double x0 = -B/(4.0*C);
      if (x0 >= -1.0 && x0 <= 1.0) {
        const double fv = 1.0 - C - (B*B)/(8.0*C);
        m = std::min(m, fv);
      }
    }
    return m;
  };
  const double Dmin = Dmin_over_domain(B_eff, C_eff);
  double pen_den = 0.0;
  if (Dmin < EPS_FLOOR) {
    const double deficit = EPS_FLOOR - Dmin;
    pen_den = LAMBDA_DEN * deficit * deficit;
  }

  // (2) amplitude box on selected physical amplitudes
  auto over2 = [&](double A){ const double v = std::fabs(A) - A_MAX_AMP; return (v>0.0)? v*v : 0.0; };

  double pen_amp = 0.0;
  if (LIMIT_A_UUCOS) pen_amp += over2(g_ge_ctx.rVA * aUUc);
  if (LIMIT_A_ULSIN) pen_amp += over2(g_ge_ctx.rVA * aUL1);
  if (LIMIT_A_ULSIN2)pen_amp += over2(g_ge_ctx.rBA * aUL2);
  if (LIMIT_A_LU)    pen_amp += over2(g_ge_ctx.rWA * aLU);
  if (LIMIT_A_LLCOS) pen_amp += over2(g_ge_ctx.rWA * aLLc);
  if (LIMIT_A_LL0)   pen_amp += over2(g_ge_ctx.rCA * aLL);
  pen_amp *= LAMBDA_AMP;

  // final objective
  f = chi2 + pen_den + pen_amp;
}

// ─────────────────────────────────────────────────────────────────────
// Build three asymmetry histograms (BSA/TSA/DSA) for one kinematic bin
// and also compute ⟨sinθγ⟩ per φ-bin for that bin.
// Returns { hALU, hAUL, hALL } and updates g_ge_ctx.{sTG_phi_mean, sTG_phi_centered}
// and the smoothed spline of m^c(φ).
// ─────────────────────────────────────────────────────────────────────
static std::tuple<TH1D*, TH1D*, TH1D*>
createHistogramForBin_GeneralExclusive(const char* histBaseName, int binIndex, const std::string& prefix) {

  // Variable range for this kinematic bin
  const double varMin = allBins[currentFits][binIndex];
  const double varMax = allBins[currentFits][binIndex + 1];

  // φ-binning
  const int nPhiBins = 9;
  const double phiMin = 0.0;
  const double phiMax = 2.0*TMath::Pi();

  TH1D* ppp = new TH1D(Form("%s_pp", histBaseName), "", nPhiBins, phiMin, phiMax);
  TH1D* ppm = new TH1D(Form("%s_pm", histBaseName), "", nPhiBins, phiMin, phiMax);
  TH1D* pmp = new TH1D(Form("%s_mp", histBaseName), "", nPhiBins, phiMin, phiMax);
  TH1D* pmm = new TH1D(Form("%s_mm", histBaseName), "", nPhiBins, phiMin, phiMax);

  // per-φ-bin accumulators for ⟨sinθγ⟩ and counts
  std::vector<double> sTG_sum(nPhiBins, 0.0);
  std::vector<int>    sTG_cnt(nPhiBins, 0);

  // other accumulators
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
  TTreeReaderValue<double> sTG       (dataReader, "sinthetagamma"); // raw sinθγ per event
  TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

  while (dataReader.Next()) {
    if (!kinematicCuts->applyCuts(currentFits, false)) continue;
    if (*currentVariable < varMin || *currentVariable >= varMax) continue;

    sumVar += *currentVariable; nEvt += 1.0;
    sumPol += *beam_pol;
    if (*target_pol > 0) { sumTargetPosPol += *target_pol; ++nPosT; }
    else if (*target_pol < 0) { sumTargetNegPol += *target_pol; ++nNegT; }

    // fill raw yields (by helicity/target)
    if (*helicity > 0 && *target_pol < 0) { ppm->Fill(*phi); }
    else if (*helicity < 0 && *target_pol > 0) { pmp->Fill(*phi); }
    if (*helicity > 0 && (*target_pol >= 0)) { ppp->Fill(*phi); }
    else if (*helicity < 0 && (*target_pol <= 0)) { pmm->Fill(*phi); }

    // accumulate per-φ-bin sinθγ (unweighted average across accepted events)
    int ib = std::max(1, std::min(nPhiBins, (int)std::floor(((*phi - phiMin)/(phiMax-phiMin))*nPhiBins) + 1));
    sTG_sum[ib-1] += *sTG;
    sTG_cnt[ib-1] += 1;
  }
  dataReader.Restart();

  const double meanVar = (nEvt>0)? sumVar/nEvt : 0.0;
  const double meanPb  = (nEvt>0)? sumPol/nEvt : 0.0;
  const double Ptp     = (nPosT>0)? (sumTargetPosPol/nPosT) : 1.0;
  const double Ptm     = (nNegT>0)? -(sumTargetNegPol/nNegT) : 1.0;

  // final asymmetry histograms
  TH1D* hALU = new TH1D(Form("%s_ALU", histBaseName), "", nPhiBins, phiMin, phiMax);
  TH1D* hAUL = new TH1D(Form("%s_AUL", histBaseName), "", nPhiBins, phiMin, phiMax);
  TH1D* hALL = new TH1D(Form("%s_ALL", histBaseName), "", nPhiBins, phiMin, phiMax);

  // per-φ weights from total charge-normalized yield (sum of four samples)
  std::vector<double> w_tot(nPhiBins, 0.0);

  for (int ib = 1; ib <= nPhiBins; ++ib) {
    // normalize to charges (guard 0)
    std::cout << "HELLO WORLD" << std::endl;
    std::cout << ppp->GetBinContent(ib) << " " << ppm->GetBinContent(ib) << " " << pmp->GetBinContent(ib) << " " << pmm->GetBinContent(ib) << std::endl;
    std::cout << "HELLO WORLD" << std::endl;
    const double Npp = ppp->GetBinContent(ib) / std::max(cpp, 1.0);
    const double Npm = ppm->GetBinContent(ib) / std::max(cpm, 1.0);
    const double Nmp = pmp->GetBinContent(ib) / std::max(cmp, 1.0);
    const double Nmm = pmm->GetBinContent(ib) / std::max(cmm, 1.0);


    // weight for centering/normalizing ⟨sinθγ⟩
    w_tot[ib-1] = std::max(0.0, Npp + Npm + Nmp + Nmm);

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

  // Store per-φ-bin ⟨sinθγ⟩ into context (means)
  g_ge_ctx.sTG_phi_mean.assign(nPhiBins, 0.0);
  for (int i=0;i<nPhiBins;++i) {
    g_ge_ctx.sTG_phi_mean[i] = (sTG_cnt[i]>0) ? (sTG_sum[i]/sTG_cnt[i]) : 0.0;
  }

  // Compute charge-weighted φ-average ⟨sTG⟩_w and weighted std σ_w
  double wsum = 0.0, ssum = 0.0;
  for (int i = 0; i < nPhiBins; ++i) {
    const double w = std::max(0.0, w_tot[i]);
    wsum += w;
    ssum += w * g_ge_ctx.sTG_phi_mean[i];
  }
  g_ge_ctx.sTG_wmean = (wsum>0.0) ? (ssum/wsum) : 0.0;

  double vsum = 0.0;
  for (int i = 0; i < nPhiBins; ++i) {
    const double w   = std::max(0.0, w_tot[i]);
    const double dif = g_ge_ctx.sTG_phi_mean[i] - g_ge_ctx.sTG_wmean;
    vsum += w * dif * dif;
  }
  g_ge_ctx.sTG_wstd = (wsum>0.0) ? std::sqrt(vsum/wsum) : 0.0;

  // Center & normalize (if std ~ 0, zero it)
  g_ge_ctx.sTG_phi_centered.assign(nPhiBins, 0.0);
  const double eps_std = 1e-4;
  if (g_ge_ctx.sTG_wstd > eps_std) {
    for (int i=0; i<nPhiBins; ++i) {
      g_ge_ctx.sTG_phi_centered[i] =
        (g_ge_ctx.sTG_phi_mean[i] - g_ge_ctx.sTG_wmean) / g_ge_ctx.sTG_wstd;
    }
  } else {
    for (int i=0; i<nPhiBins; ++i) g_ge_ctx.sTG_phi_centered[i] = 0.0;
  }

  // Build smooth cubic spline for m^c(φ) to avoid jagged model curves
  {
    if (g_ge_ctx.sTGc_spline) { delete g_ge_ctx.sTGc_spline; g_ge_ctx.sTGc_spline = nullptr; }
    if (g_ge_ctx.sTGc_graph)  { delete g_ge_ctx.sTGc_graph;  g_ge_ctx.sTGc_graph  = nullptr; }

    const int    nb   = nPhiBins;
    const double dphi = (phiMax - phiMin) / nb;

    // Wrap with two extra points (before/after) to stabilize edges
    const int npts = nb + 2;
    g_ge_ctx.sTGc_graph = new TGraph(npts);

    // point 0: wrap from last bin
    g_ge_ctx.sTGc_graph->SetPoint(0, phiMin - 0.5*dphi, g_ge_ctx.sTG_phi_centered.back());

    // interior bin centers
    for (int i=0; i<nb; ++i) {
      const double xc = phiMin + (i + 0.5) * dphi;
      const double yc = g_ge_ctx.sTG_phi_centered[i];
      g_ge_ctx.sTGc_graph->SetPoint(i+1, xc, yc);
    }

    // last point: wrap from first bin
    g_ge_ctx.sTGc_graph->SetPoint(npts-1, phiMax + 0.5*dphi, g_ge_ctx.sTG_phi_centered.front());

    // cubic spline
    g_ge_ctx.sTGc_spline = new TSpline3("sTGc_spline", g_ge_ctx.sTGc_graph, "");
  }

  g_ge_ctx.nPhiBins = nPhiBins;
  g_ge_ctx.phiMin   = phiMin;
  g_ge_ctx.phiMax   = phiMax;

  delete ppp; delete ppm; delete pmp; delete pmm;
  return {hALU, hAUL, hALL};
}

// ─────────────────────────────────────────────────────────────────────
// Plot 1×3 canvas with model and GLOBAL χ²/ndf, with uncertainties in legend
// ─────────────────────────────────────────────────────────────────────
// true  -> compact legends (hide chi2/ndf, AUU entries, and A_T-UL/A_T-LL) and use shorter box
// false -> full legends
static bool g_ge_compact_legend = true;

// ─────────────────────────────────────────────────────────────────────

static void plotHistogramAndFit_GeneralExclusive(
  TH1D* hALU, TH1D* hAUL, TH1D* hALL,
  const double par[],                 // values (structure-function ratios / leakages)
  const double err[],                 // uncertainties in same order
  int binIndex, const std::string& prefix,
  const std::string& runSuffix,
  double globalChi2, int globalNdf) {
  const double a0    = par[0],  a1    = par[1];
  const double aLU   = par[2],  aUL1  = par[3],  aUL2 = par[4];
  const double aLL   = par[5],  aLLc  = par[6];
  const double aUUc  = par[7],  aUUc2 = par[8];
  const double aTUL  = par[9];
  const double aTLL  = par[10];

  const double eLU   = err[2],  eUL1 = err[3],  eUL2 = err[4];
  const double eLL   = err[5],  eLLc = err[6];
  const double eUUc  = err[7],  eUUc2= err[8];
  const double eTUL  = err[9];
  const double eTLL  = err[10];

  // ---------------- Model (we plot asymmetries) ----------------
  auto denom = [&](double phi) {
    return 1.0
         + g_ge_ctx.rVA * aUUc  * std::cos(phi)
         + g_ge_ctx.rBA * aUUc2 * std::cos(2.0*phi);
  };
  auto yALU = [&](double phi){
    return a0 + (g_ge_ctx.rWA * aLU * std::sin(phi)) / denom(phi);
  };
  auto yAUL = [&](double phi){
    const double sTGc = GE_sTG_centered_interp(phi, hAUL); // smoothed centered ⟨sinθγ⟩
    const double num = g_ge_ctx.rVA * aUL1 * std::sin(phi)
                     + g_ge_ctx.rBA * aUL2 * std::sin(2.0*phi)
                     + aTUL * sTGc * std::sin(phi);
    return a1 + num / denom(phi);
  };
  auto yALL = [&](double phi){
    const double sTGc = GE_sTG_centered_interp(phi, hALL);
    const double num = (g_ge_ctx.rCA * aLL)
                     + (g_ge_ctx.rWA * aLLc * std::cos(phi))
                     + (aTLL * sTGc * std::sin(phi));
    return num / denom(phi);
  };

  // ---------------- Convert to asymmetry amplitudes for legend (3 d.p.) -------
  auto r3 = [](double x){
    double r = std::round(x * 1000.0) / 1000.0;
    if (std::fabs(r) < 0.0005) r = 0.0; // suppress -0.000
    return r;
  };
  const double A_LU   = r3(g_ge_ctx.rWA * aLU);    const double dA_LU   = r3(g_ge_ctx.rWA * eLU);
  const double A_UL1  = r3(g_ge_ctx.rVA * aUL1);   const double dA_UL1  = r3(g_ge_ctx.rVA * eUL1);
  const double A_UL2  = r3(g_ge_ctx.rBA * aUL2);   const double dA_UL2  = r3(g_ge_ctx.rBA * eUL2);
  const double A_LL   = r3(g_ge_ctx.rCA * aLL);    const double dA_LL   = r3(g_ge_ctx.rCA * eLL);
  const double A_LLc  = r3(g_ge_ctx.rWA * aLLc);   const double dA_LLc  = r3(g_ge_ctx.rWA * eLLc);
  const double A_UUc  = r3(g_ge_ctx.rVA * aUUc);   const double dA_UUc  = r3(g_ge_ctx.rVA * eUUc);
  const double A_UUc2 = r3(g_ge_ctx.rBA * aUUc2);  const double dA_UUc2 = r3(g_ge_ctx.rBA * eUUc2);
  const double A_TUL  = r3(aTUL);                  const double dA_TUL  = r3(eTUL);  // (no depol)
  const double A_TLL  = r3(aTLL);                  const double dA_TLL  = r3(eTLL);  // (no depol)

  // ---------------- Canvas / layout (reduced gaps) ----------------------------
  TCanvas* c = new TCanvas(Form("cGE_%d",binIndex), "", 1600, 560);
  c->Divide(3, 1, 0.002, 0.002);

  // ---------------- Helper to draw each panel ---------------------------------
  const auto addPointsAndCurve = [&](int pad, TH1D* h, auto ymodel,
                                     const char* ytitle,
                                     std::function<void(TLegend*)> fillLegend,
                                     double ylow, double yhigh)
  {
    c->cd(pad);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.06);

    // Data points → TGraphErrors
    TGraphErrors* gr = new TGraphErrors();
    const int nb = h ? h->GetNbinsX() : 0;
    int ip = 0;
    for (int i=1;i<=nb;++i){
      const double x  = h->GetBinCenter(i);
      const double y  = h->GetBinContent(i);
      const double ey = h->GetBinError(i);
      if (!std::isfinite(y) || !std::isfinite(ey)) continue;
      gr->SetPoint(ip, x, y);
      gr->SetPointError(ip, 0.0, ey);
      ++ip;
    }
    gr->SetMarkerStyle(kFullCircle);
    gr->SetMarkerSize(1.1);
    gr->SetMarkerColor(kBlack);
    gr->SetLineColor(kBlack);

    gr->GetXaxis()->SetTitle("#phi");
    gr->GetYaxis()->SetTitle(ytitle);
    gr->GetXaxis()->SetLimits(0, 2*TMath::Pi());
    gr->GetYaxis()->SetRangeUser(ylow, yhigh);

    gr->GetXaxis()->CenterTitle(true);
    gr->GetYaxis()->CenterTitle(true);
    gr->GetXaxis()->SetTitleSize(0.052);
    gr->GetYaxis()->SetTitleSize(0.052);
    gr->GetXaxis()->SetLabelSize(0.044);
    gr->GetYaxis()->SetLabelSize(0.044);
    gr->GetYaxis()->SetTitleOffset(1.85);
    gr->GetXaxis()->SetTitleOffset(1.25);

    gr->Draw("AP");

    // Smooth model curve
    const int np = 1440;
    TGraph* gm = new TGraph(np);
    for (int j=0; j<np; ++j){
      const double phi = (2.0*TMath::Pi()) * (j/(double)(np-1));
      gm->SetPoint(j, phi, ymodel(phi));
    }
    gm->SetLineColor(kRed);
    gm->SetLineWidth(2);
    gm->Draw("L same");

    // Legend
    TLegend* L = new TLegend(0.56, 0.72, 0.95, 0.94 );
    L->SetBorderSize(1);
    L->SetLineColor(kBlack);
    L->SetFillColor(kWhite);
    L->SetFillStyle(1001);
    L->SetTextSize(0.026);
    L->SetTextAlign(12);
    L->SetMargin(0.10);

    if (!g_ge_compact_legend) {
      L->AddEntry((TObject*)0,
        Form("#chi^{2}/ndf (global) = %.1f/%d = %.2f",
             globalChi2, globalNdf,
             (globalNdf>0?globalChi2/globalNdf:0.0)), "");
    }
    fillLegend(L);
    L->Draw("same");
  };

  // ---------------- Legend fillers --------------------------------------------
  auto fillBSA = [&](TLegend* L){
    L->AddEntry((TObject*)0, Form("A_{LU}^{sin#phi} = %.3f #pm %.3f", A_LU,  dA_LU ), "");
    if (!g_ge_compact_legend) {
      L->AddEntry((TObject*)0, Form("A_{UU}^{cos#phi}  = %.3f #pm %.3f", A_UUc,  dA_UUc), "");
      L->AddEntry((TObject*)0, Form("A_{UU}^{cos2#phi} = %.3f #pm %.3f", A_UUc2, dA_UUc2), "");
    }
  };
  auto fillTSA = [&](TLegend* L){
    L->AddEntry((TObject*)0, Form("A_{UL}^{sin#phi}  = %.3f #pm %.3f", A_UL1, dA_UL1), "");
    L->AddEntry((TObject*)0, Form("A_{UL}^{sin2#phi} = %.3f #pm %.3f", A_UL2, dA_UL2), "");
    if (!g_ge_compact_legend) {
      L->AddEntry((TObject*)0, Form("A_{T-UL}^{sin#phi}   = %.3f #pm %.3f", A_TUL, dA_TUL), "");
      L->AddEntry((TObject*)0, Form("A_{UU}^{cos#phi}  = %.3f #pm %.3f", A_UUc,  dA_UUc), "");
      L->AddEntry((TObject*)0, Form("A_{UU}^{cos2#phi} = %.3f #pm %.3f", A_UUc2, dA_UUc2), "");
    }
  };
  auto fillDSA = [&](TLegend* L){
    L->AddEntry((TObject*)0, Form("A_{LL} = %.3f #pm %.3f", A_LL,  dA_LL ), "");
    L->AddEntry((TObject*)0, Form("A_{LL}^{cos#phi} = %.3f #pm %.3f", A_LLc, dA_LLc), "");
    if (!g_ge_compact_legend) {
      L->AddEntry((TObject*)0, Form("A_{T-LL}^{sin#phi}   = %.3f #pm %.3f", A_TLL, dA_TLL), "");
      L->AddEntry((TObject*)0, Form("A_{UU}^{cos#phi}  = %.3f #pm %.3f", A_UUc,  dA_UUc), "");
      L->AddEntry((TObject*)0, Form("A_{UU}^{cos2#phi} = %.3f #pm %.3f", A_UUc2, dA_UUc2), "");
    }
  };

  // ---------------- Draw the three panels -------------------------------------
  addPointsAndCurve(1, hALU, yALU, "A_{LU}", fillBSA, -0.20, 0.20);
  addPointsAndCurve(2, hAUL, yAUL, "A_{UL}", fillTSA, -0.20, 0.20);
  addPointsAndCurve(3, hALL, yALL, "A_{LL}", fillDSA, -0.10, 0.60);

  // ---------------- Title and save --------------------------------------------
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

// ─────────────────────────────────────────────────────────────────────
// Driver: simultaneous fits per bin (11 parameters)
// Writes arrays, covariance/correlation, and LaTeX table with
// [min, mean, max] for Q2, x, y, z, t′
//   - LaTeX table values rounded to two decimals
//   - Column headers are just variable names (no units)
//   - Caption specifies that Q2 and t′ are in GeV^{2}
// ─────────────────────────────────────────────────────────────────────
// ─────────────────────────────────────────────────────────────────────
// Driver: simultaneous fits per bin (11 parameters)
// ─────────────────────────────────────────────────────────────────────
// ─────────────────────────────────────────────────────────────────────
// Driver: simultaneous fits per bin (11 parameters)
// ─────────────────────────────────────────────────────────────────────
void performChi2Fits_GeneralExclusive(const char* output_file,
                                      const char* kinematic_file,
                                      const char* kinematicPlot_file,
                                      const std::string& prefix) {
  // === PRINT SWITCH for the fit-results LaTeX table ============================
  // If true  -> print everything (spin-dependent + UU terms + A_T leakages)
  // If false -> print only spin-dependent SF ratios: FLU^sin, FUL^sin, FUL^sin2, ALL, ALL^cos
  static bool g_ge_write_all_results = false;  // unchanged

  // === NEW: simple on/off switch for the extra export file ====================
  static bool g_ge_write_bin_export = true;

  // Control the leakage fits here (global switches) — unchanged
  g_fit_enable_TUL = true;
  g_fit_enable_TLL = true;
  g_fit_fixed_AT_UL = 0.0;
  g_fit_fixed_AT_LL = 0.0;

  // Prepare output streams for arrays (unchanged)
  std::ostringstream sALUoff, sAULoff, sALU, sAUL, sAUL2, sAT_UL, sALL, sALLc, sAUUc, sAUUc2, sAT_LL;
  for (auto* s : {&sALUoff,&sAULoff,&sALU,&sAUL,&sAUL2,&sAT_UL,&sALL,&sALLc,&sAUUc,&sAUUc2,&sAT_LL})
    (*s) << std::fixed << std::setprecision(9);

  sALUoff << prefix << "GEchi2FitsALUoffset = {";
  sAULoff << prefix << "GEchi2FitsAULoffset = {";
  sALU    << prefix << "GEchi2FitsALUsinphi = {";
  sAUL    << prefix << "GEchi2FitsAULsinphi = {";
  sAUL2   << prefix << "GEchi2FitsAULsin2phi = {";
  sAT_UL  << prefix << "GEchi2FitsA_T_UL = {";
  sALL    << prefix << "GEchi2FitsALL = {";
  sALLc   << prefix << "GEchi2FitsALLcosphi = {";
  sAUUc   << prefix << "GEchi2FitsAUUcosphi = {";
  sAUUc2  << prefix << "GEchi2FitsAUUcos2phi = {";
  sAT_LL  << prefix << "GEchi2FitsA_T_LL = {";

  // Kinematic LaTeX table (unchanged other than sourcing -t' from tprime)
  std::ostringstream kinLatex;
  kinLatex << "\\begin{table}[h]\n\\centering\n"
           << "\\begin{tabular}{|c|c|c|c|c|c|} \\hline\n"
           << "Bin & $Q^{2}$ & $x_{B}$ & $y$ & $z$ & $-t'$ \\\\ \\hline\n";

  // Keep kinList for downstream plotting/scripts (unchanged)
  std::ostringstream kinList;
  kinList << prefix << "GEKinematics = {";

  gSystem->mkdir("output/results", kTRUE);

  // Build suffix from output_file (unchanged)
  auto deriveSuffixFromOut = [](const char* outPath)->std::string{
    std::string s = outPath ? outPath : "";
    size_t slash = s.find_last_of("/\\");
    std::string base = (slash==std::string::npos)? s : s.substr(slash+1);
    size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0,dot);
    const std::string prefixes[] = {"asymmetries_", "kinematics_", "kinematicPlots_"};
    for (const auto& pfx : prefixes) {
      if (base.rfind(pfx, 0) == 0) { base = base.substr(pfx.size()); break; }
    }
    return base;
  };
  const std::string suffix = deriveSuffixFromOut(output_file);
  const std::string covPath  = "output/results/GE_" + prefix + "_cov_"  + suffix + ".txt";
  const std::string corrPath = "output/results/GE_" + prefix + "_corr_" + suffix + ".txt";

  // === NEW: open the extra per-bin export file (header once) ==================
  std::ofstream binExport;
  std::string binExportPath;
  if (g_ge_write_bin_export) {
    binExportPath = "output/results/GE_bin_export_" + prefix + "_" + suffix + ".txt";
    binExport.open(binExportPath, std::ios::out | std::ios::trunc);
    binExport << std::fixed << std::setprecision(9);
    binExport << "# x  -t  -t'  y  Q2\n";
  }

  // Containers to write the *fit-results* LaTeX table after the loop (unchanged)
  std::vector<double> meanVars;
  std::vector<std::vector<double>> all_pvals;  // [bin][0..10]
  std::vector<std::vector<double>> all_perrs;  // [bin][0..10]

  // Parameter names/order (11 total) — unchanged
  const int npar = 11;
  const char* names[npar] = {
    "ALU_offset","AUL_offset",
    "F_LU_sin/F_UU","F_UL_sin/F_UU","F_UL_sin2/F_UU",
    "F_LL/F_UU","F_LL_cos/F_UU",
    "F_UU_cos/F_UU","F_UU_cos2/F_UU",
    "A_T_UL","A_T_LL"
  };

  const size_t numBins = allBins[currentFits].size() - 1;

  for (size_t i = 0; i < numBins; ++i) {
    std::cout << "Beginning simultaneous chi2 GE fit for " << binNames[currentFits]
              << " bin " << i << ". " << std::endl;

    // Build histograms and per-φ-bin ⟨sinθγ⟩ (centered+normalized and smoothed)
    char hname[64]; snprintf(hname, sizeof(hname), "GE_%zu", i);
    TH1D *hALU, *hAUL, *hALL;
    std::tie(hALU, hAUL, hALL) =
      createHistogramForBin_GeneralExclusive(hname, (int)i, prefix);

    // ---- Mean & range kinematics (track min/mean/max for Q2,x,y,z, and −t′) ----
    double sumQ2=0, sumW=0, sumx=0, sumy=0, sumz=0, sumt=0, sumtmin=0, sumVar=0, nEvt=0;

    // Ranges (use extreme sentinels; signs are correct)
    double q2min= 1e300, q2max=-1e300;
    double xmin = 1e300, xmax =-1e300;
    double ymin = 1e300, ymax =-1e300;
    double zmin = 1e300, zmax =-1e300;

    // −t′ (from the tprime branch)
    double mtp_min = 1e300, mtp_max = -1e300;
    double sum_mtp = 0.0;

    // Optional moments retained (do not affect behavior elsewhere)
    double sumCos = 0.0, sumCos2 = 0.0; long long nCos = 0;

    {
      TTreeReaderValue<double> Q2(dataReader,"Q2");
      TTreeReaderValue<double> W (dataReader,"W");
      TTreeReaderValue<double> x (dataReader,"x");
      TTreeReaderValue<double> y (dataReader,"y");
      TTreeReaderValue<double> z (dataReader,"z");
      TTreeReaderValue<double> t (dataReader,"t");
      TTreeReaderValue<double> tmin(dataReader,"tmin");
      TTreeReaderValue<double> tprime(dataReader,"tprime");
      TTreeReaderValue<double> phi(dataReader,"phi");
      TTreeReaderValue<double> currentVariable(dataReader, propertyNames[currentFits].c_str());

      const double vmin = allBins[currentFits][i];
      const double vmax = allBins[currentFits][i+1];

      auto upd_minmax = [](double v, double& mn, double& mx) {
        if (v < mn) mn = v;
        if (v > mx) mx = v;
      };

      while (dataReader.Next()){
        if (!kinematicCuts->applyCuts(currentFits,false)) continue;
        if (*currentVariable < vmin || *currentVariable >= vmax) continue;

        // Accumulate means
        sumQ2 += *Q2; sumW += *W; sumx += *x; sumy += *y; sumz += *z;
        sumt  += *t;  sumtmin += *tmin; sumVar += *currentVariable; nEvt += 1.0;

        // Ranges
        upd_minmax(*Q2, q2min, q2max);
        upd_minmax(*x,  xmin,  xmax);
        upd_minmax(*y,  ymin,  ymax);
        upd_minmax(*z,  zmin,  zmax);

        // −t′ directly from tprime (store as positive number)
        const double mtp = -(*tprime);
        sum_mtp += mtp;
        upd_minmax(mtp, mtp_min, mtp_max);

        // keep ⟨cosφ⟩ moment (not used in export now)
        const double c = std::cos(*phi);
        sumCos  += c;
        sumCos2 += c*c;
        ++nCos;
      }
      dataReader.Restart();
    }

    // Means (guard nEvt=0)
    const double meanVar  = (nEvt>0)? sumVar/nEvt : 0.0;
    const double meanQ2   = (nEvt>0)? sumQ2/nEvt  : 0.0;
    const double meanW    = (nEvt>0)? sumW/nEvt   : 0.0;
    const double meanx    = (nEvt>0)? sumx/nEvt   : 0.0;
    const double meany    = (nEvt>0)? sumy/nEvt   : 0.0;
    const double meanz    = (nEvt>0)? sumz/nEvt   : 0.0;
    const double meant    = (nEvt>0)? sumt/nEvt   : 0.0;
    const double meantmin = (nEvt>0)? sumtmin/nEvt: 0.0;
    const double mean_mtp = (nEvt>0)? sum_mtp/nEvt : 0.0;   // mean of (−t′)

    // If no events, make ranges [0,0,0]
    if (nEvt == 0) {
      q2min=q2max=0.0; xmin=xmax=0.0; ymin=ymax=0.0; zmin=zmax=0.0; mtp_min=mtp_max=0.0;
    }

    // Mean and SEM for ⟨cosφ⟩ (not used in export now; retained to avoid side-effects)
    double meanCos = 0.0, eMeanCos = 0.0;
    if (nCos > 0) {
      meanCos = sumCos / static_cast<double>(nCos);
      if (nCos > 1) {
        const double s2 = std::max(0.0,
          (sumCos2 - static_cast<double>(nCos)*meanCos*meanCos) / static_cast<double>(nCos - 1));
        eMeanCos = std::sqrt(s2 / static_cast<double>(nCos));
      } else {
        eMeanCos = 0.0;
      }
    }

    // Depolarization ratios (unchanged)
    double sumDepA=0, sumDepB=0, sumDepC=0, sumDepV=0, sumDepW=0;
    {
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
        sumDepA += *DepA; sumDepB += *DepB; sumDepC += *DepC; sumDepV += *DepV; sumDepW += *DepW;
      }
      dataReader.Restart();
    }
    const double depA = (nEvt>0)? (sumDepA/nEvt) : 1.0;
    const double depB = (nEvt>0)? (sumDepB/nEvt) : 1.0;
    const double depC = (nEvt>0)? (sumDepC/nEvt) : 1.0;
    const double depV = (nEvt>0)? (sumDepV/nEvt) : 1.0;
    const double depW = (nEvt>0)? (sumDepW/nEvt) : 1.0;

    // Pass to FCN (unchanged)
    g_ge_ctx.hLU = hALU;
    g_ge_ctx.hUL = hAUL;
    g_ge_ctx.hLL = hALL;
    g_ge_ctx.rVA = (depA!=0.0)? (depV/depA) : 1.0;  // <V/A>
    g_ge_ctx.rBA = (depA!=0.0)? (depB/depA) : 1.0;  // <B/A>
    g_ge_ctx.rWA = (depA!=0.0)? (depW/depA) : 1.0;
    g_ge_ctx.rCA = (depA!=0.0)? (depC/depA) : 1.0;

    // Minuit (unchanged)
    TMinuit minuit(npar);
    minuit.SetPrintLevel(-1);
    minuit.SetErrorDef(1.0);
    minuit.SetFCN(chi2Fcn_GeneralExclusive);

    // name, initial value, step, low, up (unchanged)
    minuit.DefineParameter(0,  "ALU_offset",      0.00,  0.01,  -0.1,  0.1);
    minuit.DefineParameter(1,  "AUL_offset",      0.00,  0.01,  -0.1,  0.1);
    minuit.DefineParameter(2,  "F_LU_sin/F_UU",   0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(3,  "F_UL_sin/F_UU",   0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(4,  "F_UL_sin2/F_UU",  0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(5,  "F_LL/F_UU",       0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(6,  "F_LL_cos/F_UU",   0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(7,  "F_UU_cos/F_UU",   0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(8,  "F_UU_cos2/F_UU",  0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(9,  "A_T_UL",          0.00,  0.01,  -1.0,  1.0);
    minuit.DefineParameter(10, "A_T_LL",          0.00,  0.01,  -1.0,  1.0);

    minuit.FixParameter(0);
    minuit.FixParameter(1);
    minuit.FixParameter(2);
    minuit.FixParameter(5);
    minuit.FixParameter(8);
    minuit.FixParameter(9);
    minuit.FixParameter(10);

    // Fix leakage terms if ⟨sinθγ⟩ is flat (unchanged)
    if (!g_fit_enable_TUL || g_ge_ctx.sTG_wstd <= 1e-4) {
      minuit.FixParameter(9);
      if (g_ge_ctx.sTG_wstd <= 1e-4) {
        std::cout << "  [info] Bin " << i << ": sTG_wstd ≈ 0; fixing A_T_UL=0." << std::endl;
      }
    }
    if (!g_fit_enable_TLL || g_ge_ctx.sTG_wstd <= 1e-4) {
      minuit.FixParameter(10);
      if (g_ge_ctx.sTG_wstd <= 1e-4) {
        std::cout << "  [info] Bin " << i << ": sTG_wstd ≈ 0; fixing A_T_LL=0." << std::endl;
      }
    }

    // Strategy / minimize (unchanged)
    {
      double arglist[2]; int ier=0;
      arglist[0] = 2;               minuit.mnexcm("SET STR", arglist, 1, ier);
      arglist[0] = 5000; arglist[1] = 0.1;  minuit.mnexcm("MINImize", arglist, 2, ier);
      if (ier != 0) {
        arglist[0] = 2000;          minuit.mnexcm("SIMPLEX", arglist, 1, ier);
        arglist[0] = 5000; arglist[1] = 0.1;  minuit.mnexcm("MIGRAD", arglist, 2, ier);
      }
      minuit.mnexcm("HESSE", nullptr, 0, ier);
    }

    // Fit status and results (unchanged)
    double fmin, edm, errdef; int npari, nparx, istat;
    minuit.mnstat(fmin, edm, errdef, npari, nparx, istat);

    double pval[11], perr[11];
    for (int ip=0; ip<npar; ++ip) minuit.GetParameter(ip, pval[ip], perr[ip]);

    // GLOBAL dof (= total valid points – nvar) (unchanged)
    auto count_valid_points = [](TH1D* h){
      int n=0; if (!h) return n;
      const int nb = h->GetNbinsX();
      for (int i=1;i<=nb;++i){
        const double y=h->GetBinContent(i), e=h->GetBinError(i);
        if (std::isfinite(y) && std::isfinite(e) && e>0) ++n;
      }
      return n;
    };
    const int npts_total = count_valid_points(hALU) + count_valid_points(hAUL) + count_valid_points(hALL);
    const int ndf_global = std::max(0, npts_total - npar);
    const double chi2_global = fmin;  // FCN returns χ²

    // Plot (unchanged)
    plotHistogramAndFit_GeneralExclusive(hALU, hAUL, hALL, pval, perr,
                                         (int)i, prefix, suffix,
                                         chi2_global, ndf_global);

    // Append to arrays (unchanged)
    sALUoff << "{" << meanVar << ", " << pval[0]   << ", " << perr[0]   << "}";
    sAULoff << "{" << meanVar << ", " << pval[1]   << ", " << perr[1]   << "}";
    sALU    << "{" << meanVar << ", " << pval[2]   << ", " << perr[2]   << "}";
    sAUL    << "{" << meanVar << ", " << pval[3]   << ", " << perr[3]   << "}";
    sAUL2   << "{" << meanVar << ", " << pval[4]   << ", " << perr[4]   << "}";
    sALL    << "{" << meanVar << ", " << pval[5]   << ", " << perr[5]   << "}";
    sALLc   << "{" << meanVar << ", " << pval[6]   << ", " << perr[6]   << "}";
    sAUUc   << "{" << meanVar << ", " << pval[7]   << ", " << perr[7]   << "}";
    sAUUc2  << "{" << meanVar << ", " << pval[8]   << ", " << perr[8]   << "}";
    sAT_UL  << "{" << meanVar << ", " << pval[9]   << ", " << perr[9]   << "}";
    sAT_LL  << "{" << meanVar << ", " << pval[10]  << ", " << perr[10]  << "}";

    // Collect for the *fit-results* LaTeX table (unchanged)
    meanVars.push_back(meanVar);
    {
      std::vector<double> pv(11), pe(11);
      for (int k=0;k<11;++k){ pv[k]=pval[k]; pe[k]=perr[k]; }
      all_pvals.push_back(std::move(pv));
      all_perrs.push_back(std::move(pe));
    }

    if (i < numBins - 1) {
      sALUoff << ", "; sAULoff << ", "; sALU << ", "; sAUL << ", "; sAUL2 << ", ";
      sALL << ", "; sALLc << ", "; sAUUc << ", "; sAUUc2 << ", ";
      sAT_UL << ", "; sAT_LL << ", ";
    }

    // LaTeX row: [min, mean, max] (with −t′) (unchanged except source of −t′)
    auto triple = [](double mn, double mu, double mx){
      std::ostringstream o; o.setf(std::ios::fixed); o<<std::setprecision(2)
        << "[" << mn << ", " << mu << ", " << mx << "]";
      return o.str();
    };
    kinLatex << (i+1) << " ~&~ "
             << triple(q2min,  meanQ2,  q2max)   << " ~&~ "
             << triple(xmin,   meanx,   xmax)    << " ~&~ "
             << triple(ymin,   meany,   ymax)    << " ~&~ "
             << triple(zmin,   meanz,   zmax)    << " ~&~ "
             << triple(mtp_min, mean_mtp, mtp_max)
             << " \\\\ \\hline\n";

    // Keep kinList content (means) unchanged for downstream use
    kinList << "{" << meanQ2 << ", " << meanW << ", " << meanx << ", "
            << meany << ", " << meant << ", " << meantmin << "}";
    if (i < numBins - 1) kinList << ", ";

    // Save covariance/correlation blocks (unchanged)
    std::vector<double> cov(npar*npar, 0.0);
    minuit.mnemat(cov.data(), npar);

    std::vector<double> errv(npar);
    for (int ip=0; ip<npar; ++ip) errv[ip] = std::sqrt(std::max(cov[ip*npar+ip], 0.0));

    const double vminB = allBins[currentFits][i];
    const double vmaxB = allBins[currentFits][i+1];

    {
      std::ofstream of(covPath, std::ios::out | std::ios::app);
      of << std::setprecision(9);
      of << "## Bin " << i << "  Range: [" << vminB << ", " << vmaxB << ")  Events: " << nEvt
         << "  npts=" << npts_total << "  npar=" << npar << "  ndf=" << ndf_global
         << "  chi2=" << chi2_global << "\n";
      of << std::left << std::setw(22) << "#";
      for (int c=0;c<npar;++c) of << std::setw(22) << names[c];
      of << "\n";
      for (int r=0; r<npar; ++r) {
        of << std::left << std::setw(22) << names[r];
        for (int c=0; c<npar; ++c) of << std::setw(22) << cov[r*npar + c];
        of << "\n";
      }
      of << "\n";
    }
    {
      std::ofstream of(corrPath, std::ios::out | std::ios::app);
      of << std::setprecision(9);
      of << "## Bin " << i << "  Range: [" << vminB << ", " << vmaxB << ")  Events: " << nEvt << "\n";
      for (int r=0; r<npar; ++r) {
        for (int c=0; c<npar; ++c) {
          double denom = (errv[r]>0 && errv[c]>0) ? (errv[r]*errv[c]) : 1.0;
          double rho = cov[r*npar + c] / denom;
          of << std::setw(22) << rho;
        }
        of << "\n";
      }
      of << "\n";
    }

    // === NEW: write one line into the extra export file =======================
    if (g_ge_write_bin_export && binExport.is_open()) {
      const double mean_mt = -meant;      // ⟨-t⟩
      // mean_mtp is already ⟨-t'⟩
      binExport << meanx    << " "
                << mean_mt  << " "
                << mean_mtp << " "
                << meany    << " "
                << meanQ2   << "\n";
    }

    delete hALU; delete hAUL; delete hALL;
  }

  // Close arrays and write to file (unchanged)
  sALUoff << "};"; sAULoff << "};"; sALU << "};"; sAUL << "};";
  sAUL2  << "};"; sALL << "};"; sALLc<< "};"; sAUUc << "};"; sAUUc2 << "};";
  sAT_UL << "};"; sAT_LL << "};";

  {
    std::ofstream out(output_file, std::ios::app);
    out << sALUoff.str() << "\n";
    out << sAULoff.str() << "\n";
    out << sALU.str()    << "\n";
    out << sAUL.str()    << "\n";
    out << sAUL2.str()   << "\n";
    out << sAT_UL.str()  << "\n";
    out << sALL.str()    << "\n";
    out << sALLc.str()   << "\n";
    out << sAUUc.str()   << "\n";
    out << sAUUc2.str()  << "\n";
    out << sAT_LL.str()  << "\n";
  }

  // Finish LaTeX kinematics table & kinematics list (unchanged aside from -t′ source)
  kinLatex << "\\end{tabular}\n"
           << "\\caption{Per-bin kinematics shown as [min, mean, max] for $Q^{2}$, $x_{B}$, $y$, $z$, and $-t'$. "
           << "$Q^{2}$ and $-t'$ are in GeV$^{2}$.}\n"
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

  // =================================================================
  // Write the *fit-results* LaTeX table (unchanged in content)
  // =================================================================
  {
    // Determine variable name being fit against
    std::string varName = propertyNames[currentFits];

    // Normalize to lower/compact form for detection
    auto toLower = [](std::string s){
      for (auto& ch : s) ch = (char)std::tolower((unsigned char)ch);
      return s;
    };
    auto compact = [](std::string s){
      std::string out; out.reserve(s.size());
      for (char ch : s) if (ch!=' ' && ch!='_' && ch!='-' ) out.push_back(ch);
      return out;
    };
    const std::string key = compact(toLower(varName));

    const bool is_tprime = (key=="tprime" || key=="t'" || key=="tminustmin");
    const bool is_t      = (!is_tprime && (key=="t"));

    // Build LaTeX label for the first column
    std::string varLabel;
    if (is_tprime)      varLabel = "$\\langle -t' \\rangle$";
    else if (is_t)      varLabel = "$\\langle -t \\rangle$";
    else {
      try { varLabel = "$\\langle " + formatLabelName(varName) + " \\rangle$"; }
      catch (...) { varLabel = "$\\langle " + varName + " \\rangle$"; }
    }

    // Column spec (unchanged)
    struct Col { int idx; const char* label; double sysFrac; bool showSyst; };
    std::vector<Col> cols_spin = {
      { 2, "$F_{LU}^{\\sin\\phi}/F_{UU}$",           0.06, true }, // BSA
      { 3, "$F_{UL}^{\\sin\\phi}/F_{UU}$",           0.08, true }, // TSA
      { 4, "$F_{UL}^{\\sin2\\phi}/F_{UU}$",          0.08, true }, // TSA
      { 5, "$F_{LL}/F_{UU}$",                        0.10, true }, // DSA
      { 6, "$F_{LL}^{\\cos\\phi}/F_{UU}$",           0.10, true }  // DSA
    };

    // Full set adds UU terms and leakage amplitudes (unchanged)
    std::vector<Col> cols_full = cols_spin;
    cols_full.push_back({ 7, "$F_{UU}^{\\cos\\phi}/F_{UU}$",        0.00, false });
    cols_full.push_back({ 8, "$F_{UU}^{\\cos2\\phi}/F_{UU}$",       0.00, false });
    cols_full.push_back({ 9, "$A_{T\\text{-}UL}^{\\sin\\phi}$",     0.08, false });
    cols_full.push_back({10, "$A_{T\\text{-}LL}^{\\sin\\phi}$",     0.10, false });

    const auto& cols = g_ge_write_all_results ? cols_full : cols_spin;

    // File path for the fit-results table (unchanged)
    std::string varToken = varName;
    for (auto& ch : varToken) if (ch==' ' || ch=='\'') ch = '_';
    const std::string fitOutPath = "output/results/fit_results_GE_" + varToken + "_" + suffix + ".tex";
    std::ofstream out(fitOutPath, std::ios::out | std::ios::trunc);
    if (!out) {
      std::cerr << "[performChi2Fits_GeneralExclusive] Failed to open " << fitOutPath << " for writing.\n";
      if (binExport.is_open()) binExport.close();
      return;
    }

    // Helper: value^{±stat}_{±syst} (unchanged)
    auto entry = [](double v, double eStat, double fracSys, bool showSyst) {
      std::ostringstream s; s.setf(std::ios::fixed);
      s << std::setprecision(3) << v
        << "^{\\pm " << std::setprecision(3) << eStat << "}";
      if (showSyst) {
        const double eSys = std::fabs(v) * fracSys;
        s << "_{\\pm " << std::setprecision(3) << eSys << "}";
      }
      return s.str();
    };

    // Header (unchanged)
    std::ostringstream header;
    header << "\\begin{table}[h]\n\\centering\n"
           << "\\small\n"
           << "\\begin{tabular}{|c";
    for (size_t i=0;i<cols.size();++i) header << "|c";
    header << "|} \\hline\n";

    header << varLabel;
    for (const auto& c : cols) header << " & " << c.label;
    header << " \\\\ \\hline\n";
    out << header.str();

    // Rows (unchanged)
    for (size_t ib=0; ib<meanVars.size(); ++ib) {
      const double meanDisplay = (is_t || is_tprime) ? -meanVars[ib] : meanVars[ib];

      std::ostringstream row; row.setf(std::ios::fixed);
      row << std::setprecision(3) << meanDisplay;
      for (const auto& c : cols) {
        const double v  = all_pvals[ib][c.idx];
        const double es = all_perrs[ib][c.idx];
        row << " ~&~ " << entry(v, es, c.sysFrac, c.showSyst);
      }
      row << " \\\\ \\hline\n";
      out << row.str();
    }

    // Footer (unchanged)
    out << "\\end{tabular}\n"
        << "\\caption{Fitted structure-function ratios per bin. Entries are "
           "$\\text{value}^{\\pm\\,\\text{stat}}_{\\pm\\,\\text{syst}}$."
        << "\\label{table:GE_fitresults_" << varToken << "}\n"
        << "\\end{table}\n";
    out.close();
  }

  // binExport will auto-close on scope exit (ofstream destructor).
}
// ===================== GeneralExclusive (end) =====================