#include <TTreeReader.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include "TGraphErrors.h"
#include "BaseKinematicCuts.h"
#include "KinematicCuts.h"
#include "InclusiveKinematicCuts.h"
#include "SingleHadronKinematicCuts.h"
#include "B2BDihadronKinematicCuts.h"
#include "DihadronKinematicCuts.h"
#include "histConfigs.h"

extern std::map<std::string, HistConfig> histConfigs;

#ifndef PLOT_DATA_H
#define PLOT_DATA_H

// Function declarations
void createIntegratedKinematicPlots();
void createIntegratedKinematicPlotsForBinsAndFits();
void createCorrelationPlotsforrunnum();
void createCorrelationPlots();
void createElectronMisIDRatePlots();

// Updated FillHistogram declaration
template<typename T>
void FillHistogram(TTreeReader& reader, const std::string& branchName, TH1D* hist,
  BaseKinematicCuts& kinematicCuts, int fitIndex, bool isMC);

template<typename T>
void FillHistogramForBins(TTreeReader& reader, const std::string& branchName, TH1D* hist,
  BaseKinematicCuts& kinematicCuts, int fitIndex, bool isMC,
  double varMin, double varMax);


template<typename T1, typename T2>
void createAndFillHistogram(TTreeReader& reader, TH2D* hist, const std::string& branchX, 
  const std::string& branchY, KinematicCuts& kinematicCuts);

#endif // PLOT_DATA_H
