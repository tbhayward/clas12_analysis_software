#include <TTreeReader.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include "KinematicCuts.h"
#include "inclusiveKinematicCuts.h"
#include "b2bDihadronKinematicCuts.h"
#include "dihadronKinematicCuts.h"
#include "histConfigs.h"

extern std::map<std::string, HistConfig> histConfigs;

#ifndef PLOT_DATA_H
#define PLOT_DATA_H

// Function declarations
void createIntegratedKinematicPlots();
void createIntegratedKinematicPlotsForBinsAndFits();
void createCorrelationPlotsforrunnum();
void createCorrelationPlots();

template<typename T>
void FillHistogram(TTreeReader& reader, const std::string& branchName, TH1D* hist, KinematicCuts& kinematicCuts, int fitIndex, bool isMC);

template<typename T1, typename T2>
void createAndFillHistogram(TTreeReader& reader, TH2D* hist, const std::string& branchX, const std::string& branchY, KinematicCuts& kinematicCuts);


#endif // PLOT_DATA_H
