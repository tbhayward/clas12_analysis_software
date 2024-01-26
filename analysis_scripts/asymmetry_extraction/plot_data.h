#include <TTreeReader.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#ifndef PLOT_DATA_H
#define PLOT_DATA_H

// Function declarations
void createIntegratedKinematicPlots();
void createIntegratedKinematicPlotsForBinsAndFits();
void createCorrelationPlotsforrunnum();
void createCorrelationPlots();

template<typename T>
void FillHistogram(TTreeReader& reader, const std::string& branchName, TH1D* hist, KinematicCuts& kinematicCuts, int fitIndex);

template<typename T1, typename T2>
void createAndFillHistogram(TTreeReader& reader, TH2D* hist, const std::string& branchX, const std::string& branchY, KinematicCuts& kinematicCuts);


#endif // PLOT_DATA_H
