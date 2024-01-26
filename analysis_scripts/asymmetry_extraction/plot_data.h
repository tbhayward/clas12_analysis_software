#ifndef PLOT_DATA_H
#define PLOT_DATA_H

#include <TTreeReader.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include "KinematicCuts.h"
#include "histConfigs.h"

extern std::map<std::string, HistConfig> histConfigs;

// Function declarations
void createIntegratedKinematicPlots();
void createIntegratedKinematicPlotsForBinsAndFits();
void createCorrelationPlotsforrunnum();
void createCorrelationPlots();

// FillHistogram template definition
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

template<typename T1, typename T2>
void createAndFillHistogram(TTreeReader& reader, TH2D* hist, const std::string& branchX, 
                            const std::string& branchY, KinematicCuts& kinematicCuts);

#endif // PLOT_DATA_H
