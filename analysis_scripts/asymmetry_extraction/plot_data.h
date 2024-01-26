#ifndef PLOT_DATA_H
#define PLOT_DATA_H

// Function declarations
void createIntegratedKinematicPlots();
void createIntegratedKinematicPlotsForBinsAndFits();
void createCorrelationPlotsforrunnum();
void createCorrelationPlots();

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
                            const std::string& branchY, KinematicCuts& kinematicCuts) {
    TTreeReaderValue<T1> valX(reader, branchX.c_str());
    TTreeReaderValue<T2> valY(reader, branchY.c_str());

    reader.Restart();
    while (reader.Next()) {
        if (kinematicCuts.applyCuts(0, false)) {
            hist->Fill(*valX, *valY);
        }
    }
}

#endif // PLOT_DATA_H
