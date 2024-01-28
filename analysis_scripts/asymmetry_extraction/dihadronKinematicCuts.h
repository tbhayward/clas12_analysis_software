#pragma once
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts

class DihadronKinematicCuts : public BaseKinematicCuts { // Inherit from BaseKinematicCuts
public:
    DihadronKinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC) override; // Override applyCuts method

private:
    TTreeReaderValue<int> runnum;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> target_pol;
};
