#pragma once
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"
#include "BaseKinematicCuts.h" 

class dvcsKinematicCuts : public BaseKinematicCuts { // Inherit from BaseKinematicCuts
public:
    dvcsKinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC) override; // Override applyCuts method

private:
    TTreeReaderValue<int> runnum;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> t1;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> Emiss2;
    TTreeReaderValue<double> theta_gamma_gamma;
    TTreeReaderValue<double> pTmiss;
    TTreeReaderValue<double> Mxgammasquared;
    TTreeReaderValue<double> eta2;
};