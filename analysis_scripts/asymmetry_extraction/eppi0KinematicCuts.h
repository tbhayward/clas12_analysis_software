#pragma once
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"
#include "BaseKinematicCuts.h" 

class eppi0KinematicCuts : public BaseKinematicCuts { // Inherit from BaseKinematicCuts
public:
    eppi0KinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC) override; // Override applyCuts method

private:
    TTreeReaderValue<int> runnum;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> t1;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> Emiss2;
    TTreeReaderValue<double> Mx2;
    TTreeReaderValue<double> Mx2_1;
    TTreeReaderValue<double> Mx2_2;
    TTreeReaderValue<double> theta_pi0_pi0;
    TTreeReaderValue<double> open_angle_ep2;
    TTreeReaderValue<double> pTmiss;
    TTreeReaderValue<double> eta2;
};