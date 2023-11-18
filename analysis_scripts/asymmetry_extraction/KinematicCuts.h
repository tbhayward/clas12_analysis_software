#pragma once
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"

class KinematicCuts {
public:
    KinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC);

private:
    TTreeReaderValue<double> p1_p;
    TTreeReaderValue<double> p1_phi;
    TTreeReaderValue<double> p1_theta;
    TTreeReaderValue<double> p2_p;
    TTreeReaderValue<double> p2_phi;
    TTreeReaderValue<double> p2_theta;
    TTreeReaderValue<double> p3_p;
    TTreeReaderValue<double> p3_phi;
    TTreeReaderValue<double> p3_theta;
    // TTreeReaderValue<double> Mh;
    // TTreeReaderValue<double> Mh23;
    // TTreeReaderValue<double> z23;
    // TTreeReaderValue<double> Mx1;
    // TTreeReaderValue<double> Mx23;

    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> Mx;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> z;
    TTreeReaderValue<double> pT;
    TTreeReaderValue<double> xF;
    TTreeReaderValue<double> target_pol;
};
