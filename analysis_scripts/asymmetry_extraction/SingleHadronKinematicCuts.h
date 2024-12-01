#pragma once
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts

class SingleHadronKinematicCuts : public BaseKinematicCuts {
public:
    SingleHadronKinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC) override;

private:
    TTreeReaderValue<int> runnum;
    TTreeReaderValue<int> fiducial_status;
    TTreeReaderValue<double> e_phi;
    TTreeReaderValue<double> vz_e;
    TTreeReaderValue<double> p_phi;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> Mx2;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> t;
    TTreeReaderValue<double> tmin;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> z;
    TTreeReaderValue<double> xi;
    TTreeReaderValue<double> pT;
    TTreeReaderValue<double> xF;
    TTreeReaderValue<double> phi;
    TTreeReaderValue<double> target_pol;
};