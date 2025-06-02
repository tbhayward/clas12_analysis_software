#pragma once

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"
#include "BaseKinematicCuts.h"

class SingleHadronKinematicCuts : public BaseKinematicCuts {
public:
    explicit SingleHadronKinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC) override;

private:
    // Branch readers (added e_p, e_theta, p_theta for t‐calculation)
    TTreeReaderValue<int>    runnum;
    TTreeReaderValue<int>    fiducial_status;
    TTreeReaderValue<double> e_p;       // scattered‐electron momentum
    TTreeReaderValue<double> e_theta;   // polar angle of scattered electron
    TTreeReaderValue<double> e_phi;     // azimuth of scattered electron

    TTreeReaderValue<double> p_p;       // pion momentum
    TTreeReaderValue<double> p_theta;   // polar angle of pion
    TTreeReaderValue<double> p_phi;     // azimuth of pion

    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> Mx2;
    TTreeReaderValue<double> xF;
    TTreeReaderValue<double> pT;        // pion pT
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> xi;
    TTreeReaderValue<double> phi;
    TTreeReaderValue<double> z;
    TTreeReaderValue<double> t;         // unused for “enpi+,” we recalc
    TTreeReaderValue<double> tmin;
    TTreeReaderValue<double>    target_pol;
};