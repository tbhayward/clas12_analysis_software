#pragma once
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"

class b2bDihadronKinematicCuts {
public:
    b2bDihadronKinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC);

private:

    TTreeReaderValue<int> runnum;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> target_pol;
};
