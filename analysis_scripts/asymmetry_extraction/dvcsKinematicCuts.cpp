#include "dvcsKinematicCuts.h"
#include <iostream>
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

dvcsKinematicCuts::dvcsKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), Q2(reader, "Q2"), W(reader, "W"), 
      x(reader, "x"), y(reader, "y"), t1(reader, "t1"), eta2(reader, "eta2"),
      Emiss2(reader, "Emiss2"), theta_gamma_gamma(reader, "theta_gamma_gamma"),
      pTmiss(reader, "pTmiss"), Mxgammasquared(reader, "Mxgammasquared") {}

bool dvcsKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    // goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75 && *eta2 < 0 && *t1 > -2;// DIS cuts (and radiative photon in detector)
    goodEvent = true;
    
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
      return goodEvent;
    } else {
      // return goodEvent && *target_pol!=0;
      return goodEvent;
    }
    return false;
}