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
      pTmiss(reader, "pTmiss"), open_angle_ep2(reader, "open_angle_ep2") {}

bool dvcsKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    goodEvent = *Q2 > 1 && *W > 2 && *eta2 < 0 && *t1 > -1 && *pTmiss < 0.15 && *Emiss2 < 1 && *theta_gamma_gamma < 0.7 && *open_angle_ep2 > 5;
    goodEvent = true;
    
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
      return goodEvent;
    } else {
      // return goodEvent && *target_pol!=0;
      return goodEvent;
    }
    return false;
}