#include "eppi0KinematicCuts.h"
#include <iostream>
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

eppi0KinematicCuts::eppi0KinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), Q2(reader, "Q2"), W(reader, "W"),
      Mx2(reader, "Mx2"), Mx2_1(reader, "Mx2_1"), Mx2_2(reader, "Mx2_2"), 
      x(reader, "x"), y(reader, "y"), t1(reader, "t1"), eta2(reader, "eta2"),
      Emiss2(reader, "Emiss2"), theta_pi0_pi0(reader, "theta_pi0_pi0"),
      pTmiss(reader, "pTmiss"), open_angle_ep2(reader, "open_angle_ep2") {}

bool eppi0KinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    goodEvent = *Q2 > 1 && *W > 2 && *t1 > -1 && *pTmiss < 0.09 && *Emiss2 < 0.5 && *theta_pi0_pi0 < 0.7 && *open_angle_ep2 > 5 && *Mx2 < 0.01 && *Mx2 > -0.01;
    goodEvent = true;
    
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
      return goodEvent;
    } else {
      // return goodEvent && *target_pol!=0;
      return goodEvent;
    }
    return false;
}