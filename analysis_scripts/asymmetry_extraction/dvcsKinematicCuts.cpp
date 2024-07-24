#include "B2BDihadronKinematicCuts.h"
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

dvcsKinematicCuts::dvcsKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), Q2(reader, "Q2"), W(reader, "W"), 
      x(reader, "x"), y(reader, "y"), t1(reader, "t1"), target_pol(reader, "target_pol"),
      Emiss2(reader, "Emiss2"), theta_gamma_gamma(reader, "theta_gamma_gamma"),
      pTmiss(reader, "pTmiss"), Mxgammasquared(reader, "Mxgammasquared") {}

bool B2BDihadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75; // DIS cuts
    if (property == "dvcsx") {
        goodEvent = true;
    } else {
      std::cout << "Property, " << property << ", not detected!" << std::endl;
    }
    
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
      return goodEvent;
    } else {
      // return goodEvent && *target_pol!=0;
      return goodEvent;
    }
    return false;
}