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

    goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75 && *eta2 < 0 && *t1 > -2;// DIS cuts (and radiative photon in detector)
    goodEvent = goodEvent && *theta_gamma_gamma < 0.6 && *Emiss2 < 0.5 &&
        *pTmiss < 0.125 && *Mxgammasquared < 1.25;
    if (property == "dvcsx") {
        // goodEvent = true;
    } else if (property == "dvcst") {
        // goodEvent = true;
    } else if (property == "dvcsQ2") {
        // goodEvent = true;
    } else if (property == "dvcsQ2") {
        // goodEvent = true;
    } else if (property == "dvcsQ21x1") {
      goodEvent = goodEvent && *Q2 > 1.0 && *Q2 < 1.4 && *x > 0.00 && *x < 0.13;
    } else if (property == "dvcsQ21x2") {
      goodEvent = goodEvent && *Q2 > 1.0 && *Q2 < 1.4 && *x > 0.13 && *x < 0.21;
    } else if (property == "dvcsQ21x3") {
      goodEvent = goodEvent && *Q2 > 1.0 && *Q2 < 1.4 && *x > 0.21 && *x < 1.00;
    } else if (property == "dvcsQ22x1") {
      goodEvent = goodEvent && *Q2 > 1.4 && *Q2 < 1.8 && *x > 0.00 && *x < 0.13;
    } else if (property == "dvcsQ22x2") {
      goodEvent = goodEvent && *Q2 > 1.4 && *Q2 < 1.8 && *x > 0.13 && *x < 0.21;
    } else if (property == "dvcsQ22x3") {
      goodEvent = goodEvent && *Q2 > 1.4 && *Q2 < 1.8 && *x > 0.21 && *x < 1.00;
    } else if (property == "dvcsQ23x1") {
      goodEvent = goodEvent && *Q2 > 1.8 && *Q2 < 2.4 && *x > 0.00 && *x < 0.16;
    } else if (property == "dvcsQ23x2") {
      goodEvent = goodEvent && *Q2 > 1.8 && *Q2 < 2.4 && *x > 0.16 && *x < 0.26;
    } else if (property == "dvcsQ23x3") {
      goodEvent = goodEvent && *Q2 > 1.8 && *Q2 < 2.4 && *x > 0.26 && *x < 1.00;
    } else if (property == "dvcsQ24x1") {
      goodEvent = goodEvent && *Q2 > 2.4 && *Q2 < 3.5 && *x > 0.00 && *x < 0.21;
    } else if (property == "dvcsQ24x2") {
      goodEvent = goodEvent && *Q2 > 2.4 && *Q2 < 3.5 && *x > 0.21 && *x < 0.33;
    } else if (property == "dvcsQ24x3") {
      goodEvent = goodEvent && *Q2 > 2.4 && *Q2 < 3.5 && *x > 0.33 && *x < 1.00;
    } else if (property == "dvcsQ25x1") {
      goodEvent = goodEvent && *Q2 > 3.5 && *Q2 < 5.0 && *x > 0.00 && *x < 0.33;
    } else if (property == "dvcsQ25x2") {
      goodEvent = goodEvent && *Q2 > 3.5 && *Q2 < 5.0 && *x > 0.33 && *x < 1.00;
    } else if (property == "dvcsQ26x1") {
      goodEvent = goodEvent && *Q2 > 5.0 && *Q2 < 11.00 && *x > 0.00 && *x < 0.55;
    } else if (property == "dvcsQ26x2") {
      goodEvent = goodEvent && *Q2 > 5.0 && *Q2 < 11.00 && *x > 0.55 && *x < 1.00;
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