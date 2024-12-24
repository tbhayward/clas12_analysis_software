#include "InclusiveKinematicCuts.h"
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

InclusiveKinematicCuts::InclusiveKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), Q2(reader, "Q2"), W(reader, "W"), 
      x(reader, "x"), y(reader, "y"), e_phi(reader, "e_phi"), target_pol(reader, "target_pol") {}

bool InclusiveKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];


    if (*fiducial_status != 1) return false; // fiducial cuts

    if (property == "eX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else if (property == "xeX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else {
      std::cout << "Property, " << property << ", not detected." << std::endl;
    }

    if (property == "xBsector0") { // meant to be all six sectors
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80;
      return goodEvent;
    } 
    if (property == "xBsector1") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && (*e_phi < 0.2 || *e_phi > 5.5);
      return goodEvent;
    } 
    if (property == "xBsector2") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && (*e_phi > 0.2 && *e_phi < 1.25);
      return goodEvent;
    } 
    if (property == "xBsector3") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && (*e_phi > 1.25 && *e_phi < 2.25);
      return goodEvent;
    } 
    if (property == "xBsector4") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && (*e_phi > 2.25 && *e_phi < 3.35);
      return goodEvent;
    } 
    if (property == "xBsector5") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && (*e_phi > 3.35 && *e_phi < 4.4);
      return goodEvent;
    } 
    if (property == "xBsector6") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && (*e_phi > 4.4 && *e_phi < 5.5);
      return goodEvent;
    }
    
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
      return goodEvent;
    } else {
      // return goodEvent && *target_pol!=0;
        return goodEvent;
    }
    return false;
}