#include "InclusiveKinematicCuts.h"
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

InclusiveKinematicCuts::InclusiveKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), Q2(reader, "Q2"), W(reader, "W"), 
      x(reader, "x"), y(reader, "y"), target_pol(reader, "target_pol") {}

bool InclusiveKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    if (property == "eX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else if (property == "xeX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else {
      std::cout << "Property, " << property << ", not detected." << std::endl;
    }
    
    if (isMC || *runnum < 11571) {
      return goodEvent
    } else {
      return goodEvent && *target_pol!=0;
    }
    return false;
}