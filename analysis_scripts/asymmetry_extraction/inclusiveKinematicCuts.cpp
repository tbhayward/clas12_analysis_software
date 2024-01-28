#include "inclusiveKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>

using std::string;

inclusiveKinematicCuts::inclusiveKinematicCuts(TTreeReader& reader)
    : runnum(reader, "runnum"),  Q2(reader, "Q2"), W(reader, "W"), 
      x(reader, "x"), y(reader, "y"),  target_pol(reader, "target_pol") {}

bool inclusiveKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    if (property == "eX") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else {
      std::cout << "Property not detected" << std::endl;
    }
    return goodEvent;
}