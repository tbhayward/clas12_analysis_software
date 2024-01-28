#include "DihadronKinematicCuts.h"
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

DihadronKinematicCuts::DihadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), Q2(reader, "Q2"), W(reader, "W"), 
      x(reader, "x"), y(reader, "y"), target_pol(reader, "target_pol") {}

bool DihadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    if (property == "epipppimX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else if (property == "ekpkmX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else {
        std::cout << "Property not detected" << std::endl;
    }
    return goodEvent;
}
