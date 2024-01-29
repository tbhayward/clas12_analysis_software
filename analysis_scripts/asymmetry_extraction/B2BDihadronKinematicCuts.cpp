#include "B2BDihadronKinematicCuts.h"
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

B2BDihadronKinematicCuts::B2BDihadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), Q2(reader, "Q2"), W(reader, "W"),
      x(reader, "x"), y(reader, "y"), z1(reader, "z1"), xF1(reader, "xF1"), xF2(reader, "xF2"), 
      Mx(reader, "Mx"), Mx1(reader, "Mx1"), Mx2(reader, "Mx2"), target_pol(reader, "target_pol") {}

bool B2BDihadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    if (property == "epippX") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else if (property == "xepimpX") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75 && *z1 > 0.2 && *xF1 > 0 && *xF2 < 0;
    } else if (property == "epimpX") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else if (property == "xepimpX") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75 && *z1 > 0.2 && *xF1 > 0 && *xF2 < 0;;
    } else {
      std::cout << "Property, " << property << ", not detected." << std::endl;
    }
    return goodEvent;
}
