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

    goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75; // DIS cuts
    goodEvent = goodEvent && *x > 0.05 && *x < 0.70;
    if (property == "b2bchannel") {
      goodEvent = goodEvent;
    } else if (property == "b2bchannelMxStudy") {
      goodEvent = goodEvent; // && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2;
    } else if (property == "b2bchannelz1Study") {
      goodEvent = goodEvent; // && *xF1 > 0 && *xF2 < 0 && *Mx > 0.95 && *Mx1 > 1.4 && *Mx2 > 1.8;
    } else if (property == "b2bchannelxFStudy") {
      goodEvent = goodEvent; // && *z1 > 0.2 && *Mx > 0.95 && *Mx1 > 1.4 && *Mx2 > 1.8;
    } else if (property == "b2banalysis") {
      goodEvent = goodEvent && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2 && *Mx > 0.95 && 
        *Mx1 > 1.4 && *Mx2 > 1.8;
    } else {
      std::cout << "Property, " << property << ", not detected." << std::endl;
    }
    
    if (isMC || *runnum < 11571) {
      return goodEvent;
    } else {
      return goodEvent && *target_pol!=0;
    }
    return false;
}