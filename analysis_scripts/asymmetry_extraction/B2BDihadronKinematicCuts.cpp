#include "B2BDihadronKinematicCuts.h"
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

B2BDihadronKinematicCuts::B2BDihadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), fiducial_status(reader, "fiducial_status"), 
      Q2(reader, "Q2"), W(reader, "W"), p1_p(reader, "p1_p"),
      x(reader, "x"), y(reader, "y"), z1(reader, "z1"), xF1(reader, "xF1"), xF2(reader, "xF2"), 
      Mx2(reader, "Mx2"), Mx2_1(reader, "Mx2_1"), target_pol(reader, "target_pol") {}

bool B2BDihadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    if (*Q2 < 1 || *W < 2 || *y > 0.75 || *fiducial_status!=3) return false;
    if (*p1_p < 1.2 || *xF1 < 0 || *Mx2_1 < 3.24) return false;
    return true;

    if (property == "PTPT") {
        goodEvent = goodEvent;
    } else if (property == "b2bchannel") {
      goodEvent = goodEvent;
    } else if (property == "b2bchannelMxStudy") {
      goodEvent = goodEvent; // && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2;
    } else if (property == "b2bchannelMx1Study") {
      goodEvent = goodEvent; // && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2;
    } else if (property == "b2bchannelMx2Study") {
      goodEvent = goodEvent; // && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2;
    } else if (property == "b2bchannelxF1Study") {
      goodEvent = goodEvent && *Mx > 0.95 && *Mx1 > 1.8 && *Mx2 > 1.4;
    } else if (property == "b2bchannelxF2Study") {
      goodEvent = goodEvent && *Mx > 0.95 && *Mx1 > 1.8 && *Mx2 > 1.4;
    } else if (property == "b2banalysis") {
      goodEvent = goodEvent && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2 && *Mx > 0.95 && 
        *Mx1 > 1.8 && *Mx2 > 1.4;
    } else if (property == "b2banalysisx") {
      goodEvent = goodEvent && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2 && *Mx > 0.95 && 
        *Mx1 > 1.8 && *Mx2 > 1.4;
    } else if (property == "b2banalysispTpT") {
      // goodEvent = goodEvent && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2 && *Mx > 0.95 && 
      //   *Mx1 > 1.8 && *Mx2 > 1.4;
      goodEvent = goodEvent && *xF1 > 0 && *xF2 < 0 && *z1 > 0.2 && *Mx > 1.2 && 
        *Mx1 > 1.5 && *Mx2 > 1.2;
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