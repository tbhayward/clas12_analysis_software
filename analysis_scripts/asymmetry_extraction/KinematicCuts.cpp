#include "KinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>

using std::string;

KinematicCuts::KinematicCuts(TTreeReader& reader)
    : Q2(reader, "Q2"), W(reader, "W"), Mx(reader, "Mx"), 
      x(reader, "x"), y(reader, "y"), pT(reader, "pT"), 
      xF(reader, "xF"), target_pol(reader, "target_pol") {}

bool KinematicCuts::applyCuts(int currentFits, bool isMC) {
        bool goodEvent = false;
        string property = binNames[currentFits];

        if (property == "xF") {
            goodEvent = *Q2 > 1 && *W > 2 && *Mx > 1.4 && *y < 0.75;
        }
        if (property == "Mx") {
            goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
        }
        if (property == "Q2bin") {
            goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *x>0.2 && *x<0.3 && *pT>0.25 && 
            *pT<0.35 && *xF<0;
        }
        if (property == "PTTFR" || property ==  "xTFR" || property == "zetaTFR" || 
          property == "Q2TFR" || property ==  "x") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *xF<0;
        }
        if (property == "PTCFR" || property == "xCFR" || property == "zetaCFR" ||
          property == "Q2TFR") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *xF>0;
        } 
        //
        // epi+X
        if (property == "xFpip") { 
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75;
        }
        if (property == "PTTFRpip" || property ==  "xTFRpip" || property == "zTFRpip" || 
          property == "Q2TFRpip" || property ==  "xpip") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *xF<0;
        }
        if (property == "PTCFRpip" || property == "xCFRpip" || property == "zCFRpip" ||
          property == "Q2TFRpip") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *xF>0;
        }
        //
        // epi-X
        if (property == "xFpim") { 
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75;
        }
        if (property == "PTTFRpim" || property ==  "xTFRpim" || property == "zTFRpim" || 
          property == "Q2TFRpim" || property ==  "xpim") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *xF<0;
        }
        if (property == "PTCFRpim" || property == "xCFRpim" || property == "zCFRpim" ||
          property == "Q2TFRpim") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *xF>0;
        }
        //
        // epi+pi+X, exclusive rho
        if (property == "exclusiveRhoIntegrated") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && fabs(*Mx-0.95)<0.15;
        }

        if (isMC) {
            return goodEvent;
        } else {
            return goodEvent && *target_pol != 0;
        }
    }