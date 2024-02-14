#include "SingleHadronKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>

using std::string;

SingleHadronKinematicCuts::SingleHadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Call to the BaseKinematicCuts constructor
      runnum(reader, "runnum"), p_p(reader, "p_p"), p_theta(reader, "p_theta"), 
      Q2(reader, "Q2"), W(reader, "W"), Mx(reader, "Mx"), x(reader, "x"), 
      y(reader, "y"), z(reader, "z"), pT(reader, "pT"), xF(reader, "xF"), 
      target_pol(reader, "target_pol") {}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
        bool goodEvent = false;
        string property = binNames[currentFits];

        if (property == "xF") {
            goodEvent = *Q2 > 1 && *W > 2 && *Mx > 1.4 && *y < 0.75;
        }
        else if (property == "Mx") {
            goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
        }
        else if (property == "Q2bin") {
            goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *x>0.2 && *x<0.3 && *pT>0.25 && 
            *pT<0.35 && *xF<0;
        }
        else if (property == "PTTFR" || property ==  "xTFR" || property == "zetaTFR" || 
          property == "Q2TFR" || property ==  "x") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *xF<0;
        }
        else if (property == "PTCFR" || property == "xCFR" || property == "zetaCFR" ||
          property == "Q2TFR") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *xF>0;
        } 
        //
        // epi+X
        else if (property == "xpip") { 
            goodEvent = *Q2>1 && *W>2 && *y<0.75 && *z>0.20;
        }
        else if (property == "PTTFRpip" || property ==  "xTFRpip" || property == "zTFRpip" || 
          property == "Q2TFRpip") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *z>0.20 && *xF<0;
        }
        else if (property == "PTCFRpip" || property == "xCFRpip" || property == "zCFRpip" ||
          property == "Q2TFRpip") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *z>0.20 && *xF>0;
        } 

        else if (property == "UURCy1clasdis") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 && 
            *y>0.45 && *y<0.55;
        } else if (property == "UURCy1clasdis_noLU") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.45 && *y<0.55;
        }else if (property == "UURCy1claspyth") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.45 && *y<0.55;
        } else if (property == "UURCy2clasdis") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.55 && *y<0.65;
        } else if (property == "UURCy2claspyth") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.55 && *y<0.65;
        } else if (property == "UURCy3clasdis") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.65 && *y<0.75;
        } else if (property == "UURCy3claspyth") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.65 && *y<0.75;
        } else if (property == "UURCy3z2clasdis") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.65 && *y<0.75 && *z>0.2 && *z<0.24 && *Q2>2.0 && *Q2<2.5;
        } else if (property == "UURCy3z2clasdis_noLU") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.65 && *y<0.75 && *z>0.2 && *z<0.24 && *Q2>2.0 && *Q2<2.5;
        }else if (property == "UURCy3z2claspyth") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.65 && *y<0.75 && *z>0.2 && *z<0.24 && *Q2>2.0 && *Q2<2.5;
        } else if (property == "UURCy3z5clasdis") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.65 && *y<0.75 && *z>0.40 && *z<0.73 && *Q2>2.0 && *Q2<2.5;
        } else if (property == "UURCy3z5clasdis_noLU") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.65 && *y<0.75 && *z>0.40 && *z<0.73 && *Q2>2.0 && *Q2<2.5;
        } else if (property == "UURCy3z5claspyth") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF > 0 && *Mx > 1.5 && *z>0.2 &&
            *y>0.65 && *y<0.75 && *z>0.40 && *z<0.73 && *Q2>2.0 && *Q2<2.5;
        }
        //
        // epi-X
        else if (property == "xpim") { 
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *z>0.20 && *xF>0;
        }
        else if (property == "PTTFRpim" || property ==  "xTFRpim" || property == "zTFRpim" || 
          property == "Q2TFRpim") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *z>0.20 && *xF<0;
        }
        else if (property == "PTCFRpim" || property == "xCFRpim" || property == "zCFRpim" ||
          property == "Q2TFRpim") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *z>0.20 && *xF>0;
        }
        else {
          std::cout << "Property, " << property << ", not detected." << std::endl;
        }

        if (isMC || (*runnum < 16042 || *runnum > 17811)) {
          return goodEvent;
        } else {
          return goodEvent && *target_pol!=0;
        }
        return false;
    }