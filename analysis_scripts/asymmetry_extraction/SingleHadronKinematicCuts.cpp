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
        bool checked = false;
        string property = binNames[currentFits];

        if (property == "xF") {
            goodEvent = *Q2 > 1 && *W > 2 && *Mx > 1.9 && *y < 0.75;
        }
        else if (property == "Mx") {
            goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
        }
        else if (property == "PTTFR" || property ==  "xTFR" || property == "zetaTFR" || 
          property == "Q2TFR" || property ==  "x") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *xF<0;
        }
        else if (property == "PTCFR" || property == "xCFR" || property == "zetaCFR" ||
          property == "Q2CFR") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *xF>0;
        } 

        // multidimensional analysis checks
        else if (*Q2>1 && *W>2 && *Mx>1.9  && *y<0.75) {
          checked = true;
          size_t pos = property.find("z");
          std::string prez = property.substr(0, pos);
          std::string postz = property.substr(pos); // Skip "z" itself

          // Corrected conditional checks for prez
          if (prez == "Q2y1" || prez == "Q2y2" || prez == "Q2y3"  || prez == "Q2y4") {
              goodEvent = *Q2 > 1 && *Q2 <= 2;
          } else if (prez == "Q2y5" || prez == "Q2y6" || prez == "Q2y7"  || prez == "Q2y8") {
              goodEvent = *Q2 > 2 && *Q2 <= 3;
          } else if (prez == "Q2y9" || prez == "Q2y10" || prez == "Q2y11"  || prez == "Q2y12") {
              goodEvent = *Q2 > 3 && *Q2 <= 4;
          } else if (prez == "Q2y13" || prez == "Q2y14" || prez == "Q2y15") {
              goodEvent = *Q2 > 4 && *Q2 <= 5;
          } else if (prez == "Q2y16" || prez == "Q2y17") {
              goodEvent = *Q2 > 5 && *Q2 <= 7;
          }
          if (goodEvent) {
            // Assuming goodEvent might be true, let's adjust for y ranges
            if ((prez == "Q2y1" || prez == "Q2y5" || prez == "Q2y9" || prez == "Q2y13" || prez == "Q2y16")) {
                goodEvent = *y > 0.65 && *y <= 0.75;
            } else if ((prez == "Q2y2" || prez == "Q2y6" || prez == "Q2y10" || prez == "Q2y14" || prez == "Q2y17")) {
                goodEvent = *y > 0.55 && *y <= 0.65;
            } else if ((prez == "Q2y3" || prez == "Q2y7" || prez == "Q2y11" || prez == "Q2y15")) {
                goodEvent = *y > 0.45 && *y <= 0.55;
            } else if ((prez == "Q2y4" || prez == "Q2y8" || prez == "Q2y12")) {
                goodEvent = *y > 0.30 && *y <= 0.45;
            }

            if (goodEvent) {
              // Corrected z bin checks
              if (postz == "z1") {
                  goodEvent = *z > 0.10 && *z <= 0.25;
              } else if (postz == "z2") {
                  goodEvent = *z > 0.25 && *z <= 0.35;
              } else if (postz == "z3") {
                  goodEvent = *z > 0.35 && *z <= 0.45;
              } else if (postz == "z4") {
                  goodEvent = *z > 0.45 && *z <= 0.55;
              } else if (postz == "z5") {
                  goodEvent = *z > 0.55 && *z <= 0.75;
              }
            }
          }
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
          property == "Q2CFRpip") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *z>0.20 && *xF>0;
        } 

        else if (property == "UURCQ2y1z1") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.15 && *z<0.20;
        } else if (property == "UURCQ2y1z2") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.20 && *z<0.24;
        } else if (property == "UURCQ2y1z3") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.24 && *z<0.29;
        } else if (property == "UURCQ2y1z4") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.29 && *z<0.40;
        } else if (property == "UURCQ2y1z5") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.40 && *z<0.73;
        }

        else if (property == "UURCQ2y2z1") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.18 && *z<0.23;
        } else if (property == "UURCQ2y2z2") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.23 && *z<0.26;
        } else if (property == "UURCQ2y2z3") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.26 && *z<0.31;
        } else if (property == "UURCQ2y2z4") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.31 && *z<0.38;
        } else if (property == "UURCQ2y2z5") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.38 && *z<0.50;
        } else if (property == "UURCQ2y2z6") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.50 && *z<0.74;
        }

        else if (property == "UURCQ2y4z1") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.26 && *z<0.32;
        } else if (property == "UURCQ2y4z2") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.32 && *z<0.37;
        } else if (property == "UURCQ2y4z3") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.37 && *z<0.43;
        } else if (property == "UURCQ2y4z4") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.43 && *z<0.50;
        } else if (property == "UURCQ2y4z5") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.50 && *z<0.60;
        } else if (property == "UURCQ2y4z6") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.000 && *Q2 <= 2.423 && *z>0.60 && *z<0.71;
        } 

        else if (property == "UURCQ2y5z1") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.423 && *Q2 <= 2.987 && *z>0.15 && *z<0.19;
        } else if (property == "UURCQ2y5z2") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.423 && *Q2 <= 2.987 && *z>0.19 && *z<0.24;
        } else if (property == "UURCQ2y5z3") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.423 && *Q2 <= 2.987 && *z>0.24 && *z<0.29;
        } else if (property == "UURCQ2y5z4") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.423 && *Q2 <= 2.987 && *z>0.29 && *z<0.38;
        } else if (property == "UURCQ2y5z5") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.423 && *Q2 <= 2.987 && *z>0.38 && *z<0.50;
        } else if (property == "UURCQ2y5z6") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.30 && *y<0.45 && *Q2 > 2.423 && *Q2 <= 2.987 && *z>0.50 && *z<0.73;
        } 

        else if (property == "UURCQ2y16z1") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 5.384 && *Q2 <= 9.896 && *z>0.15 && *z<0.20;
        } else if (property == "UURCQ2y16z2") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 5.384 && *Q2 <= 9.896 && *z>0.20 && *z<0.25;
        } else if (property == "UURCQ2y16z3") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 5.384 && *Q2 <= 9.896 && *z>0.25 && *z<0.32;
        } else if (property == "UURCQ2y16z4") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 5.384 && *Q2 <= 9.896 && *z>0.32 && *z<0.41;
        } else if (property == "UURCQ2y16z5") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.65 && *y<0.75 && *Q2 > 5.384 && *Q2 <= 9.896 && *z>0.41 && *z<0.71;
        }

        else if (property == "UURCQ2y17z1") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 5.384 && *Q2 <= 7.922 && *z>0.18 && *z<0.23;
        } else if (property == "UURCQ2y17z2") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 5.384 && *Q2 <= 7.922 && *z>0.23 && *z<0.30;
        } else if (property == "UURCQ2y17z3") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 5.384 && *Q2 <= 7.922 && *z>0.30 && *z<0.38;
        } else if (property == "UURCQ2y17z4") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 5.384 && *Q2 <= 7.922 && *z>0.38 && *z<0.48;
        } else if (property == "UURCQ2y17z5") {
          goodEvent = *Q2>1 && *W>2 && *y<0.75 && *xF>0 && *Mx>1.5 && *p_p > 1.25 &&
            *y>0.55 && *y<0.65 && *Q2 > 5.384 && *Q2 <= 7.922 && *z>0.48 && *z<0.72;
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






        // else if (checked == false) {
        //   std::cout << "Property, " << property << ", not detected." << std::endl;
        // }

        if (isMC || (*runnum < 16042 || *runnum > 17811)) {
          return goodEvent;
        } else {
          return goodEvent && *target_pol!=0;
        }
        return false;
    }