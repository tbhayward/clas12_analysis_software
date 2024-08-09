#include "SingleHadronKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>
#include "TMath.h"

using std::string;

SingleHadronKinematicCuts::SingleHadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Call to the BaseKinematicCuts constructor
      runnum(reader, "runnum"), 
      e_theta(reader, "e_theta"), e_phi(reader, "e_phi"),
      p_p(reader, "p_p"), p_theta(reader, "p_theta"), p_phi(reader, "p_phi"), 
      Q2(reader, "Q2"), W(reader, "W"), Mx(reader, "Mx"), x(reader, "x"), 
      y(reader, "y"), z(reader, "z"), pT(reader, "pT"), xF(reader, "xF"),
      phi(reader, "phi"), 
      target_pol(reader, "target_pol") {}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
        bool goodEvent = false;
        bool checked = false;
        string property = binNames[currentFits];


        if (property == "xF" || "x" || "PT") {
            goodEvent = *Q2 > 1 && *W > 2 && *Mx > 1.4 && *y < 0.75;
        }
        if (property == "Mx") {
            goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
        }

        // multidimensional analysis checks
        else if (*Q2>1 && *W>2 && *Mx>1.4  && *y<0.75) {
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


        if (isMC || (*runnum < 16042 || *runnum > 17811)) {
          return goodEvent;
        } else {
          // return goodEvent && *target_pol!=0;
          return goodEvent;
        }
        return false;
    }