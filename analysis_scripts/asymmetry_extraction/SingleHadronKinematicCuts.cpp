#include "SingleHadronKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>
#include "TMath.h"

using std::string;

SingleHadronKinematicCuts::SingleHadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Call to the BaseKinematicCuts constructor
      runnum(reader, "runnum"), 
      e_theta(reader, "e_theta"), e_phi(reader, "e_phi"), vz_e(reader, "vz_e"),
      p_p(reader, "p_p"), p_theta(reader, "p_theta"), p_phi(reader, "p_phi"), 
      vz_p(reader, "vz_p"), 
      Q2(reader, "Q2"), W(reader, "W"), Mx(reader, "Mx"), x(reader, "x"), 
      y(reader, "y"), z(reader, "z"), pT(reader, "pT"), xF(reader, "xF"),
      phi(reader, "phi"), 
      target_pol(reader, "target_pol") {}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    bool checked = false;
    string property = binNames[currentFits];

    if (property == "xF" || property == "x" || property == "z" || property == "PT" || property == "runnum") {
      goodEvent = *vz_e > -10 && *vz_e < 1 && *vz_p > -10 && *vz_p < 1 && *Q2 > 1 && *W > 2 && *Mx > 1.35 && *y < 0.75;
      goodEvent = goodEvent && *x > 0.06 && *x < 0.60 && *pT > 0 && *pT < 1.2 && *xF > -1 && *xF < 1;
      return goodEvent;
    } else if (property == "Mx") {
      goodEvent = *vz_e > -10 && *vz_e < 1 && *vz_p > -10 && *vz_p < 1 && *Q2 > 1 && *W > 2 && *y < 0.75 && *Mx > 0;
      goodEvent = goodEvent && *x > 0.06 && *x < 0.60 && *pT > 0 && *pT < 1.2 && *xF > -1 && *xF < 1;
      return goodEvent;
    }
    if (property == "xFall" || property == "xall" || property == "zall" || property == "PTall") {
      goodEvent = *vz_e > -10 && *vz_e < 1 && *vz_p > -10 && *vz_p < 1 && *Q2 > 1 && *W > 2 && *Mx > 0 && *y < 0.75;
      goodEvent = goodEvent && *x > 0.06 && *x < 0.60 && *pT > 0 && *pT < 1.2 && *xF > -1 && *xF < 1;
      return goodEvent;
    } 
    if (property == "Q2multi1") {
      goodEvent = *vz_e > -10 && *vz_e < 1 && *vz_p > -10 && *vz_p < 1 && *Q2 > 1 && *W > 2 && *Mx > 0.95 && *y < 0.75;
      goodEvent = goodEvent && *x > 0.06 && *x < 0.60 && *pT > 0 && *pT < 1.2 && *xF > -1 && *xF < 1;
      goodEvent = goodEvent && *x > 0.12 && *x < 0.15 && *pT > 0.30 && *pT < 0.5 && *z > 0.16 && *z < 0.24;
      return goodEvent;
    }
    if (property == "Q2multi2") {
      goodEvent = *vz_e > -10 && *vz_e < 1 && *vz_p > -10 && *vz_p < 1 && *Q2 > 1 && *W > 2 && *Mx > 0.95 && *y < 0.75;
      goodEvent = goodEvent && *x > 0.06 && *x < 0.60 && *pT > 0 && *pT < 1.2 && *xF > -1 && *xF < 1;
      goodEvent = goodEvent && *x > 0.15 && *x < 0.18 && *pT > 0.30 && *pT < 0.5 && *z > 0.16 && *z < 0.22 ;
      return goodEvent;
    }
    if (property == "Q2multi3 ") {
      goodEvent = *vz_e > -10 && *vz_e < 1 && *vz_p > -10 && *vz_p < 1 && *Q2 > 1 && *W > 2 && *Mx > 0.95 && *y < 0.75;
      goodEvent = goodEvent && *x > 0.06 && *x < 0.60 && *pT > 0 && *pT < 1.2 && *xF > -1 && *xF < 1;
      goodEvent = goodEvent && *x > 0.18 && *x < 0.21 && *pT > 0.30 && *pT < 0.5 && *z > 0.16 && *z < 0.22;
      return goodEvent;
    }

    if (*Q2 > 1 && *W > 2 && *Mx > 1.35 && *y < 0.75 && !checked) {
      goodEvent = *vz_e > -10 && *vz_e < 1 && *vz_p > -10 && *vz_p < 1;
      size_t pos = property.find("z");
      std::string prez = property.substr(0, pos);
      std::string postz = property.substr(pos);

      bool prezGood = false, zGood = false;

      // Q2-x bin checks
      if (prez == "Q2x1") {
          prezGood = *x < 0.1;
      } else if (prez == "Q2x2") {
          prezGood = *Q2 < 1.50 && *x > 0.1 && *x < 0.14;
      } else if (prez == "Q2x3") {
          prezGood = *Q2 >= 1.50 && *Q2 < 1.70 && *x > 0.1 && *x < 0.14;
      } else if (prez == "Q2x4") {
          prezGood = *Q2 >= 1.70 && *x > 0.1 && *x < 0.14;
      } else if (prez == "Q2x5") {
          prezGood = *Q2 < 1.50 && *x > 0.14 && *x < 0.21;
      } else if (prez == "Q2x6") {
          prezGood = *Q2 >= 1.50 && *Q2 < 1.70 && *x > 0.14 && *x < 0.21;
      } else if (prez == "Q2x7") {
          prezGood = *Q2 >= 1.70 && *Q2 < 2.00 && *x > 0.14 && *x < 0.21;
      } else if (prez == "Q2x8") {
          prezGood = *Q2 >= 2.00 && *x > 0.14 && *x < 0.21;
      } else if (prez == "Q2x9") {
          prezGood = *Q2 < 2.20 && *x > 0.21 && *x < 0.30;
      } else if (prez == "Q2x10") {
          prezGood = *Q2 >= 2.20 && *Q2 < 2.60 && *x > 0.21 && *x < 0.30;
      } else if (prez == "Q2x11") {
          prezGood = *Q2 >= 2.60 && *x > 0.21 && *x < 0.30;
      } else if (prez == "Q2x12") {
          prezGood = *Q2 < 3.20 && *x > 0.30 && *x < 0.42;
      } else if (prez == "Q2x13") {
          prezGood = *Q2 >= 3.20 && *x > 0.30 && *x < 0.42;
      } else if (prez == "Q2x14") {
          prezGood = *x >= 0.42;
      }

      // z bin checks
      if (postz == "z1") {
          zGood = *z > 0 && *z <= 0.19;
      } else if (postz == "z2") {
          zGood = *z > 0.19 && *z <= 0.30;
      } else if (postz == "z3") {
          zGood = *z > 0.30 && *z <= 0.42;
      } else if (postz == "z4") {
          zGood = *z > 0.42 && *z <= 1.00;
      }

      // Final check combining all criteria
      goodEvent = goodEvent && prezGood && zGood;
  }

  //   if (*Q2 > 1 && *W > 2 && *Mx > 1.35 && *y < 0.75 && !checked) {
  //     goodEvent = *vz_e > -10 && *vz_e < 1 && *vz_p > -10 && *vz_p < 1;
  //     size_t pos = property.find("z");
  //     std::string prez = property.substr(0, pos);
  //     std::string postz = property.substr(pos);

  //     bool prezGood = false, yGood = false, zGood = false;

  //     if (prez == "Q2y1" || prez == "Q2y2" || prez == "Q2y3" || prez == "Q2y4") {
  //         prezGood = *Q2 > 1 && *Q2 <= 2;
  //     } else if (prez == "Q2y5" || prez == "Q2y6" || prez == "Q2y7" || prez == "Q2y8") {
  //         prezGood = *Q2 > 2 && *Q2 <= 3;
  //     } else if (prez == "Q2y9" || prez == "Q2y10" || prez == "Q2y11" || prez == "Q2y12") {
  //         prezGood = *Q2 > 3 && *Q2 <= 4;
  //     } else if (prez == "Q2y13" || prez == "Q2y14" || prez == "Q2y15") {
  //         prezGood = *Q2 > 4 && *Q2 <= 5;
  //     } else if (prez == "Q2y16" || prez == "Q2y17") {
  //         prezGood = *Q2 > 5 && *Q2 <= 7;
  //     }

  //     if (prezGood) {
  //         if (prez == "Q2y1" || prez == "Q2y5" || prez == "Q2y9" || prez == "Q2y13" || prez == "Q2y16") {
  //             yGood = *y > 0.65 && *y <= 0.75;
  //         } else if (prez == "Q2y2" || prez == "Q2y6" || prez == "Q2y10" || prez == "Q2y14" || prez == "Q2y17") {
  //             yGood = *y > 0.55 && *y <= 0.65;
  //         } else if (prez == "Q2y3" || prez == "Q2y7" || prez == "Q2y11" || prez == "Q2y15") {
  //             yGood = *y > 0.45 && *y <= 0.55;
  //         } else if (prez == "Q2y4" || prez == "Q2y8" || prez == "Q2y12") {
  //             yGood = *y > 0.30 && *y <= 0.45;
  //         }
  //     }

  //     if (yGood) {
  //         if (postz == "z1") {
  //             zGood = *z > 0.10 && *z <= 0.25;
  //         } else if (postz == "z2") {
  //             zGood = *z > 0.25 && *z <= 0.35;
  //         } else if (postz == "z3") {
  //             zGood = *z > 0.35 && *z <= 0.45;
  //         } else if (postz == "z4") {
  //             zGood = *z > 0.45 && *z <= 0.55;
  //         } else if (postz == "z5") {
  //             zGood = *z > 0.55 && *z <= 0.75;
  //         }
  //     }

  //     goodEvent = goodEvent && prezGood && yGood && zGood;
  // }


  if (isMC || (*runnum < 16042 || *runnum > 17811)) {
    return goodEvent;
  } else {
    // return goodEvent && *target_pol!=0;
    return goodEvent;
  }
  return false;
}