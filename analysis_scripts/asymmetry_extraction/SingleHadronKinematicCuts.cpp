#include "SingleHadronKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>
#include "TMath.h"

using std::string;

SingleHadronKinematicCuts::SingleHadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Call to the BaseKinematicCuts constructor
      runnum(reader, "runnum"), fiducial_status(reader, "fiducial_status"), 
      e_phi(reader, "e_phi"), vz_e(reader, "vz_e"),
      Q2(reader, "Q2"), W(reader, "W"), Mx2(reader, "Mx2"), x(reader, "x"), 
      t(reader, "t"), tmin(reader, "tmin"), y(reader, "y"), z(reader, "z"), 
      xi(reader, "xi"), pT(reader, "pT"), xF(reader, "xF"), phi(reader, "phi"), 
      target_pol(reader, "target_pol") {}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    bool checked = false;
    string property = binNames[currentFits];

    if (*fiducial_status != 2) return false; // fiducial cuts
    // if (*Q2 < 1) return false;
    // if (*W < 2) return false;
    // if (*y > 0.8) return false;
    if (property == "W") {
      goodEvent = true;
      return goodEvent;
    }

    if (property == "integrated") {
      goodEvent = *Q2 > 1 && *W > 2 && *Mx2 > 1.8225 && *y < 0.80;
      return goodEvent;
    }
    if (property == "Mx2") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80;
      return goodEvent;
    }
    if (property == "xF") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225;
      return goodEvent;
    } 
    if (property == "xFsmallPT") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && *pT < 0.5;
      return goodEvent;
    } 
    if (property == "xFlargePT") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && *pT > 0.5;
      return goodEvent;
    } 
    if (property == "xTFR" || property == "xi" || property == "PTTFR") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && *xF < 0;
      return goodEvent;
    }
    if (property == "xTFRsmallPT" || property == "xismallPT") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && *xF < 0 && *pT < 0.5;
      return goodEvent;
    }
    if (property == "xTFRlargePT" || property == "xilargePT") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && *xF < 0 && *pT > 0.5;
      return goodEvent;
    }
    if (property == "xCFR" || property == "z" || property == "PTCFR") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && *xF > 0.0;
      return goodEvent;
    }

    if (property == "xFsector1") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*e_phi < 0.7 || *e_phi > 5.9);
      return goodEvent;
    } 
    if (property == "xFsector2") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*e_phi > 0.7 && *e_phi < 1.8);
      return goodEvent;
    } 
    if (property == "xFsector3") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*e_phi > 1.8 && *e_phi < 2.8);
      return goodEvent;
    } 
    if (property == "xFsector4") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*e_phi > 2.8 && *e_phi < 3.8);
      return goodEvent;
    } 
    if (property == "xFsector5") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*e_phi > 3.8 && *e_phi < 4.8);
      return goodEvent;
    } 
    if (property == "xFsector6") {
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*e_phi > 4.8 && *e_phi < 5.85);
      return goodEvent;
    } 

    
    if (property == "pT1xi1x1Q21") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.063210 <= *x && 0.140320 > *x && *Q2 > 1.000070 && *Q2 < 1.914260;
        return goodEvent;
    }

    if (property == "pT1xi1x1Q22") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.063210 <= *x && 0.140320 > *x && *Q2 > 1.914260 && *Q2 < 2.639900;
        return goodEvent;
    }

    if (property == "pT1xi1x1Q23") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.063210 <= *x && 0.140320 > *x && *Q2 > 2.639900 && *Q2 < 10.285940;
        return goodEvent;
    }

    if (property == "pT1xi1x2Q21") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.140320 <= *x && 0.184680 > *x && *Q2 > 1.000020 && *Q2 < 1.916020;
        return goodEvent;
    }

    if (property == "pT1xi1x2Q22") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.140320 <= *x && 0.184680 > *x && *Q2 > 1.916020 && *Q2 < 2.641230;
        return goodEvent;
    }

    if (property == "pT1xi1x2Q23") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.140320 <= *x && 0.184680 > *x && *Q2 > 2.641230 && *Q2 < 10.004930;
        return goodEvent;
    }

    if (property == "pT1xi1x3Q21") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.184680 <= *x && 0.245720 > *x && *Q2 > 1.000130 && *Q2 < 1.913030;
        return goodEvent;
    }

    if (property == "pT1xi1x3Q22") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.184680 <= *x && 0.245720 > *x && *Q2 > 1.913030 && *Q2 < 2.636450;
        return goodEvent;
    }

    if (property == "pT1xi1x3Q23") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.184680 <= *x && 0.245720 > *x && *Q2 > 2.636450 && *Q2 < 9.941760;
        return goodEvent;
    }

    if (property == "pT1xi1x4Q21") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.245720 <= *x && 0.687810 > *x && *Q2 > 1.000030 && *Q2 < 1.909520;
        return goodEvent;
    }

    if (property == "pT1xi1x4Q22") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.245720 <= *x && 0.687810 > *x && *Q2 > 1.909520 && *Q2 < 2.633690;
        return goodEvent;
    }

    if (property == "pT1xi1x4Q23") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && -0.139410 <= *xi && 0.388150 > *xi && 0.245720 <= *x && 0.687810 > *x && *Q2 > 2.633690 && *Q2 < 10.165910;
        return goodEvent;
    }

    if (property == "pT1xi2x1Q21") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.063280 <= *x && 0.140330 > *x && *Q2 > 1.000220 && *Q2 < 1.913900;
        return goodEvent;
    }

    if (property == "pT1xi2x1Q22") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.063280 <= *x && 0.140330 > *x && *Q2 > 1.913900 && *Q2 < 2.637110;
        return goodEvent;
    }

    if (property == "pT1xi2x1Q23") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.063280 <= *x && 0.140330 > *x && *Q2 > 2.637110 && *Q2 < 9.825690;
        return goodEvent;
    }

    if (property == "pT1xi2x2Q21") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.140330 <= *x && 0.184630 > *x && *Q2 > 1.000010 && *Q2 < 1.914370;
        return goodEvent;
    }

    if (property == "pT1xi2x2Q22") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.140330 <= *x && 0.184630 > *x && *Q2 > 1.914370 && *Q2 < 2.639730;
        return goodEvent;
    }

    if (property == "pT1xi2x2Q23") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.140330 <= *x && 0.184630 > *x && *Q2 > 2.639730 && *Q2 < 10.311350;
        return goodEvent;
    }

    if (property == "pT1xi2x3Q21") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.184630 <= *x && 0.245730 > *x && *Q2 > 1.000100 && *Q2 < 1.912350;
        return goodEvent;
    }

    if (property == "pT1xi2x3Q22") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.184630 <= *x && 0.245730 > *x && *Q2 > 1.912350 && *Q2 < 2.636010;
        return goodEvent;
    }

    if (property == "pT1xi2x3Q23") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.184630 <= *x && 0.245730 > *x && *Q2 > 2.636010 && *Q2 < 10.111230;
        return goodEvent;
    }

    if (property == "pT1xi2x4Q21") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.245730 <= *x && 0.695230 > *x && *Q2 > 1.000040 && *Q2 < 1.914880;
        return goodEvent;
    }

    if (property == "pT1xi2x4Q22") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.245730 <= *x && 0.695230 > *x && *Q2 > 1.914880 && *Q2 < 2.639570;
        return goodEvent;
    }

    if (property == "pT1xi2x4Q23") {
        goodEvent = 0.000150 <= *pT && 0.412050 > *pT && 0.388150 <= *xi && 0.808840 > *xi && 0.245730 <= *x && 0.695230 > *x && *Q2 > 2.639570 && *Q2 < 10.224480;
        return goodEvent;
    }

    if (property == "pT2xi1x1Q21") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.063270 <= *x && 0.140200 > *x && *Q2 > 1.000020 && *Q2 < 1.911420;
        return goodEvent;
    }

    if (property == "pT2xi1x1Q22") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.063270 <= *x && 0.140200 > *x && *Q2 > 1.911420 && *Q2 < 2.635120;
        return goodEvent;
    }

    if (property == "pT2xi1x1Q23") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.063270 <= *x && 0.140200 > *x && *Q2 > 2.635120 && *Q2 < 10.086440;
        return goodEvent;
    }

    if (property == "pT2xi1x2Q21") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.140200 <= *x && 0.184460 > *x && *Q2 > 1.000030 && *Q2 < 1.912560;
        return goodEvent;
    }

    if (property == "pT2xi1x2Q22") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.140200 <= *x && 0.184460 > *x && *Q2 > 1.912560 && *Q2 < 2.636340;
        return goodEvent;
    }

    if (property == "pT2xi1x2Q23") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.140200 <= *x && 0.184460 > *x && *Q2 > 2.636340 && *Q2 < 10.140130;
        return goodEvent;
    }

    if (property == "pT2xi1x3Q21") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.184460 <= *x && 0.245560 > *x && *Q2 > 1.000070 && *Q2 < 1.910100;
        return goodEvent;
    }

    if (property == "pT2xi1x3Q22") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.184460 <= *x && 0.245560 > *x && *Q2 > 1.910100 && *Q2 < 2.634490;
        return goodEvent;
    }

    if (property == "pT2xi1x3Q23") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.184460 <= *x && 0.245560 > *x && *Q2 > 2.634490 && *Q2 < 10.039920;
        return goodEvent;
    }

    if (property == "pT2xi1x4Q21") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.245560 <= *x && 0.686600 > *x && *Q2 > 1.000020 && *Q2 < 1.915670;
        return goodEvent;
    }

    if (property == "pT2xi1x4Q22") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.245560 <= *x && 0.686600 > *x && *Q2 > 1.915670 && *Q2 < 2.641450;
        return goodEvent;
    }

    if (property == "pT2xi1x4Q23") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && -0.143920 <= *xi && 0.388030 > *xi && 0.245560 <= *x && 0.686600 > *x && *Q2 > 2.641450 && *Q2 < 10.113240;
        return goodEvent;
    }

    if (property == "pT2xi2x1Q21") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.063220 <= *x && 0.140470 > *x && *Q2 > 1.000090 && *Q2 < 1.911970;
        return goodEvent;
    }

    if (property == "pT2xi2x1Q22") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.063220 <= *x && 0.140470 > *x && *Q2 > 1.911970 && *Q2 < 2.638360;
        return goodEvent;
    }

    if (property == "pT2xi2x1Q23") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.063220 <= *x && 0.140470 > *x && *Q2 > 2.638360 && *Q2 < 9.993850;
        return goodEvent;
    }

    if (property == "pT2xi2x2Q21") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.140470 <= *x && 0.184810 > *x && *Q2 > 1.000030 && *Q2 < 1.915860;
        return goodEvent;
    }

    if (property == "pT2xi2x2Q22") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.140470 <= *x && 0.184810 > *x && *Q2 > 1.915860 && *Q2 < 2.640410;
        return goodEvent;
    }

    if (property == "pT2xi2x2Q23") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.140470 <= *x && 0.184810 > *x && *Q2 > 2.640410 && *Q2 < 10.054770;
        return goodEvent;
    }

    if (property == "pT2xi2x3Q21") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.184810 <= *x && 0.246040 > *x && *Q2 > 1.000200 && *Q2 < 1.915710;
        return goodEvent;
    }

    if (property == "pT2xi2x3Q22") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.184810 <= *x && 0.246040 > *x && *Q2 > 1.915710 && *Q2 < 2.638880;
        return goodEvent;
    }

    if (property == "pT2xi2x3Q23") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.184810 <= *x && 0.246040 > *x && *Q2 > 2.638880 && *Q2 < 9.914960;
        return goodEvent;
    }

    if (property == "pT2xi2x4Q21") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.246040 <= *x && 0.689050 > *x && *Q2 > 1.000020 && *Q2 < 1.915620;
        return goodEvent;
    }

    if (property == "pT2xi2x4Q22") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.246040 <= *x && 0.689050 > *x && *Q2 > 1.915620 && *Q2 < 2.641890;
        return goodEvent;
    }

    if (property == "pT2xi2x4Q23") {
        goodEvent = 0.412050 <= *pT && 1.569350 > *pT && 0.388030 <= *xi && 0.804170 > *xi && 0.246040 <= *x && 0.689050 > *x && *Q2 > 2.641890 && *Q2 < 10.078500;
        return goodEvent;
    }


  if (isMC || (*runnum < 16042 || *runnum > 17811)) {
    return goodEvent;
  } else {
    // return goodEvent && *target_pol!=0;
    return goodEvent;
  }
  return false;
}