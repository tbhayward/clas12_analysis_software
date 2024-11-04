#include "SingleHadronKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>
#include "TMath.h"

using std::string;

SingleHadronKinematicCuts::SingleHadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Call to the BaseKinematicCuts constructor
      runnum(reader, "runnum"), fiducial_status(reader, "fiducial_status"), 
      e_phi(reader, "e_phi"), vz_e(reader, "vz_e"), p_p(reader, "p_p"),  vz_p(reader, "vz_p"), 
      Q2(reader, "Q2"), W(reader, "W"), Mx2(reader, "Mx2"), x(reader, "x"), 
      t(reader, "t"), tmin(reader, "tmin"), y(reader, "y"), z(reader, "z"), 
      xi(reader, "xi"), pT(reader, "pT"), xF(reader, "xF"), phi(reader, "phi"), 
      target_pol(reader, "target_pol") {}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    bool checked = false;
    string property = binNames[currentFits];

    if (*fiducial_status != 2) return false; // fiducial cuts
    if (*Q2 < 1) return false;
    if (*W < 2) return false;
    if (*y > 0.8) return false;

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
      goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && *xF > 0.2;
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
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.063210 <= *x && 0.143940 > *x && *Q2 > 1.000070 && *Q2 < 1.932240;
        return goodEvent;
    }

    if (property == "pT1xi1x1Q22") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.063210 <= *x && 0.143940 > *x && *Q2 > 1.932240 && *Q2 < 2.673890;
        return goodEvent;
    }

    if (property == "pT1xi1x1Q23") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.063210 <= *x && 0.143940 > *x && *Q2 > 2.673890 && *Q2 < 10.285940;
        return goodEvent;
    }

    if (property == "pT1xi1x2Q21") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.143940 <= *x && 0.190670 > *x && *Q2 > 1.000020 && *Q2 < 1.933920;
        return goodEvent;
    }

    if (property == "pT1xi1x2Q22") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.143940 <= *x && 0.190670 > *x && *Q2 > 1.933920 && *Q2 < 2.675150;
        return goodEvent;
    }

    if (property == "pT1xi1x2Q23") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.143940 <= *x && 0.190670 > *x && *Q2 > 2.675150 && *Q2 < 10.004930;
        return goodEvent;
    }

    if (property == "pT1xi1x3Q21") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.190670 <= *x && 0.254060 > *x && *Q2 > 1.000130 && *Q2 < 1.931030;
        return goodEvent;
    }

    if (property == "pT1xi1x3Q22") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.190670 <= *x && 0.254060 > *x && *Q2 > 1.931030 && *Q2 < 2.670410;
        return goodEvent;
    }

    if (property == "pT1xi1x3Q23") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.190670 <= *x && 0.254060 > *x && *Q2 > 2.670410 && *Q2 < 10.032720;
        return goodEvent;
    }

    if (property == "pT1xi1x4Q21") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.254060 <= *x && 0.689540 > *x && *Q2 > 1.000030 && *Q2 < 1.927430;
        return goodEvent;
    }

    if (property == "pT1xi1x4Q22") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.254060 <= *x && 0.689540 > *x && *Q2 > 1.927430 && *Q2 < 2.667740;
        return goodEvent;
    }

    if (property == "pT1xi1x4Q23") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && -0.301350 <= *xi && 0.344580 > *xi && 0.254060 <= *x && 0.689540 > *x && *Q2 > 2.667740 && *Q2 < 10.165910;
        return goodEvent;
    }

    if (property == "pT1xi2x1Q21") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.063280 <= *x && 0.143950 > *x && *Q2 > 1.000220 && *Q2 < 1.931490;
        return goodEvent;
    }

    if (property == "pT1xi2x1Q22") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.063280 <= *x && 0.143950 > *x && *Q2 > 1.931490 && *Q2 < 2.671350;
        return goodEvent;
    }

    if (property == "pT1xi2x1Q23") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.063280 <= *x && 0.143950 > *x && *Q2 > 2.671350 && *Q2 < 9.987370;
        return goodEvent;
    }

    if (property == "pT1xi2x2Q21") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.143950 <= *x && 0.190640 > *x && *Q2 > 1.000010 && *Q2 < 1.932540;
        return goodEvent;
    }

    if (property == "pT1xi2x2Q22") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.143950 <= *x && 0.190640 > *x && *Q2 > 1.932540 && *Q2 < 2.673890;
        return goodEvent;
    }

    if (property == "pT1xi2x2Q23") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.143950 <= *x && 0.190640 > *x && *Q2 > 2.673890 && *Q2 < 10.311350;
        return goodEvent;
    }

    if (property == "pT1xi2x3Q21") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.190640 <= *x && 0.254090 > *x && *Q2 > 1.000080 && *Q2 < 1.930140;
        return goodEvent;
    }

    if (property == "pT1xi2x3Q22") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.190640 <= *x && 0.254090 > *x && *Q2 > 1.930140 && *Q2 < 2.669880;
        return goodEvent;
    }

    if (property == "pT1xi2x3Q23") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.190640 <= *x && 0.254090 > *x && *Q2 > 2.669880 && *Q2 < 10.111230;
        return goodEvent;
    }

    if (property == "pT1xi2x4Q21") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.254090 <= *x && 0.698890 > *x && *Q2 > 1.000040 && *Q2 < 1.932940;
        return goodEvent;
    }

    if (property == "pT1xi2x4Q22") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.254090 <= *x && 0.698890 > *x && *Q2 > 1.932940 && *Q2 < 2.674500;
        return goodEvent;
    }

    if (property == "pT1xi2x4Q23") {
        goodEvent = 0.000000 <= *pT && 0.408820 > *pT && 0.344580 <= *xi && 0.808840 > *xi && 0.254090 <= *x && 0.698890 > *x && *Q2 > 2.674500 && *Q2 < 10.224480;
        return goodEvent;
    }

    if (property == "pT2xi1x1Q21") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.063250 <= *x && 0.143820 > *x && *Q2 > 1.000020 && *Q2 < 1.929410;
        return goodEvent;
    }

    if (property == "pT2xi1x1Q22") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.063250 <= *x && 0.143820 > *x && *Q2 > 1.929410 && *Q2 < 2.668970;
        return goodEvent;
    }

    if (property == "pT2xi1x1Q23") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.063250 <= *x && 0.143820 > *x && *Q2 > 2.668970 && *Q2 < 10.086440;
        return goodEvent;
    }

    if (property == "pT2xi1x2Q21") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.143820 <= *x && 0.190460 > *x && *Q2 > 1.000030 && *Q2 < 1.930690;
        return goodEvent;
    }

    if (property == "pT2xi1x2Q22") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.143820 <= *x && 0.190460 > *x && *Q2 > 1.930690 && *Q2 < 2.670430;
        return goodEvent;
    }

    if (property == "pT2xi1x2Q23") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.143820 <= *x && 0.190460 > *x && *Q2 > 2.670430 && *Q2 < 10.140130;
        return goodEvent;
    }

    if (property == "pT2xi1x3Q21") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.190460 <= *x && 0.253870 > *x && *Q2 > 1.000070 && *Q2 < 1.928200;
        return goodEvent;
    }

    if (property == "pT2xi1x3Q22") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.190460 <= *x && 0.253870 > *x && *Q2 > 1.928200 && *Q2 < 2.668540;
        return goodEvent;
    }

    if (property == "pT2xi1x3Q23") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.190460 <= *x && 0.253870 > *x && *Q2 > 2.668540 && *Q2 < 10.039920;
        return goodEvent;
    }

    if (property == "pT2xi1x4Q21") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.253870 <= *x && 0.693940 > *x && *Q2 > 1.000020 && *Q2 < 1.933840;
        return goodEvent;
    }

    if (property == "pT2xi1x4Q22") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.253870 <= *x && 0.693940 > *x && *Q2 > 1.933840 && *Q2 < 2.675720;
        return goodEvent;
    }

    if (property == "pT2xi1x4Q23") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && -0.305970 <= *xi && 0.344480 > *xi && 0.253870 <= *x && 0.693940 > *x && *Q2 > 2.675720 && *Q2 < 10.113240;
        return goodEvent;
    }

    if (property == "pT2xi2x1Q21") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.063220 <= *x && 0.144090 > *x && *Q2 > 1.000090 && *Q2 < 1.929750;
        return goodEvent;
    }

    if (property == "pT2xi2x1Q22") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.063220 <= *x && 0.144090 > *x && *Q2 > 1.929750 && *Q2 < 2.672190;
        return goodEvent;
    }

    if (property == "pT2xi2x1Q23") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.063220 <= *x && 0.144090 > *x && *Q2 > 2.672190 && *Q2 < 10.073840;
        return goodEvent;
    }

    if (property == "pT2xi2x2Q21") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.144090 <= *x && 0.190820 > *x && *Q2 > 1.000030 && *Q2 < 1.933660;
        return goodEvent;
    }

    if (property == "pT2xi2x2Q22") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.144090 <= *x && 0.190820 > *x && *Q2 > 1.933660 && *Q2 < 2.674890;
        return goodEvent;
    }

    if (property == "pT2xi2x2Q23") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.144090 <= *x && 0.190820 > *x && *Q2 > 2.674890 && *Q2 < 10.153080;
        return goodEvent;
    }

    if (property == "pT2xi2x3Q21") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.190820 <= *x && 0.254390 > *x && *Q2 > 1.000200 && *Q2 < 1.933500;
        return goodEvent;
    }

    if (property == "pT2xi2x3Q22") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.190820 <= *x && 0.254390 > *x && *Q2 > 1.933500 && *Q2 < 2.673150;
        return goodEvent;
    }

    if (property == "pT2xi2x3Q23") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.190820 <= *x && 0.254390 > *x && *Q2 > 2.673150 && *Q2 < 9.983770;
        return goodEvent;
    }

    if (property == "pT2xi2x4Q21") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.254390 <= *x && 0.689050 > *x && *Q2 > 1.000020 && *Q2 < 1.933970;
        return goodEvent;
    }

    if (property == "pT2xi2x4Q22") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.254390 <= *x && 0.689050 > *x && *Q2 > 1.933970 && *Q2 < 2.676190;
        return goodEvent;
    }

    if (property == "pT2xi2x4Q23") {
        goodEvent = 0.408820 <= *pT && 1.575770 > *pT && 0.344480 <= *xi && 0.804170 > *xi && 0.254390 <= *x && 0.689050 > *x && *Q2 > 2.676190 && *Q2 < 10.126660;
        return goodEvent;
    }
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