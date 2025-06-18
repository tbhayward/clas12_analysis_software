#include "SingleHadronKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>
#include "TMath.h"

using std::string;

// Physical masses (GeV)
static constexpr double m_e  = 0.000511;  // electron
static constexpr double m_pi = 0.13957;   // charged pion

//================================================================================
// Constructor: grab every branch we’ll need
//================================================================================
SingleHadronKinematicCuts::SingleHadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader),
      runnum       (reader, "runnum"),
      fiducial_status(reader, "fiducial_status"),

      // Electron‐side branches (added e_p, e_theta)
      e_p          (reader, "e_p"),
      e_theta      (reader, "e_theta"),
      e_phi        (reader, "e_phi"),

      // Pion‐side branches (added p_theta)
      p_p          (reader, "p_p"),
      p_theta      (reader, "p_theta"),
      p_phi        (reader, "p_phi"),

      // Standard DIS / hadron variables
      Q2           (reader, "Q2"),
      W            (reader, "W"),
      Mx2          (reader, "Mx2"),
      xF           (reader, "xF"),
      pT           (reader, "pT"),
      y            (reader, "y"),
      x            (reader, "x"),
      xi           (reader, "xi"),
      phi          (reader, "phi"),
      z            (reader, "z"),
      t            (reader, "t"),
      tmin         (reader, "tmin"),
      target_pol   (reader, "target_pol")
{}

//================================================================================
// beamEnergy(run): 
//    Return beam energy (GeV) based on run number.  Matches the mapping:
//      • 6616–6783   → Eb = 10.1998 (RGA Sp19, H₂ data)
//      • 16042–17065 → Eb = 10.5473 (RGC Su22)
//      • 17067–17724 → Eb = 10.5563 (RGC Fa22)
//      • 17725–17811 → Eb = 10.5593 (RGC Sp23)
//    Outside these → 0.0  (will cause t‐calc to be nonsense and fail).
//================================================================================
static double beamEnergy(int run)
{
    if (run >= 6616  && run <= 6783)   return 10.1998;
    if (run >= 16042 && run <= 17065)  return 10.5473;
    if (run >= 17067 && run <= 17724)  return 10.5563;
    if (run >= 17725 && run <= 17811)  return 10.5593;
    return 0.0;
}

//================================================================================
// compute_t(…) 
//    Given arrays of runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi (all scalars
//    for one event), return t = (q – p_pi)^2 = (ΔE)^2 – (Δp)^2.  We assume the beam
//    travels +z with energy Eb(run).  “q” is virtual photon four‐vector: p_beam – p_e'.
//    Then p_pi = (E_pi, p_vec_pi).  Finally t = (ΔE)^2 – |Δ→p|^2.
//================================================================================
static double compute_t_scalar(int run,
                               double e_p, double e_theta, double e_phi,
                               double p_p, double p_theta, double p_phi)
{
    // 1) beam energy
    double Eb = beamEnergy(run);
    if (Eb <= 0.0) return 1e6; // invalid run → force fail

    // 2) scattered electron 4‐vector
    double E_e = std::sqrt(e_p*e_p + m_e*m_e);
    double sin_e = std::sin(e_theta);
    double cos_e = std::cos(e_theta);
    double ex = e_p * sin_e * std::cos(e_phi);
    double ey = e_p * sin_e * std::sin(e_phi);
    double ez = e_p * cos_e;

    // 3) pion 4‐vector
    double E_pi = std::sqrt(p_p*p_p + m_pi*m_pi);
    double sin_p = std::sin(p_theta);
    double cos_p = std::cos(p_theta);
    double px = p_p * sin_p * std::cos(p_phi);
    double py = p_p * sin_p * std::sin(p_phi);
    double pz = p_p * cos_p;

    // 4) virtual photon q = (Eb – E_e, –ex, –ey, Eb – ez)
    double E_q = Eb - E_e;
    double qx  = -ex;
    double qy  = -ey;
    double qz  = Eb - ez;

    // 5) Δ = q – p_pi
    double dE = E_q - E_pi;
    double dx = qx  - px;
    double dy = qy  - py;
    double dz = qz  - pz;

    // 6) t = (ΔE)^2 – (dx^2 + dy^2 + dz^2)
    return (dE*dE - (dx*dx + dy*dy + dz*dz));
}

//================================================================================
// applyCuts(…)
//    We add one new “property” name:  “enpi+”, which means “apply all of the
//    usual SingleHadronKinematicCuts plus |t|<1.0 calculated from e_p, e_theta,
//    e_phi, p_p, p_theta, p_phi, runnum.”
//================================================================================
bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC)
{
    // Basic naming lookup
    string property = binNames[currentFits];

    bool goodEvent = true;
    // 1) Standard DIS/Hadron cuts (common to almost everything):
    if (*Q2 <  1.0    ) return false;
    if (*W  <  2.0    ) return false;
    if (*y  >  0.75   ) return false;
    if (*fiducial_status != 2) return false;
    return true;
    // if (*p_p < 1.2    ) return false;
    // if (*xF  < 0.0    ) return false;
    // if (*Mx2 < 3.24   ) return false;

    // if (*runnum == 16234 || *runnum == 16235 || *runnum == 16236 || *runnum == 16243 ||
    //     *runnum == 16250 || *runnum == 16251 || *runnum == 16317 || *runnum == 16325 ||
    //     *runnum == 16658 || *runnum == 16659 || *runnum == 16660 || *runnum == 16664 || 
    //     *runnum == 16665 || *runnum == 16666 || *runnum == 16671 || *runnum == 16672 || 
    //     *runuum == 16673 || *runnum == 16674 || *runnum == 16675 || *runnum == 16676 ||
    //     *runnum == 16678 || *runnum == 16679 || *runnum == 16681 || *runnum == 16682 || 
    //     *runnum == 16683 || *runnum == 16685 || *runnum == 16686 || *runnum == 16687 || 
    //     *runnum == 16688 || *runnum == 16689 || *runnum == 16690 || *runnum == 16692 ||
    //     *runnum == 16693 || *runnum == 16695 || *runnum == 16721 || *runnum == 16722 || 
    //     *runnum == 16733 || *runnum == 16734 || *runnum == 16736) { return false; }


    // 2) If the property is “enpi,” impose |t| < 1.0 as well:
    if (property == "enpi") {
        // compute t from the branches
        int    rn     = *runnum;
        double ec_p   = *e_p;
        double ec_th  = *e_theta;
        double ec_ph  = *e_phi;
        double pi_p   = *p_p;
        double pi_th  = *p_theta;
        double pi_ph  = *p_phi;

        double t_val = compute_t_scalar(rn, ec_p, ec_th, ec_ph,
                                          pi_p, pi_th, pi_ph);

        // if (std::fabs(t_val) >= 1.0 || *Mx2 < 0.75 || *Mx2 > 1.050625) {
        if (std::fabs(t_val) < 0.07 || std::fabs(t_val) > 0.7 ||
            *y > 0.65 || *z < 0.55 || *Q2 > 8 || *Mx2 > 1.1) {
            return false;
        } else  {
            return true;
        }
    }

    if (property == "Fall18xB" || property == "Fall18pT" ||
        property == "Spring18xB" || property == "Spring18pT")
    {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 && *Mx2 > 2.25);
        return goodEvent;
    }
    if (property == "W" || property == "x") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 && *Mx2 > 1.8225);
        return goodEvent;
    }
    if (property == "integrated") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *Mx2 > 1.8225 && *y < 0.80);
        return goodEvent;
    }
    if (property == "Mx2") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80);
        return goodEvent;
    }
    if (property == "xF") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 && *Mx2 > 1.8225);
        return goodEvent;
    }
    if (property == "xFsmallPT") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && *pT < 0.5);
        return goodEvent;
    }
    if (property == "xFlargePT") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && *pT > 0.5);
        return goodEvent;
    }
    if (property == "xTFR"   || property == "xi"     || property == "PTTFR") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && *xF < 0.0);
        return goodEvent;
    }
    if (property == "xTFRsmallPT" || property == "xismallPT") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && *xF < 0.0 && *pT < 0.5);
        return goodEvent;
    }
    if (property == "xTFRlargePT" || property == "xilargePT") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && *xF < 0.0 && *pT > 0.5);
        return goodEvent;
    }
    if (property == "xCFR"   || property == "z"  || property == "PTCFR") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && *xF > 0.2);
        return goodEvent;
    }
    if (property == "x") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 && *Mx2 > 1.8225);
        return goodEvent;
    }

    // “xBsector”‐style cuts:
    if (property == "xBsector0") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 && *Mx2 > 1.8225);
        return goodEvent;
    }
    if (property == "xBsector1") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*p_phi < 0.2 || *p_phi > 5.5));
        return goodEvent;
    }
    if (property == "xBsector2") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*p_phi > 0.2 && *p_phi < 1.25));
        return goodEvent;
    }
    if (property == "xBsector3") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*p_phi > 1.25 && *p_phi < 2.25));
        return goodEvent;
    }
    if (property == "xBsector4") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*p_phi > 2.25 && *p_phi < 3.35));
        return goodEvent;
    }
    if (property == "xBsector5") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*p_phi > 3.35 && *p_phi < 4.4));
        return goodEvent;
    }
    if (property == "xBsector6") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*p_phi > 4.4 && *p_phi < 5.5));
        return goodEvent;
    }

    // “pTsector”‐style cuts:
    if (property == "pTsector0") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 && *Mx2 > 1.8225);
        return goodEvent;
    }
    if (property == "pTsector1") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*e_phi < 0.2 || *e_phi > 5.5));
        return goodEvent;
    }
    if (property == "pTsector2") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*e_phi > 0.2 && *e_phi < 1.25));
        return goodEvent;
    }
    if (property == "pTsector3") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*e_phi > 1.25 && *e_phi < 2.25));
        return goodEvent;
    }
    if (property == "pTsector4") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*e_phi > 2.25 && *e_phi < 3.35));
        return goodEvent;
    }
    if (property == "pTsector5") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*e_phi > 3.35 && *e_phi < 4.4));
        return goodEvent;
    }
    if (property == "pTsector6") {
        bool goodEvent = (*Q2 > 1.0 && *W > 2.0 && *y < 0.80 &&
                          *Mx2 > 1.8225 && (*e_phi > 4.4 && *e_phi < 5.5));
        return goodEvent;
    }

    // if (property == "pTsector0") { // meant to be all six sectors
    //   goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225;
    //   return goodEvent;
    // } 
    // if (property == "pTsector1") {
    //   goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*p_phi < 0.2 || *p_phi > 5.5);
    //   return goodEvent;
    // } 
    // if (property == "pTsector2") {
    //   goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*p_phi > 0.2 && *p_phi < 1.25);
    //   return goodEvent;
    // } 
    // if (property == "pTsector3") {
    //   goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*p_phi > 1.25 && *p_phi < 2.25);
    //   return goodEvent;
    // } 
    // if (property == "pTsector4") {
    //   goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*p_phi > 2.25 && *p_phi < 3.35);
    //   return goodEvent;
    // } 
    // if (property == "pTsector5") {
    //   goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*p_phi > 3.35 && *p_phi < 4.4);
    //   return goodEvent;
    // } 
    // if (property == "pTsector6") {
    //   goodEvent = *Q2 > 1 && *W > 2 && *y < 0.80 && *Mx2 > 1.8225 && (*p_phi > 4.4 && *p_phi < 5.5);
    //   return goodEvent;
    // } 

    
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


  // if (isMC || (*runnum < 16042 || *runnum > 17811)) {
  //   return goodEvent;
  // } else {
  //   // return goodEvent && *target_pol!=0;
  //   return goodEvent;
  // }
  return false;
}