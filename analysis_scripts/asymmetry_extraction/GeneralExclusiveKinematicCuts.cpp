#include "GeneralExclusiveKinematicCuts.h"
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
GeneralExclusiveKinematicCuts::GeneralExclusiveKinematicCuts(TTreeReader& reader)
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
// applyCuts(…)
//    We add one new “property” name:  “enpi+”, which means “apply all of the
//    usual GeneralExclusiveKinematicCuts plus |t|<1.0 calculated from e_p, e_theta,
//    e_phi, p_p, p_theta, p_phi, runnum.”
//================================================================================
bool GeneralExclusiveKinematicCuts::applyCuts(int currentFits, bool isMC)
{
    // Basic naming lookup
    string property = binNames[currentFits];

    bool goodEvent = true;
    // 1) Standard DIS/Hadron cuts (common to almost everything):
    if (*Q2 <  1.0    ) return false;
    if (*W  <  2.0    ) return false;
    if (*y  >  0.75   ) return false;

    // 2) If the property is “enpi,” impose |t| < 1.0 as well:
    if (property == "enpi") {
        goodEvent = goodEvent && *fiducial_status >= 100 && 
            std::fabs(*t) <= 1.0 && std::fabs(*t) >= 0.0 &&
            *Mx2 > 0.80 && *Mx2 < 1.00;
        return goodEvent;
    }
    if (property == "enpiLowt") {
        goodEvent = goodEvent && std::fabs(*t) >= 0.00 && std::fabs(*t) <= 0.30;
        goodEvent = goodEvent && *fiducial_status >= 100 && 
            std::fabs(*t) <= 1.0 && std::fabs(*t) >= 0.0 &&
            *Mx2 > 0.80 && *Mx2 < 1.00;
        return goodEvent;
    }
    if (property == "enpiMidt") {
        goodEvent = goodEvent && std::fabs(*t) >= 0.30 && std::fabs(*t) <= 0.70;
        goodEvent = goodEvent && *fiducial_status >= 100 && 
            std::fabs(*t) <= 1.0 && std::fabs(*t) >= 0.0 &&
            *Mx2 > 0.80 && *Mx2 < 1.00;
        return goodEvent;
    }
    if (property == "enpiHight") {
        goodEvent = goodEvent && std::fabs(*t) >= 0.70;
        goodEvent = goodEvent && *fiducial_status >= 100 && 
            std::fabs(*t) <= 1.0 && std::fabs(*t) >= 0.0 &&
            *Mx2 > 0.80 && *Mx2 < 1.00;
        return goodEvent;
    }
    if (property == "enpiHarutsBin") {
        goodEvent = goodEvent && std::fabs(*t) >= 0.47 && std::fabs(*t) <= 0.87;
        goodEvent = goodEvent && *fiducial_status >= 100 && 
            std::fabs(*t) <= 1.0 && std::fabs(*t) >= 0.0 &&
            *Mx2 > 0.80 && *Mx2 < 1.00;
        return goodEvent;
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

  // if (isMC || (*runnum < 16042 || *runnum > 17811)) {
  //   return goodEvent;
  // } else {
  //   // return goodEvent && *target_pol!=0;
  //   return goodEvent;
  // }
  return false;
}