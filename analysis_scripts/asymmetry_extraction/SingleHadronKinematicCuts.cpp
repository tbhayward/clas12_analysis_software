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
    return 10.604;
}

//================================================================================
// compute_t(…) 
//    Given arrays of runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi (all scalars
//    for one event), return t = (q – p_pi)^2 = (ΔE)^2 – (Δp)^2.  We assume the beam
//    travels +z with energy Eb(run).  “q” is virtual photon four‐vector: p_beam – p_e'.
//    Then p_pi = (E_pi, p_vec_pi).  
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

    if (property == "enpi") {
        goodEvent = goodEvent && *Mx2 > 0.0 && *Mx2 < 1.07 && 
            *x > 0.1 && *x < 0.6;
        return goodEvent;
    }
    if (property == "enpiLowxB") {
        goodEvent = goodEvent && *Mx2 > 0.0 && *Mx2 < 1.07 && 
            *x > 0.1 && *x < 0.25;
        return goodEvent;
    }
    if (property == "enpiMidLowxB") {
        goodEvent = goodEvent && *Mx2 > 0.0 && *Mx2 < 1.07 && 
            *x > 0.25 && *x < 0.35;
        return goodEvent;
    }
    if (property == "enpiMidHighxB") {
        goodEvent = goodEvent && *Mx2 > 0.0 && *Mx2 < 1.07 && 
            *x > 0.35 && *x < 0.45;
        return goodEvent;
    }
    if (property == "enpiHighxB") {
        goodEvent = goodEvent && *Mx2 > 0.0 && *Mx2 < 1.07 && 
            *x > 0.45 && *x < 0.60;
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