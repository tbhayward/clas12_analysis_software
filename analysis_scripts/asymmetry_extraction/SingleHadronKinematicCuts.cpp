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
      phi(reader, "phi"), phi(reader, "phi2"), 
      target_pol(reader, "target_pol") {}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
  bool goodEvent = false;
  bool checked = false;
  string property = binNames[currentFits];

  if (-10 > *vz_p || *vz_p > 1.5 || -9 > *vz_e || *vz_e > 2) return false;

  if (property == "integrated") {
    goodEvent = *Q2 > 1 && *W > 2 && *Mx > 1.35 && *y < 0.80;
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