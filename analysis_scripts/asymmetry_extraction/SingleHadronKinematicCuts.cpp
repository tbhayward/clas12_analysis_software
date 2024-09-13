#include "SingleHadronKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>
#include <optional>
#include "TMath.h"

using std::string;

SingleHadronKinematicCuts::SingleHadronKinematicCuts(TTreeReader& reader, TTree* tree)
    : BaseKinematicCuts(reader), // Call to the BaseKinematicCuts constructor
      runnum(reader, "runnum"),
      e_theta(reader, "e_theta"), e_phi(reader, "e_phi"), vz_e(reader, "vz_e"),
      vz_p(reader, "vz_p"),
      Q2(reader, "Q2"), W(reader, "W"), Mx(reader, "Mx"), x(reader, "x"),
      y(reader, "y"), z(reader, "z"), pT(reader, "pT"), xF(reader, "xF"),
      phi(reader, "phi"), 
      target_pol(reader, "target_pol")
{
    // Initialize Mx1, Mx2, and Mx23 only if the corresponding branches exist
    if (tree->GetBranch("Mx1")) {
        Mx1.emplace(reader, "Mx1");
    }
    if (tree->GetBranch("Mx2")) {
        Mx2.emplace(reader, "Mx2");
    }
    if (tree->GetBranch("Mx23")) {
        Mx23.emplace(reader, "Mx23");
    }
}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];  // Get the property dynamically

    // Define cuts based on property
    if (property == "epiplus") {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx > 1.5 && *y < 0.80;
        return goodEvent;
    }
    if (property == "epipluspiminus" && Mx1) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx1 > 1.5 && *y < 0.80;
        return goodEvent;
    }
    if (property == "epipluspiminus_rho0_free" && Mx1) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx1 > 1.5 && *Mx > 1.05 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppiplus" && Mx2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2 > 1.5 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppiplus_rho0_free" && Mx1 && Mx2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2 > 1.5 && **Mx1 > 0.95 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus" && Mx2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2 > 1.5 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus_rho0_free_A" && Mx2 && Mx23) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2 > 1.5 && **Mx23 > 1.05 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus_rho0_free_B" && Mx1 && Mx2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2 > 1.5 && **Mx1 > 0.95 && *y < 0.80;
        return goodEvent;
    }

    // Final event selection
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
        return goodEvent;
    } else {
        return goodEvent;
    }
}