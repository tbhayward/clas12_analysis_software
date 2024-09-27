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
      Q2(reader, "Q2"), W(reader, "W"),  x(reader, "x"),
      y(reader, "y"), z(reader, "z"), pT(reader, "pT"), xF(reader, "xF"),
      target_pol(reader, "target_pol")
{
    // Initialize Mx1, Mx2, and Mx23 only if the corresponding branches exist
    if (tree->GetBranch("Mx2_1")) {
        Mx2_1.emplace(reader, "Mx2_1");
    }
    if (tree->GetBranch("Mx2_")) {
        Mx2_2.emplace(reader, "Mx2_");
    }
    if (tree->GetBranch("Mx2_23")) {
        Mx2_23.emplace(reader, "Mx2_23");
    }
}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];  // Get the property dynamically

    // Define cuts based on property
    if (property == "epiplus") {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx2 > 3.24 && *y < 0.80;
        return goodEvent;
    }
    if (property == "epipluspiminus" && Mx2_1) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2_1 > 2.25 && *y < 0.80;
        return goodEvent;
    }
    if (property == "epipluspiminus_rho0_free" && Mx2_1) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2_1 > 2.25 && *Mx2_2 > 2.25 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppiplus" && Mx2_2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2_1 > 3.24 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppiplusNoRho" && Mx2_2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2_1 > 3.24 && **Mx2_2 > 1.35 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppiplus_rho0_free" && Mx2_1 && Mx2_2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2_2 > 2.25 && **Mx2_1 > 1.35 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus" && Mx2_2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2_2 > 2.25 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus_rho0_free_A" && Mx2_2 && Mx2_23) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2_2 > 2.25 && **Mx2_23 > 2.25 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus_rho0_free_B" && Mx2_1 && Mx2_2) {
        goodEvent = *Q2 > 1 && *W > 2 && **Mx2_2 > 2.25 && **Mx2_1 > 1.35 && *y < 0.80;
        return goodEvent;
    }

    // Final event selection
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
        return goodEvent;
    } else {
        return goodEvent;
    }
}