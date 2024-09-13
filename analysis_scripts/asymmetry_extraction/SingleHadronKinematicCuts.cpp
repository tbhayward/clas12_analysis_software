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
      phi(reader, "phi"), phi2(reader, "phi2"),
      target_pol(reader, "target_pol") {}

SingleHadronKinematicCuts::~SingleHadronKinematicCuts() {
    // Clean up dynamically allocated memory
    if (Mx1) delete Mx1;
    if (Mx2) delete Mx2;
    if (Mx23) delete Mx23;
}

bool SingleHadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];  // Get the property dynamically

    // Conditionally initialize Mx1, Mx2, Mx23 based on the property
    if ((property == "epipluspiminus" || property == "epipluspiminus_rho0_free") && !Mx1) {
        Mx1 = new TTreeReaderValue<double>(reader, "Mx1");
    }
    if ((property == "eppiplus" || property == "eppiplus_rho0_free" || property == "eppipluspiminus") && !Mx2) {
        Mx2 = new TTreeReaderValue<double>(reader, "Mx2");
    }
    if (property == "eppipluspiminus_rho0_free_A" && !Mx23) {
        Mx23 = new TTreeReaderValue<double>(reader, "Mx23");
    }

    // Define cuts based on property
    if (property == "epiplus") {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx > 1.5 && *y < 0.80;
        return goodEvent;
    }
    if (property == "epipluspiminus" && Mx1) {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx1 > 1.5 && *y < 0.80;
        return goodEvent;
    }
    if (property == "epipluspiminus_rho0_free" && Mx1) {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx1 > 1.5 && *Mx > 1.05 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppiplus" && Mx2) {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx2 > 1.5 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppiplus_rho0_free" && Mx1 && Mx2) {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx2 > 1.5 && *Mx1 > 0.95 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus" && Mx2) {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx2 > 1.5 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus_rho0_free_A" && Mx2 && Mx23) {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx2 > 1.5 && *Mx23 > 1.05 && *y < 0.80;
        return goodEvent;
    }
    if (property == "eppipluspiminus_rho0_free_B" && Mx1 && Mx2) {
        goodEvent = *Q2 > 1 && *W > 2 && *Mx2 > 1.5 && *Mx1 > 0.95 && *y < 0.80;
        return goodEvent;
    }

    // Final event selection
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
        return goodEvent;
    } else {
        return goodEvent;
    }
}