#include "formatLabelName.h"

std::string formatLabelName(const std::string& original) {
    std::map<std::string, std::string> specialLabels = {
        {"Q2", "Q^{2} (GeV^{2})"},
        {"W", "W (GeV)"},
        {"Delta_eta", "#Delta#eta"},
        {"Delta_phi", "#Delta#phi"},
        {"Delta_phi12", "#Delta#phi_{12}"},
        {"Delta_phi13", "#Delta#phi_{13}"},
        {"Delta_phi23", "#Delta#phi_{23}"},
        {"eta1", "#eta_{1}"},
        {"eta2", "#eta_{2}"},
        {"eta3", "#eta_{3}"},
        {"eta12", "#eta_{12}"},
        {"eta13", "#eta_{13}"},
        {"eta23", "#eta_{23}"},
        {"pT", "P_{T} (GeV)"},
        {"pT1", "P_{1T} (GeV)"},
        {"pT2", "P_{2T} (GeV)"},
        {"pT23", "P_{3T} (GeV)"},
        {"pT12", "P_{12T} (GeV)"},
        {"pT13", "P_{13T} (GeV)"},
        {"pT23", "P_{23T} (GeV)"},
        {"pTpT", "P_{1T}P_{2T} (GeV^{2})"},
        {"Mh", "M_{h} (GeV)"},
        {"Mh12", "M_{h12} (GeV)"},
        {"Mh13", "M_{h13} (GeV)"},
        {"Mh23", "M_{h23} (GeV)"},
        {"t", "t (GeV^{2})"},
        {"phi1", "#phi_{1}"},
        {"phi2", "#phi_{2}"},
        {"phi3", "#phi_{3}"},
        {"phi12", "#phi_{12}"},
        {"phi13", "#phi_{13}"},
        {"phiR", "#phi_{R}"},
        {"tmin", "t_{min} (GeV^{2})"},
        {"e_p", "e_{p} (GeV)"},
        {"Mx", "M_{x} (GeV)"},
        {"Mx1", "M_{x1} (GeV)"},
        {"Mx2", "M_{x2} (GeV)"},
        {"Mx3", "M_{x3} (GeV)"},
        {"Mx12", "M_{x12} (GeV)"},
        {"Mx13", "M_{x13} (GeV)"},
        {"Mx23", "M_{x23} (GeV)"},
        // {"Mx2", "M_{x}^{2} (GeV)"},
        {"p_p", "p_{p} (GeV)"},
        {"p1_p", "p1_{p} (GeV)"},
        {"p2_p", "p2_{p} (GeV)"},
        {"t1", "t_{1}"},
        {"t2", "t_{2}"},
        {"t3", "t_{3}"},
        {"xF", "x_{F}"},
        {"xF1", "x_{F1}"},
        {"xF2", "x_{F2}"},
        {"x", "x_{B}"},
        {"z1", "z_{1}"},
        {"z2", "z_{2}"},
        {"z3", "z_{3}"},
        {"z12", "z_{12}"},
        {"z13", "z_{13}"},
        {"z23", "z_{23}"},
        {"zeta1", "#zeta_{1}"},
        {"zeta2", "#zeta_{2}"},
        {"zeta3", "#zeta_{3}"},
        {"zeta12", "#zeta_{12}"},
        {"zeta13", "#zeta_{13}"},
        {"zeta23", "#zeta_{23}"},
        {"vz_e", "v_{z_{e}} (cm)"},
        {"vz_p", "v_{z_{p}} (cm)"},
        {"vz_p1", "v_{z_{p1}} (cm)"},
        {"vz_p2", "v_{z_{p2}} (cm)"},
        {"vz_p3", "v_{z_{p3}} (cm)"}
    };
  
    if (specialLabels.find(original) != specialLabels.end()) {
        return specialLabels[original];
    }

    std::string formatted = original;
    size_t pos = 0;
    while ((pos = formatted.find('_', pos)) != std::string::npos) {
        formatted.replace(pos, 1, "_{");
        size_t closing = formatted.find('_', pos + 2);
        if (closing == std::string::npos) {
            closing = formatted.length();
        }
        formatted.insert(closing, "}");
        pos = closing + 1;
    }

    if (formatted.find("theta") != std::string::npos) {
        formatted.replace(formatted.find("theta"), 5, "#theta");
    }

    if (formatted.find("zeta") != std::string::npos) {
        formatted.replace(formatted.find("zeta"), 5, "#zeta");
    }

    if (formatted.find("phi") != std::string::npos) {
        formatted.replace(formatted.find("phi"), 3, "#phi");
    }

    if (formatted.find("eta") != std::string::npos && 
        formatted.find("theta") == std::string::npos && 
        formatted.find("zeta") == std::string::npos) {
        formatted.replace(formatted.find("eta"), 3, "#eta");
    }
  
    return formatted;
}