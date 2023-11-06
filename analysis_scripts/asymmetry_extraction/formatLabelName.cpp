#include "formatLabelName.h"

std::string formatLabelName(const std::string& original) {
    std::map<std::string, std::string> specialLabels = {
        {"Q2", "Q^{2} (GeV^{2})"},
        {"W", "W (GeV)"},
        {"Delta_eta", "#Delta #eta"},
        {"Delta_phi", "#Delta #phi"},
        {"eta1", "#eta_{1}"},
        {"eta2", "#eta_{2}"},
        {"pT", "P_{T} (GeV)"},
        {"pT1", "P_{1T} (GeV)"},
        {"pT2", "P_{2T} (GeV)"},
        {"pTpT", "P_{1T}P_{2T} (GeV^{2})"},
        {"Mh", "M_{h} (GeV)"},
        {"t", "t (GeV^{2})"},
        {"tmin", "t_{min} (GeV^{2})"},
        {"e_p", "e_{p} (GeV)"},
        {"Mx", "M_{x} (GeV)"},
        {"Mx2", "M_{x}^{2} (GeV)"},
        {"p_p", "p_{p} (GeV)"},
        {"xF", "x_{F}"},
        {"xF1", "x_{F1}"},
        {"xF2", "x_{F2}"},
        {"x", "x_{B}"},
        {"z1", "z_{1}"},
        {"z2", "z_{2}"},
        {"zeta1", "zeta_{1}"},
        {"zeta2", "zeta_{2}"},
        {"vz_e", "{v_{ez}"},
        {"vz_p", "{v_{pz}"},
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