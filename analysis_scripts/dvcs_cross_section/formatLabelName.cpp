#include "formatLabelName.h"
#include <map>

std::string formatLabelName(const std::string& original) {
    std::map<std::string, std::string> specialLabels = {
        {"open_angle_ep2", "#theta_{e'#gamma} (deg)"},
        {"Mx2_2", "M M_{e'#gamma}^{2} (GeV^{2})"},
        {"theta_gamma_gamma", "#theta_{#gamma(det)#gamma(rec)}"},
        {"theta_pi0_pi0", "#theta_{#pi^{0}(det)#pi^{0}(rec)}"},
        {"placeholder", ""},
        {"Emiss2", "M E_{e'p'#gamma} (GeV)"},
        {"Mx2", "M M_{e'p'#gamma}^{2} (GeV^{2})"},
        {"Mx2_1", "M M_{e'p'}^{2} (GeV^{2})"},
        {"pTmiss", "M P_{T(e'p'#gamma)} (GeV)"}
    };

    if (specialLabels.find(original) != specialLabels.end()) {
        return specialLabels[original];
    }

    return original;  // Default return if no special formatting is defined
}