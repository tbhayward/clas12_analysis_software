#include "formatLabelName.h"
#include <map>

std::string formatLabelName(const std::string& original, const std::string& analysisType) {
    // Define common labels for both DVCS and eppi0
    std::map<std::string, std::string> commonLabels = {
        {"open_angle_ep2", "#theta_{e'#gamma} (deg)"},
        {"Mx2_2", "M M_{e'#gamma}^{2} (GeV^{2})"},
        {"theta_gamma_gamma", "#theta_{#gamma(det)#gamma(rec)}"},
        {"theta_pi0_pi0", "#theta_{#pi^{0}(det)#pi^{0}(rec)}"},
        {"xF", "x_{F}"},
        {"Emiss2", "M E_{e'p'#gamma} (GeV)"},
        {"Mx2", "M M_{e'p'#gamma}^{2} (GeV^{2})"},
        {"Mx2_1", "M M_{e'p'}^{2} (GeV^{2})"},
        {"pTmiss", "M P_{T(e'p'#gamma)} (GeV)"}
    };

    // Determine if we are dealing with DVCS or eppi0
    if (commonLabels.find(original) != commonLabels.end()) {
        std::string label = commonLabels[original];
        // If analysisType is "eppi0", replace #gamma with #pi^{0}
        if (analysisType == "eppi0") {
            size_t pos;
            while ((pos = label.find("#gamma")) != std::string::npos) {
                label.replace(pos, 6, "#pi^{0}");
            }
        }
        return label;
    }

    return original;  // Default return if no special formatting is defined
}