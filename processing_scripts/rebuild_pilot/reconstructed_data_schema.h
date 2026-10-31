#ifndef RECONSTRUCTED_DATA_SCHEMA_PRIORITY1_20261031_H
#define RECONSTRUCTED_DATA_SCHEMA_PRIORITY1_20261031_H

#include <string>
#include <unordered_set>
#include <vector>

namespace reconstructed_data_schema {

inline const std::vector<std::string>& dvcs() {
    static const std::vector<std::string> names = {
        "fiducial_status", "num_pos", "num_neg", "num_neutrals", "runnum", "evnum",
        "helicity", "detector1", "detector2", "e_p", "e_theta", "e_phi", "vz_e",
        "p1_p", "p1_theta", "p1_phi", "vz_p1", "p2_p", "p2_theta", "p2_phi",
        "vz_p2", "open_angle_ep", "open_angle_ep1", "open_angle_ep2", "open_angle_p1p2",
        "Q2", "W", "Mx2", "Mx2_1", "Mx2_2", "x", "t", "t1", "t2", "tmin", "y",
        "z", "z1", "z2", "Mh", "xF", "xF1", "xF2", "pT", "pT1", "pT2", "pTpT",
        "xi", "xi1", "xi2", "eta", "eta1", "eta2", "Delta_eta", "eta1_gN", "eta2_gN",
        "phi1", "phi2", "Delta_phi", "phih", "phiR", "theta", "DepA", "DepB", "DepC",
        "DepV", "DepW", "Emiss2", "theta_gamma_gamma", "pTmiss"
    };
    return names;
}

inline const std::vector<std::string>& pi0() {
    static const std::vector<std::string> names = {
        "fiducial_status", "num_pos", "num_neg", "num_neutrals", "runnum", "evnum",
        "helicity", "detector1", "detector2", "e_p", "e_theta", "e_phi", "vz_e",
        "p1_p", "p1_theta", "p1_phi", "vz_p1", "p2_p", "p2_theta", "p2_phi",
        "vz_p2", "open_angle_ep", "open_angle_ep1", "open_angle_ep2", "open_angle_p1p2",
        "Q2", "W", "Mx2", "Mx2_1", "Mx2_2", "x", "t", "t1", "t2", "tmin", "y",
        "z", "z1", "z2", "Mh", "xF", "xF1", "xF2", "pT", "pT1", "pT2", "pTpT",
        "xi", "xi1", "xi2", "eta", "eta1", "eta2", "Delta_eta", "eta1_gN", "eta2_gN",
        "phi1", "phi2", "Delta_phi", "phih", "phiR", "theta", "DepA", "DepB", "DepC",
        "DepV", "DepW", "Mh_gammagamma", "detector_gamma1", "detector_gamma2",
        "open_angle_egamma1", "open_angle_egamma2", "gamma_phi1", "gamma_phi2", "Emiss2",
        "theta_pi0_pi0", "pTmiss"
    };
    return names;
}

inline const std::vector<std::string>* for_mode(const std::string& mode) {
    if (mode == "dvcs_data" || mode == "pi0_as_dvcs") return &dvcs();
    if (mode == "pi0_data") return &pi0();
    return nullptr;
}

inline bool definition_is_valid(const std::string& mode,
                                const std::vector<std::string>& names) {
    const std::size_t required = mode == "pi0_data" ? 77 : 70;
    if (names.size() != required) return false;
    std::unordered_set<std::string> unique;
    for (const auto& name : names) {
        if (name.empty() || !unique.insert(name).second) return false;
    }
    return unique.size() == required;
}

}  // namespace reconstructed_data_schema

#endif
