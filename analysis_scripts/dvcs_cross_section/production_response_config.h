#ifndef PRODUCTION_RESPONSE_CONFIG_H
#define PRODUCTION_RESPONSE_CONFIG_H
#include "global_cuts.h"
#include <nlohmann/json.hpp>
using json=nlohmann::json;
static GlobalCutConfig config(const json& j) {
    GlobalCutConfig c;
#define LOAD(name) c.name = j.at(#name).get<decltype(c.name)>()
    LOAD(t1_abs_max); LOAD(open_angle_min_deg);
    LOAD(enable_reconstructed_electron_p_cut); LOAD(reconstructed_electron_p_min_GeV);
    LOAD(enable_analysis_phase_space_cuts); LOAD(enable_analysis_y_min_cut);
    LOAD(analysis_y_min); LOAD(analysis_W_min_GeV); LOAD(enable_pTmiss_cut); LOAD(pTmiss_max);
    LOAD(enable_dvcsgen_ycol_cut); LOAD(dvcsgen_ycol_cut);
    LOAD(enable_auxiliary_fiducial_cuts); LOAD(auxiliary_require_distinct_fd_sectors);
    LOAD(auxiliary_e_theta_min_deg); LOAD(auxiliary_e_theta_max_deg);
    LOAD(auxiliary_fd_proton_theta_min_deg); LOAD(auxiliary_fd_proton_theta_max_deg);
    LOAD(auxiliary_fd_proton_p_min_inb_GeV); LOAD(auxiliary_fd_proton_p_min_out_GeV);
    LOAD(auxiliary_fd_photon_theta_min_deg); LOAD(auxiliary_fd_photon_theta_max_deg);
    LOAD(auxiliary_cd_proton_p_min_GeV); LOAD(auxiliary_cd_proton_theta_min_deg);
    LOAD(auxiliary_cd_proton_theta_max_deg); LOAD(auxiliary_ft_photon_p_min_GeV);
    LOAD(enable_topology_filter); LOAD(required_detector1); LOAD(required_detector2);
    LOAD(enable_electron_fd_sector_filter); LOAD(electron_fd_sector);
    LOAD(enable_proton_fd_sector_filter); LOAD(proton_fd_sector);
    LOAD(enable_proton_cd_sector_filter); LOAD(proton_cd_sector);
    LOAD(enable_photon_fd_sector_filter); LOAD(photon_fd_sector); LOAD(excluded_runs);
#undef LOAD
    return c;
}

#endif
