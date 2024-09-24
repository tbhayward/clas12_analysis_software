#include "histConfigs.h"
#include "TMath.h"

// Define the histogram binning configurations
std::map<std::string, HistConfig> histConfigs = {
    {"open_angle_ep2", {100, 0, 60}},   // Adjusted range
    {"Mx2_2", {100, 0, 4}},
    {"theta_gamma_gamma", {100, 0, 3}}, 
    {"theta_pi0_pi0", {100, 0, 5}}, 
    {"xF", {100, -0.4, 0.2}}, // Placeholder for future use
    {"Emiss2", {100, -1, 2}},
    {"Mx2", {100, -0.03, 0.03}},
    {"Mx2_1", {100, -1, 1.5}},
    {"pTmiss", {100, 0, 0.5}}
};