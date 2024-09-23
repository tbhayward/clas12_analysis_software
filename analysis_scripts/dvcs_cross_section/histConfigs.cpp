#include "histConfigs.h"
#include "TMath.h"

// Define the histogram binning configurations
std::map<std::string, HistConfig> histConfigs = {
    {"open_angle_ep2", {100, 0, 60}},   // Adjusted range
    {"Mx2_2", {100, -1, 4}},
    {"theta_gamma_gamma", {100, 0, 7}}, 
    {"placeholder", {100, 0, 0}}, // Placeholder for future use
    {"Emiss2", {100, -1, 2}},
    {"Mx2", {100, -0.06, 0.06}},
    {"Mx2_1", {100, -1, 1.5}},
    {"pTmiss", {100, 0, 0.7}}
};