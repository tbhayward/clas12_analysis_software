#include "histConfigs.h"
#include "TMath.h"

// Define the histogram binning configurations
std::map<std::string, HistConfig> histConfigs = {
    {"open_angle_ep2", {50, 0, 60}},   // Adjusted range
    {"Mx2_2", {50, -1, 4}},
    {"theta_gamma_gamma", {50, 0, 7}}, 
    {"placeholder", {50, 0, 0}}, // Placeholder for future use
    {"Emiss2", {50, -1, 2}},
    {"Mx2", {50, -0.06, 0.06}},
    {"Mx2_1", {50, -1, 1.5}},
    {"pTmiss", {50, 0, 0.7}}
};