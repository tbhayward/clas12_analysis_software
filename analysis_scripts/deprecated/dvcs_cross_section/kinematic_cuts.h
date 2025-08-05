#ifndef KINEMATIC_CUTS_H
#define KINEMATIC_CUTS_H

#include <string>

// Declaration of the kinematic cuts function
bool apply_kinematic_cuts(
    double t,
    double open_angle_ep2,
    double theta_neutral_neutral,
    double Emiss2,
    double Mx2,
    double Mx2_1,
    double Mx2_2,
    double pTmiss,
    double xF,
    const std::string& channel,      // "dvcs" or "eppi0"
    const std::string& data_type,    // "data" or "mc"
    const std::string& run_period,   // e.g., "RGA Fa18 Inb", "RGA Fa18 Out", "RGA Sp19 Inb"
    const std::string& topology      // e.g., "(FD,FD)", "(CD,FD)", "(CD,FT)"
);

#endif // KINEMATIC_CUTS_H