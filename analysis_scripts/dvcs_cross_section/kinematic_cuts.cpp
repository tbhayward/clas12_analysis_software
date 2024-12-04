#include "kinematic_cuts.h"
#include <string>
#include <iostream>

// Function that applies the kinematic cuts
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
) {
    // Applying the specified cuts

    // Universal cuts
    if (-t >= 1) return false;
    if (open_angle_ep2 <= 5) return false;

    // Variables for mu and sigma values
    double mu_theta_neutral_neutral = 0.0;
    double sigma_theta_neutral_neutral = 0.0;
    double mu_pTmiss = 0.0;
    double sigma_pTmiss = 0.0;
    double mu_xF = 0.0;
    double sigma_xF = 0.0;
    double mu_Emiss2 = 0.0;
    double sigma_Emiss2 = 0.0;
    double mu_Mx2 = 0.0;
    double sigma_Mx2 = 0.0;
    double mu_Mx2_1 = 0.0;
    double sigma_Mx2_1 = 0.0;
    double mu_Mx2_2 = 0.0;
    double sigma_Mx2_2 = 0.0;

    // Map run_period to internal period identifiers
    std::string period;
    if (run_period == "RGA Fa18 Inb" || run_period == "first period" || run_period == "Fa18Inb") {
        period = "period1";
    } else if (run_period == "RGA Fa18 Out" || run_period == "second period" || run_period == "Fa18Out") {
        period = "period2";
    } else if (run_period == "RGA Sp19 Inb" || run_period == "third period" || run_period == "Sp19Inb") {
        period = "period3";
    } else {
        std::cerr << "Unknown run period: " << run_period << std::endl;
        return false;
    }

    // Map topology to internal identifiers
    std::string topo;
    if (topology == "(FD,FD)") {
        topo = "FD_FD";
    } else if (topology == "(CD,FD)") {
        topo = "CD_FD";
    } else if (topology == "(CD,FT)") {
        topo = "CD_FT";
    } else {
        std::cerr << "Unknown topology: " << topology << std::endl;
        return false;
    }

    // Assign mu and sigma values based on the parameters
    if (channel == "dvcs") {
        if (data_type == "data") {
            if (period == "period1") { // First period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.294;
                    sigma_theta_neutral_neutral = 0.176;
                    mu_pTmiss = 0.057;
                    sigma_pTmiss = 0.035;
                    mu_xF = -0.036;
                    sigma_xF = 0.051;
                    mu_Emiss2 = 0.278;
                    sigma_Emiss2 = 0.309;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.006;
                    mu_Mx2_1 = 0.038;
                    sigma_Mx2_1 = 0.200;
                    mu_Mx2_2 = 1.288;
                    sigma_Mx2_2 = 0.458;
                } else if (topo == "CD_FT") {
                    mu_theta_neutral_neutral = 0.228;
                    sigma_theta_neutral_neutral = 0.153;
                    mu_pTmiss = 0.031;
                    sigma_pTmiss = 0.020;
                    mu_xF = -0.012;
                    sigma_xF = 0.034;
                    mu_Emiss2 = 0.150;
                    sigma_Emiss2 = 0.226;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.068;
                    sigma_Mx2_1 = 0.236;
                    mu_Mx2_2 = 1.119;
                    sigma_Mx2_2 = 0.362;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.362;
                    sigma_theta_neutral_neutral = 0.177;
                    mu_pTmiss = 0.070;
                    sigma_pTmiss = 0.039;
                    mu_xF = -0.034;
                    sigma_xF = 0.042;
                    mu_Emiss2 = 0.228;
                    sigma_Emiss2 = 0.253;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.005;
                    mu_Mx2_1 = 0.010;
                    sigma_Mx2_1 = 0.165;
                    mu_Mx2_2 = 1.145;
                    sigma_Mx2_2 = 0.301;
                }
            } else if (period == "period2") { // Second period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.323;
                    sigma_theta_neutral_neutral = 0.181;
                    mu_pTmiss = 0.057;
                    sigma_pTmiss = 0.035;
                    mu_xF = -0.036;
                    sigma_xF = 0.050;
                    mu_Emiss2 = 0.265;
                    sigma_Emiss2 = 0.289;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.006;
                    mu_Mx2_1 = 0.047;
                    sigma_Mx2_1 = 0.166;
                    mu_Mx2_2 = 1.293;
                    sigma_Mx2_2 = 0.456;
                } else if (topo == "CD_FT") {
                    mu_theta_neutral_neutral = 0.232;
                    sigma_theta_neutral_neutral = 0.156;
                    mu_pTmiss = 0.030;
                    sigma_pTmiss = 0.020;
                    mu_xF = -0.006;
                    sigma_xF = 0.036;
                    mu_Emiss2 = 0.141;
                    sigma_Emiss2 = 0.227;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.005;
                    mu_Mx2_1 = 0.126;
                    sigma_Mx2_1 = 0.251;
                    mu_Mx2_2 = 1.112;
                    sigma_Mx2_2 = 0.371;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.388;
                    sigma_theta_neutral_neutral = 0.179;
                    mu_pTmiss = 0.071;
                    sigma_pTmiss = 0.040;
                    mu_xF = -0.026;
                    sigma_xF = 0.043;
                    mu_Emiss2 = 0.184;
                    sigma_Emiss2 = 0.226;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.005;
                    mu_Mx2_1 = 0.035;
                    sigma_Mx2_1 = 0.162;
                    mu_Mx2_2 = 1.115;
                    sigma_Mx2_2 = 0.291;
                }
            } else if (period == "period3") { // Third period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.299;
                    sigma_theta_neutral_neutral = 0.175;
                    mu_pTmiss = 0.058;
                    sigma_pTmiss = 0.035;
                    mu_xF = -0.041;
                    sigma_xF = 0.048;
                    mu_Emiss2 = 0.282;
                    sigma_Emiss2 = 0.301;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.005;
                    mu_Mx2_1 = 0.007;
                    sigma_Mx2_1 = 0.159;
                    mu_Mx2_2 = 1.293;
                    sigma_Mx2_2 = 0.448;
                } else if (topo == "CD_FT") {
                    mu_theta_neutral_neutral = 0.245;
                    sigma_theta_neutral_neutral = 0.152;
                    mu_pTmiss = 0.032;
                    sigma_pTmiss = 0.020;
                    mu_xF = -0.015;
                    sigma_xF = 0.032;
                    mu_Emiss2 = 0.140;
                    sigma_Emiss2 = 0.223;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.027;
                    sigma_Mx2_1 = 0.194;
                    mu_Mx2_2 = 1.101;
                    sigma_Mx2_2 = 0.356;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.369;
                    sigma_theta_neutral_neutral = 0.177;
                    mu_pTmiss = 0.070;
                    sigma_pTmiss = 0.039;
                    mu_xF = -0.031;
                    sigma_xF = 0.043;
                    mu_Emiss2 = 0.214;
                    sigma_Emiss2 = 0.249;
                    mu_Mx2 = 0.000;
                    sigma_Mx2 = 0.005;
                    mu_Mx2_1 = 0.019;
                    sigma_Mx2_1 = 0.166;
                    mu_Mx2_2 = 1.130;
                    sigma_Mx2_2 = 0.296;
                }
            }
        } else if (data_type == "mc") {
            if (period == "period1") { // First period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.206;
                    sigma_theta_neutral_neutral = 0.148;
                    mu_pTmiss = 0.038;
                    sigma_pTmiss = 0.026;
                    mu_xF = -0.003;
                    sigma_xF = 0.043;
                    mu_Emiss2 = 0.029;
                    sigma_Emiss2 = 0.257;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.012;
                    sigma_Mx2_1 = 0.141;
                    mu_Mx2_2 = 0.933;
                    sigma_Mx2_2 = 0.363;
                } else if (topo == "CD_FT") {
                    mu_theta_neutral_neutral = 0.146;
                    sigma_theta_neutral_neutral = 0.120;
                    mu_pTmiss = 0.021;
                    sigma_pTmiss = 0.017;
                    mu_xF = -0.009;
                    sigma_xF = 0.028;
                    mu_Emiss2 = 0.106;
                    sigma_Emiss2 = 0.193;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.003;
                    mu_Mx2_1 = 0.036;
                    sigma_Mx2_1 = 0.172;
                    mu_Mx2_2 = 1.051;
                    sigma_Mx2_2 = 0.300;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.289;
                    sigma_theta_neutral_neutral = 0.164;
                    mu_pTmiss = 0.051;
                    sigma_pTmiss = 0.032;
                    mu_xF = -0.010;
                    sigma_xF = 0.039;
                    mu_Emiss2 = 0.015;
                    sigma_Emiss2 = 0.236;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = -0.042;
                    sigma_Mx2_1 = 0.146;
                    mu_Mx2_2 = 0.886;
                    sigma_Mx2_2 = 0.272;
                }
            } else if (period == "period2") { // Second period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.235;
                    sigma_theta_neutral_neutral = 0.158;
                    mu_pTmiss = 0.038;
                    sigma_pTmiss = 0.026;
                    mu_xF = -0.002;
                    sigma_xF = 0.044;
                    mu_Emiss2 = 0.023;
                    sigma_Emiss2 = 0.236;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.013;
                    sigma_Mx2_1 = 0.124;
                    mu_Mx2_2 = 0.924;
                    sigma_Mx2_2 = 0.354;
                } else if (topo == "CD_FT") {
                    mu_theta_neutral_neutral = 0.155;
                    sigma_theta_neutral_neutral = 0.125;
                    mu_pTmiss = 0.022;
                    sigma_pTmiss = 0.017;
                    mu_xF = -0.007;
                    sigma_xF = 0.032;
                    mu_Emiss2 = 0.111;
                    sigma_Emiss2 = 0.196;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.003;
                    mu_Mx2_1 = 0.066;
                    sigma_Mx2_1 = 0.204;
                    mu_Mx2_2 = 1.063;
                    sigma_Mx2_2 = 0.311;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.326;
                    sigma_theta_neutral_neutral = 0.172;
                    mu_pTmiss = 0.052;
                    sigma_pTmiss = 0.033;
                    mu_xF = -0.008;
                    sigma_xF = 0.039;
                    mu_Emiss2 = 0.012;
                    sigma_Emiss2 = 0.212;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = -0.020;
                    sigma_Mx2_1 = 0.135;
                    mu_Mx2_2 = 0.881;
                    sigma_Mx2_2 = 0.266;
                }
            } else if (period == "period3") { // Third period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.212;
                    sigma_theta_neutral_neutral = 0.149;
                    mu_pTmiss = 0.038;
                    sigma_pTmiss = 0.027;
                    mu_xF = -0.005;
                    sigma_xF = 0.044;
                    mu_Emiss2 = 0.037;
                    sigma_Emiss2 = 0.257;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.005;
                    sigma_Mx2_1 = 0.131;
                    mu_Mx2_2 = 0.943;
                    sigma_Mx2_2 = 0.366;
                } else if (topo == "CD_FT") {
                    mu_theta_neutral_neutral = 0.150;
                    sigma_theta_neutral_neutral = 0.123;
                    mu_pTmiss = 0.021;
                    sigma_pTmiss = 0.017;
                    mu_xF = -0.010;
                    sigma_xF = 0.028;
                    mu_Emiss2 = 0.104;
                    sigma_Emiss2 = 0.1990;
                    mu_Mx2 = 0.00;
                    sigma_Mx2 = 0.003;
                    mu_Mx2_1 = 0.028;
                    sigma_Mx2_1 = 0.153;
                    mu_Mx2_2 = 1.047;
                    sigma_Mx2_2 = 0.296;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.297;
                    sigma_theta_neutral_neutral = 0.166;
                    mu_pTmiss = 0.051;
                    sigma_pTmiss = 0.032;
                    mu_xF = -0.011;
                    sigma_xF = 0.040;
                    mu_Emiss2 = 0.019;
                    sigma_Emiss2 = 0.234;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = -0.040;
                    sigma_Mx2_1 = 0.141;
                    mu_Mx2_2 = 0.890;
                    sigma_Mx2_2 = 0.271;
                }
            }
        }
    } else if (channel == "eppi0") {
        // Check for unavailable topology
        if (topo == "CD_FT") {
            std::cerr << "Topology (CD,FT) not available for eppi0 channel." << std::endl;
            return false;
        }
        if (data_type == "data") {
            if (period == "period1") { // First period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.329;
                    sigma_theta_neutral_neutral = 0.174;
                    mu_pTmiss = 0.051;
                    sigma_pTmiss = 0.032;
                    mu_xF = -0.005;
                    sigma_xF = 0.049;
                    mu_Emiss2 = 0.042;
                    sigma_Emiss2 = 0.255;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.005;
                    mu_Mx2_1 = 0.032;
                    sigma_Mx2_1 = 0.126;
                    mu_Mx2_2 = 0.945;
                    sigma_Mx2_2 = 0.345;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.355;
                    sigma_theta_neutral_neutral = 0.173;
                    mu_pTmiss = 0.061;
                    sigma_pTmiss = 0.037;
                    mu_xF = -0.015;
                    sigma_xF = 0.039;
                    mu_Emiss2 = 0.062;
                    sigma_Emiss2 = 0.210;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.002;
                    sigma_Mx2_1 = 0.101;
                    mu_Mx2_2 = 0.947;
                    sigma_Mx2_2 = 0.247;
                }
            } else if (period == "period2") { // Second period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.348;
                    sigma_theta_neutral_neutral = 0.174;
                    mu_pTmiss = 0.050;
                    sigma_pTmiss = 0.032;
                    mu_xF = -0.002;
                    sigma_xF = 0.047;
                    mu_Emiss2 = 0.034;
                    sigma_Emiss2 = 0.232;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.045;
                    sigma_Mx2_1 = 0.111;
                    mu_Mx2_2 = 0.937;
                    sigma_Mx2_2 = 0.346;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.372;
                    sigma_theta_neutral_neutral = 0.172;
                    mu_pTmiss = 0.060;
                    sigma_pTmiss = 0.036;
                    mu_xF = -0.009;
                    sigma_xF = 0.038;
                    mu_Emiss2 = 0.046;
                    sigma_Emiss2 = 0.201;
                    mu_Mx2 = 0.000;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.026;
                    sigma_Mx2_1 = 0.101;
                    mu_Mx2_2 = 0.932;
                    sigma_Mx2_2 = 0.261;
                }
            } else if (period == "period3") { // Third period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.334;
                    sigma_theta_neutral_neutral = 0.172;
                    mu_pTmiss = 0.051;
                    sigma_pTmiss = 0.032;
                    mu_xF = -0.008;
                    sigma_xF = 0.047;
                    mu_Emiss2 = 0.044;
                    sigma_Emiss2 = 0.252;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.016;
                    sigma_Mx2_1 = 0.104;
                    mu_Mx2_2 = 0.946;
                    sigma_Mx2_2 = 0.345;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.359;
                    sigma_theta_neutral_neutral = 0.172;
                    mu_pTmiss = 0.060;
                    sigma_pTmiss = 0.036;
                    mu_xF = -0.012;
                    sigma_xF = 0.040;
                    mu_Emiss2 = 0.048;
                    sigma_Emiss2 = 0.210;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.007;
                    sigma_Mx2_1 = 0.100;
                    mu_Mx2_2 = 0.933;
                    sigma_Mx2_2 = 0.247;
                }
            }
        } else if (data_type == "mc") {
            if (period == "period1") { // First period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.261;
                    sigma_theta_neutral_neutral = 0.159;
                    mu_pTmiss = 0.046;
                    sigma_pTmiss = 0.031;
                    mu_xF = 0.002;
                    sigma_xF = 0.045;
                    mu_Emiss2 = -0.008;
                    sigma_Emiss2 = 0.244;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.026;
                    sigma_Mx2_1 = 0.118;
                    mu_Mx2_2 = 0.879;
                    sigma_Mx2_2 = 0.333;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.310;
                    sigma_theta_neutral_neutral = 0.168;
                    mu_pTmiss = 0.057;
                    sigma_pTmiss = 0.036;
                    mu_xF = -0.002;
                    sigma_xF = 0.038;
                    mu_Emiss2 = -0.012;
                    sigma_Emiss2 = 0.210;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.004;
                    sigma_Mx2_1 = 0.106;
                    mu_Mx2_2 = 0.852;
                    sigma_Mx2_2 = 0.245;
                }
            } else if (period == "period2") { // Second period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.283;
                    sigma_theta_neutral_neutral = 0.165;
                    mu_pTmiss = 0.046;
                    sigma_pTmiss = 0.031;
                    mu_xF = 0.005;
                    sigma_xF = 0.045;
                    mu_Emiss2 = -0.020;
                    sigma_Emiss2 = 0.222;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.003;
                    mu_Mx2_1 = 0.025;
                    sigma_Mx2_1 = 0.106;
                    mu_Mx2_2 = 0.857;
                    sigma_Mx2_2 = 0.329;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.333;
                    sigma_theta_neutral_neutral = 0.170;
                    mu_pTmiss = 0.058;
                    sigma_pTmiss = 0.036;
                    mu_xF = 0.000;
                    sigma_xF = 0.038;
                    mu_Emiss2 = -0.021;
                    sigma_Emiss2 = 0.200;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.013;
                    sigma_Mx2_1 = 0.108;
                    mu_Mx2_2 = 0.841;
                    sigma_Mx2_2 = 0.252;
                }
            } else if (period == "period3") { // Third period
                if (topo == "CD_FD") {
                    mu_theta_neutral_neutral = 0.262;
                    sigma_theta_neutral_neutral = 0.160;
                    mu_pTmiss = 0.046;
                    sigma_pTmiss = 0.031;
                    mu_xF = 0.001;
                    sigma_xF = 0.045;
                    mu_Emiss2 = -0.002;
                    sigma_Emiss2 = 0.242;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.022;
                    sigma_Mx2_1 = 0.110;
                    mu_Mx2_2 = 0.885;
                    sigma_Mx2_2 = 0.334;
                } else if (topo == "FD_FD") {
                    mu_theta_neutral_neutral = 0.319;
                    sigma_theta_neutral_neutral = 0.168;
                    mu_pTmiss = 0.057;
                    sigma_pTmiss = 0.035;
                    mu_xF = -0.002;
                    sigma_xF = 0.038;
                    mu_Emiss2 = -0.014;
                    sigma_Emiss2 = 0.207;
                    mu_Mx2 = -0.001;
                    sigma_Mx2 = 0.004;
                    mu_Mx2_1 = 0.004;
                    sigma_Mx2_1 = 0.104;
                    mu_Mx2_2 = 0.853;
                    sigma_Mx2_2 = 0.243;
                }
            }
        }
    } else {
        std::cerr << "Unknown channel: " << channel << std::endl;
        return false;
    }

    // Apply the +/- 3 sigma cuts
    if (!(theta_neutral_neutral >= (mu_theta_neutral_neutral - 3 * sigma_theta_neutral_neutral) &&
          theta_neutral_neutral <= (mu_theta_neutral_neutral + 3 * sigma_theta_neutral_neutral)))
        return false;

    if (!(pTmiss >= (mu_pTmiss - 3 * sigma_pTmiss) &&
          pTmiss <= (mu_pTmiss + 3 * sigma_pTmiss)))
        return false;

    if (!(xF >= (mu_xF - 3 * sigma_xF) &&
          xF <= (mu_xF + 3 * sigma_xF)))
        return false;

    if (!(Emiss2 >= (mu_Emiss2 - 3 * sigma_Emiss2) &&
          Emiss2 <= (mu_Emiss2 + 3 * sigma_Emiss2)))
        return false;

    if (!(Mx2 >= (mu_Mx2 - 3 * sigma_Mx2) &&
          Mx2 <= (mu_Mx2 + 3 * sigma_Mx2)))
        return false;

    if (!(Mx2_1 >= (mu_Mx2_1 - 3 * sigma_Mx2_1) &&
          Mx2_1 <= (mu_Mx2_1 + 3 * sigma_Mx2_1)))
        return false;

    if (!(Mx2_2 >= (mu_Mx2_2 - 3 * sigma_Mx2_2) &&
          Mx2_2 <= (mu_Mx2_2 + 3 * sigma_Mx2_2)))
        return false;

    return true;  // Event passes all cuts
}