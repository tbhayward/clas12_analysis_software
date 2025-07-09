#include "dilution_factor.h"
#include <cmath>
#include <iostream>
#include <TRandom3.h>

// Function to determine the appropriate 4D bin prefix
std::string determine_4D_bin(double Q2, double x, double z) {
    std::string prez;  // Stores the Q2-x bin prefix
    std::string postz; // Stores the z-bin suffix

    // Determine the Q2-x bin
    if (x < 0.1) {
        prez = "Q2x1";
    } else if (x > 0.1 && x < 0.14) {
        if (Q2 < 1.50) prez = "Q2x2";
        else if (Q2 >= 1.50 && Q2 < 1.70) prez = "Q2x3";
        else if (Q2 >= 1.70) prez = "Q2x4";
    } else if (x > 0.14 && x < 0.21) {
        if (Q2 < 1.50) prez = "Q2x5";
        else if (Q2 >= 1.50 && Q2 < 1.70) prez = "Q2x6";
        else if (Q2 >= 1.70 && Q2 < 2.00) prez = "Q2x7";
        else if (Q2 >= 2.00) prez = "Q2x8";
    } else if (x > 0.21 && x < 0.30) {
        if (Q2 < 2.20) prez = "Q2x9";
        else if (Q2 >= 2.20 && Q2 < 2.60) prez = "Q2x10";
        else if (Q2 >= 2.60) prez = "Q2x11";
    } else if (x > 0.30 && x < 0.42) {
        if (Q2 < 3.20) prez = "Q2x12";
        else if (Q2 >= 3.20) prez = "Q2x13";
    } else if (x >= 0.42) {
        prez = "Q2x14";
    }

    // Determine the z-bin
    if (z > 0 && z <= 0.19) {
        postz = "z1";
    } else if (z > 0.19 && z <= 0.30) {
        postz = "z2";
    } else if (z > 0.30 && z <= 0.42) {
        postz = "z3";
    } else if (z > 0.42 && z <= 1.00) {
        postz = "z4";
    }

    // Combine to get the full prefix for the 4D bin
    std::string full_prefix = prez + postz;
    return full_prefix;  // Return the complete 4D bin prefix
}


// Main function to calculate the dilution factor
double dilution_factor(double Q2, double x, double z, double pT, const std::string& prefix) {
    TRandom3 rand_gen;

    // if (prefix == "x") { return 0.136288+0.283353*x+-0.125745*std::pow(x,2); }
    // if (prefix == "xF") { return 0.136288+0.283353*x+-0.125745*std::pow(x,2); }
    // if (prefix == "PT") { return 0.136288+0.283353*x+-0.125745*std::pow(x,2); }

    // if (prefix == "xall") { return 0.163228+-0.0367291*x+0.804498*std::pow(x,2); }
    // if (prefix == "xFall") { return 0.163228+-0.0367291*x+0.804498*std::pow(x,2); }
    // if (prefix == "PTall") { return 0.163228+-0.0367291*x+0.804498*std::pow(x,2); }

    // Determine if the prefix is a one-dimensional case
    bool isAllPrefix = (prefix.find("all") != std::string::npos);
    std::string basePrefix = isAllPrefix ? prefix.substr(0, prefix.size() - 3) : prefix;

    if (basePrefix == "Q2" || basePrefix == "x" || basePrefix == "y" || basePrefix == "z" || basePrefix == "zeta" || basePrefix == "PT" || basePrefix == "xF" || basePrefix == "Mx") {
        std::string fourD_prefix = determine_4D_bin(Q2, x, z);  // Get the 4D bin prefix
        if (fourD_prefix.empty()) {
            std::cerr << "Error: No matching 4D bin found for the given Q2, x, and z values." << std::endl;
            return 0.19;  // Fallback value
        }

        // Call the appropriate 4D bin function
        if (isAllPrefix) {
            return all_4D_bin(fourD_prefix, rand_gen);
        } else {
            return standard_4D_bin(fourD_prefix, rand_gen);
        }
    }

    // If the prefix was already a 4D bin, call the appropriate function directly
    if (isAllPrefix) {
        return all_4D_bin(prefix, rand_gen);
    } else {
        return standard_4D_bin(prefix, rand_gen);
    }
}

// double dilution_factor(double currentVariable, const std::string& pref0

//   TRandom3 rand_gen;

//   if (prefix == "Q2") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
//   if (prefix == "x") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
//   if (prefix == "y") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
//   if (prefix == "z") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
//   if (prefix == "zeta") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
//   if (prefix == "PT") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
//   if (prefix == "xF") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }
//   if (prefix == "Mx") { return 0.136288+0.283353*currentVariable+-0.125745*std::pow(currentVariable,2); }

//   if (prefix == "Q2all") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
//   if (prefix == "xall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
//   if (prefix == "yall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
//   if (prefix == "zall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
//   if (prefix == "zetaall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
//   if (prefix == "PTall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }
//   if (prefix == "xFall") { return 0.163228+-0.0367291*currentVariable+0.804498*std::pow(currentVariable,2); }


//   if (prefix == "Q2x1z1") { double sigma = 0.00941327; return 0.149075 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x1z2") { double sigma = 0.0132642; return 0.162577 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x1z3") { double sigma = 0.0229451; return 0.141929 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x1z4") { double sigma = 0.0665136; return 0.181142 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x2z1") { double sigma = 0.00902454; return 0.141377 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x2z2") { double sigma = 0.00594617; return 0.15723 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x2z3") { double sigma = 0.0121427; return 0.164111 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x2z4") { double sigma = 0.00819397; return 0.184566 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x3z1") { double sigma = 0.0111792; return 0.159587 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x3z2") { double sigma = 0.00459304; return 0.180858 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x3z3") { double sigma = 0.0104746; return 0.17691 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x3z4") { double sigma = 0.0159035; return 0.160347 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x4z1") { double sigma = 0.00977397; return 0.175046 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x4z2") { double sigma = 0.00518098; return 0.186986 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x4z3") { double sigma = 0.0110679; return 0.182736 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x4z4") { double sigma = 0.0221351; return 0.176834 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x5z1") { double sigma = 0.146381; return 0.348494 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x5z2") { double sigma = 0.0108159; return 0.154464 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x5z3") { double sigma = 0.0123293; return 0.168112 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x5z4") { double sigma = 0.0113772; return 0.164353 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x6z1") { double sigma = 0.0219882; return 0.132161 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x6z2") { double sigma = 0.0106923; return 0.160582 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x6z3") { double sigma = 0.00607725; return 0.167237 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x6z4") { double sigma = 0.0110221; return 0.145863 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x7z1") { double sigma = 0.00527885; return 0.141097 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x7z2") { double sigma = 0.00544389; return 0.181025 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x7z3") { double sigma = 0.00507036; return 0.189089 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x7z4") { double sigma = 0.00954266; return 0.186328 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x8z1") { double sigma = 0.00882734; return 0.168493 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x8z2") { double sigma = 0.00430375; return 0.199761 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x8z3") { double sigma = 0.00440225; return 0.192723 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x8z4") { double sigma = 0.00854078; return 0.198004 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x9z2") { double sigma = 0.00738543; return 0.184486 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x9z3") { double sigma = 0.00710463; return 0.195825 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x9z4") { double sigma = 0.00505832; return 0.197091 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x10z1") { double sigma = 0.042577; return 0.0624643 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x10z2") { double sigma = 0.00775736; return 0.198887 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x10z3") { double sigma = 0.00395318; return 0.212409 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x10z4") { double sigma = 0.0023216; return 0.22681 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x11z1") { double sigma = 0.00856085; return 0.165176 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x11z2") { double sigma = 0.00353569; return 0.202314 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x11z3") { double sigma = 0.00314671; return 0.214343 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x11z4") { double sigma = 0.00615722; return 0.208119 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x12z2") { double sigma = 0.009042; return 0.207689 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x12z3") { double sigma = 0.00594662; return 0.223018 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x12z4") { double sigma = 0.00583763; return 0.242526 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x13z1") { double sigma = 0.00990338; return 0.174887 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x13z2") { double sigma = 0.00306666; return 0.218621 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x13z3") { double sigma = 0.00275945; return 0.240082 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x13z4") { double sigma = 0.00645786; return 0.241681 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x14z1") { double sigma = 0.0253505; return 0.115816 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x14z2") { double sigma = 0.0123555; return 0.246474 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x14z3") { double sigma = 0.00707888; return 0.266938 + rand_gen.Gaus(0, sigma); }

//   if (prefix == "Q2x14z4") { double sigma = 0.00480178; return 0.261037 + rand_gen.Gaus(0, sigma); }
  
//   return 0.19;
// }