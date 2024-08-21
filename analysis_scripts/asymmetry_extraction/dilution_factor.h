#pragma once
#include <string>
#include <TRandom3.h>

// Function to determine the appropriate 4D bin prefix
std::string determine_4D_bin(double Q2, double x, double z);

// Function to calculate the dilution factor
double dilution_factor(double Q2, double x, double z, double pT, const std::string& prefix);

// Function for standard 4D bin dilution factor calculation
double standard_4D_bin(const std::string& prefix, TRandom3& rand_gen);

// Function for "all" 4D bin dilution factor calculation
double all_4D_bin(const std::string& prefix, TRandom3& rand_gen);