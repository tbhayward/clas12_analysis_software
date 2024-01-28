#pragma once

#include <map>
#include <string>
#include <vector>
#include "histConfigs.h"

extern std::map<std::string, std::vector<double>> bins_map;
extern std::vector<std::vector<double>> allBins;
extern std::vector<std::string> binNames;
extern std::vector<std::string> propertyNames;
extern std::vector<std::string> variable_names;
extern double cmm;
extern double cpm;
extern double cmp;
extern double cpp;
extern double total_charge_carbon;
extern int currentFits;
extern int channel; 

extern std::map<std::string, HistConfig> histConfigs;