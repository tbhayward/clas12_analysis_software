#ifndef HISTCONFIGS_H
#define HISTCONFIGS_H

#include <map>
#include <string>

struct HistConfig {
    int bins;
    double min;
    double max;
};

// Global variable declaration
extern std::map<std::string, HistConfig> histConfigs;

#endif