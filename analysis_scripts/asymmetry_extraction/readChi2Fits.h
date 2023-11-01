#ifndef readChi2Fits_H
#define readChi2Fits_H

#include <string>
#include <vector>
#include <map>

std::map<std::string, std::vector<std::vector<double>>> readChi2Fits(const std::string& filepath);

#endif // readChi2Fits_H
