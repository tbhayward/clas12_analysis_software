#ifndef DETERMINE_EXCLUSIVITY_H
#define DETERMINE_EXCLUSIVITY_H

#include <TTreeReader.h>
#include <string>

void determine_exclusivity(const std::string& analysisType, const std::string& topology, TTreeReader& dataReader, TTreeReader& mcReader, const std::string& outputDir, const std::string& plotTitle);

#endif