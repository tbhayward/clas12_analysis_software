#ifndef DETERMINE_EXCLUSIVITY_H
#define DETERMINE_EXCLUSIVITY_H

#include <TTreeReader.h>

void determine_exclusivity(TTreeReader& dataReader, TTreeReader& mcReader, const std::string& outputDir);

#endif