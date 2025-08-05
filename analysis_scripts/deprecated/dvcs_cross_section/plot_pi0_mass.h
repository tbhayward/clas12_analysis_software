#ifndef PLOT_PI0_MASS_H
#define PLOT_PI0_MASS_H

#include <string>
#include <TTreeReader.h>

// Function to plot the invariant mass of two photons
void plot_pi0_mass(TTreeReader& dataReader1, TTreeReader& dataReader2, TTreeReader& dataReader3,
                   TTreeReader& mcReader1, TTreeReader& mcReader2, TTreeReader& mcReader3,
                   const std::string& outputDir);

#endif  // PLOT_PI0_MASS_H