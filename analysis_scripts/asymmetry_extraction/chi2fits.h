#ifndef CHI2FITS_H
#define CHI2FITS_H

#include <string>
#include <cmath>
#include <iostream>

double asymmetry_value_calculation(double currentVariable, const std::string& prefix, ...);
double asymmetry_error_calculation(double currentVariable, const std::string& prefix, ...);
double BSA_funcToFit(double* x, double* par);
double TSA_funcToFit(double* x, double* par);
double DSA_funcToFit(double* x, double* par);

#endif // CHI2FITS_H
