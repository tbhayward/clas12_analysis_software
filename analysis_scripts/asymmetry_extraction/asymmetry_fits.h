#ifndef ASYMMETRY_FITS_H
#define ASYMMETRY_FITS_H

#include <string>
#include <cmath>
#include <iostream>

double asymmetry_value_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, double Ptp, double Ptm, 
  int asymmetry_index);
double asymmetry_error_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, double Ptp, double Ptm, 
  int asymmetry_index);
double BSA_funcToFit(double* x, double* par);
double TSA_funcToFit(double* x, double* par);
double DSA_funcToFit(double* x, double* par);

#endif // ASYMMETRY_FITS_H
