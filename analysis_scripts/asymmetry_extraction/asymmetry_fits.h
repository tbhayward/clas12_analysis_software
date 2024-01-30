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
double Legendre_P(int ell, int m, float theta);
double BSA_inclusive(double* par);
double BSA_single_hadron(double* x, double* par);
double BSA_b2b_dihadron(double* x, double* y, double* par);
double BSA_dihadron(double* x, double* y, double* z, double* par);
double TSA_inclusive(double* par);
double TSA_single_hadron(double* x, double* par);
double TSA_b2b_dihadron(double* x, double* y, double* par);
double TSA_dihadron(double* x, double* y, double* z, double* par);
double DSA_single_hadron(double* x, double* par);
double DSA_inclusive(double* par);
double BSA_b2b_dihadron(double* x, double* y, double* par);
double DSA_dihadron(double* x, double* y, double* z, double* par);

#endif // ASYMMETRY_FITS_H
