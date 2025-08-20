#ifndef ASYMMETRY_FITS_H
#define ASYMMETRY_FITS_H

#include <string>
#include <cmath>
#include <iostream>

// Existing declarations…
double asymmetry_value_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, 
  double Ptp, double Ptm, int asymmetry_index);
double asymmetry_error_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, 
  double Ptp, double Ptm, int asymmetry_index);
double Legendre_P(int ell, int m, float theta);
double BSA_inclusive(double* par);
double BSA_single_hadron(double* x, double* par);
double BSA_b2b_dihadron(double* x, double* par);
double BSA_dihadron(double* x,double* par);
double BSA_dvcs(double* x, double* par);
double TSA_inclusive(double* par);
double TSA_single_hadron(double* x, double* par);
double TSA_b2b_dihadron(double* x, double* par);
double TSA_dihadron(double* x, double* par);
double DSA_single_hadron(double* x, double* par);
double DSA_inclusive(double* par);
double DSA_b2b_dihadron(double* x, double* par);
double DSA_dihadron(double* x, double* par);

// ===== GeneralExclusive model helpers =====
// Denominator: 1 + B*cos(phi) + C*cos(2phi)
double GE_den(double phi, double B_UUcos, double C_UUcos2);

// BSA(phi) = ( ALU * sin(phi) ) / GE_den(...)
double GE_model_BSA(double phi, double ALU,
                    double B_UUcos, double C_UUcos2);

// TSA(phi) = ( AUL*sin(phi) + AUL2*sin(2phi) + Atg*sTG*sin(phi) ) / GE_den(...)
double GE_model_TSA_with_sTG(double phi, double AUL, double AUL2, double Atg, double sTG,
                             double B_UUcos, double C_UUcos2);

// DSA(phi) = ( ALL + ALL2*cos(phi) ) / GE_den(...)
double GE_model_DSA(double phi, double ALL, double ALL2,
                    double B_UUcos, double C_UUcos2);

// ----- TF1 wrappers -----
// x[0] = phi
// BSA params: [0]=ALU, [1]=B_UUcos, [2]=C_UUcos2
double BSA_general_exclusive_TF1(double* x, double* par);

// TSA params: [0]=AUL, [1]=AUL2, [2]=Atg, [3]=sTG(const), [4]=B_UUcos, [5]=C_UUcos2
double TSA_general_exclusive_TF1(double* x, double* par);

// DSA params: [0]=ALL, [1]=ALL2, [2]=B_UUcos, [3]=C_UUcos2
double DSA_general_exclusive_TF1(double* x, double* par);

#endif // ASYMMETRY_FITS_H