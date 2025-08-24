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
double BSA_dihadron(double* x, double* par);
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

// TSA(phi) = ( AUL*sin(phi) + AUL2*sin(2phi) + A_T_UL*sTG*sin(phi) ) / GE_den(...)
// NOTE: for the leakage we pass a *centered* sTG per φ-bin from the fit driver.
double GE_model_TSA_with_sTG(double phi, double AUL, double AUL2, double A_T_UL, double sTG,
                             double B_UUcos, double C_UUcos2);

// DSA(phi) = ( ALL + ALLcos*cos(phi) + A_T_LL*sTG*sin(phi) ) / GE_den(...)
double GE_model_DSA_with_sTG(double phi, double ALL, double ALLcos, double A_T_LL, double sTG,
                             double B_UUcos, double C_UUcos2);

// ----- TF1 wrappers -----
// x[0] = phi
// BSA params: [0]=ALU, [1]=B_UUcos, [2]=C_UUcos2
double BSA_general_exclusive_TF1(double* x, double* par);

// TSA params: [0]=AUL, [1]=AUL2, [2]=A_T_UL, [3]=(IGNORED), [4]=B_UUcos, [5]=C_UUcos2
// The centered sTG(φ-bin) is injected inside the fit; par[3] is ignored for overlays.
double TSA_general_exclusive_TF1(double* x, double* par);

// DSA params: [0]=ALL, [1]=ALLcos, [2]=A_T_LL, [3]=(IGNORED), [4]=B_UUcos, [5]=C_UUcos2
double DSA_general_exclusive_TF1(double* x, double* par);

#endif // ASYMMETRY_FITS_H