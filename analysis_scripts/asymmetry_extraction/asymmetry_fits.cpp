#include "common_vars.h"
#include "asymmetry_fits.h"
#include "dilution_factor.h" 

// Function to fit the beam-spin asymmetry histogram
double BSA_funcToFit(double* x, double* par) {
  // Retrieve the parameters 
  double ALU_offset = par[0];
  double ALU_sinphi = par[1];
  // double AUU_cosphi = par[2];
  // double AUU_cos2phi = par[3];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALU_offset + ALU_sinphi*sin(phi);
  // return (ALU_offset + ALU_sinphi*sin(phi)) / (1 + AUU_cosphi*cos(phi) + AUU_cos2phi*cos(2*phi));
}

// Function to fit the target-spin asymmetry histogram
double TSA_funcToFit(double* x, double* par) {
  // Retrieve the parameters A
  double AUL_offset = par[0];
  double AUL_sinphi = par[1];
  double AUL_sin2phi = par[2];
  // double AUU_cosphi = par[3];
  // double AUU_cos2phi = par[4];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return AUL_offset + AUL_sinphi*sin(phi)+AUL_sin2phi*sin(2*phi);
  // return (AUL_offset + AUL_sinphi*sin(phi)+AUL_sin2phi*sin(2*phi)) /
  //   (1 + AUU_cosphi*cos(phi) + AUU_cos2phi*cos(2*phi));
}

// Function to fit the double-spin asymmetry histogram
double DSA_funcToFit(double* x, double* par) {
  // Retrieve the parameters A
  double ALL = par[0];
  double ALL_cosphi = par[1];
  // double AUU_cosphi = par[2];
  // double AUU_cos2phi = par[3];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALL+ALL_cosphi*cos(phi);
  // return (ALL+ALL_cosphi*cos(phi)) / (1 + AUU_cosphi*cos(phi) + AUU_cos2phi*cos(2*phi));
}

double asymmetry_value_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, double Ptp, double Ptm, 
  int asymmetry_index) {
  double Df = dilution_factor(currentVariable, prefix); // dilution factor
  // return the asymmetry value 
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      return (1 / meanPol) * (Ptm*(Npp-Nmp)+Ptp*(Npm-Nmm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    case 1: // target-spin asymmetry
      return (1 / Df) * ((Npp+Nmp)-(Npm+Nmm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    case 2: // double-spin asymmetry
      return (1 / (Df*meanPol)) * ((Npp-Nmp)+(Nmm-Npm)) / (Ptm*(Npp+Nmp)+Ptp*(Npm+Nmm));
    default:
      std::cout << "Invalid asymmetry_index!" << std::endl;
      return 0;
  }
}

double asymmetry_error_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, double Ptp, double Ptm, 
  int asymmetry_index) {
  double Df = dilution_factor(currentVariable, prefix); // dilution factor
  // return the asymmetry error 
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry
      return (2 / meanPol) * std::sqrt(
        ((cmm*cpm*cpp*Nmp*std::pow(Ptm,2)*std::pow(Npp*Ptm+Npm*Ptp,2))+
        (cmp*cpm*cpp*Nmm*std::pow(Ptp,2)*std::pow(Npp*Ptm+Npm*Ptp,2))+
        (cmm*cmp*std::pow(Nmp*Ptm+Nmm*Ptp,2)*(cpm*Npp*std::pow(Ptm,2)+cpp*Npm*std::pow(Ptp,2))))/
        (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
    case 1: // target-spin asymmetry
      return (1 / Df) * std::sqrt(
        (((cmp*cpm*cpp*Nmm*std::pow(Nmp+Npp,2)+cmm*cmp*cpp*Npm*std::pow(Nmp+Npp,2)+
        cmm*cpm*std::pow(Nmm+Npm,2)*(cpp*Nmp+cmp*Npp))*std::pow(Ptm+Ptp,2))) /
        (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
    case 2: // double-spin asymmetry
      return (1 / (Df*meanPol)) * std::sqrt(
        (cmp*cpm*cpp*Nmm*std::pow((Nmp+Npp)*Ptm+(Nmp+2*Npm-Npp)*Ptp,2) + 
        cmm*cmp*cpp*Npm*std::pow(Nmp*(Ptm-Ptp)+2*Nmm*Ptp+Npp*(Ptm+Ptp),2) +
        cmm*cpm*(cmp*Npp*std::pow((-Nmm+2*Nmp+Npm)*Ptm+(Nmm+Npm)*Ptp,2) +
        cpp*Nmp*std::pow((Nmm-Npm+2*Npp)*Ptm+(Nmm+Npm)*Ptp,2))) / 
        (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
    default:
      std::cout << "Invalid asymmetry_index!" << std::endl;
      return 0;
  }
}