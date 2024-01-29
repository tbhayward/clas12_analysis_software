#include "common_vars.h"
#include "asymmetry_fits.h"
#include "dilution_factor.h" 

/******** Legendre Polynomial ********/
double Legendre_P(int ell, int m, float theta) {
  if (ell == 0 && m == 0) return 1;
  if (ell == 1 && m == 0) return cos(theta);
  if ((ell == 1 && m == 1) || (ell == 1 && m == -1)) return sin(theta);
  if (ell == 2 && m == 0) return 0.5*(3*cos(theta)*cos(theta)-1);
  if ((ell == 2 && m == 1) || (ell == 2 && m == -1)) return sin(2*theta);
  if ((ell == 2 && m == 2) || (ell == 2 && m == -2)) return sin(theta)*sin(theta);
  if (ell == 3 && m == 0) return 0.5*cos(theta)*(5*cos(theta)*cos(theta)-3);
  if ((ell == 3 && m == 1) || (ell == 3 && m == -1)) return (5*cos(theta)*cos(theta)-1)*sin(theta);
  if ((ell == 3 && m == 2) || (ell == 3 && m == -2)) return sin(2*theta)*sin(theta);
  if ((ell == 3 && m == 3) || (ell == 3 && m == -3)) return sin(theta)*sin(theta)*sin(theta);
  return 0;
}



/******** BEAM-SPIN ASYMMETRY ********/
// Function to fit the beam-spin asymmetry histogram

double BSA_inclusive(double* par) {
  // Retrieve the parameters 
  double ALU_offset = par[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALU_offset;
}

double BSA_single_hadron(double* x, double* par) {
  // Retrieve the parameters 
  double ALU_offset = par[0];
  double ALU_sinphi = par[1];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALU_offset + ALU_sinphi*sin(phi);
}

double BSA_b2b_dihadron(double* x, double* y, double* par) {
  // Retrieve the parameters 
  double ALU_offset = par[0];
  double ALU_sinphi1 = par[1];
  double ALU_sinphi2 = par[2];
  double ALU_sinDeltaphi = par[3];
  // Retrieve the phi variables from the input arrays
  double phi1 = x[0];
  double phi2 = y[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALU_offset + ALU_sinphi1*sin(phi1) + ALU_sinphi2*sin(phi2) +
    ALU_sinDeltaphi*sin(phi1 - phi2);
}

double BSA_dihadron(double* x, double* y, double* z, double* par) {
  // Retrieve the parameters
  double ALU_offset = par[0];
  double ALU_W_ell0_m0 = par[1];
  double ALU_C_ell1_m1 = par[2];
  double ALU_W_ell1_mn1 = par[3];
  double ALU_W_ell1_m0 = par[4];
  double ALU_W_ell1_m1 = par[5];
  double ALU_C_ell2_m1 = par[6];
  double ALU_C_ell2_m2 = par[7];
  double ALU_W_ell2_mn2 = par[8];
  double ALU_W_ell2_mn1 = par[9];
  double ALU_W_ell2_m0 = par[10];
  double ALU_W_ell2_m1 = par[11];
  double ALU_W_ell2_m2 = par[12];

  // Retrieve the angle variables from the input arrays
  double phih = x[0];
  double phiR = y[0];
  double theta = z[0];
  // Calculate and return the value of the function for the given phi and parameters
  return ALU_offset + 
    ALU_W_ell0_m0*Legendre_P(0,0,theta)*sin(phih) +             // tw3, ell=0, m=0
    ALU_C_ell1_m1*2*Legendre_P(1,1,theta)*sin(phih - phiR) +    // tw2, ell=1, m=1
    ALU_W_ell1_mn1*Legendre_P(1,-1,theta)*sin(2*phih-phiR) +    // tw3, ell=1, m=-1
    ALU_W_ell1_m0*Legendre_P(1,0,theta)*sin(phih) +             // tw3, ell=1, m=0
    ALU_W_ell1_m1*Legendre_P(1,1,theta)*sin(phiR) +             // tw3, ell=1, m=1
    ALU_C_ell2_m1*2*Legendre_P(2,1,theta)*sin(phih-phiR) +      // tw2, ell=2, m=1
    ALU_C_ell2_m2*2*Legendre_P(2,2,theta)*sin(2*phih-2*phiR) +  // tw2, ell=2, m=2
    ALU_W_ell2_mn2*Legendre_P(2,-2,theta)*sin(3*phih-2*phiR) +  // tw3, ell=2, m=-2
    ALU_W_ell2_mn1*Legendre_P(2,-1,theta)*sin(2*phih-phiR) +    // tw3, ell=2, m=-1
    ALU_W_ell2_m0*Legendre_P(2,0,theta)*sin(phih) +             // tw3, ell=2, m=0
    ALU_W_ell2_m1*Legendre_P(2,1,theta)*sin(phiR) +             // tw3, ell=2, m=1
    ALU_W_ell2_m2*Legendre_P(2,2,theta)*sin(-phih+2*phiR);      // tw3, ell=2, m=2
}

/******** TARGET-SPIN ASYMMETRY ********/

double TSA_inclusive(double* par) {
  // Retrieve the parameters A
  double AUL_offset = par[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return AUL_offset;
}

double TSA_single_hadron(double* x, double* par) {
  // Retrieve the parameters A
  double AUL_offset = par[0];
  double AUL_sinphi = par[1];
  double AUL_sin2phi = par[2];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return AUL_offset + AUL_sinphi*sin(phi)+AUL_sin2phi*sin(2*phi);
}

double TSA_b2b_dihadron(double* x, double* y, double* par) {
  // Retrieve the parameters 
  double AUL_offset = par[0];
  double AUL_sinphi1 = par[1];
  double AUL_sinphi2 = par[2];
  double AUL_sin2phi1 = par[3];
  double AUL_sin2phi2 = par[4];
  double AUL_sinDeltaphi = par[5];
  double AUL_sinSumphi = par[6];
  // Retrieve the phi variable from the input x array
  double phi1 = x[0];
  double phi2 = y[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return AUL_offset + AUL_sinphi1*sin(phi1) + AUL_sinphi2*sin(phi2) +
    AUL_sin2phi1*sin(2*phi1) + AUL_sin2phi2*sin(2*phi2) +
    + AUL_sinDeltaphi*sin(phi1 - phi2) + AUL_sinSumphi*sin(phi1 + phi2);
}

double TSA_dihadron(double* x, double* y, double* z, double* par) {
  // Retrieve the parameters
  double ALU_offset = par[0];
  double ALU_B_ell0_m0 = par[1];
  double ALU_V_ell0_m0 = par[2];
  double ALU_A_ell1_m1 = par[3];
  double ALU_B_ell1_mn1 = par[4];
  double ALU_B_ell1_m0 = par[5];
  double ALU_B_ell1_m1 = par[6];
  double ALU_V_ell1_mn1 = par[7];
  double ALU_V_ell1_m0 = par[8];
  double ALU_V_ell1_m1 = par[9];
  double ALU_A_ell2_m1 = par[10];
  double ALU_A_ell2_m2 = par[11];
  double ALU_B_ell2_mn2 = par[12];
  double ALU_B_ell2_mn1 = par[13];
  double ALU_B_ell2_m0 = par[14];
  double ALU_B_ell2_m1 = par[15];
  double ALU_B_ell2_m2 = par[16];
  double ALU_V_ell2_mn2 = par[17];
  double ALU_V_ell2_mn1 = par[18];
  double ALU_V_ell2_m0 = par[19];
  double ALU_V_ell2_m1 = par[20]; 
  double ALU_V_ell2_m2 = par[21];

  // Retrieve the angle variables from the input arrays
  double phih = x[0];
  double phiR = y[0];
  double theta = z[0];
  // Calculate and return the value of the function for the given phi and parameters
  return ALU_offset + 
    ALU_B_ell0_m0*Legendre_P(0,0,theta)*sin(2*phih) +             // tw2, B, ell=0, m=0
    ALU_V_ell0_m0*Legendre_P(0,0,theta)*sin(phih) +               // tw3, V, ell=0, m=0
    ALU_A_ell1_m1*Legendre_P(1,1,theta)*sin(-phih+phiR) +         // tw2, A, ell=1, m=1
    ALU_B_ell1_mn1*Legendre_P(1,-1,theta)*sin(3*phih-phiR) +      // tw2, B, ell=1, m=-1
    ALU_B_ell1_m0*Legendre_P(1,0,theta)*sin(2*phih) +             // tw2, B, ell=1, m=0
    ALU_B_ell1_m1*Legendre_P(1,1,theta)*sin(phih+phiR) +          // tw2, B, ell=1, m=1
    ALU_V_ell1_mn1*Legendre_P(1,-1,theta)*sin(2*phih-phiR) +      // tw3, V, ell=1, m=-1
    ALU_V_ell1_m0*Legendre_P(1,0,theta)*sin(phih) +               // tw3, V, ell=1, m=0
    ALU_V_ell1_m1*Legendre_P(1,1,theta)*sin(phiR) +               // tw3, V, ell=1, m=1
    ALU_A_ell2_m1*Legendre_P(2,1,theta)*sin(-phih+phiR) +         // tw2, A, ell=2, m=1
    ALU_A_ell2_m2*Legendre_P(2,2,theta)*sin(-2*phih+2*phiR) +     // tw2, A, ell=2, m=2
    ALU_B_ell2_mn2*Legendre_P(2,-2,theta)*sin(4*phih-2*phiR) +    // tw2, B, ell=2, m=-2
    ALU_B_ell2_mn1*Legendre_P(2,-1,theta)*sin(3*phih-phiR) +      // tw2, B, ell=2, m=-1
    ALU_B_ell2_m0*Legendre_P(2,0,theta)*sin(2*phih) +             // tw2, B, ell=2, m=0
    ALU_B_ell2_m1*Legendre_P(2,1,theta)*sin(phih+phiR) +          // tw2, B, ell=2, m=1
    ALU_B_ell2_m2*Legendre_P(2,2,theta)*sin(2*phiR) +             // tw2, B, ell=2, m=2
    ALU_V_ell2_mn2*Legendre_P(2,-2,theta)*sin(3*phih-2*phiR) +    // tw3, V, ell=2, m=-2
    ALU_V_ell2_mn1*Legendre_P(2,-1,theta)*sin(2*phih-phiR) +      // tw3, V, ell=2, m=-1
    ALU_V_ell2_m0*Legendre_P(2,0,theta)*sin(phih) +               // tw3, V, ell=2, m=0
    ALU_V_ell2_m1*Legendre_P(2,1,theta)*sin(phiR) +               // tw3, V, ell=2, m=1
    ALU_V_ell2_m2*Legendre_P(2,2,theta)*sin(-phih+2*phiR);        // tw3, V, ell=2, m=2
}

/******** DOUBLE-SPIN ASYMMETRY ********/

double DSA_inclusive(double* par) {
  // Retrieve the parameters A
  double ALL = par[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALL;
}

double DSA_single_hadron(double* x, double* par) {
  // Retrieve the parameters A
  double ALL = par[0];
  double ALL_cosphi = par[1];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALL + ALL_cosphi*cos(phi);
}

double DSA_single_hadron(double* x, double* y, double* par) {
  // Retrieve the parameters A
  double ALL = par[0];
  double ALL_cosphi1 = par[1];
  double ALL_cosphi2 = par[2];
  // Retrieve the phi variable from the input x array
  double phi1 = x[0];
  double phi2 = y[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALL + ALL_cosphi1*cos(phi1) + ALL_cosphi2*cos(phi2);
}

double DSA_dihadron(double* x, double* y, double* z, double* par) {
  // Retrieve the parameters
  double ALU_offset = par[0];
  double ALU_C_ell0_m0 = par[1];
  double ALU_W_ell0_m0 = par[2];
  double ALU_C_ell1_m0 = par[3];
  double ALU_C_ell1_m1 = par[4];
  double ALU_W_ell1_mn1 = par[5];
  double ALU_W_ell1_m0 = par[6];
  double ALU_W_ell1_m1 = par[7];
  double ALU_C_ell2_m0 = par[8];
  double ALU_C_ell2_m1 = par[9];
  double ALU_C_ell2_m2 = par[10];
  double ALU_W_ell2_mn2 = par[11];
  double ALU_W_ell2_mn1 = par[12];
  double ALU_W_ell2_m0 = par[13];
  double ALU_W_ell2_m1 = par[14];
  double ALU_W_ell2_m2 = par[15];

  // Retrieve the angle variables from the input arrays
  double phih = x[0];
  double phiR = y[0];
  double theta = z[0];
  // Calculate and return the value of the function for the given phi and parameters
  return ALU_offset + 
    ALU_C_ell0_m0*2*Legendre_P(0,0,theta)*cos(0) +            // tw2, ell=0, m=0
    ALU_W_ell0_m0*Legendre_P(0,0,theta)*cos(phih) +           // tw3, ell=0, m=0
    ALU_C_ell1_m0*2*Legendre_P(1,0,theta)*cos(0) +            // tw2, ell=1, m=0
    ALU_C_ell1_m1*4*Legendre_P(1,1,theta)*cos(phih-phiR) +    // tw2, ell=1, m=1
    ALU_W_ell1_mn1*Legendre_P(1,-1,theta)*cos(2*phih-phiR) +  // tw3, ell=1, m=-1
    ALU_W_ell1_m0*Legendre_P(1,0,theta)*cos(phih) +           // tw3, ell=1, m=0
    ALU_W_ell1_m1*Legendre_P(1,1,theta)*cos(phiR) +           // tw3, ell=1, m=1
    ALU_C_ell2_m0*2*Legendre_P(2,0,theta)*cos(0) +            // tw2, ell=2, m=0
    ALU_C_ell2_m1*4*Legendre_P(2,1,theta)*cos(phih-phiR) +    // tw2, ell=2, m=1
    ALU_C_ell2_m2*4*Legendre_P(2,2,theta)*cos(2*phih-2*phiR) +// tw2, ell=2, m=2
    ALU_W_ell2_mn2*Legendre_P(2,-2,theta)*cos(3*phih-2*phiR) +// tw3, ell=2, m=-2
    ALU_W_ell2_mn1*Legendre_P(2,-1,theta)*cos(2*phih-phiR) +  // tw3, ell=2, m=-1
    ALU_W_ell2_m0*Legendre_P(2,0,theta)*cos(phih) +           // tw3, ell=2, m=0
    ALU_W_ell2_m1*Legendre_P(2,1,theta)*cos(phiR) +           // tw3, ell=2, m=1
    ALU_W_ell2_m2*Legendre_P(2,2,theta)*cos(-phih+2*phiR);    // tw3, ell=2, m=2
}


/******** VALUE CALCULATIONS ********/

double asymmetry_value_calculation(double currentVariable, const std::string& prefix, 
  double Npp, double Npm, double Nmp, double Nmm, double meanPol, double Ptp, double Ptm, 
  int asymmetry_index) {
  double Df = dilution_factor(currentVariable, prefix); // dilution factor
  std::cout << Df << std::endl;
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