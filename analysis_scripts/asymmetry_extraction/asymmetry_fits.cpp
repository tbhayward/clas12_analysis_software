#include "common_vars.h"
#include "asymmetry_fits.h"
#include "dilution_factor.h" 
#include <cmath>

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

double BSA_dvcs(double* x, double* par) {
  // Retrieve the parameters 
  double ALU_offset = par[0];
  double ALU_sinphi = par[1];
  double AUU_cosphi = par[2];
  // Retrieve the phi variable from the input x array
  double phi = x[0];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALU_offset + ALU_sinphi*sin(phi)/(1+AUU_cosphi*cos(phi));
}

double BSA_b2b_dihadron(double* x, double* par) {
  // Retrieve the parameters 
  double ALU_offset = par[0];
  double ALU_sinphi1 = par[1];
  double ALU_sinphi2 = par[2];
  double ALU_sinDeltaphi = par[3];
  double ALU_sin2Deltaphi = par[4];
  // Retrieve the phi variables from the input arrays
  double phi1 = x[0];
  double phi2 = x[1];
  // Calculate and return the value of the function for the given phi and parameters 
  return ALU_offset + ALU_sinphi1*sin(phi1) + ALU_sinphi2*sin(phi2) +
    ALU_sinDeltaphi*sin(phi1 - phi2) + ALU_sin2Deltaphi*sin(2*phi1 - 2*phi2);
}

double BSA_dihadron(double* x, double* par) {
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
  double phiR = x[1];
  double theta = x[2];
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

double TSA_b2b_dihadron(double* x, double* par) {
  // Retrieve the parameters 
  double AUL_offset = par[0];
  double AUL_sinphi1 = par[1];
  double AUL_sinphi2 = par[2];
  double AUL_sin2phi1 = par[3];
  double AUL_sin2phi2 = par[4];
  double AUL_sinDeltaphi = par[5];
  double AUL_sin2Deltaphi = par[6];
  double AUL_sinSumphi = par[7];
  // Retrieve the phi variable from the input x array
  double phi1 = x[0];
  double phi2 = x[1];
  // Calculate and return the value of the function for the given phi and parameters 
  return AUL_offset + AUL_sinphi1*sin(phi1) + AUL_sinphi2*sin(phi2) +
    AUL_sin2phi1*sin(2*phi1) + AUL_sin2phi2*sin(2*phi2) +
    + AUL_sinDeltaphi*sin(phi1 - phi2) + AUL_sin2Deltaphi*sin(2*phi1 - 2*phi2) +
    AUL_sinSumphi*sin(phi1 + phi2);
}

double TSA_dihadron(double* x, double* par) {
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
  double phiR = x[1];
  double theta = x[2];
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
  // return ALL;
}

double DSA_b2b_dihadron(double* x, double* par) {
  // Retrieve the parameters A
  double ALL = par[0];
  double ALL_cosphi1 = par[1];
  double ALL_cosphi2 = par[2];
  // Retrieve the phi variable from the input x array
  double phi1 = x[0];
  double phi2 = x[1];
  // Calculate and return the value of the function for the given phi and parameters 
  // return ALL;
  return ALL + ALL_cosphi1*cos(phi1) + ALL_cosphi2*cos(phi2);
}

double DSA_dihadron(double* x, double* par) {
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
  double phiR = x[1];
  double theta = x[2];
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

double asymmetry_value_calculation(double currentVariable, const std::pair<double, double>& dilutionFactor,
  const std::string& prefix, double Npp, double Npm, double Nmp, double Nmm,
  double meanPol, double Ptp, double Ptm, 
  int asymmetry_index) {
  // double Df = dilution_factor(Q2, xB, z, pT, prefix); // dilution factor
  double Df = dilutionFactor.first;

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

double asymmetry_error_calculation(double currentVariable, const std::pair<double, double>& dilutionFactor, 
  const std::string& prefix, double npp, double npm, double nmp, double nmm,
  double Pb, double Ptp, double Ptm, 
  int asymmetry_index) {
  // double Df = dilution_factor(Q2, xB, z, pT, prefix); // dilution factor
  double Df = dilutionFactor.first;
  double sigmaDf = dilutionFactor.second;
  double sigmaPb = 0.015;
  double sigmaPtp = 0.025;
  double sigmaPtm = 0.025;

  // return the asymmetry error 
  switch (asymmetry_index) {
    case 0: // beam-spin asymmetry

      // return (2 / meanPol) * std::sqrt(
      //   ((cmm*cpm*cpp*Nmp*std::pow(Ptm,2)*std::pow(Npp*Ptm+Npm*Ptp,2))+
      //   (cmp*cpm*cpp*Nmm*std::pow(Ptp,2)*std::pow(Npp*Ptm+Npm*Ptp,2))+
      //   (cmm*cmp*std::pow(Nmp*Ptm+Nmm*Ptp,2)*(cpm*Npp*std::pow(Ptm,2)+cpp*Npm*std::pow(Ptp,2))))/
      //   (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
      return sqrt(
    (4 * cmm * pow(cmp, 3) * pow(cpm, 4) * pow(cpp, 4) * pow(nmm, 3) * nmp * Ptm * pow(Ptp, 3) * pow(sigmaPb, 2) +
     pow(cmp, 4) * pow(cpm, 4) * pow(cpp, 4) * pow(nmm, 4) * pow(Ptp, 4) * pow(sigmaPb, 2) -
     4 * pow(cmm, 3) * cmp * pow(cpm, 2) * pow(cpp, 2) * nmm * nmp *
       (-pow(cpm, 2) * pow(cpp, 2) * pow(nmp, 2) * pow(Ptm, 3) * Ptp * pow(sigmaPb, 2) +
        pow(cmp, 2) * (Ptp *
           (-2 * pow(Pb, 2) * Ptm * (pow(cpm, 2) * npp * pow(Ptm, 2) + pow(cpp, 2) * npm * pow(Ptp, 2)) +
            Ptm * pow((cpm * npp * Ptm + cpp * npm * Ptp), 2) * pow(sigmaPb, 2) +
            2 * cpm * cpp * npm * npp * pow(Pb, 2) * Ptp * pow(sigmaPtm, 2)) +
           2 * cpm * cpp * npm * npp * pow(Pb, 2) * pow(Ptm, 2) * pow(sigmaPtp, 2))) +
     pow(cmm, 4) *
       (pow(cpm, 4) * pow(cpp, 4) * pow(nmp, 4) * pow(Ptm, 4) * pow(sigmaPb, 2) +
        pow(cmp, 4) * pow((cpm * npp * Ptm + cpp * npm * Ptp), 4) * pow(sigmaPb, 2) +
        2 * pow(cmp, 2) * pow(cpm, 2) * pow(cpp, 2) * nmp *
          (2 * cpm * cpp * npm * npp * pow(Ptm, 3) * Ptp * (2 * pow(Pb, 2) - nmp * pow(sigmaPb, 2)) +
           pow(cpm, 2) * npp * pow(Ptm, 4) * (2 * (nmp + npp) * pow(Pb, 2) - nmp * npp * pow(sigmaPb, 2)) +
           pow(cpp, 2) * npm *
             (pow(Ptp, 2) *
                (pow(Ptm, 2) * (2 * (nmp + npm) * pow(Pb, 2) - nmp * npm * pow(sigmaPb, 2)) +
                 2 * nmp * npm * pow(Pb, 2) * pow(sigmaPtm, 2)) +
              2 * nmp * npm * pow(Pb, 2) * pow(Ptm, 2) * pow(sigmaPtp, 2)))) +
     2 * pow(cmm, 2) * pow(cmp, 2) * pow(cpm, 2) * pow(cpp, 2) * nmm *
       (3 * pow(cpm, 2) * pow(cpp, 2) * nmm * pow(nmp, 2) * pow(Ptm, 2) * pow(Ptp, 2) * pow(sigmaPb, 2) +
        pow(cmp, 2) *
          (2 * cpm * cpp * npm * npp * Ptm * pow(Ptp, 3) * (2 * pow(Pb, 2) - nmm * pow(sigmaPb, 2)) +
           pow(cpp, 2) * npm * pow(Ptp, 4) * (2 * (nmm + npm) * pow(Pb, 2) - nmm * npm * pow(sigmaPb, 2)) +
           pow(cpm, 2) * npp *
             (pow(Ptp, 2) *
                (pow(Ptm, 2) * (2 * (nmm + npp) * pow(Pb, 2) - nmm * npp * pow(sigmaPb, 2)) +
                 2 * nmm * npp * pow(Pb, 2) * pow(sigmaPtm, 2)) +
              2 * nmm * npp * pow(Pb, 2) * pow(Ptm, 2) * pow(sigmaPtp, 2))))) /
    pow((cmm * cpm * (cpp * nmp + cmp * npp) * Pb * Ptm + cmp * cpp * (cpm * nmm + cmm * npm) * Pb * Ptp), 4));
    case 1: // target-spin asymmetry
      // return (1 / Df) * std::sqrt(
      //   (((cmp*cpm*cpp*Nmm*std::pow(Nmp+Npp,2)+cmm*cmp*cpp*Npm*std::pow(Nmp+Npp,2)+
      //   cmm*cpm*std::pow(Nmm+Npm,2)*(cpp*Nmp+cmp*Npp))*std::pow(Ptm+Ptp,2))) /
      //   (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
      return sqrt(
    (pow(cmm, 2) * pow(cmp, 2) * pow(cpm, 2) * pow(cpp, 4) * pow(Df, 2) * nmp * 
     pow((cpm * nmm + cmm * npm), 2) * pow((Ptm + Ptp), 2) +
     pow(cmm, 2) * pow(cmp, 4) * pow(cpm, 2) * pow(cpp, 2) * pow(Df, 2) * 
     pow((cpm * nmm + cmm * npm), 2) * npp * pow((Ptm + Ptp), 2) +
     pow(cmm, 2) * pow(cmp, 2) * pow(cpm, 4) * pow(cpp, 2) * pow(Df, 2) * 
     nmm * pow((cpp * nmp + cmp * npp), 2) * pow((Ptm + Ptp), 2) +
     pow(cmm, 4) * pow(cmp, 2) * pow(cpm, 2) * pow(cpp, 2) * pow(Df, 2) * 
     npm * pow((cpp * nmp + cmp * npp), 2) * pow((Ptm + Ptp), 2) +
     pow((cmm * cpm * cpp * nmp - cmp * cpp * (cpm * nmm + cmm * npm) + 
          cmm * cmp * cpm * npp), 2) * 
     pow((cmm * cpm * (cpp * nmp + cmp * npp) * Ptm + 
          cmp * cpp * (cpm * nmm + cmm * npm) * Ptp), 2) * pow(sigmaDf, 2) +
     pow(cmm, 2) * pow(cpm, 2) * pow(Df, 2) * pow((cpp * nmp + cmp * npp), 2) * 
     pow((cmm * cpm * cpp * nmp - cmp * cpp * (cpm * nmm + cmm * npm) + 
          cmm * cmp * cpm * npp), 2) * pow(sigmaPtm, 2) +
     pow(cmp, 2) * pow(cpp, 2) * pow(Df, 2) * pow((cpm * nmm + cmm * npm), 2) * 
     pow((cmm * cpm * cpp * nmp - cmp * cpp * (cpm * nmm + cmm * npm) + 
          cmm * cmp * cpm * npp), 2) * pow(sigmaPtp, 2)) /
    (pow(Df, 4) * pow((cmm * cpm * (cpp * nmp + cmp * npp) * Ptm + 
                       cmp * cpp * (cpm * nmm + cmm * npm) * Ptp), 4)));
    case 2: // double-spin asymmetry
      // return (1 / (Df*meanPol)) * std::sqrt(
      //   (cmp*cpm*cpp*Nmm*std::pow((Nmp+Npp)*Ptm+(Nmp+2*Npm-Npp)*Ptp,2) + 
      //   cmm*cmp*cpp*Npm*std::pow(Nmp*(Ptm-Ptp)+2*Nmm*Ptp+Npp*(Ptm+Ptp),2) +
      //   cmm*cpm*(cmp*Npp*std::pow((-Nmm+2*Nmp+Npm)*Ptm+(Nmm+Npm)*Ptp,2) +
      //   cpp*Nmp*std::pow((Nmm-Npm+2*Npp)*Ptm+(Nmm+Npm)*Ptp,2))) / 
      //   (cmm*cmp*cpm*cpp*std::pow((Nmp+Npp)*Ptm+(Nmm+Npm)*Ptp,4)));
      return sqrt(
    (pow(cmp, 2) * pow(cpp, 2) * pow(Df, 2) * npm * pow(Pb, 2) * 
     pow((cmm * pow(cpm, 2) * (cpp * nmp + cmp * npp) * Ptm + 
          cmp * pow(cpm, 2) * cpp * nmm * Ptp + 
          pow(cmm, 2) * (cmp * cpp * nmm - cpm * cpp * nmp + 
          cmp * cpm * npp) * Ptp), 2) +
     pow(cmm, 2) * pow(cmp, 2) * pow(cpm, 2) * pow(cpp, 2) * pow(Df, 2) * 
     npp * pow(Pb, 2) * pow((cmp * cpm * (npm * Ptm + nmm * Ptp) + 
          cmm * (-cmp * nmm * Ptm + 2 * cpm * nmp * Ptm + 
          cmp * npm * Ptp)), 2) +
     pow(cmm, 2) * pow(cmp, 2) * pow(cpm, 2) * pow(cpp, 2) * pow(Df, 2) * 
     nmp * pow(Pb, 2) * pow((cpm * cpp * (-npm * Ptm + nmm * Ptp) + 
          cmm * (cpp * nmm * Ptm + 2 * cpm * npp * Ptm + 
          cpp * npm * Ptp)), 2) +
     pow(cmp, 2) * pow(cpp, 2) * pow(Df, 2) * nmm * pow(Pb, 2) * 
     pow((cmp * pow(cpm, 2) * cpp * npm * Ptp + 
          cmm * pow(cpm, 2) * (cpp * nmp - cmp * npp) * Ptp + 
          pow(cmm, 2) * (cpm * (cpp * nmp + cmp * npp) * Ptm + 
          cmp * cpp * npm * Ptp)), 2) +
     pow((cmp * cpm * cpp * npm - cmm * (cmp * cpp * nmm - 
          cpm * cpp * nmp + cmp * cpm * npp)), 2) * 
     pow(Pb, 2) * pow((cmm * cpm * (cpp * nmp + cmp * npp) * Ptm + 
          cmp * cpp * (cpm * nmm + cmm * npm) * Ptp), 2) * pow(sigmaDf, 2) +
     pow(Df, 2) * pow((cmp * cpm * cpp * npm - 
          cmm * (cmp * cpp * nmm - cpm * cpp * nmp + 
          cmp * cpm * npp)), 2) * 
     pow((cmm * cpm * (cpp * nmp + cmp * npp) * Ptm + 
          cmp * cpp * (cpm * nmm + cmm * npm) * Ptp), 2) * pow(sigmaPb, 2) +
     pow(cmm, 2) * pow(cpm, 2) * pow(Df, 2) * pow((cpp * nmp + cmp * npp), 2) * 
     pow((cmp * cpm * cpp * npm - 
          cmm * (cmp * cpp * nmm - cpm * cpp * nmp + 
          cmp * cpm * npp)), 2) * pow(Pb, 2) * pow(sigmaPtm, 2) +
     pow(cmp, 2) * pow(cpp, 2) * pow(Df, 2) * pow((cpm * nmm + cmm * npm), 2) * 
     pow((cmp * cpm * cpp * npm - 
          cmm * (cmp * cpp * nmm - cpm * cpp * nmp + 
          cmp * cpm * npp)), 2) * pow(Pb, 2) * pow(sigmaPtp, 2)) /
    (pow(Df, 4) * pow(Pb, 4) * 
     pow((cmm * cpm * (cpp * nmp + cmp * npp) * Ptm + 
          cmp * cpp * (cpm * nmm + cmm * npm) * Ptp), 4)));
    default:
      std::cout << "Invalid asymmetry_index!" << std::endl;
      return 0;
  }
}