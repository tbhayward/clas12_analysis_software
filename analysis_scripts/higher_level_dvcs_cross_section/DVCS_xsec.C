#include <cmath>
#include <iostream>
#include "TMath.h"

// -------------------------------------------------------------------------------------------------
//  Global constants and forward declarations
// -------------------------------------------------------------------------------------------------

// Fundamental constants
double PI     = TMath::Pi();             // π
double alpha  = 1.0/137.036;             // fine-structure constant alpha
double alpha3 = TMath::Power(alpha, 3);  // alpha³
double hbarc2 = 0.38938;                 // (ℏc)² in GeV²*mbarn units

// Particle masses and magnetic moment
double m   = 0.000511;   // electron mass (GeV)
double M   = 0.93827;    // proton mass   (GeV)
double muP = 2.79285;    // proton magnetic moment

// Forward declarations for nucleon form factors (Dirac/Pauli & Sachs)
double GetF1(double T);
double GetF2(double T);
double GetGMP(double tau);
double GetGEP(double tau);

// Forward declarations for Compton Form Factors (CFFs)
// Imaginary parts
double GetImH(double xi, double t);
double GetImHt(double xi, double t);
double GetImE(double xi, double t);
double GetImEt(double xi, double t);
// Real parts
double GetReH(double xi, double t);
double GetReHt(double xi, double t);
double GetReE(double xi, double t);
double GetReEt(double xi, double t);

// Flags to switch on/off particular GPD contributions
bool hasH  = false;  // include H ?
bool hasHt = false;  // include Htilde ?
bool hasE  = false;  // include E ?
bool hasEt = false;  // include Etilde ?

// -------------------------------------------------------------------------------------------------
//  BMK_DVCS class: encapsulates DVCS/BH kinematics & cross sections
// -------------------------------------------------------------------------------------------------
class BMK_DVCS {
public:
    // Primary inputs
    double q_beam;     // beam charge  (+1 or -1)
    double L_beam;     // beam helicity (+1, 0, -1)
    double L_target;   // target helicity
    double EB;         // beam energy (GeV)
    double xB;         // Bjorken x
    double Q2;         // virtuality (GeV²)
    double t;          // momentum transfer (negative)
    double phi;        // azimuthal angle (degrees)
    double theta_Tpol, // target polarization polar angle (rad)
           phi_Tpol;   // target polarization azimuthal angle (rad)

    // Secondary, derived variables
    double xi;     // skewness
    double nu;     // nu = Q² / (2 M xB)
    double y;      // inelasticity
    double eps;    // epsilon parameter (depolarization, ratio of long/tan virtual photon)
    double eps2;   // epsilon squared
    double phi_BMK;// BMK convention: pi(1−phi/180°)
    double t_min;  // minimal kinematically allowed t
    double K2, K;  // BH kinematic factor
    double J;      // Jacobian-like combination
    double Ktild2, Ktilda; // another BH factor

    double Jacob; // Jacobian from (xB,y)->(xB,Q2)

    // Electromagnetic form factors
    double F1, F2;
    double FF_comb1, FF_comb2, FF_comb3; 

    // Compton form factors (CFFs)
    double ImH,  ImHt,  ImE, ImEt;
    double ReH,  ReHt,  ReE, ReEt;

    bool VERB;   // verbosity flag

    // Constructor: set primaries + compute secondaries
    BMK_DVCS(double rq_beam, double rL_beam, double rL_target, double rEB, double rxB,      
            double rQ2, double rt, double rphi, 
            double rtheta_Tpol = 0, double rphi_Tpol = 0);

    // Re-set primaries after construction, then re-compute all derived
    void setPrimaryVars(double rq_beam,  double rL_beam,  double rL_target,
                        double rEB,      double rxB,     double rQ2,
                        double rt,       double rphi,
                        double rtheta_Tpol = 0, double rphi_Tpol = 0);

    // Compute all secondary variables & form factors/CFFs
    void setSecondaryVars();

    // Cross sections & asymmetries
    double CrossSection();      // unpolarized 
    double TPolCrossSection();  // transversely polarized

    // Beam-spin asymmetry (BSA) variants
    double BSA();   // q_beam = -1, electron beam
    double pBSA();  // q_beam = +1, positron beam

    // Target-spin asymmetries
    double TLSA();    // Longitudinal target spin asymmetry
    double TLLSA();   // Double longitudinal (beam+target)
    double TTSAx();   // Transverse target, x-direction
    double TTSAy();   // Transverse target, y-direction
    double TTSSAx();  // twist-three combos
    double TTSSAy();  // twist-three combos

    // Charge asymmetries
    double BCA();   // Beam charge asymmetry
    double BCSA();  // Beam & spin combined
    double BC0SA(); // special combination

    // More exotic asymmetries
    double BCLA();   double BCLLA();
    double BCTxA();  double BCTyA();

    // Underlying squared amplitudes
    double T2();       // BH² + DVCS² − q_beam * 2 Re(BH·DVCS*)
    double BH2();      // Bethe–Heitler² term
    double DVCS2();    // DVCS² term
    double BHDVCS();   // interference term (drives ALU DVCS asymmetry)

    // Harmonically decomposed coefficients for BH, I, DVCS
    double c0_BH();  double c1_BH();  double c2_BH();
    double c0_BH_LP();  double c1_BH_LP();
    double c0_BH_TP();  double c1_BH_TP();  double s1_BH_TP();
    double BHP1();      double BHP2();

    double c0_I();   double c1_I();  double s1_I();
    double c0_I_LP(); double c1_I_LP(); double s1_I_LP();
    double c0_I_TP(); double c1_I_TP(); double s1_I_TP();

    double c0_DVCS();     double c0_DVCS_LP();   double c0_DVCS_TP();
};

// -------------------------------------------------------------------------------------------------
//  Implementation of BMK_DVCS methods
// -------------------------------------------------------------------------------------------------

// Constructor: assign inputs, convert angles, then compute secondaries
BMK_DVCS::BMK_DVCS(double rq_beam, double rL_beam, double rL_target,
                   double rEB,    double rxB,     double rQ2,
                   double rt,     double rphi,
                   double rtheta_Tpol, double rphi_Tpol)
  : q_beam(rq_beam),
    L_beam(rL_beam),
    L_target(rL_target),
    EB(rEB),
    xB(rxB),
    Q2(rQ2),
    t(-TMath::Abs(rt)),                  // ensure t ≤ 0
    phi(rphi),
    // convert degree inputs->radians
    theta_Tpol(PI * rtheta_Tpol/180.0),
    phi_Tpol( PI * rphi_Tpol/180.0),
    VERB(false)
{
    setSecondaryVars();  // build all derived kinematics & form factors
}

// Re-assign primaries after construction
void BMK_DVCS::setPrimaryVars(double rq_beam, double rL_beam, double rL_target,
                              double rEB,     double rxB,     double rQ2,
                              double rt,      double rphi,
                              double rtheta_Tpol, double rphi_Tpol)
{
    q_beam      = rq_beam;
    L_beam      = rL_beam;
    L_target    = rL_target;
    EB          = rEB;
    xB          = rxB;
    Q2          = rQ2;
    t           = -TMath::Abs(rt);
    phi         = rphi;
    theta_Tpol  = PI * rtheta_Tpol/180.0;
    phi_Tpol    = PI * rphi_Tpol/180.0;

    setSecondaryVars();  // recompute everything
}

// Build all secondary variables, form factors, and CFFs
void BMK_DVCS::setSecondaryVars()
{
    // --- Skewness xi = xB (1 + t/(2Q²)) / (2 − xB + xB t/Q²)
    xi = xB * (1 + 0.5*t/Q2) / (2 - xB + xB*t/Q2);

    // BMK azimuth: phi_BMK = pi (1 − phi_trento/180)
    phi_BMK = PI*(1 - phi/180.0);

    // Energy transfer nu = Q² / (2 M xB)
    nu = Q2 / (2.0 * M * xB);

    // Inelasticity y = nu / EB
    y  = nu / EB;

    // Jacobian factor Dy/dQ² = y/Q²
    Jacob = y/Q2;

    // epsilon = 2 xB M / sqrt(Q²)
    eps  = 2 * xB * M / TMath::Sqrt(Q2);
    eps2 = eps*eps;

    // t_min: minimal kinematic t
    t_min = -Q2 * ( 2*(1 - xB)*(1 - TMath::Sqrt(1 + eps2)) + eps2 )
            / (4*xB*(1 - xB) + eps2);

    // Compute BH kinematic factors K², J, etc.
    if(t_min < t) {
        K2 = 0;
    } else {
        // see Eq. (2.10) in 1005.5209
        K2 = -(t - t_min)/Q2 * (1 - xB) * (1 - y - 0.25*y*y*eps2)
             * ( TMath::Sqrt(1 + eps2)
                 + (4*xB*(1 - xB) + eps2)/(4*(1 - xB)) * (t - t_min)/Q2 );
    }
    K = TMath::Sqrt(K2);

    // J factor in BH denominator
    J = (1 - y - 0.5*y*eps2)*(1 + t/Q2) - (1 - xB)*(2 - y)*t/Q2;

    // Another BH factor Ktilde²
    if(t_min < t) {
        Ktild2 = 0;
    } else {
        Ktild2 = (t_min - t) * (
                  (1 - xB)*(1 + eps2)
                  + (t_min - t)*(eps2 + 4*(1 - xB)*xB)/(4*Q2)
                 );
    }
    Ktilda = TMath::Sqrt(Ktild2);

    // --- Electromagnetic form factors
    F1 = GetF1(t);
    F2 = GetF2(t);
    // Useful combinations in BH numerator
    FF_comb1 = F1*F1 - t*F2*F2/(4*M*M);
    FF_comb2 = TMath::Power(F1 + F2, 2);
    FF_comb3 = F1 + t*F2/(4*M*M);

    // --- Compton form factors (CFFs)
    ImH  = GetImH(xi, t);
    ImHt = GetImHt(xi, t);
    ImE  = GetImE(xi, t);
    ImEt = GetImEt(xi, t);

    ReH  = GetReH(xi, t);
    ReHt = GetReHt(xi, t);
    ReE  = GetReE(xi, t);
    ReEt = GetReEt(xi, t);

    // Verbose debug printout
    if(VERB){
        std::cout
            << "Primary kine:  EB=" << EB
            << ", xB=" << xB
            << ", Q2=" << Q2
            << ", t=" << t
            << ", phi=" << phi << "\n"
            << " Derived: xi=" << xi
            << ", y=" << y
            << ", epsilon squared=" << eps2
            << ", F1=" << F1
            << ", F2=" << F2 << "\n"
            << " t_min=" << t_min
            << ", K²=" << K2
            << ", J=" << J << "\n";
    }
}

// -------------------------------------------------------------------------------------------------
//   Cross sections and asymmetries
// -------------------------------------------------------------------------------------------------

// Unpolarized differential cross section (nb/GeV⁴)
double BMK_DVCS::CrossSection()
{
    if(VERB){
        std::cout << "Prefactor ="
                  << (1e9 * hbarc2 * 2*PI * Jacob * alpha3 * xB * y)
                  / (16*PI*PI * Q2 * TMath::Sqrt(1 + eps2)) << std::endl;
    }
    // Master formula × T2 amplitude
    return 1e6 * hbarc2 * alpha3 * xB * y / (16*PI*PI * Q2 * TMath::Sqrt(1 + eps2)) * T2();
}

// Transversely polarized target cross section (nb/GeV⁴ dphi_Tpol)
double BMK_DVCS::TPolCrossSection()
{
    // see Eq. (1) in 1212.6674
    return 1e6 * hbarc2 * alpha3 * xB * y*y / (16*PI*PI * Q2*Q2 * TMath::Sqrt(1 + eps2)) * T2();
}

// Beam spin asymmetry BSA = (sigma(+) − sigma(−)) / (sigma(+) + sigma(−))
// Here we flip L_beam and q_beam appropriately
double BMK_DVCS::BSA()
{
    q_beam  = -1;    // electron beam
    L_target = 0;    // unpolarized
    L_beam   = 1;    // positive helicity
    double X1 = CrossSection();
    L_beam   = -1;   // negative helicity
    double X2 = CrossSection();
    return (X1 - X2)/(X1 + X2);
}

// Same but for positron beam
double BMK_DVCS::pBSA()
{
    q_beam   = +1;  
    L_target = 0;
    L_beam   = 1;
    double X1 = CrossSection();
    L_beam   = -1;
    double X2 = CrossSection();
    return (X1 - X2)/(X1 + X2);
}

// Longitudinal target-spin asymmetry (beam unpolarized, target helicity flip)
double BMK_DVCS::TLSA()
{
    q_beam  = -1;
    L_beam  = 0;
    L_target = 1;    // positive target helicity
    theta_Tpol = 0;  // polarization along z
    double X1 = CrossSection();
    L_target = -1;   // flip target
    double X2 = CrossSection();
    return (X1 - X2)/(X1 + X2);
}

// Double-spin LL asymmetry
double BMK_DVCS::TLLSA()
{
    q_beam    = -1;
    L_beam    = +1;
    L_target  = +1; theta_Tpol = 0;
    double X1 = CrossSection();
    L_beam    = -1; L_target = -1;
    double X2 = CrossSection();
    L_beam    = -1; L_target = +1;
    double X3 = CrossSection();
    L_beam    = +1; L_target = -1;
    double X4 = CrossSection();
    return (X1 + X2 - X3 - X4)/(X1 + X2 + X3 + X4);
}

// Transverse target spin asymmetry, x-component
double BMK_DVCS::TTSAx()
{
    q_beam    = -1;
    L_beam    = 0;
    L_target  = +1;
    theta_Tpol = PI/2; // in-plane
    phi_Tpol   = 0;
    double X1 = CrossSection();
    L_target = -1;
    double X2 = CrossSection();
    return (X1 - X2)/(X1 + X2);
}

// Transverse target spin asymmetry, y-component
double BMK_DVCS::TTSAy()
{
    q_beam    = -1;
    L_beam    = 0;
    L_target  = +1;
    theta_Tpol = PI/2; // in-plane
    phi_Tpol   = PI/2; // rotated 90°
    double X1 = CrossSection();
    L_target = -1;
    double X2 = CrossSection();
    return (X1 - X2)/(X1 + X2);
}

//    Double transverse spin asymmetry along x (twist-three):
//    A_x^{TTSSA} = (sigma(↑↑) + sigma(↓↓) − sigma(↑↓) − sigma(↓↑)) / 
//          (sigma(↑↑) + sigma(↑↓) + sigma(↓↑) + sigma(↓↓))
double BMK_DVCS::TTSSAx() {
    // set beam charge, beam helicity, target helicity, and transverse-polar angle phi_Tpol=0
    q_beam    = -1;          // electron beam
    L_beam    = +1;          // beam helicity +1
    L_target  = +1;          // target helicity +1
    theta_Tpol = PI/2;       // polarization in x–y plane
    phi_Tpol   = 0;          // along x

    // sigma(↑↑)
    double xsec1 = CrossSection();

    // flip target helicity: sigma(↑↓)
    L_target = -1;
    double xsec2 = CrossSection();

    // flip beam helicity: sigma(↓↑)
    L_beam   = -1;
    L_target = +1;
    double xsec3 = CrossSection();

    // flip only target again: sigma(↓↓)
    L_target = -1;
    double xsec4 = CrossSection();

    // build TTSSA_x
    return (xsec1 + xsec4 - xsec2 - xsec3)
         / (xsec1 + xsec2 + xsec3 + xsec4);
}

//    Double transverse spin asymmetry along y (twist-three):
//    same structure but phi_Tpol = 90° -> y-direction
double BMK_DVCS::TTSSAy() {
    q_beam    = -1;
    L_beam    = +1;
    L_target  = +1;
    theta_Tpol = PI/2;       // still transverse
    phi_Tpol   = PI/2;       // rotated into y

    double xsec1 = CrossSection();  // sigma(↑↑)
    L_target = -1;
    double xsec2 = CrossSection();  // sigma(↑↓)

    L_beam   = -1;
    L_target = +1;
    double xsec3 = CrossSection();  // sigma(↓↑)
    L_target = -1;
    double xsec4 = CrossSection();  // sigma(↓↓)

    return (xsec1 + xsec4 - xsec2 - xsec3)
         / (xsec1 + xsec2 + xsec3 + xsec4);
}

//    Beam charge asymmetry (unpolarized):
//    BCA = (sigma(e⁺) − sigma(e⁻)) / (sigma(e⁺) + sigma(e⁻))
double BMK_DVCS::BCA() {
    L_beam   = 0;   // no beam helicity
    L_target = 0;   // unpolarized target
    q_beam   = +1;  // positron
    double xsec1 = CrossSection();
    q_beam   = -1;  // electron
    double xsec2 = CrossSection();
    return (xsec1 - xsec2) / (xsec1 + xsec2);
}

//    Beam & spin combined asymmetry (BCSA):
//    (sigma(+,+) − sigma(−,+) − sigma(+,−) + sigma(−,−)) / total
double BMK_DVCS::BCSA() {
    L_target = 0;    // unpolarized target
    q_beam   = +1;
    L_beam   = +1;
    double xsec1 = CrossSection();  // (+,+)
    
    L_beam   = -1;
    double xsec2 = CrossSection();  // (−,+)
    
    q_beam   = -1;
    L_beam   = +1;
    double xsec3 = CrossSection();  // (+,−)
    
    L_beam   = -1;
    double xsec4 = CrossSection();  // (−,−)
    
    return ( xsec1 - xsec2 - xsec3 + xsec4 )
         / ( xsec1 + xsec2 + xsec3 + xsec4 );
}

//    Zero-helicity beam charge & spin asymmetry (BC0SA):
//    (sigma(+,+) − sigma(−,+) + sigma(+,−) − sigma(−,−)) / total
double BMK_DVCS::BC0SA() {
    L_target = 0;
    q_beam   = +1;
    L_beam   = +1;
    double xsec1 = CrossSection();  // (+,+)
    
    L_beam   = -1;
    double xsec2 = CrossSection();  // (−,+)
    
    q_beam   = -1;
    L_beam   = +1;
    double xsec3 = CrossSection();  // (+,−)
    
    L_beam   = -1;
    double xsec4 = CrossSection();  // (−,−)
    
    return ( xsec1 - xsec2 + xsec3 - xsec4 )
         / ( xsec1 + xsec2 + xsec3 + xsec4 );
}

//    Beam–charge longitudinal asymmetry (BCLA):
//    similar pattern but with L_target != 0, L_beam = 0
double BMK_DVCS::BCLA() {
    q_beam   = +1;
    L_beam   = 0;    // no beam helicity
    L_target = +1;   // polarized target
    double xsec1 = CrossSection();  // e⁺, ↑
    
    L_target = -1;
    double xsec2 = CrossSection();  // e⁺, ↓
    
    q_beam   = -1;
    L_target = +1;
    double xsec3 = CrossSection();  // e⁻, ↑
    
    L_target = -1;
    double xsec4 = CrossSection();  // e⁻, ↓
    
    return ( xsec1 - xsec2 - xsec3 + xsec4 )
         / ( xsec1 + xsec2 + xsec3 + xsec4 );
}

//    Double spin charge–longitudinal asymmetry (BCLLA):
//    eight combinations over q_beam, L_beam, L_target
double BMK_DVCS::BCLLA() {
    // first block: L_beam=+1
    q_beam   = +1;
    L_beam   = +1;
    L_target = +1;
    double xsec1 = CrossSection();  // +,+,↑
    
    L_target = -1;
    double xsec2 = CrossSection();  // +,+,↓
    
    q_beam   = -1;
    L_target = +1;
    double xsec3 = CrossSection();  // −,+,↑
    
    L_target = -1;
    double xsec4 = CrossSection();  // −,+,↓
    
    // second block: L_beam = −1
    q_beam   = +1;
    L_beam   = -1;
    L_target = +1;
    double xsec5 = CrossSection();  // +,−,↑
    
    L_target = -1;
    double xsec6 = CrossSection();  // +,−,↓
    
    q_beam   = -1;
    L_target = +1;
    double xsec7 = CrossSection();  // −,−,↑
    
    L_target = -1;
    double xsec8 = CrossSection();  // −,−,↓
    
    return ( xsec1 - xsec2 + xsec3 - xsec4
           - xsec5 + xsec6 - xsec7 + xsec8 )
         / ( xsec1 + xsec2 + xsec3 + xsec4
           + xsec5 + xsec6 + xsec7 + xsec8 );
}

// Beam–charge transverse–x asymmetry (BCTxA):
double BMK_DVCS::BCTxA() {
    q_beam    = +1;
    L_beam    = 0;
    theta_Tpol = PI/2;
    phi_Tpol   = 0;
    double xsec1 = CrossSection();  // e⁺, ↑_x
    
    L_target = -1;
    double xsec2 = CrossSection();  // e⁺, ↓_x
    
    q_beam   = -1;
    L_target = +1;
    double xsec3 = CrossSection();  // e⁻, ↑_x
    
    L_target = -1;
    double xsec4 = CrossSection();  // e⁻, ↓_x
    
    return ( xsec1 + xsec4 - xsec2 - xsec3 )
         / ( xsec1 + xsec2 + xsec3 + xsec4 );
}

// Beam–charge transverse–y asymmetry (BCTyA):
double BMK_DVCS::BCTyA() {
    q_beam    = +1;
    L_beam    = 0;
    theta_Tpol = PI/2;
    phi_Tpol   = PI/2;
    double xsec1 = CrossSection();  // e⁺, ↑_y
    
    L_target = -1;
    double xsec2 = CrossSection();  // e⁺, ↓_y
    
    q_beam   = -1;
    L_target = +1;
    double xsec3 = CrossSection();  // e⁻, ↑_y
    
    L_target = -1;
    double xsec4 = CrossSection();  // e⁻, ↓_y
    
    return ( xsec1 + xsec4 - xsec2 - xsec3 )
         / ( xsec1 + xsec2 + xsec3 + xsec4 );
}
// T2 = |BH|² + |DVCS|² − q_beam·2Re(BH·DVCS*)
double BMK_DVCS::T2()
{
    return BH2() + DVCS2() - q_beam * BHDVCS();
}

// BH² term: harmonic decomposition in phi_BMK
double BMK_DVCS::BH2()
{
    // Denominator ~ xB² y² t BHP1 BHP2 (1+epsilon²)²
    double DENOM = xB*xB * y*y * t * BHP1() * BHP2() * TMath::Power(1 + eps2, 2);

    // Build c0, c1, c2 + sin harmonic if target polarized
    double res_c0 = c0_BH();
    if(L_target != 0) {
        // add longitudinal + transverse contributions
        res_c0 += L_target * (
                   TMath::Cos(theta_Tpol)*c0_BH_LP()
                 + TMath::Sin(theta_Tpol)*c0_BH_TP()
                 );
    }
    double res_c1 = c1_BH();
    if(L_target != 0) {
        res_c1 += L_target * (
                   TMath::Cos(theta_Tpol)*c1_BH_LP()
                 + TMath::Sin(theta_Tpol)*c1_BH_TP()
                 );
    }
    double res_c2 = c2_BH();

    // sin phi term only for transverse target
    double res_s1 = 0;
    if(L_target != 0) {
        res_s1 = L_target * TMath::Sin(theta_Tpol) * s1_BH_TP();
    }

    // Combine harmonics
    double H = res_c0
             + res_c1 * TMath::Cos(phi_BMK)
             + res_c2 * TMath::Cos(2*phi_BMK)
             + res_s1 * TMath::Sin(phi_BMK);

    if(VERB) {
        std::cout << "BH harmonics: c0=" << res_c0
                  << " c1=" << res_c1
                  << " c2=" << res_c2
                  << " s1=" << res_s1
                  << " -> H=" << H
                  << " / DEN=" << DENOM
                  << std::endl;
    }
    return H / DENOM;
}

// Pure DVCS² term
double BMK_DVCS::DVCS2()
{
    double DENOM = y*y * Q2;
    double res_c0 = c0_DVCS();
    if(L_target != 0) {
        res_c0 += L_target * (
                  TMath::Cos(theta_Tpol) * c0_DVCS_LP()
                + TMath::Sin(theta_Tpol) * c0_DVCS_TP()
                );
    }
    return res_c0 / DENOM;
}

// BH–DVCS interference (I) term
double BMK_DVCS::BHDVCS()
{
    double DENOM = xB * y*y*y * t * BHP1() * BHP2();

    // c0, c1, s1 from interference
    double r0 = c0_I();
    if(L_target != 0) {
        r0 += L_target * (
              TMath::Cos(theta_Tpol) * c0_I_LP()
            + TMath::Sin(theta_Tpol) * c0_I_TP()
            );
    }
    double r1 = c1_I();
    if(L_target != 0) {
        r1 += L_target * (
              TMath::Cos(theta_Tpol) * c1_I_LP()
            + TMath::Sin(theta_Tpol) * c1_I_TP()
            );
    }
    double s1 = s1_I();
    if(L_target != 0) {
        s1 += L_target * (
              TMath::Cos(theta_Tpol) * s1_I_LP()
            + TMath::Sin(theta_Tpol) * s1_I_TP()
            );
    }

    double H = r0
             + r1 * TMath::Cos(phi_BMK)
             + s1 * TMath::Sin(phi_BMK);
    return H / DENOM;
}

// -------------------------------------------------------------------------------------------------
//   Bethe–Heitler & interference coefficients (c0_BH, c1_BH, … c1_I_TP, … c0_DVCS, … )
//   Each one corresponds to a specific equation in the BMK formalism.
// -------------------------------------------------------------------------------------------------

double BMK_DVCS::c0_BH(void){
    // (35) from 0112108
    double c0_BH_unp = 0;
    c0_BH_unp  = (2+eps2)
               * ( 4*xB*xB*M*M*TMath::Power(1+t/Q2,2)/t
                   + 4*(1-xB)*(1+xB*t/Q2) )
               * FF_comb1;
    c0_BH_unp += 4*xB*xB
               * ( xB + (1-xB+0.5*eps2)*TMath::Power(1-t/Q2,2)
                   - xB*(1-2*xB)*t*t/(Q2*Q2) )
               * FF_comb2;
    c0_BH_unp *= TMath::Power(2-y,2.);
    c0_BH_unp += 8*K2
               * ( (2+3*eps2)*Q2*FF_comb1/t + 2*xB*xB*FF_comb2 );
    c0_BH_unp += 8*(1+eps2)*(1-y-eps2*y*y/4)
               * (  2*eps2*(1-t/(4*M*M))*FF_comb1
                   - xB*xB*TMath::Power(1-t/Q2,2)*FF_comb2 );
    return c0_BH_unp;
}

double BMK_DVCS::c0_BH_LP(void){
    // (38) from 0112108
    double c0_BH_LP = 0;
    c0_BH_LP  = (1 - (1-xB)*t/Q2)
               * ( xB*xB*M*M/t * TMath::Power(1+t/Q2,2)
                   + (1-xB)*(1+xB*t/Q2) )
               * FF_comb3;
    c0_BH_LP += 0.5*(0.5*xB*(1-t/Q2)-0.25*t/(M*M))
               * ( 2-xB
                   -2*TMath::Power(1-xB,2)*t/Q2
                   + eps2*(1-t/Q2)
                   - xB*(1-2*xB)*t*t/(Q2*Q2) )
               * (F1+F2);
    c0_BH_LP  = 8 * L_beam * xB * (2-y) * y
              * TMath::Sqrt(1+eps2)/(1-0.25*t/(M*M))
              * (F1+F2) * c0_BH_LP;
    return c0_BH_LP;
}

double BMK_DVCS::c0_BH_TP(void){
    // (40) from 0112108
    double c0_BH_TP = xB*xB*xB *M*M/Q2 * (1 - t/Q2) * (F1+F2)
                   + (1-(1-xB)*t/Q2)
                   * ( xB*xB*M*M/t *(1-t/Q2)*F1 + 0.5*xB*F2 );
    c0_BH_TP = -8 * L_beam * TMath::Cos(phi_Tpol)
             * (2-y) * y * TMath::Sqrt(Q2)/M
             * TMath::Sqrt(1+eps2) * K
             / TMath::Sqrt(1-y-eps2*y*y*0.25)
             * (F1+F2) * c0_BH_TP;
    return c0_BH_TP;
}

double BMK_DVCS::c1_BH(void){
    // (36) from 0112108
    return 8*K*(2.-y)
         * ( (4*xB*xB*M*M/t - 2*xB - eps2)*FF_comb1
             + 2*xB*xB*(1-(1-2*xB)*t/Q2)*FF_comb2 );
}

double BMK_DVCS::c1_BH_LP(void){
    // (39) from 0112108
    double c1_BH_LP = (1+xB
                       - (3-2*xB)*(1+xB*t/Q2)
                       - 4*xB*xB*M*M/t*(1+t*t/(Q2*Q2)))
                    * FF_comb3;
    c1_BH_LP += (0.5*t/(M*M) - xB*(1-t/Q2))
              * (1-xB+xB*t/Q2) * (F1+F2);
    c1_BH_LP  = -8 * L_beam * xB * y * K
               * TMath::Sqrt(1+eps2)/(1-0.25*t/(M*M))
               * (F1+F2) * c1_BH_LP;
    return c1_BH_LP;
}

double BMK_DVCS::c1_BH_TP(void){
    // (41) from 0112108
    double c1_BH_TP = 2*K2*Q2
                    /(t*(1-y-eps2*y*y*0.25))
                    * ( xB*(1-t/Q2)*F1
                        + 0.25*t/(M*M)*F2 )
                    + (1+eps2)*xB*(1-t/Q2)*FF_comb3;
    c1_BH_TP = -16 * L_beam * TMath::Cos(phi_Tpol)
             * xB * y
             * TMath::Sqrt(1-y-eps2*y*y*0.25)
             * M/TMath::Sqrt(Q2)
             * TMath::Sqrt(1+eps2)
             * (F1+F2)
             * c1_BH_TP;
    return c1_BH_TP;
}

double BMK_DVCS::s1_BH_TP(void){
    // (42) from 0112108
    return -16 * L_beam * TMath::Sin(phi_Tpol)
         * y * xB*xB
         * TMath::Sqrt(1-y-0.25*eps2*y*y)
         * M/TMath::Sqrt(Q2)
         * TMath::Power(1+eps2,1.5)
         * (1-t/Q2)
         * (F1+F2)
         * FF_comb3;
}

double BMK_DVCS::c2_BH(void){
    // (37) from 0112108
    return 8*xB*xB*K2
         * (4*M*M*FF_comb1/t + 2*FF_comb2);
}

double BMK_DVCS::BHP1(void){
    // BH propagator factor 1
    return - ( J + 2*TMath::Sqrt(K2) * TMath::Cos(phi_BMK) )
           / ( y*(1.+eps2) );
}

double BMK_DVCS::BHP2(void){
    // BH propagator factor 2
    return 1 + t/Q2 - BHP1();
}

double BMK_DVCS::c0_I(void){
    // (A.1) from 1005.5209 and (B.1) from 1212.6674 (unpolarized interference)
    double C_pp0_unp0 = Ktild2/Q2
                      * TMath::Power(2-y,2) / TMath::Sqrt(1+eps2);
    C_pp0_unp0 += t/Q2 * (1-y-0.25*eps2*y*y) * (2-xB)
                * (1 + (2*xB*(2-xB+0.5*(TMath::Sqrt(1+eps2)-1)+0.5*eps2/xB)*t/Q2+eps2)
                    /((2-xB)*(1+TMath::Sqrt(1+eps2))));
    C_pp0_unp0 *= -4*(2-y)*(1+TMath::Sqrt(1+eps2)) / TMath::Power(1+eps2,2);

    double C_pp0_unpV = (1-y-0.25*y*y*eps2) * 0.5*(1+TMath::Sqrt(1+eps2)) * (1+t/Q2)
                      * (1 + (TMath::Sqrt(1+eps2)-1+2*xB) /(1+TMath::Sqrt(1+eps2))*t/Q2);
    C_pp0_unpV += TMath::Power(2-y,2)*Ktild2 /(TMath::Sqrt(1+eps2)*Q2);
    C_pp0_unpV *= 8*(2-y)*xB*t /(TMath::Power(1+eps2,2)*Q2);

    double C_pp0_unpA = 0.5*(1+TMath::Sqrt(1+eps2)) * (1+TMath::Sqrt(1+eps2)-xB
                         + (TMath::Sqrt(1+eps2)-1 + xB*(3+TMath::Sqrt(1+eps2)-2*xB)
                             /(1+TMath::Sqrt(1+eps2)))*t/Q2) - 2*Ktild2/Q2;
    C_pp0_unpA = 8*(2-y)/TMath::Power(1+eps2,2) * t/Q2 * (TMath::Power(2-y,2)
                   *Ktild2/(TMath::Sqrt(1+eps2)*Q2) *0.5*(1+TMath::Sqrt(1+eps2)-2*xB)
                   + (1-y-0.25*eps2*y*y)*C_pp0_unpA);

    // Compton form factor combinations
    double CFF_unp0 = F1*ReH - t/(4*M*M)*F2*ReE + xB*(F1+F2)*ReHt/(2-xB + xB*t/Q2);
    double CFF_unpV = xB*(F1+F2)*(ReH+ReE) /(2-xB + xB*t/Q2);
    double CFF_unpA = xB*(F1+F2)*ReHt /(2-xB + xB*t/Q2);

    return C_pp0_unp0*CFF_unp0
         + C_pp0_unpV*CFF_unpV
         + C_pp0_unpA*CFF_unpA;
}

double BMK_DVCS::c0_I_LP(void){
    // (A.5) from 1005.5209 (polarized interference)
    double C_pp0_LP0 = TMath::Power(2-y,2)*Ktild2/Q2 + (1-y-0.25*eps2*y*y)
                       *(xB*t/Q2 - (1-t/Q2)*eps2*0.5) *(1 + (TMath::Sqrt(1+eps2)-1+2*xB)
                          /(1+TMath::Sqrt(1+eps2))*t/Q2);
    C_pp0_LP0 = -4*L_beam*y*(1+TMath::Sqrt(1+eps2)) /TMath::Power(1+eps2,2.5)*C_pp0_LP0;

    double C_pp0_LPV = (2-xB+1.5*eps2) *(1 + (4*(1-xB)*xB+eps2)/(4-2*xB+3*eps2)*t/Q2)
                     *(1 + (TMath::Sqrt(1+eps2)-1+2*xB) /(1+TMath::Sqrt(1+eps2))*t/Q2);
    C_pp0_LPV = TMath::Power(2-y,2) *(1+TMath::Sqrt(1+eps2)+2*xB)
               /(1+TMath::Sqrt(1+eps2))*Ktild2/Q2 + (1-y-0.25*eps2*y*y)*C_pp0_LPV;
    C_pp0_LPV = 4*L_beam*y*(1+TMath::Sqrt(1+eps2)) /TMath::Power(1+eps2,2.5)*t/Q2*C_pp0_LPV;

    double C_pp0_LPA = (1+TMath::Sqrt(1+eps2)) *(1-(1-2*xB)*t/Q2)
                     *(1 + (TMath::Sqrt(1+eps2)-1+2*xB) /(1+TMath::Sqrt(1+eps2))*t/Q2);
    C_pp0_LPA = 4*L_beam*y/TMath::Power(1+eps2,2.5)
               *xB*t/Q2 *(2*TMath::Power(2-y,2)*Ktild2/Q2 + (1-y-0.25*eps2*y*y)*C_pp0_LPA);

    double CFF_LP0 = xB/(2-xB+xB*t/Q2)*(F1+F2) *(ReH + 0.5*xB*(1-t/Q2)*ReE
                     -(1-2*xB)*t*ReHt/Q2 - 0.25*t*ReEt/(M*M)) + 2/(2-xB+xB*t/Q2)*F1
                     *((1-xB)*(1+xB*t/Q2) + 0.5*xB + xB*xB*M*M*(3+t/Q2)/Q2)
                     *ReHt + 0.5*t/(M*M)*xB *(0.25*t/(M*M)-0.5*xB*(1-t/Q2)) *ReEt;

    double CFF_LPV = xB/(2-xB+xB*t/Q2)*(F1+F2) *(ReH + 0.5*xB*(1-t/Q2)*ReE);
    double CFF_LPA = xB/(2-xB+xB*t/Q2)*(F1+F2) *(ReHt + 2*xB*M*M*ReHt/Q2 + 0.5*xB*ReEt);

    return C_pp0_LP0*CFF_LP0
         + C_pp0_LPV*CFF_LPV
         + C_pp0_LPA*CFF_LPA;
}

double BMK_DVCS::c0_I_TP(void){
    // (61) from 0112108 (transverse polarized interference)
    double CI_TPM = (xB*xB*F1 - (1-xB)*t*F2/(M*M))/(2-xB) * ImH
                  + (0.25*t/(M*M)
                     *((2-xB)*F1 + xB*xB*F2/(2-xB))
                     + xB*xB*F1/(2-xB)) * ImE
                  - xB*xB/(2-xB)*(F1+F2)
                    *(ImHt + 0.25*t/(M*M)*ImEt);

    double CI_TPP = (F1+F2)
                  *( xB*xB/(2-xB)*(ReH+0.5*xB*ReE)
                    + 0.25*xB*t/(M*M)*ReE )
                  - xB*xB*F1/(2-xB)
                    *(ReHt+0.5*xB*ReEt)
                  + 0.25*t/(M*M)
                    *(4*(1-xB)*F2*ReHt/(2-xB)
                      - (xB*F1 + xB*xB*F2/(2-xB))*ReEt);

    double res = (2-y)*TMath::Sin(phi_Tpol)
               * TMath::Power(2-y,2)/(1-y) * CI_TPM
               - L_beam*y*TMath::Cos(phi_Tpol)
               * (TMath::Power(2-y,2)/(1-y) + 2) * CI_TPP;
    return 8*M*TMath::Sqrt(1-y)*K/TMath::Sqrt(Q2) * res;
}

double BMK_DVCS::c1_I(void){
    // (A.1) from 1005.5209 (unpolarized 1st harmonic interference)
    double C_pp1_unp0 = 1
                      - (1-3*xB)*t/Q2
                      + (1 - TMath::Sqrt(1+eps2) + 3*eps2)
                        /(1+TMath::Sqrt(1+eps2)-eps2)*xB*t/Q2;
    C_pp1_unp0 = -4*K*(2-2*y+y*y+0.5*eps2*y*y)
                *(1+TMath::Sqrt(1+eps2)-eps2)
                /TMath::Power(1+eps2,2.5)*C_pp1_unp0;
    C_pp1_unp0 = -16*K*(1-y-0.25*y*y*eps2)
                /TMath::Power(1+eps2,2.5)
                * ((1 + (1-xB)*0.5*(TMath::Sqrt(1+eps2)-1)/xB)
                   * xB*t/Q2 - 0.75*eps2)
                + C_pp1_unp0;

    double C_pp1_unpV = 16*K*xB*t
                      /(Q2*TMath::Power(1+eps2,2.5))
                      * ( TMath::Power(2-y,2)
                          *(1-(1-2*xB)*t/Q2)
                          + (1-y-0.25*eps2*y*y)
                            *0.5*(1+TMath::Sqrt(1+eps2)-2*xB)
                            *(t-t_min)/Q2 );

    double C_pp1_unpA = -TMath::Power(2-y,2)
                      * (1-0.5*xB + 0.25*(1+TMath::Sqrt(1+eps2)-2*xB)
                         *(1-t/Q2)
                         + (4*xB*(1-xB)+eps2)
                           *0.5/TMath::Sqrt(1+eps2)*(t-t_min)/Q2 );
    C_pp1_unpA = -16*K*t
                /(Q2*TMath::Power(1+eps2,2))
                * ((1-y-0.25*eps2*y*y)
                   * (1-(1-2*xB)*t/Q2
                      + (4*xB*(1-xB)+eps2)
                        *0.25/TMath::Sqrt(1+eps2)*(t-t_min)/Q2)
                   + C_pp1_unpA);

    double CFF_unp0 = F1*ReH
                    - t/(4*M*M)*F2*ReE
                    + xB*(F1+F2)*ReHt/(2-xB + xB*t/Q2);
    double CFF_unpV = xB*(F1+F2)*(ReH+ReE)
                    /(2-xB + xB*t/Q2);
    double CFF_unpA = xB*(F1+F2)*ReHt
                    /(2-xB + xB*t/Q2);

    return C_pp1_unp0*CFF_unp0
         + C_pp1_unpV*CFF_unpV
         + C_pp1_unpA*CFF_unpA;
}

double BMK_DVCS::c1_I_LP(void){
    // (A.5) from 1005.5209 (polarized 1st harmonic interference)
    double C_pp1_LP0 = 1
                      - (1 - 2*xB*(2+TMath::Sqrt(1+eps2))
                        /(1+TMath::Sqrt(1+eps2)-eps2))*t/Q2;
    C_pp1_LP0 = -4*L_beam*K*y*(2-y)
               /TMath::Power(1+eps2,2.5)
               * (1+TMath::Sqrt(1+eps2)-eps2)
               * C_pp1_LP0;

    double C_pp1_LPV = 1
                      - (1 + (1-eps2)/TMath::Sqrt(1+eps2)
                        -2*xB*(1+4*(1-xB)/TMath::Sqrt(1+eps2)))
                        /(2*(TMath::Sqrt(1+eps2)+2*(1-xB)))
                        * (t-t_min)/Q2;
    C_pp1_LPV = 8*L_beam*K*y*(2-y)
               /TMath::Power(1+eps2,2)
               * (TMath::Sqrt(1+eps2)+2*(1-xB))
               * t/Q2 * C_pp1_LPV;

    double C_pp1_LPA = 16*L_beam*K*y*(2-y)
                     /TMath::Power(1+eps2,2.5)
                     * xB*t/Q2 * (1-(1-2*xB)*t/Q2);

    double CFF_LP0 = xB/(2-xB+xB*t/Q2)*(F1+F2)
                   *(ReH + 0.5*xB*(1-t/Q2)*ReE)
                   + (1 + M*M/Q2*xB*xB*(3+t/Q2)/(2-xB+xB*t/Q2))
                     *F1*ReHt
                   - t/Q2*(2*xB*(1-2*xB)/(2-xB+xB*t/Q2))*F2*ReHt
                   - xB/(2-xB+xB*t/Q2)
                     *(0.5*xB*(1-t/Q2)*F1 + 0.25*t*F2/(M*M))
                     *ReEt;

    double CFF_LPV = xB/(2-xB+xB*t/Q2)*(F1+F2)
                   *(ReH + 0.5*xB*(1-t/Q2)*ReE);
    double CFF_LPA = xB/(2-xB+xB*t/Q2)*(F1+F2)
                   *(ReHt + 2*xB*M*M*ReHt/Q2 + 0.5*xB*ReEt);

    return C_pp1_LP0*CFF_LP0
         + C_pp1_LPV*CFF_LPV
         + C_pp1_LPA*CFF_LPA;
}

double BMK_DVCS::c1_I_TP(void){
    // (62) from 0112108 (transverse 1st harmonic interference)
    double CI_TPM = (xB*xB*F1 - (1-xB)*t*F2/(M*M))/(2-xB) * ImH
                  + (0.25*t/(M*M)
                     *((2-xB)*F1 + xB*xB*F2/(2-xB))
                     + xB*xB*F1/(2-xB)) * ImE
                  - xB*xB/(2-xB)*(F1+F2)
                    *(ImHt + 0.25*t/(M*M)*ImEt);

    double CI_TPP = (F1+F2)
                  *(xB*xB/(2-xB)*(ReH+0.5*xB*ReE)
                    + 0.25*xB*t/(M*M)*ReE)
                  - xB*xB*F1/(2-xB)*(ReHt+0.5*xB*ReEt)
                  + 0.25*t/(M*M)
                    *(4*(1-xB)*F2*ReHt/(2-xB)
                      - (xB*F1 + xB*xB*F2/(2-xB))*ReEt);

    double res = TMath::Sin(phi_Tpol)*(2-2*y+y*y)*CI_TPM
               - L_beam*y*(2-y)*TMath::Cos(phi_Tpol)*CI_TPP;
    return 8*M*TMath::Sqrt(1-y)/TMath::Sqrt(Q2) * res;
}

double BMK_DVCS::s1_I(void){
    // (A.1) from 1005.5209 (unpolarized sin phi interference)
    double S_pp1_unp0 = 8*L_beam*K*(2-y)*y/(1+eps2) * (1 + (1-xB + 0.5*(TMath::Sqrt(1+eps2)-1))
                         /(1+eps2)*(t-t_min)/Q2);
    double S_pp1_unpV = -8*L_beam*K*(2-y)*y/TMath::Power(1+eps2,2) * xB*t/Q2
                      * (TMath::Sqrt(1+eps2)-1 + (1+TMath::Sqrt(1+eps2)-2*xB)*t/Q2);
    double S_pp1_unpA = 8*L_beam*K*(2-y)*y/(1+eps2) * t/Q2
                      * (1 - (1-2*xB)*(1+TMath::Sqrt(1+eps2)-2*xB)
                         *0.5/TMath::Sqrt(1+eps2)*(t-t_min)/Q2);

    double CFF_unp0 = F1*ImH - t/(4*M*M)*F2*ImE + xB*(F1+F2)*ImHt/(2-xB+xB*t/Q2);
    double CFF_unpV = xB*(F1+F2)*(ImH+ImE) /(2-xB+xB*t/Q2);
    double CFF_unpA = xB*(F1+F2)*ImHt /(2-xB+xB*t/Q2);

    // (2.28-2.30) from 1005.5209 for the CFF_unp0, CFF_unpV, CFF_unpA terms
    return S_pp1_unp0*CFF_unp0
         + S_pp1_unpV*CFF_unpV
         + S_pp1_unpA*CFF_unpA;
}

double BMK_DVCS::s1_I_LP(void){
    // (A.5) from 1005.5209 (polarized sin phi interference)
    double S_pp1_LP0 = 4*K*(2-2*y+y*y+0.5*eps2*y*y)
                     /TMath::Power(1+eps2,3)
                     *(1+TMath::Sqrt(1+eps2))
                     *(2*TMath::Sqrt(1+eps2)-1
                       +(1+TMath::Sqrt(1+eps2)-2*xB)
                        /(1+TMath::Sqrt(1+eps2))*t/Q2)
                     + 8*K*(2-2*y+0.5*eps2*y*y)
                     /TMath::Power(1+eps2,3)
                     *(3*eps2/2
                       +(1-TMath::Sqrt(1+eps2)-0.5*eps2-xB*(3-TMath::Sqrt(1+eps2)))
                        *t/Q2);
    double S_pp1_LPV = (1 - (1 - TMath::Sqrt(1+eps2)
                           + 0.5*eps2
                           -2*xB*(3*(1-xB)-TMath::Sqrt(1+eps2)))
                           /(4 - xB*(TMath::Sqrt(1+eps2)+3)
                             +2.5*eps2)
                           *t/Q2)
                     *32*K*(1-y+0.25*eps2*y*y)
                     /TMath::Power(1+eps2,3)
                     *(1 - 0.25*(3+TMath::Sqrt(1+eps2))*xB
                       +5*eps2/8)*t/Q2
                     + 8*K*(2-2*y+y*y+0.5*eps2*y*y)
                     /TMath::Power(1+eps2,2)
                     *t/Q2
                     *(1 - (1-2*xB)
                        *(1+TMath::Sqrt(1+eps2)-2*xB)
                        /(2*(1+eps2))*(t-t_min)/Q2);
    double S_pp1_LPA  = -8*K*(2-2*y+y*y+0.5*eps2*y*y)
                      /TMath::Power(1+eps2,3)
                      *xB*t/Q2
                      *(TMath::Sqrt(1+eps2)-1
                        +(1+TMath::Sqrt(1+eps2)-2*xB)
                         *t/Q2)
                      + 8*K*(1-y+0.25*eps2*y*y)
                      /TMath::Power(1+eps2,3)
                      *(3+TMath::Sqrt(1+eps2))
                      *xB*t/Q2
                      *(1 - (3-TMath::Sqrt(1+eps2)-6*xB)
                         /(3+TMath::Sqrt(1+eps2))*t/Q2);

    double CFF_LP0 = xB/(2-xB+xB*t/Q2)*(F1+F2)
                   *(ImH + 0.5*xB*(1-t/Q2)*ImE)
                   + (1 + M*M/Q2*xB*xB*(3+t/Q2)/(2-xB+xB*t/Q2))
                     *F1*ImHt
                   - t/Q2*(2*xB*(1-2*xB)
                            /(2-xB+xB*t/Q2))*F2*ImHt
                   - xB/(2-xB+xB*t/Q2)
                     *(0.5*xB*(1-t/Q2)*F1+0.25*t*F2/(M*M))*ImEt;
    double CFF_LPV = xB/(2-xB+xB*t/Q2)*(F1+F2)
                   *(ImH + 0.5*xB*(1-t/Q2)*ImE);
    double CFF_LPA = xB/(2-xB+xB*t/Q2)*(F1+F2)
                   *(ImHt + 2*xB*M*M*ImHt/Q2 + 0.5*xB*ImEt);

    return S_pp1_LP0*CFF_LP0
         + S_pp1_LPV*CFF_LPV
         + S_pp1_LPA*CFF_LPA;
}

double BMK_DVCS::s1_I_TP(void){
    // (71) from 0112108 (transverse sin phi interference)
    double CI_TPM = (xB*xB*F1 - (1-xB)*t*F2/(M*M))
                  /(2-xB)*ReH
                  + (0.25*t/(M*M)
                     *((2-xB)*F1 + xB*xB*F2/(2-xB))
                     + xB*xB*F1/(2-xB))*ReE
                  - xB*xB/(2-xB)*(F1+F2)
                    *(ReHt + 0.25*t/(M*M)*ReEt);

    double CI_TPP = (F1+F2)
                  *( xB*xB/(2-xB)*(ImH+0.5*xB*ImE)
                    + 0.25*xB*t/(M*M)*ImE )
                  - xB*xB*F1/(2-xB)*(ImHt+0.5*xB*ImEt)
                  + 0.25*t/(M*M)
                    *(4*(1-xB)*F2*ImHt/(2-xB)
                      - (xB*F1 + xB*xB*F2/(2-xB))*ImEt);

    double res = TMath::Cos(phi_Tpol)*(2-2*y+y*y)*CI_TPP
               + L_beam*y*(2-y)*TMath::Sin(phi_Tpol)*CI_TPM;
    return 8*M*TMath::Sqrt(1-y)/TMath::Sqrt(Q2) * res;
}

double BMK_DVCS::c0_DVCS(void){
    // (2.18) from 1005.5209 & (37) from 1212.6674 (pure DVCS)
    double C_0_unp = 2*(2-2*y+y*y+0.5*eps2*y*y)/(1+eps2);

    double C_DVCS_unp = 4*(1-xB)*(1+xB*t/Q2)
                     /TMath::Power(2-xB+xB*t/Q2,2)
                     * (ReH*ReH + ImH*ImH + ReHt*ReHt + ImHt*ImHt)
                     + (2+t/Q2)*eps2
                     /TMath::Power(2-xB+xB*t/Q2,2)
                     * (ReHt*ReHt + ImHt*ImHt)
                     - 0.25*t/(M*M)
                       * (ReE*ReE + ImE*ImE)
                     - xB*xB
                     /TMath::Power(2-xB+xB*t/Q2,2)
                     * ( TMath::Power(1+t/Q2,2)
                         * (2*ReH*ImE + 2*ReE*ImH + ReE*ReE + ImE*ImE)
                         + 2*(ReHt*ImEt + ReEt*ImHt)
                         + 0.25*t/(M*M)*(ReE*ReE + ImE*ImE) );

    return C_0_unp * C_DVCS_unp;
}

double BMK_DVCS::c0_DVCS_LP(void){
    // (2.20) from 1005.5209 & (40) from 1212.6674 (polarized DVCS)
    double C_0_LP = 2 * L_beam * y * (2-y) / TMath::Sqrt(1+eps2);

    double C_DVCS_LP = (4*(1-xB)*(1+xB*t/Q2)
                      + 2*(1-xB+0.5*(1+t/Q2))*eps2)
                     /TMath::Power(2-xB+xB*t/Q2,2)
                     * (2*ReH*ImHt + 2*ReHt*ImH)
                     - (xB*xB*(1+xB*t/Q2-(1-xB)*t/Q2)*t/Q2)
                       /TMath::Power(2-xB+xB*t/Q2,2)
                       * (2*ReH*ImEt + 2*ReEt*ImH + 2*ReHt*ImE + 2*ReE*ImHt)
                     - 0.5*(4*xB*(1-xB)*(1+xB*t/Q2)*t/Q2
                           + xB*TMath::Power(1+t/Q2,2)*eps2)
                       /TMath::Power(2-xB+xB*t/Q2,2)
                       * (2*ReHt*ImE + 2*ReE*ImHt)
                     - xB*(0.5*xB*xB*TMath::Power(1+t/Q2,2)
                            /TMath::Power(2-xB+xB*t/Q2,2)
                            + 0.25*t/(M*M))
                       /TMath::Power(2-xB+xB*t/Q2,2)
                       * (2*ReE*ImEt + 2*ReEt*ImE);

    return C_0_LP * C_DVCS_LP;
}

double BMK_DVCS::c0_DVCS_TP(void){
    // (43), (50), (51) from 1212.6674 (transverse polarized DVCS)
    double C_0_TPP = (2-y)/(1+eps2) * Ktilda/M
                   * L_beam * TMath::Cos(phi_Tpol);

    double C_0_TPM = -(2-y)/(1+eps2) * Ktilda/M
                   * TMath::Sin(phi_Tpol)
                   * (2-2*y+y*y+0.5*eps2*y*y)/(2-y);

    double C_DVCS_TPP = xB*(2*ReH*ImEt+2*ReEt*ImH)
                     + 4*xB*(1-2*xB)*M*M/Q2*(2*ReH*ImHt+2*ReHt*ImH)
                     - (2-xB+xB*t/Q2+0.5*(3+t/Q2)*eps2)
                       *(2*ReHt*ImE+2*ReE*ImHt)
                     + 0.5*xB*xB*(1-t/Q2)
                       *(2*ReE*ImEt+2*ReEt*ImE);
    C_DVCS_TPP *= 2/TMath::Power(2-xB+xB*t/Q2,2);

    double C_DVCS_TPM = (2*(2*ImH*ReE - 2*ReH*ImE)
                       - 2*xB*(2*ImHt*ReEt - 2*ReHt*ImEt))
                       /TMath::Power(2-xB+xB*t/Q2,2);

    return C_0_TPP*C_DVCS_TPP + C_0_TPM*C_DVCS_TPM;
}

// -------------------------------------------------------------------------------------------------
//  Proton Sachs form factors via 12-parameter z-expansion fits
//  —————————————————————————————————————————————————————————————————————————————————————————————
//  GetGMP: magnetic form factor G_M^p(tau) multiplied by proton magnetic moment μ_p
//  GetGEP: electric form factor G_E^p(tau)
//  See: Arrington et al., “Hard two‐photon exchange in elastic e±p scattering,” EPJA 54 (2018).
// -------------------------------------------------------------------------------------------------

double GetGMP(double tau) {
    // tau = Q² / (4 M²)
    double Q2   = 4 * M * M * tau;
    double tcut = 4 * 0.13957 * 0.13957; // four-pion threshold (GeV²)
    double t0   = -0.7;                 // expansion parameter

    // Conformal mapping to z-variable
    double z = ( TMath::Sqrt(tcut + Q2) - TMath::Sqrt(tcut - t0) )
             / ( TMath::Sqrt(tcut + Q2) + TMath::Sqrt(tcut - t0) );

    // 12-coefficient expansion in z
    double res = 0;
    res +=  0.264142994136;
    res += -1.095306122120 * z;
    res +=  1.218553781780 * TMath::Power(z, 2);
    res +=  0.661136493537 * TMath::Power(z, 3);
    res += -1.405678925030 * TMath::Power(z, 4);
    res += -1.356418438880 * TMath::Power(z, 5);
    res +=  1.447029155340 * TMath::Power(z, 6);
    res +=  4.235669735900 * TMath::Power(z, 7);
    res += -5.334045653410 * TMath::Power(z, 8);
    res += -2.916300520960 * TMath::Power(z, 9);
    res +=  8.707403067570 * TMath::Power(z, 10);
    res += -5.706999943750 * TMath::Power(z, 11);
    res +=  1.280814375890 * TMath::Power(z, 12);

    // Return mu_p × expansion result
    return muP * res;
}

double GetGEP(double tau) {
    // tau = Q² / (4 M²)
    double Q2   = 4 * M * M * tau;
    double tcut = 4 * 0.13957 * 0.13957;
    double t0   = -0.7;

    // Same z mapping as for G_M
    double z = ( TMath::Sqrt(tcut + Q2) - TMath::Sqrt(tcut - t0) )
             / ( TMath::Sqrt(tcut + Q2) + TMath::Sqrt(tcut - t0) );

    // 12-coefficient expansion for G_E^p
    double res = 0;
    res +=  0.239163298067;
    res += -1.109858574410 * z;
    res +=  1.444380813060 * TMath::Power(z, 2);
    res +=  0.479569465603 * TMath::Power(z, 3);
    res += -2.286894741870 * TMath::Power(z, 4);
    res +=  1.126632984980 * TMath::Power(z, 5);
    res +=  1.250619843540 * TMath::Power(z, 6);
    res += -3.631020471590 * TMath::Power(z, 7);
    res +=  4.082217023790 * TMath::Power(z, 8);
    res +=  0.504097346499 * TMath::Power(z, 9);
    res += -5.085120460510 * TMath::Power(z, 10);
    res +=  3.967742543950 * TMath::Power(z, 11);
    res += -0.981529071103 * TMath::Power(z, 12);

    return res;
}

// Dirac & Pauli form factors
double GetF1(double T) {
    double t = TMath::Abs(T);
    double tau = t / (4*M*M);
    return (GetGEP(tau) + tau*GetGMP(tau)) / (1 + tau);
}

double GetF2(double T) {
    double t = TMath::Abs(T);
    double tau = t / (4*M*M);
    return (GetGMP(tau) - GetGEP(tau)) / (1 + tau);
}

// -------------------------------------------------------------------------------------------------
//   Compton Form Factor (CFF) models—imaginary parts
// -------------------------------------------------------------------------------------------------
double renormImag = 1.0;
double renormReal = 1.0;

// -----------------------------------------------------------------------------

// GPD‐H defaults (original values, BKM? VGG? not entirely clear)
double r_H      = 0.9;
double alpha0_H = 0.43;
double alpha1_H = 0.85;
double n_H      = 1.35;
double b_H      = 0.4;
double Mm2_H    = 0.64;
double P_H      = 1.0;
double GetImH(double xi, double t) {
    if(!hasH) return 0.0;
    // build H‐ansatz
    double alphaH   = alpha0_H + alpha1_H * t;
    double pref = TMath::Pi()*5.0/9.0 * n_H * r_H / (1 + xi);
    double xfac = TMath::Power(2*xi/(1+xi), -alphaH);
    double yfac = TMath::Power((1 - xi)/(1+xi), b_H);
    double tfac = TMath::Power(1 - ((1 - xi)/(1+xi))*t/Mm2_H, -P_H);
    return renormImag * pref * xfac * yfac * tfac * 2.0; // *2 correction from VGG
}

// GPD‐Htilde defaults
double r_Ht      = 7.0;
double alpha0_Ht = 0.43;
double alpha1_Ht = 0.85;
double n_Ht      = 0.6;
double b_Ht      = 2.0;
double Mm2_Ht    = 0.8;
double P_Ht      = 1.0;
double GetImHt(double xi, double t) {
    if(!hasHt) return 0.0;
    double alphaHt   = alpha0_Ht + alpha1_Ht * t;
    double pref  = TMath::Pi()*5.0/9.0 * n_Ht * r_Ht / (1 + xi);
    double xfac  = TMath::Power(2*xi/(1+xi), -alphaHt);
    double yfac  = TMath::Power((1 - xi)/(1+xi), b_Ht);
    double tfac  = TMath::Power(1 - ((1 - xi)/(1+xi))*t/Mm2_Ht, -P_Ht);
    return renormImag * pref * xfac * yfac * tfac * 0.4; // *0.4 correction from VGG
}

// GPD‐E defaults (same as H for valence)
double r_E      = 0.9;
double alpha0_E = 0.43;
double alpha1_E = 0.85;
double n_E      = 1.35;
double b_E      = 0.4;
double Mm2_E    = 0.64;
double P_E      = 1.0;
double GetImE(double xi, double t) {
    // from HERMES CFF paper 1301.1230 argue ImE = 0.5 x ImH (the factor of 2 is in ImH)
    if(!hasE) return 0.0;
    double alphaE    = alpha0_E + alpha1_E * t;
    double pref  = TMath::Pi()*5.0/9.0 * n_E * r_E / (1 + xi);
    double xfac  = TMath::Power(2*xi/(1+xi), -alphaE);
    double yfac  = TMath::Power((1 - xi)/(1+xi), b_E);
    double tfac  = TMath::Power(1 - ((1 - xi)/(1+xi))*t/Mm2_E, -P_E);
    return renormImag * pref * xfac * yfac * tfac;
}

// GPD‐Etilde defaults
double r_Et      = 1;
double alpha0_Et = 0.0;
double alpha1_Et = 0.0;
double n_Et      = 0.0; // this sets ImEt to 0 in default model
double b_Et      = 0.0;
double Mm2_Et    = 0.0;
double P_Et      = 0.0;
double GetImEt(double xi, double t) {
    if(!hasEt) return 0.0;
    double alphaEt    = alpha0_Et + alpha1_Et * t;
    double pref  = TMath::Pi()*5.0/9.0 * n_Et * r_Et / (1 + xi);
    double xfac  = TMath::Power(2*xi/(1+xi), -alphaEt);
    double yfac  = TMath::Power((1 - xi)/(1+xi), b_Et);
    double tfac  = TMath::Power(1 - ((1 - xi)/(1+xi))*t/Mm2_Et, -P_Et);
    return renormImag * pref * xfac * yfac * tfac;
}

// -------------------------------------------------------------------------------------------------
//   Compton Form Factor (CFF) models—real parts
// -------------------------------------------------------------------------------------------------

double GetReH(double xi, double t) {
    if(!hasH) return 0.0;
    // simple polynomial ansatz
    double res = -12*xi*TMath::Power(1 - xi, 2) * TMath::Sqrt(TMath::Abs(t))
               / TMath::Power(1 - t/0.7, 2);
    res += -3 * TMath::Power(1 - xi, 4) / TMath::Power(1 - t/1.1, 2);
    return renormReal * res / (1 + TMath::Power(t/0.8, 4));
}

double GetReHt(double xi, double t) {
    if(!hasHt) return 0.0;
    double res = -12*xi*TMath::Power(1 - xi, 2) / TMath::Power(1 - t/1.5, 2);
    return renormReal * res;
}

double GetReE(double xi, double t) {
    if(!hasE) return 0.0;
    double res = -7 * xi*TMath::Power(1 - xi, 2) * TMath::Sqrt(TMath::Abs(t))
               / TMath::Power(1 - t/0.7, 2);
    res += -3 * TMath::Power(1 - xi, 2) / TMath::Power(1 - t/1.2, 2);
    return renormReal * res / (1 + TMath::Power(t, 4));
}

double GetReEt(double xi, double t) {
    if(!hasEt) return 0.0;
    // small-t behavior ~1/t
    double res = 10.0/t * 1.0/TMath::Power(1 + TMath::Power(3*xi, 4), 1);
    return renormReal * res;
}