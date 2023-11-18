#include "KinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>

using std::string;

// KinematicCuts::KinematicCuts(TTreeReader& reader)
//     : p1_p(reader, "p1_p"), p1_phi(reader, "p1_phi"), p1_theta(reader, "p1_theta"), 
//       p2_p(reader, "p2_p"), p2_phi(reader, "p2_phi"), p2_theta(reader, "p2_theta"), 
//       Q2(reader, "Q2"), W(reader, "W"), Mx(reader, "Mx"), Mx1(reader, "Mx1"),
//       Mx23(reader, "Mx23"), Mh23(reader, "Mh23"), 
//       x(reader, "x"), y(reader, "y"), z23(reader, "z23"), target_pol(reader, "target_pol") {}

KinematicCuts::KinematicCuts(TTreeReader& reader)
    : p1_p(reader, "p1_p"), p1_phi(reader, "p1_phi"), p1_theta(reader, "p1_theta"), 
      p2_p(reader, "p2_p"), p2_phi(reader, "p2_phi"), p2_theta(reader, "p2_theta"), 
      p3_p(reader, "p3_p"), p3_phi(reader, "p3_phi"), p3_theta(reader, "p3_theta"), 
      Mx1(reader, "Mx1"), Mx23(reader, "Mx23"), Mh23(reader, "Mh23"), z23(reader, "z23"),
      Q2(reader, "Q2"), W(reader, "W"), Mx(reader, "Mx"), 
      x(reader, "x"), y(reader, "y"), z(reader, "z"), pT(reader, "pT"), 
      xF(reader, "xF"), target_pol(reader, "target_pol") {}

// KinematicCuts::KinematicCuts(TTreeReader& reader)
//     : Q2(reader, "Q2"), W(reader, "W"), Mx(reader, "Mx"), 
//       x(reader, "x"), y(reader, "y"), z(reader, "z"), pT(reader, "pT"), 
//       xF(reader, "xF"), target_pol(reader, "target_pol") {}

// Function to convert spherical to Cartesian coordinates
void SphericalToCartesian(double p, double phi, double theta, double &x, double &y, double &z) {
    x = p * sin(theta) * cos(phi);
    y = p * sin(theta) * sin(phi);
    z = p * cos(theta);
}

bool KinematicCuts::applyCuts(int currentFits, bool isMC) {
        bool goodEvent = false;
        string property = binNames[currentFits];

        if (property == "xF") {
            goodEvent = *Q2 > 1 && *W > 2 && *Mx > 1.4 && *y < 0.75;
        }
        if (property == "Mx") {
            goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
        }
        if (property == "Q2bin") {
            goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *x>0.2 && *x<0.3 && *pT>0.25 && 
            *pT<0.35 && *xF<0;
        }
        if (property == "PTTFR" || property ==  "xTFR" || property == "zetaTFR" || 
          property == "Q2TFR" || property ==  "x") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *xF<0;
        }
        if (property == "PTCFR" || property == "xCFR" || property == "zetaCFR" ||
          property == "Q2TFR") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.4 && *y<0.75 && *xF>0;
        } 
        //
        // epi+X
        if (property == "xFpip") { 
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75;
        }
        if (property == "PTTFRpip" || property ==  "xTFRpip" || property == "zTFRpip" || 
          property == "Q2TFRpip" || property ==  "xpip") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *xF<0;
        }
        if (property == "PTCFRpip" || property == "xCFRpip" || property == "zCFRpip" ||
          property == "Q2TFRpip") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *xF>0;
        }
        //
        // epi-X
        if (property == "xFpim") { 
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75;
        }
        if (property == "PTTFRpim" || property ==  "xTFRpim" || property == "zTFRpim" || 
          property == "Q2TFRpim" || property ==  "xpim") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *xF<0;
        }
        if (property == "PTCFRpim" || property == "xCFRpim" || property == "zCFRpim" ||
          property == "Q2TFRpim") {
          goodEvent = *Q2>1 && *W>2 && *Mx>1.5 && *y<0.75 && *xF>0;
        }
        
        // // epi+pi+X, exclusive rho
        // if (property == "exclusiveRhoIntegrated" || property == "exclusiveRhoIntegratedx" ||
        //     property == "exclusiveRhoIntegratedt") {

        //   goodEvent = *Q2>1 && *W>2 && *y<0.75 && fabs(*Mx1-0.775)<0.10 && 
        //     fabs(*Mx23-0.938)<0.10 && fabs(*Mh23-0.775)<0.15 && *z23>0.80;
        // }
        // if (property == "exclusiveRhoTransversex" || property == "exclusiveRhoTransverset") {
        //   // Convert spherical coordinates to Cartesian for p1
        //   double p2_x, p2_y, p2_z;
        //   SphericalToCartesian(*p2_p, *p2_phi, *p2_theta, p2_x, p2_y, p2_z);

        //   // Convert spherical coordinates to Cartesian for p2
        //   double p3_x, p3_y, p3_z;
        //   SphericalToCartesian(*p3_p, *p3_phi, *p3_theta, p3_x, p3_y, p3_z);

        //   // Calculate the difference in components
        //   double delta_x = p2_x - p3_x;
        //   double delta_y = p2_y - p3_y;
        //   double delta_z = p2_z - p3_z;

        //   // Calculate the magnitude of the vector difference
        //   double magnitude = sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);

        //   goodEvent = *Q2>1 && *W>2 && *y<0.75 && fabs(*Mx1-0.775)<0.10 && 
        //     fabs(*Mx23-0.938)<0.10 && fabs(*Mh23-0.775)<0.15 && *z23>0.80 && magnitude<1.5;
        // }
        // if (property == "exclusiveRhoLongitudinalx" || property == "exclusiveRhoLongitudinalt") {
        //   // Convert spherical coordinates to Cartesian for p1
        //   double p2_x, p2_y, p2_z;
        //   SphericalToCartesian(*p2_p, *p2_phi, *p2_theta, p2_x, p2_y, p2_z);

        //   // Convert spherical coordinates to Cartesian for p2
        //   double p3_x, p3_y, p3_z;
        //   SphericalToCartesian(*p3_p, *p3_phi, *p3_theta, p3_x, p3_y, p3_z);

        //   // Calculate the difference in components
        //   double delta_x = p2_x - p3_x;
        //   double delta_y = p2_y - p3_y;
        //   double delta_z = p2_z - p3_z;

        //   // Calculate the magnitude of the vector difference
        //   double magnitude = sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);

        //   goodEvent = *Q2>1 && *W>2 && *y<0.75 && fabs(*Mx1-0.775)<0.10 && 
        //     fabs(*Mx23-0.938)<0.10 && fabs(*Mh23-0.775)<0.15 && *z23>0.80 && magnitude>1.5;
        // }

        if (isMC) {
            return goodEvent;
        } else {
            return goodEvent && *target_pol != 0;
        }
    }