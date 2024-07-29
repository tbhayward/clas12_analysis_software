/**
 *
 * @author Timothy B. Hayward
 */

package analyzers;

import org.jlab.clas.physics.*;

public class kinematic_variables {
    
    double particle_mass (int pid) {
	if (pid==11||pid==-11) { // electron is pid=11, positron is pid=-11
            return 0.0005109989461;
        } else if (pid == 22) {
            return 0; 
        } else if (pid == 111) {
            return 0.1349768;
	} else if (pid==211||pid==-211) { // pions
            return 0.139570;
	} else if (pid==321||pid==-321) { // kaons
            return 0.493677;
	} else if (pid==2212||pid==-2212) { // protons
            return 0.938272;
//            return 1.875;
	} else if (pid==113) { // rho0
            return 0.7754;
        }
            return -1;
    }
    
    /*~~~~~~~~~~~~~~~~~ DIS ~~~~~~~~~~~~~~~~~*/
    
    double Q2(LorentzVector lv_e) {
        return -lv_e.mass2();
    }
    
    double nu(LorentzVector lv_beam, LorentzVector lv_e) {
        return lv_beam.e()-lv_e.e();
    }
    
    double x(double Q2, double nu) { // x-Bjorken
        return Q2 / (2 * particle_mass(2212) * nu);
    }
    
    double W(double Q2, double nu) {
        return Math.pow(Math.pow(particle_mass(2212),2)+2*particle_mass(2212)*nu - Q2, 0.5);
    }
    
    double y(double nu, LorentzVector lv_beam) {
        return nu/lv_beam.e();
    }
    
    double gamma(double Q2, double x) {
        return 2*particle_mass(2212)*x/Math.pow(Q2, 0.5);
    }
    
    double tmin(double x) {
        return -Math.pow((particle_mass(2212)*x),2)/(1-x);
    }
    
    double t(double p, double theta) {
        double mp = particle_mass(2212); // proton mass in GeV
        double E = mp; // target proton energy (written as E to help checking calculation but at rest)
        return 2*mp*(E - p) - 2*Math.sqrt(mp*mp + E*E)*Math.sqrt(mp*mp + p*p) +
            2*Math.sqrt(mp*mp + E*E)*Math.sqrt(mp*mp + p*p)*Math.cos(theta);
    }
 
    double Depolarization_A(double gamma, double y) {
        return 1/(1+gamma*gamma)*(1-y+y*y/2+y*y*gamma*gamma/4);
    }
    
    double Depolarization_B(double gamma, double y) {
        return 1/(1+gamma*gamma)*(1-y-y*y*gamma*gamma/4);
    }
    
    double Depolarization_C(double gamma, double y) {
        return (y/Math.pow(1+gamma*gamma, 0.5))*(1-y/2);
    }
    
    double Depolarization_V(double gamma, double y) {
        return (2-y)/(1+gamma*gamma)*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
    }
    
    double Depolarization_W(double gamma, double y) {
        return y/(Math.pow(1+gamma*gamma, 0.5))*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
    }
        
    double Mx(LorentzVector lv_q, LorentzVector lv_target) {
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target);  
        return lv_Mx.mass();
    }
    
    double Mx2(LorentzVector lv_q, LorentzVector lv_target) {
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target);  
        return lv_Mx.mass2();
    }
    
    double Mx(LorentzVector lv_q, LorentzVector lv_target, LorentzVector lv_p) {
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target); lv_Mx.sub(lv_p);  
        return lv_Mx.mass();
    }
    
    double Mx2(LorentzVector lv_q, LorentzVector lv_target, LorentzVector lv_p) {
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target); lv_Mx.sub(lv_p);
        return lv_Mx.mass2();
    }
    
    /*~~~~~~~~~~~~~~~~~ Single Hadron ~~~~~~~~~~~~~~~~~*/
    double z(LorentzVector lv_p, LorentzVector lv_q) {
        return lv_p.e()/lv_q.e();
    }
    
    /*~~~~~~~~~~~~~~~~~ Lorentz Boosts ~~~~~~~~~~~~~~~~~*/
    
    // boost to gamma*-nucleon center of mass frame
    LorentzVector lv_boost_gN(LorentzVector lv_target, LorentzVector lv_q, LorentzVector lv_p_gN) {
        // set up boost to gamma*-nucleon center of mass frame
        LorentzVector gN = new LorentzVector(lv_q);
	gN.add(lv_target);
	Vector3 gNBoost = gN.boostVector();
	gNBoost.negative();
        lv_p_gN.boost(gNBoost);
        return lv_p_gN;
    }
    
    /*~~~~~~~~~~~~~~~~~ Exclusivity ~~~~~~~~~~~~~~~~~*/
    
    double Emiss0(LorentzVector lv_beam, LorentzVector lv_target, LorentzVector lv_e) {
        return lv_beam.e() + lv_target.e() - lv_e.e();
    }
    
    double Emiss1(LorentzVector lv_beam, LorentzVector lv_target, LorentzVector lv_e,
            LorentzVector lv_p) {
        return lv_beam.e() + lv_target.e() - lv_e.e() - lv_p.e();
    }
    
    double Emiss2(LorentzVector lv_beam, LorentzVector lv_target, LorentzVector lv_e,
            LorentzVector lv_p1, LorentzVector lv_p2) {
        return lv_beam.e() + lv_target.e() - lv_e.e() - lv_p1.e() - lv_p2.e();
    }
    
    double Emiss3(LorentzVector lv_beam, LorentzVector lv_target, LorentzVector lv_e,
            LorentzVector lv_p1, LorentzVector lv_p2, LorentzVector lv_p3) {
        return lv_beam.e() + lv_target.e() - lv_e.e() - lv_p1.e() - lv_p2.e() - lv_p3.e();
    }
    
    double theta_gamma_gamma(LorentzVector lv_beam, LorentzVector lv_target, LorentzVector lv_e,
            LorentzVector lv_p, LorentzVector lv_gamma) {
        
        // Calculate the missing momentum (expected photon momentum)
        LorentzVector lv_expected_gamma = new LorentzVector(lv_beam);
        lv_expected_gamma.sub(lv_e); lv_expected_gamma.sub(lv_p); lv_expected_gamma.add(lv_target);

        // Get the momentum vectors
        Vector3 p_expected_gamma = lv_expected_gamma.vect();
        Vector3 p_detected_gamma = lv_gamma.vect();
        // Calculate the dot product of the two vectors
        double dot_product = p_expected_gamma.dot(p_detected_gamma);
        // Calculate the magnitudes of the two vectors
        double mag_expected_gamma = p_expected_gamma.mag();
        double mag_detected_gamma = p_detected_gamma.mag();
        // Calculate the cosine of the angle between the two vectors
        double cosTheta = dot_product / (mag_expected_gamma * mag_detected_gamma);
        // Convert the cosine to an angle in radians
        double theta = Math.acos(cosTheta);
        // Convert radians to degrees 
        double theta_degrees = Math.toDegrees(theta);

        return theta_degrees;  // Return the angle in degrees
    }
    
    double pTmiss(LorentzVector lv_beam, LorentzVector lv_target, LorentzVector lv_e, 
            LorentzVector lv_p, LorentzVector lv_gamma) {
        // Calculate the initial momentum (beam + target)
        LorentzVector lv_initial = new LorentzVector(lv_beam);
        lv_initial.add(lv_target);
        // Calculate the final momentum (scattered electron + detected proton + detected photon)
        LorentzVector lv_final = new LorentzVector(lv_e);
        lv_final.add(lv_p);
        lv_final.add(lv_gamma);
        // Calculate the missing momentum (initial - final)
        LorentzVector lv_missing = new LorentzVector(lv_initial);
        lv_missing.sub(lv_final);
        // Extract the transverse components of the missing momentum
        double pTmiss_x = lv_missing.px();
        double pTmiss_y = lv_missing.py();
        // Calculate the magnitude of the transverse missing momentum
        double pTmiss = Math.sqrt(pTmiss_x * pTmiss_x + pTmiss_y * pTmiss_y);

        return pTmiss;  // Return the missing transverse momentum
    }
}
