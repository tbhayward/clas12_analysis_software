/**
 *
 * @author Timothy B. Hayward
 */

package analyzers;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

import org.jlab.clas.physics.*;

public class kinematic_variables {
    
    double particle_mass (int pid) {
	if (pid==11||pid==-11) { // electron is pid=11, positron is pid=-11
            return 0.0005109989461;
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
    
}
