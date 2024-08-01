/**
 *
 * @author Timothy B. Hayward
 */

package extended_kinematic_fitters;

import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.pdg.PDGParticle;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;

import org.jlab.clas.physics.*;

public class parent_hadron_creation {
    
    Particle pi0_check(PhysicsEvent physEvent, int current_p1, int current_p2) {
        PDGDatabase PDGDatabase = new PDGDatabase();
        // Check if the particle with PID -111 exists
        if (!PDGDatabase.isValidId(-111)) {
            // Create a new particle
            PDGParticle particle = new PDGParticle("pi0background", -111, 0.1349766, 0);
            // Add the particle to the database
            PDGDatabase.addParticle(particle);
        }
        
        Particle gamma_1 = physEvent.getParticle("[22,"+current_p1+"]");
        LorentzVector lv_gamma_1 = new LorentzVector();
        lv_gamma_1.setPxPyPzM(gamma_1.px(), gamma_1.py(), gamma_1.pz(), 0);
                    
        Particle gamma_2 = physEvent.getParticle("[22,"+current_p2+"]");
        LorentzVector lv_gamma_2 = new LorentzVector();
        lv_gamma_2.setPxPyPzM(gamma_2.px(), gamma_2.py(), gamma_2.pz(), 0);
                    
        LorentzVector lv_pi0_candidate = new LorentzVector(lv_gamma_1); lv_pi0_candidate.add(lv_gamma_2);
        
        // candidate is in the mass range for pi0
        if (0.11 < lv_pi0_candidate.mass() && lv_pi0_candidate.mass() < 0.16) {
            Particle part = new Particle(111, lv_pi0_candidate.px(), lv_pi0_candidate.py(), 
                lv_pi0_candidate.pz(), 0, 0, 0);
            return part;
        }
        
        // candidate is in the sideband of the pi0 mass spectrum, assigned to pid == -111 (negative of pi0 pid)
        if (0.22 < lv_pi0_candidate.mass() && lv_pi0_candidate.mass() < 0.45) {
            Particle part = new Particle(-111, lv_pi0_candidate.px(), lv_pi0_candidate.py(), 
                lv_pi0_candidate.pz(), 0, 0, 0);
            return part;
        }
        
        return null;
    }
    
    Particle rho0_check(PhysicsEvent physEvent, int current_p1, int current_p2) {
        Particle pi_1 = physEvent.getParticle("[211,"+current_p1+"]");
        LorentzVector lv_pi_1 = new LorentzVector();
        lv_pi_1.setPxPyPzM(pi_1.px(), pi_1.py(), pi_1.pz(), 0);
                    
        Particle pi_2 = physEvent.getParticle("[-211,"+current_p1+"]");
        LorentzVector lv_pi_2 = new LorentzVector();
        lv_pi_2.setPxPyPzM(pi_2.px(), pi_2.py(), pi_2.pz(), 0);
                    
        LorentzVector lv_rho_candidate = new LorentzVector(lv_pi_1); lv_rho_candidate.add(lv_pi_2);
        
        // candidate is in the mass range for pi0
        if (0.65 < lv_rho_candidate.mass() && lv_rho_candidate.mass() < 0.90) {
            // add rho0 to event
            Particle part = new Particle(113, lv_rho_candidate.px(), lv_rho_candidate.py(), 
                lv_rho_candidate.pz(), 0, 0, 0);
            return part;
        }
        
        return null;
    }
    
    Particle rhop_check(PhysicsEvent physEvent, int current_p1, int current_p2) {
        Particle pi_1 = physEvent.getParticle("[211,"+current_p1+"]");
        LorentzVector lv_pi_1 = new LorentzVector();
        lv_pi_1.setPxPyPzM(pi_1.px(), pi_1.py(), pi_1.pz(), 0);
                    
        Particle pi_2 = physEvent.getParticle("[111,"+current_p1+"]");
        LorentzVector lv_pi_2 = new LorentzVector();
        lv_pi_2.setPxPyPzM(pi_2.px(), pi_2.py(), pi_2.pz(), 0);
                    
        LorentzVector lv_rho_candidate = new LorentzVector(lv_pi_1); lv_rho_candidate.add(lv_pi_2);
        
        // candidate is in the mass range for pi0
        if (0.65 < lv_rho_candidate.mass() && lv_rho_candidate.mass() < 0.90) {
            // add rho+ to event
            Particle part = new Particle(213, lv_rho_candidate.px(), lv_rho_candidate.py(), 
                lv_rho_candidate.pz(), 0, 0, 0);
            return part;
        }
        
        return null;
    }
    
    Particle rhom_check(PhysicsEvent physEvent, int current_p1, int current_p2) {
        PDGDatabase PDGDatabase = new PDGDatabase();
        if (!PDGDatabase.isValidId(-213)) {
            // Create a new particle
            PDGParticle particle = new PDGParticle("rho-", -213, 0.7754, 0);
            // Add the particle to the database
            PDGDatabase.addParticle(particle);
        }
        
        Particle pi_1 = physEvent.getParticle("[-211,"+current_p1+"]");
        LorentzVector lv_pi_1 = new LorentzVector();
        lv_pi_1.setPxPyPzM(pi_1.px(), pi_1.py(), pi_1.pz(), 0);
                    
        Particle pi_2 = physEvent.getParticle("[111,"+current_p1+"]");
        LorentzVector lv_pi_2 = new LorentzVector();
        lv_pi_2.setPxPyPzM(pi_2.px(), pi_2.py(), pi_2.pz(), 0);
                    
        LorentzVector lv_rho_candidate = new LorentzVector(lv_pi_1); lv_rho_candidate.add(lv_pi_2);
        
        // candidate is in the mass range for pi0
        if (0.65 < lv_rho_candidate.mass() && lv_rho_candidate.mass() < 0.90) {
            // add rho+ to event
            Particle part = new Particle(-213, lv_rho_candidate.px(), lv_rho_candidate.py(), 
                lv_rho_candidate.pz(), 0, 0, 0);
            return part;
        }
        
        return null;
    }
}
