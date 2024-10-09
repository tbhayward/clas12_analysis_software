/**
 *
 * @author Timothy B. Hayward
 */

package extended_kinematic_fitters;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

public class monte_carlo_fitter extends GenericKinematicFitter {

    protected final Double mybeam;

    public monte_carlo_fitter(double beam) {
        super(beam);
        mybeam = beam;
    } 
    /**
     * Returns PhysicsEvent object with reconstructed particles.
     *
     * @param event - DataEvent object
     * @return PhysicsEvent : event containing particles.
     */
    @Override
    public PhysicsEvent getPhysicsEvent(DataEvent event) {
        
        boolean banks_test = true; // check to see if the event has all of the banks present
        if (!(event.hasBank("MC::Particle"))) {
            banks_test = false;
        } 
        if (banks_test) {
            PhysicsEvent physEvent = new PhysicsEvent();
            HipoDataBank mcpbBank = (HipoDataBank) event.getBank("MC::Particle"); // load MC bank
            for (int current_part = 0; current_part < mcpbBank.rows(); current_part++) {
                int pid = mcpbBank.getInt("pid", current_part);
                if (pid!=0) { // do not load tracks that the EventBuilder did not assign a pid to 
                    float vx = mcpbBank.getFloat("vx",current_part);
                    float vy = mcpbBank.getFloat("vy",current_part);
                    float vz = mcpbBank.getFloat("vz",current_part);
                    float px = mcpbBank.getFloat("px", current_part);
                    float py = mcpbBank.getFloat("py", current_part);
                    float pz = mcpbBank.getFloat("pz", current_part);
                    Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
                    physEvent.addParticle(part);
                }
            }
            
            int num_gamma = physEvent.countByPid(22); // number of photons in event

            parent_hadron_creation parent_hadron_creation = new parent_hadron_creation();

            for (int current_p1 = 0; current_p1 < num_gamma; current_p1++) {
                for (int current_p2 = 0; current_p2 < num_gamma; current_p2++) {
                    if (current_p1 == current_p2) {
                        continue;
                    }
                    Particle part = parent_hadron_creation.pi0_check(physEvent, current_p1, current_p2);

                    // Check if a valid Particle was returned before adding it to the event
                    if (part != null) {
                        physEvent.addParticle(part);
                    }
                }
            }
            
            int num_pip = physEvent.countByPid(211); // number of pi+ in event
            int num_pim = physEvent.countByPid(-211); // number of pi- in event
            int num_pi0 = physEvent.countByPid(111); // number of pi0 in event

            for (int current_p1 = 0; current_p1 < num_pip; current_p1++) {
                for (int current_p2 = 0; current_p2 < num_pim; current_p2++) {
                    if (current_p1 == current_p2) {
                        continue;
                    }
                    Particle part = parent_hadron_creation.rho0_check(physEvent, current_p1, current_p2);

                    // Check if a valid Particle was returned before adding it to the event
                    if (part != null) {
                        physEvent.addParticle(part);
                    }
                }
            }

            for (int current_p1 = 0; current_p1 < num_pip; current_p1++) {
                for (int current_p2 = 0; current_p2 < num_pi0; current_p2++) {
                    if (current_p1 == current_p2) {
                        continue;
                    }
                    Particle part = parent_hadron_creation.rhop_check(physEvent, current_p1, current_p2);

                    // Check if a valid Particle was returned before adding it to the event
                    if (part != null) {
                        physEvent.addParticle(part);
                    }
                }
            }

            for (int current_p1 = 0; current_p1 < num_pim; current_p1++) {
                for (int current_p2 = 0; current_p2 < num_pi0; current_p2++) {
                    if (current_p1 == current_p2) {
                        continue;
                    }
                    Particle part = parent_hadron_creation.rhom_check(physEvent, current_p1, current_p2);

                    // Check if a valid Particle was returned before adding it to the event
                    if (part != null) {
                        physEvent.addParticle(part);
                    }
                }
            }
            
            return physEvent;
        }
    return new PhysicsEvent(this.mybeam);
    }
}
    
