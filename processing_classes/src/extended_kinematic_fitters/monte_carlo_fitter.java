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
            return physEvent;
        }
    return new PhysicsEvent(this.mybeam);
    }
}
    
