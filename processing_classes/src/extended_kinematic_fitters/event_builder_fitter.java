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

public class event_builder_fitter extends GenericKinematicFitter {

    protected final Double mybeam;

    public event_builder_fitter(double beam) {
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
        if (!(event.hasBank("REC::Particle"))) {
            banks_test = false;
        } 
        if (banks_test) {
            PhysicsEvent physEvent = new PhysicsEvent();
            HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Particle"); // load particle bank
            for (int current_part = 0; current_part < eventBank.rows(); current_part++) {
                int pid = eventBank.getInt("pid", current_part);
                if (pid!=0) { // do not load tracks that the EventBuilder did not assign a pid to 
                    float vx = eventBank.getFloat("vx",current_part);
                    float vy = eventBank.getFloat("vy",current_part);
                    float vz = eventBank.getFloat("vz",current_part);
                    float px = eventBank.getFloat("px", current_part);
                    float py = eventBank.getFloat("py", current_part);
                    float pz = eventBank.getFloat("pz", current_part);
                    Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
                    physEvent.addParticle(part);
                }
            }
            return physEvent;
        }
    return new PhysicsEvent(this.mybeam);
    }
}
    
