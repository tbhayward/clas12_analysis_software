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

import org.jlab.clas.physics.*;

public class generic_tests {
    
    public boolean banks_test(DataEvent event) {
        String[] bankNames = 
            {"RUN::config","REC::Particle","REC::Calorimeter","REC::Track","REC::Traj","REC::Cherenkov"};
        for (String bankName : bankNames) {
            if (!event.hasBank(bankName)) { return false; }
        }
        return true;
    }
    
}
