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
    
    double Q2(LorentzVector lv) {
        return -lv.mass2();
    }
    
}
