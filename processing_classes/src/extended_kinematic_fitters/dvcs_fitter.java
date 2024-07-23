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


public class dvcs_fitter extends GenericKinematicFitter {

    protected final Double mybeam;

    public dvcs_fitter(double beam) {
        super(beam);
        mybeam = beam;
    } 



}