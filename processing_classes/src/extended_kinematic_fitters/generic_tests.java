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
    
    public boolean forward_detector_cut(int particle_Index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", particle_Index);
        return (Math.abs(status)<4000 || Math.abs(status)>4999) && Math.abs(status)>1999;
    }
    
    public boolean nphe_cut(int particle_Index, HipoDataBank cc_Bank) {
        for (int current_Row = 0; current_Row < cc_Bank.rows(); current_Row++) {
            if (cc_Bank.getInt("pindex", current_Row)==particle_Index) {
                return cc_Bank.getFloat("nphe", current_Row) > 2;
            }
        }
        return false; 
    }
    
    public int rich_detector_pid(int particle_Index, HipoDataBank rich_Bank) {
        // Iterate through the rows of the data bank
        for (int current_Row = 0; current_Row < rich_Bank.rows(); current_Row++) {
            // Get the pindex for the current row
            int pindex = rich_Bank.getInt("pindex", current_Row);
            
            // Check if the pindex value matches the specified particle
            if (pindex == particle_Index) {
                return rich_Bank.getInt("best_PID", current_Row);
            }
        }
        return 0;
    }
    
}
