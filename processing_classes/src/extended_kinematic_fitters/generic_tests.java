/**
 *
 * @author Timothy B. Hayward
 */

package extended_kinematic_fitters;

import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

public class generic_tests {
    
    public boolean banks_test(DataEvent event) {
        String[] bankNames = 
            {"RUN::config","REC::Event","REC::Particle","REC::Calorimeter","REC::Track","REC::Traj","REC::Cherenkov"};
        for (String bankName : bankNames) {
            if (!event.hasBank(bankName)) { return false; }
        }
        return true;
    }
    
//    public boolean forward_detector_cut(int particle_Index, HipoDataBank rec_Bank) {
//        int status = rec_Bank.getInt("status", particle_Index);
//        return (Math.abs(status)<4000 || Math.abs(status)>4999) && Math.abs(status)>1999;
//    }
    
    public boolean forward_detector_cut(int particle_Index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", particle_Index);
        return (Math.abs(status)>=2000 && Math.abs(status)<4000);
    }
    
    public boolean central_detector_cut(int particle_Index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", particle_Index);
        return (Math.abs(status)>=4000 && Math.abs(status)<5000);
    }
    
    public boolean forward_tagger_cut(int particle_Index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", particle_Index);
        return (Math.abs(status)>=1000 && Math.abs(status)<2000);
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
    
    public boolean vertex_cut(int particle_Index, HipoDataBank rec_Bank, 
        HipoDataBank run_Bank) {
        // pass2 derived vertex cuts
        int pid = rec_Bank.getInt("pid", particle_Index);
        float vz = rec_Bank.getFloat("vz", particle_Index);
        float vz_e = rec_Bank.getFloat("vz", 0);
        double Delta_vz = vz_e - vz;
        float polarity = run_Bank.getFloat("torus", 0);
        
        if (pid == 11) {
            if (polarity < 0) { return vz > -6 && vz < 1; }
            if (polarity > 0) { return vz > -7 && vz < 0; }
        } else if (pid > 0) { // positive hadrons
            if (polarity < 0) { return Delta_vz > 1.15 - 3*2.44 && Delta_vz < 1.15 + 3*2.44; }
            if (polarity > 0) { return Delta_vz > -0.86 - 3*2.24 && Delta_vz < -0.86 + 3*2.24; }
        } else if (pid < 0) { // negative hadrons
            if (polarity < 0) { return Delta_vz > 0.03 - 3*2.24 && Delta_vz < 0.03 + 3*2.24; }
            if (polarity > 0) { return Delta_vz > 0.38 - 3*2.55 && Delta_vz < 0.38 + 3*2.55; }
        }
        return true;  // track didn't match any pid?
    }
    
    public boolean theta_cut(int particle_Index, HipoDataBank rec_Bank) {
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double r = Math.pow(px*px + py*py + pz*pz, 0.5);
        double theta = (180/Math.PI)*Math.acos(pz/r);
        
        return true;
//        return theta>0 && theta<90;
    }
}
