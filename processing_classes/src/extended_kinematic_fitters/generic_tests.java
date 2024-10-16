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
            {"RUN::config","REC::Event","REC::Particle","REC::Calorimeter","REC::Traj","REC::Cherenkov"};
        for (String bankName : bankNames) {
            if (!event.hasBank(bankName)) { System.out.println(bankName); return false; }
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
    
    public boolean vertex_cut(int particle_Index, HipoDataBank rec_Bank, HipoDataBank run_Bank) {
        // pass2 derived vertex cuts
        int charge = rec_Bank.getInt("charge", particle_Index);
        float vz = rec_Bank.getFloat("vz", particle_Index);
        
        int runnum = run_Bank.getInt("run", 0);
        
        if (runnum == 11) {
            if (charge > 0) {
                return -10 < vz && vz < 1.5;
            } else if (charge < 0) {
                return -9 < vz && vz < 2;
            }
        } else if (runnum >= 6616 && runnum <= 6783) { // RGC Su22
            if (charge > 0) {
               return -10 < vz && vz < 1.5; 
            } else if (charge < 0) {
               return -9 < vz && vz < 2;
            }
        }
        
        return -11 < vz && vz < 3;
    }
    
}
