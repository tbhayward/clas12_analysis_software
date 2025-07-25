/**
 *
 * @author Timothy B. Hayward
 */
package extended_kinematic_fitters;

import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

public class generic_tests {

    public boolean banks_test(DataEvent event) {
        String[] bankNames
                = {"RUN::config", "REC::Event", "REC::Particle", "REC::Calorimeter", "REC::Traj", "REC::Track", "REC::Cherenkov"};
        for (String bankName : bankNames) {
            if (!event.hasBank(bankName)) {
                return false;
            }
        }
        return true;
    }

    public static double p_calculation(double x, double y, double z) {
        return Math.sqrt(x * x + y * y + z * z);
    }

    public static double phi_calculation(double x, double y) {
        // tracks are given with Cartesian values and so must be converted to cylindrical
        double phi = Math.toDegrees(Math.atan2(x, y));
        phi = phi - 90;
        if (phi < 0) {
            phi = 360 + phi;
        }
        phi = 360 - phi;
        return phi;
    }

    public static double theta_calculation(double x, double y, double z) {
        // convert cartesian coordinates to polar angle
        double r = Math.pow(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2), 0.5);
        return (double) (180 / Math.PI) * Math.acos(z / r);
    }

    public static double x_calculation(double p, double theta, double phi) {
        // Convert angles from degrees to radians
        theta = Math.toRadians(theta);
        phi = Math.toRadians(phi);

        // Calculate x using spherical to cartesian conversion
        return p * Math.sin(theta) * Math.cos(phi);
    }

    public static double y_calculation(double p, double theta, double phi) {
        // Convert angles from degrees to radians
        theta = Math.toRadians(theta);
        phi = Math.toRadians(phi);

        // Calculate y using spherical to cartesian conversion
        return p * Math.sin(theta) * Math.sin(phi);
    }

    public static double z_calculation(double p, double theta) {
        // Convert angle from degrees to radians
        theta = Math.toRadians(theta);

        // Calculate z using spherical to cartesian conversion
        return p * Math.cos(theta);
    }

    public boolean forward_detector_cut(int particle_Index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", particle_Index);
        return (Math.abs(status) >= 2000 && Math.abs(status) < 4000);
    }

    public boolean central_detector_cut(int particle_Index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", particle_Index);
        return (Math.abs(status) >= 4000 && Math.abs(status) < 5000);
    }

    public boolean forward_tagger_cut(int particle_Index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", particle_Index);
        return (Math.abs(status) >= 1000 && Math.abs(status) < 2000);
    }

    public int sector(int particle_Index, HipoDataBank track_Bank) {
        for (int current_Row = 0; current_Row < track_Bank.rows(); current_Row++) {
            // Get the pindex and layer values for the current row
            int pindex = track_Bank.getInt("pindex", current_Row);
            if (pindex == particle_Index) {
                return track_Bank.getInt("sector", current_Row);
            }
        }
        return -1; // no match found?
    }

    public boolean nphe_cut(int particle_Index, HipoDataBank cc_Bank) {
        for (int current_Row = 0; current_Row < cc_Bank.rows(); current_Row++) {
            if (cc_Bank.getInt("pindex", current_Row) == particle_Index) {
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
        } else if (runnum >= 4003 && runnum <= 4325) { // RGA Sp18 Inb
            if (charge > 0) {
                return -7.879 < vz && vz < 1.515;
            } else if (charge < 0) {
                return -6.0606 < vz && vz < 1.8182;
            }
        } else if ((runnum >= 3173 && runnum <= 3293) || (runnum >= 3863 && runnum <= 3987)) { // RGA Sp18 Out
            if (charge > 0) {
                return -6.6667 < vz && vz < 2.7273;
            } else if (charge < 0) {
                return -7.273 < vz && vz < 0.9091;
            }
        } else if ((runnum >= 5032 && runnum <= 5419) || (runnum >= 6616 && runnum <= 6783)) { // RGA Fa18 Inb and RGA Sp19 Inb
            if (charge > 0) {
                return -8.485 < vz && vz < 0.606;
            } else if (charge < 0) {
                return -6.364 < vz && vz < 1.515;
            }
        } else if (runnum >= 5422 && runnum <= 5666) { // RGA Fa18 Out
            if (charge > 0) {
                return -6.970 < vz && vz < 1.818;
            } else if (charge < 0) {
                return -7.879 < vz && vz < 0.303;
            }
        } else if (runnum >= 16043 && runnum <= 16772) { // RGC Su22
            if (charge > 0) {
                return -9.394 < vz && vz < -0.606;
            } else if (charge < 0) {
                return -7.576 < vz && vz < 0.303;
            }
        } else if (runnum >= 16843 && runnum <= 17811) { // RGC Fa22 & RGC Sp23
            if (charge > 0) {
                return -8.788 < vz && vz < 0.303;
            } else if (charge < 0) {
                return -5.758 < vz && vz < 1.515;
            }
        }

        return -9 < vz && vz < 2;
    }

}
