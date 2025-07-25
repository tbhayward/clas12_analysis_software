package analyzers;

/**
 *
 * @author tbhayward
 */
import extended_kinematic_fitters.fiducial_cuts;
import extended_kinematic_fitters.generic_tests;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.clas.physics.*;

public class Inclusive {

    protected byte helicity;
    protected int runnum;

    protected int fiducial_status = -1;

    protected int num_elec, num_piplus, num_piminus, num_kplus, num_kminus, num_protons, num_particles;
    protected int num_pos, num_neg, num_neutrals;
    protected int num_positrons, num_antiprotons;

    protected double Q2, W, gamma, nu, x, y, t, tmin, Mx, Mx2;

    protected double e_px, e_py, e_pz, e_p, e_e, e_theta, e_phi, vz_e; // electron kinematics

    // depolarization vectors defining the polarization lost during the transfer from beam to 
    // the virtual photon. 
    // in ALU BSAs the twist 2 terms are modified by C/A and the twist 3 terms by W/A
    // B and V come in AUL
    protected double Depolarization_A;
    protected double Depolarization_B;
    protected double Depolarization_C;
    protected double Depolarization_V;
    protected double Depolarization_W;

    public static boolean channel_test(Inclusive variables) {
//        if (variables.helicity == 0 && variables.runnum != 11) {
//            return false;
//        }
        if (variables.Q2() < 1) {
            return false;
        }
        if (variables.W() < 2) {
            return false;
        }
        if (variables.y() > 0.80) {
            return false;
        }
        return true;
    }

    public Inclusive(DataEvent event, PhysicsEvent recEvent, double Eb) {
        // provide the PDG PID of the two hadrons

        kinematic_variables kinematic_variables = new kinematic_variables();

        // load banks
        HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
        HipoDataBank configBank = (HipoDataBank) event.getBank("RUN::config");
        HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle");
        HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
        HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");

        helicity = eventBank.getByte("helicity", 0);
        runnum = configBank.getInt("run", 0); // used for beam energy and polarization

        num_elec = recEvent.countByPid(11); // returns number of electrons
        num_positrons = recEvent.countByPid(-11); // returns number of positrons
        num_piplus = recEvent.countByPid(211);
        num_piminus = recEvent.countByPid(-211);
        num_kplus = recEvent.countByPid(321);
        num_kminus = recEvent.countByPid(-321);
        num_protons = recEvent.countByPid(2212);
        num_antiprotons = recEvent.countByPid(-2212);
        num_particles = num_elec + num_piplus + num_piminus + num_kplus + num_kminus + num_protons;
        num_pos = num_positrons + num_piplus + num_kplus + num_protons;
        num_neg = num_elec + num_piminus + num_kminus + num_antiprotons;
        num_neutrals = recEvent.countByPid(22) + recEvent.countByPid(2112);

        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();

        boolean electron_pcal_fiducial = fiducial_cuts.pcal_fiducial_cut(0, 1, configBank, rec_Bank, cal_Bank);
        boolean electron_fd_fiducial = fiducial_cuts.dc_fiducial_cut(0, rec_Bank, traj_Bank, configBank);
        boolean e_fiducial_check = electron_pcal_fiducial && electron_fd_fiducial;

        // Check if all checks pass
        if (e_fiducial_check) {
            fiducial_status = 1; // Set to 1 if electron checks pass
        }

        // Set up Lorentz vectors
        // beam electron
        LorentzVector lv_beam = new LorentzVector();
        lv_beam.setPxPyPzM(0, 0, Math.pow(Eb * Eb - kinematic_variables.particle_mass(11) * kinematic_variables.particle_mass(11), 0.5),
                kinematic_variables.particle_mass(11));
        LorentzVector lv_target = new LorentzVector();
        // target, proton for RGA... what mass to use for RGB (deuterium target)?
        double target_pz = 0;
//        double target_pz = kinematic_variables.rand_Gaussian(0,0.033);
//        double target_pz = kinematic_variables.rand_table(0,0.033);
        lv_target.setPxPyPzM(0, 0, target_pz, kinematic_variables.particle_mass(2212));
        // pull from rec banks for outgoing particles
        // electron
        String electron_index = "[11,0]"; // highest p, kinematic fitter should require FD etc
        Particle scattered_electron = recEvent.getParticle(electron_index); //
        LorentzVector lv_e = new LorentzVector();
        lv_e.setPxPyPzM(scattered_electron.px(), scattered_electron.py(),
                scattered_electron.pz(), kinematic_variables.particle_mass(11));
        // hadrons set up below (to allow for iteration over more than two hadrons in an event)

        // kinematics of electron
        e_px = lv_e.px();
        e_py = lv_e.py();
        e_pz = lv_e.pz();
        e_p = lv_e.p();
        e_e = lv_e.e();
        e_theta = scattered_electron.theta();
        e_phi = scattered_electron.phi();
        if (e_phi < 0) {
            e_phi = 2 * Math.PI + e_phi;
        }

        // DIS variables
        LorentzVector lv_q = new LorentzVector(lv_beam);
        lv_q.sub(lv_e);
        Q2 = kinematic_variables.Q2(lv_q);
        nu = kinematic_variables.nu(lv_beam, lv_e);
        x = kinematic_variables.x(Q2, nu);
        W = kinematic_variables.W(Q2, nu);
        y = kinematic_variables.y(nu, lv_beam);
        gamma = kinematic_variables.gamma(Q2, x);
        t = kinematic_variables.t(lv_e.p(), e_theta);
        tmin = kinematic_variables.tmin(x);

        // Depolarization variables
        Depolarization_A = kinematic_variables.Depolarization_A(gamma, y);
        Depolarization_B = kinematic_variables.Depolarization_B(gamma, y);
        Depolarization_C = kinematic_variables.Depolarization_C(gamma, y);
        Depolarization_V = kinematic_variables.Depolarization_V(gamma, y);
        Depolarization_W = kinematic_variables.Depolarization_W(gamma, y);

        vz_e = scattered_electron.vz();

        // missing mass calculations
        Mx = kinematic_variables.Mx(lv_q, lv_target);
        Mx2 = kinematic_variables.Mx2(lv_q, lv_target);  // missing mass squared
    }

    public int get_helicity() { // -1, 0, or 1. 0 equals unassigned by EventBuilder
        if (runnum >= 4326 && runnum <= 5666) {
            return -1 * helicity;
        } else if (runnum >= 6616 && runnum <= 6783) {
            return -1 * helicity;
        } else if (runnum >= 6120 && runnum <= 6604) {
            return -1 * helicity;
        } else if (runnum >= 11093 && runnum <= 11283) {
            return helicity;
        } else if (runnum >= 11284 && runnum < 11300) {
            return -1 * helicity;
        } else if (runnum >= 11323 && runnum < 11571) {
            return helicity;
        }
        return helicity;
    }

    public int get_runnum() {
        return runnum;
    }

    public int get_num_pos() {
        return num_pos;
    }

    public int get_num_neg() {
        return num_neg;
    }

    public int get_num_neutrals() {
        return num_neutrals;
    }

    public int get_fiducial_status() {
        return fiducial_status;
    }

    public int num_elec() {
        return num_elec;
    } // returns number of electrons

    public int num_piplus() {
        return num_piplus;
    } // returns number of piplus

    public int num_piminus() {
        return num_piminus;
    } // returns number of piminus

    public int num_kplus() {
        return num_kplus;
    }// returns number of kplus

    public int num_kminus() {
        return num_kminus;
    } // returns number of kminus

    public int num_protons() {
        return num_protons;
    } // returns number of protons

    public double Q2() {
        return Double.valueOf(Math.round(Q2 * 100000)) / 100000;
    } // returns Q2

    public double W() {
        return Double.valueOf(Math.round(W * 100000)) / 100000;
    }// returns W

    public double gamma() {
        return Double.valueOf(Math.round(gamma * 100000)) / 100000;
    } // returns gamma

    public double nu() {
        return Double.valueOf(Math.round(nu * 100000)) / 100000;
    }// returns nu

    public double x() {
        return Double.valueOf(Math.round(x * 100000)) / 100000;
    }// returns x

    public double y() {
        return Double.valueOf(Math.round(y * 100000)) / 100000;
    }// returns y

    public double t() {
        return Double.valueOf(Math.round(t * 100000)) / 100000;
    }// returns t

    public double tmin() {
        return Double.valueOf(Math.round(tmin * 100000)) / 100000;
    }// returns t

    public double Mx() {
        return Double.valueOf(Math.round(Mx * 100000)) / 100000;
    }// returns Mx(ep1p2)

    public double Mx2() {
        return ((int) (Mx2 * 100000)) / 100000.0;
    }

    public double e_px() {
        return Double.valueOf(Math.round(e_px * 100000)) / 100000;
    }// returns electron lab frame px

    public double e_py() {
        return Double.valueOf(Math.round(e_py * 100000)) / 100000;
    }// returns electron lab frame py

    public double e_pz() {
        return Double.valueOf(Math.round(e_pz * 100000)) / 100000;
    }// returns electron lab frame pz

    public double e_p() {
        return Double.valueOf(Math.round(e_p * 100000)) / 100000;
    }// returns electron lab frame p

    public double e_e() {
        return Double.valueOf(Math.round(e_e * 100000)) / 100000;
    }// returns electron lab frame energy

    public double e_theta() {
        return Double.valueOf(Math.round(e_theta * 100000)) / 100000;
    } // returns electron lab 
    // frame polar angle

    public double e_phi() {
        return Double.valueOf(Math.round(e_phi * 100000)) / 100000;
    } // returns electron lab 
    // frame polar angle

    public double Depolarization_A() {
        return ((int) (Depolarization_A * 100000)) / 100000.0;
    }

    public double Depolarization_B() {
        return ((int) (Depolarization_B * 100000)) / 100000.0;
    }

    public double Depolarization_C() {
        return ((int) (Depolarization_C * 100000)) / 100000.0;
    }

    public double Depolarization_V() {
        return ((int) (Depolarization_V * 100000)) / 100000.0;
    }

    public double Depolarization_W() {
        return ((int) (Depolarization_W * 100000)) / 100000.0;
    }

    public double vz_e() {
        return Double.valueOf(Math.round(vz_e * 100000)) / 100000;
    }// returns electron z vertex

}
