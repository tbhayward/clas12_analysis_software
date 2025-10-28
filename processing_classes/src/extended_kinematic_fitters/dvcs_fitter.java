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

    ////////////////////////////////////////////////////////////////////////////////////////////////
    public boolean electron_test(int particle_Index, double p,
            HipoDataBank rec_Bank, HipoDataBank cal_Bank,
            HipoDataBank traj_Bank, HipoDataBank run_Bank, HipoDataBank cc_Bank) {

        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();

        return true
                && p > 2.0
                && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
                && generic_tests.vertex_cut(particle_Index, rec_Bank, run_Bank)
//                && pid_cuts.calorimeter_energy_cut(particle_Index, cal_Bank, run_Bank)
                && pid_cuts.calorimeter_sampling_fraction_cut(particle_Index, p, run_Bank, cal_Bank)
                && pid_cuts.calorimeter_diagonal_cut(particle_Index, p, cal_Bank, run_Bank)
                && fiducial_cuts.pcal_fiducial_cut(particle_Index, 2, run_Bank, rec_Bank, cal_Bank)
                && fiducial_cuts.dc_fiducial_cut(particle_Index, rec_Bank, traj_Bank, run_Bank)
                ;
    }

    public boolean proton_test(int particle_Index, int pid, float vz, double trigger_electron_vz,
            HipoDataBank rec_Bank, HipoDataBank cal_Bank,
            HipoDataBank traj_Bank, HipoDataBank run_Bank) {

        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();

        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(px * px + py * py + pz * pz);

        float torus = run_Bank.getFloat("torus", 0);

        boolean passesForwardDetector = generic_tests.forward_detector_cut(particle_Index, rec_Bank);
        boolean passesCentralDetector = generic_tests.central_detector_cut(particle_Index, rec_Bank);

        return true
                && (passesCentralDetector ? p > 0.3 : true)
                && (passesForwardDetector && (torus > 0) ? p > 0.42 : true)
                && (passesForwardDetector && (torus < 0) ? p > 0.50 : true)
                && p < 1.2 // this bound is enforced at p < 1.14 by -t < 1, done here to speed up processing
                && generic_tests.theta_calculation(px, py, pz) < 64.23
                && generic_tests.vertex_cut(particle_Index, rec_Bank, run_Bank)
                && (passesForwardDetector
                        ? fiducial_cuts.dc_fiducial_cut(particle_Index, rec_Bank, traj_Bank, run_Bank)
                        : true)
                && (passesCentralDetector
                        ? fiducial_cuts.cvt_fiducial_cut(particle_Index, rec_Bank, traj_Bank, 2)
                        : true)
//                && (passesForwardDetector // dedicated PID cuts for forward
//                        ? pid_cuts.charged_hadron_chi2pid_cut(particle_Index, rec_Bank, run_Bank)
//                        : true)
//                && (passesCentralDetector
//                        ? pid_cuts.charged_hadron_chi2pid_cut(particle_Index, rec_Bank, run_Bank)
//                        : true) //            && charged_hadron_chi2pid_cut(particle_Index, rec_Bank)
                ;
    }

    public boolean photon_test(int particle_Index, HipoDataBank run_Bank, HipoDataBank rec_Bank,
            HipoDataBank cal_Bank, HipoDataBank ft_Bank, LorentzVector lv_e) {

        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();

        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px, 2) + Math.pow(py, 2) + Math.pow(pz, 2));
        LorentzVector lv_gamma = new LorentzVector();
        lv_gamma.setPxPyPzM(px, py, pz, 0.0);

        boolean passesForwardDetector = generic_tests.forward_detector_cut(particle_Index, rec_Bank);
        boolean passesForwardTagger = generic_tests.forward_tagger_cut(particle_Index, rec_Bank);

        return true
                && p > 2.00
                && (passesForwardDetector || passesForwardTagger)
                && (passesForwardDetector
                        ? fiducial_cuts.pcal_fiducial_cut(particle_Index, 1, run_Bank, rec_Bank, cal_Bank)
                        : fiducial_cuts.forward_tagger_fiducial_cut(particle_Index, rec_Bank, ft_Bank))
                && pid_cuts.beta_cut(particle_Index, rec_Bank)
                ;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    @Override
    public PhysicsEvent getPhysicsEvent(DataEvent event) {

        generic_tests generic_tests = new generic_tests();
        if (generic_tests.banks_test(event)) {
            PhysicsEvent physEvent = new PhysicsEvent();
            // load the hipo banks
            // assumption is we are using trains which would require all of these banks to exist
            HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle");
            HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
            HipoDataBank cc_Bank = (HipoDataBank) event.getBank("REC::Cherenkov");
            HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");
            HipoDataBank track_Bank = (HipoDataBank) event.getBank("REC::Track");
            HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config");
            HipoDataBank ft_Bank = null;
            if (event.hasBank("REC::ForwardTagger")) {
                ft_Bank = (HipoDataBank) event.getBank("REC::ForwardTagger");
            }
            double vz_e = -999;

            LorentzVector lv_e = new LorentzVector();
            if (rec_Bank.getInt("pid", 0) == 11) {
                // trigger particle was an electron
                // highest momentum electron listed first (used for DIS calculations)
                float px = rec_Bank.getFloat("px", 0);
                float py = rec_Bank.getFloat("py", 0);
                float pz = rec_Bank.getFloat("pz", 0);
                double p = Math.sqrt(px * px + py * py + pz * pz);
                lv_e.setPxPyPzM(px, py, pz, 0.0005109989461);
                vz_e = rec_Bank.getFloat("vz", 0);

            } else {
                return physEvent;
            } // trigger particle was not an electron

            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
                int pid = rec_Bank.getInt("pid", particle_Index);
                float px = rec_Bank.getFloat("px", particle_Index);
                float py = rec_Bank.getFloat("py", particle_Index);
                float pz = rec_Bank.getFloat("pz", particle_Index);
                float vx = rec_Bank.getFloat("vx", particle_Index);
                float vy = rec_Bank.getFloat("vy", particle_Index);
                float vz = rec_Bank.getFloat("vz", particle_Index);
                double p = Math.sqrt(px * px + py * py + pz * pz);

                int sector = generic_tests.sector(particle_Index, track_Bank); // 0 FT/CD, 1-6 FD

                int runnum = run_Bank.getInt("run", 0);
                int runPeriod = -1;
                if (runnum >= 4763 && runnum <= 5666) {
                    runPeriod = 3; // this is how Richard Capobianco labeled them in his momentum corrections
                } // RGA Fa18
                else if (runnum >= 6616 && runnum <= 6783) {
                    runPeriod = 2;
                } // RGA Sp19 
                boolean inbending = false;
                boolean outbending = false;
                if (run_Bank.getFloat("torus", 0) == 1) {
                    outbending = true;
                } else {
                    inbending = true;
                }

                energy_loss_corrections energy_loss_corrections = new energy_loss_corrections();
                momentum_corrections momentum_corrections = new momentum_corrections();

                if (pid == 11 && electron_test(particle_Index, p, rec_Bank, cal_Bank,
                        traj_Bank, run_Bank, cc_Bank)) {

                    float[] momentum = {px, py, pz};
//                    if (inbending) {
//                        momentum_corrections.inbending_momentum_corrections(momentum, sector, 0, runPeriod, runPeriod, 0, 0);
//                    } else if (outbending) {
//                        momentum_corrections.outbending_momentum_corrections(momentum, sector, 0, runPeriod, runPeriod, 0, 0);
//                    }
                    px = momentum[0];
                    py = momentum[1];
                    pz = momentum[2];

                    Particle electron = new Particle(pid, px, py, pz, vx, vy, vz_e);
                    physEvent.addParticle(electron);
                }

                if (pid == 2212 && proton_test(particle_Index, pid, vz, vz_e, rec_Bank, cal_Bank,
                        traj_Bank, run_Bank)) {

                    float[] momentum = {px, py, pz};
                    energy_loss_corrections.proton_energy_loss_corrections(particle_Index, momentum, rec_Bank, run_Bank);
//                    energy_loss_corrections.krishna_energy_loss_corrections(particle_Index, momentum, rec_Bank, run_Bank);
//                    energy_loss_corrections.mariana_proton_energy_loss_corrections(particle_Index, momentum, rec_Bank, run_Bank);

                    px = momentum[0];
                    py = momentum[1];
                    pz = momentum[2];
                    Particle part = new Particle(pid, px, py, pz, vx, vy, vz);
                    physEvent.addParticle(part);
                }

                if (pid == 22 && photon_test(particle_Index, run_Bank, rec_Bank, cal_Bank, ft_Bank, lv_e)) {

                    float[] momentum = {px, py, pz};
                    energy_loss_corrections.sebastian_photon_energy_loss_corrections(particle_Index, momentum, rec_Bank, run_Bank);
//
                    px = momentum[0];
                    py = momentum[1];
                    pz = momentum[2];
                    Particle part = new Particle(pid, px, py, pz, vx, vy, vz);
                    physEvent.addParticle(part);
                }
            }

            int num_gamma = physEvent.countByPid(22);
            for (int current_p1 = 0; current_p1 < num_gamma; current_p1++) {
                for (int current_p2 = 0; current_p2 < num_gamma; current_p2++) {
                    parent_hadron_creation parent_hadron_creation = new parent_hadron_creation();
                    Particle part = parent_hadron_creation.pi0_check(physEvent, current_p1, current_p2);

                    // Check if a valid Particle was returned before adding it to the event
                    if (part != null) {
                        physEvent.addParticle(part);
                    }
                }
            }

            return physEvent;
        }
        return new PhysicsEvent(this.mybeam);
    }
}
