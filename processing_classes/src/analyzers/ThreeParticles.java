package analyzers;

/**
 *
 * @author tbhayward
 */
import extended_kinematic_fitters.fiducial_cuts;
import extended_kinematic_fitters.generic_tests;
import extended_kinematic_fitters.pid_cuts;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.clas.physics.*;

public class ThreeParticles {

    protected byte helicity;
    protected int runnum;

    protected int e_pid_cuts = 0;
    protected int fiducial_status = -1;
    protected int pid_status = -1;
    protected int detector1 = -1;
    protected int detector2 = -1;

    protected int num_elec, num_piplus, num_piminus, num_kplus, num_kminus, num_protons, num_particles;
    protected int num_pos, num_neg, num_neutrals;
    protected int num_positrons, num_antiprotons;

    // labels are unnumbered if they refer to the dihadron (perhaps a meson) and numbered for individual
    // hadrons. Convention is ordered by mass, then charge. For example in pi+pi- pi+ is hadron 1
    // in proton+pi+ the proton is p1, in k+pi- the kaon is p1.
    protected double Q2, W, gamma, nu, x, y, t, t1, t2, tmin, z, z1, z2;
    protected double Mx2, Mx2_1, Mx2_2; // Mx is the Mx(ep1p2), Mx1 is Mx(ep1), Mx2 is Mx(ep2), Mx3 is Mx(e)
    protected double Mh, pT, pT1, pT2, xF, xF1, xF2, zeta, zeta1, zeta2, xi, xi1, xi2;
    protected double eta, eta1, eta2, eta_gN, eta1_gN, eta2_gN;
    // eta is the rapidity, preferred by theorists in the Breit frame (e.g. eta1 is in Breit) 
    // eta_gN is the rapidity in the gamma*-nucleon COM frame
    // the difference between two rapidities is Lorentz invariant, i.e.
    // eta1-eta2 = eta1_COM - eta2_COM

    // theta is defined as the angle between the hadron p1 in the pair center-of-mass frame 
    // and the direction of the pair, Ph, in the photon-target rest frame.
    // phih and phiR are defined from the vectors Ph = p1+p2 and 2R = p1-p2, see 
    // hep-ex:2101.04842
    protected double theta, phih, phiR, phi1, phi2, Delta_phi;

    // exclusivity
    protected double Emiss0, Emiss1, Emiss2, Emiss3;
    protected double theta_gamma_gamma;
    protected double pTmiss;

    // depolarization vectors defining the polarization lost during the transfer from beam to 
    // the virtual photon. 
    // in ALU BSAs the twist 2 terms are modified by C/A and the twist 3 terms by W/A
    // B and V come in AUL
    protected double Depolarization_A;
    protected double Depolarization_B;
    protected double Depolarization_C;
    protected double Depolarization_V;
    protected double Depolarization_W;

    protected double e_px, e_py, e_pz, e_p, e_e, e_theta, e_phi; // electron kinematics
    protected double p_px, p_py, p_pz, p_p, p_e; // dihadron kinematics, mass is Mh
    protected double p1_px, p1_py, p1_pz, p1_p, p1_e, p1_theta, p1_phi; // p1 kinematics
    protected double p2_px, p2_py, p2_pz, p2_p, p2_e, p2_theta, p2_phi; // p2 kinematics
    protected double vx_e, vx_p1, vx_p2, vy_e, vy_p1, vy_p2, vz_e, vz_p1, vz_p2;
    protected double open_angle_ep, open_angle_ep1, open_angle_ep2, open_angle_p1p2;

    protected double p_Breit_pz, p1_Breit_pz, p2_Breit_pz, p_gN_pz, p1_gN_pz, p2_gN_pz;

    protected float e_nphe, p1_nphe, p2_nphe; // number of photoelectrons in cherenkov counter
    protected float e_chi2pid, p1_chi2pid, p2_chi2pid; // chi2pid of the CLAS12 EventBuilder
    // defined as the sampling fraction for electron candidates and the distance away from the 
    // mean beta(p) value from the TOF for hadron candidates

    protected double p1_COM_phi, p1_COM_theta, p2_COM_phi, p2_COM_theta, COM_Delta_phi, COM_Delta_theta;
    protected double COM_open_angle;

    protected double gN_angle_p1_p2, gN_angle_p1_X, gN_angle_p2_X;

    // RICH variables
    protected int emilay1 = -9999;
    protected int emico1 = -9999;
    protected int emqua1 = -9999;
    protected int best_PID1 = -9999;
    protected float RQ1 = -9999;
    protected float ReQ1 = -9999;
    protected float el_logl1 = -9999;
    protected float pi_logl1 = -9999;
    protected float k_logl1 = -9999;
    protected float pr_logl1 = -9999;
    protected float best_ch1 = -9999;
    protected float best_c21 = -9999;
    protected float best_RL1 = -9999;
    protected float best_ntot1 = -9999;

    protected int emilay2 = -9999;
    protected int emico2 = -9999;
    protected int emqua2 = -9999;
    protected int best_PID2 = -9999;
    protected float RQ2 = -9999;
    protected float ReQ2 = -9999;
    protected float el_logl2 = -9999;
    protected float pi_logl2 = -9999;
    protected float k_logl2 = -9999;
    protected float pr_logl2 = -9999;
    protected float best_ch2 = -9999;
    protected float best_c22 = -9999;
    protected float best_RL2 = -9999;
    protected float best_ntot2 = -9999;

    public static boolean channel_test(ThreeParticles variables) {
//        if (variables.helicity == 0 && variables.runnum != 11) {
//            return false;
//        }
//        if (variables.Q2() < 1.00) {
//            return false;
//        } 
//        else if (variables.W() < 2) {
//            return false;
//        } 
//        else if (variables.y() > 0.80) {
//            return false;
//        }
        return true;
    }

    public static int getIndex(HipoDataBank rec_Bank, int input_pid, int input_index) {
        int index = -1;
        for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
            int pid = rec_Bank.getInt("pid", particle_Index);
            if (pid == input_pid) {
                index++;
            }
            if (index == input_index) {
                return particle_Index;
            }
        }
        return -1;
    }

    public ThreeParticles(DataEvent event, PhysicsEvent recEvent, int p1PID, int p1Index, int p2PID, int p2Index,
            double Eb) {
        // provide the PDG PID of the two hadrons

        kinematic_variables kinematic_variables = new kinematic_variables();

        // load banks
        HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
        HipoDataBank configBank = (HipoDataBank) event.getBank("RUN::config");
        HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle");
        HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
        HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");
        int p1_rec_index = getIndex(rec_Bank, p1PID, p1Index);
        int p2_rec_index = getIndex(rec_Bank, p2PID, p2Index);
        HipoDataBank rich_Bank = null;
        if (event.hasBank("RICH::Particle")) {
            rich_Bank = (HipoDataBank) event.getBank("RICH::Particle");

            for (int row = 0; row < rich_Bank.rows(); row++) {
                int pindex = rich_Bank.getInt("pindex", row);

                // proceed if this row is for p1 or p2
                if (pindex == p1_rec_index || pindex == p2_rec_index) {
                    // choose suffix: 1 for p1, 2 for p2
                    boolean isFirst = (pindex == p1_rec_index);

                    //  assign into either the “1” variables or the “2” variables
                    if (isFirst) {
                        emilay1 = rich_Bank.getByte("emilay", row);
                        emico1 = rich_Bank.getByte("emico", row);
                        emqua1 = rich_Bank.getShort("emqua", row);
                        best_PID1 = rich_Bank.getShort("best_PID", row);
                        RQ1 = rich_Bank.getFloat("RQ", row);
                        ReQ1 = rich_Bank.getFloat("ReQ", row);
                        el_logl1 = rich_Bank.getFloat("el_logl", row);
                        pi_logl1 = rich_Bank.getFloat("pi_logl", row);
                        k_logl1 = rich_Bank.getFloat("k_logl", row);
                        pr_logl1 = rich_Bank.getFloat("pr_logl", row);
                        best_ch1 = rich_Bank.getFloat("best_ch", row);
                        best_c21 = rich_Bank.getFloat("best_c2", row);
                        best_RL1 = rich_Bank.getFloat("best_RL", row);
                        best_ntot1 = rich_Bank.getFloat("best_ntot", row);
                    } else {
                        emilay2 = rich_Bank.getByte("emilay", row);
                        emico2 = rich_Bank.getByte("emico", row);
                        emqua2 = rich_Bank.getShort("emqua", row);
                        best_PID2 = rich_Bank.getShort("best_PID", row);
                        RQ2 = rich_Bank.getFloat("RQ", row);
                        ReQ2 = rich_Bank.getFloat("ReQ", row);
                        el_logl2 = rich_Bank.getFloat("el_logl", row);
                        pi_logl2 = rich_Bank.getFloat("pi_logl", row);
                        k_logl2 = rich_Bank.getFloat("k_logl", row);
                        pr_logl2 = rich_Bank.getFloat("pr_logl", row);
                        best_ch2 = rich_Bank.getFloat("best_ch", row);
                        best_c22 = rich_Bank.getFloat("best_c2", row);
                        best_RL2 = rich_Bank.getFloat("best_RL", row);
                        best_ntot2 = rich_Bank.getFloat("best_ntot", row);
                    }
                }
            }
        }

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
        pid_cuts pid_cuts = new pid_cuts();

        boolean electron_pcal_fiducial = fiducial_cuts.pcal_fiducial_cut(0, 1, configBank, rec_Bank, cal_Bank);
        boolean electron_fd_fiducial = fiducial_cuts.dc_fiducial_cut(0, rec_Bank, traj_Bank, configBank);
        boolean e_fiducial_check = electron_pcal_fiducial && electron_fd_fiducial;

        boolean passesForwardDetector_1 = generic_tests.forward_detector_cut(p1_rec_index, rec_Bank)
                ? fiducial_cuts.dc_fiducial_cut(p1_rec_index, rec_Bank, traj_Bank, configBank) : true;
        boolean passesCentralDetector_1 = generic_tests.central_detector_cut(p1_rec_index, rec_Bank)
                ? fiducial_cuts.cvt_fiducial_cut(p1_rec_index, rec_Bank, traj_Bank, 1) : true;
        boolean passesForwardTagger_1 = generic_tests.forward_tagger_cut(p1_rec_index, rec_Bank)
                ? fiducial_cuts.forward_tagger_fiducial_cut(p1_rec_index, rec_Bank, cal_Bank) : true;
        boolean p1_fiducial_check = passesForwardTagger_1 && passesForwardDetector_1 && passesCentralDetector_1;

        boolean passesForwardDetector_2 = generic_tests.forward_detector_cut(p2_rec_index, rec_Bank)
                ? fiducial_cuts.dc_fiducial_cut(p2_rec_index, rec_Bank, traj_Bank, configBank) : true;
        boolean passesCentralDetector_2 = generic_tests.central_detector_cut(p2_rec_index, rec_Bank)
                ? fiducial_cuts.cvt_fiducial_cut(p2_rec_index, rec_Bank, traj_Bank, 1) : true;
        boolean passesForwardTagger_2 = generic_tests.forward_tagger_cut(p2_rec_index, rec_Bank)
                ? fiducial_cuts.forward_tagger_fiducial_cut(p2_rec_index, rec_Bank, cal_Bank) : true;
        boolean p2_fiducial_check = passesForwardTagger_2 && passesForwardDetector_2 && passesCentralDetector_2;

        // Check if all checks pass
        if (e_fiducial_check && p1_fiducial_check && p2_fiducial_check) {
            fiducial_status = 3; // Set to 3 if all checks pass
        } else {
            // Now check for specific cases where only one is false
            if (!e_fiducial_check && p1_fiducial_check && p2_fiducial_check) {
                fiducial_status = 0; // Set to 0 if only electron check is false
            } else if (e_fiducial_check && !p1_fiducial_check && p2_fiducial_check) {
                fiducial_status = 1; // Set to 1 if only p1 check is false
            } else if (e_fiducial_check && p1_fiducial_check && !p2_fiducial_check) {
                fiducial_status = 2; // Set to 2 if only p2 check is false (same status as p1)
            }
            // If more than one is false, fiducial_status remains -1 (default)
        }

        if (generic_tests.forward_tagger_cut(p1_rec_index, rec_Bank)) {
            detector1 = 0; // Forward Tagger
        } else if (generic_tests.forward_detector_cut(p1_rec_index, rec_Bank)) {
            detector1 = 1; // Forward Detector
        } else if (generic_tests.central_detector_cut(p1_rec_index, rec_Bank)) {
            detector1 = 2; // Central Detector
        }

        if (generic_tests.forward_tagger_cut(p2_rec_index, rec_Bank)) {
            detector2 = 0; // Forward Tagger
        } else if (generic_tests.forward_detector_cut(p2_rec_index, rec_Bank)) {
            detector2 = 1; // Forward Detector
        } else if (generic_tests.central_detector_cut(p2_rec_index, rec_Bank)) {
            detector2 = 2; // Central Detector
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

        if (pid_cuts.calorimeter_energy_cut(0, cal_Bank, 1)
                && pid_cuts.calorimeter_sampling_fraction_cut(0, lv_e.p(), configBank, cal_Bank)
                && pid_cuts.calorimeter_diagonal_cut(0, lv_e.p(), cal_Bank)) {
            e_pid_cuts = 1;
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

        // Depolarization variables
        Depolarization_A = kinematic_variables.Depolarization_A(gamma, y);
        Depolarization_B = kinematic_variables.Depolarization_B(gamma, y);
        Depolarization_C = kinematic_variables.Depolarization_C(gamma, y);
        Depolarization_V = kinematic_variables.Depolarization_V(gamma, y);
        Depolarization_W = kinematic_variables.Depolarization_W(gamma, y);

        // set up boost to gamma*-nucleon center of mass frame
        LorentzVector gN = new LorentzVector(lv_q);
        gN.add(lv_target);
        Vector3 gNBoost = gN.boostVector();
        gNBoost.negative();

        // set up boost to Breit frame, this needs to be cross checked
        LorentzVector Breit = new LorentzVector(lv_q);
        LorentzVector Breit_target = new LorentzVector();
        Breit_target.setPxPyPzM(0, 0, 0, 2 * x * kinematic_variables.particle_mass(2212));
        Breit.add(Breit_target);
        Vector3 BreitBoost = Breit.boostVector();
        BreitBoost.negative();

        // set up hadrons and dihadron
        String p1Index_string = "[" + p1PID + "," + p1Index + "]";
        Particle p1 = recEvent.getParticle(p1Index_string);
        String p2Index_string = "[" + p2PID + "," + p2Index + "]";
        Particle p2 = recEvent.getParticle(p2Index_string);
        String combined_index_string = "[" + p1PID + "," + p1Index + "]+[" + p2PID + "," + p2Index + "]";
        Particle dihadron = recEvent.getParticle(combined_index_string);
        Mh = dihadron.mass();

        vx_e = scattered_electron.vx();
        vx_p1 = p1.vx();
        vx_p2 = p2.vx();
        vy_e = scattered_electron.vy();
        vy_p1 = p1.vy();
        vy_p2 = p2.vy();
        vz_e = scattered_electron.vz();
        vz_p1 = p1.vz();
        vz_p2 = p2.vz();

        LorentzVector lv_p = new LorentzVector();
        lv_p.setPxPyPzM(dihadron.px(), dihadron.py(), dihadron.pz(), dihadron.mass());
        LorentzVector lv_p1 = new LorentzVector();
        lv_p1.setPxPyPzM(p1.px(), p1.py(), p1.pz(), p1.mass());
        LorentzVector lv_p2 = new LorentzVector();
        lv_p2.setPxPyPzM(p2.px(), p2.py(), p2.pz(), p2.mass());

        t = kinematic_variables.t(lv_p.p(), lv_p.theta());
        t1 = kinematic_variables.t(lv_p1.p(), lv_p1.theta());
        t2 = kinematic_variables.t(lv_p2.p(), lv_p2.theta());
        tmin = kinematic_variables.tmin(x);

        open_angle_ep = kinematic_variables.open_angle(lv_e, lv_p);
        open_angle_ep1 = kinematic_variables.open_angle(lv_e, lv_p1);
        open_angle_ep2 = kinematic_variables.open_angle(lv_e, lv_p2);
        open_angle_p1p2 = kinematic_variables.open_angle(lv_p1, lv_p2);

        Emiss2 = kinematic_variables.Emiss2(lv_beam, lv_target, lv_e, lv_p1, lv_p2);
        theta_gamma_gamma = kinematic_variables.theta_gamma_gamma(lv_beam, lv_target, lv_e, lv_p1, lv_p2);
        pTmiss = kinematic_variables.pTmiss(lv_beam, lv_target, lv_e, lv_p1, lv_p2);

        // kinematics of hadrons
        p1_px = lv_p1.px();
        p1_py = lv_p1.py();
        p1_pz = lv_p1.pz();
        p1_p = lv_p1.p();
        p1_e = p1.e();
        p1_theta = p1.theta();
        p1_phi = p1.phi();
        if (p1_phi < 0) {
            p1_phi = 2 * Math.PI + p1_phi;
        }

        p2_px = lv_p2.px();
        p2_py = lv_p2.py();
        p2_pz = lv_p2.pz();
        p2_p = lv_p2.p();
        p2_e = p2.e();
        p2_theta = p2.theta();
        p2_phi = p2.phi();
        if (p2_phi < 0) {
            p2_phi = 2 * Math.PI + p2_phi;
        }
        p_px = lv_p.px();
        p_py = lv_p.py();
        p_pz = lv_p.pz();
        p_p = lv_p.p();
        p_e = lv_p.e();

        z = kinematic_variables.z(lv_p, lv_q);
        z1 = kinematic_variables.z(lv_p1, lv_q);
        z2 = kinematic_variables.z(lv_p2, lv_q);

        // missing mass calculations
        Mx2 = kinematic_variables.Mx2(lv_q, lv_target, lv_p1, lv_p2);
        Mx2_1 = kinematic_variables.Mx2(lv_q, lv_target, lv_p1);
        Mx2_2 = kinematic_variables.Mx2(lv_q, lv_target, lv_p2);

        // boost to gamma*-nucleon center of mass frame
        LorentzVector lv_p_gN = new LorentzVector(lv_p);
        lv_p_gN = kinematic_variables.lv_boost_gN(lv_target, lv_q, lv_p_gN);
        LorentzVector lv_p1_gN = new LorentzVector(lv_p1);
        lv_p1_gN = kinematic_variables.lv_boost_gN(lv_target, lv_q, lv_p1_gN);
        LorentzVector lv_p2_gN = new LorentzVector(lv_p2);
        lv_p2_gN = kinematic_variables.lv_boost_gN(lv_target, lv_q, lv_p2_gN);
        LorentzVector lv_e_gN = new LorentzVector(lv_e);
        lv_e_gN = kinematic_variables.lv_boost_gN(lv_target, lv_q, lv_e_gN);
        Vector3 lv_e_gN_unit = new Vector3();
        lv_e_gN_unit.setMagThetaPhi(1, lv_e_gN.theta(), lv_e_gN.phi());
        LorentzVector lv_target_gN = new LorentzVector(lv_target);
        lv_target_gN = kinematic_variables.lv_boost_gN(lv_target, lv_q, lv_target_gN);
        LorentzVector lv_q_gN = new LorentzVector(lv_q);
        lv_q_gN = kinematic_variables.lv_boost_gN(lv_target, lv_q, lv_q_gN);
        Vector3 lv_q_gN_unit = new Vector3();
        lv_q_gN_unit.setMagThetaPhi(1, lv_q_gN.theta(), lv_q_gN.phi());
        // in gamma*-nucleon frame the z axis is along gamma* and the x axis is in the 
        // e-e' plane in the direction of e. the y axis is then the cross product of these two

        // boost to Breit infinite momentum frame
        LorentzVector lv_p_Breit = new LorentzVector(lv_p);
        lv_p_Breit.boost(BreitBoost);
        LorentzVector lv_p1_Breit = new LorentzVector(lv_p1);
        lv_p1_Breit.boost(BreitBoost);
        LorentzVector lv_p2_Breit = new LorentzVector(lv_p2);
        lv_p2_Breit.boost(BreitBoost);
        LorentzVector lv_e_Breit = new LorentzVector(lv_e);
        lv_e_Breit.boost(BreitBoost);
        Vector3 lv_e_Breit_unit = new Vector3();
        lv_e_Breit_unit.setMagThetaPhi(1, lv_e_Breit.theta(), lv_e_Breit.phi());
        LorentzVector lv_q_Breit = new LorentzVector(lv_q);
        lv_q_Breit.boost(BreitBoost);
        Vector3 lv_q_Breit_unit = new Vector3();
        lv_q_Breit_unit.setMagThetaPhi(1, lv_q_Breit.theta(), lv_q_Breit.phi());
        // note that in the Breit frame +z is antialigned to the direction of q

        // set up boost to dihadron rest frame (p1 and p2 center of mass frame)
        // just used to calculate theta
        Vector3 pCOMBoost = lv_p.boostVector();
        pCOMBoost.negative();
        LorentzVector lv_p_COM = new LorentzVector(lv_p);
        lv_p_COM.boost(pCOMBoost);
        LorentzVector lv_p1_COM = new LorentzVector(lv_p1);
        lv_p1_COM.boost(pCOMBoost);
        LorentzVector lv_p2_COM = new LorentzVector(lv_p2);
        lv_p2_COM.boost(pCOMBoost);
        theta = Math.acos(lv_p1_COM.vect().dot(lv_p.vect())
                / (lv_p1_COM.vect().mag() * lv_p.vect().mag()));

        LorentzVector lv_p_boost_rho = new LorentzVector();
        lv_p_boost_rho.setPxPyPzM(lv_p.px(), lv_p.py(), lv_p.pz(), kinematic_variables.particle_mass(113));
        Vector3 pCOMBoostrho = lv_p_boost_rho.boostVector();
        pCOMBoostrho.negative();
        LorentzVector lv_p_COM_rho = new LorentzVector(lv_p);
        lv_p_COM_rho.boost(pCOMBoostrho);
        LorentzVector lv_p1_COM_rho = new LorentzVector(lv_p1);
        lv_p1_COM_rho.boost(pCOMBoostrho);
        LorentzVector lv_p2_COM_rho = new LorentzVector(lv_p2);
        lv_p2_COM_rho.boost(pCOMBoostrho);
        p1_COM_phi = lv_p1_COM_rho.phi();
        if (p1_COM_phi < 0) {
            p1_COM_phi += 2 * Math.PI;
        }
        p1_COM_theta = lv_p1_COM_rho.theta();
        if (p1_COM_theta < 0) {
            p1_COM_theta += 2 * Math.PI;
        }
        p2_COM_phi = lv_p2_COM_rho.phi();
        if (p2_COM_phi < 0) {
            p2_COM_phi += 2 * Math.PI;
        }
        p2_COM_theta = lv_p2_COM_rho.theta();
        if (p2_COM_phi < 0) {
            p2_COM_phi += 2 * Math.PI;
        }
        COM_Delta_phi = p1_COM_phi - p2_COM_phi;
        if (COM_Delta_phi < 0) {
            COM_Delta_phi += 2 * Math.PI;
        }
        COM_Delta_theta = p1_COM_theta - p2_COM_theta;
        if (COM_Delta_theta < 0) {
            COM_Delta_theta += 2 * Math.PI;
        }
//        if (COM_Delta_theta > Math.PI) {COM_Delta_theta=Math.PI-COM_Delta_theta;}
        COM_Delta_phi = 180 / Math.PI * COM_Delta_phi;
        COM_Delta_theta = 180 / Math.PI * COM_Delta_theta;
        COM_open_angle = 180 / Math.PI * Math.acos(lv_p1_COM_rho.vect().dot(lv_p2_COM_rho.vect())
                / (lv_p1_COM_rho.vect().mag() * lv_p2_COM_rho.vect().mag()));

        Delta_phi = phi2 - phi1;
        if (Delta_phi < 0) {
            Delta_phi += 2 * Math.PI;
        }

        pT = lv_q_gN_unit.cross(lv_p_gN.vect()).mag();
        pT1 = lv_q_gN_unit.cross(lv_p1_gN.vect()).mag();
        pT2 = lv_q_gN_unit.cross(lv_p2_gN.vect()).mag();

        xF = 2 * (lv_p_gN.vect().dot(lv_q_gN.vect())) / (lv_q_gN.vect().mag() * W);
        xF1 = 2 * (lv_p1_gN.vect().dot(lv_q_gN.vect())) / (lv_q_gN.vect().mag() * W);
        xF2 = 2 * (lv_p2_gN.vect().dot(lv_q_gN.vect())) / (lv_q_gN.vect().mag() * W);

        zeta = lv_p_gN.e() / lv_target_gN.e();
        zeta1 = lv_p1_gN.e() / lv_target_gN.e();
        zeta2 = lv_p2_gN.e() / lv_target_gN.e();

        xi = kinematic_variables.Lorentz_vector_inner_product(lv_p_gN, lv_q_gN)
                / kinematic_variables.Lorentz_vector_inner_product(lv_target_gN, lv_q_gN);
        xi1 = kinematic_variables.Lorentz_vector_inner_product(lv_p1_gN, lv_q_gN)
                / kinematic_variables.Lorentz_vector_inner_product(lv_target_gN, lv_q_gN);
        xi2 = kinematic_variables.Lorentz_vector_inner_product(lv_p2_gN, lv_q_gN)
                / kinematic_variables.Lorentz_vector_inner_product(lv_target_gN, lv_q_gN);

        p_gN_pz = lv_p_gN.vect().dot(lv_q_gN.vect()) / lv_q_gN.vect().mag();
        p1_gN_pz = lv_p1_gN.vect().dot(lv_q_gN.vect()) / lv_q_gN.vect().mag();
        p2_gN_pz = lv_p2_gN.vect().dot(lv_q_gN.vect()) / lv_q_gN.vect().mag();
        p_Breit_pz = lv_p_Breit.vect().dot(lv_q_Breit.vect()) / lv_q_Breit.vect().mag();
        p1_Breit_pz = lv_p1_Breit.vect().dot(lv_q_Breit.vect()) / lv_q_Breit.vect().mag();
        p2_Breit_pz = lv_p2_Breit.vect().dot(lv_q_Breit.vect()) / lv_q_Breit.vect().mag();

        // Breit frame rapidity
        eta = -0.5 * Math.log((lv_p_Breit.e() + p_Breit_pz) / (lv_p_Breit.e() - p_Breit_pz));
        eta1 = -0.5 * Math.log((lv_p1_Breit.e() + p1_Breit_pz) / (lv_p1_Breit.e() - p1_Breit_pz));
        eta2 = -0.5 * Math.log((lv_p2_Breit.e() + p2_Breit_pz) / (lv_p2_Breit.e() - p2_Breit_pz));

        // gamma*-nucleon frame rapidity
        eta_gN = 0.5 * Math.log((lv_p_gN.e() + p_gN_pz) / (lv_p_gN.e() - p_gN_pz));
        eta1_gN = 0.5 * Math.log((lv_p1_gN.e() + p1_gN_pz) / (lv_p1_gN.e() - p1_gN_pz));
        eta2_gN = 0.5 * Math.log((lv_p2_gN.e() + p2_gN_pz) / (lv_p2_gN.e() - p2_gN_pz));

        Vector3 vecH = new Vector3();
        vecH.setMagThetaPhi(lv_p2_gN.vect().mag() / z2, lv_p2_gN.vect().theta(), lv_p2_gN.vect().phi());
        Vector3 vecR = new Vector3(lv_p1_gN.vect()); // not R yet
        vecR.setMagThetaPhi(lv_p1_gN.vect().mag() / z1, lv_p1_gN.vect().theta(), lv_p1_gN.vect().phi());
        vecR.sub(vecH); // this is really R now that the subtraction is done

        Vector3 vectRt = new Vector3();
        Vector3 R_Q = new Vector3();

        R_Q.setMagThetaPhi(vecR.dot(lv_q_gN_unit), lv_q_gN_unit.theta(), lv_q_gN_unit.phi());
        vectRt = vecR;
        vectRt.sub(R_Q);

        Vector3 vectPh = new Vector3(lv_p_gN.vect());
        Vector3 vectPh1 = new Vector3(lv_p1_gN.vect());
        Vector3 vectPh2 = new Vector3(lv_p2_gN.vect());
        Vector3 Pt_Q = new Vector3();
        Pt_Q.setMagThetaPhi(vecR.dot(lv_q_gN_unit), lv_q_gN_unit.theta(), lv_q_gN_unit.phi());
        Vector3 vectPhT = new Vector3(vectPh);
        vectPhT.sub(Pt_Q);
        Vector3 vectPhT1 = new Vector3(vectPh1);
        vectPhT1.sub(Pt_Q);
        Vector3 vectPhT2 = new Vector3(vectPh2);
        vectPhT2.sub(Pt_Q);

        Vector3 vT = new Vector3(lv_q_gN_unit.cross(lv_e_gN_unit));
        vT.unit();
        Vector3 vTR = new Vector3(lv_q_gN_unit.cross(vectRt));
        vTR.unit();
        Vector3 vTH = new Vector3(lv_q_gN_unit.cross(vectPhT));
        vTH.unit();
        Vector3 vTH1 = new Vector3(lv_q_gN_unit.cross(vectPhT1));
        vTH1.unit();
        Vector3 vTH2 = new Vector3(lv_q_gN_unit.cross(vectPhT2));
        vTH2.unit();

        double cosPhiR = vT.dot(vTR);
        double sinPhiR = lv_e_gN.vect().cross(vectRt).dot(lv_q_gN_unit);
        double cosPhiH = vT.dot(vTH);
        double sinPhiH = lv_e_gN.vect().cross(vectPhT).dot(lv_q_gN_unit);
        double cosPhiH1 = vT.dot(vTH1);
        double sinPhiH1 = lv_e_gN.vect().cross(vectPhT1).dot(lv_q_gN_unit);
        double cosPhiH2 = vT.dot(vTH2);
        double sinPhiH2 = lv_e_gN.vect().cross(vectPhT2).dot(lv_q_gN_unit);

        // scaling
        double rScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag() * lv_q_gN_unit.cross(vectRt).mag();
        sinPhiR = sinPhiR / rScale;
        double hScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag() * lv_q_gN_unit.cross(vectPh).mag();
        sinPhiH = sinPhiH / hScale;
        double h1Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag() * lv_q_gN_unit.cross(vectPh1).mag();
        sinPhiH1 = sinPhiH1 / h1Scale;
        double h2Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag() * lv_q_gN_unit.cross(vectPh2).mag();
        sinPhiH2 = sinPhiH2 / h2Scale;

        phih = Math.acos(cosPhiH);
        phi1 = Math.acos(cosPhiH1);
        phi2 = Math.acos(cosPhiH2);
        phiR = Math.acos(cosPhiR);

        if (sinPhiR < 0.0) {
            phiR = 2 * Math.PI - phiR;
        }
        if (sinPhiH < 0.0) {
            phih = 2 * Math.PI - phih;
        }
        if (sinPhiH1 < 0.0) {
            phi1 = 2 * Math.PI - phi1;
        }
        if (sinPhiH2 < 0.0) {
            phi2 = 2 * Math.PI - phi2;
        }

        Delta_phi = phi2 - phi1;
        if (Delta_phi < 0) {
            Delta_phi += 2 * Math.PI;
        }

    }

    public int get_helicity() { // -1, 0, or 1. 0 equals unassigned by EventBuilder
        if (runnum <= 5666) {
            return -1 * helicity;
        } else if (runnum >= 6616 && runnum <= 6783) {
            return -1 * helicity;
        } else if (runnum >= 6120 && runnum <= 6604) {
            return -1 * helicity;
        } else if (runnum >= 11093 && runnum <= 11283) {
            return helicity;
        } else if (runnum >= 11284 && runnum < 11300) {
            return -1 * helicity;
        } else if (runnum >= 11323 && runnum <= 11571) {
            return helicity;
        }
//        System.out.println("runnum not found, assigning helicity flip");
        return helicity;
    }

    public int get_runnum() {
        return runnum;
    }

    public int get_detector1() {
        return detector1;
    }

    public int get_detector2() {
        return detector2;
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
    
    public int e_pid_cuts() {
        return e_pid_cuts;
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
    }// returns tmin

    public double t1() {
        return Double.valueOf(Math.round(t1 * 100000)) / 100000;
    }// returns t1

    public double t2() {
        return Double.valueOf(Math.round(t2 * 100000)) / 100000;
    }// returns t2

    public double z() {
        return Double.valueOf(Math.round(z * 100000)) / 100000;
    }// returns z

    public double z1() {
        return Double.valueOf(Math.round(z1 * 100000)) / 100000;
    }// returns z1

    public double z2() {
        return Double.valueOf(Math.round(z2 * 100000)) / 100000;
    }// returns z2

    public double Mx2() {
        return Double.valueOf(Math.round(Mx2 * 100000)) / 100000;
    }// returns Mx(ep1p2)

    public double Mx2_1() {
        return Double.valueOf(Math.round(Mx2_1 * 100000)) / 100000;
    }// returns Mx(ep1)

    public double Mx2_2() {
        return Double.valueOf(Math.round(Mx2_2 * 100000)) / 100000;
    }// returns Mx(ep2)

    public double Mh() {
        return Double.valueOf(Math.round(Mh * 100000)) / 100000;
    }// returns Mh

    public double pT() {
        return Double.valueOf(Math.round(pT * 100000)) / 100000;
    }// returns pT

    public double pT1() {
        return Double.valueOf(Math.round(pT1 * 100000)) / 100000;
    }// returns pT1

    public double pT2() {
        return Double.valueOf(Math.round(pT2 * 100000)) / 100000;
    }// returns pT2

    public double pTpT() {
        return Double.valueOf(Math.round(pT1 * pT2 * 100000)) / 100000;
    } // returns ptpt

    public double xF() {
        return Double.valueOf(Math.round(xF * 100000)) / 100000;
    }// returns xF

    public double xF1() {
        return Double.valueOf(Math.round(xF1 * 100000)) / 100000;
    }// returns xF1

    public double xF2() {
        return Double.valueOf(Math.round(xF2 * 100000)) / 100000;
    }// returns xF2

    public double zeta() {
        return Double.valueOf(Math.round(zeta * 100000)) / 100000;
    }// returns zeta

    public double zeta1() {
        return Double.valueOf(Math.round(zeta1 * 100000)) / 100000;
    }// returns zeta1

    public double zeta2() {
        return Double.valueOf(Math.round(zeta2 * 100000)) / 100000;
    }// returns zeta2

    public double xi() {
        return Double.valueOf(Math.round(xi * 100000)) / 100000;
    }

    public double xi1() {
        return Double.valueOf(Math.round(xi1 * 100000)) / 100000;
    }

    public double xi2() {
        return Double.valueOf(Math.round(xi2 * 100000)) / 100000;
    }

    public double p1_Breit_pz() {
        return Double.valueOf(Math.round(p1_Breit_pz * 100000)) / 100000;
    }
    // returns p1 pz in Breit

    public double p2_Breit_pz() {
        return Double.valueOf(Math.round(p2_Breit_pz * 100000)) / 100000;
    }
    // returns p2 pz in Breit

    public double p1_gN_pz() {
        return Double.valueOf(Math.round(p1_gN_pz * 100000)) / 100000;
    } // returns p1 pz in gN

    public double p2_gN_pz() {
        return Double.valueOf(Math.round(p2_gN_pz * 100000)) / 100000;
    } // returns p2 pz in gN

    public double eta() {
        return Double.valueOf(Math.round(eta * 100000)) / 100000;
    }// returns eta in the Breit frame

    public double eta1() {
        return Double.valueOf(Math.round(eta1 * 100000)) / 100000;
    }// returns eta1 in the Breit frame

    public double eta2() {
        return Double.valueOf(Math.round(eta2 * 100000)) / 100000;
    }// returns eta2 in the Breit frame

    public double eta_gN() {
        return Double.valueOf(Math.round(eta_gN * 100000)) / 100000;
    }// returns eta_gN

    public double eta1_gN() {
        return Double.valueOf(Math.round(eta1_gN * 100000)) / 100000;
    }// returns eta1_gN

    public double eta2_gN() {
        return Double.valueOf(Math.round(eta2_gN * 100000)) / 100000;
    }// returns eta2_gN

    public double Delta_eta() {
        return Double.valueOf(Math.round((eta2 - eta1) * 100000)) / 100000;
    }// returns Delta_eta, 
    // should be Lorentz invariant
    // maybe set a check for difference of Delta_eta in Breit and gN frames in future...
    // make sure the sign is what you really want

    public double theta() {
        return Double.valueOf(Math.round(theta * 100000)) / 100000;
    }
    // returns theta, the dihadron "decay angle"

    public double phih() {
        return Double.valueOf(Math.round(phih * 100000)) / 100000;
    }// returns phih

    public double phiR() {
        return Double.valueOf(Math.round(phiR * 100000)) / 100000;
    }// returns phiR

    public double phi1() {
        return Double.valueOf(Math.round(phi1 * 100000)) / 100000;
    }
    // returns phi1 (gamma*-nucleon frame)

    public double phi2() {
        return Double.valueOf(Math.round(phi2 * 100000)) / 100000;
    }
    // returns phi2 (gamma*-nucleon frame)

    public double Delta_phi() {
        return Double.valueOf(Math.round(Delta_phi * 100000)) / 100000;
    }//returns Delta_phi (p1-p2)

    public double Depolarization_A() {
        return Double.valueOf(Math.round(Depolarization_A * 100000)) / 100000;
    }
    // returns Depolarization_A

    public double Depolarization_B() {
        return Double.valueOf(Math.round(Depolarization_B * 100000)) / 100000;
    }
    // returns Depolarization_B

    public double Depolarization_C() {
        return Double.valueOf(Math.round(Depolarization_C * 100000)) / 100000;
    }
    // returns Depolarization_C

    public double Depolarization_V() {
        return Double.valueOf(Math.round(Depolarization_V * 100000)) / 100000;
    }
    // returns Depolarization_V

    public double Depolarization_W() {
        return Double.valueOf(Math.round(Depolarization_W * 100000)) / 100000;
    }
    // returns Depolarization_W

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

    public double p_px() {
        return Double.valueOf(Math.round(p_px * 100000)) / 100000;
    }// returns dihadron lab frame px

    public double p_py() {
        return Double.valueOf(Math.round(p_py * 100000)) / 100000;
    }// returns dihadron lab frame py

    public double p_pz() {
        return Double.valueOf(Math.round(p_pz * 100000)) / 100000;
    }// returns dihadron lab frame pz

    public double p_p() {
        return Double.valueOf(Math.round(p_p * 100000)) / 100000;
    }// returns dihadron lab frame p

    public double p_e() {
        return Double.valueOf(Math.round(p_e * 100000)) / 100000;
    }// returns dihadron lab frame energy

    public double p1_px() {
        return Double.valueOf(Math.round(p1_px * 100000)) / 100000;
    }// returns hadron 1 lab frame px

    public double p1_py() {
        return Double.valueOf(Math.round(p1_py * 100000)) / 100000;
    }// returns hadron 1 lab frame py

    public double p1_pz() {
        return Double.valueOf(Math.round(p1_pz * 100000)) / 100000;
    }// returns hadron 1 lab frame pz

    public double p1_p() {
        return Double.valueOf(Math.round(p1_p * 100000)) / 100000;
    }// returns hadron 1 lab frame p

    public double p1_e() {
        return Double.valueOf(Math.round(p1_e * 100000)) / 100000;
    }// returns hadron 1 lab frame energy

    public double p1_theta() {
        return Double.valueOf(Math.round(p1_theta * 100000)) / 100000;
    } // returns p1 lab 
    // frame polar angle

    public double p1_phi() {
        return Double.valueOf(Math.round(p1_phi * 100000)) / 100000;
    } // returns p1 lab 
    // frame polar angle

    public double p2_px() {
        return Double.valueOf(Math.round(p2_px * 100000)) / 100000;
    }// returns hadron 2 lab frame px

    public double p2_py() {
        return Double.valueOf(Math.round(p2_py * 100000)) / 100000;
    }// returns hadron 2 lab frame py

    public double p2_pz() {
        return Double.valueOf(Math.round(p2_pz * 100000)) / 100000;
    }// returns hadron 2 lab frame pz

    public double p2_p() {
        return Double.valueOf(Math.round(p2_p * 100000)) / 100000;
    }// returns hadron 2 lab frame p

    public double p2_e() {
        return Double.valueOf(Math.round(p2_e * 100000)) / 100000;
    }// returns hadron 2 lab frame energy

    public double p2_theta() {
        return Double.valueOf(Math.round(p2_theta * 100000)) / 100000;
    } // returns p2 lab frame polar angle

    public double p2_phi() {
        return Double.valueOf(Math.round(p2_phi * 100000)) / 100000;
    } // returns p2 lab frame polar angle

    public double vx_e() {
        return Double.valueOf(Math.round(vx_e * 100000)) / 100000;
    }// returns electron x vertex

    public double vx_p1() {
        return Double.valueOf(Math.round(vx_p1 * 100000)) / 100000;
    }// returns p1 x vertex

    public double vx_p2() {
        return Double.valueOf(Math.round(vx_p2 * 100000)) / 100000;
    }// returns p2 x vertex

    public double vy_e() {
        return Double.valueOf(Math.round(vy_e * 100000)) / 100000;
    }// returns electron y vertex

    public double vy_p1() {
        return Double.valueOf(Math.round(vy_p1 * 100000)) / 100000;
    }// returns p1 y vertex

    public double vy_p2() {
        return Double.valueOf(Math.round(vy_p2 * 100000)) / 100000;
    }// returns p2 y vertex

    public double vz_e() {
        return Double.valueOf(Math.round(vz_e * 100000)) / 100000;
    }// returns electron z vertex

    public double vz_p1() {
        return Double.valueOf(Math.round(vz_p1 * 100000)) / 100000;
    }// returns p1 z vertex

    public double vz_p2() {
        return Double.valueOf(Math.round(vz_p2 * 100000)) / 100000;
    }// returns p2 z vertex

    public double open_angle_ep() {
        return Double.valueOf(Math.round(open_angle_ep * 100000)) / 100000;
    }

    public double open_angle_ep1() {
        return Double.valueOf(Math.round(open_angle_ep1 * 100000)) / 100000;
    }

    public double open_angle_ep2() {
        return Double.valueOf(Math.round(open_angle_ep2 * 100000)) / 100000;
    }

    public double open_angle_p1p2() {
        return Double.valueOf(Math.round(open_angle_p1p2 * 100000)) / 100000;
    }

    public double p1_COM_phi() {
        return Double.valueOf(Math.round(p1_COM_phi * 100000)) / 100000;
    }// p1 COM phi

    public double p1_COM_theta() {
        return Double.valueOf(Math.round(p1_COM_theta * 100000)) / 100000;
    }// p1_COM_theta

    public double p2_COM_phi() {
        return Double.valueOf(Math.round(p2_COM_phi * 100000)) / 100000;
    }// p2_COM_phi

    public double p2_COM_theta() {
        return Double.valueOf(Math.round(p2_COM_theta * 100000)) / 100000;
    }// p2_COM_theta

    public double COM_Delta_phi() {
        return Double.valueOf(Math.round(COM_Delta_phi * 100000)) / 100000;
    }// COM_Delta_phi

    public double COM_Delta_theta() {
        return Double.valueOf(Math.round(COM_Delta_theta * 100000)) / 100000;
    }// COM_Delta_theta

    public double COM_open_angle() {
        return Double.valueOf(Math.round(COM_open_angle * 100000)) / 100000;
    }// COM_open_angle

    public double gN_angle_p1_p2() {
        return Double.valueOf(Math.round(gN_angle_p1_p2 * 100000)) / 100000;
    }// gN_angle_p1_p2

    public double gN_angle_p1_X() {
        return Double.valueOf(Math.round(gN_angle_p1_X * 100000)) / 100000;
    }// gN_angle_p1_X

    public double gN_angle_p2_X() {
        return Double.valueOf(Math.round(gN_angle_p2_X * 100000)) / 100000;
    }// gN_angle_p1_X

    public double Emiss2() {
        return Double.valueOf(Math.round(Emiss2 * 100000)) / 100000;
    }// returns Emiss2

    public double theta_gamma_gamma() {
        return Double.valueOf(Math.round(theta_gamma_gamma * 100000)) / 100000;
    }// returns theta_gamma_gamma

    public double pTmiss() {
        return Double.valueOf(Math.round(pTmiss * 100000)) / 100000;
    }// returns pTmiss

    public int emilay1() {
        return emilay1;
    }

    public int emico1() {
        return emico1;
    }

    public int emqua1() {
        return emqua1;
    }

    public int best_PID1() {
        return best_PID1;
    }

    public float RQ1() {
        return RQ1;
    }

    public float ReQ1() {
        return ReQ1;
    }

    public float el_logl1() {
        return el_logl1;
    }

    public float pi_logl1() {
        return pi_logl1;
    }

    public float k_logl1() {
        return k_logl1;
    }

    public float pr_logl1() {
        return pr_logl1;
    }

    public float best_ch1() {
        return best_ch1;
    }

    public float best_c21() {
        return best_c21;
    }

    public float best_RL1() {
        return best_RL1;
    }

    public float best_ntot1() {
        return best_ntot1;
    }

    public int emilay2() {
        return emilay2;
    }

    public int emico2() {
        return emico2;
    }

    public int emqua2() {
        return emqua2;
    }

    public int best_PID2() {
        return best_PID2;
    }

    public float RQ2() {
        return RQ2;
    }

    public float ReQ2() {
        return ReQ2;
    }

    public float el_logl2() {
        return el_logl2;
    }

    public float pi_logl2() {
        return pi_logl2;
    }

    public float k_logl2() {
        return k_logl2;
    }

    public float pr_logl2() {
        return pr_logl2;
    }

    public float best_ch2() {
        return best_ch2;
    }

    public float best_c22() {
        return best_c22;
    }

    public float best_RL2() {
        return best_RL2;
    }

    public float best_ntot2() {
        return best_ntot2;
    }
}
