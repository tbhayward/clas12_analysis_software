package analyzers;

/**
 *
 * @author tbhayward
 */
import extended_kinematic_fitters.fiducial_cuts;
import extended_kinematic_fitters.generic_tests;
import extended_kinematic_fitters.momentum_corrections;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.clas.physics.*;

public class TwoParticles {

    protected byte helicity;
    protected int runnum;

    protected int fiducial_status = -1;
    protected int detector = -1;

    protected int num_elec, num_piplus, num_piminus, num_kplus, num_kminus, num_protons, num_particles;
    protected int num_pos, num_neg, num_neutrals;
    protected int num_positrons, num_antiprotons;

    protected double open_angle;
    protected double Q2, W, gamma, nu, x, y, z, t, tmin;
    protected double Mx, Mx2;
    protected double Mh, pT, xF, zeta, xi;
    protected double eta, eta_gN;

    protected double phi;

    protected double Depolarization_A;
    protected double Depolarization_B;
    protected double Depolarization_C;
    protected double Depolarization_V;
    protected double Depolarization_W;

    protected double e_px, e_py, e_pz, e_p, e_e, e_theta, e_phi; // electron kinematics
    protected double p_px, p_py, p_pz, p_p, p_e, p_theta, p_phi; // hadron kinematics
    protected double vz_e, vz_p;

    protected double p_Breit_pz, p_gN_pz;

    protected int RICH_pid;
    protected double chi2pid, beta, RQ_prob, el_prob, pi_prob, k_prob, pr_prob;

    // --- NEW: one-shot inverse-ISR photon to subtract from q ---
    private static boolean useInverseISRNext = false;
    private static double nextEgammaGeV = 0.0;
    private static double nextThetaRad = 0.0;
    private static double nextPhiRad   = 0.0;

    /** Arm a one-shot inverse-ISR correction for the next constructed TwoParticles. */
    public static void setNextInverseISRPhoton(double EgammaGeV, double thetaRad, double phiRad) {
        useInverseISRNext = true;
        nextEgammaGeV = EgammaGeV;
        nextThetaRad = thetaRad;
        nextPhiRad = phiRad;
    }

    public static boolean channel_test(TwoParticles variables) {
        if (variables.helicity == 0 && variables.runnum != 11) {
            return false;
        } else if (variables.Q2() < 1.00) {
            return false;
        } else if (variables.W() < 2) {
            return false;
        } else if (variables.y() > 0.80) {
            return false;
        } 
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

    public TwoParticles(DataEvent event, PhysicsEvent recEvent, int pPID, int pIndex, double Eb) {
        kinematic_variables kinematic_variables = new kinematic_variables();

        // banks
        HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
        HipoDataBank configBank = (HipoDataBank) event.getBank("RUN::config");
        HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle");
        HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
        HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");

        helicity = eventBank.getByte("helicity", 0);
        runnum = configBank.getInt("run", 0);

        num_elec = recEvent.countByPid(11);
        num_positrons = recEvent.countByPid(-11);
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

        int p_rec_index = getIndex(rec_Bank, pPID, pIndex);
        if (generic_tests.forward_tagger_cut(p_rec_index, rec_Bank)) {
            detector = 0;
        } else if (generic_tests.forward_detector_cut(p_rec_index, rec_Bank)) {
            detector = 1;
        } else if (generic_tests.central_detector_cut(p_rec_index, rec_Bank)) {
            detector = 2;
        }

        boolean passesForwardDetector_1 = generic_tests.forward_detector_cut(p_rec_index, rec_Bank)
                ? fiducial_cuts.dc_fiducial_cut(p_rec_index, rec_Bank, traj_Bank, configBank) : true;
        boolean passesCentralDetector_1 = generic_tests.central_detector_cut(p_rec_index, rec_Bank)
                ? fiducial_cuts.cvt_fiducial_cut(p_rec_index, rec_Bank, traj_Bank, 1) : true;
        boolean passesForwardTagger_1 = generic_tests.forward_tagger_cut(p_rec_index, rec_Bank)
                ? fiducial_cuts.forward_tagger_fiducial_cut(p_rec_index, rec_Bank, cal_Bank) : true;
        boolean p_fiducial_check = passesForwardTagger_1 && passesForwardDetector_1 && passesCentralDetector_1;

        fiducial_status = 0;
        if (electron_pcal_fiducial) fiducial_status += 1;
        if (electron_fd_fiducial)   fiducial_status += 10;
        if (p_fiducial_check)       fiducial_status += 100;

        // ==== Lorentz vectors ====
        // beam electron (no tilt)
        LorentzVector lv_beam = new LorentzVector();
        double me = kinematic_variables.particle_mass(11);
        double pBeam = Math.sqrt(Math.max(0.0, Eb * Eb - me * me));
        lv_beam.setPxPyPzM(0.0, 0.0, pBeam, me);

        // scattered electron
        String electron_index = "[11,0]";
        Particle scattered_electron = recEvent.getParticle(electron_index);
        LorentzVector lv_e = new LorentzVector();
        lv_e.setPxPyPzM(scattered_electron.px(), scattered_electron.py(),
                        scattered_electron.pz(), kinematic_variables.particle_mass(11));

        String pIndex_string = "[" + pPID + "," + pIndex + "]";
        Particle hadron = recEvent.getParticle(pIndex_string);

        vz_e = scattered_electron.vz();
        vz_p = hadron.vz();

        LorentzVector lv_p = new LorentzVector();
        lv_p.setPxPyPzM(hadron.px(), hadron.py(), hadron.pz(), hadron.mass());

        // --- Build baseline q = k - k' ---
        LorentzVector lv_q = new LorentzVector(lv_beam);
        lv_q.sub(lv_e);

        // --- Optional inverse-ISR: subtract R from q (consume one-shot) ---
        if (useInverseISRNext) {
            useInverseISRNext = false; // consume

            // guard: Egamma <= 0.999 * nu_baseline
            double nu_baseline = lv_q.e();
            double Egamma_eff = Math.min(nextEgammaGeV, Math.max(0.0, 0.999 * nu_baseline));

            if (Egamma_eff > 0.0) {
                Vector3 nhat = new Vector3();
                nhat.setMagThetaPhi(1.0, nextThetaRad, nextPhiRad);

                LorentzVector lv_R = new LorentzVector();
                lv_R.setPxPyPzM(nhat.x() * Egamma_eff,
                                nhat.y() * Egamma_eff,
                                nhat.z() * Egamma_eff,
                                0.0);
                lv_q.sub(lv_R);
            }
        }

        // target and missing mass (use corrected q)
        LorentzVector lv_target = new LorentzVector();
        momentum_corrections momentum_corrections = new momentum_corrections();
        lv_target.setPxPyPzM(0, 0, 0, kinematic_variables.particle_mass(2212));

        Mx  = kinematic_variables.Mx (lv_q, lv_target, lv_p);
        Mx2 = kinematic_variables.Mx2(lv_q, lv_target, lv_p);
        
//        /* TOGGLE ON OR OFF IF FERMI MOTION DESIRED */
//        // Simulate Fermi motion
//        org.jlab.clas.physics.Vector3 fermiP = momentum_corrections.sampleFermiMomentum(Mx2);
//        lv_target.setPxPyPzM(fermiP.x(),fermiP.y(),fermiP.z(),kinematic_variables.particle_mass(2212));
//        Mx  = kinematic_variables.Mx (lv_q, lv_target, lv_p);
//        Mx2 = kinematic_variables.Mx2(lv_q, lv_target, lv_p);


        // electron kinematics (for output)
        e_px = lv_e.px();
        e_py = lv_e.py();
        e_pz = lv_e.pz();
        e_p  = lv_e.p();
        e_e  = lv_e.e();
        e_theta = scattered_electron.theta();
        e_phi   = scattered_electron.phi();
        if (e_phi < 0) e_phi = 2 * Math.PI + e_phi;

        // ---- DIS variables (fully invariant, supports moving target) ----

        Q2 = kinematic_variables.Q2(lv_q);
        // x_B = Q2 / (2 p·q)
        x  = kinematic_variables.x(lv_q, lv_target);
        // W = sqrt( (p + q)^2 )
        W  = kinematic_variables.W(lv_q, lv_target);
        // y = (p·q) / (p·k)
        y  = kinematic_variables.y(lv_beam, lv_q, lv_target);
        // gamma from invariant x and Q2
        gamma = kinematic_variables.gamma(Q2, x);

        // Depolarization from corrected (gamma, y)
        Depolarization_A = kinematic_variables.Depolarization_A(gamma, y);
        Depolarization_B = kinematic_variables.Depolarization_B(gamma, y);
        Depolarization_C = kinematic_variables.Depolarization_C(gamma, y);
        Depolarization_V = kinematic_variables.Depolarization_V(gamma, y);
        Depolarization_W = kinematic_variables.Depolarization_W(gamma, y);

        // g*-N COM frame (from corrected q)
        LorentzVector gN = new LorentzVector(lv_q);
        gN.add(lv_target);
        Vector3 gNBoost = gN.boostVector();
        gNBoost.negative();

        // Breit frame (from corrected q)
        LorentzVector Breit = new LorentzVector(lv_q);
        LorentzVector Breit_target = new LorentzVector();
        Breit_target.setPxPyPzM(0, 0, 0, 2 * x * kinematic_variables.particle_mass(2212));
        Breit.add(Breit_target);
        Vector3 BreitBoost = Breit.boostVector();
        BreitBoost.negative();

        open_angle = kinematic_variables.open_angle(lv_e, lv_p);
//        t = kinematic_variables.t(lv_p.p(), lv_p.theta());
        t = kinematic_variables.t(lv_target, lv_p);
        tmin = kinematic_variables.tmin(x, Q2);

        // hadron kinematics
        p_px = lv_p.px();
        p_py = lv_p.py();
        p_pz = lv_p.pz();
        p_p  = lv_p.p();
        p_e  = hadron.e();
        p_theta = hadron.theta();
        p_phi   = hadron.phi();
        if (p_phi < 0) p_phi = 2 * Math.PI + p_phi;

        // z from corrected q
        z = kinematic_variables.z(lv_p, lv_q);

        // boosts with corrected q
        LorentzVector lv_p_gN = new LorentzVector(lv_p);
        lv_p_gN = kinematic_variables.lv_boost_gN(lv_target, lv_q, lv_p_gN);
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

        // Breit
        LorentzVector lv_p_Breit = new LorentzVector(lv_p);
        lv_p_Breit.boost(BreitBoost);
        LorentzVector lv_e_Breit = new LorentzVector(lv_e);
        lv_e_Breit.boost(BreitBoost);
        Vector3 lv_e_Breit_unit = new Vector3();
        lv_e_Breit_unit.setMagThetaPhi(1, lv_e_Breit.theta(), lv_e_Breit.phi());
        LorentzVector lv_q_Breit = new LorentzVector(lv_q);
        lv_q_Breit.boost(BreitBoost);
        Vector3 lv_q_Breit_unit = new Vector3();
        lv_q_Breit_unit.setMagThetaPhi(1, lv_q_Breit.theta(), lv_q_Breit.phi());

        pT = lv_q_gN_unit.cross(lv_p_gN.vect()).mag();
        xF = 2 * (lv_p_gN.vect().dot(lv_q_gN.vect())) / (lv_q_gN.vect().mag() * W);

        zeta = lv_p_gN.e() / lv_target_gN.e();

        LightConeKinematics lck = new LightConeKinematics();
        xi = lck.xi_h(lv_p_gN, lv_q_gN, lv_target_gN);

        p_gN_pz = lv_p_gN.vect().dot(lv_q_gN.vect()) / lv_q_gN.vect().mag();
        p_Breit_pz = lv_p_Breit.vect().dot(lv_q_Breit.vect()) / lv_q_Breit.vect().mag();

        // rapidities
        eta = -0.5 * Math.log((lv_p_Breit.e() + p_Breit_pz) / (lv_p_Breit.e() - p_Breit_pz));
        eta_gN = 0.5 * Math.log((lv_p_gN.e() + p_gN_pz) / (lv_p_gN.e() - p_gN_pz));

        // Trento phi using corrected q
        Vector3 vecH = new Vector3();
        vecH.setMagThetaPhi(lv_p_gN.vect().mag() / z, lv_p_gN.vect().theta(), lv_p_gN.vect().phi());
        Vector3 vecR = new Vector3(vecH);
        vecR.negative();

        double dotProductRQ = vecR.dot(lv_q_gN_unit);
        Vector3 R_Q = new Vector3(lv_q_gN_unit);
        R_Q.setMagThetaPhi(dotProductRQ, lv_q_gN_unit.theta(), lv_q_gN_unit.phi());

        Vector3 vectPhT = new Vector3(lv_p_gN.vect());
        vectPhT.sub(R_Q);

        Vector3 vT = lv_q_gN_unit.cross(lv_e_gN_unit);
        vT.unit();
        Vector3 vTH = lv_q_gN_unit.cross(vectPhT);
        vTH.unit();

        double cosPhiH = vT.dot(vTH);
        double sinPhiH = lv_e_gN.vect().cross(vectPhT).dot(lv_q_gN_unit);

        double hScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag() * lv_q_gN_unit.cross(vecH).mag();
        sinPhiH = sinPhiH / hScale;

        phi = Math.acos(cosPhiH);
        if (sinPhiH < 0.0) {
            phi = 2 * Math.PI - phi;
        }
        // end
    }

    public int get_helicity() {
        if (runnum >= 4326 && runnum <= 5666) return -1 * helicity;
        else if (runnum >= 6616 && runnum <= 6783) return -1 * helicity;
        else if (runnum >= 6120 && runnum <= 6604) return -1 * helicity;
        else if (runnum >= 11093 && runnum <= 11283) return helicity;
        else if (runnum >= 11284 && runnum < 11300) return -1 * helicity;
        else if (runnum >= 11323 && runnum < 11571) return helicity;
        return helicity;
    }

    public int get_runnum() { return runnum; }
    public int get_detector() { return detector; }
    public int get_num_pos() { return num_pos; }
    public int get_num_neg() { return num_neg; }
    public int get_num_neutrals() { return num_neutrals; }
    public int get_fiducial_status() { return fiducial_status; }

    public int num_elec() { return num_elec; }
    public int num_piplus() { return num_piplus; }
    public int num_piminus() { return num_piminus; }
    public int num_kplus() { return num_kplus; }
    public int num_kminus() { return num_kminus; }
    public int num_protons() { return num_protons; }

    public double Q2() { return ((int) (Q2 * 100000)) / 100000.0; }
    public double W()  { return ((int) (W  * 100000)) / 100000.0; }
    public double gamma() { return ((int) (gamma * 100000)) / 100000.0; }
    public double nu() { return ((int) (nu * 100000)) / 100000.0; }
    public double x()  { return ((int) (x  * 100000)) / 100000.0; }
    public double y()  { return ((int) (y  * 100000)) / 100000.0; }

    public double t()    { return Double.valueOf(Math.round(t    * 100000)) / 100000; }
    public double tmin() { return Double.valueOf(Math.round(tmin * 100000)) / 100000; }

    public double z()   { return ((int) (z   * 100000)) / 100000.0; }
    public double Mx()  { return ((int) (Mx  * 100000)) / 100000.0; }
    public double Mx2() { return ((int) (Mx2 * 100000)) / 100000.0; }
    public double pT()  { return ((int) (pT  * 100000)) / 100000.0; }
    public double xF()  { return ((int) (xF  * 100000)) / 100000.0; }
    public double zeta(){ return ((int) (zeta* 100000)) / 100000.0; }
    public double xi()  { return ((int) (xi  * 100000)) / 100000.0; }

    public double p_Breit_pz() { return ((int) (p_Breit_pz * 100000)) / 100000.0; }
    public double p_gN_pz()    { return ((int) (p_gN_pz    * 100000)) / 100000.0; }
    public double eta()        { return ((int) (eta        * 100000)) / 100000.0; }
    public double eta_gN()     { return ((int) (eta_gN     * 100000)) / 100000.0; }
    public double phi()        { return ((int) (phi        * 100000)) / 100000.0; }

    public double Depolarization_A() { return ((int) (Depolarization_A * 100000)) / 100000.0; }
    public double Depolarization_B() { return ((int) (Depolarization_B * 100000)) / 100000.0; }
    public double Depolarization_C() { return ((int) (Depolarization_C * 100000)) / 100000.0; }
    public double Depolarization_V() { return ((int) (Depolarization_V * 100000)) / 100000.0; }
    public double Depolarization_W() { return ((int) (Depolarization_W * 100000)) / 100000.0; }

    public double e_px() { return ((int) (e_px * 100000)) / 100000.0; }
    public double e_py() { return ((int) (e_py * 100000)) / 100000.0; }
    public double e_pz() { return ((int) (e_pz * 100000)) / 100000.0; }
    public double e_p()  { return ((int) (e_p  * 100000)) / 100000.0; }
    public double e_e()  { return ((int) (e_e  * 100000)) / 100000.0; }
    public double e_theta() { return ((int) (e_theta * 100000)) / 100000.0; }
    public double e_phi()   { return ((int) (e_phi   * 100000)) / 100000.0; }

    public double p_px() { return ((int) (p_px * 100000)) / 100000.0; }
    public double p_py() { return ((int) (p_py * 100000)) / 100000.0; }
    public double p_pz() { return ((int) (p_pz * 100000)) / 100000.0; }
    public double p_p()  { return ((int) (p_p  * 100000)) / 100000.0; }
    public double p_e()  { return ((int) (p_e  * 100000)) / 100000.0; }
    public double p_theta() { return ((int) (p_theta * 100000)) / 100000.0; }
    public double p_phi()   { return ((int) (p_phi   * 100000)) / 100000.0; }

    public double vz_e() { return ((int) (vz_e * 100000)) / 100000.0; }
    public double vz_p() { return ((int) (vz_p * 100000)) / 100000.0; }
    public double open_angle() { return ((int) (open_angle * 100000)) / 100000.0; }
}