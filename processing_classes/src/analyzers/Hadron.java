/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analyzers;

/**
 *
 * @author tbhayward
 */

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.clas.physics.*;

public class Hadron {
    
    protected byte helicity;
    protected int runnum;
    
    protected double test;
    
    protected int num_elec, num_piplus, num_piminus, num_kplus, num_kminus, num_protons, num_particles;
    
    // labels are unnumbered if they refer to the dihadron (perhaps a meson) and numbered for individual
    // hadrons. Convention is ordered by mass, then charge. For example in pi+pi- pi+ is hadron 1
    // in proton+pi+ the proton is p1, in k+pi- the kaon is p1.
    protected double Q2, W, gamma, nu, x, y, z;
    protected double Mx, Mx2; 
    protected double Mh, pT, xF, zeta;
    protected double eta, eta_gN;
    // eta is the rapidity, preferred by theorists in the Breit frame (e.g. eta1 is in Breit) 
    // eta_gN is the rapidity in the gamma*-nucleon COM frame
    // the difference between two rapidities is Lorentz invariant, i.e.
    // eta1-eta2 = eta1_COM - eta2_COM

    protected double phi;
    
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
    protected double p_px, p_py, p_pz, p_p, p_e, p_theta, p_phi; // hadron kinematics
    protected double vz_e, vz_p;
    
    protected double p_Breit_pz, p_gN_pz;
    
    protected int RICH_pid;
    protected double chi2pid, beta, RQ_prob, el_prob, pi_prob, k_prob, pr_prob;
    
    public static boolean channel_test(Hadron variables) {
//        if (variables.helicity==0){ 
//            System.out.println("You're returning false because helicity = 0. Is this data or MC?");
//            return false; }
        if (variables.Q2()<1) { return false; } 
        if (variables.W()<2) { return false; } 
//        if (variables.xF()<0.0) { return false; } 
        else if (variables.y()>0.75) { return false; } 
//        else if (variables.z()<0.2 || variables.z() > 0.95) { return false; } 
//        else if (variables.p_p()<1.25) { return false; } 
//        else if (variables.Mx()<1.4) { return false; } 
        
        
        else if (variables.e_p > 12 || variables.p_p > 12 || variables.Q2 > 12 || variables.W > 6 || 
                variables.z > 2 || variables.zeta > 2 || 
                variables.Mx < - 4 || variables.pT > 2 || variables.xF < -2 || variables.xF > 2 || 
                variables.Depolarization_A > 5 || variables.Depolarization_B > 5 ||
                variables.Depolarization_C > 5 || variables.Depolarization_V > 5 ||
                variables.Depolarization_W > 5) 
                { return false; }
	return true;
    }
    
    public Hadron(DataEvent event, PhysicsEvent recEvent, int pPID, int pIndex, double Eb) {
        // provide the PDG PID of the two hadrons
        
        // load banks
        HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
        HipoDataBank configBank = (HipoDataBank) event.getBank("RUN::config");
        
        helicity = eventBank.getByte("helicity", 0);
//        helicity = eventBank.getByte("helicityRaw", 0);
        runnum = configBank.getInt("run",0); // used for beam energy and polarization
//    
//        // default beam energy set to rga fall 2018
//        double Eb = 10.6041; // rga fall 2018
//        if (runnum >= 5032 && runnum <= 5666) { Eb = 10.6041; } // RGA Fall 2018
//        else if (runnum >= 5875 && runnum <= 5995) { Eb = 6.535; } // RGK Fall 2018
//        else if (runnum >= 5674 && runnum <= 5870) { Eb = 7.546; } // RGK Fall 2018
//        else if (runnum >= 6616 && runnum <= 6783) { Eb = 10.1998; } // RGA Spring 2019
//        else if (runnum >= 6120 && runnum <= 6399) { Eb = 10.5986; } // RGB Spring 2019
//        else if (runnum >= 6409 && runnum <= 6604) { Eb = 10.1998; } // RGB Spring 2019
//        else if (runnum >= 11093 && runnum <= 11283) { Eb = 10.4096; } // RGB Fall 2019
//        else if (runnum >= 11284 && runnum <= 11300) { Eb = 4.17179; } // RGB Fall 2019
//        else if (runnum >= 11323 && runnum <= 11571) { Eb = 10.3894; } // RGB Spring 2020
//        else if (runnum >= 16290 && runnum <= 16353) { Eb = 10.5473; } // RGC preliminary summer 2022
        
        num_elec = recEvent.countByPid(11); // returns number of electrons
	num_piplus = recEvent.countByPid(211); 
	num_piminus = recEvent.countByPid(-211);
	num_kplus = recEvent.countByPid(321);
	num_kminus = recEvent.countByPid(-321);
        num_protons = recEvent.countByPid(2212);
        num_particles = num_elec+num_piplus+num_piminus+num_kplus+num_kminus+num_protons;
        
        // Set up Lorentz vectors
        // beam electron
        LorentzVector lv_beam = new LorentzVector();
	lv_beam.setPxPyPzM(0, 0, Math.pow(Eb*Eb-particle_mass(11)*particle_mass(11),0.5), 
                particle_mass(11));
        LorentzVector lv_target = new LorentzVector();
        // target, proton for RGA... what mass to use for RGB (deuterium target)?
	lv_target.setPxPyPzM(0, 0, 0, particle_mass(2212));
        // pull from rec banks for outgoing particles
        // electron
        String electron_index = "[11,0]"; // highest p, kinematic fitter should require FD etc
	Particle scattered_electron = recEvent.getParticle(electron_index); //
        LorentzVector lv_e = new LorentzVector();
	lv_e.setPxPyPzM(scattered_electron.px(), scattered_electron.py(), 
            scattered_electron.pz(), particle_mass(11));
        // hadrons set up below (to allow for iteration over more than two hadrons in an event)
        
        // kinematics of electron
        e_px = lv_e.px(); e_py = lv_e.py(); e_pz = lv_e.pz(); e_p = lv_e.p(); e_e = lv_e.e(); 
        e_theta = scattered_electron.theta();
        e_phi = scattered_electron.phi();
        if (e_phi < 0) { e_phi = 2*Math.PI + e_phi; }
                
        // DIS variables
        LorentzVector lv_q = new LorentzVector(lv_beam); lv_q.sub(lv_e);
	Q2 = -lv_q.mass2();
	nu = lv_beam.e()-lv_e.e();
	x  = Q2 / (2 * particle_mass(2212) * nu);
	W  = Math.pow(Math.pow(particle_mass(2212),2)+2*particle_mass(2212)*nu - Q2, 0.5);
	y = nu/lv_beam.e();
        gamma = 2*particle_mass(2212)*x/Math.pow(Q2, 0.5);
        
        // Depolarization variables
        Depolarization_A = 1/(1+gamma*gamma)*(1-y+y*y/2+y*y*gamma*gamma/4);
        Depolarization_B = 1/(1+gamma*gamma)*(1-y-y*y*gamma*gamma/4);
        Depolarization_C = (y/Math.pow(1+gamma*gamma, 0.5))*(1-y/2);
        Depolarization_V = (2-y)/(1+gamma*gamma)*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
        Depolarization_W = y/(Math.pow(1+gamma*gamma, 0.5))*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
    
        // set up boost to gamma*-nucleon center of mass frame
        LorentzVector gN = new LorentzVector(lv_q);
	gN.add(lv_target);
	Vector3 gNBoost = gN.boostVector();
	gNBoost.negative();
        
        // set up boost to Breit frame, this needs to be cross checked
        LorentzVector Breit = new LorentzVector(lv_q);
        LorentzVector Breit_target = new LorentzVector();
        Breit_target.setPxPyPzM(0, 0, 0, 2*x*particle_mass(2212));
        Breit.add(Breit_target);
        Vector3 BreitBoost = Breit.boostVector();
        BreitBoost.negative();
        
        // set up hadrons 
        String pIndex_string = "["+pPID+","+pIndex+"]";
        Particle hadron = recEvent.getParticle(pIndex_string);
        
        vz_e = scattered_electron.vz();
        vz_p = hadron.vz();
        
        LorentzVector lv_p = new LorentzVector();
        lv_p.setPxPyPzM(hadron.px(), hadron.py(), hadron.pz(), hadron.mass());
        
        // kinematics of hadrons
        p_px = lv_p.px(); p_py = lv_p.py(); p_pz = lv_p.pz(); p_p = lv_p.p(); p_e = hadron.e(); 
        p_theta = hadron.theta();
        p_phi = hadron.phi();
        if (p_phi < 0) { p_phi = 2*Math.PI + p_phi; }
    
        z = lv_p.e()/lv_q.e();
        
        // missing mass calculations
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target); lv_Mx.sub(lv_p); 
        Mx = lv_Mx.mass();
        Mx2 = lv_Mx.mass2(); // missing mass squared
        
        // boost to gamma*-nucleon center of mass frame
        LorentzVector lv_p_gN = new LorentzVector(lv_p); lv_p_gN.boost(gNBoost);
        LorentzVector lv_e_gN = new LorentzVector(lv_e); lv_e_gN.boost(gNBoost);
        LorentzVector lv_target_gN = new LorentzVector(lv_target); lv_target_gN.boost(gNBoost);
        Vector3 lv_e_gN_unit = new Vector3();
        lv_e_gN_unit.setMagThetaPhi(1, lv_e_gN.theta(), lv_e_gN.phi());
        LorentzVector lv_q_gN = new LorentzVector(lv_q); lv_q_gN.boost(gNBoost);
        Vector3 lv_q_gN_unit = new Vector3();
        lv_q_gN_unit.setMagThetaPhi(1, lv_q_gN.theta(), lv_q_gN.phi());
        Vector3 lv_y_gN = new Vector3(); // I guess this isn't really a lorenzvector so the lv name is wrong 
        lv_y_gN = lv_q_gN.vect().cross(lv_e_gN.vect()); 
        // in gamma*-nucleon frame the z axis is along gamma* and the x axis is in the 
        // e-e' plane in the direction of e. the y axis is then the cross product of these two
        
        // boost to Breit infinite momentum frame
        LorentzVector lv_p_Breit = new LorentzVector(lv_p); lv_p_Breit.boost(BreitBoost);
        LorentzVector lv_e_Breit = new LorentzVector(lv_e); lv_e_Breit.boost(BreitBoost);
        Vector3 lv_e_Breit_unit = new Vector3();
        lv_e_Breit_unit.setMagThetaPhi(1, lv_e_Breit.theta(), lv_e_Breit.phi());
        LorentzVector lv_q_Breit = new LorentzVector(lv_q); lv_q_Breit.boost(BreitBoost);
        Vector3 lv_q_Breit_unit = new Vector3();
        lv_q_Breit_unit.setMagThetaPhi(1, lv_q_Breit.theta(), lv_q_Breit.phi());
        // note that in the Breit frame +z is antialigned to the direction of q
  
        
        pT = lv_q_gN_unit.cross(lv_p_gN.vect()).mag();
    
        xF =  2*(lv_p_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        
        zeta = lv_p_gN.e()/lv_target_gN.e(); // only really applicable when p1 is a proton
    
        p_gN_pz = lv_p_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p_Breit_pz = lv_p_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        
        // Breit frame rapidity
        eta = -0.5*Math.log((lv_p_Breit.e()+p_Breit_pz) /  (lv_p_Breit.e()-p_Breit_pz));
        
        // gamma*-nucleon frame rapidity
        eta_gN = 0.5*Math.log((lv_p_gN.e()+p_gN_pz) / (lv_p_gN.e()-p_gN_pz));
        
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

        // scaling
        double hScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag() * lv_q_gN_unit.cross(vecH).mag();
        sinPhiH = sinPhiH / hScale;

        phi = Math.acos(cosPhiH);

        if (sinPhiH < 0.0) {
            phi = 2 * Math.PI - phi;
        }

        
        // see trento conventions: https://arxiv.org/pdf/hep-ph/0410050.pdf        
    }
    
    
    public int get_helicity() { // -1, 0, or 1. 0 equals unassigned by EventBuilder
        if (runnum <= 5666) {
            return -1*helicity;
        } else if ( runnum >= 6616 && runnum <= 6783) {
            return -1*helicity;
        } else if ( runnum >= 6120 && runnum <= 6604) { 
            return -1*helicity;
        } else if ( runnum >= 11093 && runnum <= 11283) {
            return helicity;
        } else if ( runnum >= 11284 && runnum < 11300) {
            return -1*helicity;
        } else if ( runnum >= 11323 && runnum < 11571) {
            return helicity;
        }
        return helicity;
    }
    
    public int get_runnum() { return runnum; }
    public int num_elec() { return num_elec; }
    public int num_piplus() { return num_piplus; }
    public int num_piminus() { return num_piminus; }
    public int num_kplus() { return num_kplus; }
    public int num_kminus() { return num_kminus; }
    public int num_protons() { return num_protons; }
    public double test() { return test; }

    public double Q2() { return ((int) (Q2 * 100000)) / 100000.0; }
    public double W() { return ((int) (W * 100000)) / 100000.0; }
    public double gamma() { return ((int) (gamma * 100000)) / 100000.0; }
    public double nu() { return ((int) (nu * 100000)) / 100000.0; }
    public double x() { return ((int) (x * 100000)) / 100000.0; }
    public double y() { return ((int) (y * 100000)) / 100000.0; }
    public double z() { return ((int) (z * 100000)) / 100000.0; }
    public double Mx() { return ((int) (Mx * 100000)) / 100000.0; }
    public double Mx2() { return ((int) (Mx2 * 100000)) / 100000.0; }
    public double pT() { return ((int) (pT * 100000)) / 100000.0; }
    public double xF() { return ((int) (xF * 100000)) / 100000.0; }
    public double zeta() { return ((int) (zeta * 100000)) / 100000.0; }
    public double p_Breit_pz() { return ((int) (p_Breit_pz * 100000)) / 100000.0; }
    public double p_gN_pz() { return ((int) (p_gN_pz * 100000)) / 100000.0; }
    public double eta() { return ((int) (eta * 100000)) / 100000.0; }
    public double eta_gN() { return ((int) (eta_gN * 100000)) / 100000.0; }
    public double phi() { return ((int) (phi * 100000)) / 100000.0; }
    public double Depolarization_A() { return ((int) (Depolarization_A * 100000)) / 100000.0; }
    public double Depolarization_B() { return ((int) (Depolarization_B * 100000)) / 100000.0; }
    public double Depolarization_C() { return ((int) (Depolarization_C * 100000)) / 100000.0; }
    public double Depolarization_V() { return ((int) (Depolarization_V * 100000)) / 100000.0; }
    public double Depolarization_W() { return ((int) (Depolarization_W * 100000)) / 100000.0; }
    public double e_px() { return ((int) (e_px * 100000)) / 100000.0; }
    public double e_py() { return ((int) (e_py * 100000)) / 100000.0; }
    public double e_pz() { return ((int) (e_pz * 100000)) / 100000.0; }
    public double e_p() { return ((int) (e_p * 100000)) / 100000.0; }
    public double e_e() { return ((int) (e_e * 100000)) / 100000.0; }
    public double e_theta() { return ((int) (e_theta * 100000)) / 100000.0; }
    public double e_phi() { return ((int) (e_phi * 100000)) / 100000.0; }
    public double p_px() { return ((int) (p_px * 100000)) / 100000.0; }
    public double p_py() { return ((int) (p_py * 100000)) / 100000.0; }
    public double p_pz() { return ((int) (p_pz * 100000)) / 100000.0; }
    public double p_p() { return ((int) (p_p * 100000)) / 100000.0; }
    public double p_e() { return ((int) (p_e * 100000)) / 100000.0; }
    public double p_theta() { return ((int) (p_theta * 100000)) / 100000.0; }
    public double p_phi() { return ((int) (p_phi * 100000)) / 100000.0; }
    public double vz_e() { return ((int) (vz_e * 100000)) / 100000.0; }
    public double vz_p() { return ((int) (vz_p * 100000)) / 100000.0; }
    
    
    private static double particle_mass (int pid) {
	if (pid==11||pid==-11) { // electron is pid=11, positron is pid=-11
            return 0.0005109989461;
        } else if (pid == 111) {
            return 0.1349768;
	} else if (pid==211||pid==-211) { // pions
            return 0.139570;
	} else if (pid==321||pid==-321) { // kaons
            return 0.493677;
	} else if (pid==2212||pid==-2212) { // protons
            return 0.938272;
//            return 1.875;
	} else if (pid==113) { // rho0
            return 0.7754;
        }
            return -1;
    }
}
