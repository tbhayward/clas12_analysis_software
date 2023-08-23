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


public class Dihadrons {
    
    protected byte helicity;
    protected int runnum;
    
    protected double test;
    
    protected int num_elec, num_piplus, num_piminus, num_kplus, num_kminus, num_protons, num_particles;
    
    // labels are unnumbered if they refer to the dihadron (perhaps a meson) and numbered for individual
    // hadrons. Convention is ordered by mass, then charge. For example in pi+pi- pi+ is hadron 1
    // in proton+pi+ the proton is p1, in k+pi- the kaon is p1.
    protected double Q2, W, gamma, nu, x, y, z, z1, z2;
    protected double Mx, Mx1, Mx2, Mx3; // Mx is the Mx(ep1p2), Mx1 is Mx(ep1), Mx2 is Mx(ep2), Mx3 is Mx(e)
    protected double Mh, pT, pT1, pT2, xF, xF1, xF2, zeta, zeta1, zeta2;
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
    protected double p1_px_unc, p1_py_unc, p1_pz_unc, p1_p_unc, mass2, mass2_unc;
    protected double p2_px, p2_py, p2_pz, p2_p, p2_e, p2_theta, p2_phi; // p2 kinematics
    protected double vx_e, vx_p1, vx_p2, vy_e, vy_p1, vy_p2, vz_e, vz_p1, vz_p2;
    
    protected double p_Breit_pz, p1_Breit_pz, p2_Breit_pz, p_gN_pz, p1_gN_pz, p2_gN_pz;
    
    protected float e_nphe, p1_nphe, p2_nphe; // number of photoelectrons in cherenkov counter
    protected float e_chi2pid, p1_chi2pid, p2_chi2pid; // chi2pid of the CLAS12 EventBuilder
    // defined as the sampling fraction for electron candidates and the distance away from the 
    // mean beta(p) value from the TOF for hadron candidates
    
    protected int p1_RICH_pid;
    protected double p1_beta, p1_RQ_prob, p1_el_prob, p1_pi_prob, p1_k_prob, p1_pr_prob;
    protected double p2_beta, p2_RQ_prob, p2_el_prob, p2_pi_prob, p2_k_prob, p2_pr_prob;
    
    protected double p1_COM_phi, p1_COM_theta, p2_COM_phi, p2_COM_theta, COM_Delta_phi, COM_Delta_theta;
    protected double COM_open_angle;
    
    protected double gN_angle_p1_p2, gN_angle_p1_X, gN_angle_p2_X;
    
    public static boolean channel_test(Dihadrons variables) {
//        if (variables.helicity==0){ return false; }
        if (variables.Q2()<1) { return false; } 
        else if (variables.W()<2) { return false; } 
        else if (variables.y()>0.75) { return false; } 
//        else if (variables.Mx() < 0.80 || variables.Mx() > 1.00) { return false; }
//        else if (variables.Mx1()<1.50) { return false; } 
//        else if (variables.Mx2()<1.50) { return false; } 
        else if (variables.Mh()<0.60 || variables.Mh() > 0.90) { return false; }
//        else if (variables.xF1()>0.0) { return false; } 
//        else if (variables.xF2()<0.0) { return false; } 
//        else if (variables.z2()<0.2) { return false; } 
//        else if (variables.vz_e() < -12 || variables.vz_e() > 6) { return false; }
//        else if (variables.vz_p1() < -12 || variables.vz_p1() > 6) { return false; }
//        else if (variables.vz_p2() < -12 || variables.vz_p2() > 6) { return false; }
//        else if (variables.vz_e() < -10 || variables.vz_e() > 2.5) { return false; }
//        else if (variables.vz_p1() < -8 || variables.vz_p1() > 3) { return false; }
//        else if (variables.vz_p2() < -8 || variables.vz_p2() > 3) { return false; }
//
	return true;
    }
    
    public Dihadrons(DataEvent event, PhysicsEvent recEvent, int p1PID, int p1Index, int p2PID, int p2Index) {
        // provide the PDG PID of the two hadrons
        
        // load banks
        HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
        HipoDataBank configBank = (HipoDataBank) event.getBank("RUN::config");
        HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle"); // load particle bank
        
        helicity = eventBank.getByte("helicity", 0);
//        helicity = eventBank.getByte("helicityRaw", 0);
        runnum = configBank.getInt("run",0); // used for beam energy and polarization
        
    
        // default beam energy set to rga fall 2018
        double Eb = 10.6041; // rga fall 2018
        if (runnum >= 5032 && runnum <= 5666) { Eb = 10.6041; } // RGA Fall 2018
        else if (runnum >= 5875 && runnum <= 5995) { Eb = 6.535; } // RGK Fall 2018
        else if (runnum >= 5674 && runnum <= 5870) { Eb = 7.546; } // RGK Fall 2018
        else if (runnum >= 6616 && runnum <= 6783) { Eb = 10.1998; } // RGA Spring 2019
        else if (runnum >= 6120 && runnum <= 6399) { Eb = 10.5986; } // RGB Spring 2019
        else if (runnum >= 6409 && runnum <= 6604) { Eb = 10.1998; } // RGB Spring 2019
        else if (runnum >= 11093 && runnum <= 11283) { Eb = 10.4096; } // RGB Fall 2019
        else if (runnum >= 11284 && runnum <= 11300) { Eb = 4.17179; } // RGB Fall 2019
        else if (runnum >= 11323 && runnum <= 11571) { Eb = 10.3894; } // RGB Spring 2020
        else if (runnum >= 16290 && runnum <= 16353) { Eb = 10.5473; } // RGC preliminary summer 2022
//        Eb = 22.00;
        
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
    
        // calculate the Mx3, the missing mass of just Mx(e)
        LorentzVector lv_Mx3 = new LorentzVector(lv_q); lv_Mx3.add(lv_target); 
        Mx3 = lv_Mx3.mass();
        
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
        
        // set up hadrons and dihadron
        String p1Index_string = "["+p1PID+","+p1Index+"]";
        Particle p1 = recEvent.getParticle(p1Index_string);
        String p2Index_string = "["+p2PID+","+p2Index+"]";
        Particle p2 = recEvent.getParticle(p2Index_string);
        String combined_index_string = "["+p1PID+","+p1Index+"]+["+p2PID+","+p2Index+"]";
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
        
        
        // missing mass calculations
        // uncorrected
        LorentzVector lv_Mx_unc = new LorentzVector(lv_q); lv_Mx_unc.add(lv_target); lv_Mx_unc.sub(lv_p1); 
        lv_Mx_unc.sub(lv_p2);
        mass2_unc = lv_Mx_unc.mass2();
        
        // kinematics of hadrons
        p1_px = lv_p1.px(); p1_py = lv_p1.py(); p1_pz = lv_p1.pz(); p1_p = lv_p1.p(); p1_e = p1.e();
        p1_theta = p1.theta();
        p1_phi = p1.phi();
        if (p1_phi < 0) { p1_phi = 2*Math.PI + p1_phi; }
        
        p2_px = lv_p2.px(); p2_py = lv_p2.py(); p2_pz = lv_p2.pz(); p2_p = lv_p2.p(); p2_e = p2.e();
        p2_theta = p2.theta();
        p2_phi = p2.phi();
        if (p2_phi < 0) { p2_phi = 2*Math.PI + p2_phi; }
        p_px = lv_p.px(); p_py = lv_p.py(); p_pz = lv_p.pz(); p_p = lv_p.p(); p_e = lv_p.e();

        z = lv_p.e()/lv_q.e();
        z1 = lv_p1.e()/lv_q.e();
	z2 = lv_p2.e()/lv_q.e();
        
        // missing mass calculations
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target); lv_Mx.sub(lv_p1); lv_Mx.sub(lv_p2);
        Mx = lv_Mx.mass(); mass2 = lv_Mx.mass2(); 
        LorentzVector lv_Mx1 = new LorentzVector(lv_q); lv_Mx1.add(lv_target); lv_Mx1.sub(lv_p1); 
        Mx1 = lv_Mx1.mass(); // missing mass with p1 observed
        LorentzVector lv_Mx2 = new LorentzVector(lv_q); lv_Mx2.add(lv_target); lv_Mx2.sub(lv_p2);
        Mx2 = lv_Mx2.mass(); // missing mass with p2 observed
        
        // boost to gamma*-nucleon center of mass frame
        LorentzVector lv_p_gN = new LorentzVector(lv_p); lv_p_gN.boost(gNBoost);
        LorentzVector lv_p1_gN = new LorentzVector(lv_p1); lv_p1_gN.boost(gNBoost);
        LorentzVector lv_p2_gN = new LorentzVector(lv_p2); lv_p2_gN.boost(gNBoost);
        LorentzVector lv_e_gN = new LorentzVector(lv_e); lv_e_gN.boost(gNBoost);
        LorentzVector lv_target_gN = new LorentzVector(lv_target); lv_target_gN.boost(gNBoost);
        LorentzVector lv_Mx_gN = new LorentzVector(lv_Mx); lv_Mx_gN.boost(gNBoost);
        Vector3 lv_e_gN_unit = new Vector3();
        lv_e_gN_unit.setMagThetaPhi(1, lv_e_gN.theta(), lv_e_gN.phi());
        LorentzVector lv_q_gN = new LorentzVector(lv_q); lv_q_gN.boost(gNBoost);
        Vector3 lv_q_gN_unit = new Vector3();
        lv_q_gN_unit.setMagThetaPhi(1, lv_q_gN.theta(), lv_q_gN.phi());
        Vector3 lv_y_gN = new Vector3(); // I guess this isn't really a lorenzvector so the lv name is wrong 
        lv_y_gN = lv_q_gN.vect().cross(lv_e_gN.vect()); 
        // in gamma*-nucleon frame the z axis is along gamma* and the x axis is in the 
        // e-e' plane in the direction of e. the y axis is then the cross product of these two
        // calculate open angles between p1.p2, p1.X and p2.X
        gN_angle_p1_p2 = 180/Math.PI*Math.acos(lv_p1_gN.vect().dot(lv_p2_gN.vect())/(lv_p1_gN.vect().mag()*lv_p2_gN.vect().mag()));
        gN_angle_p1_X = 180/Math.PI*Math.acos(lv_p1_gN.vect().dot(lv_Mx_gN.vect())/(lv_p1_gN.vect().mag()*lv_Mx_gN.vect().mag()));
        gN_angle_p2_X = 180/Math.PI*Math.acos(lv_p2_gN.vect().dot(lv_Mx_gN.vect())/(lv_p2_gN.vect().mag()*lv_Mx_gN.vect().mag()));
        
        // boost to Breit infinite momentum frame
        LorentzVector lv_p_Breit = new LorentzVector(lv_p); lv_p_Breit.boost(BreitBoost);
        LorentzVector lv_p1_Breit = new LorentzVector(lv_p1); lv_p1_Breit.boost(BreitBoost);
        LorentzVector lv_p2_Breit = new LorentzVector(lv_p2); lv_p2_Breit.boost(BreitBoost);
        LorentzVector lv_e_Breit = new LorentzVector(lv_e); lv_e_Breit.boost(BreitBoost);
        Vector3 lv_e_Breit_unit = new Vector3();
        lv_e_Breit_unit.setMagThetaPhi(1, lv_e_Breit.theta(), lv_e_Breit.phi());
        LorentzVector lv_q_Breit = new LorentzVector(lv_q); lv_q_Breit.boost(BreitBoost);
        Vector3 lv_q_Breit_unit = new Vector3();
        lv_q_Breit_unit.setMagThetaPhi(1, lv_q_Breit.theta(), lv_q_Breit.phi());
        // note that in the Breit frame +z is antialigned to the direction of q
        
  
        // set up boost to dihadron rest frame (p1 and p2 center of mass frame)
        // just used to calculate theta
        Vector3 pCOMBoost = lv_p.boostVector(); pCOMBoost.negative();
        LorentzVector lv_p_COM = new LorentzVector(lv_p); lv_p_COM.boost(pCOMBoost);
        LorentzVector lv_p1_COM = new LorentzVector(lv_p1); lv_p1_COM.boost(pCOMBoost);
        LorentzVector lv_p2_COM = new LorentzVector(lv_p2); lv_p2_COM.boost(pCOMBoost);
//        theta = Math.acos(lv_p1_COM.vect().dot(lv_p_gN.vect()) / 
//                (lv_p1_COM.vect().mag()*lv_p_gN.vect().mag()));
        theta = Math.acos(lv_p1_COM.vect().dot(lv_p.vect()) / 
                (lv_p1_COM.vect().mag()*lv_p.vect().mag()));
        
        
        LorentzVector lv_p_boost_rho = new LorentzVector();
	lv_p_boost_rho.setPxPyPzM(lv_p.px(), lv_p.py(), lv_p.pz(), particle_mass(113));
        Vector3 pCOMBoostrho = lv_p_boost_rho.boostVector(); pCOMBoostrho.negative();
        LorentzVector lv_p_COM_rho = new LorentzVector(lv_p); lv_p_COM_rho.boost(pCOMBoostrho); 
        LorentzVector lv_p1_COM_rho = new LorentzVector(lv_p1); lv_p1_COM_rho.boost(pCOMBoostrho);
        LorentzVector lv_p2_COM_rho = new LorentzVector(lv_p2); lv_p2_COM_rho.boost(pCOMBoostrho);
        p1_COM_phi = lv_p1_COM_rho.phi(); if (p1_COM_phi < 0) { p1_COM_phi+=2*Math.PI; }
        p1_COM_theta = lv_p1_COM_rho.theta(); if (p1_COM_theta < 0) { p1_COM_theta+=2*Math.PI; }
        p2_COM_phi = lv_p2_COM_rho.phi(); if (p2_COM_phi < 0) { p2_COM_phi+=2*Math.PI; }
        p2_COM_theta = lv_p2_COM_rho.theta(); if (p2_COM_phi < 0) { p2_COM_phi+=2*Math.PI; }
        COM_Delta_phi = p1_COM_phi-p2_COM_phi; if (COM_Delta_phi < 0) { COM_Delta_phi+=2*Math.PI; }
        COM_Delta_theta = p1_COM_theta-p2_COM_theta; if (COM_Delta_theta < 0) { COM_Delta_theta+=2*Math.PI; }
//        if (COM_Delta_theta > Math.PI) {COM_Delta_theta=Math.PI-COM_Delta_theta;}
        COM_Delta_phi = 180/Math.PI*COM_Delta_phi;
        COM_Delta_theta = 180/Math.PI*COM_Delta_theta;
        COM_open_angle = 180/Math.PI*Math.acos(lv_p1_COM_rho.vect().dot(lv_p2_COM_rho.vect())/
                (lv_p1_COM_rho.vect().mag()*lv_p2_COM_rho.vect().mag()));
        
        
        
        Delta_phi = phi2 - phi1;
        if (Delta_phi < 0) { Delta_phi+=2*Math.PI; }
        
        pT = lv_q_gN_unit.cross(lv_p_gN.vect()).mag();
        pT1 = lv_q_gN_unit.cross(lv_p1_gN.vect()).mag();
        pT2 = lv_q_gN_unit.cross(lv_p2_gN.vect()).mag();
    
        xF =  2*(lv_p_gN.vect() .dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        xF1 = 2*(lv_p1_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        xF2 = 2*(lv_p2_gN.vect().dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        
        zeta = lv_p_gN.e()/lv_target_gN.e(); 
        zeta1 = lv_p1_gN.e()/lv_target_gN.e(); 
        zeta2 = lv_p2_gN.e()/lv_target_gN.e(); 
        
        p_gN_pz = lv_p_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p1_gN_pz = lv_p1_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p2_gN_pz = lv_p2_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p_Breit_pz = lv_p_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        p1_Breit_pz = lv_p1_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        p2_Breit_pz = lv_p2_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        
        // Breit frame rapidity
        eta = -0.5*Math.log((lv_p_Breit.e()+p_Breit_pz) /  (lv_p_Breit.e()-p_Breit_pz));
        eta1 = -0.5*Math.log((lv_p1_Breit.e()+p1_Breit_pz) /  (lv_p1_Breit.e()-p1_Breit_pz));
        eta2 = -0.5*Math.log((lv_p2_Breit.e()+p2_Breit_pz) /   (lv_p2_Breit.e()-p2_Breit_pz));
        
        // gamma*-nucleon frame rapidity
        eta_gN = 0.5*Math.log((lv_p_gN.e()+p_gN_pz) / (lv_p_gN.e()-p_gN_pz));
        eta1_gN = 0.5*Math.log((lv_p1_gN.e()+p1_gN_pz) /  (lv_p1_gN.e()-p1_gN_pz));
        eta2_gN = 0.5*Math.log((lv_p2_gN.e()+p2_gN_pz) /  (lv_p2_gN.e()-p2_gN_pz));
        
        Vector3 vecH = new Vector3();
        vecH.setMagThetaPhi(lv_p2_gN.vect().mag()/z2, lv_p2_gN.vect().theta(), lv_p2_gN.vect().phi());
        Vector3 vecR = new Vector3(lv_p1_gN.vect()); // not R yet
        vecR.setMagThetaPhi(lv_p1_gN.vect().mag()/z1, lv_p1_gN.vect().theta(), lv_p1_gN.vect().phi());
        vecR.sub(vecH); // this is really R now that the subtraction is done
        
	Vector3 vectRt = new Vector3();
	Vector3 R_Q = new Vector3();
							
	R_Q.setMagThetaPhi(vecR.dot(lv_q_gN_unit), lv_q_gN_unit.theta(), lv_q_gN_unit.phi());
	vectRt=vecR;
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
        
        Vector3 vT = new Vector3(lv_q_gN_unit.cross(lv_e_gN_unit)); vT.unit();
	Vector3 vTR = new Vector3(lv_q_gN_unit.cross(vectRt)); vTR.unit();
	Vector3 vTH = new Vector3(lv_q_gN_unit.cross(vectPhT)); vTH.unit();
        Vector3 vTH1 = new Vector3(lv_q_gN_unit.cross(vectPhT1)); vTH1.unit();
        Vector3 vTH2 = new Vector3(lv_q_gN_unit.cross(vectPhT2)); vTH2.unit();
        
        double cosPhiR = vT.dot(vTR);
	double sinPhiR = lv_e_gN.vect().cross(vectRt).dot(lv_q_gN_unit);
	double cosPhiH = vT.dot(vTH);
	double sinPhiH = lv_e_gN.vect().cross(vectPhT).dot(lv_q_gN_unit);
        double cosPhiH1 = vT.dot(vTH1);
	double sinPhiH1 = lv_e_gN.vect().cross(vectPhT1).dot(lv_q_gN_unit);
        double cosPhiH2 = vT.dot(vTH2);
	double sinPhiH2 = lv_e_gN.vect().cross(vectPhT2).dot(lv_q_gN_unit);

        // scaling
	double rScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectRt).mag();
	sinPhiR = sinPhiR/rScale;
	double hScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh).mag();
	sinPhiH = sinPhiH/hScale;
        double h1Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh1).mag();
	sinPhiH1 = sinPhiH1/h1Scale;
        double h2Scale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh2).mag();
	sinPhiH2 = sinPhiH2/h2Scale;
        
        phih = Math.acos(cosPhiH); 
        phi1 = Math.acos(cosPhiH1); 
        phi2 = Math.acos(cosPhiH2);
	phiR = Math.acos(cosPhiR); 
        
        if (sinPhiR < 0.0) {  phiR = 2*Math.PI - phiR; }
	if (sinPhiH < 0.0) { phih = 2*Math.PI - phih; }
        if (sinPhiH1 < 0.0) { phi1 = 2*Math.PI - phi1; }
        if (sinPhiH2 < 0.0) { phi2 = 2*Math.PI - phi2; }
        
        Delta_phi = phi2 - phi1;
        if (Delta_phi < 0) { Delta_phi+=2*Math.PI; }
         
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
        } else if ( runnum >= 11323 && runnum <= 11571) {
            return helicity;
        }
//        System.out.println("runnum not found, assigning helicity flip");
        return helicity;
    }
    
    public int get_runnum() { return runnum; }; // returns run number for polarizations and energy
    
    public int num_elec() { return num_elec; } // returns number of electrons
    
    public int num_piplus() { return num_piplus; } // returns number of piplus
    
    public int num_piminus() { return num_piminus; } // returns number of piminus
    
    public int num_kplus() { return num_kplus; }// returns number of kplus
    
    public int num_kminus() { return num_kminus; } // returns number of kminus
    
    public int num_protons() { return num_protons; } // returns number of protons
    
    public double test() { return test; } // returns test var

    public double Q2() { return Double.valueOf(Math.round(Q2*100000))/100000; } // returns Q2
    
    public double W() { return Double.valueOf(Math.round(W*100000))/100000; }// returns W
    
    public double gamma() { return Double.valueOf(Math.round(gamma*100000))/100000; } // returns gamma
    
    public double nu() { return Double.valueOf(Math.round(nu*100000))/100000; }// returns nu
    
    public double x() { return Double.valueOf(Math.round(x*100000))/100000; }// returns x
    
    public double y() { return Double.valueOf(Math.round(y*100000))/100000; }// returns y
    
    public double z() { return Double.valueOf(Math.round(z*100000))/100000; }// returns z
    
    public double z1() { return Double.valueOf(Math.round(z1*100000))/100000; }// returns z1
    
    public double z2() { return Double.valueOf(Math.round(z2*100000))/100000; }// returns z2

    public double Mx() { return Double.valueOf(Math.round(Mx*100000))/100000; }// returns Mx(ep1p2)
    
    public double Mx1() { return Double.valueOf(Math.round(Mx1*100000))/100000; }// returns Mx(ep1)
    
    public double Mx2() { return Double.valueOf(Math.round(Mx2*100000))/100000; }// returns Mx(ep2)
    
    public double Mx3() { return Double.valueOf(Math.round(Mx3*100000))/100000; }// returns Mx(e)
    
    public double Mh() { return Double.valueOf(Math.round(Mh*100000))/100000; }// returns Mh
    
    public double pT() { return Double.valueOf(Math.round(pT*100000))/100000; }// returns pT
    
    public double pT1() { return Double.valueOf(Math.round(pT1*100000))/100000; }// returns pT1
    
    public double pT2() { return Double.valueOf(Math.round(pT2*100000))/100000; }// returns pT2
    
    public double pTpT() { return Double.valueOf(Math.round(pT1*pT2*100000))/100000; } // returns ptpt
    
    public double xF() { return Double.valueOf(Math.round(xF*100000))/100000; }// returns xF
    
    public double xF1() { return Double.valueOf(Math.round(xF1*100000))/100000; }// returns xF1
    
    public double xF2() { return Double.valueOf(Math.round(xF2*100000))/100000; }// returns xF2
    
    public double zeta() { return Double.valueOf(Math.round(zeta*100000))/100000; }// returns zeta
    
    public double zeta1() { return Double.valueOf(Math.round(zeta1*100000))/100000; }// returns zeta1
    
    public double zeta2() { return Double.valueOf(Math.round(zeta2*100000))/100000; }// returns zeta2
    
    public double p1_Breit_pz() { return Double.valueOf(Math.round(p1_Breit_pz*100000))/100000; } 
    // returns p1 pz in Breit
    
    public double p2_Breit_pz() { return Double.valueOf(Math.round(p2_Breit_pz*100000))/100000; } 
    // returns p2 pz in Breit
    
    public double p1_gN_pz() { return Double.valueOf(Math.round(p1_gN_pz*100000))/100000; } // returns p1 pz in gN
    
    public double p2_gN_pz() { return Double.valueOf(Math.round(p2_gN_pz*100000))/100000; } // returns p2 pz in gN
    
    public double eta() { return Double.valueOf(Math.round(eta*100000))/100000; }// returns eta in the Breit frame
    
    public double eta1() { return Double.valueOf(Math.round(eta1*100000))/100000; }// returns eta1 in the Breit frame
    
    public double eta2() { return Double.valueOf(Math.round(eta2*100000))/100000; }// returns eta2 in the Breit frame
    
    public double eta_gN() { return Double.valueOf(Math.round(eta_gN*100000))/100000; }// returns eta_gN
    
    public double eta1_gN() { return Double.valueOf(Math.round(eta1_gN*100000))/100000; }// returns eta1_gN
    
    public double eta2_gN() { return Double.valueOf(Math.round(eta2_gN*100000))/100000; }// returns eta2_gN
    
    public double Delta_eta() { return Double.valueOf(Math.round((eta2-eta1)*100000))/100000; }// returns Delta_eta, 
    // should be Lorentz invariant
    // maybe set a check for difference of Delta_eta in Breit and gN frames in future...
    // make sure the sign is what you really want
    
    public double theta() { return Double.valueOf(Math.round(theta*100000))/100000; }
    // returns theta, the dihadron "decay angle"
    
    public double phih() { return Double.valueOf(Math.round(phih*100000))/100000; }// returns phih
    
    public double phiR() { return Double.valueOf(Math.round(phiR*100000))/100000; }// returns phiR
    
    public double phi1() { return Double.valueOf(Math.round(phi1*100000))/100000; }
    // returns phi1 (gamma*-nucleon frame)
    
    public double phi2() { return Double.valueOf(Math.round(phi2*100000))/100000; }
    // returns phi2 (gamma*-nucleon frame)
    
    public double Delta_phi() {return Double.valueOf(Math.round(Delta_phi*100000))/100000; }//returns Delta_phi (p1-p2)
     
    public double Depolarization_A() { return Double.valueOf(Math.round(Depolarization_A*100000))/100000; }
    // returns Depolarization_A
    
    public double Depolarization_B() { return Double.valueOf(Math.round(Depolarization_B*100000))/100000; }
    // returns Depolarization_B
    
    public double Depolarization_C() { return Double.valueOf(Math.round(Depolarization_C*100000))/100000; }
    // returns Depolarization_C
    
    public double Depolarization_V() { return Double.valueOf(Math.round(Depolarization_V*100000))/100000; }
    // returns Depolarization_V
    
    public double Depolarization_W() { return Double.valueOf(Math.round(Depolarization_W*100000))/100000; }
    // returns Depolarization_W
    
    public double e_px() { return Double.valueOf(Math.round(e_px*100000))/100000; }// returns electron lab frame px
    
    public double e_py() { return Double.valueOf(Math.round(e_py*100000))/100000; }// returns electron lab frame py
    
    public double e_pz() { return Double.valueOf(Math.round(e_pz*100000))/100000; }// returns electron lab frame pz
    
    public double e_p() { return Double.valueOf(Math.round(e_p*100000))/100000; }// returns electron lab frame p
    
    public double e_e() { return Double.valueOf(Math.round(e_e*100000))/100000; }// returns electron lab frame energy
    
    public double e_theta() { return Double.valueOf(Math.round(e_theta*100000))/100000; } // returns electron lab 
    // frame polar angle
    
    public double e_phi() { return Double.valueOf(Math.round(e_phi*100000))/100000; } // returns electron lab 
    // frame polar angle
    
    public double p_px() { return Double.valueOf(Math.round(p_px*100000))/100000; }// returns dihadron lab frame px
    
    public double p_py() { return Double.valueOf(Math.round(p_py*100000))/100000; }// returns dihadron lab frame py
    
    public double p_pz() { return Double.valueOf(Math.round(p_pz*100000))/100000; }// returns dihadron lab frame pz
    
    public double p_p() { return Double.valueOf(Math.round(p_p*100000))/100000; }// returns dihadron lab frame p
    
    public double p_e() { return Double.valueOf(Math.round(p_e*100000))/100000; }// returns dihadron lab frame energy
    
    public double p1_px() { return Double.valueOf(Math.round(p1_px*100000))/100000; }// returns hadron 1 lab frame px
    
    public double p1_py() { return Double.valueOf(Math.round(p1_py*100000))/100000; }// returns hadron 1 lab frame py
    
    public double p1_pz() { return Double.valueOf(Math.round(p1_pz*100000))/100000; }// returns hadron 1 lab frame pz
    
    public double p1_p() { return Double.valueOf(Math.round(p1_p*100000))/100000; }// returns hadron 1 lab frame p
    
    public double p1_px_unc() { return Double.valueOf(Math.round(p1_px_unc*100000))/100000; }// returns hadron 1 lab frame px
    
    public double p1_py_unc() { return Double.valueOf(Math.round(p1_py_unc*100000))/100000; }// returns hadron 1 lab frame py
    
    public double p1_pz_unc() { return Double.valueOf(Math.round(p1_pz_unc*100000))/100000; }// returns hadron 1 lab frame pz
    
    public double p1_p_unc() { return Double.valueOf(Math.round(p1_p_unc*100000))/100000; }// returns hadron 1 lab frame p
    
    public double mass2() { return Double.valueOf(Math.round(mass2*100000))/100000; }// returns mass2
    
    public double mass2_unc() { return Double.valueOf(Math.round(mass2_unc*100000))/100000; }// returns mass2
    
    public double p1_e() { return Double.valueOf(Math.round(p1_e*100000))/100000; }// returns hadron 1 lab frame energy
    
    public double p1_theta() { return Double.valueOf(Math.round(p1_theta*100000))/100000; } // returns p1 lab 
    // frame polar angle
    
    public double p1_phi() { return Double.valueOf(Math.round(p1_phi*100000))/100000; } // returns p1 lab 
    // frame polar angle
    
    public double p2_px() { return Double.valueOf(Math.round(p2_px*100000))/100000; }// returns hadron 2 lab frame px
    
    public double p2_py() { return Double.valueOf(Math.round(p2_py*100000))/100000; }// returns hadron 2 lab frame py
    
    public double p2_pz() { return Double.valueOf(Math.round(p2_pz*100000))/100000; }// returns hadron 2 lab frame pz
    
    public double p2_p() { return Double.valueOf(Math.round(p2_p*100000))/100000; }// returns hadron 2 lab frame p
    
    public double p2_e() { return Double.valueOf(Math.round(p2_e*100000))/100000; }// returns hadron 2 lab frame energy
    
    public double p2_theta() { return Double.valueOf(Math.round(p2_theta*100000))/100000; } // returns p2 lab frame polar angle
    
    public double p2_phi() { return Double.valueOf(Math.round(p2_phi*100000))/100000; } // returns p2 lab frame polar angle
    
    public double vx_e() { return Double.valueOf(Math.round(vx_e*100000))/100000; }// returns electron x vertex
    
    public double vx_p1() { return Double.valueOf(Math.round(vx_p1*100000))/100000; }// returns p1 x vertex
    
    public double vx_p2() { return Double.valueOf(Math.round(vx_p2*100000))/100000; }// returns p2 x vertex
    
    public double vy_e() { return Double.valueOf(Math.round(vy_e*100000))/100000; }// returns electron y vertex
    
    public double vy_p1() { return Double.valueOf(Math.round(vy_p1*100000))/100000; }// returns p1 y vertex
    
    public double vy_p2() { return Double.valueOf(Math.round(vy_p2*100000))/100000; }// returns p2 y vertex
    
    public double vz_e() { return Double.valueOf(Math.round(vz_e*100000))/100000; }// returns electron z vertex
    
    public double vz_p1() { return Double.valueOf(Math.round(vz_p1*100000))/100000; }// returns p1 z vertex
    
    public double vz_p2() { return Double.valueOf(Math.round(vz_p2*100000))/100000; }// returns p2 z vertex
    
    public double e_nphe() { return Double.valueOf(Math.round(e_nphe*100000))/100000; }// electron nphe
    
    public double p1_nphe() { return Double.valueOf(Math.round(p1_nphe*100000))/100000; }// p1 nphe
    
    public double p2_nphe() { return Double.valueOf(Math.round(p2_nphe*100000))/100000; }// p2 nphe
    
    public double e_chi2pid() { return Double.valueOf(Math.round(e_chi2pid*100000))/100000; }// electron chi2pid
    
    public double p1_chi2pid() { return Double.valueOf(Math.round(p1_chi2pid*100000))/100000; }// p1 chi2pid
    
    public double p2_chi2pid() { return Double.valueOf(Math.round(p2_chi2pid*100000))/100000; }// p2 chi2pid
   
    public double p1_COM_phi() { return Double.valueOf(Math.round(p1_COM_phi*100000))/100000; }// p1 COM phi
    
    public double p1_COM_theta() { return Double.valueOf(Math.round(p1_COM_theta*100000))/100000; }// p1_COM_theta
    
    public double p2_COM_phi() { return Double.valueOf(Math.round(p2_COM_phi*100000))/100000; }// p2_COM_phi
    
    public double p2_COM_theta() { return Double.valueOf(Math.round(p2_COM_theta*100000))/100000; }// p2_COM_theta
    
    public double COM_Delta_phi() { return Double.valueOf(Math.round(COM_Delta_phi*100000))/100000; }// COM_Delta_phi
    
    public double COM_Delta_theta() { return Double.valueOf(Math.round(COM_Delta_theta*100000))/100000; }// COM_Delta_theta
    
    public double COM_open_angle() { return Double.valueOf(Math.round(COM_open_angle*100000))/100000; }// COM_open_angle
    
    public double gN_angle_p1_p2() { return Double.valueOf(Math.round(gN_angle_p1_p2*100000))/100000; }// gN_angle_p1_p2
    
    public double gN_angle_p1_X() { return Double.valueOf(Math.round(gN_angle_p1_X*100000))/100000; }// gN_angle_p1_X
    
    public double gN_angle_p2_X() { return Double.valueOf(Math.round(gN_angle_p2_X*100000))/100000; }// gN_angle_p1_X
    
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
