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

//import org.jlab.detector.scalers.DaqScalersSequence;

public class analysis_fitter extends GenericKinematicFitter {

    protected final Double mybeam;

    public analysis_fitter(double beam) {
        super(beam);
        mybeam = beam;
    } 
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    public boolean electron_test(int particle_Index, double p,
            HipoDataBank rec_Bank, HipoDataBank cal_Bank,  HipoDataBank track_Bank, 
            HipoDataBank traj_Bank, HipoDataBank run_Bank, HipoDataBank cc_Bank) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();
        
        return true
//            && p > 2.2 // higher cut ultimately enforced when we cut on y, this speeds processing
            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && pid_cuts.calorimeter_energy_cut(particle_Index, cal_Bank) 
            && pid_cuts.calorimeter_sampling_fraction_cut(particle_Index, p, run_Bank, cal_Bank)
            && pid_cuts.calorimeter_diagonal_cut(particle_Index, p, cal_Bank)
            && generic_tests.vertex_cut(particle_Index, rec_Bank, run_Bank)    
            && fiducial_cuts.pcal_fiducial_cut(particle_Index, rec_Bank, cal_Bank)
            && fiducial_cuts.pass1_dc_fiducial_cut(particle_Index,rec_Bank,track_Bank,traj_Bank,run_Bank)
            ;
    }
    
    public boolean pion_test(int particle_Index, int pid, float vz, double trigger_electron_vz, HipoDataBank rec_Bank, 
            HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();
        
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
//            && p > 1.25
//            && p < 5.00 // this wasn't used in the dihadron publication but was used in the submitted single pion
            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && generic_tests.vertex_cut(particle_Index, rec_Bank, run_Bank) 
            && pid_cuts.charged_hadron_pass2_chi2pid_cut(particle_Index, rec_Bank)
//            && pid_cuts.charged_pion_generic_chi2pid_cut(particle_Index, rec_Bank)
//            && pid_cuts.charged_hadron_chi2pid_cut(particle_Index, rec_Bank)
            && fiducial_cuts.pass1_dc_fiducial_cut(particle_Index,rec_Bank,track_Bank,traj_Bank,run_Bank)
              ;
    }
    
    public boolean kaon_test(int particle_Index, int pid, float vz, double trigger_electron_vz, HipoDataBank rec_Bank, 
            HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();
        
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
//            && p > 1.00
//            && p < 3.5 
            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && generic_tests.vertex_cut(particle_Index, rec_Bank, run_Bank) 
            && pid_cuts.charged_hadron_pass2_chi2pid_cut(particle_Index, rec_Bank)
//            && pid_cuts.charged_hadron_chi2pid_cut(particle_Index, rec_Bank)
            && fiducial_cuts.pass1_dc_fiducial_cut(particle_Index,rec_Bank,track_Bank,traj_Bank,run_Bank)
              ;
    }
    
    public boolean proton_test(int particle_Index, int pid, float vz, double trigger_electron_vz, 
            HipoDataBank rec_Bank, HipoDataBank cal_Bank, HipoDataBank track_Bank, 
            HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();
        
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
//            && p > 0.4
            && generic_tests.vertex_cut(particle_Index, rec_Bank, run_Bank) 
            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && fiducial_cuts.pass1_dc_fiducial_cut(particle_Index,rec_Bank,track_Bank,traj_Bank,run_Bank)
            && pid_cuts.charged_hadron_pass2_chi2pid_cut(particle_Index, rec_Bank)
//            && charged_hadron_chi2pid_cut(particle_Index, rec_Bank)
              ;
    }
    
    public boolean photon_test(int particle_Index, HipoDataBank rec_Bank, HipoDataBank cal_Bank, 
            LorentzVector lv_e) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();
        
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        LorentzVector lv_gamma = new LorentzVector();
        lv_gamma.setPxPyPzM(px, py, pz, 0.0);
        
        return true
            && p > 0.50
//            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && fiducial_cuts.pcal_fiducial_cut(particle_Index, rec_Bank, cal_Bank)
            && pid_cuts.beta_cut(particle_Index, rec_Bank)
//            && pid_cuts.e_gamma_open_angle_cut(lv_e, lv_gamma)
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
            HipoDataBank track_Bank = (HipoDataBank) event.getBank("REC::Track");
            HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");
            HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config");
            
            double vz_e = -999;
        
            LorentzVector lv_e = new LorentzVector(); 
            if (rec_Bank.getInt("pid", 0) == 11) { 
                // trigger particle was an electron
                // highest momentum electron listed first (used for DIS calculations)
                float px = rec_Bank.getFloat("px", 0);
                float py = rec_Bank.getFloat("py", 0);
                float pz = rec_Bank.getFloat("pz", 0);
                double p = Math.sqrt(px*px+py*py+pz*pz);
                lv_e.setPxPyPzM(px, py, pz, 0.0005109989461);
                vz_e = rec_Bank.getFloat("vz",0);
                
            } else { return physEvent; } // trigger particle was not an electron
            
            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
                int pid = rec_Bank.getInt("pid", particle_Index);
                float px = rec_Bank.getFloat("px", particle_Index);
                float py = rec_Bank.getFloat("py", particle_Index);
                float pz = rec_Bank.getFloat("pz", particle_Index);
                float vx = rec_Bank.getFloat("vx",particle_Index);
                float vy = rec_Bank.getFloat("vy",particle_Index);
                float vz = rec_Bank.getFloat("vz",particle_Index);
                double p = Math.sqrt(px*px+py*py+pz*pz);
                
                energy_loss energy_loss = new energy_loss();
                
                if (pid == 11 && electron_test(particle_Index, p, rec_Bank, cal_Bank, track_Bank, 
                        traj_Bank, run_Bank, cc_Bank)) {
                    // this checks all of the PID requirements, if it passes all of them the electron is 
                    // added to the event below
//                  double fe = energy_loss.pass2_fd_energy_loss(run_Bank.getFloat("torus", 0), pid, px, py, pz);
                    double fe = 1;
                    Particle part = new Particle(pid,fe*px,fe*py,fe*pz,vx,vy,vz_e);
                    physEvent.addParticle(part);
                }
                
                if (Math.abs(pid)==211 && pion_test(particle_Index, pid, vz, vz_e, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                    // check for pion PID
                    
//                   double fe = EnergyLoss(run_Bank.getFloat("torus", 0), pid, px, py, pz);
                   double fe = 1;
                   Particle part = new Particle(pid,fe*px,fe*py,fe*pz,vx,vy,vz);
                   physEvent.addParticle(part);   
                }
                
                if (Math.abs(pid)==321 && kaon_test(particle_Index, pid, vz, vz_e, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                    // check for pion PID
                   
//                   double fe = energy_loss.pass2_fd_energy_loss(run_Bank.getFloat("torus", 0), pid, px, py, pz);
                   double fe = 1;
                   Particle part = new Particle(pid,fe*px,fe*py,fe*pz,vx,vy,vz);
                   physEvent.addParticle(part);   
                }
                
                if (pid==2212 && proton_test(particle_Index, pid, vz, vz_e, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                   
//                   double fe = energy_loss.pass2_energy_loss(run_Bank.getFloat("torus", 0), pid, px, py, pz);
                   double fe = 1;
                   Particle part = new Particle(pid,fe*px,fe*py,fe*pz,vx,vy,vz);
                   physEvent.addParticle(part);   
               }
                
                if (pid==22 && photon_test(particle_Index, rec_Bank, cal_Bank, lv_e)) {
                   
                   Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
                   physEvent.addParticle(part);   
               }
            }

            
//            // check for >= 2 photons (in order to reconstruct pi0)
//            int num_EB_photons = 0;
//            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
//                int pid = rec_Bank.getInt("pid", particle_Index);
//                if (pid == 22) num_EB_photons ++;
//            }
//            
//            int num_analysis_photons = 0; // number of photons meeting cuts
//            LorentzVector lv_gamma_1 = new LorentzVector();
//            LorentzVector lv_gamma_2 = new LorentzVector();
//            LorentzVector lv_gamma_3 = new LorentzVector();
//            LorentzVector lv_gamma_4 = new LorentzVector();
//            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
//                int pid = rec_Bank.getInt("pid", particle_Index);
//                
//                if (pid!=22 || trigger_electron_vz == -9 || num_EB_photons <2) { continue; }
//                // requires the particle to be photon by EventBuilder & for an electron to have been assigned to event
//                
//                
//                // load momenta and vertices
//                float px = rec_Bank.getFloat("px", particle_Index);
//                float py = rec_Bank.getFloat("py", particle_Index);
//                float pz = rec_Bank.getFloat("pz", particle_Index);
//                float vx = rec_Bank.getFloat("vx",particle_Index);
//                float vy = rec_Bank.getFloat("vy",particle_Index);
//                float vz = rec_Bank.getFloat("vz",particle_Index);
//                LorentzVector lv_gamma = new LorentzVector(); lv_gamma.setPxPyPzM(px, py, pz, 0);
//                if (photon_test(particle_Index, rec_Bank, cal_Bank, lv_e, lv_gamma)) {
//                   
//                   Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
//                   physEvent.addParticle(part);  
//                   num_analysis_photons++;
//                   if (num_analysis_photons==1) { lv_gamma_1.setPxPyPzM(px,py,pz,0); }
//                   else if (num_analysis_photons==2) { lv_gamma_2.setPxPyPzM(px,py,pz,0); }
//                   else if (num_analysis_photons==3) { lv_gamma_3.setPxPyPzM(px,py,pz,0); }
//                   else if (num_analysis_photons==4) { lv_gamma_4.setPxPyPzM(px,py,pz,0); }
//               }
//            }
//            
//            double lower_mgg = 0.11;
//            double upper_mgg = 0.16; // signal
////            double lower_mgg = 0.22;
////            double upper_mgg = 0.45; // background
//            if (num_analysis_photons>1) {
//                LorentzVector lv_gamma_1_2 = new LorentzVector(lv_gamma_1); lv_gamma_1_2.add(lv_gamma_2);
//                if (lv_gamma_1_2.mass() > lower_mgg && lv_gamma_1_2.mass() < upper_mgg) {
//                    Particle part = new Particle(111,lv_gamma_1_2.px(),lv_gamma_1_2.py(),lv_gamma_1_2.pz(),0,0,0);
//                    physEvent.addParticle(part);
//                }
//                if (num_analysis_photons>2) {
//                    LorentzVector lv_gamma_1_3 = new LorentzVector(lv_gamma_1); lv_gamma_1_3.add(lv_gamma_3);
//                    LorentzVector lv_gamma_2_3 = new LorentzVector(lv_gamma_2); lv_gamma_2_3.add(lv_gamma_3);
//                    if (lv_gamma_1_3.mass() > lower_mgg && lv_gamma_1_3.mass() < upper_mgg) {
//                        Particle part = new Particle(111,lv_gamma_1_3.px(),lv_gamma_1_3.py(),lv_gamma_1_3.pz(),0,0,0);
//                        physEvent.addParticle(part);
//                    }
//                    if (lv_gamma_2_3.mass() > lower_mgg && lv_gamma_2_3.mass() < upper_mgg) {
//                        Particle part = new Particle(111,lv_gamma_2_3.px(),lv_gamma_2_3.py(),lv_gamma_2_3.pz(),0,0,0);
//                        physEvent.addParticle(part);
//                    }
//                    if (num_analysis_photons>3) {
//                        LorentzVector lv_gamma_1_4 = new LorentzVector(lv_gamma_1); lv_gamma_1_4.add(lv_gamma_4);
//                        LorentzVector lv_gamma_2_4 = new LorentzVector(lv_gamma_2); lv_gamma_2_3.add(lv_gamma_4);
//                        LorentzVector lv_gamma_3_4 = new LorentzVector(lv_gamma_3); lv_gamma_2_3.add(lv_gamma_4);
//                        if (lv_gamma_1_4.mass() > lower_mgg && lv_gamma_1_4.mass() < upper_mgg) {
//                            Particle part = new Particle(111,lv_gamma_1_4.px(),lv_gamma_1_4.py(),lv_gamma_1_4.pz(),
//                                    0,0,0);
//                            physEvent.addParticle(part);
//                        }
//                        if (lv_gamma_2_4.mass() > lower_mgg && lv_gamma_2_4.mass() < upper_mgg) {
//                            Particle part = new Particle(111,lv_gamma_2_4.px(),lv_gamma_2_4.py(),lv_gamma_2_4.pz(),
//                                    0,0,0);
//                            physEvent.addParticle(part);
//                        }
//                        if (lv_gamma_3_4.mass() > lower_mgg && lv_gamma_3_4.mass() < upper_mgg) {
//                            Particle part = new Particle(111,lv_gamma_3_4.px(),lv_gamma_3_4.py(),lv_gamma_3_4.pz(),
//                                    0,0,0);
//                            physEvent.addParticle(part);
//                        }
//                    }
//                }
//            }
            
            
            
            return physEvent;
        }
    return new PhysicsEvent(this.mybeam);
    }
}