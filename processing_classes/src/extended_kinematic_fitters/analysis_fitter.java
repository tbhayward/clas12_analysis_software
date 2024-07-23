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
    
    public double EnergyLoss(float polarity, int pid, double px, double py, double pz) {
        double dp = 0; // scale size
        double r = Math.pow(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2),0.5);
        double theta = (180/Math.PI)*Math.acos(pz/r);
        double p = Math.sqrt(px*px+py*py+pz*pz);
        
        if (polarity < 0) {
            if (theta < 27) {
                switch (pid) { 
                    case 11:
                        dp = 0.00730953+Math.exp(-14.0303+0.957154*p);
                        break;
                        
                    case 211:
                        dp = -36.3456+Math.exp(3.59313-0.0000065*p);
                        break;

                    case -211:
                        dp = 0.00396596+Math.exp(6.14785-16.8218*p);
                        break;

                    case 321:
                        dp = -34.5647+Math.exp(3.54297-0.0000136581*p);
                        break;

                    case -321:
                        dp = 0.00340915+Math.exp(-0.567643-6.31624*p);
                        break;
                        
                    case 2212:
                        dp = 0.00428181+Math.exp(-2.62675-3.39989*p);
                        break;
                }
            } else if (theta >= 27) {
                switch (pid) { 
                    case 11:
                        dp = 0.00552679+Math.exp(-4.2381-0.660486*p);
                        break;
                        
                    case 211:
                        dp = -45.9594+Math.exp(3.82795-0.0000358352*p);
                        break;

                    case -211:
                        dp = -0.000914991+Math.exp(-4.29887-0.555948*p);
                        break;

                    case 321:
                        dp = 0.00514143+Math.exp(-3.65383-2.19434*p);
                        break;

                    case -321:
                        dp = 0.00106601+Math.exp(-3.87516-1.04774*p);
                        break;
                        
                    case 2212:
                        dp = 0.00597654+Math.exp(-2.4691-2.37887*p);
                        break;
                }
            }
        } else if (polarity > 0) {
            if (theta < 27) {
                switch (pid) { 
                    case 11:
                        dp = -30.2+Math.exp(3.40808-0.0000048*p);
                        break;
                    
                    case 211:
                        dp = 0.00392301+Math.exp(3.83445-13.9503*p);
                        break;

                    case -211:
                        dp = -39.0688+Math.exp(3.66543-0.000006173*p);
                        break;

                    case 321:
                        dp = 0.00377599+Math.exp(-0.386989-7.33292*p);
                        break;

                    case -321:
                        dp = 0.00377308+Math.exp(-2.11549-6.31415*p);
                        break;
                        
                    case 2212:
                        dp = 0.00443825+Math.exp(-0.883026-5.07364*p);
                        break;
                }
            } else if (theta >= 27) {
                switch (pid) { 
                    case 11:
                        dp = -0.0154673+Math.exp(-3.42119-0.096196*p);
                        break;
                        
                    case 211:
                        dp = 0.00185544+Math.exp(-4.41637-0.878962*p);
                        break;

                    case -211:
                        dp = -91.7796+Math.exp(4.51949-0.0000206923*p);
                        break;

                    case 321:
                        dp = 0.00279201+Math.exp(-3.64049-1.65387*p);
                        break;

                    case -321:
                        dp = 0.00430273+Math.exp(-3.86056-1.6304*p);
                        break;
                        
                    case 2212:
                        dp = 0.00475655+Math.exp(-1.99834-2.81217*p);
                        break;
                }
            }
        }
        
        double fe = (dp+p)/p;
        return fe;
        
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////  
//    public boolean current_test(DataEvent event) {
//        HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config");
//        
//        DaqScalersSequence chargeSeq = DaqScalersSequence.readSequence(inputList);
//        // get beam current
//        double current = chargeSeq.getInterval(timeStamp).getBeamCurrent();
//    }


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // PID enhancements

    public boolean calorimeter_energy_cut(int particle_Index, HipoDataBank cal_Bank) {
        // Iterate through the rows of the data bank
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            // Get the pindex and layer values for the current row
            int pindex = cal_Bank.getInt("pindex", current_Row);
            int layer = cal_Bank.getInt("layer", current_Row);

            // Check if the pindex and layer values match the specified particle and layer
            if (pindex == particle_Index && layer == 1) {
                // Get the energy value for the current row
                float energy = cal_Bank.getFloat("energy", current_Row);

                // Check if the energy is greater than the threshold value
                return energy > 0.07;
            }
        }

        // If no matching rows are found, return false
        return false;
    }

    public boolean calorimeter_sampling_fraction_cut(int particle_Index, double p, HipoDataBank run_Bank, 
            HipoDataBank cal_Bank) {
        double scale = 3.5; // how many std away from mean to cut on
        int sector = -1;
        double cal_energy = 0;
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==particle_Index)  {
                sector = cal_Bank.getInt("sector", current_Row) - 1; // subtract one to start at index of 0
                cal_energy+= cal_Bank.getFloat("energy", current_Row);
            }
        }
        
        int runnum = run_Bank.getInt("run", 0);
        
        double e_cal_sampl_mu_RGA[][] = {{0.2531, 0.2550, 0.2514, 0.2494, 0.2528, 0.2521},
            {-0.6502, -0.7472, -0.7674, -0.4913, -0.3988, -0.703},
            {4.939, 5.350, 5.102, 6.440, 6.149, 4.957}};
        
        double e_cal_sampl_sigma_RGA[][] = {{0.002726, 0.004157, 0.00522, 0.005398, 0.008453, 0.006553}, 
            {1.062, 0.859, 0.5564, 0.6576, 0.3242, 0.4423}, 
            {-4.089, -3.318, -2.078, -2.565, -0.8223, -1.274}};
        
        double e_cal_sampl_mu_RGBSp19[][] = {{0.2520, 0.2520, 0.2479, 0.2444, 0.2463, 0.2478},
            {-0.8615, -0.8524, -0.6848, -0.5521, -0.5775, -0.7327},
            {5.596, 6.522, 5.752, 5.278, 6.430, 5.795}};
        
        double e_cal_sampl_sigma_RGBSp19[][] = {{-0.02963, -0.1058, -0.05087, -0.04524, -0.02951, -0.01769}, 
            {20.4, 129.3, 0.6191, 0.6817, 20.84, 8.44}, 
            {-41.44, -101.6, -2.673, -2.606, -42.67, -21.73}};
        
        double e_cal_sampl_mu_RGBW20[][] = {{0.2433, 0.2421, 0.2415, 0.2486, 0.2419, 0.2447},
            {-0.8052, -1.0495, -1.1747, -0.5170, -0.6840, -0.9022},
            {5.2750, 4.4886, 4.4935, 5.9044, 5.6716, 4.9288}};
        
        double e_cal_sampl_sigmaRGBW20[][] = {{0.0120, 0.0164, 0.0120, 0.0108, 0.0147, 0.0077}, 
            {0.1794, 0.1519, 0.1379, 0.1838, 0.0494, 0.3509}, 
            {-0.0695, 0.1553, 0.3300, 0.4330, 1.1032, -0.7996}};
            
            
        // sampling fraction is cal_energy/p
        if( (runnum>=11323 && runnum<=11571) || (runnum>=11093 && runnum<=11300) ) { 
        // RGB winter 2020 // (also using this for RGB fall 2019, but it should be updated! TODO)
        
            double mean = e_cal_sampl_mu_RGBW20[0][sector]+(e_cal_sampl_mu_RGBW20[1][sector]/1000)*
                (p-e_cal_sampl_mu_RGBW20[2][sector])*(p-e_cal_sampl_mu_RGBW20[2][sector]);
        
            double std = e_cal_sampl_sigmaRGBW20[0][sector] + e_cal_sampl_sigmaRGBW20[1][sector] / 
                (10 * (p-e_cal_sampl_sigmaRGBW20[2][sector]));
        
            return ((cal_energy/p) > (mean-scale*std)) && ((cal_energy/p) < (mean+scale*std));
        } else if ( runnum>=6120 && runnum<=6604 ) { // RGB Sp19
        
            double mean = e_cal_sampl_mu_RGBSp19[0][sector]+(e_cal_sampl_mu_RGBSp19[1][sector]/1000)*
                (p-e_cal_sampl_mu_RGBSp19[2][sector])*(p-e_cal_sampl_mu_RGBSp19[2][sector]);
        
            double std = e_cal_sampl_sigma_RGBSp19[0][sector] + e_cal_sampl_sigma_RGBSp19[1][sector] / 
                (10 * (p-e_cal_sampl_sigma_RGBSp19[2][sector]));
       
            
            return ((cal_energy/p) > (mean-scale*std)) && ((cal_energy/p) < (mean+scale*std));
        } else {
            double mean = e_cal_sampl_mu_RGA[0][sector]+(e_cal_sampl_mu_RGA[1][sector]/1000)*
                (p-e_cal_sampl_mu_RGA[2][sector])*(p-e_cal_sampl_mu_RGA[2][sector]);
        
            double std = e_cal_sampl_sigma_RGA[0][sector] + e_cal_sampl_sigma_RGA[1][sector] / 
                (10 * (p-e_cal_sampl_sigma_RGA[2][sector]));
        
            return ((cal_energy/p) > (mean-scale*std)) && ((cal_energy/p) < (mean+scale*std));
        }
    }
    
    public boolean calorimeter_diagonal_cut(int particle_Index, double p, HipoDataBank cal_Bank) {
        // Only apply diagonal cut above 4.5 GeV
        if (p < 4.5) { return true; }

        // Initialize the sum of energy for PCAL and ECAL inner layers
        double pcal_plus_ecal_inner = 0;

        // Iterate through the rows of the data bank
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            // Check if the current row corresponds to the particle index and the required layer (1 or 4)
            int pindex = cal_Bank.getInt("pindex", current_Row);
            int layer = cal_Bank.getInt("layer", current_Row);

            if (pindex == particle_Index && (layer == 1 || layer == 4)) {
                // Add the energy value to the sum
                pcal_plus_ecal_inner += cal_Bank.getFloat("energy", current_Row);
            }
        }

        // Check if the energy ratio is above the threshold
        return 0.2 < pcal_plus_ecal_inner / p;
    }
    
    public boolean hadron_pass2_cut(int particle_Index, HipoDataBank rec_Bank) {
        // Retrieve values from the data bank
        float chi2pid = rec_Bank.getFloat("chi2pid", particle_Index);
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);

        // Calculate the momentum
        double p = Math.sqrt(px * px + py * py + pz * pz);
        // Determine the constant C based on the particle ID
        int pid = rec_Bank.getInt("pid", particle_Index);
        double mu = 1;
        double sigma =1 ;
        
        if (pid == -211) {
            mu = -0.063;
            sigma = 0.947;
            if (p < 3) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else {
                return chi2pid > mu-3*sigma && chi2pid < chi2pid_cut(-0.7, mu+3*sigma+0.7, 0.9, 3, p);
            }
        }
        
        if (pid == -321) {
            mu = 0.115;
            sigma = 0.958;
            if (p < 2) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else {
                return chi2pid > chi2pid_cut(1.7, mu-3*sigma-1.7, 1.0, 2.0, p) && chi2pid < mu+3*sigma;
            }
        }
        
        if (pid == 211) {
            mu = -0.067;
            sigma = 0.956;
            if (p < 3.5) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else {
                return chi2pid > mu-3*sigma && chi2pid < chi2pid_cut(-0.55, mu+3*sigma+0.55, 0.55, 3.5, p);
            }
        }
        
        if (pid == 321) {
            mu = 0.082;
            sigma = 0.985;
            if (p < 2) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else if (p > 2 && p < 2.5) {
                return chi2pid > chi2pid_cut(1.2, mu-3*sigma-1.2, 0.6, 2, p) && chi2pid < mu+3*sigma;
            } else {
                return chi2pid > chi2pid_cut(1.2, mu-3*sigma-1.2, 0.6, 2, p) && chi2pid < chi2pid_cut(2.6, mu+3*sigma-2.6, 0.3, 2.5, p);
            }
        }
        
        if (pid == 2212) {
            mu = 0.372;
            sigma = 1.192;
            if (p < 2) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else {
                return chi2pid > chi2pid_cut(2, mu-3*sigma-2, 0.9, 2, p) && chi2pid < mu+3*sigma;
            }
        }
        
        return false; // not a hadron? 
    }
    
    public double chi2pid_cut(double a, double b, double c, double x0, double p) {
        return a + b * Math.exp(-(p-x0)/c);
    }
    
    public boolean pion_chi2pid_cut(int particle_Index, HipoDataBank rec_Bank) {
        // Retrieve values from the data bank
        float chi2pid = rec_Bank.getFloat("chi2pid", particle_Index);
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);

        // Calculate the momentum
        double p = Math.sqrt(px * px + py * py + pz * pz);

        // Determine the constant C based on the particle ID
        int pid = rec_Bank.getInt("pid", particle_Index);
        double C = (pid == 211) ? 0.88 : 0.93; 

        // Check if the momentum is less than the threshold
        if (p < 2.44) {
            // Simplify the range check using Math.abs
            return Math.abs(chi2pid) < 3 * C;
        } else {
            // Calculate the upper limit of the range check for chi2pid
            double upperLimit = C * (0.00869 + 14.98587 * Math.exp(-p / 1.18236) + 1.81751 * Math.exp(-p / 4.86394));
            // Check if chi2pid is within the range
            return chi2pid < upperLimit && chi2pid > -3 * C;
        }
    }

    public boolean hadron_chi2pid_cut(int particle_Index, HipoDataBank rec_Bank) {
        float chi2pid = rec_Bank.getFloat("chi2pid", particle_Index);
        
        return Math.abs(chi2pid) < 5;
    }
    
    public boolean random_photon_cut(int particle_Index, HipoDataBank rec_Bank) {
        // remove random photons
        float beta = rec_Bank.getFloat("beta", particle_Index);
        return (beta > 0.90 && beta < 1.10);
    }
    
    public boolean open_angle_cut(LorentzVector lv_e, LorentzVector lv_gamma) {
        return 180/Math.PI*Math.acos(lv_e.vect().dot(lv_gamma.vect())/
                (lv_e.vect().mag()*lv_gamma.vect().mag())) > 8; 
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    public boolean electron_test(int particle_Index, double p, float vz, double trigger_electron_vz, 
            HipoDataBank rec_Bank, HipoDataBank cal_Bank,  HipoDataBank track_Bank, 
            HipoDataBank traj_Bank, HipoDataBank run_Bank, HipoDataBank cc_Bank) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        
        return true
//            && p > 2.2 // higher cut ultimately enforced when we cut on y, this speeds processing
            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && calorimeter_energy_cut(particle_Index, cal_Bank) 
            && calorimeter_sampling_fraction_cut(particle_Index, p, run_Bank, cal_Bank)
            && calorimeter_diagonal_cut(particle_Index, p, cal_Bank)
            && generic_tests.vertex_cut(particle_Index, trigger_electron_vz, rec_Bank, run_Bank)    
            && fiducial_cuts.pcal_fiducial_cut(particle_Index, rec_Bank, cal_Bank)
            && fiducial_cuts.pass1_dc_fiducial_cut(particle_Index,rec_Bank,track_Bank,traj_Bank,run_Bank)
            ;
    }
    
    public boolean pion_test(int particle_Index, int pid, float vz, double trigger_electron_vz, HipoDataBank rec_Bank, 
            HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
//            && p > 1.25
//            && p < 5.00 // this wasn't used in the dihadron publication but was used in the submitted single pion
            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && generic_tests.vertex_cut(particle_Index, trigger_electron_vz, rec_Bank, run_Bank) 
            && hadron_pass2_cut(particle_Index, rec_Bank)
//            && pion_chi2pid_cut(particle_Index, rec_Bank)
//            && hadron_chi2pid_cut(particle_Index, rec_Bank)
            && fiducial_cuts.pass1_dc_fiducial_cut(particle_Index,rec_Bank,track_Bank,traj_Bank,run_Bank)
              ;
    }
    
    public boolean kaon_test(int particle_Index, int pid, float vz, double trigger_electron_vz, HipoDataBank rec_Bank, 
            HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
//            && p > 1.00
//            && p < 3.5 
            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && generic_tests.vertex_cut(particle_Index, trigger_electron_vz, rec_Bank, run_Bank) 
            && hadron_pass2_cut(particle_Index, rec_Bank)
//            && hadron_chi2pid_cut(particle_Index, rec_Bank)
            && fiducial_cuts.pass1_dc_fiducial_cut(particle_Index,rec_Bank,track_Bank,traj_Bank,run_Bank)
              ;
    }
    
    public boolean proton_test(int particle_Index, int pid, float vz, double trigger_electron_vz, 
            HipoDataBank rec_Bank, HipoDataBank cal_Bank, HipoDataBank track_Bank, 
            HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
//            && p > 0.4
            && generic_tests.vertex_cut(particle_Index, trigger_electron_vz, rec_Bank, run_Bank) 
            && generic_tests.forward_detector_cut(particle_Index, rec_Bank)
            && fiducial_cuts.pass1_dc_fiducial_cut(particle_Index,rec_Bank,track_Bank,traj_Bank,run_Bank)
            && hadron_pass2_cut(particle_Index, rec_Bank)
//            && hadron_chi2pid_cut(particle_Index, rec_Bank)
              ;
    }
    
    public boolean photon_test(int particle_Index, HipoDataBank rec_Bank, HipoDataBank cal_Bank, 
            LorentzVector lv_e, LorentzVector lv_gamma) {
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
            && p > 0.50
//            && forward_detector_cut(particle_Index, rec_Bank)
//            && pcal_fiducial_cut(particle_Index, rec_Bank, cal_Bank)
//            && random_photon_cut(particle_Index, rec_Bank)
            && open_angle_cut(lv_e, lv_gamma);
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
            
            double trigger_electron_vz = -99;
            double pion_vz = -99;
            double p_max = 0;
            int p_max_index = -99; // find the index of the highest energy electron before FD cut
            // cycle over the particles in recBank and investigate electron and pion IDs
            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
                int pid = rec_Bank.getInt("pid", particle_Index);
                // find electron candidates assigned by EventBuilder
                if (pid!=11) { continue; } // we are only investigating electrons at this point
                float px = rec_Bank.getFloat("px", particle_Index);
                float py = rec_Bank.getFloat("py", particle_Index);
                float pz = rec_Bank.getFloat("pz", particle_Index);
                double p = Math.sqrt(px*px+py*py+pz*pz);

                float vz = rec_Bank.getFloat("vz", particle_Index);

                // searching for the highest momentum electron to use as the "trigger electron"
                if (p > p_max) {
                    p_max = p;
                    p_max_index = particle_Index;
                    trigger_electron_vz = vz;
                }
            }
            

            LorentzVector lv_e = new LorentzVector(); 
            if (p_max_index >= 0) { 
                int pid = rec_Bank.getInt("pid", p_max_index);
                // require that the highest momentum electron be in the forward detector 
                // THIS MAY BE MODIFIED IN A FUTURE ANALYSIS
                float px = rec_Bank.getFloat("px", p_max_index);
                float py = rec_Bank.getFloat("py", p_max_index);
                float pz = rec_Bank.getFloat("pz", p_max_index);
                double p = Math.sqrt(px*px+py*py+pz*pz);
                float vx = rec_Bank.getFloat("vx",p_max_index);
                float vy = rec_Bank.getFloat("vy",p_max_index);
                float vz = rec_Bank.getFloat("vz",p_max_index);
                if (electron_test(p_max_index, p, vz, trigger_electron_vz, rec_Bank, 
                    cal_Bank, track_Bank, traj_Bank,  run_Bank, cc_Bank)) {
                    // this checks all of the PID requirements, if it passes all of them the electron is 
                    // added to the event below
//                    double fe = EnergyLoss(run_Bank.getFloat("torus", 0), pid, px, py, pz);
                   double fe = 1;
                    Particle part = new Particle(pid,fe*px,fe*py,fe*pz,vx,vy,vz);
                    physEvent.addParticle(part);
                    lv_e.setPxPyPzM(px, py, pz, 0.0005109989461);
                }
            }

            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
                int pid = rec_Bank.getInt("pid", particle_Index);
                
                if ((Math.abs(pid)!=211) || trigger_electron_vz == -99) { continue; }
                // requires the particle to be pion by EventBuilder and for an electron to have been assigned to event
                // if no electron was assigned we just skip
                
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", particle_Index);
                float py = rec_Bank.getFloat("py", particle_Index);
                float pz = rec_Bank.getFloat("pz", particle_Index);
                float vx = rec_Bank.getFloat("vx",particle_Index);
                float vy = rec_Bank.getFloat("vy",particle_Index);
                float vz = rec_Bank.getFloat("vz",particle_Index);
                pion_vz = vz;
                if (pion_test(particle_Index, pid, vz, trigger_electron_vz, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                    // check for pion PID
                    
//                   double fe = EnergyLoss(run_Bank.getFloat("torus", 0), pid, px, py, pz);
                   double fe = 1;
                   Particle part = new Particle(pid,fe*px,fe*py,fe*pz,vx,vy,vz);
                   physEvent.addParticle(part);   
                }
//                if (event.hasBank("RICH::Particle")) {
//                    HipoDataBank rich_Bank = (HipoDataBank) event.getBank("RICH::Particle");
//                    pid = rich_detector_pid(particle_Index, rich_Bank);
//                    if (pid != 0) {
//                        Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
//                        physEvent.addParticle(part);
//                    } 
//                }
            }
            
            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
                int pid = rec_Bank.getInt("pid", particle_Index);
                
                if ((Math.abs(pid)!=321) || trigger_electron_vz == -99) { continue; }
                // requires the particle to be pion by EventBuilder and 
                // for an electron to have been assigned to event
                // if no electron was assigned we just skip
                
//                if (pid == 321) { pid = 211; }
//                if (pid == -321) { pid = 211; }
                
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", particle_Index);
                float py = rec_Bank.getFloat("py", particle_Index);
                float pz = rec_Bank.getFloat("pz", particle_Index);
                float vx = rec_Bank.getFloat("vx",particle_Index);
                float vy = rec_Bank.getFloat("vy",particle_Index);
                float vz = rec_Bank.getFloat("vz",particle_Index);
                pion_vz = vz;
                if (kaon_test(particle_Index, pid, vz, trigger_electron_vz, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                    // check for pion PID
                   
                   double fe = EnergyLoss(run_Bank.getFloat("torus", 0), pid, px, py, pz);
//                   double fe = 1;
                   Particle part = new Particle(pid,fe*px,fe*py,fe*pz,vx,vy,vz);
                   physEvent.addParticle(part);   
                }
//                if (event.hasBank("RICH::Particle")) {
//                    HipoDataBank rich_Bank = (HipoDataBank) event.getBank("RICH::Particle");
//                    pid = rich_detector_pid(particle_Index, rich_Bank);
//                    if (pid != 0) {
//                        Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
//                        physEvent.addParticle(part);
//                    } 
//                }
            }
            
            
            for (int particle_Index = 0; particle_Index < rec_Bank.rows(); particle_Index++) {
                int pid = rec_Bank.getInt("pid", particle_Index);
                
                if (pid!=2212 || trigger_electron_vz == -99) { continue; }
                // requires the particle to be proton by EventBuilder & for an electron to have been assigned to event
                
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", particle_Index);
                float py = rec_Bank.getFloat("py", particle_Index);
                float pz = rec_Bank.getFloat("pz", particle_Index);
                double p = Math.sqrt(px*px+py*py+pz*pz);
                float vx = rec_Bank.getFloat("vx",particle_Index);
                float vy = rec_Bank.getFloat("vy",particle_Index);
                float vz = rec_Bank.getFloat("vz",particle_Index);
                
//                // Fermi motion compensation
//                int runnum = run_Bank.getInt("run",0);
//                double dp = 0; // scale size for fermi motion (or previously energy loss)
//                if (runnum > 16000) {
//                    dp = 0.002*p*p*Math.exp(-p*p/16000); 
//                }
////                dp = 0;
//                double fe = (dp+p)/p;

                
                if (proton_test(particle_Index, pid, vz, trigger_electron_vz, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                   
//                   double fe = EnergyLoss(run_Bank.getFloat("torus", 0), pid, px, py, pz);
                   double fe = 1;
                   Particle part = new Particle(pid,fe*px,fe*py,fe*pz,vx,vy,vz);
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