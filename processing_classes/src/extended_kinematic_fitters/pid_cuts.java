/**
 *
 * @author Timothy B. Hayward
 */

package extended_kinematic_fitters;

import org.jlab.io.hipo.HipoDataBank;

import org.jlab.clas.physics.*;


public class pid_cuts {

    /*~~~~~~~~~~~~~~~~~ Electrons ~~~~~~~~~~~~~~~~~*/
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
    
    
    /*~~~~~~~~~~~~~~~~~ Charged Hadrons ~~~~~~~~~~~~~~~~~*/
    public double const_plus_exponential(double a, double b, double c, double x0, double p) {
        return a + b * Math.exp(-(p-x0)/c);
    }
    
    public boolean charged_pion_generic_chi2pid_cut(int particle_Index, HipoDataBank rec_Bank) {
        // pass1 chi2pid id for pions, derived by S. Diehl
        
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
    
    public boolean charged_hadron_chi2pid_cut(int particle_Index, HipoDataBank rec_Bank) {
        float chi2pid = rec_Bank.getFloat("chi2pid", particle_Index);
        
        return Math.abs(chi2pid) < 5;
    }
    
    public boolean charged_hadron_pass2_chi2pid_cut(int particle_Index, HipoDataBank rec_Bank) {
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
                return chi2pid > mu-3*sigma && chi2pid < const_plus_exponential(-0.7, mu+3*sigma+0.7, 0.9, 3, p);
            }
        }
        
        if (pid == -321) {
            mu = 0.115;
            sigma = 0.958;
            if (p < 2) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else {
                return chi2pid > const_plus_exponential(1.7, mu-3*sigma-1.7, 1.0, 2.0, p) && chi2pid < mu+3*sigma;
            }
        }
        
        if (pid == 211) {
            mu = -0.067;
            sigma = 0.956;
            if (p < 3.5) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else {
                return chi2pid > mu-3*sigma && chi2pid < const_plus_exponential(-0.55, mu+3*sigma+0.55, 0.55, 3.5, p);
            }
        }
        
        if (pid == 321) {
            mu = 0.082;
            sigma = 0.985;
            if (p < 2) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else if (p > 2 && p < 2.5) {
                return chi2pid > const_plus_exponential(1.2, mu-3*sigma-1.2, 0.6, 2, p) && chi2pid < mu+3*sigma;
            } else {
                return chi2pid > const_plus_exponential(1.2, mu-3*sigma-1.2, 0.6, 2, p) && chi2pid < const_plus_exponential(2.6, mu+3*sigma-2.6, 0.3, 2.5, p);
            }
        }
        
        if (pid == 2212) {
            mu = 0.372;
            sigma = 1.192;
            if (p < 2) {
                return chi2pid > mu-3*sigma && chi2pid < mu+3*sigma; 
            } else {
                return chi2pid > const_plus_exponential(2, mu-3*sigma-2, 0.9, 2, p) && chi2pid < mu+3*sigma;
            }
        }
        
        return false; // not a charged hadron? 
    }
    
    /*~~~~~~~~~~~~~~~~~ Photons ~~~~~~~~~~~~~~~~~*/
    
    public boolean beta_cut(int particle_Index, HipoDataBank rec_Bank) {
        // remove random photons
        float beta = rec_Bank.getFloat("beta", particle_Index);
        return (beta > 0.90 && beta < 1.10);
    }
    
    public boolean e_gamma_open_angle_cut(LorentzVector lv_e, LorentzVector lv_gamma) {
        return 180/Math.PI*Math.acos(lv_e.vect().dot(lv_gamma.vect())/
                (lv_e.vect().mag()*lv_gamma.vect().mag())) > 8; 
    }
}

