/**
 *
 * @author Timothy B. Hayward
 */
package extended_kinematic_fitters;

import org.jlab.io.hipo.HipoDataBank;

import org.jlab.clas.physics.*;

public class pid_cuts {

    /*~~~~~~~~~~~~~~~~~ Electrons ~~~~~~~~~~~~~~~~~*/
    public boolean calorimeter_energy_cut(int particle_Index, HipoDataBank cal_Bank, HipoDataBank run_Bank) {

        int runnum = run_Bank.getInt("run", 0);

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
                if (runnum >= 16043 && runnum <= 17811) {
                    return energy > 0.15;
                }
                return energy > 0.07;
            }
        }

        // If no matching rows are found, return false
        return false;
    }

    public boolean calorimeter_sampling_fraction_cut(int particle_Index, double p, HipoDataBank run_Bank,
            HipoDataBank cal_Bank) {
        double cal_energy = 0;
        int sector = -1;
        double max_energy = 0;

        // Sum calorimeter energy & find sector
        for (int i = 0; i < cal_Bank.rows(); i++) {
            if (cal_Bank.getInt("pindex", i) == particle_Index) {
                double e = cal_Bank.getFloat("energy", i);
                cal_energy += e;
                if (e > max_energy) {
                    max_energy = e;
                    sector = cal_Bank.getInt("sector", i);
                }
            }
        }

        double sf = cal_energy / p;
        int runnum = run_Bank.getInt("run", 0);

        // Sector must be valid
        if (sector < 1 || sector > 6) {
            return false;
        }
        int idx = sector - 1;

        // Prepare coefficient arrays
        double[][] lowerCoeffs = new double[6][3];
        double[][] upperCoeffs = new double[6][3];

        // RGA Fall 2018 Inbending
        if (runnum >= 4991 && runnum <= 5419) {
            lowerCoeffs[0] = new double[]{0.167070, 0.010930, -0.001085};
            upperCoeffs[0] = new double[]{0.302129, -0.000735, -0.000409};
            lowerCoeffs[1] = new double[]{0.146344, 0.023547, -0.002660};
            upperCoeffs[1] = new double[]{0.313502, -0.004386, 0.000036};
            lowerCoeffs[2] = new double[]{0.136665, 0.031226, -0.003705};
            upperCoeffs[2] = new double[]{0.310201, -0.001985, -0.000242};
            lowerCoeffs[3] = new double[]{0.152493, 0.018234, -0.001604};
            upperCoeffs[3] = new double[]{0.315967, -0.004037, 0.000058};
            lowerCoeffs[4] = new double[]{0.143086, 0.023846, -0.002510};
            upperCoeffs[4] = new double[]{0.321858, -0.009117, 0.000657};
            lowerCoeffs[5] = new double[]{0.152500, 0.021085, -0.002384};
            upperCoeffs[5] = new double[]{0.313124, -0.003053, -0.000193};

            // RGA Fall 2018 Outbending
        } else if (runnum >= 5423 && runnum <= 5666) {
            lowerCoeffs[0] = new double[]{0.168531, 0.013221, -0.001070};
            upperCoeffs[0] = new double[]{0.305273, -0.002978, -0.000101};
            lowerCoeffs[1] = new double[]{0.166549, 0.013526, -0.001129};
            upperCoeffs[1] = new double[]{0.298919, -0.000799, -0.000253};
            lowerCoeffs[2] = new double[]{0.165926, 0.016098, -0.001352};
            upperCoeffs[2] = new double[]{0.301111, -0.000298, -0.000414};
            lowerCoeffs[3] = new double[]{0.172671, 0.009819, -0.000641};
            upperCoeffs[3] = new double[]{0.303462, -0.002267, -0.000161};
            lowerCoeffs[4] = new double[]{0.166557, 0.010846, -0.000709};
            upperCoeffs[4] = new double[]{0.297885, -0.001656, -0.000110};
            lowerCoeffs[5] = new double[]{0.168393, 0.012321, -0.000936};
            upperCoeffs[5] = new double[]{0.298590, -0.000390, -0.000383};

            // RGA Spring 2019 Inbending
        } else if (runnum >= 6616 && runnum <= 6783) {
            lowerCoeffs[0] = new double[]{0.159840, 0.016237, -0.001594};
            upperCoeffs[0] = new double[]{0.309811, -0.003301, -0.000171};
            lowerCoeffs[1] = new double[]{0.138572, 0.026177, -0.003200};
            upperCoeffs[1] = new double[]{0.306967, -0.000238, -0.000421};
            lowerCoeffs[2] = new double[]{0.123413, 0.038137, -0.004472};
            upperCoeffs[2] = new double[]{0.312705, -0.003361, -0.000150};
            lowerCoeffs[3] = new double[]{0.142108, 0.023657, -0.002314};
            upperCoeffs[3] = new double[]{0.312121, -0.004128, 0.000035};
            lowerCoeffs[4] = new double[]{0.134147, 0.029226, -0.003312};
            upperCoeffs[4] = new double[]{0.307159, -0.001868, -0.000248};
            lowerCoeffs[5] = new double[]{0.144287, 0.024536, -0.002842};
            upperCoeffs[5] = new double[]{0.311596, -0.002988, -0.000228};

            // Monte Carlo
        } else if (runnum == 11) {
            lowerCoeffs[0] = new double[]{0.182342, 0.010612, -0.000989};
            upperCoeffs[0] = new double[]{0.316417, -0.007887, 0.000497};
            lowerCoeffs[1] = new double[]{0.174139, 0.015019, -0.001926};
            upperCoeffs[1] = new double[]{0.319061, -0.009477, 0.000800};
            lowerCoeffs[2] = new double[]{0.168963, 0.017387, -0.002206};
            upperCoeffs[2] = new double[]{0.318202, -0.008970, 0.000741};
            lowerCoeffs[3] = new double[]{0.176921, 0.012012, -0.001191};
            upperCoeffs[3] = new double[]{0.315149, -0.008103, 0.000559};
            lowerCoeffs[4] = new double[]{0.177276, 0.013250, -0.001440};
            upperCoeffs[4] = new double[]{0.316444, -0.008119, 0.000571};
            lowerCoeffs[5] = new double[]{0.180294, 0.011205, -0.001183};
            upperCoeffs[5] = new double[]{0.317384, -0.008433, 0.000605};

            // New Su22 (16043–16772)
        } else if (runnum >= 16043 && runnum <= 16772) {
            lowerCoeffs[0] = new double[]{0.193134, 0.005244, -0.000566};
            upperCoeffs[0] = new double[]{0.302522, -0.005450, 0.000329};
            lowerCoeffs[1] = new double[]{0.189430, 0.006233, -0.000609};
            upperCoeffs[1] = new double[]{0.303987, -0.000174, -0.000490};
            lowerCoeffs[2] = new double[]{0.189653, 0.007051, -0.000707};
            upperCoeffs[2] = new double[]{0.302645, 0.002085, -0.000874};
            lowerCoeffs[3] = new double[]{0.178332, 0.010998, -0.001015};
            upperCoeffs[3] = new double[]{0.307951, -0.001060, -0.000424};
            lowerCoeffs[4] = new double[]{0.181131, 0.008697, -0.000784};
            upperCoeffs[4] = new double[]{0.310068, -0.003615, -0.000236};
            lowerCoeffs[5] = new double[]{0.185641, 0.008087, -0.000799};
            upperCoeffs[5] = new double[]{0.304876, 0.000195, -0.000638};

            // New Fa22 (16843–17408)
        } else if (runnum >= 16843 && runnum <= 17408) {
            lowerCoeffs[0] = new double[]{0.176355, 0.011013, -0.001184};
            upperCoeffs[0] = new double[]{0.321309, -0.012014, 0.001052};
            lowerCoeffs[1] = new double[]{0.183261, 0.007594, -0.000720};
            upperCoeffs[1] = new double[]{0.310003, -0.003690, -0.000113};
            lowerCoeffs[2] = new double[]{0.193621, 0.004625, -0.000448};
            upperCoeffs[2] = new double[]{0.309858, -0.002373, -0.000340};
            lowerCoeffs[3] = new double[]{0.185285, 0.007387, -0.000686};
            upperCoeffs[3] = new double[]{0.314824, -0.002479, -0.000235};
            lowerCoeffs[4] = new double[]{0.184573, 0.006317, -0.000532};
            upperCoeffs[4] = new double[]{0.312326, -0.005413, 0.000004};
            lowerCoeffs[5] = new double[]{0.182386, 0.008435, -0.000777};
            upperCoeffs[5] = new double[]{0.316822, -0.006715, 0.000075};

            // New Sp23 (17477–17811)
        } else if (runnum >= 17477 && runnum <= 17811) {
            lowerCoeffs[0] = new double[]{0.190334, 0.004279, -0.000448};
            upperCoeffs[0] = new double[]{0.307456, -0.006925, 0.000473};
            lowerCoeffs[1] = new double[]{0.185558, 0.007610, -0.000723};
            upperCoeffs[1] = new double[]{0.309016, -0.002613, -0.000219};
            lowerCoeffs[2] = new double[]{0.187881, 0.007244, -0.000712};
            upperCoeffs[2] = new double[]{0.309075, -0.004073, -0.000147};
            lowerCoeffs[3] = new double[]{0.179524, 0.010042, -0.000902};
            upperCoeffs[3] = new double[]{0.314158, -0.003706, -0.000116};
            lowerCoeffs[4] = new double[]{0.182893, 0.006639, -0.000552};
            upperCoeffs[4] = new double[]{0.311758, -0.006891, 0.000130};
            lowerCoeffs[5] = new double[]{0.188799, 0.004953, -0.000454};
            upperCoeffs[5] = new double[]{0.309367, -0.004730, -0.000132};

        } else {
            return sf > 0.19;
        }

        // Compute the quadratic limits
        double[] lo = lowerCoeffs[idx];
        double[] hi = upperCoeffs[idx];
        double lower_limit = lo[0] + lo[1] * p + lo[2] * p * p;
        double upper_limit = hi[0] + hi[1] * p + hi[2] * p * p;

        // Su22 sector 1 extra cap at 0.28
        if (runnum >= 16043 && runnum <= 16772 && idx == 0) {
            upper_limit = Math.min(upper_limit, 0.28);
        }

        // Final cut
        return sf > lower_limit && sf < upper_limit;
    }

    public boolean pass_1_calorimeter_sampling_fraction_cut(int particle_Index, double p, HipoDataBank run_Bank,
            HipoDataBank cal_Bank) {
        double scale = 3.5; // how many std away from mean to cut on
        int sector = -1;
        double cal_energy = 0;
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row) == particle_Index) {
                sector = cal_Bank.getInt("sector", current_Row) - 1; // subtract one to start at index of 0
                cal_energy += cal_Bank.getFloat("energy", current_Row);
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
        if ((runnum >= 11323 && runnum <= 11571) || (runnum >= 11093 && runnum <= 11300)) {
            // RGB winter 2020 // (also using this for RGB fall 2019, but it should be updated! TODO)

            double mean = e_cal_sampl_mu_RGBW20[0][sector] + (e_cal_sampl_mu_RGBW20[1][sector] / 1000)
                    * (p - e_cal_sampl_mu_RGBW20[2][sector]) * (p - e_cal_sampl_mu_RGBW20[2][sector]);

            double std = e_cal_sampl_sigmaRGBW20[0][sector] + e_cal_sampl_sigmaRGBW20[1][sector]
                    / (10 * (p - e_cal_sampl_sigmaRGBW20[2][sector]));

            return ((cal_energy / p) > (mean - scale * std)) && ((cal_energy / p) < (mean + scale * std));
        } else if (runnum >= 6120 && runnum <= 6604) { // RGB Sp19

            double mean = e_cal_sampl_mu_RGBSp19[0][sector] + (e_cal_sampl_mu_RGBSp19[1][sector] / 1000)
                    * (p - e_cal_sampl_mu_RGBSp19[2][sector]) * (p - e_cal_sampl_mu_RGBSp19[2][sector]);

            double std = e_cal_sampl_sigma_RGBSp19[0][sector] + e_cal_sampl_sigma_RGBSp19[1][sector]
                    / (10 * (p - e_cal_sampl_sigma_RGBSp19[2][sector]));

            return ((cal_energy / p) > (mean - scale * std)) && ((cal_energy / p) < (mean + scale * std));
        } else {
            double mean = e_cal_sampl_mu_RGA[0][sector] + (e_cal_sampl_mu_RGA[1][sector] / 1000)
                    * (p - e_cal_sampl_mu_RGA[2][sector]) * (p - e_cal_sampl_mu_RGA[2][sector]);

            double std = e_cal_sampl_sigma_RGA[0][sector] + e_cal_sampl_sigma_RGA[1][sector]
                    / (10 * (p - e_cal_sampl_sigma_RGA[2][sector]));

            return ((cal_energy / p) > (mean - scale * std)) && ((cal_energy / p) < (mean + scale * std));
        }
    }

    public boolean calorimeter_diagonal_cut(int particle_Index, double p, HipoDataBank cal_Bank, HipoDataBank run_Bank) {
        // Only apply diagonal cut above 4.5 GeV
        if (p < 4.9) {
            return true;
        }

        int runnum = run_Bank.getInt("run", 0);
        if (!(runnum >= 16043 && runnum <= 17811)) {

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
            return 0.19 < pcal_plus_ecal_inner / p;
        } else {
            double pcal_energy = 0;
            double ecin_energy = 0;
            // Iterate through the rows of the data bank
            for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
                int pindex = cal_Bank.getInt("pindex", current_Row);
                int layer = cal_Bank.getInt("layer", current_Row);
                if (pindex == particle_Index && layer == 1) {
                    pcal_energy += cal_Bank.getFloat("energy", current_Row);
                } else if (pindex == particle_Index && layer == 4) {
                    ecin_energy += cal_Bank.getFloat("energy", current_Row);
                }
            }
            return (ecin_energy / p) >= -0.625 * (pcal_energy / p) + 0.15;
        }
    }

    /*~~~~~~~~~~~~~~~~~ Charged Hadrons ~~~~~~~~~~~~~~~~~*/
    public double const_plus_exponential(double a, double b, double c, double x0, double p) {
        return a + b * Math.exp(-(p - x0) / c);
    }

    public boolean pass_1_pion_generic_chi2pid_cut(int particle_Index, HipoDataBank rec_Bank) {
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

    public boolean charged_hadron_chi2pid_cut(int particle_Index, HipoDataBank rec_Bank, HipoDataBank run_Bank) {

        int pid = rec_Bank.getInt("pid", particle_Index);
        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(px * px + py * py + pz * pz);
        float chi2pid = rec_Bank.getFloat("chi2pid", particle_Index);

        int runnum = run_Bank.getInt("run", 0);

        generic_tests generic_tests = new generic_tests();
        boolean isForwardDetector = generic_tests.forward_detector_cut(particle_Index, rec_Bank);
        boolean isCentralDetector = generic_tests.central_detector_cut(particle_Index, rec_Bank);

        if (isForwardDetector) {
            if (p >= 5.0) {
                return false;
            }
            if (runnum >= 6616 && runnum <= 6783) {
                if (pid == 211) {
                    return -0.05 - 1.08 * 3.0 < chi2pid && chi2pid < -0.05 + 1.08 * 3.0;
                }
                if (pid == -211) {
                    return -0.02 - 1.08 * 3.0 < chi2pid && chi2pid < -0.02 + 1.08 * 3.0;
                }
                if (pid == 321) {
                    return 0.03 - 1.08 * 3.0 < chi2pid && chi2pid < 0.03 + 1.08 * 3.0;
                }
                if (pid == -321) {
                    return -0.14 - 1.30 * 3.0 < chi2pid && chi2pid < -0.14 + 1.30 * 3.0;
                }
                if (pid == 2212) {
                    return 0.36 - 1.36 * 3.5 < chi2pid && chi2pid < 0.36 + 1.36 * 3.5;
                }
            } else if (runnum == 11) { // MC
                if (pid == 211) {
                    return 0.06 - 1.23 * 3.0 < chi2pid && chi2pid < 0.06 + 1.23 * 3.0;
                }
                if (pid == -211) {
                    return -0.03 - 1.26 * 3.0 < chi2pid && chi2pid < -0.03 + 1.26 * 3.0;
                }
                if (pid == 321) {
                    return -0.02 - 1.31 * 3.0 < chi2pid && chi2pid < -0.02 + 1.31 * 3.0;
                }
                if (pid == -321) {
                    return -0.04 - 1.04 * 3.0 < chi2pid && chi2pid < -0.04 + 1.04 * 3.0;
                }
                if (pid == 2212) {
                    return 0.07 - 1.35 * 3.5 < chi2pid && chi2pid < 0.07 + 1.35 * 3.5;
                }
            }
        }
        if (isCentralDetector) {
            if (pid == 2212 && p >= 1.5) {
                return false;
            }
            if (p < 0.25) {
                return false;
            }
            if (runnum >= 6616 && runnum <= 6783) {
                if (pid == 211) {
                    return -0.18 - 1.76 * 3.0 < chi2pid && chi2pid < -0.18 + 1.76 * 3.0;
                }
                if (pid == -211) {
                    return -0.21 - 1.68 * 3.0 < chi2pid && chi2pid < -0.21 + 1.68 * 3.0;
                }
                if (pid == 321) {
                    return 0.59 - 2.06 * 3.0 < chi2pid && chi2pid < 0.59 + 2.06 * 3.0;
                }
                if (pid == -321) {
                    return -0.20 - 1.72 * 3.0 < chi2pid && chi2pid < -0.20 + 1.72 * 3.0;
                }
                if (pid == 2212) {
                    return 0.75 - 2.15 * 3.5 < chi2pid && chi2pid < 0.75 + 2.15 * 3.5;
                }
            } else if (runnum == 11) { // MC
                if (pid == 211) {
                    return 0.09 - 1.57 * 3.0 < chi2pid && chi2pid < 0.09 + 1.57 * 3.0;
                }
                if (pid == -211) {
                    return -0.02 - 1.52 * 3.0 < chi2pid && chi2pid < -0.02 + 1.52 * 3.0;
                }
                if (pid == 321) {
                    return -0.14 - 1.75 * 3.0 < chi2pid && chi2pid < -0.14 + 1.75 * 3.0;
                }
                if (pid == -321) {
                    return -0.37 - 1.44 * 3.0 < chi2pid && chi2pid < -0.37 + 1.44 * 3.0;
                }
                if (pid == 2212) {
                    return 0.45 - 1.91 * 3.5 < chi2pid && chi2pid < 0.45 + 1.91 * 3.5;
                }
            }
        }

        return Math.abs(chi2pid) < 4.0;
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
        double sigma = 1;

        if (pid == -211) {
            mu = -0.063;
            sigma = 0.947;
            if (p < 3) {
                return chi2pid > mu - 3 * sigma && chi2pid < mu + 3 * sigma;
            } else {
                return chi2pid > mu - 3 * sigma && chi2pid < const_plus_exponential(-0.7, mu + 3 * sigma + 0.7, 0.9, 3, p);
            }
        }

        if (pid == -321) {
            mu = 0.115;
            sigma = 0.958;
            if (p < 2) {
                return chi2pid > mu - 3 * sigma && chi2pid < mu + 3 * sigma;
            } else {
                return chi2pid > const_plus_exponential(1.7, mu - 3 * sigma - 1.7, 1.0, 2.0, p) && chi2pid < mu + 3 * sigma;
            }
        }

        if (pid == 211) {
            mu = -0.067;
            sigma = 0.956;
            if (p < 3.5) {
                return chi2pid > mu - 3 * sigma && chi2pid < mu + 3 * sigma;
            } else {
                return chi2pid > mu - 3 * sigma && chi2pid < const_plus_exponential(-0.55, mu + 3 * sigma + 0.55, 0.55, 3.5, p);
            }
        }

        if (pid == 321) {
            mu = 0.082;
            sigma = 0.985;
            if (p < 2) {
                return chi2pid > mu - 3 * sigma && chi2pid < mu + 3 * sigma;
            } else if (p > 2 && p < 2.5) {
                return chi2pid > const_plus_exponential(1.2, mu - 3 * sigma - 1.2, 0.6, 2, p) && chi2pid < mu + 3 * sigma;
            } else {
                return chi2pid > const_plus_exponential(1.2, mu - 3 * sigma - 1.2, 0.6, 2, p) && chi2pid < const_plus_exponential(2.6, mu + 3 * sigma - 2.6, 0.3, 2.5, p);
            }
        }

        if (pid == 2212) {
            mu = 0.372;
            sigma = 1.192;
            if (p < 2) {
                return chi2pid > mu - 3 * sigma && chi2pid < mu + 3 * sigma;
            } else {
                return chi2pid > const_plus_exponential(2, mu - 3 * sigma - 2, 0.9, 2, p) && chi2pid < mu + 3 * sigma;
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

}
