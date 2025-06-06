package analyzers;

import org.jlab.clas.physics.PhysicsEvent;
import java.util.List;
import java.util.ArrayList;

/**
 *
 * @author tbhayward
 */
public class BeamEnergy {

    protected double Eb;
    protected double[] beam_percentage = new double[]{0.99995, 0.98985,
        0.97975, 0.96965, 0.95955, 0.94945, 0.93935,
        0.92925, 0.91915, 0.90905, 0.89895, 0.88885, 0.87875, 0.86865,
        0.85855, 0.84845, 0.83835, 0.82825, 0.81815, 0.80805, 0.79795,
        0.78785, 0.77775, 0.76765, 0.75755, 0.74745, 0.73735, 0.72725,
        0.71715, 0.70705, 0.69695, 0.68685, 0.67675, 0.66665, 0.65655,
        0.64645, 0.63635, 0.62625, 0.61615, 0.60605, 0.59595, 0.58585,
        0.57575, 0.56565, 0.55555, 0.54545, 0.53535, 0.52525, 0.51515,
        0.50505, 0.49495, 0.48485, 0.47475, 0.46465, 0.45455, 0.44445,
        0.43435, 0.42425, 0.41415, 0.40405, 0.39395, 0.38385, 0.37375,
        0.36365, 0.35355, 0.34345, 0.33335, 0.32325, 0.31315, 0.30305,
        0.29295, 0.28285, 0.27275, 0.26265, 0.25255, 0.24245, 0.23235,
        0.22225, 0.21215, 0.20205, 0.19195, 0.18185, 0.17175, 0.16165,
        0.15155, 0.14145, 0.13135, 0.12125, 0.11115, 0.10105, 0.09095,
        0.08085, 0.07075, 0.06065, 0.05055, 0.04045, 0.03035, 0.02025,
        0.01015, 0.0000499999};
    protected double[] beam_likelihood = new double[]{0.791947, 0.808145,
        0.825926, 0.838389, 0.847606, 0.855139,
        0.861519, 0.867214, 0.872296, 0.876636, 0.88062, 0.884296, 0.887679,
        0.890929, 0.893896, 0.896875, 0.899571, 0.902198, 0.904849, 0.907113,
        0.909285, 0.911526, 0.913519, 0.915516, 0.917553, 0.919402, 0.921154,
        0.923104, 0.924949, 0.926656, 0.928376, 0.929996, 0.931685, 0.933225,
        0.934817, 0.936259, 0.937816, 0.939315, 0.940875, 0.942379, 0.943863,
        0.945302, 0.946724, 0.948114, 0.949478, 0.950765, 0.952125, 0.953384,
        0.954781, 0.956184, 0.957518, 0.958847, 0.960215, 0.961672, 0.963001,
        0.964353, 0.96572, 0.967056, 0.968406, 0.969727, 0.971283, 0.972702,
        0.974054, 0.975527, 0.976977, 0.97845, 0.980008, 0.981609, 0.983213,
        0.984729, 0.986328, 0.987938, 0.989692, 0.991332, 0.992895, 0.994582,
        0.996264, 0.997781, 0.999169, 0.999838, 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

    public BeamEnergy(PhysicsEvent recEvent, int runnum, boolean isRadiative) {
        // default beam energy set to rga fall 2018
        Eb = 10.6041; // RGA Fall 2018
//        Eb = 6.535; // RGK Fall 2018 6 GeV
//        Eb = 7.546; // RGK Fall 2018 7 GeV
        if (runnum >= 5032 && runnum <= 5666) {
            Eb = 10.6041;
        } // RGA Fall 2018
        else if (runnum >= 2365 && runnum <= 2598) {
            Eb = 2.22193;
        } // "RGA" Spring 2018 engineering
        else if (runnum >= 3030 && runnum <= 3106) {
            Eb = 6.42313;
        } // RGA Spring 2018 6 GeV
        else if (runnum >= 3819 && runnum <= 3862) {
            Eb = 6.42313;
        } // RGA Spring 2018 6 GeV
        else if (runnum >= 3172 && runnum <= 3817) {
            Eb = 10.5940;
        } // RGA Spring 2018 10 GeV
        else if (runnum >= 3863 && runnum <= 4326) {
            Eb = 10.5940;
        } // RGA Spring 2018 10 GeV
        else if (runnum >= 5875 && runnum <= 6000) {
            Eb = 6.535;
        } // RGK Fall 2018 6 GeV
        else if (runnum >= 5674 && runnum <= 5870) {
            Eb = 7.546;
        } // RGK Fall 2018 7 GeV
        else if (runnum >= 6616 && runnum <= 6783) {
            Eb = 10.1998;
        } // RGA Spring 2019
        else if (runnum >= 6120 && runnum <= 6399) {
            Eb = 10.5986;
        } // RGB Spring 2019
        else if (runnum >= 6409 && runnum <= 6604) {
            Eb = 10.1998;
        } // RGB Spring 2019
        else if (runnum >= 11093 && runnum <= 11283) {
            Eb = 10.4096;
        } // RGB Fall 2019
        else if (runnum >= 11284 && runnum <= 11300) {
            Eb = 4.17179;
        } // RGB Fall 2019
        else if (runnum >= 11323 && runnum <= 11571) {
            Eb = 10.3894;
        } // RGB Spring 2020
        else if (runnum >= 16042 && runnum <= 17065) {
            Eb = 10.5473;
        } // RGC Summer 2022
        else if (runnum >= 17067 && runnum <= 17724) {
            Eb = 10.5563;
        } // RGC Fa 2022
        else if (runnum >= 17725 && runnum <= 17811) {
            Eb = 10.5593;
        } // RGC Fa 2022 / Sp 2023

        // determine a minimum energy the electron must have in order to create the particles in the event
        int num_elec = recEvent.countByPid(11); // returns number of electrons
        int num_positrons = recEvent.countByPid(-11); // returns number of positrons
        int num_piplus = recEvent.countByPid(211);
        int num_piminus = recEvent.countByPid(-211);
        int num_pi0 = recEvent.countByPid(111);
        int num_kplus = recEvent.countByPid(321);
        int num_kminus = recEvent.countByPid(-321);
        int num_protons = recEvent.countByPid(2212);
        int num_antiprotons = recEvent.countByPid(-2212);

//        if (isRadiative) {
//            double double_random = Math.random();
//            for (int i = 0; i < beam_likelihood.length; i++) {
//                if (double_random < beam_likelihood[i]) {
//                    Eb = Eb * beam_percentage[i];
//                    break;
//                }
//            }
//        }
        if (isRadiative) {
            // Calculate the total rest mass of the final state particles
            double total_mass = (num_elec-1) * 0.000511 + num_positrons * 0.000511
                    + num_piplus * 0.13957 + num_piminus * 0.13957
                    + num_pi0 * 0.1349766 + num_kplus * 0.493677
                    + num_kminus * 0.493677 + (num_protons-1) * 0.938272
                    + num_antiprotons * 0.938272;

            // Determine the minimum beam percentage allowed
            double min_beam_percentage = total_mass / Eb;

            // Create new arrays for beam percentages and likelihoods that satisfy energy conservation
            List<Double> new_beam_percentage = new ArrayList<>();
            List<Double> new_beam_likelihood = new ArrayList<>();

            for (int i = 0; i < beam_percentage.length; i++) {
                if (beam_percentage[i] >= min_beam_percentage) {
                    new_beam_percentage.add(beam_percentage[i]);
                    new_beam_likelihood.add(beam_likelihood[i]);
                }
            }

            // Check if there are valid beam percentages available
            if (!new_beam_percentage.isEmpty()) {
                // Normalize the new beam likelihoods
                double max_likelihood = new_beam_likelihood.get(new_beam_likelihood.size() - 1);
                for (int i = 0; i < new_beam_likelihood.size(); i++) {
                    new_beam_likelihood.set(i, new_beam_likelihood.get(i) / max_likelihood);
                }

                // Update the beam energy using the adjusted arrays
                double double_random = Math.random();
                for (int i = 0; i < new_beam_likelihood.size(); i++) {
                    if (double_random < new_beam_likelihood.get(i)) {
                        Eb = Eb * new_beam_percentage.get(i);
                        break;
                    }
                }
            } else {
                // If no valid beam percentages, set Eb to total_mass to conserve energy
                Eb = total_mass;
            }
        }
    }

    public double Eb() {
        return Eb;
    }

}
