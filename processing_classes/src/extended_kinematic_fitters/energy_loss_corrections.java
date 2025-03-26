package extended_kinematic_fitters;

import org.jlab.io.hipo.HipoDataBank;

/**
 *
 * @author tbhayward
 */
public class energy_loss_corrections {

    public static double p_calculation(double x, double y, double z) {
        return Math.sqrt(x * x + y * y + z * z);
    }

    public static double phi_calculation(double x, double y) {
        // tracks are given with Cartesian values and so must be converted to cylindrical
        double phi = Math.toDegrees(Math.atan2(x, y));
        phi = phi - 90;
        if (phi < 0) {
            phi = 360 + phi;
        }
        phi = 360 - phi;
        return phi;
    }

    public static double theta_calculation(double x, double y, double z) {
        // convert cartesian coordinates to polar angle
        double r = Math.pow(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2), 0.5);
        return (double) (180 / Math.PI) * Math.acos(z / r);
    }

    public static double x_calculation(double p, double theta, double phi) {
        // Convert angles from degrees to radians
        theta = Math.toRadians(theta);
        phi = Math.toRadians(phi);

        // Calculate x using spherical to cartesian conversion
        return p * Math.sin(theta) * Math.cos(phi);
    }

    public static double y_calculation(double p, double theta, double phi) {
        // Convert angles from degrees to radians
        theta = Math.toRadians(theta);
        phi = Math.toRadians(phi);

        // Calculate y using spherical to cartesian conversion
        return p * Math.sin(theta) * Math.sin(phi);
    }

    public static double z_calculation(double p, double theta) {
        // Convert angle from degrees to radians
        theta = Math.toRadians(theta);

        // Calculate z using spherical to cartesian conversion
        return p * Math.cos(theta);
    }

    public void proton_energy_loss_corrections(int particle_Index, float[] p_array,
            HipoDataBank rec_Bank, HipoDataBank run_Bank) {

        double px = p_array[0];
        double py = p_array[1];
        double pz = p_array[2];

        double p = p_calculation(px, py, pz);
        double theta = theta_calculation(px, py, pz);
        double phi = phi_calculation(px, py);

        generic_tests generic_tests = new generic_tests();
        boolean isForwardDetector = generic_tests.forward_detector_cut(particle_Index, rec_Bank);
        boolean isCentralDetector = generic_tests.central_detector_cut(particle_Index, rec_Bank);

        int runnum = run_Bank.getInt("run", 0);

        double dp = 0;
        double A_p = 0;
        double B_p = 0;
        double C_p = 0;
        double A_theta = 0;
        double B_theta = 0;
        double C_theta = 0;
        double A_phi = 0;
        double B_phi = 0;
        double C_phi = 0;
        if (runnum >= 4763 && runnum <= 5419) { // RGA Fa18 Inb

            if (isForwardDetector && !isCentralDetector) {

                A_p = 0.0099626 - 0.0002414 * theta - 0.0000020 * theta * theta;
                B_p = -0.01428267 + 0.00042833 * theta + 0.00001081 * theta * theta;
                C_p = 0.01197102 - 0.00055673 * theta + 0.00000785 * theta * theta;

                A_theta = 0.0683831 - 0.0083821 * theta + 0.0001670 * theta * theta;
                B_theta = -0.15834256 + 0.02630760 * theta - 0.00064126 * theta * theta;
                C_theta = 0.11587509 - 0.01679559 * theta + 0.00038915 * theta * theta;

                A_phi = 0.0416510 - 0.0064212 * theta + 0.0000622 * theta * theta;
                B_phi = 0.28414191 - 0.00047647 * theta + 0.00010357 * theta * theta;
                C_phi = -0.25690893 + 0.00886707 * theta - 0.00016081 * theta * theta;

            } else if (!isForwardDetector && isCentralDetector) {

                A_p = -0.2383991 + 0.0124992 * theta - 0.0001646 * theta * theta;
                B_p = 0.60123885 - 0.03128464 * theta + 0.00041314 * theta * theta;
                C_p = -0.44080146 + 0.02209857 * theta - 0.00028224 * theta * theta;

                A_theta = 0.1000890 - 0.0039222 * theta + 0.0000359 * theta * theta;
                B_theta = -0.0130680 + 0.0004545 * theta - 0.0000026 * theta * theta;
                C_theta = 0;

                A_phi = 0.0776934 - 0.0059632 * theta + 0.0000749 * theta * theta;
                B_phi = -0.31582008 + 0.01649220 * theta - 0.00018505 * theta * theta;
                C_phi = 0.10909746 - 0.00530642 * theta + 0.00005627 * theta * theta;
            }

        } else if (runnum >= 5423 && runnum <= 5666) { // RGA Fa18 Out

            if (isForwardDetector && !isCentralDetector) {

                A_p = 0.0135790 - 0.0005303 * theta;
                B_p = -0.02165929 + 0.00121123 * theta;
                C_p = 0.0;

                // A_theta, B_theta, C_theta
                A_theta = -0.3715486 + 0.0272810 * theta - 0.0006278 * theta * theta + 0.0000040 * theta * theta * theta;
                B_theta = 2.00009939 - 0.20781779 * theta + 0.00721092 * theta * theta - 0.00008343 * theta * theta * theta;
                C_theta = 0;

                // A_phi, B_phi, C_phi
                A_phi = -0.9701486 + 0.1213124 * theta - 0.0049215 * theta * theta + 0.0000640 * theta * theta * theta;
                B_phi = 2.85034691 - 0.34405076 * theta + 0.01347377 * theta * theta - 0.00016663 * theta * theta * theta;
                C_phi = 0;

            } else if (!isForwardDetector && isCentralDetector) {

                // A_p, B_p, C_p
                A_p = -0.1927861 + 0.0099546 * theta - 0.0001299 * theta * theta;
                B_p = 0.44307822 - 0.02309469 * theta + 0.00030784 * theta * theta;
                C_p = -0.32938000 + 0.01648659 * theta - 0.00021181 * theta * theta;

                // A_theta, B_theta, C_theta
                A_theta = 0.0581473 - 0.0021818 * theta + 0.0000181 * theta * theta;
                B_theta = 0.00915748 - 0.00040748 * theta + 0.00000562 * theta * theta;
                C_theta = 0; // No C_theta for rga_fa18_out

                // A_phi, B_phi, C_phi
                A_phi = -0.0733814 + 0.0010335 * theta - 0.0000044 * theta * theta;
                B_phi = -0.06127800 + 0.00492239 * theta - 0.00005683 * theta * theta;
                C_phi = 0.02586507 - 0.00160176 * theta + 0.00001642 * theta * theta;

            }

        } else if (runnum >= 6616 && runnum <= 6783) { // RGA Sp19 Inb
            if (isForwardDetector && !isCentralDetector) {

                // A_p, B_p, C_p
                A_p = 0.0095205 - 0.0001914 * theta - 0.0000031 * theta * theta;
                B_p = -0.01365658 + 0.00036322 * theta + 0.00001217 * theta * theta;
                C_p = 0.01175256 - 0.00053407 * theta + 0.00000742 * theta * theta;

                // A_theta, B_theta, C_theta
                A_theta = 0.0723069 - 0.0085078 * theta + 0.0001702 * theta * theta;
                B_theta = -0.16048057 + 0.02561073 * theta - 0.00062158 * theta * theta;
                C_theta = 0.10954630 - 0.01566605 * theta + 0.00036132 * theta * theta;

                // A_phi, B_phi, C_phi
                A_phi = 0.0486986 - 0.0067579 * theta + 0.0000638 * theta * theta;
                B_phi = 0.26803189 + 0.00016245 * theta + 0.00010433 * theta * theta;
                C_phi = -0.24522460 + 0.00826646 * theta - 0.00015640 * theta * theta;

            } else if (!isForwardDetector && isCentralDetector) {

                // A_p, B_p (no C_p for rga_fa18_out and no theta^2 term)
                A_p = -0.2716918 + 0.0142491 * theta - 0.0001862 * theta * theta;
                B_p = 0.65945101 - 0.03431360 * theta + 0.00045036 * theta * theta;
                C_p = -0.46602726 + 0.02335623 * theta - 0.00029720 * theta * theta;

                // A_theta, B_theta, C_theta
                A_theta = 0.2550377 - 0.0107983 * theta + 0.0001116 * theta * theta;
                B_theta = -0.14022533 + 0.00596067 * theta - 0.00006172 * theta * theta;
                C_theta = 0;

                // A_phi, B_phi, C_phi
                A_phi = -0.5459156 + 0.0219868 * theta - 0.0002349 * theta * theta;
                B_phi = 0.74223687 - 0.03037065 * theta + 0.00032761 * theta * theta;
                C_phi = -0.29798258 + 0.01246744 * theta - 0.00013525 * theta * theta;

            }

        } else if (runnum >= 16089 && runnum <= 16786) { // RGC Su22 Inb

            if (isForwardDetector && !isCentralDetector) {

                A_p = 0.0109317 - 0.0000194 * theta - 0.0000117 * theta * theta;
                B_p = -0.00910576 - 0.00035154 * theta + 0.00003905 * theta * theta;
                C_p = 0.01225782 - 0.00012805 * theta - 0.00000820 * theta * theta;

                A_theta = 0.0644813 - 0.0079393 * theta + 0.0001566 * theta * theta;
                B_theta = -0.13787609 + 0.02395150 * theta - 0.00058811 * theta * theta;
                C_theta = 0.10551548 - 0.01569699 * theta + 0.00036501 * theta * theta;

                A_phi = 0.0787287 - 0.0075095 * theta + 0.0000669 * theta * theta;
                B_phi = 0.03705727 + 0.01332536 * theta - 0.00009908 * theta * theta;
                C_phi = -0.10680417 - 0.00141926 * theta + 0.00001672 * theta * theta;

            } else if (!isForwardDetector && isCentralDetector) {

                A_p = -0.3951652 + 0.0202840 * theta - 0.0002660 * theta * theta;
                B_p = 0.93238668 - 0.04803619 * theta + 0.00063215 * theta * theta;
                C_p = -0.59146847 + 0.02997697 * theta - 0.00038773 * theta * theta;

                // A_theta, B_theta, C_theta
                A_theta = 0.0644813 - 0.0079393 * theta + 0.0001566 * theta * theta;
                B_theta = -0.13787609 + 0.02395150 * theta - 0.00058811 * theta * theta;
                C_theta = 0.10551548 - 0.01569699 * theta + 0.00036501 * theta * theta;

                // A_phi, B_phi, C_phi
                A_phi = 0.0787287 - 0.0075095 * theta + 0.0000669 * theta * theta;
                B_phi = 0.03705727 + 0.01332536 * theta - 0.00009908 * theta * theta;
                C_phi = -0.10680417 - 0.00141926 * theta + 0.00001672 * theta * theta;
            }

        }

        if (isForwardDetector && !isCentralDetector) {
            p += A_p + B_p / p + C_p / (p * p);
            theta += A_theta + B_theta / theta + C_theta / (theta * theta);
            phi += A_phi + B_phi / phi + C_phi / (phi * phi);
        } else if (!isForwardDetector && isCentralDetector) {
            p += A_p + B_p * p + C_p * p * p;
            theta += A_theta + B_theta / theta + C_theta / (theta * theta);
            phi += A_phi + B_phi / phi + C_phi / (phi * phi);
        }

        // Update the px, py, pz values
        p_array[0] = (float) x_calculation(p, theta, phi);
        p_array[1] = (float) y_calculation(p, theta, phi);
        p_array[2] = (float) z_calculation(p, theta);
    }
    
    public static void electron_energy_loss_corrections(int particle_Index, float[] p_array,
            HipoDataBank rec_Bank, HipoDataBank run_Bank) {

        double px = p_array[0];
        double py = p_array[1];
        double pz = p_array[2];

        double p = Math.sqrt(px * px + py * py + pz * pz);
        double theta = Math.acos(pz / p); // in radians
        double phi = Math.atan2(py, px); // in radians

        generic_tests genericTests = new generic_tests();
        boolean isForwardTagger = genericTests.forward_tagger_cut(particle_Index, rec_Bank);

        int runnum = run_Bank.getInt("run", 0);

        // Check if run is in allowed ranges
        boolean validRun = (runnum >= 4763 && runnum <= 5419) || 
                          (runnum >= 5423 && runnum <= 5666) || 
                          (runnum >= 6616 && runnum <= 6783);
        if (!validRun) {
            return;
        }

        double mass = 0.00511; // Electron mass in GeV
        double E = Math.sqrt(p * p + mass * mass);
        double E_new = E;
        double theta_new = theta;

        if (isForwardTagger) {
            // Apply spring19 fit correction for forward tagger
            E_new = E + 0.085643 - 0.0288063 * E + 0.00894691 * E * E - 0.000725449 * E * E * E;
        } else {
            // Apply Richard Tyson's correction for other cases
            double[] corOut_PDep = {-6.520e-02, 7.099e-03, -5.929e-05, 2.145e-01, -1.153e-01};
            double[] corIn_PDep = {-9.538e-03,6.661e-03,-3.333e-04,-6.136e-03,3.611e-02};
            double delta = corOut_PDep[0] + corOut_PDep[1] * E + corOut_PDep[2] * E * E 
                         + corOut_PDep[3] / E + corOut_PDep[4] / (E * E);
            E_new = E + delta * E;

            // Theta correction
            double p0 = -0.08783;
            double p1 = 73.43;
            double p2 = 0.03231;

            double phi_deg = Math.toDegrees(phi);
            double theta_correction_deg = p0 * Math.sin(Math.toRadians(p1 + phi_deg)) + p2;
            theta_new += Math.toRadians(theta_correction_deg);
        }

        // Calculate new momentum magnitude (approximated as E_new)
        double new_p = (E_new);

        // Calculate new components based on theta_new and original phi
        double px_new = new_p * Math.sin(theta_new) * Math.cos(phi);
        double py_new = new_p * Math.sin(theta_new) * Math.sin(phi);
        double pz_new = new_p * Math.cos(theta_new);

        // Update the momentum array
        p_array[0] = (float) px_new;
        p_array[1] = (float) py_new;
        p_array[2] = (float) pz_new;
    }

    public static void photon_energy_loss_corrections(int particle_Index, float[] p_array,
            HipoDataBank rec_Bank, HipoDataBank run_Bank) {

        double px = p_array[0];
        double py = p_array[1];
        double pz = p_array[2];

        double p = Math.sqrt(px * px + py * py + pz * pz);
        double theta = Math.acos(pz / p); // in radians
        double theta_deg = Math.toDegrees(theta);

        generic_tests genericTests = new generic_tests();
        boolean isForwardTagger = genericTests.forward_tagger_cut(particle_Index, rec_Bank);

        int runnum = run_Bank.getInt("run", 0);

        // Check if run is in allowed ranges
        boolean validRun = (runnum >= 4763 && runnum <= 5419) || 
                          (runnum >= 5423 && runnum <= 5666) || 
                          (runnum >= 6616 && runnum <= 6783);
        if (!validRun) {
            return;
        }

        double new_p = p;

        if (isForwardTagger) {
            // Apply spring19 fit correction for forward tagger
            new_p += 0.085643 - 0.0288063 * p + 0.00894691 * p * p - 0.000725449 * p * p * p;

            double x0 = 117.969;
            double a = -5.26113;
            double b = 0.0133069;
            double c = 0.000141808;
            double d = -3.30662e-06;

            double delta = a + b * (theta_deg - x0) + c * Math.pow(theta_deg - x0, 2) 
                         + d * Math.pow(theta_deg - x0, 3);
            new_p -= delta;
        } else {
            if (theta_deg < 20) {
                double x0 = 27.8139;
                double a = -0.028155;
                double b = -0.00496108;
                double c = -0.000369676;
                double d = -5.68592e-06;
                new_p -= (a + b * (theta_deg - x0) + c * Math.pow(theta_deg - x0, 2) 
                        + d * Math.pow(theta_deg - x0, 3));
            } else {
                double x0 = 17.4002;
                double a = -0.0338089;
                double b = 0.00296576;
                double c = -0.000212105;
                double d = -2.19736e-05;
                new_p -= (a + b * (theta_deg - x0) + c * Math.pow(theta_deg - x0, 2) 
                        + d * Math.pow(theta_deg - x0, 3));
            }
        }

        // Scale the momentum components to preserve direction
        if (p != 0) {
            double scale = (new_p) / p;
            p_array[0] = (float) (px * scale);
            p_array[1] = (float) (py * scale);
            p_array[2] = (float) (pz * scale);
        }
    }

}
