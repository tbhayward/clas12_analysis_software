/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package extended_kinematic_fitters;

import org.jlab.io.hipo.HipoDataBank;

/**
 *
 * @author tbhayward
 */
public class energy_loss_corrections {

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

    public double proton_energy_loss_corrections(int particle_Index, double px, double py, double pz,
            HipoDataBank rec_Bank, HipoDataBank run_Bank) {

        double p = Math.sqrt(px * px + py * py + pz * pz);
        double theta = theta_calculation(px, py, pz);
//        double phi = phi_calculation(px, py);

        generic_tests generic_tests = new generic_tests();
        boolean isForwardDetector = generic_tests.forward_detector_cut(particle_Index, rec_Bank);
        boolean isCentralDetector = generic_tests.forward_detector_cut(particle_Index, rec_Bank);

        int runnum = run_Bank.getInt("run", particle_Index);
        System.out.println(runnum);

        double dp = 0;
        double A_p = 0;
        double B_p = 0;
        double C_p = 0;
//        double A_theta, B_theta, C_theta;
//        double A_phi, B_phi, C_phi;
        if (runnum >= 4763 && runnum <= 5419) { // RGA Fa18 Inb

            if (isForwardDetector && !isCentralDetector) {
                A_p = 0.0099626 - 0.0002414 * theta - 0.0000020 * theta * theta;
                B_p = -0.01428267 + 0.00042833 * theta + 0.00001081 * theta * theta;
                C_p = 0.01197102 - 0.00055673 * theta + 0.00000785 * theta * theta;
            } else if (!isForwardDetector && isCentralDetector) {
                A_p = -0.2383991 + 0.0124992 * theta - 0.0001646 * theta * theta;
                B_p = 0.60123885 - 0.03128464 * theta + 0.00041314 * theta * theta;
                C_p = -0.44080146 + 0.02209857 * theta - 0.00028224 * theta * theta;
            }

        } else if (runnum >= 5423 && runnum <= 5666) { // RGA Fa18 Out

            if (isForwardDetector && !isCentralDetector) {
                A_p = 0.0135790 - 0.0005303 * theta;
                B_p = -0.02165929 + 0.00121123 * theta;
                C_p = 0.0;
            } else if (!isForwardDetector && isCentralDetector) {
                A_p = -0.1927861 + 0.0099546 * theta - 0.0001299 * theta * theta;
                B_p = 0.44307822 - 0.02309469 * theta + 0.00030784 * theta * theta;
                C_p = -0.32938000 + 0.01648659 * theta - 0.00021181 * theta * theta;
            }

        } else if (runnum >= 6616 && runnum <= 6783) { // RGA Sp19 Inb

            if (isForwardDetector && !isCentralDetector) {
                A_p = 0.0095205 - 0.0001914 * theta - 0.0000031 * theta * theta;
                B_p = -0.01365658 + 0.00036322 * theta + 0.00001217 * theta * theta;
                C_p = 0.01175256 - 0.00053407 * theta + 0.00000742 * theta * theta;
            } else if (!isForwardDetector && isCentralDetector) {
                A_p = -0.2716918 + 0.0142491 * theta - 0.0001862 * theta * theta;
                B_p = 0.65945101 - 0.03431360 * theta + 0.00045036 * theta * theta;
                C_p = -0.46602726 + 0.02335623 * theta - 0.00029720 * theta * theta;
            }

        } else if (runnum >= 16089 && runnum <= 16786) { // RGC Su22 Inb

            if (isForwardDetector && !isCentralDetector) {
                A_p = 0.0109317 - 0.0000194 * theta - 0.0000117 * theta * theta;
                B_p = -0.00910576 - 0.00035154 * theta + 0.00003905 * theta * theta;
                C_p = 0.01225782 - 0.00012805 * theta - 0.00000820 * theta * theta;
            } else if (!isForwardDetector && isCentralDetector) {
                A_p = -0.3951652 + 0.0202840 * theta - 0.0002660 * theta * theta;
                B_p = 0.93238668 - 0.04803619 * theta + 0.00063215 * theta * theta;
                C_p = -0.59146847 + 0.02997697 * theta - 0.00038773 * theta * theta;
            }

        }

        if (isForwardDetector && !isCentralDetector) {
            dp = A_p + B_p / p + C_p / (p * p);
        } else if (!isForwardDetector && isCentralDetector) {
            dp = A_p + B_p * p + C_p * p * p;
        }
        return (dp + p) / p; // fe
    }

}
