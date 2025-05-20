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

    public void stefan_piplus_energy_loss_corrections(int particle_Index, float[] p_array,
            HipoDataBank rec_Bank, HipoDataBank run_Bank, HipoDataBank track_Bank) {

        double dp = 0;
        double px = p_array[0];
        double py = p_array[1];
        double pz = p_array[2];

        double p = p_calculation(px, py, pz);
        double theta = theta_calculation(px, py, pz);

        generic_tests generic_tests = new generic_tests();
        boolean isForwardDetector = generic_tests.forward_detector_cut(particle_Index, rec_Bank);
        boolean isCentralDetector = generic_tests.central_detector_cut(particle_Index, rec_Bank);

        int runnum = run_Bank.getInt("run", 0);
        int runPeriod = -1;
        if (runnum >= 4763 && runnum <= 5419) {
            runPeriod = 1;
        } // RGA Fa18 Inb
        else if (runnum >= 5423 && runnum <= 5666) {
            runPeriod = 2;
        } // RGA Fa18 Out
        else if (runnum >= 6616 && runnum <= 6783) {
            runPeriod = 3;
        } // RGA Sp19 Inb

        if (runPeriod == 1 || runPeriod == 3) { // RGA Fa18 Inb and RGA Sp19 Inb
            if (isForwardDetector) {
                if (theta < 27 && p < 2.5) {
                    dp = 0.00342646 + (-0.00282934) * p + (0.00205983) * Math.pow(p, 2) + (-0.00043158) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta < 27 && p >= 2.5) {
                    dp = 0.00342646 + (-0.00282934) * 2.5 + (0.00205983) * Math.pow(2.5, 2) + (-0.00043158) * Math.pow(2.5, 3) + (0) * Math.pow(2.5, 4);
                } else if (theta > 27 && theta < 28 && p < 1.83) {
                    dp = 0.00328565 + (-0.00376042) * p + (0.00433886) * Math.pow(p, 2) + (-0.00141614) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 27 && theta < 28 && p >= 1.83) {
                    dp = 0.00328565 + (-0.00376042) * 1.83 + (0.00433886) * Math.pow(1.83, 2) + (-0.00141614) * Math.pow(1.83, 3) + (0) * Math.pow(1.83, 4);
                } else if (theta > 28 && theta < 29 && p < 2) {
                    dp = 0.00328579 + (-0.00281121) * p + (0.00342749) * Math.pow(p, 2) + (-0.000932614) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 28 && theta < 29 && p >= 2) {
                    dp = 0.00328579 + (-0.00281121) * 2 + (0.00342749) * Math.pow(2, 2) + (-0.000932614) * Math.pow(2, 3) + (0) * Math.pow(2, 4);
                } else if (theta > 29 && theta < 30 && p < 1.9) {
                    dp = 0.00167358 + (0.00441871) * p + (-0.000834667) * Math.pow(p, 2) + (-0.000137968) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 29 && theta < 30 && p >= 1.9) {
                    dp = 0.00167358 + (0.00441871) * 1.9 + (-0.000834667) * Math.pow(1.9, 2) + (-0.000137968) * Math.pow(1.9, 3) + (0) * Math.pow(1.9, 4);
                } else if (theta > 30 && theta < 31 && p < 1.9) {
                    dp = 0.00274159 + (0.00635686) * p + (-0.00380977) * Math.pow(p, 2) + (0.00071627) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 30 && theta < 31 && p >= 1.9) {
                    dp = 0.00274159 + (0.00635686) * 1.9 + (-0.00380977) * Math.pow(1.9, 2) + (0.00071627) * Math.pow(1.9, 3) + (0) * Math.pow(1.9, 4);
                } else if (theta > 31 && theta < 32 && p < 1.8) {
                    dp = 0.00450241 + (0.00248969) * p + (-0.00336795) * Math.pow(p, 2) + (0.00111193) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 31 && theta < 32 && p >= 1.8) {
                    dp = 0.00450241 + (0.00248969) * 1.8 + (-0.00336795) * Math.pow(1.8, 2) + (0.00111193) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                } else if (theta > 32 && theta < 33 && p < 1.8) {
                    dp = 0.00505593 + (-0.00246203) * p + (0.00172984) * Math.pow(p, 2) + (-0.000406701) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 32 && theta < 33 && p >= 1.8) {
                    dp = 0.00505593 + (-0.00246203) * 1.8 + (0.00172984) * Math.pow(1.8, 2) + (-0.000406701) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                } else if (theta > 33 && theta < 34 && p < 1.8) {
                    dp = 0.00273402 + (0.00440449) * p + (-0.00373488) * Math.pow(p, 2) + (0.000996612) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 33 && theta < 34 && p >= 1.8) {
                    dp = 0.00273402 + (0.00440449) * 1.8 + (-0.00373488) * Math.pow(1.8, 2) + (0.000996612) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                } else if (theta > 34 && theta < 35 && p < 1.8) {
                    dp = 0.00333542 + (0.00439874) * p + (-0.00397776) * Math.pow(p, 2) + (0.00105586) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 34 && theta < 35 && p >= 1.8) {
                    dp = 0.00333542 + (0.00439874) * 1.8 + (-0.00397776) * Math.pow(1.8, 2) + (0.00105586) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                } else if (theta > 35 && theta < 36 && p < 1.8) {
                    dp = 0.00354663 + (0.00565397) * p + (-0.00513503) * Math.pow(p, 2) + (0.00153346) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 35 && theta < 36 && p >= 1.8) {
                    dp = 0.00354663 + (0.00565397) * 1.8 + (-0.00513503) * Math.pow(1.8, 2) + (0.00153346) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                } else if (theta > 36 && theta < 37 && p < 1.8) {
                    dp = 0.00333909 + (0.00842367) * p + (-0.0077321) * Math.pow(p, 2) + (0.0022489) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 36 && theta < 37 && p >= 1.8) {
                    dp = 0.00333909 + (0.00842367) * 1.8 + (-0.0077321) * Math.pow(1.8, 2) + (0.0022489) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                } else if (theta > 37 && theta < 38 && p < 1.4) {
                    dp = 0.00358828 + (0.0112108) * p + (-0.0133854) * Math.pow(p, 2) + (0.00486924) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 37 && theta < 38 && p >= 1.4) {
                    dp = 0.00358828 + (0.0112108) * 1.4 + (-0.0133854) * Math.pow(1.4, 2) + (0.00486924) * Math.pow(1.4, 3) + (0) * Math.pow(1.4, 4);
                } else if (theta > 38 && theta < 39 && p < 1.3) {
                    dp = 0.00354343 + (0.0117121) * p + (-0.0129649) * Math.pow(p, 2) + (0.00455602) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 38 && theta < 39 && p >= 1.3) {
                    dp = 0.00354343 + (0.0117121) * 1.3 + (-0.0129649) * Math.pow(1.3, 2) + (0.00455602) * Math.pow(1.3, 3) + (0) * Math.pow(1.3, 4);
                } else if (theta > 39 && theta < 40 && p < 0.9) {
                    dp = -0.00194951 + (0.0409713) * p + (-0.0595861) * Math.pow(p, 2) + (0.0281588) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 39 && theta < 40 && p >= 0.9) {
                    dp = -0.00194951 + (0.0409713) * 0.9 + (-0.0595861) * Math.pow(0.9, 2) + (0.0281588) * Math.pow(0.9, 3) + (0) * Math.pow(0.9, 4);
                } else if (theta > 40 && theta < 41 && p < 0.75) {
                    dp = -0.0099217 + (0.0808096) * p + (-0.119836) * Math.pow(p, 2) + (0.0559553) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 40 && theta < 41 && p >= 0.75) {
                    dp = -0.0099217 + (0.0808096) * 0.75 + (-0.119836) * Math.pow(0.75, 2) + (0.0559553) * Math.pow(0.75, 3) + (0) * Math.pow(0.75, 4);
                } else if (theta > 41 && theta < 42 && p < 0.65) {
                    dp = 0.00854898 + (0.00025037) * p + (-0.0113992) * Math.pow(p, 2) + (0.0145178) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 41 && theta < 42 && p >= 0.65) {
                    dp = 0.00854898 + (0.00025037) * 0.65 + (-0.0113992) * Math.pow(0.65, 2) + (0.0145178) * Math.pow(0.65, 3) + (0) * Math.pow(0.65, 4);
                } else if (theta > 42 && p < 0.65) {
                    dp = 0.00564818 + (0.00706606) * p + (0.0042602) * Math.pow(p, 2) + (-0.01141) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 42 && p >= 0.65) {
                    dp = 0.00564818 + (0.00706606) * 0.65 + (0.0042602) * Math.pow(0.65, 2) + (-0.01141) * Math.pow(0.65, 3) + (0) * Math.pow(0.65, 4);
                }
            } else if (isCentralDetector) {
                if (theta < 39 && p < 0.7) {
                    dp = -0.045 + (-0.102652) + (0.455589) * p + (-0.671635) * Math.pow(p, 2) + (0.303814) * Math.pow(p, 3);
                } else if (theta < 39 && p >= 0.7) {
                    dp = -0.045 + (-0.102652) + (0.455589) * 0.7 + (-0.671635) * Math.pow(0.7, 2) + (0.303814) * Math.pow(0.7, 3);
                } else if (theta > 39 && theta < 40 && p < 1.4) {
                    dp = 0.0684552 + (-0.766492) * p + (1.73092) * Math.pow(p, 2) + (-1.46215) * Math.pow(p, 3) + (0.420127) * Math.pow(p, 4);
                } else if (theta > 39 && theta < 40 && p >= 1.4) {
                    dp = 0.0684552 + (-0.766492) * 1.4 + (1.73092) * Math.pow(1.4, 2) + (-1.46215) * Math.pow(1.4, 3) + (0.420127) * Math.pow(1.4, 4);
                } else if (theta > 40 && theta < 41 && p < 1.45) {
                    dp = 0.751549 + (-7.4593) * p + (26.8037) * Math.pow(p, 2) + (-47.1576) * Math.pow(p, 3) + (43.8527) * Math.pow(p, 4) + (-20.7039) * Math.pow(p, 5) + (3.90931) * Math.pow(p, 6);
                } else if (theta > 40 && theta < 41 && p >= 1.45) {
                    dp = 0.751549 + (-7.4593) * 1.45 + (26.8037) * Math.pow(1.45, 2) + (-47.1576) * Math.pow(1.45, 3) + (43.8527) * Math.pow(1.45, 4) + (-20.7039) * Math.pow(1.45, 5) + (3.90931) * Math.pow(1.45, 6);
                } else if (theta > 41 && theta < 42 && p < 1.2) {
                    dp = -1.35043 + (10.0788) * p + (-30.4829) * Math.pow(p, 2) + (47.7792) * Math.pow(p, 3) + (-40.996) * Math.pow(p, 4) + (18.2662) * Math.pow(p, 5) + (-3.30449) * Math.pow(p, 6);
                } else if (theta > 41 && theta < 42 && p >= 1.2) {
                    dp = -1.35043 + (10.0788) * 1.2 + (-30.4829) * Math.pow(1.2, 2) + (47.7792) * Math.pow(1.2, 3) + (-40.996) * Math.pow(1.2, 4) + (18.2662) * Math.pow(1.2, 5) + (-3.30449) * Math.pow(1.2, 6);
                } else if (theta > 42 && theta < 43 && p < 1.3) {
                    dp = -0.0231195 + (0.0744589) * p + (-0.0807029) * Math.pow(p, 2) + (0.0264266) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 42 && theta < 43 && p >= 1.3) {
                    dp = -0.0231195 + (0.0744589) * 1.3 + (-0.0807029) * Math.pow(1.3, 2) + (0.0264266) * Math.pow(1.3, 3) + (0) * Math.pow(1.3, 4);
                } else if (theta > 43 && theta < 44 && p < 1.1) {
                    dp = -0.00979928 + (0.0351043) * p + (-0.0365865) * Math.pow(p, 2) + (0.00977218) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 43 && theta < 44 && p >= 1.1) {
                    dp = -0.00979928 + (0.0351043) * 1.1 + (-0.0365865) * Math.pow(1.1, 2) + (0.00977218) * Math.pow(1.1, 3) + (0) * Math.pow(1.1, 4);
                } else if (theta > 44 && theta < 45 && p < 1.1) {
                    dp = 0.00108491 + (-0.00924885) * p + (0.0216431) * Math.pow(p, 2) + (-0.0137762) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 44 && theta < 45 && p >= 1.1) {
                    dp = 0.00108491 + (-0.00924885) * 1.1 + (0.0216431) * Math.pow(1.1, 2) + (-0.0137762) * Math.pow(1.1, 3) + (0) * Math.pow(1.1, 4);
                } else if (theta > 45 && theta < 55 && p < 1.3) {
                    dp = 0.0092263 + (-0.0676178) * p + (0.168778) * Math.pow(p, 2) + (-0.167463) * Math.pow(p, 3) + (0.05661) * Math.pow(p, 4);
                } else if (theta > 45 && theta < 55 && p >= 1.3) {
                    dp = 0.0092263 + (-0.0676178) * 1.3 + (0.168778) * Math.pow(1.3, 2) + (-0.167463) * Math.pow(1.3, 3) + (0.05661) * Math.pow(1.3, 4);
                } else if (theta > 55 && theta < 65 && p < 1.05) {
                    dp = 0.00805642 + (-0.0670962) * p + (0.188536) * Math.pow(p, 2) + (-0.20571) * Math.pow(p, 3) + (0.0765) * Math.pow(p, 4);
                } else if (theta > 55 && theta < 65 && p >= 1.05) {
                    dp = 0.00805642 + (-0.0670962) * 1.05 + (0.188536) * Math.pow(1.05, 2) + (-0.20571) * Math.pow(1.05, 3) + (0.0765) * Math.pow(1.05, 4);
                } else if (theta > 65 && theta < 75 && p < 0.75) {
                    dp = 0.00312202 + (-0.0269717) * p + (0.0715236) * Math.pow(p, 2) + (-0.0545622) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 65 && theta < 75 && p >= 0.75) {
                    dp = 0.00312202 + (-0.0269717) * 0.75 + (0.0715236) * Math.pow(0.75, 2) + (-0.0545622) * Math.pow(0.75, 3) + (0) * Math.pow(0.75, 4);
                } else if (theta > 75 && theta < 85 && p < 0.65) {
                    dp = 0.00424971 + (-0.0367683) * p + (0.10417) * Math.pow(p, 2) + (-0.0899651) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 75 && theta < 85 && p >= 0.65) {
                    dp = 0.00424971 + (-0.0367683) * 0.65 + (0.10417) * Math.pow(0.65, 2) + (-0.0899651) * Math.pow(0.65, 3) + (0) * Math.pow(0.65, 4);
                } else if (theta > 85 && theta < 95 && p < 0.5) {
                    dp = 0.00654123 + (-0.0517915) * p + (0.147888) * Math.pow(p, 2) + (-0.14253) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 85 && theta < 95 && p >= 0.5) {
                    dp = 0.00654123 + (-0.0517915) * 0.5 + (0.147888) * Math.pow(0.5, 2) + (-0.14253) * Math.pow(0.5, 3) + (0) * Math.pow(0.5, 4);
                } else if (theta > 95 && theta < 105 && p < 0.45) {
                    dp = -0.00111721 + (0.00478119) * p + (0.0158753) * Math.pow(p, 2) + (-0.052902) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 95 && theta < 105 && p >= 0.45) {
                    dp = -0.00111721 + (0.00478119) * 0.45 + (0.0158753) * Math.pow(0.45, 2) + (-0.052902) * Math.pow(0.45, 3) + (0) * Math.pow(0.45, 4);
                } else if (theta > 105 && theta < 115 && p < 0.35) {
                    dp = -0.00239839 + (0.00790738) * p + (0.0311713) * Math.pow(p, 2) + (-0.104157) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 105 && theta < 115 && p >= 0.35) {
                    dp = -0.00239839 + (0.00790738) * 0.35 + (0.0311713) * Math.pow(0.35, 2) + (-0.104157) * Math.pow(0.35, 3) + (0) * Math.pow(0.35, 4);
                } else if (theta > 115 && theta < 125 && p < 0.35) {
                    dp = -0.00778793 + (0.0256774) * p + (0.0932503) * Math.pow(p, 2) + (-0.32771) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 115 && theta < 125 && p >= 0.35) {
                    dp = -0.00778793 + (0.0256774) * 0.35 + (0.0932503) * Math.pow(0.35, 2) + (-0.32771) * Math.pow(0.35, 3) + (0) * Math.pow(0.35, 4);
                } else if (theta > 125 && theta < 135 && p < 0.35) {
                    dp = -0.00292778 + (-0.00536697) * p + (-0.00414351) * Math.pow(p, 2) + (0.0196431) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 125 && theta < 135 && p >= 0.35) {
                    dp = -0.00292778 + (-0.00536697) * 0.35 + (-0.00414351) * Math.pow(0.35, 2) + (0.0196431) * Math.pow(0.35, 3) + (0) * Math.pow(0.35, 4);
                };
            }
        } else if (runPeriod == 2) { // RGA Fa18 Out
            if (isForwardDetector) {
                if (theta < 27 && p < 2.3) {
                    dp = 0.00389945 + (-0.004062) * p + (0.00321842) * Math.pow(p, 2) + (-0.000698299) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta < 27 && p >= 2.3) {
                    dp = 0.00389945 + (-0.004062) * 2.3 + (0.00321842) * Math.pow(2.3, 2) + (-0.000698299) * Math.pow(2.3, 3) + (0) * Math.pow(2.3, 4);
                } else if (theta > 27 && theta < 28 && p < 1.7) {
                    dp = 0.00727132 + (-0.0117989) * p + (0.00962999) * Math.pow(p, 2) + (-0.00267005) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 27 && theta < 28 && p >= 1.7) {
                    dp = 0.00727132 + (-0.0117989) * 1.7 + (0.00962999) * Math.pow(1.7, 2) + (-0.00267005) * Math.pow(1.7, 3) + (0) * Math.pow(1.7, 4);
                } else if (theta > 28 && theta < 29 && p < 2) {
                    dp = 0.00844551 + (-0.0128097) * p + (0.00945956) * Math.pow(p, 2) + (-0.00237992) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 28 && theta < 29 && p >= 2) {
                    dp = 0.00844551 + (-0.0128097) * 2 + (0.00945956) * Math.pow(2, 2) + (-0.00237992) * Math.pow(2, 3) + (0) * Math.pow(2, 4);
                } else if (theta > 29 && theta < 30 && p < 1.9) {
                    dp = 0.00959007 + (-0.0139218) * p + (0.0122966) * Math.pow(p, 2) + (-0.0034012) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 29 && theta < 30 && p >= 1.9) {
                    dp = 0.00959007 + (-0.0139218) * 1.9 + (0.0122966) * Math.pow(1.9, 2) + (-0.0034012) * Math.pow(1.9, 3) + (0) * Math.pow(1.9, 4);
                } else if (theta > 30 && theta < 31 && p < 1.9) {
                    dp = 0.00542816 + (-5.10739e-05) * p + (0.000572038) * Math.pow(p, 2) + (-0.000488883) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 30 && theta < 31 && p >= 1.9) {
                    dp = 0.00542816 + (-5.10739e-05) * 1.9 + (0.000572038) * Math.pow(1.9, 2) + (-0.000488883) * Math.pow(1.9, 3) + (0) * Math.pow(1.9, 4);
                } else if (theta > 31 && theta < 32 && p < 1.8) {
                    dp = 0.0060391 + (-0.000516936) * p + (-0.00286595) * Math.pow(p, 2) + (0.00136604) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 31 && theta < 32 && p >= 1.8) {
                    dp = 0.0060391 + (-0.000516936) * 1.8 + (-0.00286595) * Math.pow(1.8, 2) + (0.00136604) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                } else if (theta > 32 && theta < 33 && p < 1.6) {
                    dp = 0.0140305 + (-0.0285832) * p + (0.0248799) * Math.pow(p, 2) + (-0.00701311) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 32 && theta < 33 && p >= 1.6) {
                    dp = 0.0140305 + (-0.0285832) * 1.6 + (0.0248799) * Math.pow(1.6, 2) + (-0.00701311) * Math.pow(1.6, 3) + (0) * Math.pow(1.6, 4);
                } else if (theta > 33 && theta < 34 && p < 1.5) {
                    dp = 0.010815 + (-0.0194244) * p + (0.0174474) * Math.pow(p, 2) + (-0.0049764) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 33 && theta < 34 && p >= 1.5) {
                    dp = 0.010815 + (-0.0194244) * 1.5 + (0.0174474) * Math.pow(1.5, 2) + (-0.0049764) * Math.pow(1.5, 3) + (0) * Math.pow(1.5, 4);
                } else if (theta > 34 && theta < 35 && p < 1.6) {
                    dp = 0.0105522 + (-0.0176248) * p + (0.0161142) * Math.pow(p, 2) + (-0.00472288) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 34 && theta < 35 && p >= 1.6) {
                    dp = 0.0105522 + (-0.0176248) * 1.6 + (0.0161142) * Math.pow(1.6, 2) + (-0.00472288) * Math.pow(1.6, 3) + (0) * Math.pow(1.6, 4);
                } else if (theta > 35 && theta < 36 && p < 1.5) {
                    dp = 0.0103938 + (-0.0164003) * p + (0.0164045) * Math.pow(p, 2) + (-0.00517012) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 35 && theta < 36 && p >= 1.5) {
                    dp = 0.0103938 + (-0.0164003) * 1.5 + (0.0164045) * Math.pow(1.5, 2) + (-0.00517012) * Math.pow(1.5, 3) + (0) * Math.pow(1.5, 4);
                } else if (theta > 36 && theta < 37 && p < 1.8) {
                    dp = 0.0441471 + (-0.183937) * p + (0.338784) * Math.pow(p, 2) + (-0.298985) * Math.pow(p, 3) + (0.126905) * Math.pow(p, 4) + (-0.0208286) * Math.pow(p, 5);
                } else if (theta > 36 && theta < 37 && p >= 1.8) {
                    dp = 0.0441471 + (-0.183937) * 1.8 + (0.338784) * Math.pow(1.8, 2) + (-0.298985) * Math.pow(1.8, 3) + (0.126905) * Math.pow(1.8, 4) + (-0.0208286) * Math.pow(1.8, 5);
                } else if (theta > 37 && theta < 38 && p < 1.7) {
                    dp = 0.0726119 + (-0.345004) * p + (0.697789) * Math.pow(p, 2) + (-0.685948) * Math.pow(p, 3) + (0.327195) * Math.pow(p, 4) + (-0.0605621) * Math.pow(p, 5);
                } else if (theta > 37 && theta < 38 && p >= 1.7) {
                    dp = 0.0726119 + (-0.345004) * 1.7 + (0.697789) * Math.pow(1.7, 2) + (-0.685948) * Math.pow(1.7, 3) + (0.327195) * Math.pow(1.7, 4) + (-0.0605621) * Math.pow(1.7, 5);
                } else if (theta > 38 && theta < 39 && p < 1.6) {
                    dp = 0.0247648 + (-0.0797376) * p + (0.126535) * Math.pow(p, 2) + (-0.086545) * Math.pow(p, 3) + (0.0219304) * Math.pow(p, 4);
                } else if (theta > 38 && theta < 39 && p >= 1.6) {
                    dp = 0.0247648 + (-0.0797376) * 1.6 + (0.126535) * Math.pow(1.6, 2) + (-0.086545) * Math.pow(1.6, 3) + (0.0219304) * Math.pow(1.6, 4);
                } else if (theta > 39 && theta < 40 && p < 1.2) {
                    dp = 0.0208867 + (-0.0492068) * p + (0.0543187) * Math.pow(p, 2) + (-0.0183393) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 39 && theta < 40 && p >= 1.2) {
                    dp = 0.0208867 + (-0.0492068) * 1.2 + (0.0543187) * Math.pow(1.2, 2) + (-0.0183393) * Math.pow(1.2, 3) + (0) * Math.pow(1.2, 4);
                } else if (theta > 40 && theta < 41 && p < 1.0) {
                    dp = 0.0148655 + (-0.0203483) * p + (0.00835867) * Math.pow(p, 2) + (0.00697134) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 40 && theta < 41 && p >= 1.0) {
                    dp = 0.0148655 + (-0.0203483) * 1.0 + (0.00835867) * Math.pow(1.0, 2) + (0.00697134) * Math.pow(1.0, 3) + (0) * Math.pow(1.0, 4);
                } else if (theta > 41 && theta < 42 && p < 0.7) {
                    dp = 0.0223585 + (-0.0365262) * p + (-0.0150027) * Math.pow(p, 2) + (0.0854164) * Math.pow(p, 3) + (-0.0462718) * Math.pow(p, 4);
                } else if (theta > 41 && theta < 42 && p >= 0.7) {
                    dp = 0.007617;
                } else if (theta > 42 && p < 0.75) {
                    dp = 0.0152373 + (-0.0106377) * p + (-0.0257573) * Math.pow(p, 2) + (0.0344851) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                } else if (theta > 42 && p >= 0.75) {
                    dp = 0.0152373 + (-0.0106377) * 0.75 + (-0.0257573) * Math.pow(0.75, 2) + (0.0344851) * Math.pow(0.75, 3) + (0) * Math.pow(0.75, 4);
                }
            } else if (isCentralDetector) {
                if (theta < 39 && p < 0.8) {
                    dp = -0.05 + (-0.0758897) + (0.362231) * p + (-0.542404) * Math.pow(p, 2) + (0.241344) * Math.pow(p, 3);
                }
                if (theta < 39 && p >= 0.8) {
                    dp = -0.05 + (-0.0758897) + (0.362231) * 0.8 + (-0.542404) * Math.pow(0.8, 2) + (0.241344) * Math.pow(0.8, 3);
                }
                if (theta > 39 && theta < 40 && p < 1.35) {
                    dp = 0.0355259 + (-0.589712) * p + (1.4206) * Math.pow(p, 2) + (-1.24179) * Math.pow(p, 3) + (0.365524) * Math.pow(p, 4);
                }
                if (theta > 39 && theta < 40 && p >= 1.35) {
                    dp = 0.0355259 + (-0.589712) * 1.35 + (1.4206) * Math.pow(1.35, 2) + (-1.24179) * Math.pow(1.35, 3) + (0.365524) * Math.pow(1.35, 4);
                }
                if (theta > 40 && theta < 41 && p < 1.4) {
                    dp = -0.252336 + (1.02032) * p + (-1.51461) * Math.pow(p, 2) + (0.967772) * Math.pow(p, 3) + (-0.226028) * Math.pow(p, 4);
                }
                if (theta > 40 && theta < 41 && p >= 1.4) {
                    dp = -0.252336 + (1.02032) * 1.4 + (-1.51461) * Math.pow(1.4, 2) + (0.967772) * Math.pow(1.4, 3) + (-0.226028) * Math.pow(1.4, 4);
                }
                if (theta > 41 && theta < 42 && p < 1.2) {
                    dp = -0.710129 + (4.49613) * p + (-11.01) * Math.pow(p, 2) + (12.9945) * Math.pow(p, 3) + (-7.41641) * Math.pow(p, 4) + (1.63923) * Math.pow(p, 5);
                }
                if (theta > 41 && theta < 42 && p >= 1.2) {
                    dp = -0.710129 + (4.49613) * 1.2 + (-11.01) * Math.pow(1.2, 2) + (12.9945) * Math.pow(1.2, 3) + (-7.41641) * Math.pow(1.2, 4) + (1.63923) * Math.pow(1.2, 5);
                }
                if (theta > 42 && theta < 43 && p < 1.2) {
                    dp = -0.0254912 + (0.0851432) * p + (-0.0968583) * Math.pow(p, 2) + (0.0350334) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 42 && theta < 43 && p >= 1.2) {
                    dp = -0.0254912 + (0.0851432) * 1.2 + (-0.0968583) * Math.pow(1.2, 2) + (0.0350334) * Math.pow(1.2, 3) + (0) * Math.pow(1.2, 4);
                }
                if (theta > 43 && theta < 44 && p < 1.4) {
                    dp = -0.0115965 + (0.0438726) * p + (-0.0500474) * Math.pow(p, 2) + (0.0163627) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 43 && theta < 44 && p >= 1.4) {
                    dp = -0.0115965 + (0.0438726) * 1.4 + (-0.0500474) * Math.pow(1.4, 2) + (0.0163627) * Math.pow(1.4, 3) + (0) * Math.pow(1.4, 4);
                }
                if (theta > 44 && theta < 45 && p < 1) {
                    dp = 0.00273414 + (-0.01851) * p + (0.0377032) * Math.pow(p, 2) + (-0.0226696) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 44 && theta < 45 && p >= 1) {
                    dp = 0.00273414 + (-0.01851) * 1 + (0.0377032) * Math.pow(1, 2) + (-0.0226696) * Math.pow(1, 3) + (0) * Math.pow(1, 4);
                }
                if (theta > 45 && theta < 55 && p < 1.4) {
                    dp = 0.0271952 + (-0.25981) * p + (0.960051) * Math.pow(p, 2) + (-1.76651) * Math.pow(p, 3) + (1.72872) * Math.pow(p, 4) + (-0.856946) * Math.pow(p, 5) + (0.167564) * Math.pow(p, 6);
                }
                if (theta > 45 && theta < 55 && p >= 1.4) {
                    dp = 0.0271952 + (-0.25981) * 1.4 + (0.960051) * Math.pow(1.4, 2) + (-1.76651) * Math.pow(1.4, 3) + (1.72872) * Math.pow(1.4, 4) + (-0.856946) * Math.pow(1.4, 5) + (0.167564) * Math.pow(1.4, 6);
                }
                if (theta > 55 && theta < 65 && p < 1.2) {
                    dp = 0.00734975 + (-0.0598841) * p + (0.161495) * Math.pow(p, 2) + (-0.1629) * Math.pow(p, 3) + (0.0530098) * Math.pow(p, 4);
                }
                if (theta > 55 && theta < 65 && p >= 1.2) {
                    dp = 0.00734975 + (-0.0598841) * 1.2 + (0.161495) * Math.pow(1.2, 2) + (-0.1629) * Math.pow(1.2, 3) + (0.0530098) * Math.pow(1.2, 4);
                }
                if (theta > 65 && theta < 75 && p < 0.95) {
                    dp = 0.00321351 + (-0.0289322) * p + (0.0786484) * Math.pow(p, 2) + (-0.0607041) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 65 && theta < 75 && p >= 0.95) {
                    dp = 0.00321351 + (-0.0289322) * 0.95 + (0.0786484) * Math.pow(0.95, 2) + (-0.0607041) * Math.pow(0.95, 3) + (0) * Math.pow(0.95, 4);
                }
                if (theta > 75 && theta < 85 && p < 0.7) {
                    dp = 0.00644253 + (-0.0543896) * p + (0.148933) * Math.pow(p, 2) + (-0.1256) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 75 && theta < 85 && p >= 0.7) {
                    dp = 0.00644253 + (-0.0543896) * 0.7 + (0.148933) * Math.pow(0.7, 2) + (-0.1256) * Math.pow(0.7, 3) + (0) * Math.pow(0.7, 4);
                }
                if (theta > 85 && theta < 95 && p < 0.65) {
                    dp = 0.00671152 + (-0.0537269) * p + (0.154509) * Math.pow(p, 2) + (-0.147667) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 85 && theta < 95 && p >= 0.65) {
                    dp = 0.00671152 + (-0.0537269) * 0.65 + (0.154509) * Math.pow(0.65, 2) + (-0.147667) * Math.pow(0.65, 3) + (0) * Math.pow(0.65, 4);
                }
                if (theta > 95 && theta < 105 && p < 0.45) {
                    dp = -0.000709077 + (0.00331818) * p + (0.0109241) * Math.pow(p, 2) + (-0.0351682) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 95 && theta < 105 && p >= 0.45) {
                    dp = -0.000709077 + (0.00331818) * 0.45 + (0.0109241) * Math.pow(0.45, 2) + (-0.0351682) * Math.pow(0.45, 3) + (0) * Math.pow(0.45, 4);
                }
                if (theta > 105 && theta < 115 && p < 0.45) {
                    dp = -0.00260164 + (0.00846919) * p + (0.0315497) * Math.pow(p, 2) + (-0.105756) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 105 && theta < 115 && p >= 0.45) {
                    dp = -0.00260164 + (0.00846919) * 0.45 + (0.0315497) * Math.pow(0.45, 2) + (-0.105756) * Math.pow(0.45, 3) + (0) * Math.pow(0.45, 4);
                }
                if (theta > 115 && theta < 125 && p < 0.45) {
                    dp = -0.00544336 + (0.018256) * p + (0.0664618) * Math.pow(p, 2) + (-0.240312) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 115 && theta < 125 && p >= 0.45) {
                    dp = -0.00544336 + (0.018256) * 0.45 + (0.0664618) * Math.pow(0.45, 2) + (-0.240312) * Math.pow(0.45, 3) + (0) * Math.pow(0.45, 4);
                }
                if (theta > 125 && theta < 135 && p < 0.35) {
                    dp = -0.00281073 + (-0.00495863) * p + (-0.00362356) * Math.pow(p, 2) + (0.0178764) * Math.pow(p, 3) + (0) * Math.pow(p, 4);
                }
                if (theta > 125 && theta < 135 && p >= 0.35) {
                    dp = -0.00281073 + (-0.00495863) * 0.35 + (-0.00362356) * Math.pow(0.35, 2) + (0.0178764) * Math.pow(0.35, 3) + (0) * Math.pow(0.35, 4);
                }
            }
        }

        if (p != 0) {
            double scale = (p + dp) / p;
            p_array[0] = (float) (px * scale);
            p_array[1] = (float) (py * scale);
            p_array[2] = (float) (pz * scale);
        }
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

    public void mariana_proton_energy_loss_corrections(int particle_Index, float[] p_array,
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

        if (isForwardDetector) {
            if (runnum >= 4763 && runnum <= 6783) {
                double[][] ag = {
                    {5.00, 7.00, 9.00, 11.00, 13.00, 15.00, 17.00, 19.00, 21.00, 23.00, 25.00, 27.00, 29.00, 31.00, 33.00},
                    {0.1039E-02, -0.6952E-02, -0.9509E-02, -0.9879E-02, -0.1279E-01, -0.1157E-01, -0.1018E-01, -0.9222E-02, -0.1355E-01, -0.1207E-01, -0.9474E-02, -0.2216E-01, -0.2105E-01, -0.2118E-01, -0.2360E-01},
                    {-0.6922E-03, 0.9763E-03, 0.1482E-02, 0.1530E-02, 0.2187E-02, 0.1953E-02, 0.1688E-02, 0.1668E-02, 0.2849E-02, 0.2495E-02, 0.1508E-02, 0.4215E-02, 0.3911E-02, 0.3948E-02, 0.4634E-02},
                    {0.9806E-03, 0.1157E-01, 0.1485E-01, 0.1588E-01, 0.1945E-01, 0.1736E-01, 0.1551E-01, 0.1383E-01, 0.1926E-01, 0.1720E-01, 0.1464E-01, 0.3250E-01, 0.3231E-01, 0.3296E-01, 0.3608E-01},
                    {-0.8024E-02, -0.1035E-01, -0.1240E-01, -0.1361E-01, -0.1518E-01, -0.1432E-01, -0.1341E-01, -0.1255E-01, -0.1462E-01, -0.1388E-01, -0.1574E-01, -0.2646E-01, -0.2820E-01, -0.3000E-01, -0.3259E-01}
                };
                if (theta * Math.PI / 180 <= 5) {
                    p = p - p * (ag[1][0] + ag[2][0] * p + ag[3][0] / p + ag[4][0] / (p * p));
                } else if (theta * Math.PI / 180 >= 33) {
                    p = p - p * (ag[1][14] + ag[2][14] * p + ag[3][14] / p + ag[4][14]  / (p * p));
                } else {
                    for (int i = 0; i < 10; i++) {
                        p = p - p * (ag[1][i] + ag[2][i] * p + ag[3][i] / p + ag[4][i] / (p * p));
                        break;
                    }
                }
            }
        } else if (isCentralDetector) {
            // FD corrections only
        }

        // Update the px, py, pz values
        p_array[0] = (float) x_calculation(p, theta, phi);
        p_array[1] = (float) y_calculation(p, theta, phi);
        p_array[2] = (float) z_calculation(p, theta);
    }

    public void krishna_energy_loss_corrections(int particle_Index, float[] p_array,
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
        int pid = rec_Bank.getInt("pid", particle_Index);

        int runnum = run_Bank.getInt("run", 0);

        double dp = 0;

        if (pid == 2212) { // proton, from Krishna
            if (isForwardDetector) {
                if (theta < 27 * Math.PI / 180) {
                    if (p < 2.4) {
                        p = p + (0.001046) * Math.pow(p, 4) + (-0.010446) * Math.pow(p, 3) + (0.036945) * Math.pow(p, 2)
                                + (-0.055368) * p + 0.034539;
                    } else if (p >= 2.4) {
                        p = p + 0.004741;
                    }
                } else if (theta >= 27 * Math.PI / 180) {
                    if (p < 2.4) {
                        p = p + (0.005519) * Math.pow(p, 4) + (-0.046289) * Math.pow(p, 3)
                                + (0.137504) * Math.pow(p, 2) + (-0.177027) * p + 0.094555;
                    } else if (p >= 2.4) {
                        p = p + 0.004899;
                    }
                }
            } else if (isCentralDetector) {
                // no central detector proton correction
            }
        }

        if (pid == -211) { // pi minus, from Krishna
            if (isForwardDetector) {
                if (theta < 27 * Math.PI / 180) {
                    if (p < 2.4) {
                        p = p + (0.00046571 + 0.00322164);
                    }
                } else if (theta >= 27 * Math.PI / 180) {
                    if (p < 1.7) {
                        p = p + (-0.0024313) * Math.pow(p, 3) + (0.0094416) * Math.pow(p, 2)
                                + (-0.01257967) * Math.pow(p, 1) + 0.0122432;
                    } else if (p >= 1.7) {
                        p = p + 0.006199071;
                    }
                }
            } else if (isCentralDetector) {
                // no central detector pi minus correction
            }
        }

        double pion_p = p;
        double dp_pion = 0;
        double pip_theta = theta * Math.PI / 180;
        if (pid == 211) { // pi plus, from Stefan
            if (isForwardDetector) {
                if ((runnum >= 4763 && runnum <= 5419) || (runnum >= 6616 && runnum <= 6783)) { // RGA inbending
                    if (pip_theta < 27) {
                        dp_pion = 0.00342646 + (-0.00282934) * pion_p + (0.00205983) * Math.pow(pion_p, 2) + (-0.00043158) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta < 27 && pion_p >= 2.5) {
                        dp_pion = 0.00342646 + (-0.00282934) * 2.5 + (0.00205983) * Math.pow(2.5, 2) + (-0.00043158) * Math.pow(2.5, 3) + (0) * Math.pow(2.5, 4);
                    }
                    if (pip_theta > 27 && pip_theta < 28) {
                        dp_pion = 0.00328565 + (-0.00376042) * pion_p + (0.00433886) * Math.pow(pion_p, 2) + (-0.00141614) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 27 && pip_theta < 28 && pion_p >= 1.83) {
                        dp_pion = 0.00328565 + (-0.00376042) * 1.83 + (0.00433886) * Math.pow(1.83, 2) + (-0.00141614) * Math.pow(1.83, 3) + (0) * Math.pow(1.83, 4);
                    }
                    if (pip_theta > 28 && pip_theta < 29) {
                        dp_pion = 0.00328579 + (-0.00281121) * pion_p + (0.00342749) * Math.pow(pion_p, 2) + (-0.000932614) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 28 && pip_theta < 29 && pion_p >= 2) {
                        dp_pion = 0.00328579 + (-0.00281121) * 2 + (0.00342749) * Math.pow(2, 2) + (-0.000932614) * Math.pow(2, 3) + (0) * Math.pow(2, 4);
                    }
                    if (pip_theta > 29 && pip_theta < 30) {
                        dp_pion = 0.00167358 + (0.00441871) * pion_p + (-0.000834667) * Math.pow(pion_p, 2) + (-0.000137968) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 29 && pip_theta < 30 && pion_p >= 1.9) {
                        dp_pion = 0.00167358 + (0.00441871) * 1.9 + (-0.000834667) * Math.pow(1.9, 2) + (-0.000137968) * Math.pow(1.9, 3) + (0) * Math.pow(1.9, 4);
                    }
                    if (pip_theta > 30 && pip_theta < 31) {
                        dp_pion = 0.00274159 + (0.00635686) * pion_p + (-0.00380977) * Math.pow(pion_p, 2) + (0.00071627) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 30 && pip_theta < 31 && pion_p >= 1.9) {
                        dp_pion = 0.00274159 + (0.00635686) * 1.9 + (-0.00380977) * Math.pow(1.9, 2) + (0.00071627) * Math.pow(1.9, 3) + (0) * Math.pow(1.9, 4);
                    }
                    if (pip_theta > 31 && pip_theta < 32) {
                        dp_pion = 0.00450241 + (0.00248969) * pion_p + (-0.00336795) * Math.pow(pion_p, 2) + (0.00111193) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 31 && pip_theta < 32 && pion_p >= 1.8) {
                        dp_pion = 0.00450241 + (0.00248969) * 1.8 + (-0.00336795) * Math.pow(1.8, 2) + (0.00111193) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                    }
                    if (pip_theta > 32 && pip_theta < 33) {
                        dp_pion = 0.00505593 + (-0.00246203) * pion_p + (0.00172984) * Math.pow(pion_p, 2) + (-0.000406701) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 32 && pip_theta < 33 && pion_p >= 1.8) {
                        dp_pion = 0.00505593 + (-0.00246203) * 1.8 + (0.00172984) * Math.pow(1.8, 2) + (-0.000406701) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                    }
                    if (pip_theta > 33 && pip_theta < 34) {
                        dp_pion = 0.00273402 + (0.00440449) * pion_p + (-0.00373488) * Math.pow(pion_p, 2) + (0.000996612) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 33 && pip_theta < 34 && pion_p >= 1.8) {
                        dp_pion = 0.00273402 + (0.00440449) * 1.8 + (-0.00373488) * Math.pow(1.8, 2) + (0.000996612) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                    }
                    if (pip_theta > 34 && pip_theta < 35) {
                        dp_pion = 0.00333542 + (0.00439874) * pion_p + (-0.00397776) * Math.pow(pion_p, 2) + (0.00105586) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 34 && pip_theta < 35 && pion_p >= 1.8) {
                        dp_pion = 0.00333542 + (0.00439874) * 1.8 + (-0.00397776) * Math.pow(1.8, 2) + (0.00105586) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                    }
                    if (pip_theta > 35 && pip_theta < 36) {
                        dp_pion = 0.00354663 + (0.00565397) * pion_p + (-0.00513503) * Math.pow(pion_p, 2) + (0.00153346) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 35 && pip_theta < 36 && pion_p >= 1.8) {
                        dp_pion = 0.00354663 + (0.00565397) * 1.8 + (-0.00513503) * Math.pow(1.8, 2) + (0.00153346) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                    }
                    if (pip_theta > 36 && pip_theta < 37) {
                        dp_pion = 0.00333909 + (0.00842367) * pion_p + (-0.0077321) * Math.pow(pion_p, 2) + (0.0022489) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 36 && pip_theta < 37 && pion_p >= 1.8) {
                        dp_pion = 0.00333909 + (0.00842367) * 1.8 + (-0.0077321) * Math.pow(1.8, 2) + (0.0022489) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                    }
                    if (pip_theta > 37 && pip_theta < 38) {
                        dp_pion = 0.00358828 + (0.0112108) * pion_p + (-0.0133854) * Math.pow(pion_p, 2) + (0.00486924) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 37 && pip_theta < 38 && pion_p >= 1.4) {
                        dp_pion = 0.00358828 + (0.0112108) * 1.4 + (-0.0133854) * Math.pow(1.4, 2) + (0.00486924) * Math.pow(1.4, 3) + (0) * Math.pow(1.4, 4);
                    }
                    if (pip_theta > 38 && pip_theta < 39) {
                        dp_pion = 0.00354343 + (0.0117121) * pion_p + (-0.0129649) * Math.pow(pion_p, 2) + (0.00455602) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 38 && pip_theta < 39 && pion_p >= 1.3) {
                        dp_pion = 0.00354343 + (0.0117121) * 1.3 + (-0.0129649) * Math.pow(1.3, 2) + (0.00455602) * Math.pow(1.3, 3) + (0) * Math.pow(1.3, 4);
                    }
                    if (pip_theta > 39 && pip_theta < 40) {
                        dp_pion = -0.00194951 + (0.0409713) * pion_p + (-0.0595861) * Math.pow(pion_p, 2) + (0.0281588) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 39 && pip_theta < 40 && pion_p >= 0.9) {
                        dp_pion = -0.00194951 + (0.0409713) * 0.9 + (-0.0595861) * Math.pow(0.9, 2) + (0.0281588) * Math.pow(0.9, 3) + (0) * Math.pow(0.9, 4);
                    }
                    if (pip_theta > 40 && pip_theta < 41) {
                        dp_pion = -0.0099217 + (0.0808096) * pion_p + (-0.119836) * Math.pow(pion_p, 2) + (0.0559553) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 40 && pip_theta < 41 && pion_p >= 0.75) {
                        dp_pion = -0.0099217 + (0.0808096) * 0.75 + (-0.119836) * Math.pow(0.75, 2) + (0.0559553) * Math.pow(0.75, 3) + (0) * Math.pow(0.75, 4);
                    }
                    if (pip_theta > 41 && pip_theta < 42) {
                        dp_pion = 0.00854898 + (0.00025037) * pion_p + (-0.0113992) * Math.pow(pion_p, 2) + (0.0145178) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 41 && pip_theta < 42 && pion_p >= 0.65) {
                        dp_pion = 0.00854898 + (0.00025037) * 0.65 + (-0.0113992) * Math.pow(0.65, 2) + (0.0145178) * Math.pow(0.65, 3) + (0) * Math.pow(0.65, 4);
                    }
                    if (pip_theta > 42) {
                        dp_pion = 0.00564818 + (0.00706606) * pion_p + (0.0042602) * Math.pow(pion_p, 2) + (-0.01141) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 42 && pion_p >= 0.65) {
                        dp_pion = 0.00564818 + (0.00706606) * 0.65 + (0.0042602) * Math.pow(0.65, 2) + (-0.01141) * Math.pow(0.65, 3) + (0) * Math.pow(0.65, 4);
                    }
                } else if (runnum >= 5423 && runnum <= 5666) { // RGA outbending
                    if (pip_theta < 27) {
                        dp_pion = 0.00389945 + (-0.004062) * pion_p + (0.00321842) * Math.pow(pion_p, 2) + (-0.000698299) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta < 27 && pion_p >= 2.3) {
                        dp_pion = 0.00389945 + (-0.004062) * 2.3 + (0.00321842) * Math.pow(2.3, 2) + (-0.000698299) * Math.pow(2.3, 3) + (0) * Math.pow(2.3, 4);
                    }
                    if (pip_theta > 27 && pip_theta < 28) {
                        dp_pion = 0.00727132 + (-0.0117989) * pion_p + (0.00962999) * Math.pow(pion_p, 2) + (-0.00267005) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 27 && pip_theta < 28 && pion_p >= 1.7) {
                        dp_pion = 0.00727132 + (-0.0117989) * 1.7 + (0.00962999) * Math.pow(1.7, 2) + (-0.00267005) * Math.pow(1.7, 3) + (0) * Math.pow(1.7, 4);
                    }
                    if (pip_theta > 28 && pip_theta < 29) {
                        dp_pion = 0.00844551 + (-0.0128097) * pion_p + (0.00945956) * Math.pow(pion_p, 2) + (-0.00237992) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 28 && pip_theta < 29 && pion_p >= 2) {
                        dp_pion = 0.00844551 + (-0.0128097) * 2 + (0.00945956) * Math.pow(2, 2) + (-0.00237992) * Math.pow(2, 3) + (0) * Math.pow(2, 4);
                    }
                    if (pip_theta > 29 && pip_theta < 30) {
                        dp_pion = 0.00959007 + (-0.0139218) * pion_p + (0.0122966) * Math.pow(pion_p, 2) + (-0.0034012) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 29 && pip_theta < 30 && pion_p >= 1.9) {
                        dp_pion = 0.00959007 + (-0.0139218) * 1.9 + (0.0122966) * Math.pow(1.9, 2) + (-0.0034012) * Math.pow(1.9, 3) + (0) * Math.pow(1.9, 4);
                    }
                    if (pip_theta > 30 && pip_theta < 31) {
                        dp_pion = 0.00542816 + (-5.10739e-05) * pion_p + (0.000572038) * Math.pow(pion_p, 2) + (-0.000488883) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 30 && pip_theta < 31 && pion_p >= 1.9) {
                        dp_pion = 0.00542816 + (-5.10739e-05) * 1.9 + (0.000572038) * Math.pow(1.9, 2) + (-0.000488883) * Math.pow(1.9, 3) + (0) * Math.pow(1.9, 4);
                    }
                    if (pip_theta > 31 && pip_theta < 32) {
                        dp_pion = 0.0060391 + (-0.000516936) * pion_p + (-0.00286595) * Math.pow(pion_p, 2) + (0.00136604) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 31 && pip_theta < 32 && pion_p >= 1.8) {
                        dp_pion = 0.0060391 + (-0.000516936) * 1.8 + (-0.00286595) * Math.pow(1.8, 2) + (0.00136604) * Math.pow(1.8, 3) + (0) * Math.pow(1.8, 4);
                    }
                    if (pip_theta > 32 && pip_theta < 33) {
                        dp_pion = 0.0140305 + (-0.0285832) * pion_p + (0.0248799) * Math.pow(pion_p, 2) + (-0.00701311) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 32 && pip_theta < 33 && pion_p >= 1.6) {
                        dp_pion = 0.0140305 + (-0.0285832) * 1.6 + (0.0248799) * Math.pow(1.6, 2) + (-0.00701311) * Math.pow(1.6, 3) + (0) * Math.pow(1.6, 4);
                    }
                    if (pip_theta > 33 && pip_theta < 34) {
                        dp_pion = 0.010815 + (-0.0194244) * pion_p + (0.0174474) * Math.pow(pion_p, 2) + (-0.0049764) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 33 && pip_theta < 34 && pion_p >= 1.5) {
                        dp_pion = 0.010815 + (-0.0194244) * 1.5 + (0.0174474) * Math.pow(1.5, 2) + (-0.0049764) * Math.pow(1.5, 3) + (0) * Math.pow(1.5, 4);
                    }
                    if (pip_theta > 34 && pip_theta < 35) {
                        dp_pion = 0.0105522 + (-0.0176248) * pion_p + (0.0161142) * Math.pow(pion_p, 2) + (-0.00472288) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 34 && pip_theta < 35 && pion_p >= 1.6) {
                        dp_pion = 0.0105522 + (-0.0176248) * 1.6 + (0.0161142) * Math.pow(1.6, 2) + (-0.00472288) * Math.pow(1.6, 3) + (0) * Math.pow(1.6, 4);
                    }
                    if (pip_theta > 35 && pip_theta < 36) {
                        dp_pion = 0.0103938 + (-0.0164003) * pion_p + (0.0164045) * Math.pow(pion_p, 2) + (-0.00517012) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 35 && pip_theta < 36 && pion_p >= 1.5) {
                        dp_pion = 0.0103938 + (-0.0164003) * 1.5 + (0.0164045) * Math.pow(1.5, 2) + (-0.00517012) * Math.pow(1.5, 3) + (0) * Math.pow(1.5, 4);
                    }
                    if (pip_theta > 36 && pip_theta < 37) {
                        dp_pion = 0.0441471 + (-0.183937) * pion_p + (0.338784) * Math.pow(pion_p, 2) + (-0.298985) * Math.pow(pion_p, 3) + (0.126905) * Math.pow(pion_p, 4) + (-0.0208286) * Math.pow(pion_p, 5);
                    }
                    if (pip_theta > 36 && pip_theta < 37 && pion_p >= 1.8) {
                        dp_pion = 0.0441471 + (-0.183937) * 1.8 + (0.338784) * Math.pow(1.8, 2) + (-0.298985) * Math.pow(1.8, 3) + (0.126905) * Math.pow(1.8, 4) + (-0.0208286) * Math.pow(1.8, 5);
                    }
                    if (pip_theta > 37 && pip_theta < 38) {
                        dp_pion = 0.0726119 + (-0.345004) * pion_p + (0.697789) * Math.pow(pion_p, 2) + (-0.685948) * Math.pow(pion_p, 3) + (0.327195) * Math.pow(pion_p, 4) + (-0.0605621) * Math.pow(pion_p, 5);
                    }
                    if (pip_theta > 37 && pip_theta < 38 && pion_p >= 1.7) {
                        dp_pion = 0.0726119 + (-0.345004) * 1.7 + (0.697789) * Math.pow(1.7, 2) + (-0.685948) * Math.pow(1.7, 3) + (0.327195) * Math.pow(1.7, 4) + (-0.0605621) * Math.pow(1.7, 5);
                    }
                    if (pip_theta > 38 && pip_theta < 39) {
                        dp_pion = 0.0247648 + (-0.0797376) * pion_p + (0.126535) * Math.pow(pion_p, 2) + (-0.086545) * Math.pow(pion_p, 3) + (0.0219304) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 38 && pip_theta < 39 && pion_p >= 1.6) {
                        dp_pion = 0.0247648 + (-0.0797376) * 1.6 + (0.126535) * Math.pow(1.6, 2) + (-0.086545) * Math.pow(1.6, 3) + (0.0219304) * Math.pow(1.6, 4);
                    }
                    if (pip_theta > 39 && pip_theta < 40) {
                        dp_pion = 0.0208867 + (-0.0492068) * pion_p + (0.0543187) * Math.pow(pion_p, 2) + (-0.0183393) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 39 && pip_theta < 40 && pion_p >= 1.2) {
                        dp_pion = 0.0208867 + (-0.0492068) * 1.2 + (0.0543187) * Math.pow(1.2, 2) + (-0.0183393) * Math.pow(1.2, 3) + (0) * Math.pow(1.2, 4);
                    }
                    if (pip_theta > 40 && pip_theta < 41) {
                        dp_pion = 0.0148655 + (-0.0203483) * pion_p + (0.00835867) * Math.pow(pion_p, 2) + (0.00697134) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 40 && pip_theta < 41 && pion_p >= 1.0) {
                        dp_pion = 0.0148655 + (-0.0203483) * 1.0 + (0.00835867) * Math.pow(1.0, 2) + (0.00697134) * Math.pow(1.0, 3) + (0) * Math.pow(1.0, 4);
                    }
                    if (pip_theta > 41 && pip_theta < 42) {
                        dp_pion = 0.0223585 + (-0.0365262) * pion_p + (-0.0150027) * Math.pow(pion_p, 2) + (0.0854164) * Math.pow(pion_p, 3) + (-0.0462718) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 41 && pip_theta < 42 && pion_p >= 0.7) {
                        dp_pion = 0.007617;
                    }
                    if (pip_theta > 42) {
                        dp_pion = 0.0152373 + (-0.0106377) * pion_p + (-0.0257573) * Math.pow(pion_p, 2) + (0.0344851) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 42 && pion_p >= 0.75) {
                        dp_pion = 0.0152373 + (-0.0106377) * 0.75 + (-0.0257573) * Math.pow(0.75, 2) + (0.0344851) * Math.pow(0.75, 3) + (0) * Math.pow(0.75, 4);
                    }
                }

            } else if (isCentralDetector) {
                if ((runnum >= 4763 && runnum <= 5419) || (runnum >= 6616 && runnum <= 6783)) { // RGA inbending
                    if (pip_theta < 39) {
                        dp_pion = -0.045 + (-0.102652) + (0.455589) * pion_p + (-0.671635) * Math.pow(pion_p, 2) + (0.303814) * Math.pow(pion_p, 3);
                    }
                    if (pip_theta < 39 && pion_p >= 0.7) {
                        dp_pion = -0.045 + (-0.102652) + (0.455589) * 0.7 + (-0.671635) * Math.pow(0.7, 2) + (0.303814) * Math.pow(0.7, 3);
                    }
                    if (pip_theta > 39 && pip_theta < 40) {
                        dp_pion = 0.0684552 + (-0.766492) * pion_p + (1.73092) * Math.pow(pion_p, 2) + (-1.46215) * Math.pow(pion_p, 3) + (0.420127) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 39 && pip_theta < 40 && pion_p >= 1.4) {
                        dp_pion = 0.0684552 + (-0.766492) * 1.4 + (1.73092) * Math.pow(1.4, 2) + (-1.46215) * Math.pow(1.4, 3) + (0.420127) * Math.pow(1.4, 4);
                    }
                    if (pip_theta > 40 && pip_theta < 41) {
                        dp_pion = 0.751549 + (-7.4593) * pion_p + (26.8037) * Math.pow(pion_p, 2) + (-47.1576) * Math.pow(pion_p, 3) + (43.8527) * Math.pow(pion_p, 4) + (-20.7039) * Math.pow(pion_p, 5) + (3.90931) * Math.pow(pion_p, 6);
                    }
                    if (pip_theta > 40 && pip_theta < 41 && pion_p >= 1.45) {
                        dp_pion = 0.751549 + (-7.4593) * 1.45 + (26.8037) * Math.pow(1.45, 2) + (-47.1576) * Math.pow(1.45, 3) + (43.8527) * Math.pow(1.45, 4) + (-20.7039) * Math.pow(1.45, 5) + (3.90931) * Math.pow(1.45, 6);
                    }
                    if (pip_theta > 41 && pip_theta < 42) {
                        dp_pion = -1.35043 + (10.0788) * pion_p + (-30.4829) * Math.pow(pion_p, 2) + (47.7792) * Math.pow(pion_p, 3) + (-40.996) * Math.pow(pion_p, 4) + (18.2662) * Math.pow(pion_p, 5) + (-3.30449) * Math.pow(pion_p, 6);
                    }
                    if (pip_theta > 41 && pip_theta < 42 && pion_p >= 1.2) {
                        dp_pion = -1.35043 + (10.0788) * 1.2 + (-30.4829) * Math.pow(1.2, 2) + (47.7792) * Math.pow(1.2, 3) + (-40.996) * Math.pow(1.2, 4) + (18.2662) * Math.pow(1.2, 5) + (-3.30449) * Math.pow(1.2, 6);
                    }
                    if (pip_theta > 42 && pip_theta < 43) {
                        dp_pion = -0.0231195 + (0.0744589) * pion_p + (-0.0807029) * Math.pow(pion_p, 2) + (0.0264266) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 42 && pip_theta < 43 && pion_p >= 1.3) {
                        dp_pion = -0.0231195 + (0.0744589) * 1.3 + (-0.0807029) * Math.pow(1.3, 2) + (0.0264266) * Math.pow(1.3, 3) + (0) * Math.pow(1.3, 4);
                    }
                    if (pip_theta > 43 && pip_theta < 44) {
                        dp_pion = -0.00979928 + (0.0351043) * pion_p + (-0.0365865) * Math.pow(pion_p, 2) + (0.00977218) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 43 && pip_theta < 44 && pion_p >= 1.1) {
                        dp_pion = -0.00979928 + (0.0351043) * 1.1 + (-0.0365865) * Math.pow(1.1, 2) + (0.00977218) * Math.pow(1.1, 3) + (0) * Math.pow(1.1, 4);
                    }
                    if (pip_theta > 44 && pip_theta < 45) {
                        dp_pion = 0.00108491 + (-0.00924885) * pion_p + (0.0216431) * Math.pow(pion_p, 2) + (-0.0137762) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 44 && pip_theta < 45 && pion_p >= 1.1) {
                        dp_pion = 0.00108491 + (-0.00924885) * 1.1 + (0.0216431) * Math.pow(1.1, 2) + (-0.0137762) * Math.pow(1.1, 3) + (0) * Math.pow(1.1, 4);
                    }
                    if (pip_theta > 45 && pip_theta < 55) {
                        dp_pion = 0.0092263 + (-0.0676178) * pion_p + (0.168778) * Math.pow(pion_p, 2) + (-0.167463) * Math.pow(pion_p, 3) + (0.05661) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 45 && pip_theta < 55 && pion_p >= 1.3) {
                        dp_pion = 0.0092263 + (-0.0676178) * 1.3 + (0.168778) * Math.pow(1.3, 2) + (-0.167463) * Math.pow(1.3, 3) + (0.05661) * Math.pow(1.3, 4);
                    }
                    if (pip_theta > 55 && pip_theta < 65) {
                        dp_pion = 0.00805642 + (-0.0670962) * pion_p + (0.188536) * Math.pow(pion_p, 2) + (-0.20571) * Math.pow(pion_p, 3) + (0.0765) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 55 && pip_theta < 65 && pion_p >= 1.05) {
                        dp_pion = 0.00805642 + (-0.0670962) * 1.05 + (0.188536) * Math.pow(1.05, 2) + (-0.20571) * Math.pow(1.05, 3) + (0.0765) * Math.pow(1.05, 4);
                    }
                    if (pip_theta > 65 && pip_theta < 75) {
                        dp_pion = 0.00312202 + (-0.0269717) * pion_p + (0.0715236) * Math.pow(pion_p, 2) + (-0.0545622) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 65 && pip_theta < 75 && pion_p >= 0.75) {
                        dp_pion = 0.00312202 + (-0.0269717) * 0.75 + (0.0715236) * Math.pow(0.75, 2) + (-0.0545622) * Math.pow(0.75, 3) + (0) * Math.pow(0.75, 4);
                    }
                    if (pip_theta > 75 && pip_theta < 85) {
                        dp_pion = 0.00424971 + (-0.0367683) * pion_p + (0.10417) * Math.pow(pion_p, 2) + (-0.0899651) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 75 && pip_theta < 85 && pion_p >= 0.65) {
                        dp_pion = 0.00424971 + (-0.0367683) * 0.65 + (0.10417) * Math.pow(0.65, 2) + (-0.0899651) * Math.pow(0.65, 3) + (0) * Math.pow(0.65, 4);
                    }
                    if (pip_theta > 85 && pip_theta < 95) {
                        dp_pion = 0.00654123 + (-0.0517915) * pion_p + (0.147888) * Math.pow(pion_p, 2) + (-0.14253) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 85 && pip_theta < 95 && pion_p >= 0.5) {
                        dp_pion = 0.00654123 + (-0.0517915) * 0.5 + (0.147888) * Math.pow(0.5, 2) + (-0.14253) * Math.pow(0.5, 3) + (0) * Math.pow(0.5, 4);
                    }
                    if (pip_theta > 95 && pip_theta < 105) {
                        dp_pion = -0.00111721 + (0.00478119) * pion_p + (0.0158753) * Math.pow(pion_p, 2) + (-0.052902) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 95 && pip_theta < 105 && pion_p >= 0.45) {
                        dp_pion = -0.00111721 + (0.00478119) * 0.45 + (0.0158753) * Math.pow(0.45, 2) + (-0.052902) * Math.pow(0.45, 3) + (0) * Math.pow(0.45, 4);
                    }
                    if (pip_theta > 105 && pip_theta < 115) {
                        dp_pion = -0.00239839 + (0.00790738) * pion_p + (0.0311713) * Math.pow(pion_p, 2) + (-0.104157) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 105 && pip_theta < 115 && pion_p >= 0.35) {
                        dp_pion = -0.00239839 + (0.00790738) * 0.35 + (0.0311713) * Math.pow(0.35, 2) + (-0.104157) * Math.pow(0.35, 3) + (0) * Math.pow(0.35, 4);
                    }
                    if (pip_theta > 115 && pip_theta < 125) {
                        dp_pion = -0.00778793 + (0.0256774) * pion_p + (0.0932503) * Math.pow(pion_p, 2) + (-0.32771) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 115 && pip_theta < 125 && pion_p >= 0.35) {
                        dp_pion = -0.00778793 + (0.0256774) * 0.35 + (0.0932503) * Math.pow(0.35, 2) + (-0.32771) * Math.pow(0.35, 3) + (0) * Math.pow(0.35, 4);
                    }
                    if (pip_theta > 125 && pip_theta < 135) {
                        dp_pion = -0.00292778 + (-0.00536697) * pion_p + (-0.00414351) * Math.pow(pion_p, 2) + (0.0196431) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                    }
                    if (pip_theta > 125 && pip_theta < 135 && pion_p >= 0.35) {
                        dp_pion = -0.00292778 + (-0.00536697) * 0.35 + (-0.00414351) * Math.pow(0.35, 2) + (0.0196431) * Math.pow(0.35, 3) + (0) * Math.pow(0.35, 4);
                    }
                }
            } else if (runnum >= 5423 && runnum <= 5666) { // RGA outbending
                if (pip_theta < 39) {
                    dp_pion = -0.05 + (-0.0758897) + (0.362231) * pion_p + (-0.542404) * Math.pow(pion_p, 2) + (0.241344) * Math.pow(pion_p, 3);
                }
                if (pip_theta < 39 && pion_p >= 0.8) {
                    dp_pion = -0.05 + (-0.0758897) + (0.362231) * 0.8 + (-0.542404) * Math.pow(0.8, 2) + (0.241344) * Math.pow(0.8, 3);
                }
                if (pip_theta > 39 && pip_theta < 40) {
                    dp_pion = 0.0355259 + (-0.589712) * pion_p + (1.4206) * Math.pow(pion_p, 2) + (-1.24179) * Math.pow(pion_p, 3) + (0.365524) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 39 && pip_theta < 40 && pion_p >= 1.35) {
                    dp_pion = 0.0355259 + (-0.589712) * 1.35 + (1.4206) * Math.pow(1.35, 2) + (-1.24179) * Math.pow(1.35, 3) + (0.365524) * Math.pow(1.35, 4);
                }
                if (pip_theta > 40 && pip_theta < 41) {
                    dp_pion = -0.252336 + (1.02032) * pion_p + (-1.51461) * Math.pow(pion_p, 2) + (0.967772) * Math.pow(pion_p, 3) + (-0.226028) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 40 && pip_theta < 41 && pion_p >= 1.4) {
                    dp_pion = -0.252336 + (1.02032) * 1.4 + (-1.51461) * Math.pow(1.4, 2) + (0.967772) * Math.pow(1.4, 3) + (-0.226028) * Math.pow(1.4, 4);
                }
                if (pip_theta > 41 && pip_theta < 42) {
                    dp_pion = -0.710129 + (4.49613) * pion_p + (-11.01) * Math.pow(pion_p, 2) + (12.9945) * Math.pow(pion_p, 3) + (-7.41641) * Math.pow(pion_p, 4) + (1.63923) * Math.pow(pion_p, 5);
                }
                if (pip_theta > 41 && pip_theta < 42 && pion_p >= 1.2) {
                    dp_pion = -0.710129 + (4.49613) * 1.2 + (-11.01) * Math.pow(1.2, 2) + (12.9945) * Math.pow(1.2, 3) + (-7.41641) * Math.pow(1.2, 4) + (1.63923) * Math.pow(1.2, 5);
                }
                if (pip_theta > 42 && pip_theta < 43) {
                    dp_pion = -0.0254912 + (0.0851432) * pion_p + (-0.0968583) * Math.pow(pion_p, 2) + (0.0350334) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 42 && pip_theta < 43 && pion_p >= 1.2) {
                    dp_pion = -0.0254912 + (0.0851432) * 1.2 + (-0.0968583) * Math.pow(1.2, 2) + (0.0350334) * Math.pow(1.2, 3) + (0) * Math.pow(1.2, 4);
                }
                if (pip_theta > 43 && pip_theta < 44) {
                    dp_pion = -0.0115965 + (0.0438726) * pion_p + (-0.0500474) * Math.pow(pion_p, 2) + (0.0163627) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 43 && pip_theta < 44 && pion_p >= 1.4) {
                    dp_pion = -0.0115965 + (0.0438726) * 1.4 + (-0.0500474) * Math.pow(1.4, 2) + (0.0163627) * Math.pow(1.4, 3) + (0) * Math.pow(1.4, 4);
                }
                if (pip_theta > 44 && pip_theta < 45) {
                    dp_pion = 0.00273414 + (-0.01851) * pion_p + (0.0377032) * Math.pow(pion_p, 2) + (-0.0226696) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 44 && pip_theta < 45 && pion_p >= 1) {
                    dp_pion = 0.00273414 + (-0.01851) * 1 + (0.0377032) * Math.pow(1, 2) + (-0.0226696) * Math.pow(1, 3) + (0) * Math.pow(1, 4);
                }
                if (pip_theta > 45 && pip_theta < 55) {
                    dp_pion = 0.0271952 + (-0.25981) * pion_p + (0.960051) * Math.pow(pion_p, 2) + (-1.76651) * Math.pow(pion_p, 3) + (1.72872) * Math.pow(pion_p, 4) + (-0.856946) * Math.pow(pion_p, 5) + (0.167564) * Math.pow(pion_p, 6);
                }
                if (pip_theta > 45 && pip_theta < 55 && pion_p >= 1.4) {
                    dp_pion = 0.0271952 + (-0.25981) * 1.4 + (0.960051) * Math.pow(1.4, 2) + (-1.76651) * Math.pow(1.4, 3) + (1.72872) * Math.pow(1.4, 4) + (-0.856946) * Math.pow(1.4, 5) + (0.167564) * Math.pow(1.4, 6);
                }
                if (pip_theta > 55 && pip_theta < 65) {
                    dp_pion = 0.00734975 + (-0.0598841) * pion_p + (0.161495) * Math.pow(pion_p, 2) + (-0.1629) * Math.pow(pion_p, 3) + (0.0530098) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 55 && pip_theta < 65 && pion_p >= 1.2) {
                    dp_pion = 0.00734975 + (-0.0598841) * 1.2 + (0.161495) * Math.pow(1.2, 2) + (-0.1629) * Math.pow(1.2, 3) + (0.0530098) * Math.pow(1.2, 4);
                }
                if (pip_theta > 65 && pip_theta < 75) {
                    dp_pion = 0.00321351 + (-0.0289322) * pion_p + (0.0786484) * Math.pow(pion_p, 2) + (-0.0607041) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 65 && pip_theta < 75 && pion_p >= 0.95) {
                    dp_pion = 0.00321351 + (-0.0289322) * 0.95 + (0.0786484) * Math.pow(0.95, 2) + (-0.0607041) * Math.pow(0.95, 3) + (0) * Math.pow(0.95, 4);
                }
                if (pip_theta > 75 && pip_theta < 85) {
                    dp_pion = 0.00644253 + (-0.0543896) * pion_p + (0.148933) * Math.pow(pion_p, 2) + (-0.1256) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 75 && pip_theta < 85 && pion_p >= 0.7) {
                    dp_pion = 0.00644253 + (-0.0543896) * 0.7 + (0.148933) * Math.pow(0.7, 2) + (-0.1256) * Math.pow(0.7, 3) + (0) * Math.pow(0.7, 4);
                }
                if (pip_theta > 85 && pip_theta < 95) {
                    dp_pion = 0.00671152 + (-0.0537269) * pion_p + (0.154509) * Math.pow(pion_p, 2) + (-0.147667) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 85 && pip_theta < 95 && pion_p >= 0.65) {
                    dp_pion = 0.00671152 + (-0.0537269) * 0.65 + (0.154509) * Math.pow(0.65, 2) + (-0.147667) * Math.pow(0.65, 3) + (0) * Math.pow(0.65, 4);
                }
                if (pip_theta > 95 && pip_theta < 105) {
                    dp_pion = -0.000709077 + (0.00331818) * pion_p + (0.0109241) * Math.pow(pion_p, 2) + (-0.0351682) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 95 && pip_theta < 105 && pion_p >= 0.45) {
                    dp_pion = -0.000709077 + (0.00331818) * 0.45 + (0.0109241) * Math.pow(0.45, 2) + (-0.0351682) * Math.pow(0.45, 3) + (0) * Math.pow(0.45, 4);
                }
                if (pip_theta > 105 && pip_theta < 115) {
                    dp_pion = -0.00260164 + (0.00846919) * pion_p + (0.0315497) * Math.pow(pion_p, 2) + (-0.105756) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 105 && pip_theta < 115 && pion_p >= 0.45) {
                    dp_pion = -0.00260164 + (0.00846919) * 0.45 + (0.0315497) * Math.pow(0.45, 2) + (-0.105756) * Math.pow(0.45, 3) + (0) * Math.pow(0.45, 4);
                }
                if (pip_theta > 115 && pip_theta < 125) {
                    dp_pion = -0.00544336 + (0.018256) * pion_p + (0.0664618) * Math.pow(pion_p, 2) + (-0.240312) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 115 && pip_theta < 125 && pion_p >= 0.45) {
                    dp_pion = -0.00544336 + (0.018256) * 0.45 + (0.0664618) * Math.pow(0.45, 2) + (-0.240312) * Math.pow(0.45, 3) + (0) * Math.pow(0.45, 4);
                }
                if (pip_theta > 125 && pip_theta < 135) {
                    dp_pion = -0.00281073 + (-0.00495863) * pion_p + (-0.00362356) * Math.pow(pion_p, 2) + (0.0178764) * Math.pow(pion_p, 3) + (0) * Math.pow(pion_p, 4);
                }
                if (pip_theta > 125 && pip_theta < 135 && pion_p >= 0.35) {
                    dp_pion = -0.00281073 + (-0.00495863) * 0.35 + (-0.00362356) * Math.pow(0.35, 2) + (0.0178764) * Math.pow(0.35, 3) + (0) * Math.pow(0.35, 4);
                }
            }
        }

        // Update the px, py, pz values
        p_array[0] = (float) x_calculation(p, theta, phi);
        p_array[1] = (float) y_calculation(p, theta, phi);
        p_array[2] = (float) z_calculation(p, theta);
    }

    
    public static void mariana_electron_energy_loss_corrections(int particle_Index, float[] p_array,
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
        boolean validRun = (runnum >= 4763 && runnum <= 5419)
                || (runnum >= 5423 && runnum <= 5666)
                || (runnum >= 6616 && runnum <= 6783);
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
            
            if ((runnum >= 4763 && runnum <= 5419) || (runnum >= 5423 && runnum <= 5666)) {
                E_new = E + 0.0208922 + 0.050158 * Math.pow(E,1) - 0.0181107 * Math.pow(E,2) + 
                        0.00305671 * Math.pow(E,3) - 0.000178235 * Math.pow(E,4);
            } else if (runnum >= 6616 && runnum <= 6783) {
                E_new = E + 0.085643 - 0.0288063 * Math.pow(E,1) + 0.00894691 * Math.pow(E,2) - 
                        0.000725449 * Math.pow(E,3);
            }
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
    
    // needs further validation, seems to not improve anything
    public static void sebastian_electron_energy_loss_corrections(int particle_Index, float[] p_array,
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
        boolean validRun = (runnum >= 4763 && runnum <= 5419)
                || (runnum >= 5423 && runnum <= 5666)
                || (runnum >= 6616 && runnum <= 6783);
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
            double[] corIn_PDep = {-9.538e-03, 6.661e-03, -3.333e-04, -6.136e-03, 3.611e-02};
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

    // needs further validation, seems to not improve anything
    public static void sebastian_photon_energy_loss_corrections(int particle_Index, float[] p_array,
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
        boolean validRun = (runnum >= 4763 && runnum <= 5419)
                || (runnum >= 5423 && runnum <= 5666)
                || (runnum >= 6616 && runnum <= 6783);
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
