package extended_kinematic_fitters;

import static extended_kinematic_fitters.energy_loss_corrections.p_calculation;
import static extended_kinematic_fitters.energy_loss_corrections.phi_calculation;
import static extended_kinematic_fitters.energy_loss_corrections.theta_calculation;
import static extended_kinematic_fitters.energy_loss_corrections.x_calculation;
import static extended_kinematic_fitters.energy_loss_corrections.y_calculation;
import static extended_kinematic_fitters.energy_loss_corrections.z_calculation;
import org.jlab.io.hipo.HipoDataBank;
import java.util.Arrays;
import java.util.Random;
import org.jlab.clas.physics.Vector3;

/**
 *
 * @author tbhayward (translated form Richard Capobianco's C++)
 */
public class momentum_corrections {

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

    public void inbending_momentum_corrections(float[] p_array, int sec, int ivec,
            int corEl, int corPip, int corPim, int corPro) {
        // 'Px'/'Py'/'Pz'   ==> Corresponds to the Cartesian Components of the particle momentum being corrected
        // 'sec'            ==> Corresponds to the Forward Detector Sectors where the given particle is detected (6 total)
        // 'ivec'           ==> Corresponds to the particle being corrected (See below)    
        // (*) ivec = 0 --> Electron Corrections
        // (*) ivec = 1 --> Pi+ Corrections
        // (*) ivec = 2 --> Pi- Corrections
        // (*) ivec = 3 --> Proton Corrections
        // 'corEl'/'corPip'/'corPim'/'corPro' ==> Controls which version of the particle correction is used
        // Includes:
        // (*) Correction On/Off
        // (*) Pass Version
        // (*) Data Set (Fall 2018 or Spring 2019)
        // 'corEl'         ==> Controls the ELECTRON Corrections
        // corEl == 0  --> No Correction (Off)
        // corEl == 1  --> Fall  2018 - Pass 1
        // corEl == 2  --> Sping 2019 - Pass 2
        // corEl == 3  --> Fall  2018 - Pass 2
        // 'corPip'        ==> Controls the π+ PION Corrections
        // corPip == 0 --> No Correction
        // corPip == 1 --> Fall  2018 - Pass 1
        // corPip == 2 --> Sping 2019 - Pass 2
        // corPip == 3 --> Fall  2018 - Pass 2
        // 'corPim'        ==> Controls the π- PION Corrections
        // corPim == 0 --> No Correction
        // corPim == 1 --> Fall  2018 - Pass 1 (Created by Nick Trotta)
        // 'corPro'        ==> Controls the PROTON Corrections (Momentum)
        // corPro == 0 --> No Correction
        // corPro == 1 --> Fall  2018 - Pass 1

        double Px = p_array[0];
        double Py = p_array[1];
        double Pz = p_array[2];
        // Momentum Magnitude
        double pp = Math.sqrt(Px * Px + Py * Py + Pz * Pz);

        // Initializing the correction factor
        double dp = 0;

        // Defining Phi Angle
        double Phi = (180 / 3.1415926) * Math.atan2(Py, Px);

        // Central Detector Corrections Not Included (Yet)
        // (Initial) Shift of the Phi Angle (done to realign sectors whose data is separated when plotted from ±180˚)
        if (((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)) {
            Phi += 360;
        }

        // Getting Local Phi Angle
        double PhiLocal = Phi - (sec - 1) * 60;

        // Applying Shift Functions to Phi Angles (local shifted phi = phi)
        double phi = PhiLocal;

        // For Electron Shift
        if (ivec == 0) {
            phi = PhiLocal - 30 / pp;
        }

        // For π+ Pion/Proton Shift
        if (ivec == 1 || ivec == 3) {
            phi = PhiLocal + (32 / (pp - 0.05));
        }

        // For π- Pion Shift
        if (ivec == 2) {
            phi = PhiLocal - (32 / (pp - 0.05));
        }

        //===============//===============//     No Corrections     //===============//===============//
        if (corEl == 0 && ivec == 0) { // No Electron Correction
            dp = 0;
        }
        if (corPip == 0 && ivec == 1) { // No π+ Pion Correction
            dp = 0;
        }
        if (corPim == 0 && ivec == 2) { // No π- Pion Correction
            dp = 0;
        }
        if (corPro == 0 && ivec == 3) { // No Proton Correction
            dp = 0;
        }
        //==============//==============//     No Corrections (End)     //==============//==============//

        //==============================//     Electron Corrections     //==============================//
        if (corEl != 0 && ivec == 0) {
            if (corEl == 1) { // Fall 2018 - Pass 1 Corrections
                if (sec == 1) {
                    dp = ((-4.3303e-06) * phi * phi + (1.1006e-04) * phi + (-5.7235e-04)) * pp * pp + ((3.2555e-05) * phi * phi + (-0.0014559) * phi + (0.0014878)) * pp + ((-1.9577e-05) * phi * phi + (0.0017996) * phi + (0.025963));
                }
                if (sec == 2) {
                    dp = ((-9.8045e-07) * phi * phi + (6.7395e-05) * phi + (-4.6757e-05)) * pp * pp + ((-1.4958e-05) * phi * phi + (-0.0011191) * phi + (-0.0025143)) * pp + ((1.2699e-04) * phi * phi + (0.0033121) * phi + (0.020819));
                }
                if (sec == 3) {
                    dp = ((-5.9459e-07) * phi * phi + (-2.8289e-05) * phi + (-4.3541e-04)) * pp * pp + ((-1.5025e-05) * phi * phi + (5.7730e-04) * phi + (-0.0077582)) * pp + ((7.3348e-05) * phi * phi + (-0.001102) * phi + (0.057052));
                }
                if (sec == 4) {
                    dp = ((-2.2714e-06) * phi * phi + (-3.0360e-05) * phi + (-8.9322e-04)) * pp * pp + ((2.9737e-05) * phi * phi + (5.1142e-04) * phi + (0.0045641)) * pp + ((-1.0582e-04) * phi * phi + (-5.6852e-04) * phi + (0.027506));
                }
                if (sec == 5) {
                    dp = ((-1.1490e-06) * phi * phi + (-6.2147e-06) * phi + (-4.7235e-04)) * pp * pp + ((3.7039e-06) * phi * phi + (-1.5943e-04) * phi + (-8.5238e-04)) * pp + ((4.4069e-05) * phi * phi + (0.0014152) * phi + (0.031933));
                }
                if (sec == 6) {
                    dp = ((1.1076e-06) * phi * phi + (4.0156e-05) * phi + (-1.6341e-04)) * pp * pp + ((-2.8613e-05) * phi * phi + (-5.1861e-04) * phi + (-0.0056437)) * pp + ((1.2419e-04) * phi * phi + (4.9084e-04) * phi + (0.049976));
                }
            }

            if (corEl == 2) { // Spring 2019 - Pass 2 Corrections
                if (sec == 1) {
                    dp = ((-1.4215599999999998e-06) * phi * phi + (4.91084e-06) * phi + (-0.00012995999999999998)) * pp * pp + ((1.6952059999999994e-05) * phi * phi + (-0.00033224299999999997) * phi + (-0.0018080400000000003)) * pp + ((-3.1853499999999996e-05) * phi * phi + (0.0016001439999999997) * phi + (0.03187985));
                }
                if (sec == 2) {
                    dp = ((-5.4471e-06) * phi * phi + (-4.69579e-05) * phi + (0.000462807)) * pp * pp + ((5.0258819999999995e-05) * phi * phi + (0.00023192399999999994) * phi + (-0.01118006)) * pp + ((-8.5754e-05) * phi * phi + (0.00017097299999999994) * phi + (0.05324023));
                }
                if (sec == 3) {
                    dp = ((-3.4392460000000005e-06) * phi * phi + (9.860100000000002e-06) * phi + (-3.8414000000000015e-05)) * pp * pp + ((1.7492300000000002e-05) * phi * phi + (-4.111499999999996e-05) * phi + (-0.0052975509999999984)) * pp + ((1.0045499999999984e-05) * phi * phi + (-2.1412000000000004e-05) * phi + (0.03514576));
                }
                if (sec == 4) {
                    dp = ((2.4865599999999998e-06) * phi * phi + (2.9090599999999996e-05) * phi + (0.00016154500000000003)) * pp * pp + ((-2.7148730000000002e-05) * phi * phi + (-0.000136352) * phi + (-0.00543832)) * pp + ((4.917660000000001e-05) * phi * phi + (-0.0001558459999999999) * phi + (0.04322285000000001));
                }
                if (sec == 5) {
                    dp = ((-5.340280000000001e-06) * phi * phi + (1.355319e-05) * phi + (0.001362661)) * pp * pp + ((5.858976999999999e-05) * phi * phi + (-0.00024119909999999995) * phi + (-0.02025752)) * pp + ((-0.0001475504) * phi * phi + (0.0005707250000000001) * phi + (0.07970399));
                }
                if (sec == 6) {
                    dp = ((-3.0325500000000003e-06) * phi * phi + (-4.7810870999999994e-05) * phi + (0.001092504)) * pp * pp + ((2.4123071999999996e-05) * phi * phi + (0.00047091400000000007) * phi + (-0.01504266)) * pp + ((-9.523899999999999e-06) * phi * phi + (-0.0008819019999999998) * phi + (0.048088700000000005));
                }
            }

            if (corEl == 3) { // Fall 2018 - Pass 2 Corrections
                if (sec == 1) {
                    dp = ((-9.82416e-06) * phi * phi + (-2.29956e-05) * phi + (0.00029664199999999996)) * pp * pp + ((0.0001113414) * phi * phi + (-2.041300000000001e-05) * phi + (-0.00862226)) * pp + ((-0.000281738) * phi * phi + (0.00058712) * phi + (0.0652737));
                    if (pp < 7) {
                        dp = dp + ((-3.4001e-06) * phi * phi + (-2.2885e-05) * phi + (9.9705e-04)) * pp * pp + ((2.1840e-05) * phi * phi + (2.4238e-04) * phi + (-0.0091904)) * pp + ((-2.9180e-05) * phi * phi + (-6.4496e-04) * phi + (0.022505));
                    } else {
                        dp = dp + ((-6.3656e-05) * phi * phi + (1.7266e-04) * phi + (-0.0017909)) * pp * pp + ((0.00104) * phi * phi + (-0.0028401) * phi + (0.02981)) * pp + ((-0.0041995) * phi * phi + (0.011537) * phi + (-0.1196));
                    }
                    dp = dp + ((3.2780000000000006e-07) * phi * phi + (6.7084e-07) * phi + (-4.390000000000004e-05)) * pp * pp + ((-7.230999999999999e-06) * phi * phi + (-2.37482e-05) * phi + (0.0004909000000000007)) * pp + ((3.285299999999999e-05) * phi * phi + (9.63723e-05) * phi + (-0.00115));
                }
                if (sec == 2) {
                    dp = ((-7.741952e-06) * phi * phi + (-2.2402167000000004e-05) * phi + (-0.00042652900000000004)) * pp * pp + ((7.54079e-05) * phi * phi + (-1.3333999999999984e-05) * phi + (0.0002420100000000004)) * pp + ((-0.000147876) * phi * phi + (0.00057905) * phi + (0.0253551));
                    if (pp < 7) {
                        dp = dp + ((5.3611e-06) * phi * phi + (8.1979e-06) * phi + (5.9789e-04)) * pp * pp + ((-4.8185e-05) * phi * phi + (-1.5188e-04) * phi + (-0.0084675)) * pp + ((9.2324e-05) * phi * phi + (6.4420e-04) * phi + (0.026792));
                    } else {
                        dp = dp + ((-6.1139e-05) * phi * phi + (5.4087e-06) * phi + (-0.0021284)) * pp * pp + ((0.0010007) * phi * phi + (9.3492e-05) * phi + (0.039813)) * pp + ((-0.0040434) * phi * phi + (-0.0010953) * phi + (-0.18112));
                    }
                    dp = dp + ((6.221217e-07) * phi * phi + (1.9596000000000003e-06) * phi + (-9.826e-05)) * pp * pp + ((-1.28576e-05) * phi * phi + (-4.36589e-05) * phi + (0.00130342)) * pp + ((5.80399e-05) * phi * phi + (0.000215388) * phi + (-0.0040414000000000005));
                }
                if (sec == 3) {
                    dp = ((-5.115364000000001e-06) * phi * phi + (-1.1983000000000004e-05) * phi + (-0.0006832899999999999)) * pp * pp + ((4.52287e-05) * phi * phi + (0.00020855000000000003) * phi + (0.0034986999999999996)) * pp + ((-9.044610000000001e-05) * phi * phi + (-0.00106657) * phi + (0.017954199999999997));
                    if (pp < 7) {
                        dp = dp + ((9.9281e-07) * phi * phi + (3.4879e-06) * phi + (0.0011673)) * pp * pp + ((-2.0071e-05) * phi * phi + (-3.1362e-05) * phi + (-0.012329)) * pp + ((6.9463e-05) * phi * phi + (3.5102e-05) * phi + (0.037505));
                    } else {
                        dp = dp + ((-3.2178e-06) * phi * phi + (4.0630e-05) * phi + (-0.005209)) * pp * pp + ((2.0884e-05) * phi * phi + (-6.8800e-04) * phi + (0.086513)) * pp + ((3.9530e-05) * phi * phi + (0.0029306) * phi + (-0.3507));
                    }
                    dp = dp + ((-4.045999999999999e-07) * phi * phi + (-1.3115999999999994e-06) * phi + (3.9510000000000006e-05)) * pp * pp + ((5.521e-06) * phi * phi + (2.4436999999999997e-05) * phi + (-0.0016887)) * pp + ((-1.0962999999999997e-05) * phi * phi + (-0.000151944) * phi + (0.009313599999999998));
                }
                if (sec == 4) {
                    dp = ((-3.9278116999999996e-06) * phi * phi + (2.2289300000000004e-05) * phi + (0.00012665000000000002)) * pp * pp + ((4.8649299999999995e-05) * phi * phi + (-0.00012554) * phi + (-0.005955500000000001)) * pp + ((-0.00014617199999999997) * phi * phi + (-0.00028571) * phi + (0.0606998));
                    if (pp < 7) {
                        dp = dp + ((-4.8455e-06) * phi * phi + (-1.2074e-05) * phi + (0.0013221)) * pp * pp + ((3.2207e-05) * phi * phi + (1.3144e-04) * phi + (-0.010451)) * pp + ((-3.7365e-05) * phi * phi + (-4.2344e-04) * phi + (0.019952));
                    } else {
                        dp = dp + ((-3.9554e-05) * phi * phi + (5.5496e-06) * phi + (-0.0058293)) * pp * pp + ((6.5077e-04) * phi * phi + (2.6735e-05) * phi + (0.095025)) * pp + ((-0.0026457) * phi * phi + (-6.1394e-04) * phi + (-0.3793));
                    }
                    dp = dp + ((-4.593089e-07) * phi * phi + (1.40673e-05) * phi + (6.69e-05)) * pp * pp + ((4.0239e-06) * phi * phi + (-0.000180863) * phi + (-0.0008272199999999999)) * pp + ((-5.1310000000000005e-06) * phi * phi + (0.00049748) * phi + (0.00255231));
                }
                if (sec == 5) {
                    dp = ((8.036599999999999e-07) * phi * phi + (2.58072e-05) * phi + (0.000360217)) * pp * pp + ((-9.932400000000002e-06) * phi * phi + (-0.0005168531) * phi + (-0.010904)) * pp + ((1.8516299999999998e-05) * phi * phi + (0.0015570900000000001) * phi + (0.066493));
                    if (pp < 7) {
                        dp = dp + ((7.7156e-07) * phi * phi + (-3.9566e-05) * phi + (-2.3589e-04)) * pp * pp + ((-9.8309e-06) * phi * phi + (3.7353e-04) * phi + (0.0020382)) * pp + ((2.9506e-05) * phi * phi + (-8.0409e-04) * phi + (-0.0045615));
                    } else {
                        dp = dp + ((-3.2410e-05) * phi * phi + (-4.3301e-05) * phi + (-0.0028742)) * pp * pp + ((5.3787e-04) * phi * phi + (6.8921e-04) * phi + (0.049578)) * pp + ((-0.0021955) * phi * phi + (-0.0027698) * phi + (-0.21142));
                    }
                    dp = dp + ((-1.2151e-06) * phi * phi + (-8.5087e-06) * phi + (4.968e-05)) * pp * pp + ((1.46998e-05) * phi * phi + (0.000115047) * phi + (-0.00039269)) * pp + ((-4.0368600000000005e-05) * phi * phi + (-0.00037078) * phi + (0.00073998));
                }
                if (sec == 6) {
                    dp = ((-1.9552099999999998e-06) * phi * phi + (8.042199999999997e-06) * phi + (-2.1324000000000028e-05)) * pp * pp + ((1.6969399999999997e-05) * phi * phi + (-6.306600000000001e-05) * phi + (-0.00485568)) * pp + ((-2.7723e-05) * phi * phi + (-6.828400000000003e-05) * phi + (0.0447535));
                    if (pp < 7) {
                        dp = dp + ((-8.2535e-07) * phi * phi + (9.1433e-06) * phi + (3.5395e-04)) * pp * pp + ((-3.4272e-06) * phi * phi + (-1.3012e-04) * phi + (-0.0030724)) * pp + ((4.9211e-05) * phi * phi + (4.5807e-04) * phi + (0.0058932));
                    } else {
                        dp = dp + ((-4.9760e-05) * phi * phi + (-7.2903e-05) * phi + (-0.0020453)) * pp * pp + ((8.0918e-04) * phi * phi + (0.0011688) * phi + (0.037042)) * pp + ((-0.0032504) * phi * phi + (-0.0046169) * phi + (-0.16331));
                    }
                    dp = dp + ((-7.153000000000002e-07) * phi * phi + (1.62859e-05) * phi + (8.129e-05)) * pp * pp + ((7.2249999999999994e-06) * phi * phi + (-0.000178946) * phi + (-0.0009485399999999999)) * pp + ((-1.3018000000000003e-05) * phi * phi + (0.00046643000000000005) * phi + (0.00266508));
                }
            }
        }
        //==============================//  Electron Corrections (End)  //==============================//

        //==============================//        π+ Corrections        //==============================//
        if (corPip != 0 && ivec == 1) {
            if (corPip == 1) { // Fall 2018 - Pass 1 Corrections
                if (sec == 1) {
                    dp = ((-5.4904e-07) * phi * phi + (-1.4436e-05) * phi + (3.1534e-04)) * pp * pp + ((3.8231e-06) * phi * phi + (3.6582e-04) * phi + (-0.0046759)) * pp + ((-5.4913e-06) * phi * phi + (-4.0157e-04) * phi + (0.010767));
                    dp = dp + ((6.1103e-07) * phi * phi + (5.5291e-06) * phi + (-1.9120e-04)) * pp * pp + ((-3.2300e-06) * phi * phi + (1.5377e-05) * phi + (7.5279e-04)) * pp + ((2.1434e-06) * phi * phi + (-6.9572e-06) * phi + (-7.9333e-05));
                    dp = dp + ((-1.3049e-06) * phi * phi + (1.1295e-05) * phi + (4.5797e-04)) * pp * pp + ((9.3122e-06) * phi * phi + (-5.1074e-05) * phi + (-0.0030757)) * pp + ((-1.3102e-05) * phi * phi + (2.2153e-05) * phi + (0.0040938));
                }
                if (sec == 2) {
                    dp = ((-1.0087e-06) * phi * phi + (2.1319e-05) * phi + (7.8641e-04)) * pp * pp + ((6.7485e-06) * phi * phi + (7.3716e-05) * phi + (-0.0094591)) * pp + ((-1.1820e-05) * phi * phi + (-3.8103e-04) * phi + (0.018936));
                    dp = dp + ((8.8155e-07) * phi * phi + (-2.8257e-06) * phi + (-2.6729e-04)) * pp * pp + ((-5.4499e-06) * phi * phi + (3.8397e-05) * phi + (0.0015914)) * pp + ((6.8926e-06) * phi * phi + (-5.9386e-05) * phi + (-0.0021749));
                    dp = dp + ((-2.0147e-07) * phi * phi + (1.1061e-05) * phi + (3.8827e-04)) * pp * pp + ((4.9294e-07) * phi * phi + (-6.0257e-05) * phi + (-0.0022087)) * pp + ((9.8548e-07) * phi * phi + (5.9047e-05) * phi + (0.0022905));
                }
                if (sec == 3) {
                    dp = ((8.6722e-08) * phi * phi + (-1.7975e-05) * phi + (4.8118e-05)) * pp * pp + ((2.6273e-06) * phi * phi + (3.1453e-05) * phi + (-0.0015943)) * pp + ((-6.4463e-06) * phi * phi + (-5.8990e-05) * phi + (0.0041703));
                    dp = dp + ((9.6317e-07) * phi * phi + (-1.7659e-06) * phi + (-8.8318e-05)) * pp * pp + ((-5.1346e-06) * phi * phi + (8.3318e-06) * phi + (3.7723e-04)) * pp + ((3.9548e-06) * phi * phi + (-6.9614e-05) * phi + (2.1393e-04));
                    dp = dp + ((5.6438e-07) * phi * phi + (8.1678e-06) * phi + (-9.4406e-05)) * pp * pp + ((-3.9074e-06) * phi * phi + (-6.5174e-05) * phi + (5.4218e-04)) * pp + ((6.3198e-06) * phi * phi + (1.0611e-04) * phi + (-4.5749e-04));
                }
                if (sec == 4) {
                    dp = ((4.3406e-07) * phi * phi + (-4.9036e-06) * phi + (2.3064e-04)) * pp * pp + ((1.3624e-06) * phi * phi + (3.2907e-05) * phi + (-0.0034872)) * pp + ((-5.1017e-06) * phi * phi + (2.4593e-05) * phi + (0.0092479));
                    dp = dp + ((6.0218e-07) * phi * phi + (-1.4383e-05) * phi + (-3.1999e-05)) * pp * pp + ((-1.1243e-06) * phi * phi + (9.3884e-05) * phi + (-4.1985e-04)) * pp + ((-1.8808e-06) * phi * phi + (-1.2222e-04) * phi + (0.0014037));
                    dp = dp + ((-2.5490e-07) * phi * phi + (-8.5120e-07) * phi + (7.9109e-05)) * pp * pp + ((2.5879e-06) * phi * phi + (8.6108e-06) * phi + (-5.1533e-04)) * pp + ((-4.4521e-06) * phi * phi + (-1.7012e-05) * phi + (7.4848e-04));
                }
                if (sec == 5) {
                    dp = ((2.4292e-07) * phi * phi + (8.8741e-06) * phi + (2.9482e-04)) * pp * pp + ((3.7229e-06) * phi * phi + (7.3215e-06) * phi + (-0.0050685)) * pp + ((-1.1974e-05) * phi * phi + (-1.3043e-04) * phi + (0.0078836));
                    dp = dp + ((1.0867e-06) * phi * phi + (-7.7630e-07) * phi + (-4.4930e-05)) * pp * pp + ((-5.6564e-06) * phi * phi + (-1.3417e-05) * phi + (2.5224e-04)) * pp + ((6.8460e-06) * phi * phi + (9.0495e-05) * phi + (-4.6587e-04));
                    dp = dp + ((8.5720e-07) * phi * phi + (-6.7464e-06) * phi + (-4.0944e-05)) * pp * pp + ((-4.7370e-06) * phi * phi + (5.8808e-05) * phi + (1.9047e-04)) * pp + ((5.7404e-06) * phi * phi + (-1.1105e-04) * phi + (-1.9392e-04));
                }
                if (sec == 6) {
                    dp = ((2.1191e-06) * phi * phi + (-3.3710e-05) * phi + (2.5741e-04)) * pp * pp + ((-1.2915e-05) * phi * phi + (2.3753e-04) * phi + (-2.6882e-04)) * pp + ((2.2676e-05) * phi * phi + (-2.3115e-04) * phi + (-0.001283));
                    dp = dp + ((6.0270e-07) * phi * phi + (-6.8200e-06) * phi + (1.3103e-04)) * pp * pp + ((-1.8745e-06) * phi * phi + (3.8646e-05) * phi + (-8.8056e-04)) * pp + ((2.0885e-06) * phi * phi + (-3.4932e-05) * phi + (4.5895e-04));
                    dp = dp + ((4.7349e-08) * phi * phi + (-5.7528e-06) * phi + (-3.4097e-06)) * pp * pp + ((1.7731e-06) * phi * phi + (3.5865e-05) * phi + (-5.7881e-04)) * pp + ((-9.7008e-06) * phi * phi + (-4.1836e-05) * phi + (0.0035403));
                }
            }

            if (corPip == 2) { // Spring 2019 - Pass 2 Corrections
                if (sec == 1) {
                    dp = ((1.07338e-06) * phi * phi + (0.00011237500000000001) * phi + (0.00046984999999999996)) * pp * pp + ((-2.9323999999999997e-06) * phi * phi + (-0.000777199) * phi + (-0.0061279)) * pp + ((3.7362e-06) * phi * phi + (0.00049608) * phi + (0.0156802));
                    if (pp < 3.5) {
                        dp = dp + ((-8.0699e-06) * phi * phi + (3.3838e-04) * phi + (0.0051143)) * pp * pp + ((3.0234e-05) * phi * phi + (-0.0015167) * phi + (-0.023081)) * pp + ((-1.3818e-05) * phi * phi + (0.0011894) * phi + (0.015812));
                    } else {
                        dp = dp + ((2.8904e-07) * phi * phi + (-1.0534e-04) * phi + (-0.0023996)) * pp * pp + ((2.3276e-06) * phi * phi + (0.0010502) * phi + (0.022682)) * pp + ((-1.9319e-05) * phi * phi + (-0.0025179) * phi + (-0.050285));
                    }
                }
                if (sec == 2) {
                    dp = ((2.97335e-06) * phi * phi + (7.68257e-05) * phi + (0.001132483)) * pp * pp + ((-1.86553e-05) * phi * phi + (-0.000511963) * phi + (-0.0111051)) * pp + ((2.16081e-05) * phi * phi + (0.000100984) * phi + (0.0189673));
                    if (pp < 3.5) {
                        dp = dp + ((-1.4761e-06) * phi * phi + (4.9397e-06) * phi + (0.0014986)) * pp * pp + ((6.4311e-06) * phi * phi + (-3.8570e-05) * phi + (-0.005309)) * pp + ((2.2896e-06) * phi * phi + (-1.8426e-04) * phi + (-0.0030622));
                    } else {
                        dp = dp + ((3.3302e-06) * phi * phi + (-8.4794e-05) * phi + (-0.0020262)) * pp * pp + ((-3.5962e-05) * phi * phi + (9.1367e-04) * phi + (0.019333)) * pp + ((9.5116e-05) * phi * phi + (-0.0023371) * phi + (-0.045778));
                    }
                }
                if (sec == 3) {
                    dp = ((1.9689700000000002e-07) * phi * phi + (-6.73721e-05) * phi + (0.001145664)) * pp * pp + ((-1.3357999999999998e-07) * phi * phi + (0.0004974620000000001) * phi + (-0.01087555)) * pp + ((5.23389e-06) * phi * phi + (-0.00038631399999999996) * phi + (0.012021909999999999));
                    if (pp < 3.5) {
                        dp = dp + ((-3.7071e-06) * phi * phi + (-6.7985e-05) * phi + (0.0073195)) * pp * pp + ((1.2081e-05) * phi * phi + (4.0719e-04) * phi + (-0.032716)) * pp + ((1.8109e-06) * phi * phi + (-5.6304e-04) * phi + (0.022124));
                    } else {
                        dp = dp + ((2.9228e-06) * phi * phi + (-7.4216e-07) * phi + (-0.0033922)) * pp * pp + ((-2.7026e-05) * phi * phi + (-7.5709e-06) * phi + (0.03267)) * pp + ((5.8592e-05) * phi * phi + (3.8319e-05) * phi + (-0.076661));
                    }
                }
                if (sec == 4) {
                    dp = ((5.4899e-07) * phi * phi + (-1.82236e-05) * phi + (0.0007486388)) * pp * pp + ((-1.0743e-06) * phi * phi + (0.000125103) * phi + (-0.00743795)) * pp + ((1.9187e-06) * phi * phi + (-5.0545e-05) * phi + (0.01528271));
                    if (pp < 3.5) {
                        dp = dp + ((-7.1834e-06) * phi * phi + (1.2815e-04) * phi + (0.004323)) * pp * pp + ((2.7688e-05) * phi * phi + (-4.9122e-04) * phi + (-0.020112)) * pp + ((-1.5879e-05) * phi * phi + (3.5148e-04) * phi + (0.013367));
                    } else {
                        dp = dp + ((-2.2635e-06) * phi * phi + (3.3612e-05) * phi + (-0.0024779)) * pp * pp + ((2.7765e-05) * phi * phi + (-4.4868e-04) * phi + (0.02433)) * pp + ((-7.6567e-05) * phi * phi + (0.0013553) * phi + (-0.058136));
                    }
                }
                if (sec == 5) {
                    dp = ((9.5628e-07) * phi * phi + (-1.4e-06) * phi + (0.00116279)) * pp * pp + ((-3.723047e-06) * phi * phi + (2.09447e-05) * phi + (-0.0101853)) * pp + ((9.326299999999999e-06) * phi * phi + (-0.0001111214) * phi + (0.0130134));
                    if (pp < 3.5) {
                        dp = dp + ((-8.2807e-06) * phi * phi + (-1.2620e-04) * phi + (0.0060821)) * pp * pp + ((3.8915e-05) * phi * phi + (6.3989e-04) * phi + (-0.028784)) * pp + ((-3.7765e-05) * phi * phi + (-7.0844e-04) * phi + (0.021177));
                    } else {
                        dp = dp + ((-8.7415e-08) * phi * phi + (3.5806e-05) * phi + (-0.0022065)) * pp * pp + ((5.3612e-06) * phi * phi + (-4.2740e-04) * phi + (0.022369)) * pp + ((-2.3587e-05) * phi * phi + (0.0011096) * phi + (-0.056773));
                    }
                }
                if (sec == 6) {
                    dp = ((5.86478e-07) * phi * phi + (3.5833999999999994e-06) * phi + (0.00108574)) * pp * pp + ((-4.433118e-06) * phi * phi + (-5.3565999999999995e-05) * phi + (-0.00873827)) * pp + ((2.0270600000000002e-05) * phi * phi + (-7.0902e-05) * phi + (0.0077521));
                    if (pp < 3.5) {
                        dp = dp + ((1.4952e-06) * phi * phi + (1.3858e-05) * phi + (0.0028677)) * pp * pp + ((-8.0852e-06) * phi * phi + (-1.1384e-04) * phi + (-0.015643)) * pp + ((9.5078e-06) * phi * phi + (1.3285e-04) * phi + (0.014019));
                    } else {
                        dp = dp + ((-5.7308e-07) * phi * phi + (-3.8697e-05) * phi + (-0.0030495)) * pp * pp + ((1.0905e-05) * phi * phi + (3.8288e-04) * phi + (0.030355)) * pp + ((-3.1873e-05) * phi * phi + (-9.6019e-04) * phi + (-0.074345));
                    }
                }
            }

            if (corPip == 3) { // Fall 2018 - Pass 2 Corrections
                if (sec == 1) {
                    dp = ((1.338454e-06) * phi * phi + (4.714629999999999e-05) * phi + (0.00014719)) * pp * pp + ((-2.8460000000000004e-06) * phi * phi + (-0.000406925) * phi + (-0.00367325)) * pp + ((-1.193548e-05) * phi * phi + (-0.000225083) * phi + (0.01544091));
                    if (pp < 2.5) {
                        dp = dp + ((1.0929e-05) * phi * phi + (-3.8002e-04) * phi + (-0.01412)) * pp * pp + ((-2.8491e-05) * phi * phi + (5.0952e-04) * phi + (0.037728)) * pp + ((1.6927e-05) * phi * phi + (1.8165e-04) * phi + (-0.027772));
                    } else {
                        dp = dp + ((4.3191e-07) * phi * phi + (-9.0581e-05) * phi + (-0.0011766)) * pp * pp + ((-3.6232e-06) * phi * phi + (0.0010342) * phi + (0.012454)) * pp + ((1.2235e-05) * phi * phi + (-0.0025855) * phi + (-0.035323));
                    }
                    dp = dp + ((-3.7494e-07) * phi * phi + (-1.5439e-06) * phi + (4.2760e-05)) * pp * pp + ((3.5348e-06) * phi * phi + (4.8165e-05) * phi + (-2.3799e-04)) * pp + ((-8.2116e-06) * phi * phi + (-7.1750e-05) * phi + (1.5984e-04));
                }
                if (sec == 2) {
                    dp = ((5.8222e-07) * phi * phi + (5.0666599999999994e-05) * phi + (0.00051782)) * pp * pp + ((3.3785e-06) * phi * phi + (-0.000343093) * phi + (-0.007453400000000001)) * pp + ((-2.2014899999999998e-05) * phi * phi + (-0.00027579899999999997) * phi + (0.015119099999999998));
                    if (pp < 2.5) {
                        dp = dp + ((9.2373e-06) * phi * phi + (-3.3151e-04) * phi + (-0.019254)) * pp * pp + ((-2.7546e-05) * phi * phi + (5.3915e-04) * phi + (0.052516)) * pp + ((2.5220e-05) * phi * phi + (7.5362e-05) * phi + (-0.033504));
                    } else {
                        dp = dp + ((2.2654e-08) * phi * phi + (-8.8436e-05) * phi + (-0.0013542)) * pp * pp + ((3.0630e-07) * phi * phi + (9.4319e-04) * phi + (0.0147)) * pp + ((-3.5941e-06) * phi * phi + (-0.0022473) * phi + (-0.036874));
                    }
                    dp = dp + ((4.3694e-07) * phi * phi + (1.1476e-05) * phi + (1.1123e-04)) * pp * pp + ((-2.4617e-06) * phi * phi + (-7.5353e-05) * phi + (-6.2511e-04)) * pp + ((-1.0387e-06) * phi * phi + (5.8447e-05) * phi + (6.4986e-04));
                }
                if (sec == 3) {
                    dp = ((-6.17815e-07) * phi * phi + (-1.4503600000000001e-05) * phi + (0.000584689)) * pp * pp + ((8.27871e-06) * phi * phi + (9.2796e-05) * phi + (-0.0078185692)) * pp + ((-1.6866360000000002e-05) * phi * phi + (-8.065000000000001e-05) * phi + (0.0159476));
                    if (pp < 2.5) {
                        dp = dp + ((1.8595e-06) * phi * phi + (3.6900e-04) * phi + (-0.0099622)) * pp * pp + ((8.4410e-06) * phi * phi + (-0.0010457) * phi + (0.027038)) * pp + ((-1.2191e-05) * phi * phi + (6.0203e-04) * phi + (-0.019176));
                    } else {
                        dp = dp + ((6.8265e-07) * phi * phi + (3.0246e-05) * phi + (-0.0011116)) * pp * pp + ((-4.8481e-06) * phi * phi + (-3.7082e-04) * phi + (0.011452)) * pp + ((7.2478e-06) * phi * phi + (9.9858e-04) * phi + (-0.027972));
                    }
                    dp = dp + ((1.8639e-07) * phi * phi + (4.9444e-06) * phi + (-2.9030e-05)) * pp * pp + ((-1.3752e-06) * phi * phi + (-3.3709e-05) * phi + (3.8288e-04)) * pp + ((1.0113e-06) * phi * phi + (5.1273e-05) * phi + (-6.7844e-04));
                }
                if (sec == 4) {
                    dp = ((9.379499999999998e-07) * phi * phi + (-2.8101700000000002e-05) * phi + (0.00053373)) * pp * pp + ((-1.6185199999999991e-06) * phi * phi + (0.00017444500000000001) * phi + (-0.005648269999999999)) * pp + ((-3.495700000000003e-06) * phi * phi + (-7.845739999999999e-05) * phi + (0.010768400000000001));
                    if (pp < 2.5) {
                        dp = dp + ((9.5779e-06) * phi * phi + (3.5339e-04) * phi + (-0.01054)) * pp * pp + ((-1.8077e-05) * phi * phi + (-0.0010543) * phi + (0.028379)) * pp + ((3.1773e-06) * phi * phi + (5.6223e-04) * phi + (-0.018865));
                    } else {
                        dp = dp + ((7.7000e-07) * phi * phi + (4.1000e-06) * phi + (-0.0010144)) * pp * pp + ((-8.1960e-06) * phi * phi + (-4.7753e-05) * phi + (0.010594)) * pp + ((2.0716e-05) * phi * phi + (1.2151e-04) * phi + (-0.028619));
                    }
                    dp = dp + ((4.8394e-07) * phi * phi + (3.6342e-06) * phi + (-2.0136e-04)) * pp * pp + ((-3.2757e-06) * phi * phi + (-3.5397e-05) * phi + (0.0015599)) * pp + ((3.2095e-06) * phi * phi + (7.9013e-05) * phi + (-0.002012));
                }
                if (sec == 5) {
                    dp = ((1.7566900000000006e-07) * phi * phi + (2.21337e-05) * phi + (0.0011632)) * pp * pp + ((2.812770000000001e-06) * phi * phi + (-0.00018654499999999998) * phi + (-0.011854620000000001)) * pp + ((-8.442900000000003e-06) * phi * phi + (-0.00011505800000000001) * phi + (0.0176174));
                    if (pp < 2.5) {
                        dp = dp + ((3.3685e-05) * phi * phi + (2.8972e-04) * phi + (-0.017862)) * pp * pp + ((-8.4089e-05) * phi * phi + (-9.8038e-04) * phi + (0.050405)) * pp + ((4.3478e-05) * phi * phi + (6.9924e-04) * phi + (-0.033066));
                    } else {
                        dp = dp + ((4.6106e-07) * phi * phi + (-3.6786e-05) * phi + (-0.0015894)) * pp * pp + ((-4.4217e-06) * phi * phi + (3.7321e-04) * phi + (0.015917)) * pp + ((7.5188e-06) * phi * phi + (-8.0676e-04) * phi + (-0.036944));
                    }
                    dp = dp + ((4.3113e-07) * phi * phi + (2.6869e-06) * phi + (-2.1326e-04)) * pp * pp + ((-3.1063e-06) * phi * phi + (-2.7152e-05) * phi + (0.0017964)) * pp + ((3.1946e-06) * phi * phi + (4.2059e-05) * phi + (-0.0031325));
                }
                if (sec == 6) {
                    dp = ((1.94354e-06) * phi * phi + (1.3306000000000006e-05) * phi + (0.00067634)) * pp * pp + ((-7.9584e-06) * phi * phi + (-7.949999999999998e-05) * phi + (-0.005861990000000001)) * pp + ((6.994000000000005e-07) * phi * phi + (-0.00022435) * phi + (0.0118564));
                    if (pp < 2.5) {
                        dp = dp + ((1.7381e-05) * phi * phi + (5.4630e-04) * phi + (-0.019637)) * pp * pp + ((-3.8681e-05) * phi * phi + (-0.0017358) * phi + (0.0565)) * pp + ((1.2268e-05) * phi * phi + (0.0011412) * phi + (-0.035608));
                    } else {
                        dp = dp + ((-8.9398e-08) * phi * phi + (-1.2347e-05) * phi + (-0.0018442)) * pp * pp + ((7.8164e-08) * phi * phi + (1.3063e-04) * phi + (0.01783)) * pp + ((8.2374e-06) * phi * phi + (-3.5862e-04) * phi + (-0.047011));
                    }
                    dp = dp + ((4.9123e-07) * phi * phi + (5.1828e-06) * phi + (-1.3898e-04)) * pp * pp + ((-3.4108e-06) * phi * phi + (-5.0009e-05) * phi + (0.0014879)) * pp + ((4.0320e-06) * phi * phi + (6.5853e-05) * phi + (-0.0032227));
                }

            }
        }
        //==============================//     π+ Corrections (End)     //==============================//

        //==============================//        π- Corrections        //==============================//
        if (corPim != 0 && ivec == 2) {
            if (sec == 1) { // Fall 2018 - Pass 1 Corrections (Only)
                dp = ((-4.0192658422317425e-06) * phi * phi - (2.660222128967742e-05) * phi + 0.004774434682983547) * pp * pp;
                dp = dp + ((1.9549520962477972e-05) * phi * phi - 0.0002456062756770577 * phi - 0.03787692408323466) * pp;
                dp = dp + (-2.128953094937459e-05) * phi * phi + 0.0002461708852239913 * phi + 0.08060704449822174 - 0.01;
            }
            if (sec == 2) {
                dp = ((1.193010521758372e-05) * phi * phi - (5.996221756031922e-05) * phi + 0.0009093437955814359) * pp * pp;
                dp = dp + ((-4.89113824430594e-05) * phi * phi + 0.00021676479488147118 * phi - 0.01861892053916726) * pp;
                dp = dp + (4.446394152208071e-05) * phi * phi - (3.6592784167335244e-05) * phi + 0.05498710249944096 - 0.01;
            }
            if (sec == 3) {
                dp = ((-1.6596664895992133e-07) * phi * phi + (6.317189710683516e-05) * phi + 0.0016364212312654086) * pp * pp;
                dp = dp + ((-2.898409777520318e-07) * phi * phi - 0.00014531513577533802 * phi - 0.025456145839203827) * pp;
                dp = dp + (2.6432552410603506e-06) * phi * phi + 0.00018447151306275443 * phi + 0.06442602664627255 - 0.01;
            }
            if (sec == 4) {
                dp = ((2.4035259647558634e-07) * phi * phi - (8.649647351491232e-06) * phi + 0.004558993439848128) * pp * pp;
                dp = dp + ((-5.981498144060984e-06) * phi * phi + 0.00010582131454222416 * phi - 0.033572004651981686) * pp;
                dp = dp + (8.70140266889548e-06) * phi * phi - 0.00020137414379966883 * phi + 0.07258774523336173 - 0.01;
            }
            if (sec == 5) {
                dp = ((2.5817024702834863e-06) * phi * phi + 0.00010132810066914441 * phi + 0.003397314538804711) * pp * pp;
                dp = dp + ((-1.5116941263931812e-05) * phi * phi - 0.00040679799541839254 * phi - 0.028144285760769876) * pp;
                dp = dp + (1.4701931057951464e-05) * phi * phi + 0.0002426350390593454 * phi + 0.06781682510174941 - 0.01;
            }
            if (sec == 6) {
                dp = ((-8.196823669099362e-07) * phi * phi - (5.280412421933636e-05) * phi + 0.0018457238328451137) * pp * pp;
                dp = dp + ((5.2675062282094536e-06) * phi * phi + 0.0001515803461044587 * phi - 0.02294371578470564) * pp;
                dp = dp + (-9.459454671739747e-06) * phi * phi - 0.0002389523716779765 * phi + 0.06428970810739926 - 0.01;
            }
        }
        //==============================//     π- Corrections (End)     //==============================//

        if (pp != 0) {
            double scale = (pp + dp) / pp;
            p_array[0] = (float) (Px * scale);
            p_array[1] = (float) (Py * scale);
            p_array[2] = (float) (Pz * scale);
        }
    }

    public void outbending_momentum_corrections(
            float[] p_array,
            int sec,
            int ivec,
            int corEl,
            int corPip,
            int corPim,
            int corPro) {
        /**
         * Apply outbending torus momentum corrections to a particle's momentum.
         *
         * @param p_array float[3] holding {Px, Py, Pz}; updated in place
         * @param sec Forward detector sector (1–6)
         * @param ivec particle index: 0=electron, 1=π⁺, 2=π⁻, 3=proton
         * @param corEl electron correction flag (0=off, 1=Fall18-P1, 2=Fall18-P2)
         * @param corPip π⁺ correction flag (0=off, 1=Fall18-P1, 2=Fall18-P2)
         * @param corPim π⁻ correction flag (0=off, 1=Fall18-P1)
         * @param corPro proton correction flag (0=off; outbending not available)
         */
        // Cartesian components
        double Px = p_array[0];
        double Py = p_array[1];
        double Pz = p_array[2];
        // Momentum magnitude
        double pp = Math.sqrt(Px * Px + Py * Py + Pz * Pz);
        // Correction factor
        double dp = 0;

        // Phi angle in degrees
        double Phi = (180.0 / 3.1415926) * Math.atan2(Py, Px);
        // Sector realignment
        if (((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)) {
            Phi += 360;
        }
        // Local phi
        double phi = Phi - (sec - 1) * 60;
        // Shifts
        if (ivec == 0) {
            phi -= 30 / pp;
        } else if (ivec == 1 || ivec == 3) {
            phi += 32 / (pp - 0.05);
        } else if (ivec == 2) {
            phi -= 32 / (pp - 0.05);
        }

        // Electron corrections (outbending)
        if (ivec == 0 && corEl != 0) {
            if (corEl == 1) { // Fall 2018 - Pass 1
                switch (sec) {
                    case 1:
                        dp = ((1.3189e-06) * phi * phi + (4.26057e-05) * phi + (-0.002322628)) * pp * pp
                                + ((-1.1409e-05) * phi * phi + (2.2188e-05) * phi + (0.02878927)) * pp
                                + ((2.4950e-05) * phi * phi + (1.6170e-06) * phi + (-0.061816275));
                        break;
                    case 2:
                        dp = ((-2.9240e-07) * phi * phi + (3.2448e-07) * phi + (-0.001848308)) * pp * pp
                                + ((4.4500e-07) * phi * phi + (4.76324e-04) * phi + (0.02219469)) * pp
                                + ((6.9220e-06) * phi * phi + (-0.00153517) * phi + (-0.0479058));
                        break;
                    case 3:
                        dp = ((2.71911e-06) * phi * phi + (1.657148e-05) * phi + (-0.001822211)) * pp * pp
                                + ((-4.96814e-05) * phi * phi + (-3.761117e-04) * phi + (0.02564148)) * pp
                                + ((1.97748e-04) * phi * phi + (9.58259e-04) * phi + (-0.05818292));
                        break;
                    case 4:
                        dp = ((1.90966e-06) * phi * phi + (-2.4761e-05) * phi + (-0.00231562)) * pp * pp
                                + ((-2.3927e-05) * phi * phi + (2.25262e-04) * phi + (0.0291831)) * pp
                                + ((8.0515e-05) * phi * phi + (-6.42098e-04) * phi + (-0.06159197));
                        break;
                    case 5:
                        dp = ((-3.6760323e-06) * phi * phi + (4.04398e-05) * phi + (-0.0021967515)) * pp * pp
                                + ((4.90857e-05) * phi * phi + (-4.37437e-04) * phi + (0.02494339)) * pp
                                + ((-1.08257e-04) * phi * phi + (0.00146111) * phi + (-0.0648485));
                        break;
                    case 6:
                        dp = ((-6.2488e-08) * phi * phi + (2.23173e-05) * phi + (-0.00227522)) * pp * pp
                                + ((1.8372e-05) * phi * phi + (-7.5227e-05) * phi + (0.032636)) * pp
                                + ((-6.6566e-05) * phi * phi + (-2.4450e-04) * phi + (-0.072293));
                        break;
                }
            }
            if (corEl == 2) { // Fall 2018 - Pass 2
                switch (sec) {
                    case 1:
                        dp = ((-5.74868e-06) * phi * phi + (2.2331e-06) * phi + (0.00025487)) * pp * pp
                                + ((8.18043e-05) * phi * phi + (9.8383e-05) * phi + (-0.0056062)) * pp
                                + ((-0.0002393096) * phi * phi + (-0.001175548) * phi + (0.0486792));
                        if (pp < 6.5) {
                            dp += ((-1.18269e-05) * phi * phi + (-3.33e-05) * phi + (0.0002424)) * pp * pp
                                    + ((8.541e-05) * phi * phi + (0.0002206) * phi + (0.000609)) * pp
                                    + ((-0.00012801) * phi * phi + (-0.0002079) * phi + (-0.007417));
                        } else {
                            dp += ((-1.1703e-06) * phi * phi + (-9.8669e-05) * phi + (-0.002884)) * pp * pp
                                    + ((2.591e-05) * phi * phi + (0.00169554) * phi + (0.045387)) * pp
                                    + ((-0.00012285) * phi * phi + (-0.00710292) * phi + (-0.17897));
                        }
                        break;
                    case 2:
                        dp = ((-4.29711e-06) * phi * phi + (1.36105e-05) * phi + (0.00040706)) * pp * pp
                                + ((6.264498e-05) * phi * phi + (9.8981e-05) * phi + (-0.0059589)) * pp
                                + ((-0.000175865) * phi * phi + (-0.001329597) * phi + (0.030682));
                        if (pp < 6.5) {
                            dp += (1.42031e-05) * phi * phi + (-4.3073e-05) * phi + (-0.0013709);
                            dp *= pp;
                        } else {
                            dp += (-7.3709e-06) * phi * phi + (3.0237e-05) * phi + (-0.00556253);
                            dp *= pp;
                        }
                        break;
                    case 3:
                        dp = ((-3.54616e-06) * phi * phi + (2.1382e-05) * phi + (-0.00029815)) * pp * pp
                                + ((4.09174e-05) * phi * phi + (-0.000258372) * phi + (0.0028638)) * pp
                                + ((-9.1374e-05) * phi * phi + (0.00082102) * phi + (0.00818));
                        if (pp < 4.75) {
                            dp += ((-4.1035e-05) * phi * phi + (-0.00034997) * phi + (0.0037639)) * pp * pp
                                    + ((0.00027659) * phi * phi + (0.0025883) * phi + (-0.027382)) * pp
                                    + ((-0.0004434) * phi * phi + (-0.0047757) * phi + (0.049737));
                        } else {
                            dp += (8.0625e-06) * phi * phi + (-1.5445e-05) * phi + (-0.00213057);
                            dp *= pp * pp;
                        }
                        break;
                    case 4:
                        dp = ((2.43038e-06) * phi * phi + (-6.4715e-06) * phi + (-0.00076198)) * pp * pp
                                + ((-2.606703e-05) * phi * phi + (7.1258e-05) * phi + (0.00984946)) * pp
                                + ((6.69394e-05) * phi * phi + (-0.000264164) * phi + (-0.0021348));
                        if (pp < 6.5) {
                            dp += (4.3553e-06) * phi * phi + (-5.2965e-05) * phi + (0.0005391);
                            dp *= pp * pp;
                        } else {
                            dp += (-1.6976e-05) * phi * phi + (0.00011328) * phi + (-0.0032102);
                            dp *= pp * pp;
                        }
                        break;
                    case 5:
                        dp = ((-1.7474881e-06) * phi * phi + (5.20199e-05) * phi + (0.00056946)) * pp * pp
                                + ((2.579126e-05) * phi * phi + (-0.00049884) * phi + (-0.00743872)) * pp
                                + ((-5.9489e-05) * phi * phi + (0.00070898) * phi + (0.0312296));
                        if (pp < 4.75) {
                            dp += (-9.6494e-05) * phi * phi + (-0.00051529) * phi + (0.0157187);
                            dp *= pp * pp;
                        } else {
                            dp += (1.3183e-06) * phi * phi + (5.7112e-05) * phi + (-0.002674);
                            dp *= pp * pp;
                        }
                        break;
                    case 6:
                        dp = ((-7.1587e-07) * phi * phi + (3.31516e-05) * phi + (-8.785999e-05)) * pp * pp
                                + ((1.21201e-05) * phi * phi + (-0.000229416) * phi + (0.001413)) * pp
                                + ((1.222e-05) * phi * phi + (0.000154449) * phi + (0.011577));
                        if (pp < 4.75) {
                            dp += (5.3238e-05) * phi * phi + (-0.000348957) * phi + (-0.0047347);
                            dp *= pp * pp;
                        } else {
                            dp += (-2.6885e-06) * phi * phi + (1.43847e-05) * phi + (-0.00193047);
                            dp *= pp * pp;
                        }
                        break;
                }
            }
        }

        // π- corrections (Fall 2018 - Pass 1 only)
        if (ivec == 2 && corPim != 0) {
            switch (sec) {
                case 1:
                    dp = ((2.7123584594392597e-06) * phi * phi + (-5.468601175954242e-05) * phi + (0.002313330256974031)) * pp * pp
                            + ((-8.039703360516874e-06) * phi * phi + (0.00044464879674067275) * phi + (-0.02546911446157775)) * pp
                            + ((3.5973669277966655e-06) * phi * phi + (-0.0003856844699023182) * phi + (0.05496480659602064) - 0.015);
                    break;
                case 2:
                    dp = ((1.9081500905303347e-06) * phi * phi + (3.310647986349362e-05) * phi + (-0.0003264357817968204)) * pp * pp
                            + ((-1.2306311457915714e-05) * phi * phi + (-6.404982516446639e-05) * phi + (-0.01287404671840319)) * pp
                            + ((9.746651642120768e-06) * phi * phi + (6.1503461629194e-05) * phi + (0.04249861359511857) - 0.015);
                    break;
                case 3:
                    dp = ((3.467960715633796e-06) * phi * phi + (-0.00011427345789836184) * phi + (0.004780571116355615)) * pp * pp
                            + ((-1.2639455891842017e-05) * phi * phi + (0.00044737258600913664) * phi + (-0.03827009444373719)) * pp
                            + ((5.8243648992776484e-06) * phi * phi + (-0.0004240381542174731) * phi + (0.06589846610477122) - 0.015);
                    break;
                case 4:
                    dp = ((-7.97757466039691e-06) * phi * phi + (-0.00011075801628158914) * phi + (0.006505144041475733)) * pp * pp
                            + ((3.570788801587046e-05) * phi * phi + (0.0005835525352273808) * phi + (-0.045031773715754606)) * pp
                            + ((-3.223327114068019e-05) * phi * phi + (-0.0006144362450858762) * phi + (0.07280937684254037) - 0.015);
                    break;
                case 5:
                    dp = ((1.990802625607816e-06) * phi * phi + (7.057771450607931e-05) * phi + (0.005399025205722829)) * pp * pp
                            + ((-7.670376562908147e-06) * phi * phi + (-0.00032508260930191955) * phi + (-0.044439500813069875)) * pp
                            + ((7.599354976329091e-06) * phi * phi + (0.0002562152836894338) * phi + (0.07195292224032898) - 0.015);
                    break;
                case 6:
                    dp = ((1.9247834787602347e-06) * phi * phi + (7.638857332736951e-05) * phi + (0.005271258583881754)) * pp * pp
                            + ((-2.7349724034956845e-06) * phi * phi + (-0.00016130256163798413) * phi + (-0.03668300882287307)) * pp
                            + ((7.40942843287096e-07) * phi * phi + (-5.785254680184232e-05) * phi + (0.06282320712979896) - 0.015);
                    break;
            }
        }

        // Apply correction
        if (pp != 0) {
            double scale = (pp + dp) / pp;
            p_array[0] = (float) (Px * scale);
            p_array[1] = (float) (Py * scale);
            p_array[2] = (float) (Pz * scale);
        }
    }

    public void jpsi_momentum_corrections(int particle_Index, float[] p_array,
            HipoDataBank rec_Bank, HipoDataBank run_Bank, HipoDataBank track_Bank) {

        double px = p_array[0];
        double py = p_array[1];
        double pz = p_array[2];

        double p = p_calculation(px, py, pz);
        double theta = theta_calculation(px, py, pz);
        double phi = phi_calculation(px, py);

        boolean inbending = run_Bank.getFloat("torus", 0) == -1, outbending = !inbending;

        double dp = 0;
        if (inbending) {
            dp = 0.0093796 * p + (-0.0007808) * p * p + (0.000163) * p * p * p + (-0.02029) + 0.04121 / p;
        } else {
            dp = (-0.06520) * p + 0.007099 * p * p + (-0.00005929) * p * p * p + 0.2145 + -(0.1153) / p;
        }

        p = +dp;
        // Update the px, py, pz values
        p_array[0] = (float) x_calculation(p, theta, phi);
        p_array[1] = (float) y_calculation(p, theta, phi);
        p_array[2] = (float) z_calculation(p, theta);

    }

    //––––––––––– Fermi-motion tables for N14 –––––––––––
    private static final double[] kfm_fm = {
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
        2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
        3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
        4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
        5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9,
        6.0
    };

    private static final double[] nk_n14 = {
        0.163119, 0.183235, 0.225102, 0.262783, 0.276924, 0.263379, 0.228534, 0.183451, 0.137812, 0.097877,
        0.066275, 0.043066, 0.026994, 0.016454, 0.009841, 0.005874, 0.003586, 0.002312, 0.001622, 0.001247,
        0.001036, 0.000903, 0.000805, 0.000724, 0.000652, 0.000585, 0.000523, 0.000465, 0.000412, 0.000364,
        0.000320, 0.000281, 0.000246, 0.000216, 0.000188, 0.000164, 0.000143, 0.000124, 0.000107, 0.000092,
        0.000079, 0.000067, 0.000057, 0.000048, 0.000040, 0.000034, 0.000029, 0.000024, 0.000020, 0.000017,
        0.000014, 0.000012, 0.000010, 0.000008, 0.000006, 0.000005, 0.000005, 0.000004, 0.000003, 0.000003,
        0.000002
    };

    private static final double HBARC = 0.19732;  // GeV·fm
    private static final Random rand = new Random();
    private static double[] cdfN;             // cumulative 4π k² n(k) dk
    private static boolean fermiInit = false;
    private static double fermiScale = 1.0;

//––– fit parameters for D_f(Mx2)=poly₄+Gauss from python result
    private static final double a0 = -0.0236548;
    private static final double a1 = 0.355323;
    private static final double a2 = -0.210453;
    private static final double a3 = 0.0477062;
    private static final double a4 = -0.00254719;
    private static final double GA = 0.189287;
    private static final double GM = 0.889894;
    private static final double GS = 0.0794461;

    public static void setFermiScale(double scale) {
        fermiScale = scale;
    }

    /**
     * Dilution factor D_f(Mx2) = a0 + a1 x + ... + a4 x⁴ + GA·exp[-(x–GM)²/(2·GS²)]
     */
    private static double dilutionFactor(double mx2) {
        double poly = a0 + a1 * mx2 + a2 * mx2 * mx2 + a3 * mx2 * mx2 * mx2 + a4 * Math.pow(mx2, 4);
        double gauss = GA * Math.exp(-(mx2 - GM) * (mx2 - GM) / (2 * GS * GS));
        return poly + gauss;
    }

    private static void initFermi() {
        if (fermiInit) {
            return;
        }
        int N = kfm_fm.length;
        double[] F = new double[N];
        for (int i = 0; i < N; i++) {
            F[i] = 4 * Math.PI * kfm_fm[i] * kfm_fm[i] * nk_n14[i];
        }
        cdfN = new double[N];
        cdfN[0] = 0;
        double integral = 0;
        for (int i = 0; i < N - 1; i++) {
            double dk = kfm_fm[i + 1] - kfm_fm[i];
            integral += 0.5 * (F[i] + F[i + 1]) * dk;
            cdfN[i + 1] = integral;
        }
        for (int i = 0; i < N; i++) {
            cdfN[i] /= integral;
        }
        fermiInit = true;
    }

    /**
     * @param mx2 the event’s missing‐mass squared (GeV²)
     * @return A Vector3(px,py,pz) drawn from the N14 Fermi distribution (GeV/c), with probability of smearing = 1 –
     * D_f(mx2); otherwise (0,0,0).
     */
    public static Vector3 sampleFermiMomentum(double mx2) {
        initFermi();
        // smear with probability = 1 – D_f(mx2)
        double pSmear = 1.0 - dilutionFactor(mx2);
        System.out.println(dilutionFactor(mx2));
        if (rand.nextDouble() < pSmear) {
            // leave at rest
            return new Vector3(0.0, 0.0, 0.0);
        }
        // now sample k from CDF
        double r = rand.nextDouble();
        int idx = Arrays.binarySearch(cdfN, r);
        if (idx < 0) {
            idx = -idx - 2;
        }
        idx = Math.max(0, Math.min(idx, cdfN.length - 2));
        double k1 = kfm_fm[idx], k2 = kfm_fm[idx + 1];
        double c1 = cdfN[idx], c2 = cdfN[idx + 1];
        double kfm = k1 + (r - c1) * (k2 - k1) / (c2 - c1);
        // convert → |p| and apply scale
        double pMag = kfm * HBARC * fermiScale;
        // random isotropic direction
        double cosT = 2 * rand.nextDouble() - 1;
        double sinT = Math.sqrt(1 - cosT * cosT);
        double phi = 2 * Math.PI * rand.nextDouble();
        double px = pMag * sinT * Math.cos(phi);
        double py = pMag * sinT * Math.sin(phi);
        double pz = pMag * cosT;
        return new Vector3(px, py, pz);
    }

}
