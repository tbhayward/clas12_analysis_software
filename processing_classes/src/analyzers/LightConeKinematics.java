
package analyzers;

import org.jlab.clas.physics.*;

/**
 *
 * @author tbhayward
 */
// Convert Lorentz vectors to light-cone coordinates
public class LightConeKinematics {

    // Function to calculate a^+ component
    public double aPlus(LorentzVector lv) {
        return (lv.e() + lv.pz()) / Math.sqrt(2);
    }

    // Function to calculate a^- component
    public double aMinus(LorentzVector lv) {
        return (lv.e() - lv.pz()) / Math.sqrt(2);
    }

    // Function to return the transverse components (a_⊥)
    public double[] aTransverse(LorentzVector lv) {
        return new double[]{lv.px(), lv.py()};
    }

    // Light-cone inner product for transverse components
    public double transverseInnerProduct(LorentzVector lv1, LorentzVector lv2) {
        return lv1.px() * lv2.px() + lv1.py() * lv2.py();
    }

    // Full inner product for light-cone coordinates
    public double lightConeInnerProduct(LorentzVector lv1, LorentzVector lv2) {
        double plusProduct = aPlus(lv1) * aMinus(lv2) + aMinus(lv1) * aPlus(lv2);
        double transverseProduct = transverseInnerProduct(lv1, lv2);
        return plusProduct - transverseProduct;
    }

    // Calculate xi_h in light-cone coordinates
    public double xi_h(LorentzVector lv_p_gN, LorentzVector lv_q_gN, LorentzVector lv_target_gN) {
        double numerator = lightConeInnerProduct(lv_p_gN, lv_q_gN);
        double denominator = lightConeInnerProduct(lv_target_gN, lv_q_gN);
        return numerator / denominator;
    }

}
