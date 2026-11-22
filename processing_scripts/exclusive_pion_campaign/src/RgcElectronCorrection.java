package exclusivepion;
/** RGC pi+ note, channel_selection.tex, electron correction tables.
 * p_corrected = (1 + delta(theta_degrees)) p_raw; theta clamped to [8,15].
 * Applied once, after detector/PID selection and before all kinematics.
 */
public final class RgcElectronCorrection {
 private static final double[][][] C={
  {{-.003276,3.9818e-4,0},{.004490,-1.4350e-4,0},{.004677,-4.2114e-4,0},{.002444,-9.3691e-5,0},{.001267,-2.2314e-4,0},{-.001890,6.3591e-5,0}},
  {{1.9971e-4,1.0813e-4,0},{.004853,-1.1926e-4,0},{.009825,-6.3016e-4,0},{.003173,2.0016e-4,0},{.007813,-6.7382e-4,0},{.005404,-4.6323e-4,0}},
  {{-9.8063e-4,1.0435e-4,0},{.016937,-.002505,1.0245e-4},{.008040,-6.4998e-4,0},{-.010207,.002954,-1.3456e-4},{.005466,-4.2851e-4,0},{.020615,-.003468,1.3870e-4}},
  {{-.017000,.003460,-1.6146e-4},{.009163,-4.5048e-4,0},{.014589,-.001394,0},{.004071,-3.5391e-4,0},{.003038,-5.0218e-4,0},{.003121,-3.6366e-4,0}}
 };
 public static double scale(int run,int sector,double thetaDegrees,double solenoid){
  int period;
  if(run>=16043&&run<=16772)period=0;
  else if(run>=16843&&run<=17408){
   if(!Double.isFinite(solenoid)||solenoid==0)throw new IllegalArgumentException("Missing Fa22 solenoid polarity");
   period=solenoid<0?1:2;
  }else if(run>=17477&&run<=17811)period=3;
  else if(run<16000)return 1; // RGB: no RGC calibration applied.
  else throw new IllegalArgumentException("Unmapped RGC electron correction run "+run);
  if(sector<1||sector>6||!Double.isFinite(thetaDegrees))throw new IllegalArgumentException("Invalid electron sector/theta");
  double t=Math.max(8,Math.min(15,thetaDegrees));double[] c=C[period][sector-1];return 1+c[0]+t*(c[1]+t*c[2]);
 }
}
