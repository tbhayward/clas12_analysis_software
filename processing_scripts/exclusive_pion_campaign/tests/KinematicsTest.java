package exclusivepion;
import static exclusivepion.ExclusivePionVariables.*;
public class KinematicsTest {
 static void near(double a,double b){if(Math.abs(a-b)>1e-9)throw new AssertionError(a+" != "+b);}
 public static void main(String[] args){
  double eb=10.55;V beam=V.track(0,0,Math.sqrt(eb*eb-ME*ME),ME),e=V.track(1.5,0,5,ME),q=beam.sub(e),total=q.add(new V(MN,0,0,0));double W=Math.sqrt(total.m2());
  double[] boost={-total.x/total.e,-total.y/total.e,-total.z/total.e};V qc=q.boost(boost),ec=e.boost(boost);double[] z=unit(qc.xyz()),normal=unit(cross(qc.xyz(),ec.xyz())),xx=cross(normal,z);
  double ep=(W*W+MPI*MPI-MP*MP)/(2*W),pp=Math.sqrt(ep*ep-MPI*MPI);
  for(double theta:new double[]{.01,.4,1.2})for(double az:new double[]{.3,2.,4.5,6.}){
   double[] v=new double[3];for(int i=0;i<3;i++)v[i]=pp*(Math.cos(theta)*z[i]+Math.sin(theta)*(Math.cos(az)*xx[i]+Math.sin(az)*normal[i]));
   V pic=new V(ep,v[0],v[1],v[2]),pc=V.track(-v[0],-v[1],-v[2],MP);double[] inv={-boost[0],-boost[1],-boost[2]};V pi=pic.boost(inv),p=pc.boost(inv);
   double[] row=calculate(new double[]{16137,123,-1,eb,1,2,2,1,1,0,2,3},e,p,pi,0,0,0,0,0,1,1);
   if(row==null)throw new AssertionError("Unexpected cut");String[] names=INPUT_SCHEMA.split(" ");if(row.length!=names.length)throw new AssertionError("Schema size");java.util.Map<String,Double> m=new java.util.HashMap<>();for(int i=0;i<row.length;i++)m.put(names[i],row[i]);
   near(m.get("Mx2"),0);near(m.get("Emiss"),0);near(m.get("pTmiss"),0);near(m.get("Mx2_1"),MPI*MPI);near(m.get("Mx2_2"),MP*MP);near(m.get("t"),q.sub(pi).m2());near(m.get("phi2"),az);near(m.get("tprime"),2*qc.p()*pp*(Math.cos(theta)-1));near(m.get("pT2"),pp*Math.sin(theta));near(m.get("xF2"),2*pp*Math.cos(theta)/W);
  }
  V forward=new V(ep,pp*z[0],pp*z[1],pp*z[2]);near(qc.sub(forward).m2(),forwardT(-q.m2(),W));
  near(helicity(6420,1),-1);near(helicity(11290,1),-1);near(helicity(11323,1),1);near(helicity(16137,1),1);
  near(RgcElectronCorrection.scale(16137,1,10,-1),1-.003276+3.9818e-4*10);
  near(RgcElectronCorrection.scale(16900,2,10,-1),1+.004853-1.1926e-4*10);
  near(RgcElectronCorrection.scale(17200,2,10,1),1+.016937-.002505*10+1.0245e-4*100);
  near(RgcElectronCorrection.scale(17500,1,10,-1),1-.017+.003460*10-1.6146e-4*100);
  for(int r:new int[]{16137,16900,17200,17500})for(int sector=1;sector<=6;sector++){
   double sol=r==17200?1:-1;near(RgcElectronCorrection.scale(r,sector,1,sol),RgcElectronCorrection.scale(r,sector,8,sol));near(RgcElectronCorrection.scale(r,sector,30,sol),RgcElectronCorrection.scale(r,sector,15,sol));
  }
  near(RgcElectronCorrection.scale(11323,0,10,-1),1);
  System.out.println("Electron corrections: four periods, six sectors, theta clamping and RGB exclusion passed");
  System.out.println("Kinematics: 12 physical events, phi signs, t/tmin/tprime, conservation, and schema passed");
 }
}
