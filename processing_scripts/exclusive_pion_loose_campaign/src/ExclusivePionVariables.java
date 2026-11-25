package exclusivepion;
import java.util.*;
/** GeV, radians. Stationary neutron reference; signed tprime = t-tmin. */
public final class ExclusivePionVariables {
 public static final double MN=.9395654133, MP=.9382720813, MPI=.13957039, ME=.00051099895;
 public static final String INPUT_SCHEMA="runnum evnum helicity beam_energy p1_index p2_index detector1 detector2 sector_e sector1 sector2 fiducial_status e_p e_theta e_phi vz_e p1_p p1_theta p1_phi vz_p1 p2_p p2_theta p2_phi vz_p2 chi2pid1 chi2pid2 beta1 beta2 Q2 W x y nu t tmin tprime Mx2 Mx2_1 Mx2_2 Emiss pTmiss Mh z2 xF2 pT2 phi1 phi2 Delta_phi open_angle_ep1 open_angle_ep2 DepA DepB DepC DepV DepW e_index detector_e status_e status1 status2 e_raw_px e_raw_py e_raw_pz p1_raw_px p1_raw_py p1_raw_pz p2_raw_px p2_raw_py p2_raw_pz chi2pid_e beta_e pass_e_vertex pass_e_sampling pass_e_diagonal pass_e_pcal1 pass_e_pcal2 pass_e_dc pass_p1_momentum pass_p1_theta pass_p1_vertex pass_p1_dc pass_p1_cvt1 pass_p1_cvt2 pass_p1_chi2pid pass_p2_pid pass_p2_vertex pass_p2_dc pass_p2_cvt1 pass_p2_cvt2 pass_e_legacy pass_p1_legacy pass_p2_legacy electron_correction proton_correction pass_dis pass_q2_1p5 kinematics_valid";
 public static final Set<String> INTS=new HashSet<>(Arrays.asList("runnum evnum helicity p1_index p2_index detector1 detector2 sector_e sector1 sector2 fiducial_status e_index detector_e status_e status1 status2 pass_e_vertex pass_e_sampling pass_e_diagonal pass_e_pcal1 pass_e_pcal2 pass_e_dc pass_p1_momentum pass_p1_theta pass_p1_vertex pass_p1_dc pass_p1_cvt1 pass_p1_cvt2 pass_p1_chi2pid pass_p2_pid pass_p2_vertex pass_p2_dc pass_p2_cvt1 pass_p2_cvt2 pass_e_legacy pass_p1_legacy pass_p2_legacy electron_correction proton_correction pass_dis pass_q2_1p5 kinematics_valid".split(" ")));
 public static class V {
  public final double e,x,y,z;
  public V(double e,double x,double y,double z){this.e=e;this.x=x;this.y=y;this.z=z;}
  public static V track(double x,double y,double z,double m){return new V(Math.sqrt(x*x+y*y+z*z+m*m),x,y,z);}
  public V add(V v){return new V(e+v.e,x+v.x,y+v.y,z+v.z);}
  public V sub(V v){return new V(e-v.e,x-v.x,y-v.y,z-v.z);}
  public double m2(){return e*e-x*x-y*y-z*z;}
  public double[] xyz(){return new double[]{x,y,z};}
  public double p(){return Math.sqrt(dot(xyz(),xyz()));}
  public V boost(double[] b){double b2=dot(b,b);if(b2==0)return this;double g=1/Math.sqrt(1-b2),bp=dot(b,xyz()),f=(g-1)*bp/b2+g*e;return new V(g*(e+bp),x+f*b[0],y+f*b[1],z+f*b[2]);}
 }
 static double dot(double[] a,double[] b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
 static double[] cross(double[] a,double[] b){return new double[]{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
 static double[] unit(double[] a){double p=Math.sqrt(dot(a,a));if(!(p>0))return new double[]{Double.NaN,Double.NaN,Double.NaN};return new double[]{a[0]/p,a[1]/p,a[2]/p};}
 static double clamp(double v){return Math.max(-1,Math.min(1,v));}
 static double angle(V a,V b){return Math.toDegrees(Math.acos(clamp(dot(a.xyz(),b.xyz())/(a.p()*b.p()))));}
 static double phi(V q,V electron,V h){double[] qh=unit(q.xyz()),n=unit(cross(q.xyz(),electron.xyz())),nh=unit(cross(q.xyz(),h.xyz()));double v=Math.atan2(dot(qh,cross(n,nh)),dot(n,nh));return v<0?v+2*Math.PI:v;}
 static void lab(List<Double> a,V p,double vz){a.add(p.p());a.add(Math.acos(clamp(p.z/p.p())));double f=Math.atan2(p.y,p.x);a.add(f<0?f+2*Math.PI:f);a.add(vz);}
 public static int helicity(int run,int raw){
  if((run>=4326&&run<=5666)||(run>=6616&&run<=6783)||(run>=6120&&run<=6604)||(run>=11284&&run<11300))return -raw;
  return raw;
 }
 public static double forwardT(double Q2,double W){double s=W*W,eq=(s-Q2-MN*MN)/(2*W),ep=(s+MPI*MPI-MP*MP)/(2*W);double pp=Math.sqrt(Math.max(0,ep*ep-MPI*MPI));return -Q2+MPI*MPI-2*eq*ep+2*Math.sqrt(eq*eq+Q2)*pp;}
 public static double[] calculate(double[] meta,V e,V p,V pi,double vze,double vzp,double vzpi,double c1,double c2,double b1,double b2){
  double eb=meta[3];V beam=V.track(0,0,Math.sqrt(eb*eb-ME*ME),ME),q=beam.sub(e),target=new V(MN,0,0,0),total=q.add(target);
  double Q2=-q.m2(),nu=q.e,y=nu/eb,W2=total.m2();
  
  double W=Math.sqrt(W2),x=Q2/(2*MN*nu),gamma=2*MN*x/Math.sqrt(Q2),g2=gamma*gamma;
  double t=q.sub(pi).m2(),tmin=W>=MP+MPI?forwardT(Q2,W):Double.NaN;V missing=total.sub(p).sub(pi);
  double[] boost={-total.x/total.e,-total.y/total.e,-total.z/total.e};V qc=q.boost(boost),ec=e.boost(boost),pc=p.boost(boost),pic=pi.boost(boost);
  double phi1=phi(qc,ec,pc),phi2=phi(qc,ec,pic),pt2=Math.sqrt(dot(cross(unit(qc.xyz()),pic.xyz()),cross(unit(qc.xyz()),pic.xyz())));
  double depRoot=Math.sqrt(1-y-y*y*g2/4);
  List<Double> a=new ArrayList<>();for(double v:meta)a.add(v);lab(a,e,vze);lab(a,p,vzp);lab(a,pi,vzpi);
  for(double v:new double[]{c1,c2,b1,b2,Q2,W,x,y,nu,t,tmin,t-tmin,missing.m2(),total.sub(p).m2(),total.sub(pi).m2(),missing.e,Math.hypot(missing.x,missing.y),Math.sqrt(Math.max(0,p.add(pi).m2())),pi.e/nu,2*dot(pic.xyz(),unit(qc.xyz()))/W,pt2,phi1,phi2,phi2-phi1,angle(e,p),angle(e,pi),(1-y+y*y/2+y*y*g2/4)/(1+g2),(1-y-y*y*g2/4)/(1+g2),y*(1-y/2)/Math.sqrt(1+g2),(2-y)*depRoot/(1+g2),y*depRoot/Math.sqrt(1+g2)})a.add(v);
  double[] out=new double[a.size()];for(int i=0;i<out.length;i++){out[i]=a.get(i);if(Double.isInfinite(out[i]))out[i]=Double.NaN;}return out;
 }
 public static double[] append(double[] base,double[] extra){double[] out=Arrays.copyOf(base,INPUT_SCHEMA.split(" ").length);System.arraycopy(extra,0,out,base.length,extra.length);int n=out.length;out[n-3]=base[28]>=1&&base[29]>=2&&base[31]<=.8&&base[32]>0?1:0;out[n-2]=base[28]>1.5?1:0;out[n-1]=1;for(int i=28;i<55;i++)if(!Double.isFinite(base[i]))out[n-1]=0;return out;}

}
