package exclusivepion;
import java.util.*;
/** GeV, radians. Stationary neutron reference; signed tprime = t-tmin. */
public final class SinglePionVariables {
 public static final double MN=.9395654133, MP=.9382720813, MPI=.13957039, ME=.00051099895;
 public static final String INPUT_SCHEMA="runnum evnum helicity beam_energy p_index detector sector_e sector fiducial_status e_pcal_fiducial_1 e_pcal_fiducial_2 e_dc_fiducial pion_dc_fiducial rec_proton_count e_p e_theta e_phi vz_e p_p p_theta p_phi vz_p chi2pid beta Q2 W x y nu t tmin tprime Mx Mx2 Mx2_deuteron Mx2_proton Emiss missing_px missing_py missing_pz missing_p pTmiss z xF pT phi open_angle DepA DepB DepC DepV DepW";
 public static final Set<String> INTS=new HashSet<>(Arrays.asList("runnum evnum helicity p_index detector sector_e sector fiducial_status e_pcal_fiducial_1 e_pcal_fiducial_2 e_dc_fiducial pion_dc_fiducial rec_proton_count".split(" ")));
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
 static double[] unit(double[] a){double p=Math.sqrt(dot(a,a));if(!(p>0))throw new IllegalArgumentException("Degenerate angle");return new double[]{a[0]/p,a[1]/p,a[2]/p};}
 static double clamp(double v){return Math.max(-1,Math.min(1,v));}
 static double angle(V a,V b){return Math.toDegrees(Math.acos(clamp(dot(a.xyz(),b.xyz())/(a.p()*b.p()))));}
 static double phi(V q,V electron,V h){double[] qh=unit(q.xyz()),n=unit(cross(q.xyz(),electron.xyz())),nh=unit(cross(q.xyz(),h.xyz()));double v=Math.atan2(dot(qh,cross(n,nh)),dot(n,nh));return v<0?v+2*Math.PI:v;}
 static void lab(List<Double> a,V p,double vz){a.add(p.p());a.add(Math.acos(clamp(p.z/p.p())));double f=Math.atan2(p.y,p.x);a.add(f<0?f+2*Math.PI:f);a.add(vz);}
 public static int helicity(int run,int raw){
  if((run>=4326&&run<=5666)||(run>=6616&&run<=6783)||(run>=6120&&run<=6604)||(run>=11284&&run<11300))return -raw;
  return raw;
 }
 public static double forwardT(double Q2,double W){double s=W*W,eq=(s-Q2-MN*MN)/(2*W),ep=(s+MPI*MPI-MP*MP)/(2*W);double pp=Math.sqrt(Math.max(0,ep*ep-MPI*MPI));return -Q2+MPI*MPI-2*eq*ep+2*Math.sqrt(eq*eq+Q2)*pp;}
 public static final double MD=1.87561294257;
 public static double[] calculate(double[] meta,V e,V pi,double vze,double vzpi,double chi2,double beta){
  double eb=meta[3];V beam=V.track(0,0,Math.sqrt(eb*eb-ME*ME),ME),q=beam.sub(e),total=q.add(new V(MN,0,0,0));
  double Q2=-q.m2(),nu=q.e,y=nu/eb,W2=total.m2();
  if(Q2<1||W2<4||y>.8||nu<=0||pi.p()<=0)return null;
  double W=Math.sqrt(W2),x=Q2/(2*MN*nu),gamma=2*MN*x/Math.sqrt(Q2),g2=gamma*gamma;
  double t=q.sub(pi).m2(),tmin=forwardT(Q2,W);V missing=total.sub(pi);
  double[] boost={-total.x/total.e,-total.y/total.e,-total.z/total.e};V qc=q.boost(boost),ec=e.boost(boost),pic=pi.boost(boost);
  double phi,pt;
  try{phi=phi(qc,ec,pic);pt=Math.sqrt(dot(cross(unit(qc.xyz()),pic.xyz()),cross(unit(qc.xyz()),pic.xyz())));}
  catch(IllegalArgumentException ex){return null;}
  double depRoot=Math.sqrt(Math.max(0,1-y-y*y*g2/4));
  List<Double> a=new ArrayList<>();for(double v:meta)a.add(v);lab(a,e,vze);lab(a,pi,vzpi);a.add(chi2);a.add(beta);
  for(double v:new double[]{Q2,W,x,y,nu,t,tmin,t-tmin,Math.copySign(Math.sqrt(Math.abs(missing.m2())),missing.m2()),missing.m2(),q.add(new V(MD,0,0,0)).sub(pi).m2(),q.add(new V(MP,0,0,0)).sub(pi).m2(),missing.e,missing.x,missing.y,missing.z,missing.p(),Math.hypot(missing.x,missing.y),pi.e/nu,2*dot(pic.xyz(),unit(qc.xyz()))/W,pt,phi,angle(e,pi),(1-y+y*y/2+y*y*g2/4)/(1+g2),(1-y-y*y*g2/4)/(1+g2),y*(1-y/2)/Math.sqrt(1+g2),(2-y)*depRoot/(1+g2),y*depRoot/Math.sqrt(1+g2)})a.add(v);
  double[] out=new double[a.size()];for(int i=0;i<out.length;i++){out[i]=a.get(i);if(!Double.isFinite(out[i]))throw new IllegalArgumentException("Nonfinite "+INPUT_SCHEMA.split(" ")[i]);}return out;
 }
}
