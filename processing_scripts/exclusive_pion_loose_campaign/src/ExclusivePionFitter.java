package exclusivepion;
import extended_kinematic_fitters.*;
import org.jlab.io.hipo.HipoDataBank;
import java.util.*;
import static exclusivepion.ExclusivePionVariables.*;
/** Reuses existing electron/PID/correction helpers; defaults are untouched. */
public final class ExclusivePionFitter {
 public static final String EXTRA_SCHEMA="e_index detector_e status_e status1 status2 e_raw_px e_raw_py e_raw_pz p1_raw_px p1_raw_py p1_raw_pz p2_raw_px p2_raw_py p2_raw_pz chi2pid_e beta_e pass_e_vertex pass_e_sampling pass_e_diagonal pass_e_pcal1 pass_e_pcal2 pass_e_dc pass_p1_momentum pass_p1_theta pass_p1_vertex pass_p1_dc pass_p1_cvt1 pass_p1_cvt2 pass_p1_chi2pid pass_p2_pid pass_p2_vertex pass_p2_dc pass_p2_cvt1 pass_p2_cvt2 pass_e_legacy pass_p1_legacy pass_p2_legacy electron_correction proton_correction";
 public static final class Track {public int index,detector,sector,status,correction;public V v,raw;public double vz,chi2,beta;public boolean fid;public int[] flags;}
 public static final class Selected {public int invalidTracks=0,belowElectronThreshold=0;public List<Track> electrons=new ArrayList<>(),protons=new ArrayList<>(),pions=new ArrayList<>();}
 interface Check {boolean test();}
 static int flag(boolean available,Check c){if(!available)return -1;try{return c.test()?1:0;}catch(RuntimeException ex){return -1;}}
 static boolean all(int... flags){for(int f:flags)if(f!=1)return false;return true;}
 private Track track(int i,int pid,HipoDataBank r,HipoDataBank cal,HipoDataBank tr,HipoDataBank traj,HipoDataBank cfg,HipoDataBank cc){
  Track t=new Track();t.index=i;t.status=r.getShort("status",i);int st=Math.abs(t.status);t.detector=st>=2000&&st<4000?1:st>=4000&&st<5000?2:st>=1000&&st<2000?3:0;
  generic_tests g=new generic_tests();fiducial_cuts f=new fiducial_cuts();pid_cuts pc=new pid_cuts();
  t.sector=-1;if(tr!=null){try{t.sector=g.sector(i,tr);}catch(RuntimeException ex){t.sector=-1;}}double mass=pid==11?ME:pid==2212?MP:MPI;
  t.raw=V.track(r.getFloat("px",i),r.getFloat("py",i),r.getFloat("pz",i),mass);
  if(!Double.isFinite(t.raw.p())||!(t.raw.p()>0))return null;
  if(pid==11&&!(t.raw.p()>2))return null;
  t.v=t.raw;t.vz=r.getFloat("vz",i);t.chi2=r.getFloat("chi2pid",i);t.beta=r.getFloat("beta",i);

  int vertex=flag(Double.isFinite(t.vz),()->g.vertex_cut(i,r,cfg));
  int dc=flag(t.detector==1&&traj!=null,()->f.dc_fiducial_cut(i,r,traj,cfg));
  int cvt1=flag(t.detector==2&&traj!=null,()->f.cvt_fiducial_cut(i,r,traj,1));
  int cvt2=flag(t.detector==2&&traj!=null,()->f.cvt_fiducial_cut(i,r,traj,2));
  if(pid==11){
   int sf=flag(cal!=null,()->pc.calorimeter_sampling_fraction_cut(i,t.raw.p(),cfg,cal));
   int diagonal=flag(cal!=null,()->pc.calorimeter_diagonal_cut(i,t.raw.p(),cal,cfg));
   int pcal1=flag(cal!=null,()->f.pcal_fiducial_cut(i,1,cfg,r,cal));
   int pcal2=flag(cal!=null,()->f.pcal_fiducial_cut(i,2,cfg,r,cal));
   t.flags=new int[]{vertex,sf,diagonal,pcal1,pcal2,dc};t.fid=all(pcal1,dc);
   int run=cfg.getInt("run",0);boolean rgc=run>=16000&&run<=17811;
   t.correction=rgc?-1:0;
   if(rgc&&t.detector==1&&t.sector>=1&&t.sector<=6){try{double scale=RgcElectronCorrection.scale(run,t.sector,Math.toDegrees(Math.acos(t.raw.z/t.raw.p())),cfg.getFloat("solenoid",0));if(Double.isFinite(scale)&&scale>0){t.v=V.track(t.raw.x*scale,t.raw.y*scale,t.raw.z*scale,ME);t.correction=1;}}catch(IllegalArgumentException ex){t.correction=-1;}}
  }else if(pid==2212){
   double torus=cfg.getFloat("torus",0),mom=t.raw.p();int momentum=(t.detector==2?mom>.3:t.detector==1?(torus>0?mom>.42:torus<0?mom>.50:true):false)?1:0;
   t.flags=new int[]{momentum,Math.toDegrees(Math.acos(t.raw.z/mom))<64.23?1:0,vertex,dc,cvt1,cvt2,!Double.isFinite(t.chi2)?-1:Math.abs(t.chi2)<3?1:0};t.fid=(t.detector==1?dc:cvt1)==1;
   t.correction=-1;if(t.detector==1||t.detector==2){try{float[] m={(float)t.raw.x,(float)t.raw.y,(float)t.raw.z};new energy_loss_corrections().proton_energy_loss_corrections(i,m,r,cfg);V corrected=V.track(m[0],m[1],m[2],MP);if(Double.isFinite(corrected.p())&&corrected.p()>0){t.v=corrected;t.correction=1;}}catch(RuntimeException ex){t.correction=-1;}}
  }else{t.flags=new int[]{flag(Double.isFinite(t.chi2),()->pc.charged_hadron_pass2_chi2pid_cut(i,r)),vertex,dc,cvt1,cvt2};t.fid=(t.detector==1?dc:cvt1)==1;}
  return t;
 }
 public Selected select(HipoDataBank r,HipoDataBank cal,HipoDataBank tr,HipoDataBank traj,HipoDataBank cfg,HipoDataBank cc){Selected s=new Selected();for(int i=0;i<r.rows();i++){int pid=r.getInt("pid",i);if(pid!=11&&pid!=2212&&pid!=-211)continue;Track t=track(i,pid,r,cal,tr,traj,cfg,cc);if(t==null){double mom=Math.sqrt(Math.pow(r.getFloat("px",i),2)+Math.pow(r.getFloat("py",i),2)+Math.pow(r.getFloat("pz",i),2));if(pid==11&&Double.isFinite(mom)&&mom<=2)s.belowElectronThreshold++;else s.invalidTracks++;continue;}if(pid==11)s.electrons.add(t);else if(pid==2212)s.protons.add(t);else s.pions.add(t);}return s;}
 public static double[] extras(Track e,Track p,Track pi){List<Double> a=new ArrayList<>();for(int v:new int[]{e.index,e.detector,e.status,p.status,pi.status})a.add((double)v);for(Track t:new Track[]{e,p,pi})for(double v:new double[]{t.raw.x,t.raw.y,t.raw.z})a.add(v);a.add(e.chi2);a.add(e.beta);for(Track t:new Track[]{e,p,pi})for(int v:t.flags)a.add((double)v);
  a.add(e.index==0&&e.detector==1&&all(e.flags[0],e.flags[1],e.flags[2],e.flags[4],e.flags[5])?1.:0.);
  a.add(all(p.flags[0],p.flags[1],p.flags[2],p.detector==1?p.flags[3]:p.flags[5])?1.:0.);
  a.add(pi.detector==1&&pi.flags[0]==1?1.:0.);a.add((double)e.correction);a.add((double)p.correction);
  double[] out=new double[a.size()];for(int i=0;i<out.length;i++)out[i]=a.get(i);return out;
 }
 public static int fiducial(Track e,Track p,Track pi){if(e.fid&&p.fid&&pi.fid)return 3;if(!e.fid&&p.fid&&pi.fid)return 0;if(e.fid&&!p.fid&&pi.fid)return 1;if(e.fid&&p.fid&&!pi.fid)return 2;return -1;}
}
