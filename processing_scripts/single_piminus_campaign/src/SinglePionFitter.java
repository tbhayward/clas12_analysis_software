package exclusivepion;
import extended_kinematic_fitters.*;
import org.jlab.io.hipo.HipoDataBank;
import java.util.*;
import static exclusivepion.SinglePionVariables.*;
/** Use the shared pi+ electron and pion tests. No proton requirement or veto. */
public final class SinglePionFitter {
 private final analysis_fitter reference;
 public SinglePionFitter(){this(new analysis_fitter(10.6));}
 SinglePionFitter(analysis_fitter reference){this.reference=reference;}
 private final generic_tests g=new generic_tests();
 private final fiducial_cuts f=new fiducial_cuts();
 private final pid_cuts pid=new pid_cuts();
 public static final class Track {public int index,sector,detector;public V v;public double vz,chi2,beta;public boolean dc,pcal1,pcal2;}
 public static final class Selected {public Track electron;public int protons;public List<Track> pions=new ArrayList<>();}
 public Selected select(HipoDataBank r,HipoDataBank cal,HipoDataBank tr,HipoDataBank traj,HipoDataBank cfg){
  Selected s=new Selected();if(r.rows()==0||r.getInt("pid",0)!=11)return s;
  V e=V.track(r.getFloat("px",0),r.getFloat("py",0),r.getFloat("pz",0),ME);
  if(!(reference.electron_test(0,e.p(),r,cal,traj,cfg,null)))return s;
  Track et=new Track();et.index=0;et.vz=r.getFloat("vz",0);et.sector=g.sector(0,tr);
  et.dc=f.dc_fiducial_cut(0,r,traj,cfg);et.pcal1=f.pcal_fiducial_cut(0,1,cfg,r,cal);et.pcal2=f.pcal_fiducial_cut(0,2,cfg,r,cal);
  double scale=RgcElectronCorrection.scale(cfg.getInt("run",0),et.sector,Math.toDegrees(Math.acos(e.z/e.p())),cfg.getFloat("solenoid",0));
  et.v=V.track(e.x*scale,e.y*scale,e.z*scale,ME);s.electron=et;
  for(int i=1;i<r.rows();i++){
   int p=r.getInt("pid",i);if(p==2212)s.protons++;
   if(p!=-211||!reference.pion_test(i,p,r.getFloat("vz",i),et.vz,r,cal,traj,cfg))continue;
   Track t=new Track();t.index=i;t.detector=g.forward_detector_cut(i,r)?1:(g.central_detector_cut(i,r)?2:0);t.sector=g.sector(i,tr);t.v=V.track(r.getFloat("px",i),r.getFloat("py",i),r.getFloat("pz",i),MPI);
   if(t.v.p()<=0)continue;t.vz=r.getFloat("vz",i);t.chi2=r.getFloat("chi2pid",i);t.beta=r.getFloat("beta",i);t.dc=t.detector==1&&f.dc_fiducial_cut(i,r,traj,cfg);s.pions.add(t);
  }return s;
 }
 public static int fiducial(Track e,Track pi){boolean ef=e.dc&&e.pcal1&&e.pcal2;if(ef&&pi.dc)return 3;if(!ef&&pi.dc)return 0;if(ef&&!pi.dc)return 1;return -1;}
}
