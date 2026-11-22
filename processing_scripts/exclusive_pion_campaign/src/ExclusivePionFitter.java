package exclusivepion;
import extended_kinematic_fitters.*;
import org.jlab.io.hipo.HipoDataBank;
import java.util.*;
import static exclusivepion.ExclusivePionVariables.*;
/** Reuses existing electron/PID/correction helpers; defaults are untouched. */
public final class ExclusivePionFitter extends analysis_fitter {
 public ExclusivePionFitter(){super(10.6041);}
    @Override
    public boolean proton_test(int particle_Index, int pid, float vz, double trigger_electron_vz,
            HipoDataBank rec_Bank, HipoDataBank cal_Bank,
            HipoDataBank traj_Bank, HipoDataBank run_Bank) {

        generic_tests generic_tests = new generic_tests();
        fiducial_cuts fiducial_cuts = new fiducial_cuts();
        pid_cuts pid_cuts = new pid_cuts();

        float torus = run_Bank.getFloat("torus", 0);

        float px = rec_Bank.getFloat("px", particle_Index);
        float py = rec_Bank.getFloat("py", particle_Index);
        float pz = rec_Bank.getFloat("pz", particle_Index);
        double p = Math.sqrt(Math.pow(px, 2) + Math.pow(py, 2) + Math.pow(pz, 2));
        boolean passesForwardDetector = generic_tests.forward_detector_cut(particle_Index, rec_Bank);
        boolean passesCentralDetector = generic_tests.central_detector_cut(particle_Index, rec_Bank);

        return true
                                && (passesCentralDetector ? p > 0.3 : true)
                                && (passesForwardDetector && (torus > 0) ? p > 0.42 : true)
                                && (passesForwardDetector && (torus < 0) ? p > 0.50 : true)
                                && (passesForwardDetector || passesCentralDetector) // campaign: no upper proton momentum cut
                                && generic_tests.theta_calculation(px, py, pz) < 64.23
                                && generic_tests.vertex_cut(particle_Index, rec_Bank, run_Bank)
                && (passesForwardDetector
                        ? fiducial_cuts.dc_fiducial_cut(particle_Index, rec_Bank, traj_Bank, run_Bank)
                        : true)
                && (passesCentralDetector
                        ? fiducial_cuts.cvt_fiducial_cut(particle_Index, rec_Bank, traj_Bank, 2)
                        : true) //               
//                && (passesForwardDetector // dedicated PID cuts for forward
//                                        ? pid_cuts.charged_hadron_chi2pid_cut(particle_Index, rec_Bank, run_Bank)
//                                        : true)
//                                && (passesCentralDetector
//                                        ? pid_cuts.charged_hadron_chi2pid_cut(particle_Index, rec_Bank, run_Bank)
//                                        : true) //            && charged_hadron_chi2pid_cut(particle_Index, rec_Bank)
                ;
    }


 public static final class Track {public int index,detector,sector;public V v;public double vz,chi2,beta;public boolean fid;}
 public static final class Selected {public Track electron;public List<Track> protons=new ArrayList<>(),pions=new ArrayList<>();}
 private Track track(int i,HipoDataBank r,HipoDataBank tr,HipoDataBank cfg,HipoDataBank traj,boolean proton){
  Track t=new Track();t.index=i;t.detector=new generic_tests().central_detector_cut(i,r)?2:1;t.sector=new generic_tests().sector(i,tr);
  float[] mom={r.getFloat("px",i),r.getFloat("py",i),r.getFloat("pz",i)};
  if(proton)new energy_loss_corrections().proton_energy_loss_corrections(i,mom,r,cfg);
  t.v=V.track(mom[0],mom[1],mom[2],proton?MP:MPI);t.vz=r.getFloat("vz",i);t.chi2=r.getFloat("chi2pid",i);t.beta=r.getFloat("beta",i);
  fiducial_cuts f=new fiducial_cuts();t.fid=t.detector==2?f.cvt_fiducial_cut(i,r,traj,1):f.dc_fiducial_cut(i,r,traj,cfg);return t;
 }
 public Selected select(HipoDataBank r,HipoDataBank cal,HipoDataBank tr,HipoDataBank traj,HipoDataBank cfg,HipoDataBank cc){
  Selected s=new Selected();if(r.rows()==0||r.getInt("pid",0)!=11)return s;
  double px=r.getFloat("px",0),py=r.getFloat("py",0),pz=r.getFloat("pz",0);V e=V.track(px,py,pz,ME);
  if(!electron_test(0,e.p(),r,cal,traj,cfg,cc))return s;
  Track et=new Track();et.index=0;et.v=e;et.vz=r.getFloat("vz",0);et.detector=1;et.sector=new generic_tests().sector(0,tr);fiducial_cuts f=new fiducial_cuts();et.fid=f.pcal_fiducial_cut(0,1,cfg,r,cal)&&f.dc_fiducial_cut(0,r,traj,cfg);
  double scale=RgcElectronCorrection.scale(cfg.getInt("run",0),et.sector,Math.toDegrees(Math.acos(e.z/e.p())),cfg.getFloat("solenoid",0));
  et.v=V.track(e.x*scale,e.y*scale,e.z*scale,ME);s.electron=et;
  generic_tests g=new generic_tests();
  for(int i=1;i<r.rows();i++){int pid=r.getInt("pid",i);float vz=r.getFloat("vz",i);
   if(pid==-211&&g.forward_detector_cut(i,r)&&pion_test(i,pid,vz,et.vz,r,cal,traj,cfg))s.pions.add(track(i,r,tr,cfg,traj,false));
   else if(pid==2212&&proton_test(i,pid,vz,et.vz,r,cal,traj,cfg))s.protons.add(track(i,r,tr,cfg,traj,true));
  }return s;
 }
 public static int fiducial(Track e,Track p,Track pi){if(e.fid&&p.fid&&pi.fid)return 3;if(!e.fid&&p.fid&&pi.fid)return 0;if(e.fid&&!p.fid&&pi.fid)return 1;if(e.fid&&p.fid&&!pi.fid)return 2;return -1;}
}
