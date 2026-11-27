package exclusivepion;
import extended_kinematic_fitters.analysis_fitter;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.jnp.hipo4.data.*;
public class SelectionTest {
 static HipoDataBank bank(String name,int rows,String... columns){Schema.SchemaBuilder s=new Schema.SchemaBuilder(name,100,1);for(String c:columns){String[] x=c.split(":");s.addEntry(x[0],x[1],"");}return new HipoDataBank(new Bank(s.build(),rows));}
 static void require(boolean v,String msg){if(!v)throw new AssertionError(msg);}
 public static void main(String[] args){
  // Inject only the electron decision to isolate topology and real pion PID.
  analysis_fitter ref=new analysis_fitter(10.6){public boolean electron_test(int i,double p,HipoDataBank r,HipoDataBank c,HipoDataBank t,HipoDataBank cfg,HipoDataBank cc){return p>2;}};
  HipoDataBank r=bank("REC::Particle",3,"pid:I","status:S","px:F","py:F","pz:F","vz:F","chi2pid:F","beta:F");
  for(int i=0;i<3;i++){r.setFloat("pz",i,i==0?4:1);r.setShort("status",i,(short)(i==1?4500:2500));r.setFloat("beta",i,1);}
  r.setInt("pid",0,11);r.setInt("pid",1,-211);r.setInt("pid",2,22);
  HipoDataBank cfg=bank("RUN::config",1,"run:I","torus:F","solenoid:F");cfg.setInt("run",0,6420);
  HipoDataBank tr=bank("REC::Track",0,"pindex:S","sector:B");
  HipoDataBank cal=bank("REC::Calorimeter",0,"pindex:S","layer:B","sector:B","lv:F","lw:F");
  HipoDataBank traj=bank("REC::Traj",0,"pindex:S","detector:B","layer:B","x:F","y:F","z:F");
  SinglePionFitter f=new SinglePionFitter(ref);
  SinglePionFitter.Selected s=f.select(r,cal,tr,traj,cfg);
  require(s.electron!=null&&s.pions.size()==1&&s.protons==0,"No detected proton needed");
  require(s.pions.get(0).detector==2,"CD pion retained and tagged");
  r.setInt("pid",2,2212);s=f.select(r,cal,tr,traj,cfg);
  require(s.pions.size()==1&&s.protons==1,"Extra proton does not veto event");
  r.setFloat("pz",0,2);require(f.select(r,cal,tr,traj,cfg).electron==null,"Electron reference decision enforced");
  System.out.println("Single-pion topology: no proton requirement/veto, CD tagging, electron decision passed");
 }
}
