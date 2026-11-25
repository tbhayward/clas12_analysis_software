package exclusivepion;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.jnp.hipo4.data.*;
import static exclusivepion.ExclusivePionVariables.*;
public class LooseSelectionTest {
 static HipoDataBank bank(String name,int rows,String... columns){Schema.SchemaBuilder s=new Schema.SchemaBuilder(name,100,1);for(String c:columns){String[] x=c.split(":");s.addEntry(x[0],x[1],"");}return new HipoDataBank(new Bank(s.build(),rows));}
 static void require(boolean b,String m){if(!b)throw new AssertionError(m);}
 public static void main(String[] args){
  HipoDataBank r=bank("REC::Particle",6,"pid:I","status:S","px:F","py:F","pz:F","vz:F","chi2pid:F","beta:F");
  int[] pid={11,11,2212,-211,11,11};double[] momenta={2,3,.1,.2,2.1,1.9};
  for(int i=0;i<pid.length;i++){r.setInt("pid",i,pid[i]);r.setShort("status",i,(short)(i==3?4500:2500));r.setFloat("px",i,(float)(momenta[i]*.2));r.setFloat("pz",i,(float)(momenta[i]*Math.sqrt(.96)));r.setFloat("vz",i,100);r.setFloat("chi2pid",i,50);r.setFloat("beta",i,.8f);}
  // Exact threshold uses an exactly representable longitudinal vector.
  r.setFloat("px",0,0);r.setFloat("pz",0,2);
  HipoDataBank cfg=bank("RUN::config",1,"run:I","torus:F","solenoid:F");cfg.setInt("run",0,16962);cfg.setFloat("torus",0,-1);cfg.setFloat("solenoid",0,-1);
  ExclusivePionFitter.Selected s=new ExclusivePionFitter().select(r,null,null,null,cfg,null);
  require(s.electrons.size()==2,"Only raw electrons strictly above 2 GeV");require(s.electrons.get(0).index==1,"Nontrigger electron retained");
  require(s.protons.size()==1&&s.pions.size()==1,"Loose hadrons retained despite vertex, momentum, PID and missing fiducial information");
  require(s.pions.get(0).detector==2,"CD pion tagged and retained");require(s.protons.get(0).flags[0]==0,"Low proton momentum flagged");
  require(s.electrons.get(0).flags[1]==-1,"Missing calorimeter marked unavailable");require(s.electrons.get(0).correction==-1,"Unavailable electron correction leaves raw momentum");
  double[] base=calculate(new double[]{16962,1,1,10.55,2,3,1,2,-1,-1,-1,-1},s.electrons.get(0).v,s.protons.get(0).v,s.pions.get(0).v,100,100,100,50,50,.8,.8);
  double[] row=append(base,ExclusivePionFitter.extras(s.electrons.get(0),s.protons.get(0),s.pions.get(0)));require(row.length==INPUT_SCHEMA.split(" ").length,"Complete schema");
  V high=V.track(0,0,12,ME);base=calculate(new double[]{16962,1,1,10.55,2,3,1,2,-1,-1,-1,-1},high,s.protons.get(0).v,s.pions.get(0).v,0,0,0,0,0,1,1);require(base!=null,"NonDIS candidate is retained");
  row=append(base,ExclusivePionFitter.extras(s.electrons.get(0),s.protons.get(0),s.pions.get(0)));require(row[row.length-3]==0,"NonDIS flagged");require(row[row.length-1]==0,"Undefined kinematics flagged");
  System.out.println("Loose selection: threshold, multiple electrons, CD pion, low-p proton, unavailable banks and nonDIS retention passed");
 }
}
