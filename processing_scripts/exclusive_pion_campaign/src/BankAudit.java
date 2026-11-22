package exclusivepion;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.base.DataEvent;
import org.jlab.io.base.DataBank;
import java.util.*;
/** Sample records directly: no fitter, corrections, QADB, or physics selection. */
public final class BankAudit {
 public static void main(String[] args){
  if(args.length!=2)throw new IllegalArgumentException("input.hipo max-records");
  String[] names={"RUN::config","REC::Event","REC::Particle","REC::Calorimeter","REC::Traj","REC::Track","REC::Cherenkov"};
  long[] missing=new long[names.length],empty=new long[names.length];long seen=0,complete=0,electrons=0;Map<String,Long> patterns=new LinkedHashMap<>();
  HipoDataSource reader=new HipoDataSource();reader.open(args[0]);
  try{while(reader.hasEvent()&&seen<Long.parseLong(args[1])){
   DataEvent event=reader.getNextEvent();seen++;StringJoiner absent=new StringJoiner(",");boolean all=true;
   for(int i=0;i<names.length;i++){
    if(!event.hasBank(names[i])){missing[i]++;all=false;absent.add(names[i]);}
    else if(event.getBank(names[i]).rows()==0)empty[i]++;
   }
   if(all)complete++;
   String key=all?"all required banks present":"missing "+absent;patterns.put(key,patterns.getOrDefault(key,0L)+1);
   if(event.hasBank("REC::Particle")){
    DataBank particles=event.getBank("REC::Particle");boolean electron=false;
    for(int i=0;i<particles.rows();i++)if(particles.getInt("pid",i)==11){electron=true;break;}
    if(electron)electrons++;
   }
   if(seen<=3){System.out.println("SAMPLE record "+seen+": "+key);String[] banks=event.getBankList();if(banks!=null)System.out.println("  present bank list: "+Arrays.toString(banks));}
  }}finally{reader.close();}
  System.out.println("AUDIT includes first record; sampled="+seen+" complete_bank_sets="+complete+" records_with_electron="+electrons);
  for(int i=0;i<names.length;i++)System.out.println(names[i]+": absent="+missing[i]+" present_but_empty="+empty[i]);
  for(Map.Entry<String,Long> p:patterns.entrySet())System.out.println(p.getValue()+" records: "+p.getKey());
 }
}
