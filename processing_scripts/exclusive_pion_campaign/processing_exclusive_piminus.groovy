import exclusivepion.*
import analyzers.BeamEnergy
import org.jlab.clas.physics.PhysicsEvent
import org.jlab.io.hipo.*
import clasqa.QADB
import groovy.json.JsonOutput
import java.nio.file.Files
import java.nio.file.StandardCopyOption

// Arguments: exact input-list file, output text file, required run number.
if(args.length<3||args.length>4)throw new IllegalArgumentException('Expected input-list output.txt run')
def inputs=new File(args[0]).readLines().findAll{it.trim()}
int expectedRun=Integer.parseInt(args[2])
def counts=[seen:0L,first_records:0L,missing_banks:0L,non_physics:0L,excluded_run:0L,zero_helicity:0L,qa_reject:0L,qa_pass:0L,no_electron:0L,no_proton:0L,no_pion:0L,pairs:0L,kinematic_reject:0L,accepted_events:0L,rows:0L]
def progressFile=new File(args.length==4?args[3]:args[1]+'.progress.json')
long begun=System.currentTimeMillis(),lastUpdate=0
String stage='qadb_startup'
def publish={ boolean force ->
 long now=System.currentTimeMillis()
 if(force||now-lastUpdate>=15000){
  def tmp=new File(progressFile.path+'.tmp')
  tmp.text=JsonOutput.toJson([stage:stage,updated:now/1000.0,elapsed_seconds:(now-begun)/1000.0,counts:counts])
  Files.move(tmp.toPath(),progressFile.toPath(),StandardCopyOption.REPLACE_EXISTING)
  lastUpdate=now
 }
}
publish(true)
def qa=new QADB('latest')
['TotalOutlier','TerminalOutlier','MarginalOutlier','SectorLoss','LowLiveTime','Misc','ChargeHigh','ChargeNegative','ChargeUnknown'].each{qa.checkForDefect(it)}
[6736,6737,6738,6739,6740,6741,6742,6743,6744,6746,6747,6748,6749,6750,6751,6753,6754,6755,6756,6757,16194,16089,16185,16308,16184,16307,16309,16872,16975,17763,17764,17765,17766,17767,17768,17179,17180,17181,17182,17183,17188,17189,17252].each{qa.allowMiscBit(it)}
def excluded=[5247,5345,17331,16987,17079,17190,17639,16850,16851,16852,16855,16879] as Set
def fitter=new ExclusivePionFitter()
stage='processing';publish(true)
def requiredBanks=['RUN::config','REC::Event','REC::Particle','REC::Calorimeter','REC::Traj','REC::Track','REC::Cherenkov']
requiredBanks.each{counts['missing_'+it.replace('::','_')]=0L}
Double energy=null
new File(args[1]).withWriter('UTF-8'){writer->
 writer.write('# '+ExclusivePionVariables.INPUT_SCHEMA+'\n')
 inputs.each{path->
  def reader=new HipoDataSource()
  try {
   println('Opening '+path);reader.open(path)
   // Preserve the established train-file convention: first record is metadata.
   if(reader.hasEvent()){reader.getNextEvent();counts.first_records++}
   while(reader.hasEvent()){
    def event=reader.getNextEvent();counts.seen++;if(counts.seen%1000==0)publish(false)
    boolean missing=false
    for(String bank:requiredBanks){if(!event.hasBank(bank)){counts['missing_'+bank.replace('::','_')]++;missing=true}}
    if(missing){counts.missing_banks++;continue}
    def cfg=event.getBank('RUN::config')
    if(cfg.rows()==0||event.getBank('REC::Event').rows()==0||event.getBank('REC::Particle').rows()==0){counts.non_physics++;continue}
    int run=cfg.getInt('run',0),ev=cfg.getInt('event',0)
    if(!PhysicsRecord.acceptRun(run,expectedRun)){counts.non_physics++;continue}
    if((run>16600&&run<16700)||(run>17768&&run<=17811)||excluded.contains(run)){counts.excluded_run++;continue}
    int h=event.getBank('REC::Event').getByte('helicity',0)
    if(h==0){counts.zero_helicity++;continue}
    if(!qa.pass(run,ev)){counts.qa_reject++;continue}
    counts.qa_pass++
    if(energy==null)energy=new BeamEnergy(new PhysicsEvent(),run,false).Eb()
    def s=fitter.select(event.getBank('REC::Particle'),event.getBank('REC::Calorimeter'),event.getBank('REC::Track'),event.getBank('REC::Traj'),cfg,event.getBank('REC::Cherenkov'))
    if(s.electron==null){counts.no_electron++;continue}
    if(s.protons.isEmpty()){counts.no_proton++;continue}
    if(s.pions.isEmpty()){counts.no_pion++;continue}
    def e=s.electron
    long rowsBefore=counts.rows
    s.protons.each{p->s.pions.each{pi->
     counts.pairs++
     double[] meta=[run,ev,ExclusivePionVariables.helicity(run,h),energy,p.index,pi.index,p.detector,pi.detector,e.sector,p.sector,pi.sector,ExclusivePionFitter.fiducial(e,p,pi)] as double[]
     def values=ExclusivePionVariables.calculate(meta,e.v,p.v,pi.v,e.vz,p.vz,pi.vz,p.chi2,pi.chi2,p.beta,pi.beta)
     if(values!=null){
      def line=new StringBuilder(1024)
      for(int i=0;i<values.length;i++){if(i>0)line.append(' ');line.append(Double.toString(values[i]))}
      line.append('\n');writer.write(line.toString());counts.rows++
     }else{counts.kinematic_reject++}
    }}
    if(counts.rows>rowsBefore)counts.accepted_events++
   }
  } finally {reader.close()}
 }
}
stage='processed';publish(true)
println('SUMMARY '+counts.collect{k,v->"${k}=${v}"}.join(' '))
