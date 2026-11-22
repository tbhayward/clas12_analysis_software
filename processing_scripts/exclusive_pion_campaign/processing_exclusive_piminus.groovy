import exclusivepion.*
import analyzers.BeamEnergy
import org.jlab.clas.physics.PhysicsEvent
import org.jlab.io.hipo.*
import clasqa.QADB

// Arguments: exact input-list file, output text file, required run number.
if(args.length!=3)throw new IllegalArgumentException('Expected input-list output.txt run')
def inputs=new File(args[0]).readLines().findAll{it.trim()}
int expectedRun=Integer.parseInt(args[2])
def qa=new QADB('latest')
['TotalOutlier','TerminalOutlier','MarginalOutlier','SectorLoss','LowLiveTime','Misc','ChargeHigh','ChargeNegative','ChargeUnknown'].each{qa.checkForDefect(it)}
[6736,6737,6738,6739,6740,6741,6742,6743,6744,6746,6747,6748,6749,6750,6751,6753,6754,6755,6756,6757,16194,16089,16185,16308,16184,16307,16309,16872,16975,17763,17764,17765,17766,17767,17768,17179,17180,17181,17182,17183,17188,17189,17252].each{qa.allowMiscBit(it)}
def excluded=[5247,5345,17331,16987,17079,17190,17639,16850,16851,16852,16855,16879] as Set
def fitter=new ExclusivePionFitter()
long seen=0,rows=0,firstRecords=0,missingBanks=0
Double energy=null
new File(args[1]).withWriter('UTF-8'){writer->
 writer.write('# '+ExclusivePionVariables.INPUT_SCHEMA+'\n')
 inputs.each{path->
  def reader=new HipoDataSource()
  try {
   println('Opening '+path);reader.open(path)
   // Preserve the established train-file convention: first record is metadata.
   if(reader.hasEvent()){reader.getNextEvent();firstRecords++}
   while(reader.hasEvent()){
    def event=reader.getNextEvent();seen++
    if(!event.hasBank('RUN::config')){missingBanks++;continue}
    def cfg=event.getBank('RUN::config');int run=cfg.getInt('run',0),ev=cfg.getInt('event',0)
    if(run!=expectedRun)throw new IllegalStateException("Input inventory run ${expectedRun} does not match HIPO run ${run}: ${path}")
    if((run>16600&&run<16700)||(run>17768&&run<=17811)||excluded.contains(run))continue
    def names=['REC::Event','REC::Particle','REC::Calorimeter','REC::Traj','REC::Track','REC::Cherenkov']
    if(!names.every{event.hasBank(it)}){missingBanks++;continue}
    int h=event.getBank('REC::Event').getByte('helicity',0)
    if(h==0||!qa.pass(run,ev))continue
    if(energy==null)energy=new BeamEnergy(new PhysicsEvent(),run,false).Eb()
    def s=fitter.select(event.getBank('REC::Particle'),event.getBank('REC::Calorimeter'),event.getBank('REC::Track'),event.getBank('REC::Traj'),cfg,event.getBank('REC::Cherenkov'))
    if(s.electron==null||s.protons.isEmpty()||s.pions.isEmpty())continue
    def e=s.electron
    s.protons.each{p->s.pions.each{pi->
     double[] meta=[run,ev,ExclusivePionVariables.helicity(run,h),energy,p.index,pi.index,p.detector,pi.detector,e.sector,p.sector,pi.sector,ExclusivePionFitter.fiducial(e,p,pi)] as double[]
     def values=ExclusivePionVariables.calculate(meta,e.v,p.v,pi.v,e.vz,p.vz,pi.vz,p.chi2,pi.chi2,p.beta,pi.beta)
     if(values!=null){
      def line=new StringBuilder(1024)
      for(int i=0;i<values.length;i++){if(i>0)line.append(' ');line.append(Double.toString(values[i]))}
      line.append('\n');writer.write(line.toString());rows++
     }
    }}
   }
  } finally {reader.close()}
 }
}
println("SUMMARY seen=${seen} first_records=${firstRecords} missing_banks=${missingBanks} rows=${rows}")
