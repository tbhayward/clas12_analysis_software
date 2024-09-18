// dump QADB data for a specified file; this is for QADB validation and
// testing, but may still be useful as an example; see dumpQADB.cpp for a more
// user-friendly example
// - this program does requires a HIPO file

// imports
import clasqa.QADB // access QADB
import org.jlab.io.hipo.HipoDataSource // to read HIPO files

int maxbit = 5; // max defect bit number

// print separator
def sep = { s,n -> n.times{print s}; println ""; }

// print error
def err = { s -> System.err << "ERROR: $s\n" }


// instantiate HIPO file reader for specified HIPO file
def inHipoFile
if(args.length>=1) inHipoFile = args[0]
else { System.err << "ERROR: specify hipo file\n"; return; }
def reader = new HipoDataSource()
reader.open(inHipoFile)


// instantiate QADB
QADB qa = new QADB()
// alternatively, specify run range to restrict QADB (may be more efficient)
//QADB qa = new QADB(5000,5500);


// event loop
int evnum,runnum,filenum
int filenumTmp = -1000
def defname
int chargeInt
while(reader.hasEvent()) {

  event = reader.getNextEvent()

  // get run and event numbers
  runnum = event.getBank("RUN::config").getInt('run',0)
  evnum = event.getBank("RUN::config").getInt('event',0)

  // query
  qa.query(runnum,evnum)

  // skip tag1 events
  if(runnum==0) continue

  // apply QA cut
  //if(!(qa.OkForAsymmetry(runnum,evnum))) continue

  // get filenum; only print info if we queried a different file
  filenum = qa.getFilenum()
  if(filenum==filenumTmp) continue
  filenumTmp = filenum

  sep("=",50)
  //err("test error print")

  // check run and file number accessors: make sure that they
  // are equal to what we asked for
  println "- run,evnum,evnumRange,filenum"
  println qa.getRunnum() + " " + runnum
  if(qa.getRunnum() != runnum) err("QADB::getRunnum != runnum");

  // check event number: report an error if evnum min>=max
  println evnum
  println qa.getEvnumMin() + " " + qa.getEvnumMax()
  if(qa.getEvnumMin() >= qa.getEvnumMax())
    err("GetEvnumMin() >= GetEvnumMax()");
  println filenum

  // print charge (convert to pC and truncate, for easier comparison)
  println "- charge,comment"
  println ((int)(1000*qa.getCharge()))
  qa.accumulateCharge();

  // print comment
  println "\"" + qa.getComment() + "\""


  // print overall defect info
  println "- defect"
  println qa.getDefect()
  for(int sec=1; sec<=6; sec++) 
    print " " + qa.getDefect(sec); println ""

  // print defect bit info
  for(int bit=0; bit<=maxbit; bit++) {
    // translate bit number to name; check if QADB::Bit returns correct bit
    switch(bit) {
      case 0: defname="TotalOutlier"; break;
      case 1: defname="TerminalOutlier"; break;
      case 2: defname="MarginalOutlier"; break;
      case 3: defname="SectorLoss"; break;
      case 4: defname="LowLiveTime"; break;
      case 5: defname="Misc"; break;
    };
    if(qa.Bit(defname) != bit) err("QADB::Bit problem");
    println qa.Bit(defname) + " $bit $defname"
    // print defect info
    println (qa.hasDefect(defname)?1:0)
    for(int sec=1; sec<=6; sec++) 
      print " " + (qa.hasDefect(defname,sec)?1:0); println ""
    // print bit masking
    qa.setMaskBit(defname);
    println qa.getMask() + " " + (qa.pass(runnum,evnum)?1:0)
    qa.setMaskBit(defname,false);
  };

  // print QA cuts (see above for custom cut check with mask)
  println "- cuts"
  println (qa.golden(runnum,evnum)?1:0)
  println (qa.OkForAsymmetry(runnum,evnum)?1:0)
}
sep("=",50);

// print accumulated charge
println "- charge"
println ((int)(1000*qa.getAccumulatedCharge()))

sep("=",50);
