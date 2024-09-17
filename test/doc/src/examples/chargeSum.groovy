// calculate total analyzed charge for an example event loop
// with QA cuts enabled
// - you must specify a hipo file as an argument

// imports
import org.jlab.io.hipo.HipoDataSource // to read HIPO files
import clasqa.QADB // access QADB


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


// define variables
def evCount = 0
def event
def runnum,evnum


// event loop
println "begin event loop..."
while(reader.hasEvent()) {
  if(evCount%10000==0) println "$evCount events analyzed"

  // truncate event loop (for quick testing)
  if(evCount>1e5) { println "event loop truncated!"; break; }

  event = reader.getNextEvent()

  // get run and event numbers
  runnum = event.getBank("RUN::config").getInt('run',0)
  evnum = event.getBank("RUN::config").getInt('event',0)

  // QA cuts
  if(qa.OkForAsymmetry(runnum,evnum)) {

    // accumulate charge; note that although the call to
    // QADB::accumulateCharge() charge happens for each
    // event within a DST file that passed the QA cuts, that
    // file's charge will only be accumulated once, so
    // overcounting is not possible 
    qa.accumulateCharge()

    /* continue your analysis here */

  }

  evCount++
}

// print charge
println "\ntotal accumulated charge analyzed:\nrun=$runnum  charge=" +
  qa.getAccumulatedCharge() + " nC"
