// calculate total analyzed charge for a run period,
// with specified QA cuts enabled
// note: if syncCheck.groovy errors are present in the run range,
// the final charge value might be a bit wrong...

import org.jlab.io.hipo.HipoDataSource
import clasqa.QADB


// arguments: run range, from runLB to runUB
def runLB,runUB
if(args.length==2) {
  runLB = args[0].toInteger()
  runUB = args[1].toInteger()
}
else { System.err << "ARGUMENTS: [runLB] [runUB]\n"; return; }
println "run range: $runLB to $runUB"
QADB qa = new QADB(runLB,runUB)


// loop over runs and files
int runnum
int evnum
qa.getQaTree().each{ runnumStr, runTree ->
  runnum = runnumStr.toInteger()
  println runnum
  runTree.each{ filenumStr, fileTree ->
    evnum = fileTree['evnumMin']

    // QA cut /////////////////
    //if(qa.query(runnum,evnum)) qa.accumulateCharge() // no qa cut
    //if(qa.golden(runnum,evnum)) qa.accumulateCharge()
    if(qa.OkForAsymmetry(runnum,evnum)) qa.accumulateCharge()
    ///////////////////////////

  }
}

// print charge
println "\ntotal accumulated charge: " +
  qa.getAccumulatedCharge()/1e6 + " mC"
