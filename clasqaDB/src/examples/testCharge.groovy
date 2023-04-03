// for charge cross check
// - prints charge for a selection of runs, applying specified QA criteria

// imports
import clasqa.QADB // access QADB


// instantiate QADB
QADB qa = new QADB()
// alternatively, specify run range to restrict QADB (may be more efficient)
//QADB qa = new QADB(5000,5500);
int evnum
def defname
int chargeInt


// specify run number list
def runnumList = [
  5682,
  5683,
  5695,
  5751,
  5827
]


// specify QA criteria
qa.setMaskBit("SectorLoss")
qa.setMaskBit("MarginalOutlier")


// loop over runs
runnumList.each{ runnum ->

  // loop over files
  for(int filenum=0; filenum<=qa.getMaxFilenum(runnum); filenum+=5) {

    // query by file number
    qa.queryByFilenum(runnum,filenum)
    evnum = qa.getEvnumMin() // evnum needed for QA cut methods

    // accumulate charge, if QA criteria are satisfied
    if( 
      qa.pass(runnum,evnum) &&
      !( runnum==5827 && (filenum==10 || filenum==45 || filenum==50) )
    ) {
      qa.accumulateCharge()
    }

  } // end file loop

  // print this run's charge, and reset
  println "$runnum " + qa.getAccumulatedCharge()
  qa.resetAccumulatedCharge()

} // end run loop
