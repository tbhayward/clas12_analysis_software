// print list of golden and silver runs stored in entire QADB
// - golden: no defects
// - silver: allows terminal outliers, including cases where the terminal
//           outlier occurs with a marginal outler

import clasqa.QADB

QADB qa = new QADB()
def gold,silver
def defect

// open output file writer
def outFileObj = new File(System.getenv('QADB')+"/text/listOfGoldenRuns.txt")
def outFile = outFileObj.newWriter(false)

// silver masks
def terminalBitOnly = 0x1 << qa.Bit('TerminalOutlier')
def terminalAndMarginalOnly = 
  terminalBitOnly + (0x1 << qa.Bit('MarginalOutlier'))

// run loop (sorted)
qa.getQaTree().sort{ a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
  runnum,runQA ->

  gold = true
  silver = true

  // file loop (sorted)
  runQA.sort{ c,d -> c.key.toInteger() <=> d.key.toInteger() }.each{
    filenum, fileQA ->

    defect = fileQA.defect
    if(defect>0) {
      gold = false
      if(defect!=terminalBitOnly && defect!=terminalAndMarginalOnly) {
        silver = false
        return
      }
    }
  } 

  outFile << " $runnum "
  if(gold) outFile << " gold (and silver)\n"
  else if(silver) outFile << " silver\n"
  else outFile << " defect\n"
} // end run loop

outFile.close()
