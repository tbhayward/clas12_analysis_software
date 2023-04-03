// - used by bin/MakeLatexTables2.sh to generate latex tables
// - cf. makeLatexTable.groovy

import clasqa.QADB

// arguments: run range
def runLB,runUB
if(args.length==2) {
  runLB = args[0].toInteger()
  runUB = args[1].toInteger()
}
else { System.err << "ARGUMENTS: [runLB] [runUB]\n"; return; }
println "run range: $runLB to $runUB"
QADB qa = new QADB(runLB,runUB)

def nPassQA, nTotal
def exList = []
def fieldList = []
int evnum,runnum

// output latex file
def dbDirN = System.getenv('QADB') + '/qadb'
if(dbDirN==null) {
  System.err << "ERROR: env var QADB not set; source env.sh\n\n\n"
  return
}
outfileN="${dbDirN}/qaTable.tex"
def outfileF = new File(outfileN)
def outfileW = outfileF.newWriter(false)

// loop over runs
qa.getQaTree().sort{a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
  runnumStr, runTree ->
  runnum = runnumStr.toInteger()

  // resets
  nPassQA = 0
  nTotal = 0
  exList.clear()
  fieldList.clear()

  // loop over files
  runTree.sort{a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
    filenumStr, fileTree ->
    filenum = filenumStr.toInteger()
    evnum = fileTree['evnumMin']
    if(qa.OkForAsymmetry(runnum,evnum)) nPassQA++ // QA CUT /////
    else exList << filenum
    nTotal++
  }

  // include run iff some files passed QA cut
  if(nPassQA>0) {

    // prepare latex table fields
    def frac = sprintf("%.2f",100*nPassQA/nTotal)
    def exListStr = exList.join(', ')
    fieldList << "$runnum"
    fieldList << "${frac}\\\\%"
    fieldList << "$exListStr"

    // latex output
    def texStr = fieldList.join(' & ')
    texStr += " \\\\\\\\\\\\hline\n"
    outfileW << texStr
  }
}

outfileW.close()
