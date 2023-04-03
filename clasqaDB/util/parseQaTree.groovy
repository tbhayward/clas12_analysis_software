// parses qa/qaTree.json into human readable format

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import clasqa.Tools
Tools T = new Tools()
/*
import clasqa.QADB
QADB qa = new QADB()
*/

infile="qa/qaTree.json"
outfile="qa/qaTable.dat"
if(args.length>=1) {
  infile = args[0]
  outfile = "${infile}.table"
}
println "reading $infile ..."

def outfileF = new File(outfile)
def outfileW = outfileF.newWriter(false)

def slurper = new JsonSlurper()
def jsonFile = new File(infile)
def qaTree = slurper.parse(jsonFile)
def defStr = []
qaTree.sort{a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
  run, runTree ->
  outfileW << "\nRUN: $run\n"
  runTree.sort{a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
    file, fileTree ->
    def defect = fileTree.defect
    //defStr=[run,file,defect,Integer.toBinaryString(defect)]
    defStr = [run,file]

    def getSecList = { bitNum ->
      def secList = []
      fileTree.sectorDefects.each{
        if(bitNum in it.value) secList+=it.key
      }
      return secList
    }

    if(defect>0) {
      T.bitNames.eachWithIndex { str,i ->
        if(defect >> i & 0x1) defStr += " ${i}-" + str + getSecList(i)
      }
    } else defStr += " GOLDEN"

    if(fileTree.comment!=null) {
      if(fileTree.comment.length()>0) defStr += " :: " + fileTree.comment
    }

    /*
    int runnum = run.toInteger()
    int evnum = fileTree['evnumMin'].toInteger()
    if(qa.OkForAsymmetry(runnum,evnum)) {
      defStr += " -- OK FOR ASYM"
    }
    */

    outfileW << defStr.join(' ') << "\n"
    //outfileW << fileTree.sectorDefects << "\n"
  }
}

outfileW.close()
println "produced $outfile"
