// - used by bin/MakeLatexTables.sh to generate latex tables
// - cf. makeLatexTable2.groovy

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import clasqa.Tools
Tools T = new Tools()


def dbDirN = System.getenv('QADB') + '/qadb'
if(dbDirN==null) {
  System.err << "ERROR: env var QADB not set; source env.sh\n\n\n"
  return
}

def qaTree = [:]
def slurper = new JsonSlurper()
def dbDir = new File(dbDirN)
def dbFilter = ~/^qaTree.json$/
dbDir.traverse( 
  type:groovy.io.FileType.FILES,
  maxDepth:1,
  nameFilter:dbFilter) { dbFile ->
  println "read " + dbFile.getAbsoluteFile()
  slurper.parse(dbFile).each{ obj -> qaTree << obj }
}

outfileN="${dbDirN}/qaTable.tex"

def outfileF = new File(outfileN)
def outfileW = outfileF.newWriter(false)
def nGolden, nTotal
def frac

def bits = 0..<(T.bitDefinitions.size())
def nDefect = bits.collect{0}
def defect

def writeStr = ""

def rowStr = []
def runStr
qaTree.sort{a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
  run, runTree ->

  nGolden = 0
  nTotal = 0
  nDefect = bits.collect{0}
  runTree.sort{a,b -> a.key.toInteger() <=> b.key.toInteger() }.each{
    file, fileTree ->
    defect = fileTree.defect.toInteger()
    if(defect==0) nGolden++
    else bits.each{ bit -> if((defect>>bit)&0x1) nDefect[bit]++ }
    nTotal++
  }

  runStr = "${run}"
  if( run.toInteger() in [ 5128, 5129, 5130, 5158,
               5159, 5160, 5163, 5165,
               5166, 5167, 5168, 5169,
               5180, 5181, 5182, 5183,
               5567 ]) runStr+='^*'

  frac = sprintf("%.2f",100*nGolden/nTotal)
  writeStr = "\$$runStr\$ & \$${frac}\\\\%\$"
  nDefect.each{ n ->
    frac = sprintf("%.2f",100*n/nTotal)
    writeStr += " & \$${frac}\\\\%\$"
  }
  writeStr += " \\\\\\\\\\\\hline\n"
  outfileW << writeStr
}

outfileW.close()
