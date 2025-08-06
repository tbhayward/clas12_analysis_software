#!/usr/bin/env groovy

import org.jlab.io.hipo.HipoDataSource
import clasqa.QADB

// path to your CSV
def csvPath = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

// runs for which we allow the Misc bit
def allowMiscRuns = [
     5046, 5047, 5051, /* … */ 17768
]

// control flag and current header text
boolean inBlock = false
String currentHeader = ""

// formatting helper (four decimal places)
def fmt4 = { double x -> String.format("%.4f", x) }

// print intro
println "Reading CSV: $csvPath"
println "Output: runnum,totalCharge(nC),posCharge(nC),negCharge(nC),pol,polErr"
println ""

new File(csvPath).eachLine { lineRaw ->
  String line = lineRaw.trim()

  // detect new RGC block header
  if ( line.startsWith("# RGC Su22") ||
       line.startsWith("# RGC Fa22") ||
       line.startsWith("# RGC Sp23") ) {
    inBlock = true
    currentHeader = line.substring(1).trim()    // drop the leading ‘#’
    println "\n# ${currentHeader}"
    return
  }

  // detect end of RGC blocks
  if ( line.startsWith("# RGA Fa18") ) {
    inBlock = false
    return
  }

  if (!inBlock || line.startsWith("#") || line.isEmpty()) return

  // parse the CSV
  def parts = line.split(",")
  int runnum    = parts[0].toInteger()
  String pol    = parts[4]
  String polErr = parts[5]

  // new QADB instance for this run
  QADB qa = new QADB("latest", runnum, runnum)
  [ 'TotalOutlier','TerminalOutlier','MarginalOutlier','SectorLoss',
    'Misc','ChargeHigh','ChargeNegative','ChargeUnknown','PossiblyNoBeam'
  ].each { qa.checkForDefect(it) }
  allowMiscRuns.each { qa.allowMiscBit(it) }

  // snapshot before this run
  double beforeTot = qa.getAccumulatedCharge()
  double beforeP   = qa.getAccumulatedChargeHL(1)
  double beforeM   = qa.getAccumulatedChargeHL(-1)

  // loop over files for this run
  qa.getQaTree().each { runKey, runTree ->
    runTree.each { fileKey, fileTree ->
      int evnum = fileTree['evnumMin']
      if (qa.pass(runnum, evnum)) {
        qa.accumulateCharge()
        qa.accumulateChargeHL()
      }
    }
  }

  // compute per‐run charges
  double runTot = qa.getAccumulatedCharge()     - beforeTot
  double runP   = qa.getAccumulatedChargeHL(1)  - beforeP
  double runM   = qa.getAccumulatedChargeHL(-1) - beforeM

  // print with four‐decimal rounding
  println "${runnum}," +
          "${fmt4(runTot)}," +
          "${fmt4(runP)}," +
          "${fmt4(runM)}," +
          "${pol},${polErr}"
}