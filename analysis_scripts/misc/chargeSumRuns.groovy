#!/usr/bin/env groovy

import org.jlab.io.hipo.HipoDataSource
import clasqa.QADB

// path to your CSV
def csvPath = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

// full list of runs to allow the 'Misc' QA bit
def allowMiscRuns = [
     5046, 5047, 5051, 5128, 5129, 5130, 5158, 5159,
     5160, 5163, 5165, 5166, 5167, 5168, 5169, 5180,
     5181, 5182, 5183, 5400, 5448, 5495, 5496, 5505,
     5567, 5610, 5617, 5621, 5623, 6736, 6737, 6738,
     6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
     6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
     6757,
     16194, 16089, 16185, 16308, 16184, 16307, 16309,
     16872, 16975,
     17763, 17764, 17765, 17766, 17767, 17768
]

// QA defect list
def defectList = [
  'TotalOutlier','TerminalOutlier','MarginalOutlier','SectorLoss',
  'LowLiveTime','Misc','ChargeHigh','ChargeNegative',
  'ChargeUnknown','PossiblyNoBeam'
]

// helper to format 4 decimals
def fmt4 = { double x -> String.format("%.4f", x) }

// print header
println "Reading CSV: $csvPath"
println "Output columns:"
println "  runnum, totalCharge(nC), posCharge(nC), negCharge(nC), pol, polErr"
println ""

// control flags
boolean inBlock = false
String currentHeader = ""

// process CSV
new File(csvPath).eachLine { raw ->
  String line = raw.trim()

  // start of a relevant block
  if (line.startsWith("# RGC Su22") ||
      line.startsWith("# RGC Fa22") ||
      line.startsWith("# RGC Sp23")) {
    inBlock = true
    currentHeader = line.substring(1).trim()
    println "# ${currentHeader}"
    return
  }
  // end marker
  if (line.startsWith("# RGA Fa18")) {
    inBlock = false
    return
  }
  if (!inBlock || line.startsWith("#") || line.isEmpty()) return

  // parse CSV fields
  def parts  = line.split(",")
  int runnum = parts[0].toInteger()
  String pol    = parts[4]
  String polErr = parts[5]

  // instantiate QADB for exactly this run
  QADB qa = new QADB("latest", runnum, runnum)
  defectList.each { qa.checkForDefect(it) }
  allowMiscRuns.each { qa.allowMiscBit(it) }

  // snapshot before accumulation
  double beforeTot = qa.getAccumulatedCharge()
  double beforeP   = qa.getAccumulatedChargeHL(1)
  double beforeM   = qa.getAccumulatedChargeHL(-1)

  // accumulate over this run's files
  def runTree = qa.getQaTree().get(runnum.toString())
  if (runTree) {
    runTree.each { _, fileTree ->
      int evnum = fileTree['evnumMin']
      if (qa.pass(runnum, evnum)) {
        qa.accumulateCharge()
        qa.accumulateChargeHL()
      }
    }
  } else {
    System.err.println("Warning: no QADB entries for run ${runnum}")
  }

  // compute per‐run charges
  double runTot = qa.getAccumulatedCharge()     - beforeTot
  double runP   = qa.getAccumulatedChargeHL(1)  - beforeP
  double runM   = qa.getAccumulatedChargeHL(-1) - beforeM

  // print with 4-decimal precision
  println "${runnum}," +
          "${fmt4(runTot)}," +
          "${fmt4(runP)}," +
          "${fmt4(runM)}," +
          "${pol},${polErr}"
}