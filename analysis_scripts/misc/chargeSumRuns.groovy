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

// build QADB once (no run range)
QADB qa = new QADB("latest")
// enable QA cuts
[
  'TotalOutlier','TerminalOutlier','MarginalOutlier','SectorLoss',
  'LowLiveTime','Misc','ChargeHigh','ChargeNegative',
  'ChargeUnknown','PossiblyNoBeam'
].each { qa.checkForDefect(it) }
// allow specified runs to bypass 'Misc'
allowMiscRuns.each { qa.allowMiscBit(it) }

// formatting helper: four decimal places
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

  // detect start of RGC Su22/Fa22/Sp23 block
  if (line.startsWith("# RGC Su22") ||
      line.startsWith("# RGC Fa22") ||
      line.startsWith("# RGC Sp23")) {
    inBlock = true
    currentHeader = line.substring(1).trim()
    println "# ${currentHeader}"
    return
  }
  // detect end (RGA Fa18)
  if (line.startsWith("# RGA Fa18")) {
    inBlock = false
    return
  }
  if (!inBlock || line.startsWith("#") || line.isEmpty()) return

  // parse CSV columns
  def parts = line.split(",")
  int runnum    = parts[0].toInteger()
  String pol    = parts[4]
  String polErr = parts[5]

  // snapshot counters before this run
  double beforeTot = qa.getAccumulatedCharge()
  double beforeP   = qa.getAccumulatedChargeHL(1)
  double beforeM   = qa.getAccumulatedChargeHL(-1)

  // get that run's file tree
  def runTree = qa.getQaTree().get(runnum.toString())
  if (runTree) {
    runTree.each { fileKey, fileTree ->
      int evnum = fileTree['evnumMin']
      if (qa.pass(runnum, evnum)) {
        qa.accumulateCharge()
        qa.accumulateChargeHL()
      }
    }
  } else {
    System.err.println("Warning: no QADB entries for run ${runnum}")
  }

  // compute per-run charge
  double runTot = qa.getAccumulatedCharge()     - beforeTot
  double runP   = qa.getAccumulatedChargeHL(1)  - beforeP
  double runM   = qa.getAccumulatedChargeHL(-1) - beforeM

  // print CSV line with four‐decimal rounding
  println "${runnum}," +
          "${fmt4(runTot)}," +
          "${fmt4(runP)}," +
          "${fmt4(runM)}," +
          "${pol},${polErr}"
}