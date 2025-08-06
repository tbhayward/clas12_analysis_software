#!/usr/bin/env groovy

import org.jlab.io.hipo.HipoDataSource
import clasqa.QADB

// path to your CSV
def csvPath = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

// list of runs to allow Misc bit (same as before)
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

// tell the user what we’re doing
println "Reading CSV: $csvPath"
println "Will process runs under all RGC Su22/Fa22/Sp23 headings, stopping at # RGA Fa18 Inb Supplemental"
println "Output columns: runnum,totalCharge(nC),posCharge(nC),negCharge(nC),pol,polErr"
println ""

// control flag
boolean inBlock = false

new File(csvPath).eachLine { line ->
  line = line.trim()
  // start processing when we hit any RGC Su22/Fa22/Sp23 header
  if ( line.startsWith("# RGC Su22") ||
       line.startsWith("# RGC Fa22") ||
       line.startsWith("# RGC Sp23") ) {
    inBlock = true
    return
  }
  // stop when we hit the next header
  if ( line.startsWith("# RGA Fa18") ) {
    inBlock = false
    return
  }
  if (!inBlock) return

  // skip comment lines
  if (line.startsWith("#") || line.isEmpty()) return

  // parse CSV fields
  def parts = line.split(",")
  int runnum    = parts[0].toInteger()
  String pol    = parts[4]
  String polErr = parts[5]

  // build a new QADB for this single run
  QADB qa = new QADB("latest", runnum, runnum)
  // apply your defects
  [ 'TotalOutlier','TerminalOutlier','MarginalOutlier','SectorLoss',
    'Misc','ChargeHigh','ChargeNegative','ChargeUnknown','PossiblyNoBeam'
  ].each { qa.checkForDefect(it) }
  allowMiscRuns.each { qa.allowMiscBit(it) }

  // snapshot before we start this run’s files
  double beforeTot = qa.getAccumulatedCharge()
  double beforeP   = qa.getAccumulatedChargeHL(1)
  double beforeM   = qa.getAccumulatedChargeHL(-1)

  // loop all files in this run
  qa.getQaTree().each { rStr, runTree ->
    // there’s only one entry, but we loop anyway
    runTree.each { fStr, fileTree ->
      int evnum = fileTree['evnumMin']
      if (qa.pass(runnum, evnum)) {
        qa.accumulateCharge()
        qa.accumulateChargeHL()
      }
    }
  }

  // compute per-run charge
  double runTot = qa.getAccumulatedCharge() - beforeTot
  double runP   = qa.getAccumulatedChargeHL(1) - beforeP
  double runM   = qa.getAccumulatedChargeHL(-1) - beforeM

  // print CSV line
  println "${runnum},${runTot},${runP},${runM},${pol},${polErr}"
}