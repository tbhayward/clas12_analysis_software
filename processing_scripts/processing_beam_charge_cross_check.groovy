#!/usr/bin/env groovy
/*
 * author Timothy B. Hayward
 * 
 * For each HIPO file in the given directory, accumulates beam charge
 * by helicity state via QADB and then prints:
 *   runnum,total,posHel,negHel,noHel,0,0
 */

import java.io.File
import groovy.io.FileType
import org.jlab.io.hipo.HipoDataSource
import org.jlab.io.hipo.HipoDataEvent
import clasqa.QADB

// --- parse arguments ---
if (this.args.length < 1) {
    System.err.println("ERROR: Please specify a directory containing .hipo files")
    System.exit(1)
}
File dir = new File(this.args[0])
if (!dir.isDirectory()) {
    System.err.println("ERROR: ${dir} is not a directory")
    System.exit(1)
}

// gather all .hipo files
def hipoList = []
dir.eachFileRecurse(FileType.FILES) { f ->
    if (f.name.toLowerCase().endsWith(".hipo")) {
        hipoList << f
    }
}
if (hipoList.isEmpty()) {
    System.err.println("ERROR: No .hipo files found in ${dir}")
    System.exit(1)
}

// optional second arg: max files to process
int nFiles = (this.args.length >= 2) 
    ? this.args[1].toInteger() 
    : hipoList.size()
if (nFiles <= 0 || nFiles > hipoList.size()) {
    println "WARNING: invalid max count; processing all ${hipoList.size()} files"
    nFiles = hipoList.size()
}

// helicity states in QADB
def HLstates = [-1, 0, 1]

// instantiate QADB
QADB qa = new QADB("latest")

// --- process each file ---
for (int i = 0;  i < nFiles;  i++) {
    File hipofile = hipoList[i]
    println "\n=== File ${i+1} of ${nFiles}: ${hipofile.name} ==="

    // reset per-file accumulators in QADB
    qa.resetAccumulatedChargeHL()

    HipoDataSource reader = new HipoDataSource()
    reader.open(hipofile)

    long eventCount = 0
    int fileRun = -1

    // loop over events exactly once
    while (reader.hasEvent()) {
        HipoDataEvent event = reader.getNextEvent()
        eventCount++
        if ((eventCount % 1_000_000) == 0) {
            print "  processed ${eventCount} events...\r"
        }

        // extract run & event
        int runnum = event.getBank("RUN::config").getInt("run", 0)
        int evnum   = event.getBank("RUN::config").getInt("event", 0)
        if (fileRun < 0) {
            fileRun = runnum
        }

        // feed into QADB
        qa.query(runnum, evnum)
        try {
            qa.accumulateChargeHL()
        } catch (GroovyRuntimeException e) {
            // happens if QADB has no valid charge → skip
        }
    }

    reader.close()

    // gather per-helicity totals
    double negSum  = qa.getAccumulatedChargeHL(-1) ?: 0.0
    double zeroSum = qa.getAccumulatedChargeHL(0)  ?: 0.0
    double posSum  = qa.getAccumulatedChargeHL(1)  ?: 0.0
    double total   = negSum + zeroSum + posSum

    // print CSV line
    println "${fileRun},${total},${posSum},${negSum},${zeroSum},0,0"
}

println "\nAll done."