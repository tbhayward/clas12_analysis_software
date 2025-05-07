#!/usr/bin/env groovy
/*
 * author Timothy B. Hayward
 */

import java.io.File
import org.jlab.io.hipo.HipoDataSource
import org.jlab.io.hipo.HipoDataEvent
import groovy.io.FileType
import clasqa.QADB

class processing_beamCharge {

    static void main(String[] args) {
        // ————— check args ————— //
        if (!args) {
            System.err.println("ERROR: Please enter a HIPO file directory as the first argument")
            System.exit(1)
        }

        // ————— collect .hipo files ————— //
        def hipoList = []
        new File(args[0]).eachFileRecurse(FileType.FILES) { f ->
            if (f.name.endsWith('.hipo')) hipoList << f
        }
        if (hipoList.isEmpty()) {
            System.err.println("ERROR: No .hipo files found in ${args[0]}")
            System.exit(1)
        }

        // ————— limit file count ————— //
        int nFiles = (args.length >= 2) ? args[1].toInteger() : hipoList.size()
        if (nFiles <= 0 || nFiles > hipoList.size()) {
            println("WARNING: Invalid file count; defaulting to all ${hipoList.size()} files")
            nFiles = hipoList.size()
        }

        // ————— prepare QA ————— //
        QADB qa = new QADB("latest")
        // apply the same defect cuts as your colleague
        qa.checkForDefect("SectorLoss")
        qa.checkForDefect("MarginalOutlier")

        def HLstate = [-1, 0, 1]

        // ————— loop over files ————— //
        for (int idx = 0; idx < nFiles; idx++) {
            File hipofile = hipoList[idx]
            println "\nOpening file ${idx+1}/${nFiles}: ${hipofile.name}\n"

            HipoDataSource reader = new HipoDataSource()
            reader.open(hipofile)

            qa.resetAccumulatedChargeHL()
            long eventCount = 0

            // ————— process each event exactly once ————— //
            while (reader.hasEvent()) {
                HipoDataEvent event = reader.getNextEvent()
                eventCount++
                if ((eventCount % 1_000_000) == 0) {
                    print "  processed ${eventCount} events…\r"
                }

                int runnum = event.getBank("RUN::config").getInt("run",  0)
                // skip Hall C bleed-through
                if (runnum > 16600 && runnum < 16700) {
                    println "\n  Bleed-through run ${runnum}; stopping this file."
                    break
                }
                int evnum = event.getBank("RUN::config").getInt("event",0)

                // update QA
                qa.query(runnum, evnum)

                // only accumulate if QA cuts pass
                if (qa.pass(runnum, evnum)) {
                    try {
                        qa.accumulateChargeHL()
                    } catch (GroovyRuntimeException e) {
                        // no valid charge for this event → skip it
                    }
                }
            }

            // ————— report per-helicity totals ————— //
            println "\nTotals for ${hipofile.name}:"
            HLstate.each { hel ->
                // if nothing accumulated, default to 0.0
                Double tot = qa.getAccumulatedChargeHL(hel) ?: 0.0
                println "  HL charge(${hel}) = ${tot}"
            }

            reader.close()
        }

        println "\nAll files processed."
    }
}