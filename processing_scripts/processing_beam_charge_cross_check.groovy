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
        // mirror your colleague’s defect cuts:
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

            // ————— process each event once ————— //
            while (reader.hasEvent()) {
                HipoDataEvent event = reader.getNextEvent()
                eventCount++
                if ((eventCount % 1_000_000) == 0) {
                    print "  processed ${eventCount} events…\r"
                }

                int runnum = event.getBank("RUN::config").getInt("run", 0)
                // skip known bleed-through
                if (runnum > 16600 && runnum < 16700) {
                    println "\n  Bleed-through run ${runnum} → stopping file."
                    // break
                }
                int evnum = event.getBank("RUN::config").getInt("event", 0)

                // tell QADB about this event
                qa.query(runnum, evnum)

                // only if QA cuts pass AND there is a real beam charge:
                if (qa.pass(runnum, evnum)) {
                    // attempt to fetch the “beam charge” for this helicity;
                    // if QADB had no data you’d get null, so guard it:
                    Double bc = qa.getBeamChargeHL()
                    if (bc != null) {
                        qa.accumulateChargeHL()
                    }
                }
            }

            // ————— report per-helicity totals ————— //
            println "\nTotals for file ${hipofile.name}:"
            HLstate.each { hel ->
                Double tot = qa.getAccumulatedChargeHL(hel) ?: 0.0
                println "  HL charge(${hel}) = ${tot}"
            }

            reader.close()
        }

        println "\nAll files processed."
    }
}