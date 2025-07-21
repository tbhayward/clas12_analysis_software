#!/usr/bin/env groovy
/*
 * author Timothy B. Hayward
 *
 * For each HIPO file in the given directory, accumulates beam charge
 * by helicity state via QADB and then prints out the full CSV summary
 * up to and including that file:
 *
 *   runnum,total,posHel,negHel,noHel,0,0
 */

import java.io.File
import groovy.io.FileType
import org.jlab.io.hipo.HipoDataSource
import org.jlab.io.hipo.HipoDataEvent
import clasqa.QADB

class processing_beamCharge {

    static void main(String[] args) {
        if (!args) {
            System.err.println("ERROR: Please specify a directory containing .hipo files")
            System.exit(1)
        }
        File dir = new File(args[0])
        if (!dir.isDirectory()) {
            System.err.println("ERROR: ${dir} is not a directory")
            System.exit(1)
        }

        // collect all .hipo files
        def hipoList = []
        dir.eachFileRecurse(FileType.FILES) { f ->
            if (f.name.toLowerCase().endsWith(".hipo")) hipoList << f
        }
        if (hipoList.isEmpty()) {
            System.err.println("ERROR: No .hipo files found in ${dir}")
            System.exit(1)
        }

        // optional cap on number of files
        int nFiles = (args.length >= 2) ? args[1].toInteger() : hipoList.size()
        if (nFiles <= 0 || nFiles > hipoList.size()) {
            println "WARNING: invalid max count; processing all ${hipoList.size()} files"
            nFiles = hipoList.size()
        }

        // instantiate QADB once
        QADB qa = new QADB("latest")

        // we’ll accumulate one CSV line per file here
        String beamChargeList = ""

        double negSum  = 0;
        double zeroSum = 0;
        double posSum  = 0;
        double total   = negSum + zeroSum + posSum

        // process each file in turn
        for (int i = 0; i < nFiles; i++) {
            File hipofile = hipoList[i]
            println "\n=== File ${i+1}/${nFiles}: ${hipofile.name} ==="

            // reset QADB’s internal charge accumulators
            qa.resetAccumulatedChargeHL()

            HipoDataSource reader = new HipoDataSource()
            reader.open(hipofile)

            long eventCount = 0
            int fileRun = -1    // <— will set once, on first real runnum
            runnum = fileRun;
            if (runnum > 17768) process_event = false; // outbending RGC Sp23
            if (runnum == 17331 || runnum == 16987 || runnum == 17079 || runnum == 17190 || runnum == 17639) process_event = false; // low live time
            if (runnum == 16850 || runnum == 16851 || runnum == 16852 || runnum == 16855 || runnum == 16879) process_event = false; // luminosity scans

            while (reader.hasEvent()) {
                HipoDataEvent event = reader.getNextEvent()
                eventCount++
                if ((eventCount % 1000000) == 0) {
                    print "  processed ${eventCount} events...\r"
                }

                // pull run & event from the bank (defaults to 0 if missing)
                int runnum = event.getBank("RUN::config").getInt("run", 0)
                int evnum   = event.getBank("RUN::config").getInt("event",0)

                // only set fileRun once, and only if it’s non-zero
                if (fileRun < 0 && runnum != 0) {
                    fileRun = runnum
                }

                // feed this event to QADB
                qa.query(runnum, evnum)
                try {
                    qa.accumulateChargeHL()
                } catch (GroovyRuntimeException e) {
                    // skip events where QADB had no valid charge
                }

                negSum  = qa.getAccumulatedChargeHL(-1) ?: 0.0
	            zeroSum = qa.getAccumulatedChargeHL(0)  ?: 0.0
	            posSum  = qa.getAccumulatedChargeHL(1)  ?: 0.0
	            total   = negSum + zeroSum + posSum
            }

            reader.close()

            // pull out the three helicity sums
            // double negSum  = qa.getAccumulatedChargeHL(-1) ?: 0.0
            // double zeroSum = qa.getAccumulatedChargeHL(0)  ?: 0.0
            // double posSum  = qa.getAccumulatedChargeHL(1)  ?: 0.0
            // double total   = negSum + zeroSum + posSum
            negSum  = qa.getAccumulatedChargeHL(-1) ?: 0.0
            zeroSum = qa.getAccumulatedChargeHL(0)  ?: 0.0
            posSum  = qa.getAccumulatedChargeHL(1)  ?: 0.0
            total   = negSum + zeroSum + posSum

            // if we never saw a non-zero run, fall back to 0
            if (fileRun < 0) {
                fileRun = 0
            }

            // format each float to 3 decimal places
            String line = String.format(
                "%d,%.3f,%.3f,%.3f,%.3f,0,0%n",
                fileRun,
                total,
                posSum,
                negSum,
                zeroSum
            )

            // append the new line
            beamChargeList += line

            // print *all* accumulated lines so far
            println "\n--- Accumulated summary so far ---"
            print beamChargeList
        }

        println "\nAll files processed."
    }
}