#!/usr/bin/env groovy
/*
 * author Timothy B. Hayward
 */

import java.io.File
import org.jlab.io.hipo.HipoDataSource
import org.jlab.io.hipo.HipoDataEvent
import groovy.io.FileType
import clasqa.QADB

public class processing_beamCharge {

    public static void main(String[] args) {
        // ~~~~~~~~~~~~~~~~ check input directory ~~~~~~~~~~~~~~~~ //
        if (!args) {
            System.err.println("ERROR: Please enter a HIPO file directory as the first argument")
            System.exit(1)
        }

        // ~~~~~~~~~~~~~~~~ gather .hipo files ~~~~~~~~~~~~~~~~ //
        def hipoList = []
        new File(args[0]).eachFileRecurse(FileType.FILES) { file ->
            if (file.name.endsWith('.hipo')) {
                hipoList << file
            }
        }
        if (hipoList.isEmpty()) {
            System.err.println("ERROR: No .hipo files found in directory ${args[0]}")
            System.exit(1)
        }

        // ~~~~~~~~~~~~~~~~ limit to requested number of files ~~~~~~~~~~~~~~~~ //
        int nFiles = args.length >= 2 ? Integer.parseInt(args[1]) : hipoList.size()
        if (nFiles <= 0 || nFiles > hipoList.size()) {
            println("WARNING: Invalid file count; processing all ${hipoList.size()} files")
            nFiles = hipoList.size()
        }

        // helicity states to report
        def HLstate = [-1, 0, 1]

        // initialize QA database
        QADB qa = new QADB("latest")

        // ~~~~~~~~~~~~~~~~ loop over each HIPO file ~~~~~~~~~~~~~~~~ //
        for (int fileIndex = 0; fileIndex < nFiles; fileIndex++) {
            File hipoFile = hipoList[fileIndex]
            println("\nOpening file ${fileIndex + 1} of ${nFiles}: ${hipoFile}\n")

            HipoDataSource reader = new HipoDataSource()
            reader.open(hipoFile)

            long numEvents = 0

            // process every event exactly once
            while (reader.hasEvent()) {
                HipoDataEvent event = reader.getNextEvent()
                numEvents++
                if (numEvents % 1_000_000 == 0) {
                    print("Processed ${numEvents} events... ")
                }

                int runnum = event.getBank("RUN::config").getInt("run", 0)
                // skip Hall C bleed-through
                if (runnum > 16600 && runnum < 16700) {
                    println("\nBleedthrough run ${runnum} encountered; stopping this file.")
                    break
                }
                int evnum = event.getBank("RUN::config").getInt("event", 0)

                // update QA and accumulate charge by helicity
                qa.query(runnum, evnum)
                qa.accumulateChargeHL()
            }

            // report final totals per helicity
            HLstate.each { hel ->
                Double charge = qa.getAccumulatedChargeHL(hel) ?: 0.0
                println("HL charge(${hel}) = ${charge}")
            }

            // reset for next file
            qa.resetAccumulatedChargeHL()
            reader.close()
        }

        println("\nProcessing complete.")
    }
}