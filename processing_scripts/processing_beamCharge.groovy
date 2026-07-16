/*
 * author Timothy B. Hayward
 *
 * Extract RUN::scaler / HEL::scaler charge.
 *
 * Normal operation:
 *   - Reads tag-1 events.
 *   - Calculates RUN::scaler charge from:
 *
 *         max(fcupgated) - min(fcupgated)
 *
 *   - Sums HEL::scaler fcupgated values separately for helicity
 *     +1, -1 and 0.
 *
 * Diagnostic operation:
 *   - Enable with the optional --diagnostic argument.
 *   - The normal output is still calculated from tag-1 events.
 *   - A second unrestricted scan is performed over all event tags.
 *   - At the end, the script compares tag-1 and all-tag RUN::scaler
 *     coverage for every run.
 *
 * Usage:
 *
 *   run-groovy processing_beamCharge.groovy \
 *       <hipo-directory> [number-of-files] [train|dst] [--diagnostic]
 *
 * Examples:
 *
 *   run-groovy processing_beamCharge.groovy /path/to/train -1 train
 *
 *   run-groovy processing_beamCharge.groovy \
 *       /path/to/train -1 train --diagnostic
 *
 *   run-groovy processing_beamCharge.groovy \
 *       /path/to/dst -1 dst --diagnostic
 */

import java.io.File;
import java.util.Locale;

import groovy.io.FileType;

import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.SchemaFactory;


public class processing_beamCharge {

    /*
     * Accumulates the normal tag-1 charge results for one run.
     */
    static class ChargeAccumulator {

        int runnum = 0;

        double runScalerChargeMin = Double.POSITIVE_INFINITY;
        double runScalerChargeMax = Double.NEGATIVE_INFINITY;
        boolean foundRunScalerCharge = false;

        double posHelbeamChargeTotal = 0.0d;
        double negHelbeamChargeTotal = 0.0d;
        double noHelbeamChargeTotal = 0.0d;

        int nFiles = 0;

        double getRunScalerCharge() {

            if (!foundRunScalerCharge) {
                return 0.0d;
            }

            return runScalerChargeMax - runScalerChargeMin;
        }
    }


    /*
     * Describes one RUN::scaler endpoint.
     *
     * An endpoint is either the minimum or maximum fcupgated value observed
     * within a selected event-tag stream.
     */
    static class ScalerEndpoint {

        boolean valid = false;

        double charge = Double.NaN;

        int runnum = 0;
        int tag = Integer.MIN_VALUE;

        int runConfigEvent = -1;
        long timestamp = -1L;

        String filename = "";
        int fileNumber = -1;

        int scalerRow = -1;

        long unrestrictedEventIndex = -1L;
        long scalerRecordIndex = -1L;

        String formatLocation() {

            if (!valid) {
                return "not found";
            }

            String timestampText =
                timestamp >= 0L ? Long.toString(timestamp) : "unavailable";

            String eventText =
                runConfigEvent >= 0 ?
                    Integer.toString(runConfigEvent) :
                    "unavailable";

            return String.format(
                Locale.US,
                "charge=%.9f nC, tag=%d, RUN::config.event=%s, " +
                "timestamp=%s, scaler row=%d, file #%d, " +
                "all-tag event index=%d, scaler-record index=%d%n" +
                "        file=%s",
                charge,
                tag,
                eventText,
                timestampText,
                scalerRow,
                fileNumber,
                unrestrictedEventIndex,
                scalerRecordIndex,
                filename
            );
        }
    }


    /*
     * RUN::scaler statistics for one selected event stream.
     *
     * This can represent:
     *
     *   - all event tags
     *   - tag 1 only
     *   - one specific tag
     */
    static class ScalerStreamAccumulator {

        String label = "";

        long scalerEventCount = 0L;
        long scalerRowCount = 0L;

        boolean foundScaler = false;

        double minimum = Double.POSITIVE_INFINITY;
        double maximum = Double.NEGATIVE_INFINITY;

        ScalerEndpoint minimumEndpoint = new ScalerEndpoint();
        ScalerEndpoint maximumEndpoint = new ScalerEndpoint();

        double getCharge() {

            if (!foundScaler) {
                return 0.0d;
            }

            return maximum - minimum;
        }

        void update(
            double charge,
            int runnum,
            int tag,
            int runConfigEvent,
            long timestamp,
            String filename,
            int fileNumber,
            int scalerRow,
            long unrestrictedEventIndex,
            long scalerRecordIndex
        ) {

            scalerRowCount++;
            foundScaler = true;

            if (charge < minimum) {

                minimum = charge;

                minimumEndpoint = makeEndpoint(
                    charge,
                    runnum,
                    tag,
                    runConfigEvent,
                    timestamp,
                    filename,
                    fileNumber,
                    scalerRow,
                    unrestrictedEventIndex,
                    scalerRecordIndex
                );
            }

            if (charge > maximum) {

                maximum = charge;

                maximumEndpoint = makeEndpoint(
                    charge,
                    runnum,
                    tag,
                    runConfigEvent,
                    timestamp,
                    filename,
                    fileNumber,
                    scalerRow,
                    unrestrictedEventIndex,
                    scalerRecordIndex
                );
            }
        }

        static ScalerEndpoint makeEndpoint(
            double charge,
            int runnum,
            int tag,
            int runConfigEvent,
            long timestamp,
            String filename,
            int fileNumber,
            int scalerRow,
            long unrestrictedEventIndex,
            long scalerRecordIndex
        ) {

            ScalerEndpoint endpoint = new ScalerEndpoint();

            endpoint.valid = true;
            endpoint.charge = charge;

            endpoint.runnum = runnum;
            endpoint.tag = tag;

            endpoint.runConfigEvent = runConfigEvent;
            endpoint.timestamp = timestamp;

            endpoint.filename = filename;
            endpoint.fileNumber = fileNumber;

            endpoint.scalerRow = scalerRow;

            endpoint.unrestrictedEventIndex = unrestrictedEventIndex;
            endpoint.scalerRecordIndex = scalerRecordIndex;

            return endpoint;
        }
    }


    /*
     * Holds all diagnostic information for one run.
     */
    static class DiagnosticRunAccumulator {

        int runnum = 0;

        /*
         * All RUN::scaler rows seen in the unrestricted scan.
         */
        ScalerStreamAccumulator allTags =
            new ScalerStreamAccumulator(label: "all tags");

        /*
         * Only RUN::scaler rows whose event.getEventTag() is 1.
         */
        ScalerStreamAccumulator tag1 =
            new ScalerStreamAccumulator(label: "tag 1");

        /*
         * Per-tag breakdown.
         */
        LinkedHashMap<Integer, ScalerStreamAccumulator> byTag =
            new LinkedHashMap<Integer, ScalerStreamAccumulator>();

        ScalerStreamAccumulator getTagAccumulator(int tag) {

            ScalerStreamAccumulator accumulator = byTag.get(tag);

            if (accumulator == null) {

                accumulator = new ScalerStreamAccumulator(
                    label: "tag " + tag
                );

                byTag.put(tag, accumulator);
            }

            return accumulator;
        }
    }


    /*
     * Parsed command-line configuration.
     */
    static class ProgramOptions {

        String inputDirectory = "";

        int requestedFiles = -1;
        String processingMode = "train";

        boolean diagnostic = false;
    }


    public static void main(String[] args) {

        ProgramOptions options = parseArguments(args);

        def hipoList = [];

        File inputDirectory = new File(options.inputDirectory);

        if (!inputDirectory.exists()) {

            println(
                "ERROR: Input directory does not exist: " +
                inputDirectory.getAbsolutePath()
            );

            System.exit(1);
        }

        if (!inputDirectory.isDirectory()) {

            println(
                "ERROR: Input path is not a directory: " +
                inputDirectory.getAbsolutePath()
            );

            System.exit(1);
        }

        inputDirectory.eachFileRecurse(FileType.FILES) {

            if (it.name.endsWith(".hipo")) {
                hipoList << it;
            }
        }

        hipoList = hipoList.sort {
            it.getAbsolutePath()
        };

        if (hipoList.isEmpty()) {

            println(
                "ERROR: No .hipo files were found under: " +
                inputDirectory.getAbsolutePath()
            );

            System.exit(1);
        }

        int nFiles = determineNumberOfFiles(
            hipoList,
            options.requestedFiles
        );

        int hipoTag = 1;

        println();
        println("Input directory = " + inputDirectory.getAbsolutePath());
        println("Processing mode = " + options.processingMode);
        println("Normal charge stream = tag-" + hipoTag);
        println("Number of HIPO files to process = " + nFiles);
        println("Diagnostic mode = " + options.diagnostic);
        println();

        if (options.processingMode.equals("dst")) {

            processDstMode(
                hipoList,
                nFiles,
                hipoTag
            );

        } else {

            processTrainMode(
                hipoList,
                nFiles,
                hipoTag
            );
        }

        if (options.diagnostic) {

            println();
            println(
                "============================================================"
            );
            println(
                "Beginning unrestricted RUN::scaler diagnostic scan"
            );
            println(
                "============================================================"
            );
            println();

            LinkedHashMap<Integer, DiagnosticRunAccumulator> diagnosticMap =
                runScalerDiagnosticScan(
                    hipoList,
                    nFiles
                );

            printDiagnosticReport(diagnosticMap);
        }
    }


    static ProgramOptions parseArguments(String[] args) {

        if (args == null || args.length < 1) {

            printUsage();
            System.exit(1);
        }

        ProgramOptions options = new ProgramOptions();

        List<String> positionalArguments = new ArrayList<String>();

        for (String argument : args) {

            if (argument.equalsIgnoreCase("--diagnostic")) {

                options.diagnostic = true;

            } else {

                positionalArguments.add(argument);
            }
        }

        if (positionalArguments.isEmpty()) {

            printUsage();
            System.exit(1);
        }

        options.inputDirectory = positionalArguments.get(0);

        if (positionalArguments.size() >= 2) {

            try {

                options.requestedFiles =
                    Integer.parseInt(positionalArguments.get(1));

            } catch (NumberFormatException exception) {

                println(
                    "ERROR: Invalid number-of-files argument: " +
                    positionalArguments.get(1)
                );

                printUsage();
                System.exit(1);
            }
        }

        if (positionalArguments.size() >= 3) {

            options.processingMode =
                positionalArguments.get(2).toLowerCase();
        }

        if (
            !options.processingMode.equals("train") &&
            !options.processingMode.equals("dst")
        ) {

            println(
                "ERROR: Unknown processing mode: " +
                options.processingMode
            );

            println("Allowed modes are:");
            println("  train");
            println("  dst");

            System.exit(1);
        }

        if (positionalArguments.size() > 3) {

            println(
                "ERROR: Too many positional command-line arguments."
            );

            printUsage();
            System.exit(1);
        }

        return options;
    }


    static void printUsage() {

        println();
        println(
            "Usage:"
        );

        println(
            "  run-groovy processing_beamCharge.groovy " +
            "<hipo-directory> [number-of-files] [train|dst] " +
            "[--diagnostic]"
        );

        println();
        println(
            "Examples:"
        );

        println(
            "  run-groovy processing_beamCharge.groovy " +
            "/path/to/files -1 train"
        );

        println(
            "  run-groovy processing_beamCharge.groovy " +
            "/path/to/files -1 train --diagnostic"
        );

        println(
            "  run-groovy processing_beamCharge.groovy " +
            "/path/to/files -1 dst --diagnostic"
        );

        println();
    }


    static int determineNumberOfFiles(
        def hipoList,
        int requestedFiles
    ) {

        int nFiles = hipoList.size();

        if (requestedFiles > 0) {

            if (requestedFiles <= hipoList.size()) {

                nFiles = requestedFiles;

            } else {

                println(
                    "WARNING: Requested number of files exceeds the " +
                    "number found."
                );

                println(
                    "Using all " + hipoList.size() + " files."
                );
            }

        } else {

            println(
                "Using all " + hipoList.size() + " files."
            );
        }

        return nFiles;
    }


    static void processTrainMode(
        def hipoList,
        int nFiles,
        int hipoTag
    ) {

        long numEvents = 0L;
        int currentFile = 0;

        String beamChargeList = "";

        while (currentFile < nFiles) {

            File inputFile = hipoList[currentFile];

            println();
            println(
                "Opening file " +
                Integer.toString(currentFile + 1) +
                " out of " +
                nFiles
            );

            println(
                "File: " +
                inputFile.getAbsolutePath()
            );

            println();

            double beamChargeMax = 0.0d;

            double runScalerChargeMin = Double.POSITIVE_INFINITY;
            double runScalerChargeMax = Double.NEGATIVE_INFINITY;
            double runScalerCharge = 0.0d;

            boolean foundRunScalerCharge = false;

            double posHelbeamChargeTotal = 0.0d;
            double negHelbeamChargeTotal = 0.0d;
            double noHelbeamChargeTotal = 0.0d;

            int runnum = 0;

            HipoReader reader = new HipoReader();

            /*
             * Normal operation remains restricted to tag 1.
             */
            reader.setTags(hipoTag);
            reader.open(inputFile.getAbsolutePath());

            SchemaFactory schema = reader.getSchemaFactory();

            Bank runConfigBank = new Bank(
                schema.getSchema("RUN::config")
            );

            Bank runScalerBank = new Bank(
                schema.getSchema("RUN::scaler")
            );

            Bank helScalerBank = new Bank(
                schema.getSchema("HEL::scaler")
            );

            Event event = new Event();

            while (reader.hasNext()) {

                reader.nextEvent(event);
                numEvents++;

                event.read(runConfigBank);
                event.read(runScalerBank);
                event.read(helScalerBank);

                if (
                    runConfigBank.getRows() > 0 &&
                    runnum == 0
                ) {

                    runnum = runConfigBank.getInt("run", 0);
                }

                if (runScalerBank.getRows() > 0) {

                    /*
                     * RUN::scaler fcupgated is cumulative.
                     *
                     * Process every row, even though the usual bank has
                     * one row per scaler event.
                     */
                    for (
                        int row = 0;
                        row < runScalerBank.getRows();
                        row++
                    ) {

                        double thisRunScalerCharge =
                            (double) runScalerBank.getFloat(
                                "fcupgated",
                                row
                            );

                        if (
                            thisRunScalerCharge <
                            runScalerChargeMin
                        ) {

                            runScalerChargeMin =
                                thisRunScalerCharge;
                        }

                        if (
                            thisRunScalerCharge >
                            runScalerChargeMax
                        ) {

                            runScalerChargeMax =
                                thisRunScalerCharge;
                        }

                        if (
                            thisRunScalerCharge >
                            beamChargeMax
                        ) {

                            beamChargeMax =
                                thisRunScalerCharge;
                        }

                        foundRunScalerCharge = true;
                    }

                    runScalerCharge =
                        runScalerChargeMax -
                        runScalerChargeMin;
                }

                if (helScalerBank.getRows() > 0) {

                    for (
                        int row = 0;
                        row < helScalerBank.getRows();
                        row++
                    ) {

                        int helicity =
                            helScalerBank.getInt(
                                "helicity",
                                row
                            );

                        double helBeamCharge =
                            (double) helScalerBank.getFloat(
                                "fcupgated",
                                row
                            );

                        if (helicity == 1) {

                            posHelbeamChargeTotal +=
                                helBeamCharge;

                        } else if (helicity == -1) {

                            negHelbeamChargeTotal +=
                                helBeamCharge;

                        } else if (helicity == 0) {

                            noHelbeamChargeTotal +=
                                helBeamCharge;
                        }
                    }
                }

                if (numEvents % 1000000L == 0L) {

                    print(
                        "processed: " +
                        numEvents +
                        " tag-" +
                        hipoTag +
                        " events, "
                    );

                    print(
                        "max beamCharge of current run/file = " +
                        beamChargeMax +
                        " nC, "
                    );

                    print(
                        "RUN::scaler max-min of current run/file = " +
                        runScalerCharge +
                        " nC.\n"
                    );
                }
            }

            reader.close();

            if (!foundRunScalerCharge) {
                runScalerCharge = 0.0d;
            }

            beamChargeList +=
                runnum.toString() +
                "," +
                Double.toString(runScalerCharge) +
                ",";

            beamChargeList +=
                Double.toString(posHelbeamChargeTotal) +
                ",";

            beamChargeList +=
                Double.toString(negHelbeamChargeTotal) +
                ",";

            beamChargeList +=
                Double.toString(noHelbeamChargeTotal) +
                ",0,0\n";

            println();
            println();
            println();
            print(beamChargeList);
            println();

            currentFile++;
        }
    }


    static void processDstMode(
        def hipoList,
        int nFiles,
        int hipoTag
    ) {

        long numEvents = 0L;
        int currentFile = 0;

        LinkedHashMap<Integer, ChargeAccumulator> runChargeMap =
            new LinkedHashMap<Integer, ChargeAccumulator>();

        while (currentFile < nFiles) {

            File inputFile = hipoList[currentFile];

            println();
            println(
                "Opening file " +
                Integer.toString(currentFile + 1) +
                " out of " +
                nFiles
            );

            println(
                "File: " +
                inputFile.getAbsolutePath()
            );

            println();

            int runnumForFile = 0;

            HipoReader reader = new HipoReader();

            /*
             * Normal operation remains restricted to tag 1.
             */
            reader.setTags(hipoTag);
            reader.open(inputFile.getAbsolutePath());

            SchemaFactory schema = reader.getSchemaFactory();

            Bank runConfigBank = new Bank(
                schema.getSchema("RUN::config")
            );

            Bank runScalerBank = new Bank(
                schema.getSchema("RUN::scaler")
            );

            Bank helScalerBank = new Bank(
                schema.getSchema("HEL::scaler")
            );

            Event event = new Event();

            while (reader.hasNext()) {

                reader.nextEvent(event);
                numEvents++;

                event.read(runConfigBank);
                event.read(runScalerBank);
                event.read(helScalerBank);

                if (
                    runConfigBank.getRows() > 0 &&
                    runnumForFile == 0
                ) {

                    runnumForFile =
                        runConfigBank.getInt("run", 0);
                }

                int runnumForEvent = runnumForFile;

                if (runConfigBank.getRows() > 0) {

                    runnumForEvent =
                        runConfigBank.getInt("run", 0);
                }

                if (runnumForEvent == 0) {
                    continue;
                }

                ChargeAccumulator accumulator =
                    runChargeMap.get(runnumForEvent);

                if (accumulator == null) {

                    accumulator = new ChargeAccumulator();
                    accumulator.runnum = runnumForEvent;

                    runChargeMap.put(
                        runnumForEvent,
                        accumulator
                    );
                }

                if (runScalerBank.getRows() > 0) {

                    /*
                     * RUN::scaler is cumulative across the run.
                     *
                     * Find the global minimum and maximum across every
                     * file belonging to the run.
                     */
                    for (
                        int row = 0;
                        row < runScalerBank.getRows();
                        row++
                    ) {

                        double thisRunScalerCharge =
                            (double) runScalerBank.getFloat(
                                "fcupgated",
                                row
                            );

                        if (
                            thisRunScalerCharge <
                            accumulator.runScalerChargeMin
                        ) {

                            accumulator.runScalerChargeMin =
                                thisRunScalerCharge;
                        }

                        if (
                            thisRunScalerCharge >
                            accumulator.runScalerChargeMax
                        ) {

                            accumulator.runScalerChargeMax =
                                thisRunScalerCharge;
                        }

                        accumulator.foundRunScalerCharge = true;
                    }
                }

                if (helScalerBank.getRows() > 0) {

                    for (
                        int row = 0;
                        row < helScalerBank.getRows();
                        row++
                    ) {

                        int helicity =
                            helScalerBank.getInt(
                                "helicity",
                                row
                            );

                        double helBeamCharge =
                            (double) helScalerBank.getFloat(
                                "fcupgated",
                                row
                            );

                        if (helicity == 1) {

                            accumulator.posHelbeamChargeTotal +=
                                helBeamCharge;

                        } else if (helicity == -1) {

                            accumulator.negHelbeamChargeTotal +=
                                helBeamCharge;

                        } else if (helicity == 0) {

                            accumulator.noHelbeamChargeTotal +=
                                helBeamCharge;
                        }
                    }
                }

                if (numEvents % 1000000L == 0L) {

                    print(
                        "processed: " +
                        numEvents +
                        " tag-" +
                        hipoTag +
                        " events, "
                    );

                    print(
                        "current run = " +
                        runnumForEvent +
                        ", "
                    );

                    print(
                        "RUN::scaler min = " +
                        accumulator.runScalerChargeMin +
                        " nC, "
                    );

                    print(
                        "max = " +
                        accumulator.runScalerChargeMax +
                        " nC, "
                    );

                    print(
                        "max-min = " +
                        accumulator.getRunScalerCharge() +
                        " nC.\n"
                    );
                }
            }

            reader.close();

            if (runnumForFile != 0) {

                ChargeAccumulator fileAccumulator =
                    runChargeMap.get(runnumForFile);

                if (fileAccumulator != null) {
                    fileAccumulator.nFiles++;
                }
            }

            println();
            println();
            println();
            print(formatBeamChargeList(runChargeMap));
            println();

            currentFile++;
        }

        println();
        println(
            "Final dst-mode run-combined charge list:"
        );
        println();

        print(formatBeamChargeList(runChargeMap));
        println();
    }


    static String formatBeamChargeList(
        LinkedHashMap<Integer, ChargeAccumulator> runChargeMap
    ) {

        String beamChargeList = "";

        def sortedRunNumbers =
            runChargeMap.keySet().toList().sort();

        for (runnum in sortedRunNumbers) {

            ChargeAccumulator accumulator =
                runChargeMap.get(runnum);

            double runScalerCharge =
                accumulator.getRunScalerCharge();

            beamChargeList +=
                accumulator.runnum.toString() +
                ",";

            beamChargeList +=
                Double.toString(runScalerCharge) +
                ",";

            beamChargeList +=
                Double.toString(
                    accumulator.posHelbeamChargeTotal
                ) +
                ",";

            beamChargeList +=
                Double.toString(
                    accumulator.negHelbeamChargeTotal
                ) +
                ",";

            beamChargeList +=
                Double.toString(
                    accumulator.noHelbeamChargeTotal
                ) +
                ",0,0\n";
        }

        return beamChargeList;
    }


    /*
     * Performs one unrestricted scan over every selected input file.
     *
     * No reader.setTags(...) call is made here.
     *
     * Every RUN::scaler row is entered into:
     *
     *   1. The run's all-tags accumulator.
     *   2. The accumulator for the event's specific tag.
     *   3. The run's tag-1 accumulator when event.getEventTag() == 1.
     */
    static LinkedHashMap<Integer, DiagnosticRunAccumulator>
        runScalerDiagnosticScan(
            def hipoList,
            int nFiles
        ) {

        LinkedHashMap<Integer, DiagnosticRunAccumulator> diagnosticMap =
            new LinkedHashMap<Integer, DiagnosticRunAccumulator>();

        long unrestrictedEventIndex = 0L;
        long scalerRecordIndex = 0L;

        int currentFile = 0;

        while (currentFile < nFiles) {

            File inputFile = hipoList[currentFile];

            println(
                "Diagnostic scan: opening file " +
                Integer.toString(currentFile + 1) +
                " out of " +
                nFiles
            );

            println(
                "File: " +
                inputFile.getAbsolutePath()
            );

            int runnumForFile = 0;

            HipoReader reader = new HipoReader();

            /*
             * Deliberately do not call reader.setTags(...).
             *
             * This scan visits every event tag available through the
             * unrestricted reader.
             */
            reader.open(inputFile.getAbsolutePath());

            SchemaFactory schema = reader.getSchemaFactory();

            Bank runConfigBank = new Bank(
                schema.getSchema("RUN::config")
            );

            Bank runScalerBank = new Bank(
                schema.getSchema("RUN::scaler")
            );

            Event event = new Event();

            while (reader.hasNext()) {

                reader.nextEvent(event);
                unrestrictedEventIndex++;

                int eventTag = event.getEventTag();

                event.read(runConfigBank);
                event.read(runScalerBank);

                int runnumForEvent = runnumForFile;
                int runConfigEvent = -1;
                long timestamp = -1L;

                if (runConfigBank.getRows() > 0) {

                    runnumForEvent =
                        runConfigBank.getInt("run", 0);

                    if (
                        runnumForFile == 0 &&
                        runnumForEvent != 0
                    ) {

                        runnumForFile = runnumForEvent;
                    }

                    runConfigEvent =
                        runConfigBank.getInt("event", 0);

                    timestamp =
                        runConfigBank.getLong("timestamp", 0);
                }

                if (runScalerBank.getRows() <= 0) {
                    continue;
                }

                if (runnumForEvent == 0) {

                    println(
                        "WARNING: Found RUN::scaler with no usable " +
                        "run number."
                    );

                    println(
                        "  file = " +
                        inputFile.getAbsolutePath()
                    );

                    println(
                        "  event tag = " +
                        eventTag
                    );

                    println(
                        "  unrestricted event index = " +
                        unrestrictedEventIndex
                    );

                    continue;
                }

                DiagnosticRunAccumulator runAccumulator =
                    diagnosticMap.get(runnumForEvent);

                if (runAccumulator == null) {

                    runAccumulator =
                        new DiagnosticRunAccumulator();

                    runAccumulator.runnum =
                        runnumForEvent;

                    diagnosticMap.put(
                        runnumForEvent,
                        runAccumulator
                    );
                }

                /*
                 * Count one scaler-bearing event for the unrestricted
                 * stream and for this specific event tag.
                 */
                runAccumulator.allTags.scalerEventCount++;

                ScalerStreamAccumulator tagAccumulator =
                    runAccumulator.getTagAccumulator(eventTag);

                tagAccumulator.scalerEventCount++;

                if (eventTag == 1) {
                    runAccumulator.tag1.scalerEventCount++;
                }

                for (
                    int row = 0;
                    row < runScalerBank.getRows();
                    row++
                ) {

                    scalerRecordIndex++;

                    double scalerCharge =
                        (double) runScalerBank.getFloat(
                            "fcupgated",
                            row
                        );

                    /*
                     * Update unrestricted extrema.
                     */
                    runAccumulator.allTags.update(
                        scalerCharge,
                        runnumForEvent,
                        eventTag,
                        runConfigEvent,
                        timestamp,
                        inputFile.getAbsolutePath(),
                        currentFile + 1,
                        row,
                        unrestrictedEventIndex,
                        scalerRecordIndex
                    );

                    /*
                     * Update extrema for this exact event tag.
                     */
                    tagAccumulator.update(
                        scalerCharge,
                        runnumForEvent,
                        eventTag,
                        runConfigEvent,
                        timestamp,
                        inputFile.getAbsolutePath(),
                        currentFile + 1,
                        row,
                        unrestrictedEventIndex,
                        scalerRecordIndex
                    );

                    /*
                     * Update the explicit tag-1 comparison stream.
                     */
                    if (eventTag == 1) {

                        runAccumulator.tag1.update(
                            scalerCharge,
                            runnumForEvent,
                            eventTag,
                            runConfigEvent,
                            timestamp,
                            inputFile.getAbsolutePath(),
                            currentFile + 1,
                            row,
                            unrestrictedEventIndex,
                            scalerRecordIndex
                        );
                    }
                }

                if (
                    unrestrictedEventIndex %
                    1000000L ==
                    0L
                ) {

                    println(
                        "Diagnostic scan processed " +
                        unrestrictedEventIndex +
                        " unrestricted events; " +
                        "latest RUN::scaler event tag = " +
                        eventTag +
                        ", run = " +
                        runnumForEvent
                    );
                }
            }

            reader.close();

            currentFile++;
        }

        println();
        println(
            "Unrestricted diagnostic scan complete."
        );

        println(
            "Total unrestricted events read = " +
            unrestrictedEventIndex
        );

        println(
            "Total RUN::scaler rows found = " +
            scalerRecordIndex
        );

        println();

        return diagnosticMap;
    }


    static void printDiagnosticReport(
        LinkedHashMap<Integer, DiagnosticRunAccumulator> diagnosticMap
    ) {

        println();
        println(
            "============================================================"
        );
        println(
            "RUN::scaler tag-1 versus all-tag diagnostic report"
        );
        println(
            "============================================================"
        );
        println();

        if (diagnosticMap.isEmpty()) {

            println(
                "No RUN::scaler banks were found during the " +
                "unrestricted diagnostic scan."
            );

            return;
        }

        def sortedRunNumbers =
            diagnosticMap.keySet().toList().sort();

        for (runnum in sortedRunNumbers) {

            DiagnosticRunAccumulator runAccumulator =
                diagnosticMap.get(runnum);

            ScalerStreamAccumulator tag1 =
                runAccumulator.tag1;

            ScalerStreamAccumulator allTags =
                runAccumulator.allTags;

            println(
                "------------------------------------------------------------"
            );

            println(
                "Run " +
                runnum
            );

            println(
                "------------------------------------------------------------"
            );

            println();
            println(
                "Per-tag RUN::scaler coverage:"
            );
            println();

            printf(
                Locale.US,
                "%-10s %15s %15s %20s %20s %20s%n",
                "event tag",
                "scaler events",
                "scaler rows",
                "minimum (nC)",
                "maximum (nC)",
                "max-min (nC)"
            );

            printf(
                Locale.US,
                "%-10s %15s %15s %20s %20s %20s%n",
                "---------",
                "-------------",
                "-----------",
                "------------",
                "------------",
                "------------"
            );

            def sortedTags =
                runAccumulator.byTag.keySet().toList().sort();

            for (tag in sortedTags) {

                ScalerStreamAccumulator tagAccumulator =
                    runAccumulator.byTag.get(tag);

                printf(
                    Locale.US,
                    "%-10d %15d %15d %20.9f %20.9f %20.9f%n",
                    tag,
                    tagAccumulator.scalerEventCount,
                    tagAccumulator.scalerRowCount,
                    tagAccumulator.minimum,
                    tagAccumulator.maximum,
                    tagAccumulator.getCharge()
                );
            }

            println();

            printf(
                Locale.US,
                "%-10s %15d %15d %20.9f %20.9f %20.9f%n",
                "ALL",
                allTags.scalerEventCount,
                allTags.scalerRowCount,
                allTags.minimum,
                allTags.maximum,
                allTags.getCharge()
            );

            println();
            println(
                "Tag-1 versus all-tag charge comparison:"
            );
            println();

            if (!tag1.foundScaler) {

                println(
                    "No tag-1 RUN::scaler records were found for " +
                    "this run."
                );

                println(
                    "All-tag charge = " +
                    formatCharge(allTags.getCharge()) +
                    " nC"
                );

                println();

                printEndpointSection(
                    "All-tag minimum",
                    allTags.minimumEndpoint
                );

                printEndpointSection(
                    "All-tag maximum",
                    allTags.maximumEndpoint
                );

                println();

                continue;
            }

            double tag1Charge =
                tag1.getCharge();

            double allTagCharge =
                allTags.getCharge();

            double chargeDifference =
                allTagCharge -
                tag1Charge;

            double percentRelativeToTag1 =
                tag1Charge != 0.0d ?
                    100.0d *
                    chargeDifference /
                    tag1Charge :
                    Double.NaN;

            double percentRelativeToAll =
                allTagCharge != 0.0d ?
                    100.0d *
                    chargeDifference /
                    allTagCharge :
                    Double.NaN;

            double lowerEndpointExtension =
                tag1.minimum -
                allTags.minimum;

            double upperEndpointExtension =
                allTags.maximum -
                tag1.maximum;

            double endpointExtensionSum =
                lowerEndpointExtension +
                upperEndpointExtension;

            printf(
                Locale.US,
                "Tag-1 minimum:                     %.9f nC%n",
                tag1.minimum
            );

            printf(
                Locale.US,
                "Tag-1 maximum:                     %.9f nC%n",
                tag1.maximum
            );

            printf(
                Locale.US,
                "Tag-1 max-min charge:              %.9f nC%n",
                tag1Charge
            );

            println();

            printf(
                Locale.US,
                "All-tag minimum:                   %.9f nC%n",
                allTags.minimum
            );

            printf(
                Locale.US,
                "All-tag maximum:                   %.9f nC%n",
                allTags.maximum
            );

            printf(
                Locale.US,
                "All-tag max-min charge:            %.9f nC%n",
                allTagCharge
            );

            println();

            printf(
                Locale.US,
                "All-tag minus tag-1 charge:        %.9f nC%n",
                chargeDifference
            );

            printf(
                Locale.US,
                "Difference relative to tag 1:      %.9f%%%n",
                percentRelativeToTag1
            );

            printf(
                Locale.US,
                "Difference relative to all tags:   %.9f%%%n",
                percentRelativeToAll
            );

            println();

            printf(
                Locale.US,
                "Lower-end extension:%n" +
                "  tag1_min - all_min =             %.9f nC%n",
                lowerEndpointExtension
            );

            printf(
                Locale.US,
                "Upper-end extension:%n" +
                "  all_max - tag1_max =             %.9f nC%n",
                upperEndpointExtension
            );

            printf(
                Locale.US,
                "Sum of endpoint extensions:        %.9f nC%n",
                endpointExtensionSum
            );

            printf(
                Locale.US,
                "Charge difference minus extension sum: %.12f nC%n",
                chargeDifference -
                endpointExtensionSum
            );

            println();

            if (
                approximatelyEqual(
                    chargeDifference,
                    endpointExtensionSum
                )
            ) {

                println(
                    "CHECK: The all-tag/tag-1 charge difference is " +
                    "fully explained by the different observed " +
                    "minimum and maximum endpoints."
                );

            } else {

                println(
                    "CHECK: The endpoint-extension sum does not " +
                    "numerically reproduce the reported charge " +
                    "difference within tolerance."
                );
            }

            println();

            printEndpointSection(
                "Tag-1 minimum endpoint",
                tag1.minimumEndpoint
            );

            printEndpointSection(
                "All-tag minimum endpoint",
                allTags.minimumEndpoint
            );

            printEndpointSection(
                "Tag-1 maximum endpoint",
                tag1.maximumEndpoint
            );

            printEndpointSection(
                "All-tag maximum endpoint",
                allTags.maximumEndpoint
            );

            println();

            printEndpointInterpretation(
                tag1,
                allTags
            );

            println();
        }

        println(
            "============================================================"
        );

        println(
            "End of RUN::scaler diagnostic report"
        );

        println(
            "============================================================"
        );

        println();
    }


    static void printEndpointSection(
        String heading,
        ScalerEndpoint endpoint
    ) {

        println(
            heading +
            ":"
        );

        println(
            "        " +
            endpoint.formatLocation()
        );

        println();
    }


    static void printEndpointInterpretation(
        ScalerStreamAccumulator tag1,
        ScalerStreamAccumulator allTags
    ) {

        boolean minimumIsSameRecord =
            endpointsMatch(
                tag1.minimumEndpoint,
                allTags.minimumEndpoint
            );

        boolean maximumIsSameRecord =
            endpointsMatch(
                tag1.maximumEndpoint,
                allTags.maximumEndpoint
            );

        println(
            "Endpoint interpretation:"
        );

        println(
            "  Same minimum record: " +
            minimumIsSameRecord
        );

        println(
            "  Same maximum record: " +
            maximumIsSameRecord
        );

        if (
            minimumIsSameRecord &&
            maximumIsSameRecord
        ) {

            println(
                "  The tag-1 and unrestricted streams use the " +
                "same endpoint records."
            );

            return;
        }

        if (!minimumIsSameRecord) {

            println(
                "  The unrestricted stream obtains its minimum " +
                "from a different event record."
            );

            println(
                "  That additional minimum is on event tag " +
                allTags.minimumEndpoint.tag +
                "."
            );
        }

        if (!maximumIsSameRecord) {

            println(
                "  The unrestricted stream obtains its maximum " +
                "from a different event record."
            );

            println(
                "  That additional maximum is on event tag " +
                allTags.maximumEndpoint.tag +
                "."
            );
        }
    }


    static boolean endpointsMatch(
        ScalerEndpoint first,
        ScalerEndpoint second
    ) {

        if (
            first == null ||
            second == null ||
            !first.valid ||
            !second.valid
        ) {

            return false;
        }

        return (
            first.tag ==
            second.tag &&
            first.runnum ==
            second.runnum &&
            first.runConfigEvent ==
            second.runConfigEvent &&
            first.timestamp ==
            second.timestamp &&
            first.filename.equals(second.filename) &&
            first.scalerRow ==
            second.scalerRow &&
            Double.compare(
                first.charge,
                second.charge
            ) ==
            0
        );
    }


    static boolean approximatelyEqual(
        double first,
        double second
    ) {

        double scale = Math.max(
            1.0d,
            Math.max(
                Math.abs(first),
                Math.abs(second)
            )
        );

        return (
            Math.abs(first - second) <=
            1.0e-10d * scale
        );
    }


    static String formatCharge(double charge) {

        return String.format(
            Locale.US,
            "%.9f",
            charge
        );
    }
}