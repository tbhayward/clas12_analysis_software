/*
 * author Timothy B. Hayward
 *
 * extract RUN::scaler / HEL::scaler charge.
 */

import java.io.File;
import groovy.io.FileType;

import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.SchemaFactory;

public class processing_beamCharge {

    static class ChargeAccumulator {
        int runnum = 0;

        float runScalerChargeMin = Float.POSITIVE_INFINITY;
        float runScalerChargeMax = Float.NEGATIVE_INFINITY;
        boolean foundRunScalerCharge = false;

        float posHelbeamChargeTotal = 0.0f;
        float negHelbeamChargeTotal = 0.0f;
        float noHelbeamChargeTotal = 0.0f;

        int nFiles = 0;

        float getRunScalerCharge() {
            if (!foundRunScalerCharge) {
                return 0.0f;
            }

            return runScalerChargeMax - runScalerChargeMin;
        }
    }

    public static void main(String[] args) {

        // ~~~~~~~~~~~~~~~~ set up input parameters ~~~~~~~~~~~~~~~~ //

        if (!args || args.length < 1) {
            println("ERROR: Please enter a hipo file directory as the first argument");
            System.exit(0);
        }

        def hipo_list = [];

        (args[0] as File).eachFileRecurse(FileType.FILES) {
            if (it.name.endsWith(".hipo")) {
                hipo_list << it;
            }
        }

        hipo_list = hipo_list.sort { it.getAbsolutePath() };

        int n_files = hipo_list.size();

        if (args.length >= 2) {
            int requested_files = Integer.parseInt(args[1]);

            if (requested_files > 0 && requested_files <= hipo_list.size()) {
                n_files = requested_files;
            } else {
                println("WARNING: Number of files not specified, non-positive, or number too large.");
                println("Setting # of files to be equal to number of files in the directory.");
                println("There are ${hipo_list.size()} files.");
            }
        } else {
            println("WARNING: Number of files not specified or number too large.");
            println("Setting # of files to be equal to number of files in the directory.");
            println("There are ${hipo_list.size()} files.");
        }

        String processingMode = "train";

        if (args.length >= 3) {
            processingMode = args[2].toLowerCase();
        }

        if (!(processingMode.equals("train") || processingMode.equals("dst"))) {
            println("ERROR: Unknown processing mode: " + processingMode);
            println("Allowed modes are:");
            println("  train");
            println("  dst");
            System.exit(0);
        }

        println();
        println("Processing mode = " + processingMode);
        println();

        int hipo_tag = 1;

        if (processingMode.equals("dst")) {
            processDstMode(hipo_list, n_files, hipo_tag);
        } else {
            processTrainMode(hipo_list, n_files, hipo_tag);
        }
    }

    static void processTrainMode(def hipo_list, int n_files, int hipo_tag) {

        int num_events = 0;
        int current_file = 0;
        String beamChargeList = "";

        while (current_file < n_files) {

            File inputFile = hipo_list[current_file];

            println();
            println("Opening file " + Integer.toString(current_file + 1) + " out of " + n_files);
            println("File: " + inputFile.getAbsolutePath());
            println();

            float beamChargeMax = 0.0f;

            float runScalerChargeMin = Float.POSITIVE_INFINITY;
            float runScalerChargeMax = Float.NEGATIVE_INFINITY;
            float runScalerCharge = 0.0f;
            boolean foundRunScalerCharge = false;

            float posHelbeamChargeTotal = 0.0f;
            float negHelbeamChargeTotal = 0.0f;
            float noHelbeamChargeTotal = 0.0f;

            int runnum = 0;

            HipoReader reader = new HipoReader();
            reader.setTags(hipo_tag);
            reader.open(inputFile.getAbsolutePath());

            SchemaFactory schema = reader.getSchemaFactory();
            Bank runConfigBank = new Bank(schema.getSchema("RUN::config"));
            Bank runScalerBank = new Bank(schema.getSchema("RUN::scaler"));
            Bank helScalerBank = new Bank(schema.getSchema("HEL::scaler"));

            Event event = new Event();

            while (reader.hasNext()) {

                reader.nextEvent(event);
                num_events++;

                event.read(runConfigBank);
                event.read(runScalerBank);
                event.read(helScalerBank);

                if (runConfigBank.getRows() > 0 && runnum == 0) {
                    runnum = runConfigBank.getInt("run", 0);
                }

                if (runScalerBank.getRows() > 0) {

                    /*
                     * RUN::scaler fcupgated is cumulative in the scaler stream.
                     */
                    float thisRunScalerCharge = runScalerBank.getFloat("fcupgated", 0);

                    if (thisRunScalerCharge < runScalerChargeMin) {
                        runScalerChargeMin = thisRunScalerCharge;
                    }

                    if (thisRunScalerCharge > runScalerChargeMax) {
                        runScalerChargeMax = thisRunScalerCharge;
                    }

                    if (thisRunScalerCharge > beamChargeMax) {
                        beamChargeMax = thisRunScalerCharge;
                    }

                    foundRunScalerCharge = true;
                    runScalerCharge = runScalerChargeMax - runScalerChargeMin;
                }

                if (helScalerBank.getRows() > 0) {

                    for (int row = 0; row < helScalerBank.getRows(); row++) {

                        int helicity = helScalerBank.getInt("helicity", row);
                        float helBeamCharge = helScalerBank.getFloat("fcupgated", row);

                        if (helicity == 1) {
                            posHelbeamChargeTotal += helBeamCharge;
                        } else if (helicity == -1) {
                            negHelbeamChargeTotal += helBeamCharge;
                        } else if (helicity == 0) {
                            noHelbeamChargeTotal += helBeamCharge;
                        }
                    }
                }

                if (num_events % 1000000 == 0) {
                    print("processed: " + num_events + " tag-" + hipo_tag + " events, ");
                    print("max beamCharge of current run/file = " + beamChargeMax + " nC, ");
                    print("RUN::Scaler max-min of current run/file = " + runScalerCharge + " nC.\n");
                }
            }

            reader.close();

            if (!foundRunScalerCharge) {
                runScalerCharge = 0.0f;
            }

            beamChargeList += runnum.toString() + "," + runScalerCharge.toString() + ",";
            beamChargeList += posHelbeamChargeTotal.toString() + ",";
            beamChargeList += negHelbeamChargeTotal.toString() + ",";
            beamChargeList += noHelbeamChargeTotal.toString() + ",0,0\n";

            println();
            println();
            println();
            print(beamChargeList);
            println();

            current_file++;
        }
    }

    static void processDstMode(def hipo_list, int n_files, int hipo_tag) {

        int num_events = 0;
        int current_file = 0;

        LinkedHashMap<Integer, ChargeAccumulator> runChargeMap = new LinkedHashMap<Integer, ChargeAccumulator>();

        while (current_file < n_files) {

            File inputFile = hipo_list[current_file];

            println();
            println("Opening file " + Integer.toString(current_file + 1) + " out of " + n_files);
            println("File: " + inputFile.getAbsolutePath());
            println();

            int runnumForFile = 0;

            HipoReader reader = new HipoReader();
            reader.setTags(hipo_tag);
            reader.open(inputFile.getAbsolutePath());

            SchemaFactory schema = reader.getSchemaFactory();
            Bank runConfigBank = new Bank(schema.getSchema("RUN::config"));
            Bank runScalerBank = new Bank(schema.getSchema("RUN::scaler"));
            Bank helScalerBank = new Bank(schema.getSchema("HEL::scaler"));

            Event event = new Event();

            while (reader.hasNext()) {

                reader.nextEvent(event);
                num_events++;

                event.read(runConfigBank);
                event.read(runScalerBank);
                event.read(helScalerBank);

                if (runConfigBank.getRows() > 0 && runnumForFile == 0) {
                    runnumForFile = runConfigBank.getInt("run", 0);
                }

                int runnumForEvent = runnumForFile;

                if (runConfigBank.getRows() > 0) {
                    runnumForEvent = runConfigBank.getInt("run", 0);
                }

                if (runnumForEvent == 0) {
                    continue;
                }

                ChargeAccumulator acc = runChargeMap.get(runnumForEvent);

                if (acc == null) {
                    acc = new ChargeAccumulator();
                    acc.runnum = runnumForEvent;
                    runChargeMap.put(runnumForEvent, acc);
                }

                if (runScalerBank.getRows() > 0) {

                    /*
                     * In dst mode, there can be many files per run.
                     *
                     * The RUN::scaler fcupgated value is cumulative across the run, so
                     * we need the smallest and largest values seen across every file
                     * belonging to the same run.
                     */
                    float thisRunScalerCharge = runScalerBank.getFloat("fcupgated", 0);

                    if (thisRunScalerCharge < acc.runScalerChargeMin) {
                        acc.runScalerChargeMin = thisRunScalerCharge;
                    }

                    if (thisRunScalerCharge > acc.runScalerChargeMax) {
                        acc.runScalerChargeMax = thisRunScalerCharge;
                    }

                    acc.foundRunScalerCharge = true;
                }

                if (helScalerBank.getRows() > 0) {

                    for (int row = 0; row < helScalerBank.getRows(); row++) {

                        int helicity = helScalerBank.getInt("helicity", row);
                        float helBeamCharge = helScalerBank.getFloat("fcupgated", row);

                        if (helicity == 1) {
                            acc.posHelbeamChargeTotal += helBeamCharge;
                        } else if (helicity == -1) {
                            acc.negHelbeamChargeTotal += helBeamCharge;
                        } else if (helicity == 0) {
                            acc.noHelbeamChargeTotal += helBeamCharge;
                        }
                    }
                }

                if (num_events % 1000000 == 0) {
                    print("processed: " + num_events + " tag-" + hipo_tag + " events, ");
                    print("current run = " + runnumForEvent + ", ");
                    print("RUN::Scaler min = " + acc.runScalerChargeMin + " nC, ");
                    print("max = " + acc.runScalerChargeMax + " nC, ");
                    print("max-min = " + acc.getRunScalerCharge() + " nC.\n");
                }
            }

            reader.close();

            if (runnumForFile != 0) {
                ChargeAccumulator fileAcc = runChargeMap.get(runnumForFile);

                if (fileAcc != null) {
                    fileAcc.nFiles++;
                }
            }

            println();
            println();
            println();
            print(formatBeamChargeList(runChargeMap));
            println();

            current_file++;
        }

        println();
        println("Final dst-mode run-combined charge list:");
        println();
        print(formatBeamChargeList(runChargeMap));
        println();
    }

    static String formatBeamChargeList(LinkedHashMap<Integer, ChargeAccumulator> runChargeMap) {

        String beamChargeList = "";

        def sortedRunNumbers = runChargeMap.keySet().toList().sort();

        for (runnum in sortedRunNumbers) {

            ChargeAccumulator acc = runChargeMap.get(runnum);

            float runScalerCharge = acc.getRunScalerCharge();

            beamChargeList += acc.runnum.toString() + ",";
            beamChargeList += runScalerCharge.toString() + ",";
            beamChargeList += acc.posHelbeamChargeTotal.toString() + ",";
            beamChargeList += acc.negHelbeamChargeTotal.toString() + ",";
            beamChargeList += acc.noHelbeamChargeTotal.toString() + ",0,0\n";
        }

        return beamChargeList;
    }
}