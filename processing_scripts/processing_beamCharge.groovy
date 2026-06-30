/*
 * author Timothy B. Hayward
 *
 * Read only HIPO tag-1 events and extract RUN::scaler / HEL::scaler charge.
 */

import java.io.File;
import groovy.io.FileType;

import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.SchemaFactory;

public class processing_beamCharge {

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

        int n_files = args.length < 2 || Integer.parseInt(args[1]) > hipo_list.size()
            ? hipo_list.size()
            : Integer.parseInt(args[1]);

        if (args.length < 2 || Integer.parseInt(args[1]) > hipo_list.size()) {
            println("WARNING: Number of files not specified or number too large.");
            println("Setting # of files to be equal to number of files in the directory.");
            println("There are ${hipo_list.size()} files.");
        }

        // HIPO event tag to process.
        // Tag 1 is usually the scaler / configuration stream.
        int hipo_tag = 1;

        int num_events = 0;
        int current_file = 0;

        String beamChargeList = "";

        while (current_file < n_files) {

            File inputFile = hipo_list[current_file];

            println();
            println();
            println("Opening file " + Integer.toString(current_file + 1) + " out of " + n_files);
            println("File: " + inputFile.getAbsolutePath());
            println("Reading only HIPO tag-" + hipo_tag + " events");
            println();
            println();

            float beamChargeMax = 0.0f;
            float runScalerCharge = 0.0f;
            float posHelbeamChargeTotal = 0.0f;
            float negHelbeamChargeTotal = 0.0f;
            float noHelbeamChargeTotal = 0.0f;

            int runnum = 0;

            HipoReader reader = new HipoReader();

            /*
             * This is the important line.
             *
             * It must be called before reader.open(...).
             * Only events with HIPO tag 1 will be returned by reader.nextEvent(event).
             */
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

                    if (thisRunScalerCharge > runScalerCharge) {
                        runScalerCharge = thisRunScalerCharge;
                    }

                    if (thisRunScalerCharge > beamChargeMax) {
                        beamChargeMax = thisRunScalerCharge;
                    }
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
                    print("max beamCharge of current run/file = " + beamChargeMax + " nC.\n");
                }
            }

            reader.close();

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
}