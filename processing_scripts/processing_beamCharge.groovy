/*
 * author Timothy B. Hayward
 * 
 */

import java.io.File;
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; 
import analyzers.*;

// filetype for gathering files in directory
import groovy.io.FileType;


public class processing_beamCharge {

	public static void main(String[] args) {
		// ~~~~~~~~~~~~~~~~ set up input paramaeters ~~~~~~~~~~~~~~~~ //

		// Check if an argument is provided
		if (!args) {
		    // Print an error message and exit the program if the input directory is not specified
		    println("ERROR: Please enter a hipo file directory as the first argument");
		    System.exit(0);
		}
		// If the input directory is provided, iterate through each file recursively
		def hipo_list = []
		(args[0] as File).eachFileRecurse(FileType.FILES) 
			{ if (it.name.endsWith('.hipo')) hipo_list << it }

		// Set the number of files to process based on the provided 4th argument
		// use the size of the hipo_list if no argument provided
		int n_files = args.length < 2 || Integer.parseInt(args[1]) > hipo_list.size()
		    ? hipo_list.size() : Integer.parseInt(args[1]);
		if (args.length < 2 || Integer.parseInt(args[1]) > hipo_list.size()) {
		    // Print warnings and information if the number of files is not specified or too large
		    println("WARNING: Number of files not specified or number too large.")
		    println("Setting # of files to be equal to number of files in the directory.");
		    println("There are $hipo_list.size files.");
		}

		int num_events = 0;
		int current_file = 0;
		String beamChargeList = '';
		float beamChargeMax = 0;
		float posHelbeamChargeTotal = 0;
		float negHelbeamChargeTotal = 0;
		int runnum;
		while (current_file < n_files) {
			beamChargeMax = 0;
			runnum = 0;
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files); println(); println();
			// limit to a certain number of files defined by n_files

			HipoDataSource reader = new HipoDataSource();

			reader.open(hipo_list[current_file]); // open next hipo file
			current_file++;
			HipoDataEvent event = reader.getNextEvent(); 

			while(reader.hasEvent()==true){
				num_events++; 
				if (num_events%1000000 == 0) { // not necessary, just updates output
					print("processed: "+num_events+" events, max beamCharge of current ");
					print("run = "+beamChargeMax+" nC.\n");
				}

				if (event.hasBank("RUN::scaler")) {
					float beamCharge = event.getBank("RUN::scaler").getFloat("fcupgated",0);
					if (beamCharge > beamChargeMax) { beamChargeMax = beamCharge; }
    			}

    			// if (event.hasBank("HEL::scaler") && 
    			// 	event.getBank("HEL::scaler").getFloat("clock",0) > 30000) {
    			if (event.hasBank("HEL::scaler")) {
    				if (event.getBank("HEL::scaler").getInt("helicity",0) == 1) {
    					float beamCharge = event.getBank("HEL::scaler").getFloat("fcupgated",0);
    					posHelbeamChargeTotal+=beamCharge;
					} else if (event.getBank("HEL::scaler").getInt("helicity",0) == -1) {
						float beamCharge = event.getBank("HEL::scaler").getFloat("fcupgated",0);
    					negHelbeamChargeTotal+=beamCharge;
					}
    			}

				// get run and event numbers
				event = reader.getNextEvent();
				if (runnum==0) {
					runnum = event.getBank("RUN::config").getInt('run',0);
				}

			}
			beamChargeList+=runnum.toString()+","+beamChargeMax.toString()+",";
			beamChargeList+=posHelbeamChargeTotal.toString()+","
			beamChargeList+=negHelbeamChargeTotal.toString()+",0,0\n";
			println(); println(); println();
			print(beamChargeList);
			println();
			beamChargeMax = 0;
			posHelbeamChargeTotal = 0;
			negHelbeamChargeTotal = 0;
		}

	}
}