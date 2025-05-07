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

// dilks CLAS QA analysis
import clasqa.QADB


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
		float noHelbeamChargeTotal = 0;
		def HLstate = [-1,0,1]

		// setup QA database
		QADB qa = new QADB("latest");

		while (current_file < n_files) {
			beamChargeMax = 0;
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files); println(); println();
			// limit to a certain number of files defined by n_files

			HipoDataSource reader = new HipoDataSource();

			reader.open(hipo_list[current_file]); // open next hipo file
			current_file++;
			HipoDataEvent event = reader.getNextEvent(); 

			int runnum = event.getBank("RUN::config").getInt('run', 0);
			int evnum = event.getBank("RUN::config").getInt('event', 0);

			while(reader.hasEvent()==true){
				qa.query(runnum,evnum);
				qa.accumulateChargeHL();
			}
		}
		HLstate.each{ value ->
		    println "HL charge(" + value + ")= " + qa.getAccumulatedChargeHL(value)
		}
	}
}