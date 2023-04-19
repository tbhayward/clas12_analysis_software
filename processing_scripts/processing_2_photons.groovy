/*
 * author Timothy B. Hayward
 * 
 * SIDIS dihadron 
 */

// import CLAS12 physics classes
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

public static void main(String[] args) {

	// Start time
	long startTime = System.currentTimeMillis();

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

	// Set the output file name based on the provided 3rd argument or use the default name
	String output_file = args.length < 2 ? "hadron_dummy_out.txt" : args[1];
	if (args.length < 2) 
	    println("WARNING: Specify an output file name. Set to \"hadron_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

	// Set the number of files to process based on the provided 4th argument
	// use the size of the hipo_list if no argument provided
	int n_files = args.length < 5 || Integer.parseInt(args[4]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[4]);
	if (args.length < 5 || Integer.parseInt(args[4]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified or too large
	    println("WARNING: Number of files not specified or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

	// ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

	// declare physics event variables
	double gamma_p_1, gamma_p_2;
	double e_p;
	double beta_1, beta_2;
	double lv_1, lv_2, lw_1, lw_2;
	double opening_angle_1, opening_angle_2;
	double gamma_gamma_mass, gamma_gamma_p;

	// load my kinematic fitter/PID
	// GenericKinematicFitter fitter = new analysis_fitter(10.6041); 
	// GenericKinematicFitter fitter = new monte_carlo_fitter(10.6041);
	GenericKinematicFitter fitter = new event_builder_fitter(10.6041); 
	// GenericKinematicFitter fitter = new proton_energy_loss_corrections_fitter(10.6041); 
	
	// set filter for final states
	EventFilter filter = new EventFilter("11:22:22:X+:X-:Xn"); 
	
	// setup QA database
	QADB qa = new QADB();

	// create a StringBuilder for accumulating lines
	StringBuilder batchLines = new StringBuilder();

	int num_events = 0;
	int max_lines = 1000;
	int lineCount = 0;
	for (current_file in 0..<n_files) {
		// limit to a certain number of files defined by n_files
		println("\n Opening file "+Integer.toString(current_file+1)
			+" out of "+n_files+".\n"); 
		
		HipoDataSource reader = new HipoDataSource();
		reader.open(hipo_list[current_file]); // open next hipo file
		HipoDataEvent event = reader.getNextEvent(); 

		while (reader.hasEvent()) {
		    ++num_events;
		    if (num_events % 100000 == 0) { // not necessary, just updates output
		        print("processed: " + num_events + " events. ");
		    }

		    // get run and event numbers
		    event = reader.getNextEvent();
		    // collect info for QA
		    int runnum = event.getBank("RUN::config").getInt('run', 0);
		    int evnum = event.getBank("RUN::config").getInt('event', 0);

		    PhysicsEvent research_Event = fitter.getPhysicsEvent(event);

		    // do not use the qa if it is MC (runnum = 11) 
		    // do not use the qa if the run is from RGC (until QA is produced!)
		    boolean process_event = filter.isValid(research_Event) && 
		    	(runnum == 11 || runnum >= 11571 || qa.OkForAsymmetry(runnum, evnum));

		    if (process_event) {
		    	// load particle bank
				HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle"); 
	            HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
	            HipoDataBank cc_Bank = (HipoDataBank) event.getBank("REC::Cherenkov");
	            HipoDataBank track_Bank = (HipoDataBank) event.getBank("REC::Track");
	            HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");
	            HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config");

		        // get # of particles 
		        int num_photons = research_Event.countByPid(22);

		        // cycle over all photons
		        for (int current_p1 = 0; current_p1 < num_photons; current_p1++) { 
		        	for (int current_p2 = 0; current_p2 < num_photons; current_p2++) { 
		        		if (current_p1 == current_p2) {continue; }
		        		Particle e = rec_Bank.getParticle([11,0]);
		        		LorentzVector lv_e = new LorentzVector();
				        lv_e.setPxPyPzM(e.px(), e.py(), e.pz(), e.mass());
		        		Particle p1 = rec_Bank.getParticle([22,current_p1]);
		        		LorentzVector lv_p1 = new LorentzVector();
				        lv_p1.setPxPyPzM(p1.px(), p1.py(), p1.pz(), p1.mass());
		        		Particle p2 = rec_Bank.getParticle([22,current_p2]);
				        LorentzVector lv_p2 = new LorentzVector();
				        lv_p2.setPxPyPzM(p2.px(), p2.py(), p2.pz(), p2.mass());
				        Particle p12 = rec_Bank.getParticle([22,current_p1]+[22,current_p2]);
				        LorentzVector lv_p12 = new LorentzVector();
				        lv_p12.setPxPyPzM(p12.px(), p12.py(), p12.pz(), p12.mass());

				        double open_angle_1 = 180/Math.PI*Math.acos(lv_e.vect().dot(lv_p1.vect())/
            				(lv_e.vect().mag()*lv_p1.vect().mag()));
						double open_angle_2 = 180/Math.PI*Math.acos(lv_e.vect().dot(lv_p2.vect())/
            				(lv_e.vect().mag()*lv_p2.vect().mag()));

		                // Use a StringBuilder to append all data in a single call
		                StringBuilder line = new StringBuilder();
		                line.append(lv_p1.p()).append(" ")
		                	.append(lv_p2.p()).append(" ")
		                	.append(0).append(" ")
		                	.append(0).append(" ")
		                	.append(0).append(" ")
		                	.append(0).append(" ")
		                	.append(0).append(" ")
		                	.append(0).append(" ")
		                	.append(opening_angle_1).append(" ")
		                	.append(opening_angle_2).append(" ")
		                	.append(lv_p12.mass()).append(" ")
		                	.append(lv_p12.p()).append("\n");

		                // Append the line to the batchLines StringBuilder
		                batchLines.append(line.toString());
		                lineCount++; // Increment the line count

		                // If the line count reaches 1000, write to the file and reset
		                if (lineCount >= max_lines) {
		                    file.append(batchLines.toString());
		                    batchLines.setLength(0);
		                    lineCount = 0;
		                }
			        }
		        }
		    }
		reader.close();
		}

		// Write any remaining lines in the batchLines StringBuilder to the file
		if (batchLines.length() > 0) {
		    file.append(batchLines.toString());
		    batchLines.setLength(0);
		}
	}

	writer.close();

	// End time
	long endTime = System.currentTimeMillis()
	// Calculate the elapsed time
	long elapsedTime = endTime - startTime
	// Print the elapsed time in milliseconds
	println("Elapsed time: ${elapsedTime} ms");

}