/*
 * author Timothy B. Hayward
 * 
 * SIDIS hadron 
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

def executeShellCommand(String command) {
    ProcessBuilder processBuilder = new ProcessBuilder(command.split(" "))
    processBuilder.redirectErrorStream(true)
    Process process = processBuilder.start()
    process.inputStream.text
}

public static void main(String[] args) {

	// Start time
	long startTime = System.currentTimeMillis();

	String shellName = executeShellCommand("echo \$0").trim()
	if (shellName.contains("csh") || shellName.contains("tcsh")) {
	    executeShellCommand("source clasqaDB/env.csh")
	} else {
	    executeShellCommand("source clasqaDB/env.sh")
	}
	// dilks CLAS QA analysis
	import clasqa.QADB


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

	// Set the PDG PID for p1 based on the provided 2nd argument or default to 211 (pi+)
	String p1_Str = args.length < 2 ? "211" : args[1];
	if (args.length < 2) println("WARNING: Specify a PDG PID for p1! Set to pi+ (211).");
	println("Set p1 PID = $p1_Str");
	int p1_int = p1_Str.toInteger(); // Convert p1_Str to integer

	// Set the output file name based on the provided 3rd argument or use the default name
	String output_file = args.length < 3 ? "hadron_dummy_out.txt" : args[2];
	if (args.length < 3) 
	    println("WARNING: Specify an output file name. Set to \"hadron_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

	// Set the number of files to process based on the provided 4th argument
	// use the size of the hipo_list if no argument provided
	int n_files = args.length < 4 || Integer.parseInt(args[3]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[3]);
	if (args.length < 4 || Integer.parseInt(args[3]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified or too large
	    println("WARNING: Number of files not specified or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

	// ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

	// declare physics event variables
	int helicity;
	double e_p, e_theta, e_phi, p_phi, p_p, p_theta;
	double Q2, W, y, Mx, Mx2, x, z, xF, pT, eta, zeta, phi, vz_e, vz_p;
	double Depolarization_A, Depolarization_B, Depolarization_C;
	double Depolarization_V, Depolarization_W;

	// load my kinematic fitter/PID
	GenericKinematicFitter fitter = new analysis_fitter(10.6041); 
	// GenericKinematicFitter fitter = new monte_carlo_fitter(10.6041);
	// GenericKinematicFitter fitter = new event_builder_fitter(10.6041); 
	// GenericKinematicFitter fitter = new proton_energy_loss_corrections_fitter(10.6041); 
	
	// set filter for final states
	EventFilter filter = new EventFilter("11:"+p1_Str+":X+:X-:Xn"); 
	
	// setup QA database
	QADB qa = new QADB();

	// create a StringBuilder for accumulating lines
	StringBuilder batchLines = new StringBuilder();

	int num_events = 0;
	int max_lines = 1000;
	int lineCount = 0;
	// for (int current_file; current_file<n_files; current_file++) {
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

		        // get # of particles w/ pid1
		        int num_p1 = research_Event.countByPid(p1_Str.toInteger()); 

		        // cycle over all hadrons
		        for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { 

		            Hadron variables = new Hadron(event, research_Event,
		                    p1_int, current_p1);
		            // this is my class for defining all relevant kinematic variables

		            if (variables.channel_test(variables)) {
		                helicity = variables.get_helicity(); // helicity of event

		                // lab kinematics
		                e_p = variables.e_p(); // lab frame momentum
		                e_theta = variables.e_theta(); // lab polar angle
		                e_phi = variables.e_phi(); // lab azimuthal angle
		                p_phi = variables.p_phi(); // lab azimuthal angle
		                p_p = variables.p_p(); // lab momentum
		                p_theta = variables.p_theta(); // lab polar angle

		                // DIS variables
		                Q2 = variables.Q2(); // exchanged virtual photon energy
		                W = variables.W(); // hadronic mass
		                y = variables.y(); // E_scat/E_beam
		                Mx = variables.Mx(); // missing mass
		                Mx2 = variables.Mx2(); // missing mass square

		                // SIDIS variables
		                x = variables.x(); // Bjorken-x
		                z = variables.z(); // fractional hadron energy wrt virtual photon
		                xF = variables.xF(); // Feynman-x
		                pT = variables.pT(); // transverse momentum of hadron
		                eta = variables.eta(); // rapidity
		                zeta = variables.zeta(); // fractional longitudinal momentum of hadron

		                // angles
		                phi = variables.phi(); // trento phi of the hadron

		                // vertices
		                vz_e = variables.vz_e();
		                vz_p = variables.vz_p();

		                // depolarization factors
		                Depolarization_A = variables.Depolarization_A();
		                Depolarization_B = variables.Depolarization_B();
		                Depolarization_C = variables.Depolarization_C();
		                Depolarization_V = variables.Depolarization_V();
				    	Depolarization_W = variables.Depolarization_W();

		                // Use a StringBuilder to append all data in a single call
		                StringBuilder line = new StringBuilder();
		                line.append(runnum).append(" ")
		                	.append(evnum).append(" ")
		                	.append(helicity).append(" ")
		                	.append(e_p).append(" ")
		                	.append(e_theta).append(" ")
		                	.append(e_phi).append(" ")
		                	.append(vz_e).append(" ")
		                	.append(p_p).append(" ")
		                	.append(p_theta).append(" ")
		                	.append(p_phi).append(" ")
		                	.append(vz_p).append(" ")
		                	.append(Q2).append(" ")
		                	.append(W).append(" ")
		                	.append(Mx).append(" ")
		                	.append(Mx2).append(" ")
		                	.append(x).append(" ")
		                	.append(y).append(" ")
		                	.append(z).append(" ")
		                	.append(xF).append(" ")
		                	.append(pT).append(" ")
		                	.append(zeta).append(" ")
		                	.append(eta).append(" ")
		                	.append(phi).append(" ")
		                	.append(Depolarization_A).append(" ")
		                    .append(Depolarization_B).append(" ")
		                    .append(Depolarization_C).append(" ")
		                    .append(Depolarization_V).append(" ")
		                    .append(Depolarization_W).append("\n");

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

		println("\n1:runnum, 2:evnum, 3:helicity, 4:e_p, 5:e_theta, 6:e_phi, 7:vz_e,"+
		"8:p_p, 9:p_theta, 10:p_phi, 11:vz_p, 12:Q2, 13:W, 14:Mx, 15: Mx2, 16:x, 17:y, 18:z,"+
		"19:xF, 20:pT, 21:zeta, 22:eta, 23:phi (trento), "+
		"24:DepA, 25:DepB, 26:DepC, 27:DepV, 28:DepW\n");

		println("Set p1 PID = $p1_Str");
		println("output file is: $file");
	}

	writer.close();

	// End time
	long endTime = System.currentTimeMillis()
	// Calculate the elapsed time
	long elapsedTime = endTime - startTime
	// Print the elapsed time in milliseconds
	println("Elapsed time: ${elapsedTime} ms");

}