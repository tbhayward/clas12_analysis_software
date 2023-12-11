/*
 * author Timothy B. Hayward
 * 
 * SIDIS hadron mc (generated and reconstructed variables saved)
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

public static double phi_calculation (double x, double y) {
	// tracks are given with Cartesian values and so must be converted to cylindrical
	double phi = Math.toDegrees(Math.atan2(x,y));
	phi = phi - 90;
	if (phi < 0) {
		phi = 360 + phi;
	}
	phi = 360 - phi;
	return phi;	
}

public static double theta_calculation (double x, double y, double z) {
	// convert cartesian coordinates to polar angle
	double r = Math.pow(Math.pow(x,2)+Math.pow(y,2)+Math.pow(z,2),0.5);
	return (double) (180/Math.PI)*Math.acos(z/r);
}

public static void main(String[] args) {

	double scale = 3.0;

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

	int hadron_pair_counts = 0;
	GenericKinematicFitter research_fitter = new analysis_fitter(10.6041);
	GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6041);
	EventFilter filter = new EventFilter("11:"+p1_Str+":"+":X+:X-:Xn");

	// create a StringBuilder for accumulating lines
	StringBuilder batchLines = new StringBuilder();

	int num_events = 0;
	int max_lines = 100;
	int lineCount = 0;
	for (current_file in 0..<n_files) {
		// limit to a certain number of files defined by n_files
		println("\n Opening file "+Integer.toString(current_file+1)
			+" out of "+n_files+".\n"); 
		
		HipoDataSource reader = new HipoDataSource();
		reader.open(hipo_list[current_file]); // open next hipo file
		HipoDataEvent event = reader.getNextEvent(); 

		while(reader.hasEvent()){
			++num_events; 
			if (num_events%100000 == 0) { 
				print("processed: "+num_events+" events. ");
			}

			// get run and event numbers
			event = reader.getNextEvent();
		    int runnum = event.getBank("RUN::config").getInt('run',0);
		    int evnum = event.getBank("RUN::config").getInt('event',0);

		    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);
		    PhysicsEvent mc_Event = mc_fitter.getPhysicsEvent(event);

			if (filter.isValid(research_Event)) {

				HipoDataBank recBank = (HipoDataBank) event.getBank("REC::Event");
				HipoDataBank particleBank = (HipoDataBank) event.getBank("REC::Particle");
				HipoDataBank lundBank = (HipoDataBank) event.getBank("MC::Lund");
				HipoDataBank mcBank = (HipoDataBank) event.getBank("MC::Particle");


				int num_p1 = research_Event.countByPid(p1_Str.toInteger()); 

				for (int current_p1 = 0; current_p1 < num_p1; current_p1++) {

					Particle exp_e = research_Event.getParticleByPid(11,0);
					Particle exp_p1 = research_Event.getParticleByPid(p1_Str.toInteger(),current_p1);

					BeamEnergy Eb = new BeamEnergy(runnum, false);
					Hadron variables = new Hadron(event, research_Event, 
						p1_Str.toInteger(), current_p1, Eb.Eb());
					Hadron mc_variables = new Hadron(event, mc_Event, 
						p1_Str.toInteger(), current_p1, Eb.Eb());

					if (variables.channel_test(variables)) {

						// lab kinematics data
						double p_p = variables.p_p();
						double p_theta = variables.p_theta();
						double p_phi = variables.p_phi();
						double vz_p = variables.vz_p();

						// lab kinematics MC
						double mc_p_p = mc_variables.p_p();
						double mc_p_theta = mc_variables.p_theta();
						double mc_p_phi = mc_variables.p_phi();
						double mc_vz_p = mc_variables.vz_p();

						matching_p1 = false;

						matching_p1_pid = 0;
						mc_p1_parent_index = 0;
						for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
							int pid = mcBank.getInt("pid", current_part);
							if (matching_p1) { continue; }
							double mc_px = mcBank.getFloat("px", current_part);
							double mc_py = mcBank.getFloat("py", current_part);
							double mc_pz = mcBank.getFloat("pz", current_part);

							double mc_phi = phi_calculation(mc_px, mc_py);
							double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

							double exp_phi = phi_calculation(exp_p1.px(), exp_p1.py());
							double exp_theta = theta_calculation(exp_p1.px(), exp_p1.py(), 
								exp_p1.pz());

							matching_p1 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
								Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
							if (matching_p1) {
								matching_p1_pid = pid;
							}
						}

						// int rich_pid = 0;
						float beta, chi2pid;
						int particle_Index = -1;
						// if (event.hasBank("RICH::Particle")) {
							// HipoDataBank richBank = (HipoDataBank) event.getBank("RICH::Particle");
							for(int current_part = 0; current_part < particleBank.rows(); current_part++) {
								if (particleBank.getInt("vz", current_part) == vz_p) {
									beta = particleBank.getFloat("beta", current_Part);
				            		chi2pid = particleBank.getFloat("chi2pid", current_Part);
				            		particle_Index = current_part;
								}
							}
							// for (int current_Row = 0; current_Row < rich_Bank.rows(); current_Row++) {
				            // // Get the pindex for the current row
				            // int pindex = rich_Bank.getInt("pindex", current_Row);
				            // // Check if the pindex value matches the specified particle
				            // if (pindex == particle_Index) {
				            //     rich_pid = rich_Bank.getInt("best_PID", current_Row);
				            // }
						// }

						// Use a StringBuilder to append all data in a single call
		                StringBuilder line = new StringBuilder();
						line.append(p_p).append(" ")
							.append(mc_p_p).append(" ")
							.append(p_theta).append(" ")
							.append(mc_p_theta).append(" ")
							.append(p_phi).append(" ")
							.append(mc_p_phi).append(" ")
							.append(matching_p1_pid).append(" ")
							.append(beta).append(" ")
							.append(chi2pid).append("\n");
							// .append(rich_pid).append("\n");

						// Append the line to the batchLines StringBuilder
		                batchLines.append(line.toString());
		                lineCount++; // Increment the line count

		                // If the line count reaches 100, write to the file and reset
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

		println(); println();
		print("1: p_p, 2: mc_p_p, 3: p_theta, 4: mc_p_theta, 5: p_phi, 6: mc_p_phi, ");
		print("7: matching_p_pid, 8: beta, 9: chi2pid.\n");

		println(); println();
		println("Set p1 PID = "+p1_Str+"\n");
		println("output file is: "+file);
		println();
	}

}