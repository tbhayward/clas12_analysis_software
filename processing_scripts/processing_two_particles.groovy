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

// dilks CLAS QA analysis
import clasqa.QADB // access QADB

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

	// Start time
	long startTime = System.currentTimeMillis();

	// ~~~~~~~~~~~~~~~~ set up input parameters ~~~~~~~~~~~~~~~~ //

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
	// If the argument is "0", default to the full list size
	int n_files = args.length < 4 || Integer.parseInt(args[3]) == 0 || Integer.parseInt(args[3]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[3]);
	if (args.length < 4 || Integer.parseInt(args[3]) == 0 || Integer.parseInt(args[3]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified, set to 0, or too large
	    println("WARNING: Number of files not specified, set to 0, or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

	// Set the beam energy based on the provided 5th argument or default to 10.6
	double beam_energy = args.length < 5 ? 10.6 : Double.parseDouble(args[4]);
	if (args.length < 5) {
	    println("No beam energy provided, defaulting to 10.6 GeV.");
	    println("All MC will use 10.604 GeV. You must manually enter a beam energy to change this.")
	}

	// Set the user-provided run number if available
	Integer userProvidedRun = null
	if (args.length < 6) {
	    println("Run number not provided, will pull from hipo files.")
	    println("Think carefully about this if you are processing MC.")
	} else {
		userProvidedRun = Integer.parseInt(args[5]);
	}

	// ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

	// declare physics event variables
	int helicity, detector;
	double e_p, e_theta, e_phi, p_phi, p_p, p_theta, open_angle;
	double Q2, W, y, Mx2, x, t, tmin, z, xF, pT, eta, xi, phi, vz_e, vz_p;
	double Depolarization_A, Depolarization_B, Depolarization_C;
	double Depolarization_V, Depolarization_W;

	// load my kinematic fitter/PID
	GenericKinematicFitter fitter = new analysis_fitter(10.6041); 
	// GenericKinematicFitter fitter = new monte_carlo_fitter(10.6041);
	// GenericKinematicFitter fitter = new event_builder_fitter(10.6041);  
	
	// set filter for final states
	EventFilter filter = new EventFilter("11:"+p1_Str+":X+:X-:Xn");
	// EventFilter filter = new EventFilter("11:"+p1_Str+":Xn");
	
	// instantiate QADB
	QADB qa = new QADB()
	qa.checkForDefect('TotalOutlier')    
	qa.checkForDefect('TerminalOutlier')
	qa.checkForDefect('MarginalOutlier')
	qa.checkForDefect('SectorLoss')
	qa.checkForDefect('LowLiveTime')
	qa.checkForDefect('Misc')
	qa.checkForDefect('ChargeHigh')
	qa.checkForDefect('ChargeNegative')
	qa.checkForDefect('ChargeUnknown')
	qa.checkForDefect('PossiblyNoBeam')
	[ // list of runs with `Misc` that should be allowed, 
	 // for RGC generally empty target etc for dilution factor calculations
	 	5046, 5047, 5051, 5128, 5129, 5130, 5158, 5159,
  		5160, 5163, 5165, 5166, 5167, 5168, 5169, 5180,
  		5181, 5182, 5183, 5400, 5448, 5495, 5496, 5505,
  		5567, 5610, 5617, 5621, 5623, 6736, 6737, 6738,
  		6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
  		6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
  		6757, 16194, 16089, 16185, 16308, 16184, 16307, 16309
	].each{ run -> qa.allowMiscBit(run) }

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
			// if (num_events==10000) break;
		    ++num_events;
		    if (num_events % 500000 == 0) { // not necessary, just updates output
		        print("processed: " + num_events + " events. ");
		    }

		    // get run and event numbers
		    event = reader.getNextEvent();
		    // collect info for QA
		    int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0);
		    if (runnum > 16600 && runnum < 16700) break; // Hall C bleedthrough
		    int evnum = event.getBank("RUN::config").getInt('event', 0); 

		    PhysicsEvent research_Event = fitter.getPhysicsEvent(event);

		    HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");

		    // do not use the qa if it is MC (runnum = 11) 
		    // boolean process_event = filter.isValid(research_Event) && 
		    // 	(runnum == 11 || runnum < 5020 || runnum >= 11571 || 
		    // 	qa.OkForAsymmetry(runnum, evnum));
		    // boolean process_event = filter.isValid(research_Event) && 
		    // 	(runnum == 11 || runnum == 16194 || runnum == 16089 || runnum == 16185 ||
	    	// 	runnum == 16308 || runnum == 16184 || runnum == 16307 || runnum == 16309 ||
	    	// 	qa.OkForAsymmetry(runnum, evnum));
	    	boolean process_event = filter.isValid(research_Event) && (runnum == 11 || runnum < 5020 || runnum > 16772 ||
	    		qa.pass(runnum, evnum));
	    	if (runnum > 17768) process_event = false; // outbending RGC Sp23
	    	if (runnum == 17331) process_event = false; // low live time

		    if (process_event) {

		        // get # of particles w/ pid1
		        int num_p1 = research_Event.countByPid(p1_Str.toInteger()); 

		        // cycle over all hadrons
		        for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { 

		        	// supply runnum and boolean for radiative simulation or not
					BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false);
					// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
					double energy = (runnum == 11) ? beam_energy : Eb.Eb();
		            TwoParticles variables = new TwoParticles(event, research_Event,
		                    p1_int, current_p1, energy);
		            // this is my class for defining all relevant kinematic variables

		            if (variables.channel_test(variables)) {
		                fiducial_status = variables.get_fiducial_status(); // fiducial_status of track
		                helicity = variables.get_helicity(); // helicity of event
		                detector = variables.get_detector();
		                num_pos = variables.get_num_pos();
		                num_neg = variables.get_num_neg();
		                num_neutrals = variables.get_num_neutrals();

		                // lab kinematics
		                e_p = variables.e_p(); // lab frame momentum
		                e_theta = variables.e_theta(); // lab polar angle
		                e_phi = variables.e_phi(); // lab azimuthal angle
		                p_phi = variables.p_phi(); // lab azimuthal angle
		                p_p = variables.p_p(); // lab momentum
		                p_theta = variables.p_theta(); // lab polar angle
		                open_angle = variables.open_angle(); // angle between electron and second particle

		                // DIS variables
		                Q2 = variables.Q2(); // exchanged virtual photon energy
		                W = variables.W(); // hadronic mass
		                x = variables.x(); // Bjorken-x
		                t = variables.t(); // t
		                tmin = variables.tmin(); // tmin
		                y = variables.y(); // E_scat/E_beam
		                Mx2 = variables.Mx2(); // missing mass square

		                // SIDIS variables
		                z = variables.z(); // fractional hadron energy wrt virtual photon
		                xF = variables.xF(); // Feynman-x
		                pT = variables.pT(); // transverse momentum of hadron
		                eta = variables.eta(); // rapidity
		                xi = variables.xi(); // fractional longitudinal momentum of hadron

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
		                line.append(fiducial_status).append(" ")
							.append(num_pos).append(" ")
							.append(num_neg).append(" ")
							.append(num_neutrals).append(" ")
							.append(runnum).append(" ")
		                	.append(evnum).append(" ")
		                	.append(helicity).append(" ")
		                	.append(detector).append(" ")
		                	.append(e_p).append(" ")
		                	.append(e_theta).append(" ")
		                	.append(e_phi).append(" ")
		                	.append(vz_e).append(" ")
		                	.append(p_p).append(" ")
		                	.append(p_theta).append(" ")
		                	.append(p_phi).append(" ")
		                	.append(vz_p).append(" ")
		                	.append(open_angle).append(" ")
		                	.append(Q2).append(" ")
		                	.append(W).append(" ")
		                	.append(Mx2).append(" ")
		                	.append(x).append(" ")
		                	.append(t).append(" ")
		                	.append(tmin).append(" ")
		                	.append(y).append(" ")
		                	.append(z).append(" ")
		                	.append(xF).append(" ")
		                	.append(pT).append(" ")
		                	.append(xi).append(" ")
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

		println("1: fiducial_status, 2: num_pos, 3: num_neg, 4: num_neutrals, " +
		"5: runnum, 6: evnum, 7: helicity, 8: detector, 9: e_p, 10: e_theta, 11: e_phi, 12: vz_e, " +
		"13: open_angle, 14: p_p, 15: p_theta, 16: p_phi, 17: vz_p, 18: Q2, 19: W, 20: Mx2, " +
		"21: x, 22: t, 23: tmin, 24: y, 25: z, 26: xF, 27: pT, 28: xi, 29: eta, 30: phi (trento), " +
		"31: DepA, 32: DepB, 33: DepC, 34: DepV, 35: DepW");

		println("Set p1 PID = $p1_Str");
		println("output text file is: $file");
	}

	writer.close();

	// End time
	long endTime = System.currentTimeMillis()
	// Calculate the elapsed time
	long elapsedTime = endTime - startTime
	// Print the elapsed time in milliseconds
	println("Elapsed time: ${elapsedTime} ms");

}