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

	// Set the output file name based on the provided 2nd argument or use the default name
	String output_file = args.length < 2 ? "hadron_dummy_out.txt" : args[1];
	if (args.length < 2) 
	    println("WARNING: Specify an output file name. Set to \"hadron_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

	// Set the number of files to process based on the provided 3rd argument
	// If the argument is "0", default to the full list size
	int n_files = args.length < 3 || Integer.parseInt(args[2]) == 0 || Integer.parseInt(args[2]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[2]);
	if (args.length < 3 || Integer.parseInt(args[2]) == 0 || Integer.parseInt(args[2]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified, set to 0, or too large
	    println("WARNING: Number of files not specified, set to 0, or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

	// Set the beam energy based on the provided 4th argument or default to 10.6
	double beam_energy = args.length < 4 ? 10.6 : Double.parseDouble(args[3]);
	if (args.length < 4) {
	    println("No beam energy provided, defaulting to 10.6 GeV.");
	}

	// Set the user-provided run number if available
	Integer userProvidedRun = null
	if (args.length < 5) {
	    println("Run number not provided, will pull from hipo files.")
	    println("Think carefully about this if you are processing MC.")
	} else {
		userProvidedRun = Integer.parseInt(args[4]);
	}

	// ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

	// declare physics event variables
	int helicity, detector1, detector2, detector_gamma1, detector_gamma2;
	double e_p, e_theta, e_phi, p1_phi, p1_p, p1_theta, p2_phi, p2_p, p2_theta; 
	double vz_e, vz_p1, vz_p2;
	double open_angle_ep, open_angle_ep1, open_angle_ep2, open_angle_p1p2;
	double open_angle_egamma1, open_angle_egamma2, gamma_phi1, gamma_phi2;
	double Q2, W, y, Mx2, Mx2_1, Mx2_2; 
	double x, t, t1, t2, tmin, z, xF, pT, eta, eta_gN, xi;
	double z1, z2, xF1, xF2, Mh, Mh_gammagamma;
	double pT1, pT2, pTpT, eta1, eta2, Delta_eta, eta1_gN, eta2_gN;
	double phi1, phi2, Delta_phi, phih, phiR, theta;
	double Depolarization_A, Depolarization_B, Depolarization_C;
	double Depolarization_V, Depolarization_W;
	double Emiss2, theta_pi0_pi0, pTmiss;

	// load my kinematic fitter/PID
	// Initialize the fitter
	GenericKinematicFitter fitter;

	// Uncomment the desired fitter
	fitter = new analysis_fitter(10.6041); 
	// fitter = new monte_carlo_fitter(10.6041);
	// fitter = new event_builder_fitter(10.6041);  

	// Set filter for final states based on fitter type
	EventFilter filter;

	if (fitter instanceof analysis_fitter) {
	    // Use filter with 111 for analysis fitter
	    filter = new EventFilter("11:2212:111:22:22:Xn");
	} else if (fitter instanceof monte_carlo_fitter) {
	    // Use filter without 111 for monte_carlo_fitter
	    filter = new EventFilter("11:2212:22:22:Xn");
	} else {
	    // Default filter or another fitter type
	    filter = new EventFilter("11:2212:111:22:22:Xn");
	}
	
	// setup QA database
	QADB qa = new QADB();
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
	[ // list of runs with `Misc` that should be allowed, generally empty target etc for dilution factor calculations
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
		    ++num_events;
		    if (num_events % 1000000 == 0) { // not necessary, just updates output
		        print("processed: " + num_events + " events. ");
		    }

		    // get run and event numbers
		    event = reader.getNextEvent();
		    // collect info for QA
		    int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0);
		    int evnum = event.getBank("RUN::config").getInt('event', 0);

		    PhysicsEvent research_Event = fitter.getPhysicsEvent(event);

		    // do not use the qa if it is MC (runnum = 11) 
		    // do not use the qa if the run is from RGC (until QA is produced!)
		    // boolean process_event = filter.isValid(research_Event);
		    // boolean process_event = filter.isValid(research_Event) && 
		    // 	(runnum == 11 || runnum == 16194 || runnum == 16089 || runnum == 16185 ||
	    	// 	runnum == 16308 || runnum == 16184 || runnum == 16307 || runnum == 16309 ||
	    	// 	qa.OkForAsymmetry(runnum, evnum));
	    	boolean process_event = filter.isValid(research_Event) && (runnum == 11 || runnum < 5020 ||
	    	qa.pass(runnum, evnum));
	    	if (runnum > 17768) process_event = false; // outbending RGC Sp23
	    	if (runnum == 17331 || runnum == 16987 || runnum == 17079) process_event = false; // low live time
	    	
		    if (process_event) {


		    	int num_photons = research_Event.countByPid(22);
		    	for (int current_gamma1 = 0; current_gamma1 < num_photons; current_gamma1++) {
		    		for (int current_gamma2 = 0; current_gamma2 < num_photons; current_gamma2++) {
		    			if (current_gamma1 == current_gamma2) continue;

		    			// supply runnum and boolean for radiative simulation or not
						BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false);
						// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
						double energy = (runnum == 11) ? beam_energy : Eb.Eb();
			            ThreeParticles variables = new ThreeParticles(event, research_Event, 
							22, current_gamma1, 22, current_gamma2, energy);

			            Mh_gammagamma = variables.Mh();
			            if (Mh_gammagamma < 0.11 || Mh_gammagamma > 0.16) continue;
			            if (detector_gamma1 == 2 || detector_gamma2 == 2) continue;
			            detector_gamma1 = variables.get_detector1();
			            detector_gamma2 = variables.get_detector2();
			            open_angle_egamma1 = variables.open_angle_ep1();
			            open_angle_egamma2 = variables.open_angle_ep2();
			            gamma_phi1 = variables.phi1();
			            gamma_phi2 = variables.phi2();
		    		}
		    	}
		    	if (Mh_gammagamma < 0.11 || Mh_gammagamma > 0.16) continue;
		    	if (detector_gamma1 == 2 || detector_gamma2 == 2) continue;

        		// supply runnum and boolean for radiative simulation or not
				BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false);
				// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
				double energy = (runnum == 11) ? beam_energy : Eb.Eb();
	            ThreeParticles variables = new ThreeParticles(event, research_Event, 
					2212, 0, 111, 0, energy);
	            // this is my class for defining all relevant kinematic variables

	            if (variables.channel_test(variables)) {

	                fiducial_status = variables.get_fiducial_status(); // fiducial_status of track
	                helicity = variables.get_helicity(); // helicity of event
	                detector1 = variables.get_detector1(); 
	                detector2 = detector_gamma1; // the pi0 isn't actually detected in detector
	                // so put in the location of the photons it's formed from 
	                num_pos = variables.get_num_pos();
	                num_neg = variables.get_num_neg();
	                num_neutrals = variables.get_num_neutrals();

	                // lab kinematics
	                e_p = variables.e_p(); // lab frame momentum
	                e_theta = variables.e_theta(); // lab polar angle
	                e_phi = variables.e_phi(); // lab azimuthal angle
	                p1_phi = variables.p1_phi(); // lab azimuthal angle
	                p1_p = variables.p1_p(); // lab momentum
	                p1_theta = variables.p1_theta(); // lab polar angle
	                p2_phi = variables.p2_phi(); // lab azimuthal angle
	                p2_p = variables.p2_p(); // lab momentum
	                p2_theta = variables.p2_theta(); // lab polar angle
	                open_angle_ep = variables.open_angle_ep(); 
	                open_angle_ep1 = variables.open_angle_ep1(); 
	                open_angle_ep2 = variables.open_angle_ep2(); 
	                open_angle_p1p2 = variables.open_angle_p1p2(); 

	                // vertices
	                vz_e = variables.vz_e();
	                vz_p1 = variables.vz_p1();
	                vz_p2 = variables.vz_p2();

	                // DIS variables
	                Q2 = variables.Q2(); // exchanged virtual photon energy
	                W = variables.W(); // hadronic mass
	                x = variables.x(); // Bjorken-x
	                t = variables.t();
	                t1 = variables.t1();
	                t2 = variables.t2();
	                tmin = variables.tmin();
	                y = variables.y(); // E_scat/E_beam
	                Mx2 = variables.Mx2(); // missing mass squared
	                Mx2_1 = variables.Mx2_1(); // missing mass squared calculated with p1
	                Mx2_2 = variables.Mx2_2(); // missing mass squared calculated with p2

	                // SIDIS variables
	                z = variables.z(); // fractional hadron energy wrt virtual photon
	                xF = variables.xF(); // Feynman-x
	                pT = variables.pT(); // transverse momentum of hadron
	                eta = variables.eta(); // rapidity
	                eta_gN = variables.eta_gN();
	                xi = variables.xi(); // longitudinal momentum of hadron (fracture functions)

	                // SIDIS dihadron variables
					z1 = variables.z1();
					z2 = variables.z2();
					xF1 = variables.xF1();
					xF2 = variables.xF2();
					xi1 = variables.xi1();
					xi2 = variables.xi2(); 
					Mh = variables.Mh();
					pT1 = variables.pT1();
					pT2 = variables.pT2();
					pTpT = variables.pTpT();
					eta1 = variables.eta1();
					eta2 = variables.eta2();
					Delta_eta = variables.Delta_eta();
					eta1_gN = variables.eta1_gN();
					eta2_gN = variables.eta2_gN();

	                // angles
					phi1 = variables.phi1(); // trento phi of the p1
					phi2 = variables.phi2(); // trento phi of the p2
					Delta_phi = variables.Delta_phi();
					phih = variables.phih(); // trento phi of the dihadron
					phiR = variables.phiR(); // second azimuthal angle of dihadron
					theta = variables.theta(); // decay angle of dihadron

	                // depolarization factors
	                Depolarization_A = variables.Depolarization_A();
	                Depolarization_B = variables.Depolarization_B();
	                Depolarization_C = variables.Depolarization_C();
	                Depolarization_V = variables.Depolarization_V();
			    	Depolarization_W = variables.Depolarization_W();

			    	// exclusivity variables
			    	Emiss2 = variables.Emiss2();
			    	theta_pi0_pi0 = variables.theta_gamma_gamma();
			    	pTmiss = variables.pTmiss();

	                // Use a StringBuilder to append all data in a single call
	                StringBuilder line = new StringBuilder();
	                line.append(fiducial_status).append(" ")
						.append(num_pos).append(" ")
						.append(num_neg).append(" ")
						.append(num_neutrals).append(" ")
						.append(runnum).append(" ")
	                	.append(evnum).append(" ")
	                	.append(helicity).append(" ")
	                	.append(detector1).append(" ")
	                	.append(detector2).append(" ")
	                	.append(e_p).append(" ")
	                	.append(e_theta).append(" ")
	                	.append(e_phi).append(" ")
	                	.append(vz_e).append(" ")
	                	.append(p1_p).append(" ")
	                	.append(p1_theta).append(" ")
	                	.append(p1_phi).append(" ")
	                	.append(vz_p1).append(" ")
	                	.append(p2_p).append(" ")
	                	.append(p2_theta).append(" ")
	                	.append(p2_phi).append(" ")
	                	.append(vz_p2).append(" ")
	                	.append(open_angle_ep).append(" ")
	                	.append(open_angle_ep1).append(" ")
	                	.append(open_angle_ep2).append(" ")
	                	.append(open_angle_p1p2).append(" ")
	                	.append(Q2).append(" ")
	                	.append(W).append(" ")
	                	.append(Mx2).append(" ")
	                	.append(Mx2_1).append(" ")
	                	.append(Mx2_2).append(" ")
	                	.append(x).append(" ")
	                	.append(t).append(" ")
	                	.append(t1).append(" ")
	                	.append(t2).append(" ")
	                	.append(tmin).append(" ")
	                	.append(y).append(" ")
	                	.append(z).append(" ")
	                	.append(z1).append(" ")
	                	.append(z2).append(" ")
	                	.append(Mh).append(" ")
	                	.append(xF).append(" ")
	                	.append(xF1).append(" ")
	                	.append(xF2).append(" ")
	                	.append(pT).append(" ")
	                	.append(pT1).append(" ")
	                	.append(pT2).append(" ")
	                	.append(pTpT).append(" ")
	                	.append(xi).append(" ")
	                	.append(xi1).append(" ")
	                	.append(xi2).append(" ")
	                	.append(eta).append(" ")
	                	.append(eta1).append(" ")
	                	.append(eta2).append(" ")
	                	.append(Delta_eta).append(" ")
	                	.append(eta1_gN).append(" ")
	                	.append(eta2_gN).append(" ")
	                	.append(phi1).append(" ")
	                	.append(phi2).append(" ")
	                	.append(Delta_phi).append(" ")
	                	.append(phih).append(" ")
	                	.append(phiR).append(" ")
	                	.append(theta).append(" ")
	                	.append(Depolarization_A).append(" ")
	                    .append(Depolarization_B).append(" ")
	                    .append(Depolarization_C).append(" ")
	                    .append(Depolarization_V).append(" ")
	                    .append(Depolarization_W).append(" ")
	                    .append(Mh_gammagamma).append(" ")
	                    .append(detector_gamma1).append(" ")
	                    .append(detector_gamma2).append(" ")
	                    .append(open_angle_egamma1).append(" ")
	                    .append(open_angle_egamma2).append(" ")
	                    .append(gamma_phi1).append(" ")
	                    .append(gamma_phi2).append(" ")
	                    .append(Emiss2).append(" ")
	                    .append(theta_pi0_pi0).append(" ")
	                    .append(pTmiss).append("\n");

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
		reader.close();
		}

		// Write any remaining lines in the batchLines StringBuilder to the file
		if (batchLines.length() > 0) {
		    file.append(batchLines.toString());
		    batchLines.setLength(0);
		}

		println("1: fiducial_status, 2: num_pos, 3: num_neg, 4: num_neutrals, " +
	    "5: runnum, 6: evnum, 7: helicity, 8: detector1, 9: detector2, 10: e_p, 11: e_theta, 12: e_phi, 13: vz_e, " +
	    "14: p1_p, 15: p1_theta, 16: p1_phi, 17: vz_p1, 18: p2_p, 19: p2_theta, 20: p2_phi, 21: vz_p2, " +
	    "22: open_angle_ep, 23: open_angle_ep1, 24: open_angle_ep2, 25: open_angle_p1p2, " +
	    "26: Q2, 27: W, 28: Mx2, 29: Mx2_1, 30: Mx2_2, 31: x, 32: t, 33: t1, 34: t2, 35: tmin, 36: y, 37: z, " +
	    "38: z1, 39: z2, 40: Mh, 41: xF, 42: xF1, 43: xF2, 44: pT, 45: pT1, 46: pT2, 47: pTpT, " +
	    "48: xi, 49: xi1, 50: xi2, 51: eta, 52: eta1, 53: eta2, 54: Delta_eta, 55: eta1_gN, 56: eta2_gN, " +
	    "57: phi1, 58: phi2, 59: Delta_phi, 60: phih, 61: phiR, 62: theta, " +
	    "63: DepA, 64: DepB, 65: DepC, 66: DepV, 67: DepW, 68: Mh_gammagamma, " +
	    "69: detector_gamma1, 70: detector_gamma2, 71: open_angle_egamma1, 72: open_angle_egamma2, 73: gamma_phi1, 74: gamma_phi2" +
	    "75: Emiss2, 76: theta_pi0_pi0, 77: pTmiss.");

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