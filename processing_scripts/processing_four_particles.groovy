/*
 * author Timothy B. Hayward
 * 
 * SIDIS trihadron 
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

	// Set the PDG PID for p1 based on the provided 2nd argument or default to 211 (pi+)
	String p1_Str = args.length < 2 ? "211" : args[1];
	if (args.length < 2) println("WARNING: Specify a PDG PID for p1! Set to pi+ (211).");
	println("Set p1 PID = $p1_Str");
	int p1_int = p1_Str.toInteger(); // Convert p1_Str to integer

	// Set the PDG PID for p2 based on the provided 3rd argument or default to -211 (pi-)
	String p2_Str = args.length < 3 ? "-211" : args[2];
	if (args.length < 3) println("WARNING: Specify a PDG PID for p2! Set to pi- (-211).");
	println("Set p2 PID = $p2_Str");
	int p2_int = p2_Str.toInteger(); // Convert p2_Str to integer

	// Set the PDG PID for p3 based on the provided 4th argument or default to 2212 (p)
	String p3_Str = args.length < 4 ? "2212" : args[3];
	if (args.length < 4) println("WARNING: Specify a PDG PID for p3! Set to proton (2212).");
	println("Set p3 PID = $p3_Str");
	int p3_int = p3_Str.toInteger(); // Convert p3_Str to integer

	// Set the output file name based on the provided 5th argument or use the default name
	String output_file = args.length < 5 ? "hadron_dummy_out.txt" : args[4];
	if (args.length < 5) 
	    println("WARNING: Specify an output file name. Set to \"trihadron_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

	// Set the number of files to process based on the provided 6th argument
	// If the argument is "0", default to the full list size
	int n_files = args.length < 6 || Integer.parseInt(args[5]) == 0 || Integer.parseInt(args[5]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[5]);
	if (args.length < 6 || Integer.parseInt(args[5]) == 0 || Integer.parseInt(args[5]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified, set to 0, or too large
	    println("WARNING: Number of files not specified, set to 0, or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

	// Set the beam energy based on the provided 7th argument or default to 10.6
	double beam_energy = args.length < 7 ? 10.6 : Double.parseDouble(args[6]);
	if (args.length < 7) {
	    println("No beam energy provided, defaulting to 10.6 GeV.");
	}

	// Set the user-provided run number if available
	Integer userProvidedRun = null
	if (args.length < 8) {
	    println("Run number not provided, will pull from hipo files.")
	    println("Think carefully about this if you are processing MC.")
	} else {
		userProvidedRun = Integer.parseInt(args[7]);
	}

	// ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

	// declare physics event variables
	int helicity, detector1, detector2, detector3;
	double e_p, e_theta, e_phi; 
	double p1_phi, p1_p, p1_theta, p2_phi, p2_p, p2_theta, p3_phi, p3_p, p3_theta; 
	double vz_e, vz_p1, vz_p2, vz_p3;
	double open_angle_ep, open_angle_ep1, open_angle_ep2, open_angle_ep3;
	double open_angle_p1p2, open_angle_p1p3, open_angle_p2p3;
	double Q2, W, y, Mx2, Mx2_1, Mx2_2, Mx2_3, Mx2_12, Mx2_13, Mx2_23; 
	double x, t, t1, t2, t3, t12, t13, t23, tmin, z, xF, pT, eta, eta_gN, xi;
	double z1, z2, z3, z12, z13, z23;
	double xF1, xF2, xF3, xF12, xF13, xF23;
	double Mh, Mh12, Mh13, Mh23;
	double pT1, pT2, pT3, pT12, pT13, pT23, pTpT;
	double eta1, eta2, eta3, eta12, eta13, eta23;
	double phi1, phi2, phi3, phi12, phi13, phi23, Delta_phi12, Delta_phi13, Deltaphi23;
	double phih, phiR, theta;
	double Depolarization_A, Depolarization_B, Depolarization_C;
	double Depolarization_V, Depolarization_W;

	// load my kinematic fitter/PID
	GenericKinematicFitter fitter = new analysis_fitter(10.6041); 
	// GenericKinematicFitter fitter = new monte_carlo_fitter(10.6041);
	// GenericKinematicFitter fitter = new event_builder_fitter(10.6041); 
	
	// set filter for final states
	EventFilter filter = new EventFilter("11:"+p1_Str+":"+p2_Str+":"+p3_Str+":X+:X-:Xn");
	// EventFilter filter = new EventFilter("11:"+p1_Str+":"+p2_Str+":"+p3_Str+":Xn");
	
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
		    if (runnum > 16600 && runnum < 16700) break; // Hall C bleedthrough
		    int evnum = event.getBank("RUN::config").getInt('event', 0);

		    PhysicsEvent research_Event = fitter.getPhysicsEvent(event);

		    // do not use the qa if it is MC (runnum = 11) 
		    // do not use the qa if the run is from RGC (until QA is produced!)
		    // boolean process_event = filter.isValid(research_Event)
		    // boolean process_event = filter.isValid(research_Event) && 
		    // 	(runnum == 11 || runnum == 16194 || runnum == 16089 || runnum == 16185 ||
	    	// 	runnum == 16308 || runnum == 16184 || runnum == 16307 || runnum == 16309 ||
	    	// 	qa.OkForAsymmetry(runnum, evnum));
	    	boolean process_event = filter.isValid(research_Event) && (runnum == 11 || runnum < 5020 ||
	    	qa.pass(runnum, evnum));
	    	if (runnum > 17768) process_event = false; // outbending RGC Sp23

		    if (process_event) {

		        // get # of particles 
		        int num_p1 = research_Event.countByPid(p1_Str.toInteger());
		        int num_p2 = research_Event.countByPid(p2_Str.toInteger()); 
		        int num_p3 = research_Event.countByPid(p3_Str.toInteger()); 

		        // cycle over all hadrons
		        for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { 
		        	// cycle over all combinations
					for (int current_p2 = 0; current_p2 < num_p2; current_p2++) {
						for (int current_p3 = 0; current_p3 < num_p3; current_p3++) { 
			        		if (current_p1 == current_p2 && p1_Str.toInteger() == p2_Str.toInteger()) {continue; }
							if (current_p1 == current_p3 && p1_Str.toInteger() == p3_Str.toInteger()) {continue; }
							if (current_p2 == current_p3 && p2_Str.toInteger() == p3_Str.toInteger()) {continue; }

							// supply runnum and boolean for radiative simulation or not
							BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false);
							// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
							double energy = (runnum == 11) ? beam_energy : Eb.Eb();
				            FourParticles variables = new FourParticles(event, research_Event, p1_int, 
				            	current_p1, p2_int, current_p2, p3_int, current_p3, energy);
				            // this is my class for defining all relevant kinematic variables

				            if (variables.channel_test(variables)) {
				                fiducial_status = variables.get_fiducial_status(); // fiducial_status of track
				                helicity = variables.get_helicity(); // helicity of event
				                detector1 = variables.get_detector1();
				                detector2 = variables.get_detector2();
				                detector3 = variables.get_detector3();
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
				                p3_phi = variables.p3_phi(); // lab azimuthal angle
				                p3_p = variables.p3_p(); // lab momentum
				                p3_theta = variables.p3_theta(); // lab polar angle
				                open_angle_ep = variables.open_angle_ep;
				                open_angle_ep1 = variables.open_angle_ep1;
				                open_angle_ep2 = variables.open_angle_ep2;
				                open_angle_ep3 = variables.open_angle_ep3;
				                open_angle_p1p2 = variables.open_angle_p1p2;
				                open_angle_p1p3 = variables.open_angle_p1p3;
				                open_angle_p2p3 = variables.open_angle_p2p3;

				                // vertices
				                vz_e = variables.vz_e();
				                vz_p1 = variables.vz_p1();
				                vz_p2 = variables.vz_p2();
				                vz_p3 = variables.vz_p3();

				                // DIS variables
				                Q2 = variables.Q2(); // exchanged virtual photon energy
				                W = variables.W(); // hadronic mass
				                x = variables.x(); // Bjorken-x
				                t = variables.t();
				                t1 = variables.t1();
				                t2 = variables.t2();
				                t3 = variables.t3();
				                t12 = variables.t12();
				                t13 = variables.t13();
				                t23 = variables.t23();
				                tmin = variables.tmin();
				                y = variables.y(); // E_scat/E_beam
				                Mx2 = variables.Mx2(); // missing mass squared
				                Mx2_1 = variables.Mx2_1(); // missing mass squared calculated with p1
				                Mx2_2 = variables.Mx2_2(); // missing mass squared calculated with p2
				                Mx2_3 = variables.Mx2_3(); // missing mass squared calculated with p3
				                Mx2_12 = variables.Mx2_12(); // missing mass squared calculated with p1 and p2
				                Mx2_13 = variables.Mx2_13(); // missing mass squared calculated with p1 and p3
				                Mx2_23 = variables.Mx2_23(); // missing mass squared calculated with p2 and p3
				                
				                // SIDIS variables
				                z = variables.z(); // fractional hadron energy wrt virtual photon
				                xF = variables.xF(); // Feynman-x
				                pT = variables.pT(); // transverse momentum of hadron
				                eta = variables.eta(); // rapidity
				                eta_gN = variables.eta_gN();
				                xi = variables.xi(); // longitudinal momentum of hadron (fracture functions)

				                // SIDIS trihadron variables
								z1 = variables.z1();
								z2 = variables.z2();
								z3 = variables.z3();
								z12 = variables.z12();
								z13 = variables.z13();
								z23 = variables.z23();

								xF1 = variables.xF1();
								xF2 = variables.xF2();
								xF3 = variables.xF3();
								xF12 = variables.xF12();
								xF13 = variables.xF13();
								xF23 = variables.xF23();

								xi1 = variables.xi1();
								xi2 = variables.xi2();
								xi3 = variables.xi3(); 
								xi12 = variables.xi12();
								xi13 = variables.xi13();
								xi23 = variables.xi23();

								Mh = variables.Mh();
								Mh12 = variables.Mh12();
								Mh13 = variables.Mh13();
								Mh23 = variables.Mh23();

								pT1 = variables.pT1();
								pT2 = variables.pT2();
								pT3 = variables.pT3();
								pT12 = variables.pT12();
								pT13 = variables.pT13();
								pT23 = variables.pT23();
								pTpT = variables.pTpT();

								eta1 = variables.eta1();
								eta2 = variables.eta2();
								eta3 = variables.eta3();
								eta12 = variables.eta12();
								eta13 = variables.eta13();
								eta23 = variables.eta23();

				                // angles
								phi1 = variables.phi1(); // trento phi of the p1
								phi2 = variables.phi2(); // trento phi of the p2
								phi3 = variables.phi3(); // trento phi of the p3
								phi12 = variables.phi12(); // trento phi of the p12
								phi13 = variables.phi13(); // trento phi of the p13
								phi23 = variables.phi23(); // trento phi of the p23
								Delta_phi12 = variables.Delta_phi12();
								Delta_phi13 = variables.Delta_phi13();
								Delta_phi23 = variables.Delta_phi23();
								phih = variables.phih(); // trento phi of the dihadron p1p2
								phiR = variables.phiR(); // second azimuthal angle of dihadron p1p2
								theta = variables.theta(); // decay angle of dihadron

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
									.append(detector1).append(" ")
									.append(detector2).append(" ")
									.append(detector3).append(" ")
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
									.append(p3_p).append(" ")
									.append(p3_theta).append(" ")
									.append(p3_phi).append(" ")
									.append(vz_p3).append(" ")
									.append(open_angle_ep).append(" ")
									.append(open_angle_ep1).append(" ")
									.append(open_angle_ep2).append(" ")
									.append(open_angle_ep3).append(" ")
									.append(open_angle_p1p2).append(" ")
									.append(open_angle_p1p3).append(" ")
									.append(open_angle_p2p3).append(" ")
									.append(Q2).append(" ")
									.append(W).append(" ")
									.append(Mx2).append(" ")
									.append(Mx2_1).append(" ")
									.append(Mx2_2).append(" ")
									.append(Mx2_3).append(" ")
									.append(Mx2_12).append(" ")
									.append(Mx2_13).append(" ")
									.append(Mx2_23).append(" ")
									.append(x).append(" ")
									.append(t).append(" ")
				                	.append(t1).append(" ")
				                	.append(t2).append(" ")
				                	.append(t3).append(" ")
				                	.append(t12).append(" ")
				                	.append(t13).append(" ")
				                	.append(t23).append(" ")
				                	.append(tmin).append(" ")
									.append(y).append(" ")
									.append(z).append(" ")
									.append(z1).append(" ")
									.append(z2).append(" ")
									.append(z3).append(" ")
									.append(z12).append(" ")
									.append(z13).append(" ")
									.append(z23).append(" ")
									.append(xi).append(" ")
									.append(xi1).append(" ")
									.append(xi2).append(" ")
									.append(xi3).append(" ")
									.append(xi12).append(" ")
									.append(xi13).append(" ")
									.append(xi23).append(" ")
									.append(pT).append(" ")
									.append(pT1).append(" ")
									.append(pT2).append(" ")
									.append(pT3).append(" ")
									.append(pT12).append(" ")
									.append(pT13).append(" ")
									.append(pT23).append(" ")
									.append(Mh).append(" ")
									.append(Mh12).append(" ")
									.append(Mh13).append(" ")
									.append(Mh23).append(" ")
									.append(xF).append(" ")
									.append(xF1).append(" ")
									.append(xF2).append(" ")
									.append(xF3).append(" ")
									.append(xF12).append(" ")
									.append(xF13).append(" ")
									.append(xF23).append(" ")
									.append(eta).append(" ")
									.append(eta1).append(" ")
									.append(eta2).append(" ")
									.append(eta3).append(" ")
									.append(eta12).append(" ")
									.append(eta13).append(" ")
									.append(eta23).append(" ")
									.append(phi1).append(" ")
									.append(phi2).append(" ")
									.append(phi3).append(" ")
									.append(phi12).append(" ")
									.append(phi13).append(" ")
									.append(phi23).append(" ")
									.append(phih).append(" ")
									.append(phiR).append(" ")
									.append(theta).append(" ")
									.append(Delta_phi12).append(" ")
									.append(Delta_phi13).append(" ")
									.append(Delta_phi23).append(" ")
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
		        }
		    }
		reader.close();
		}

		// Write any remaining lines in the batchLines StringBuilder to the file
		if (batchLines.length() > 0) {
		    file.append(batchLines.toString());
		    batchLines.setLength(0);
		}

		println();
		print("1: fiducial_status, 2: num_pos, 3: num_neg, 4: num_neutrals, ");
		print("5: runnum, 6: evnum, 7: helicity, ");
		print("8: detector1, 9: detector2, 10: detector3, ");
		print("11: e_p, 12: e_theta, 13: e_phi, 14: vz_e, ");
		print("15: p1_p, 16: p1_theta, 17: p1_phi, 18: vz_p1, ");
		print("19: p2_p, 20: p2_theta, 21: p2_phi, 22: vz_p2, ");
		print("23: p3_p, 24: p3_theta, 25: p3_phi, 26: vz_p3, ");
		print("27: open_angle_ep, 28: open_angle_ep1, 29: open_angle_ep2, 30: open_angle_ep3, ");
		print("31: open_angle_p1p2, 32: open_angle_p1p3, 33: open_angle_p2p3, ");
		print("34: Q2, 35: W, ");
		print("36: Mx2, 37: Mx2_1, 38: Mx2_2, 39: Mx2_3, 40: Mx2_12, 41: Mx2_13, 42: Mx2_23, ");
		print("43: x, 44: t, 45: t1, 46: t2, 47: t3, 48: t12, 49: t13, 50: t23, 51: tmin, ");
		print("52: y, 53: z, 54: z1, 55: z2, 56: z3, 57: z12, 58: z13, 59: z23, ");
		print("60: xi, 61: xi1, 62: xi2, 63: xi3, 64: xi12, 65: xi13, 66: xi23, ");
		print("67: pT, 68: pT1, 69: pT2, 70: pT3, 71: pT12, 72: pT13, 73: pT23, ");
		print("74: Mh, 75: Mh12, 76: Mh13, 77: Mh23, ");
		print("78: xF, 79: xF1, 80: xF2, 81: xF3, 82: xF12, 83: xF13, 84: xF23, ");
		print("85: eta, 86: eta1, 87: eta2, 88: eta3, 89: eta12, 90: eta13, 91: eta23, ");
		print("92: phi1, 93: phi2, 94: phi3, 95: phi12, 96: phi13, 97: phi23, 98: phih, 99: phiR, 100: theta, ");
		print("101: Delta_phi12, 102: Delta_phi13, 103: Delta_phi23, ");
		print("104: DepA, 105: DepB, 106: DepC, 107: DepV, 108: DepW.");
		println();
		println("Set p1 PID = $p1_Str");
		println("Set p2 PID = $p2_Str");
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