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

	// Set the PDG PID for p2 based on the provided 3nd argument or default to -211 (pi-)
	String p2_Str = args.length < 3 ? "-211" : args[2];
	if (args.length < 3) println("WARNING: Specify a PDG PID for p2! Set to pi- (-211).");
	println("Set p2 PID = $p2_Str");
	int p2_int = p2_Str.toInteger(); // Convert p2_Str to integer

	// Set the PDG PID for p3 based on the provided 4nd argument or default to 2212 (p)
	String p3_Str = args.length < 4 ? "2212" : args[3];
	if (args.length < 4) println("WARNING: Specify a PDG PID for p2! Set to proton (2212).");
	println("Set p3 PID = $p3_Str");
	int p3_int = p3_Str.toInteger(); // Convert p3_Str to integer

	// Set the output file name based on the provided 4th argument or use the default name
	String output_file = args.length < 5 ? "hadron_dummy_out.txt" : args[4];
	if (args.length < 5) 
	    println("WARNING: Specify an output file name. Set to \"trihadron_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

	// Set the number of files to process based on the provided 5tth argument
	// use the size of the hipo_list if no argument provided
	int n_files = args.length < 6 || Integer.parseInt(args[5]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[5]);
	if (args.length < 6 || Integer.parseInt(args[5]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified or too large
	    println("WARNING: Number of files not specified or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

	// ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

	// declare physics event variables
	int helicity;
	double e_p, e_theta, e_phi; 
	double p1_phi, p1_p, p1_theta, p2_phi, p2_p, p2_theta, p3_phi, p3_p, p3_theta; 
	double vz_e, vz_p1, vz_p2, vz_p3;
	double Q2, W, y, Mx, Mx1, Mx2, Mx3; 
	double x, z, xF, pT, eta, eta_gN, zeta;
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

				            Trihadrons variables = new Trihadrons(event, research_Event, 
								p1_int, current_p1, p2_int, current_p2, p3_int, current_p3);
				            // this is my class for defining all relevant kinematic variables

				            if (variables.channel_test(variables)) {
				                helicity = variables.get_helicity(); // helicity of event

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

				                // vertices
				                vz_e = variables.vz_e();
				                vz_p1 = variables.vz_p1();
				                vz_p2 = variables.vz_p2();
				                vz_p3 = variables.vz_p3();

				                // DIS variables
				                Q2 = variables.Q2(); // exchanged virtual photon energy
				                W = variables.W(); // hadronic mass
				                x = variables.x(); // Bjorken-x
				                y = variables.y(); // E_scat/E_beam
				                Mx = variables.Mx(); // missing mass
				                Mx1 = variables.Mx1(); // missing mass calculated with p1
				                Mx2 = variables.Mx2(); // missing mass calculated with p2
				                Mx3 = variables.Mx3(); // missing mass calculated with p3
				                Mx12 = variables.Mx12(); // missing mass calculated with p1 and p2
				                Mx13 = variables.Mx13(); // missing mass calculated with p1 and p3
				                Mx23 = variables.Mx23(); // missing mass calculated with p2 and p3
				                println(Mx23);
				                // SIDIS variables
				                z = variables.z(); // fractional hadron energy wrt virtual photon
				                xF = variables.xF(); // Feynman-x
				                pT = variables.pT(); // transverse momentum of hadron
				                eta = variables.eta(); // rapidity
				                eta_gN = variables.eta_gN();
				                zeta = variables.zeta(); // longitudinal momentum of hadron (fracture functions)

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

								zeta1 = variables.zeta1();
								zeta2 = variables.zeta2();
								zeta3 = variables.zeta3(); 
								zeta12 = variables.zeta12();
								zeta13 = variables.zeta13();
								zeta23 = variables.zeta23();

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
								line.append(runnum).append(" ")
									.append(evnum).append(" ")
									.append(helicity).append(" ")
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
									.append(Q2).append(" ")
									.append(W).append(" ")
									.append(Mx).append(" ")
									.append(Mx1).append(" ")
									.append(Mx2).append(" ")
									.append(Mx3).append(" ")
									.append(Mx12).append(" ")
									.append(Mx13).append(" ")
									.append(Mx23).append(" ")
									.append(x).append(" ")
									.append(y).append(" ")
									.append(z).append(" ")
									.append(z1).append(" ")
									.append(z2).append(" ")
									.append(z3).append(" ")
									.append(z12).append(" ")
									.append(z13).append(" ")
									.append(z23).append(" ")
									.append(zeta).append(" ")
									.append(zeta1).append(" ")
									.append(zeta2).append(" ")
									.append(zeta3).append(" ")
									.append(zeta12).append(" ")
									.append(zeta13).append(" ")
									.append(zeta23).append(" ")
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
		print("1:runnum, 2:evnum, 3:helicity, ");
		print("4: e_p, 5: e_theta, 6: e_phi, 7: vz_e, ");
		print("8: p1_p, 9: p1_theta, 10: p1_phi, 11: vz_p1, ");
		print("12: p2_p, 13: p2_theta, 14: p2_phi, 15: vz_p2, ");
		print("16: p3_p, 17: p3_theta, 18: p3_phi, 19: vz_p3 ");
		print("20: Q2, 21: W");
		print("22: Mx, 23: Mx1, 24: Mx2, 25: Mx3, 26: Mx12, 27: Mx13, 28: Mx23, ");
		print("29: x, 30: y, ");
		print("31: z, 32: z1, 33: z2, 34: z3, 35: z12, 36: z13, 37: z23, ");
		print("38: zeta, 39: zeta1, 40: zeta2, 41: zeta3, 42: zeta12, 43: zeta13, 44: zeta23, ");
		print("45: pT, 46: pT1, 47: pT2, 48: pT3, 49: pT12, 50: pT13, 51: pT23, ");
		print("52: Mh, 53: Mh12, 54: Mh13, 55: Mh23, ");
		print("56: xF, 57: xF1, 58: xF2, 59: xF3, 60: xF12, 61: xF13, 62: xF23, ");
		print("63: eta, 64: eta1, 65: eta2, 66: eta3, 67: eta12, 68: eta13, 69: eta23, ");
		print("70: phi1, 71: phi2, 72: phi3, 73: phi12, 74: phi13, 75: phi23, 76: phih, 77: phiR, 78: theta, ");
		print("79: Delta_phi12, 80: Delta_phi13, 81: Delta_phi23, ")
		print("82: DepA, 83: DepB, 84: DepC, 85: DepV, 86: DepW.");
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