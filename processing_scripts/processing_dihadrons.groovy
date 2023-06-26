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

	// Set the output file name based on the provided 3rd argument or use the default name
	String output_file = args.length < 4 ? "hadron_dummy_out.txt" : args[3];
	if (args.length < 4) 
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
	int helicity;
	double e_p, e_theta, e_phi, p1_phi, p1_p, p1_theta, p2_phi, p2_p, p2_theta; 
	double vz_e, vz_p1, vz_p2;
	double Q2, W, y, Mx, Mx1, Mx2; 
	double x, z, xF, pT, eta, eta_gN, zeta;
	double z1, z2, xF1, xF2, Mh, pT1, pT2, pTpT, eta1, eta2, Delta_eta, eta1_gN, eta2_gN;
	double phi1, phi2, Delta_phi, phih, phiR, theta;
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
		    // boolean process_event = filter.isValid(research_Event) && 
		    // 	(runnum == 11 || runnum >= 11571 || qa.OkForAsymmetry(runnum, evnum));
		    boolean process_event = filter.isValid(research_Event);

		    if (process_event) {

		        // get # of particles 
		        int num_p1 = research_Event.countByPid(p1_Str.toInteger());
		        int num_p2 = research_Event.countByPid(p2_Str.toInteger()); 

		        // cycle over all hadrons
		        for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { 
		        	for (int current_p2 = 0; current_p2 < num_p2; current_p2++) { 
		        		if (current_p1 == current_p2 && p1_int == p2_int) {continue; }

			            Dihadrons variables = new Dihadrons(event, research_Event, 
							p1_int, current_p1, p2_int, current_p2);
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

			                // vertices
			                vz_e = variables.vz_e();
			                vz_p1 = variables.vz_p1();
			                vz_p2 = variables.vz_p2();

			                // DIS variables
			                Q2 = variables.Q2(); // exchanged virtual photon energy
			                W = variables.W(); // hadronic mass
			                x = variables.x(); // Bjorken-x
			                y = variables.y(); // E_scat/E_beam
			                Mx = variables.Mx(); // missing mass
			                Mx1 = variables.Mx1(); // missing mass calculated with p1
			                Mx2 = variables.Mx2(); // missing mass calculated with p2

			                // SIDIS variables
			                z = variables.z(); // fractional hadron energy wrt virtual photon
			                xF = variables.xF(); // Feynman-x
			                pT = variables.pT(); // transverse momentum of hadron
			                eta = variables.eta(); // rapidity
			                eta_gN = variables.eta_gN();
			                zeta = variables.zeta(); // longitudinal momentum of hadron (fracture functions)

			                // SIDIS dihadron variables
							z1 = variables.z1();
							z2 = variables.z2();
							xF1 = variables.xF1();
							xF2 = variables.xF2();
							zeta1 = variables.zeta1();
							zeta2 = variables.zeta2(); 
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
			                	.append(Q2).append(" ")
			                	.append(W).append(" ")
			                	.append(Mx).append(" ")
			                	.append(Mx1).append(" ")
			                	.append(Mx2).append(" ")
			                	.append(x).append(" ")
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
			                	.append(zeta).append(" ")
			                	.append(zeta1).append(" ")
			                	.append(zeta2).append(" ")
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
		reader.close();
		}

		// Write any remaining lines in the batchLines StringBuilder to the file
		if (batchLines.length() > 0) {
		    file.append(batchLines.toString());
		    batchLines.setLength(0);
		}

		println("\n1:runnum, 2:evnum, 3:helicity, 4:e_p, 5:e_theta, 6:e_phi, 7:vz_e,"+
		"8:p1_p, 9:p1_theta, 10:p1_phi, 11:vz_p1, 12:p2_p, 13:p2_theta, 14:p2_phi, 15:vz_p2"+
		"16:Q2, 17:W, 18:Mx, 19: Mx1, 20: Mx2, 21:x, 22:y, 23:z,"+
		"24:z1, 25:z2, 26:Mh, 27:xF, 28:xF1, 29:xF2, 30:pT, 31:pT1, 32:pT2, 33:pTpT"+
		"34:zeta, 35:zeta1, 36:zeta2, 37:eta, 38:eta1, 39:eta2, 40:Delta_eta, 41:eta1_gN, 42:eta2_gN"+
		"43:phi1, 44:phi2, 45:Delta_phi, 46:phih, 47:phiR, 48:theta"+
		"49:DepA, 50:DepB, 51:DepC, 52:DepV, 53:DepW\n");

		println("Set p1 PID = $p1_Str");
		println("Set p2 PID = $p2_Str");
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