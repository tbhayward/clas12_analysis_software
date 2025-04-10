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
	String output_file = args.length < 2 ? "dvcs_dummy_out.txt" : args[1];
	if (args.length < 2) 
	    println("WARNING: Specify an output file name. Set to \"dvcs_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

	// Set the number of files to process based on the provided 3rd argument or list size
	// If the argument is "0", default to the full list size
	int n_files = args.length < 3 || Integer.parseInt(args[2]) == 0 || Integer.parseInt(args[2]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[2]);
	if (args.length < 3 || Integer.parseInt(args[2]) == 0 || Integer.parseInt(args[2]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified, set to 0, or too large
	    println("WARNING: Number of files not specified, set to 0 or number too large.")
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

	println("Script will analyze generated banks and any reconstructed banks and save both.")

	// ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //

	// declare physics event variables
	int reconstructed; // variable to declare whether the generated event has reconstructed particles
	int fiducial_status, num_pos, num_neg, num_neutrals, helicity, detector1, detector2;
	double e_p, e_theta, e_phi, p1_phi, p1_p, p1_theta, p2_phi, p2_p, p2_theta; 
	double vz_e, vz_p1, vz_p2;
	double open_angle_ep, open_angle_ep1, open_angle_ep2, open_angle_p1p2;
	double Q2, W, y, Mx2, Mx2_1, Mx2_2; 
	double x, t, t1, t2, tmin, z, xF, pT, eta, eta_gN, xi;
	double z1, z2, xF1, xF2, Mh, pT1, pT2, pTpT, eta1, eta2, Delta_eta, eta1_gN, eta2_gN;
	double phi1, phi2, Delta_phi, phih, phiR, theta;
	double Depolarization_A, Depolarization_B, Depolarization_C;
	double Depolarization_V, Depolarization_W;
	double Emiss2, theta_gamma_gamma, pTmiss;

	double gen_e_p, gen_e_theta, gen_e_phi, gen_p1_phi, gen_p1_p, gen_p1_theta, gen_p2_phi, gen_p2_p, gen_p2_theta; 
	double gen_vz_e, gen_vz_p1, gen_vz_p2;
	double gen_open_angle_ep, gen_open_angle_ep1, gen_open_angle_ep2, gen_open_angle_p1p2;
	double gen_Q2, gen_W, gen_y, gen_Mx2, gen_Mx2_1, gen_Mx2_2; 
	double gen_x, gen_t, gen_t1, gen_t2, gen_tmin, gen_z, gen_xF, gen_pT, gen_eta, gen_eta_gN, gen_xi;
	double gen_z1, gen_z2, gen_xF1, gen_xF2, gen_Mh, gen_pT1, gen_pT2, gen_pTpT, gen_eta1, gen_eta2, gen_Delta_eta, gen_eta1_gN, gen_eta2_gN;
	double gen_phi1, gen_phi2, gen_Delta_phi, gen_phih, gen_phiR, gen_theta;
	double gen_Depolarization_A, gen_Depolarization_B, gen_Depolarization_C;
	double gen_Depolarization_V, gen_Depolarization_W;
	double gen_Emiss2, gen_theta_gamma_gamma, gen_pTmiss;

	// load my kinematic fitter/PID
	GenericKinematicFitter rec_fitter = new dvcs_fitter(10.6041); 
	GenericKinematicFitter gen_fitter = new monte_carlo_fitter(10.6041);
	// GenericKinematicFitter fitter = new event_builder_fitter(10.6041); 
	
	// set filter for final states
	EventFilter filter = new EventFilter("11:2212:22:Xn"); 

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

			reconstructed = -999; // if you intend to use –999 as a flag (though normally this would be a boolean)
			fiducial_status = -999;
			num_pos = -999;
			num_neg = -999;
			num_neutrals = -999;
			helicity = -999;
			detector1 = -999;
			detector2 = -999;

			e_p = -999;
			e_theta = -999;
			e_phi = -999;
			p1_phi = -999;
			p1_p = -999;
			p1_theta = -999;
			p2_phi = -999;
			p2_p = -999;
			p2_theta = -999;

			vz_e = -999;
			vz_p1 = -999;
			vz_p2 = -999;

			open_angle_ep = -999;
			open_angle_ep1 = -999;
			open_angle_ep2 = -999;
			open_angle_p1p2 = -999;

			Q2 = -999;
			W = -999;
			y = -999;
			Mx2 = -999;
			Mx2_1 = -999;
			Mx2_2 = -999;

			x = -999;
			t = -999;
			t1 = -999;
			t2 = -999;
			tmin = -999;
			z = -999;
			xF = -999;
			pT = -999;
			eta = -999;
			eta_gN = -999;
			xi = -999;

			z1 = -999;
			z2 = -999;
			xF1 = -999;
			xF2 = -999;
			Mh = -999;
			pT1 = -999;
			pT2 = -999;
			pTpT = -999;
			eta1 = -999;
			eta2 = -999;
			Delta_eta = -999;
			eta1_gN = -999;
			eta2_gN = -999;

			phi1 = -999;
			phi2 = -999;
			Delta_phi = -999;
			phih = -999;
			phiR = -999;
			theta = -999;

			Depolarization_A = -999;
			Depolarization_B = -999;
			Depolarization_C = -999;

			Depolarization_V = -999;
			Depolarization_W = -999;

			Emiss2 = -999;
			theta_gamma_gamma = -999;
			pTmiss = -999;

			// And similarly for the gen_ variables:
			gen_e_p = -999;
			gen_e_theta = -999;
			gen_e_phi = -999;
			gen_p1_phi = -999;
			gen_p1_p = -999;
			gen_p1_theta = -999;
			gen_p2_phi = -999;
			gen_p2_p = -999;
			gen_p2_theta = -999;

			gen_vz_e = -999;
			gen_vz_p1 = -999;
			gen_vz_p2 = -999;

			gen_open_angle_ep = -999;
			gen_open_angle_ep1 = -999;
			gen_open_angle_ep2 = -999;
			gen_open_angle_p1p2 = -999;

			gen_Q2 = -999;
			gen_W = -999;
			gen_y = -999;
			gen_Mx2 = -999;
			gen_Mx2_1 = -999;
			gen_Mx2_2 = -999;

			gen_x = -999;
			gen_t = -999;
			gen_t1 = -999;
			gen_t2 = -999;
			gen_tmin = -999;
			gen_z = -999;
			gen_xF = -999;
			gen_pT = -999;
			gen_eta = -999;
			gen_eta_gN = -999;
			gen_xi = -999;

			gen_z1 = -999;
			gen_z2 = -999;
			gen_xF1 = -999;
			gen_xF2 = -999;
			gen_Mh = -999;
			gen_pT1 = -999;
			gen_pT2 = -999;
			gen_pTpT = -999;
			gen_eta1 = -999;
			gen_eta2 = -999;
			gen_Delta_eta = -999;
			gen_eta1_gN = -999;
			gen_eta2_gN = -999;

			gen_phi1 = -999;
			gen_phi2 = -999;
			gen_Delta_phi = -999;
			gen_phih = -999;
			gen_phiR = -999;
			gen_theta = -999;

			gen_Depolarization_A = -999;
			gen_Depolarization_B = -999;
			gen_Depolarization_C = -999;

			gen_Depolarization_V = -999;
			gen_Depolarization_W = -999;

			gen_Emiss2 = -999;
			gen_theta_gamma_gamma = -999;
			gen_pTmiss = -999;
		    ++num_events;
		    if (num_events % 500000 == 0) { // not necessary, just updates output
		        print("processed: " + num_events + " events. ");
		    }

		    // get run and event numbers
		    event = reader.getNextEvent();
		    // collect info for QA
		    int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0);
		    int evnum = event.getBank("RUN::config").getInt('event', 0);

		    PhysicsEvent rec_Event = rec_fitter.getPhysicsEvent(event);
		    PhysicsEvent gen_Event = gen_fitter.getPhysicsEvent(event);

	    	if (runnum != 11) {
			    throw new IllegalArgumentException("error: runnum != 11; this script is intended for use on MC (runnum ==11).")
			}
			boolean process_event = filter.isValid(gen_Event)
			reconstructed = filter.isValid(rec_Event) ? 1 : 0;
			println(reconstructed);

	    	if (process_event) {

		        // get # of particles 
		        int num_p1 = gen_Event.countByPid(2212);
		        int num_p2 = gen_Event.countByPid(22); 

        		// supply runnum and boolean for radiative simulation or not
				BeamEnergy Eb = new BeamEnergy(gen_Event, runnum, false);
				// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
				double energy = (runnum == 11) ? beam_energy : Eb.Eb();
	            ThreeParticles variables = new ThreeParticles(event, gen_Event, 
					2212, 0, 22, 0, energy);
	            // this is my class for defining all relevant kinematic variables
	            if (variables.channel_test(variables)) {
	                // lab kinematics
	                gen_e_p = variables.e_p(); // lab frame momentum
	                gen_e_theta = variables.e_theta(); // lab polar angle
	                gen_e_phi = variables.e_phi(); // lab azimuthal angle
	                gen_p1_phi = variables.p1_phi(); // lab azimuthal angle
	                gen_p1_p = variables.p1_p(); // lab momentum
	                gen_p1_theta = variables.p1_theta(); // lab polar angle
	                gen_p2_phi = variables.p2_phi(); // lab azimuthal angle
	                gen_p2_p = variables.p2_p(); // lab momentum
	                gen_p2_theta = variables.p2_theta(); // lab polar angle
	                gen_open_angle_ep = variables.open_angle_ep;
	                gen_open_angle_ep1 = variables.open_angle_ep1;
	                gen_open_angle_ep2 = variables.open_angle_ep2;
	                gen_open_angle_p1p2 = variables.open_angle_p1p2;

	                // vertices
	                gen_vz_e = variables.vz_e();
	                gen_vz_p1 = variables.vz_p1();
	                gen_vz_p2 = variables.vz_p2();

	                // DIS variables
	                gen_Q2 = variables.Q2(); // exchanged virtual photon energy
	                gen_W = variables.W(); // hadronic mass
	                gen_x = variables.x(); // Bjorken-x
	                gen_t = variables.t();
	                gen_t1 = variables.t1();
	                gen_t2 = variables.t2();
	                gen_tmin = variables.tmin();
	                gen_y = variables.y(); // E_scat/E_beam
	                gen_Mx2 = variables.Mx2(); // missing mass
	                gen_Mx2_1 = variables.Mx2_1(); // missing mass calculated with p1
	                gen_Mx2_2 = variables.Mx2_2(); // missing mass squared

	                // SIDIS variables
	                gen_z = variables.z(); // fractional hadron energy wrt virtual photon
	                gen_xF = variables.xF(); // Feynman-x
	                gen_pT = variables.pT(); // transverse momentum of hadron
	                

	                // SIDIS dihadron variables
					gen_z1 = variables.z1();
					gen_z2 = variables.z2();
					gen_xF1 = variables.xF1();
					gen_xF2 = variables.xF2();
					gen_Mh = variables.Mh();
					gen_pT1 = variables.pT1();
					gen_pT2 = variables.pT2();

	                // angles
					gen_phi1 = variables.phi1(); // trento phi of the p1
					gen_phi2 = variables.phi2(); // trento phi of the p2
					gen_Delta_phi = variables.Delta_phi();

	                // depolarization factors
	                gen_Depolarization_A = variables.Depolarization_A();
	                gen_Depolarization_B = variables.Depolarization_B();
	                gen_Depolarization_C = variables.Depolarization_C();
	                gen_Depolarization_V = variables.Depolarization_V();
			    	gen_Depolarization_W = variables.Depolarization_W();

			    	// exclusivity variables
			    	gen_Emiss2 = variables.Emiss2();
			    	gen_theta_gamma_gamma = variables.theta_gamma_gamma();
			    	gen_pTmiss = variables.pTmiss();
	            }  
		    }

		    if (process_event && reconstructed) {

		        // get # of particles 
		        int num_p1 = rec_Event.countByPid(2212);
		        int num_p2 = rec_Event.countByPid(22); 

        		// supply runnum and boolean for radiative simulation or not
				BeamEnergy Eb = new BeamEnergy(rec_Event, runnum, false);
				// Use the input beam energy if runnum == 11, otherwise use Eb.Eb()
				double energy = (runnum == 11) ? beam_energy : Eb.Eb();
	            ThreeParticles variables = new ThreeParticles(event, rec_Event, 
					2212, 0, 22, 0, energy);
	            // this is my class for defining all relevant kinematic variables
	            if (variables.channel_test(variables)) {
	            	fiducial_status = variables.get_fiducial_status();
	                helicity = variables.get_helicity(); // helicity of event
	                detector1 = variables.get_detector1();
	                detector2 = variables.get_detector2();
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
	                open_angle_ep = variables.open_angle_ep;
	                open_angle_ep1 = variables.open_angle_ep1;
	                open_angle_ep2 = variables.open_angle_ep2;
	                open_angle_p1p2 = variables.open_angle_p1p2;

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
	                Mx2 = variables.Mx2(); // missing mass
	                Mx2_1 = variables.Mx2_1(); // missing mass calculated with p1
	                Mx2_2 = variables.Mx2_2(); // missing mass squared

	                // SIDIS variables
	                z = variables.z(); // fractional hadron energy wrt virtual photon
	                xF = variables.xF(); // Feynman-x
	                pT = variables.pT(); // transverse momentum of hadron

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

	                // angles
					phi1 = variables.phi1(); // trento phi of the p1
					phi2 = variables.phi2(); // trento phi of the p2
					Delta_phi = variables.Delta_phi();

	                // depolarization factors
	                Depolarization_A = variables.Depolarization_A();
	                Depolarization_B = variables.Depolarization_B();
	                Depolarization_C = variables.Depolarization_C();
	                Depolarization_V = variables.Depolarization_V();
			    	Depolarization_W = variables.Depolarization_W();

			    	// exclusivity variables
			    	Emiss2 = variables.Emiss2();
			    	theta_gamma_gamma = variables.theta_gamma_gamma();
			    	pTmiss = variables.pTmiss();
	            }  
		    }
		    // Use a StringBuilder to append all data in a single call
	                StringBuilder line = new StringBuilder();
	                // first the generated variables
	                line.append(gen_e_p).append(" ")
	                	.append(gen_e_theta).append(" ")
	                	.append(gen_e_phi).append(" ")
	                	.append(gen_vz_e).append(" ")
	                	.append(gen_p1_p).append(" ")
	                	.append(gen_p1_theta).append(" ")
	                	.append(gen_p1_phi).append(" ")
	                	.append(gen_vz_p1).append(" ")
	                	.append(gen_p2_p).append(" ")
	                	.append(gen_p2_theta).append(" ")
	                	.append(gen_p2_phi).append(" ")
	                	.append(gen_vz_p2).append(" ")
	                	.append(gen_open_angle_ep).append(" ")
	                	.append(gen_open_angle_ep1).append(" ")
	                	.append(gen_open_angle_ep2).append(" ")
	                	.append(gen_open_angle_p1p2).append(" ")
	                	.append(gen_Q2).append(" ")
	                	.append(gen_W).append(" ")
	                	.append(gen_Mx2).append(" ")
	                	.append(gen_Mx2_1).append(" ")
	                	.append(gen_Mx2_2).append(" ")
	                	.append(gen_x).append(" ")
	                	.append(gen_t).append(" ")
	                	.append(gen_t1).append(" ")
	                	.append(gen_t2).append(" ")
	                	.append(gen_tmin).append(" ")
	                	.append(gen_y).append(" ")
	                	.append(gen_z).append(" ")
	                	.append(gen_z1).append(" ")
	                	.append(gen_z2).append(" ")
	                	.append(gen_Mh).append(" ")
	                	.append(gen_xF).append(" ")
	                	.append(gen_xF1).append(" ")
	                	.append(gen_xF2).append(" ")
	                	.append(gen_pT).append(" ")
	                	.append(gen_pT1).append(" ")
	                	.append(gen_pT2).append(" ")
	                	.append(gen_phi1).append(" ")
	                	.append(gen_phi2).append(" ")
	                	.append(gen_Delta_phi).append(" ")
	                	.append(gen_Depolarization_A).append(" ")
	                    .append(gen_Depolarization_B).append(" ")
	                    .append(gen_Depolarization_C).append(" ")
	                    .append(gen_Depolarization_V).append(" ")
	                    .append(gen_Depolarization_W).append(" ")
	                    .append(gen_Emiss2).append(" ")
	                    .append(gen_theta_gamma_gamma).append(" ")
	                    .append(gen_pTmiss).append(" ")
	                    // now reconstructed variables
	                    .append(fiducial_status).append(" ")
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
	                	.append(phi1).append(" ")
	                	.append(phi2).append(" ")
	                	.append(Delta_phi).append(" ")
	                	.append(Depolarization_A).append(" ")
	                    .append(Depolarization_B).append(" ")
	                    .append(Depolarization_C).append(" ")
	                    .append(Depolarization_V).append(" ")
	                    .append(Depolarization_W).append(" ")
	                    .append(Emiss2).append(" ")
	                    .append(theta_gamma_gamma).append(" ")
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
		reader.close();

		// Write any remaining lines in the batchLines StringBuilder to the file
		if (batchLines.length() > 0) {
		    file.append(batchLines.toString());
		    batchLines.setLength(0);
		}

		println("1: gen_e_p, 2: gen_e_theta, 3: gen_e_phi, 4: gen_vz_e, " +
		    "5: gen_p1_p, 6: gen_p1_theta, 7: gen_p1_phi, 8: gen_vz_p1, " +
		    "9: gen_p2_p, 10: gen_p2_theta, 11: gen_p2_phi, 12: gen_vz_p2, " +
		    "13: gen_open_angle_ep, 14: gen_open_angle_ep1, 15: gen_open_angle_ep2, 16: gen_open_angle_p1p2, " +
		    "17: gen_Q2, 18: gen_W, 19: gen_Mx2, 20: gen_Mx2_1, 21: gen_Mx2_2, " +
		    "22: gen_x, 23: gen_t, 24: gen_t1, 25: gen_t2, 26: gen_tmin, " +
		    "27: gen_y, 28: gen_z, " +
		    "29: gen_z1, 30: gen_z2, 31: gen_Mh, " +
		    "32: gen_xF, 33: gen_xF1, 34: gen_xF2, " +
		    "35: gen_pT, 36: gen_pT1, 37: gen_pT2, " +
		    "38: gen_phi1, 39: gen_phi2, 40: gen_Delta_phi, 41: gen_phih, " +
		    "42: gen_phiR, 43: gen_theta, " +
		    "44: gen_DepA, 45: gen_DepB, 46: gen_DepC, 47: gen_DepV, 48: gen_DepW, " +
		    "49: gen_Emiss2, 50: gen_theta_gamma_gamma, " +
		    "51: gen_pTmiss");

		println("52: fiducial_status, 53: num_pos, 54: num_neg, 55: num_neutrals, " +
		    "56: runnum, 57: evnum, 58: helicity, 59: detector1, 60: detector2, " +
		    "61: e_p, 62: e_theta, 63: e_phi, 64: vz_e, " +
		    "65: p1_p, 66: p1_theta, 67: p1_phi, 68: vz_p1, " +
		    "69: p2_p, 70: p2_theta, 71: p2_phi, 72: vz_p2, " +
		    "73: open_angle_ep, 74: open_angle_ep1, 75: open_angle_ep2, 76: open_angle_p1p2, " +
		    "77: Q2, 78: W, 79: Mx2, 80: Mx2_1, 81: Mx2_2, " +
		    "82: x, 83: t, 84: t1, 85: t2, 86: tmin, 87: y, 88: z, " +
		    "89: z1, 90: z2, 91: Mh, " +
		    "92: xF, 93: xF1, 94: xF2, " +
		    "95: pT, 96: pT1, 97: pT2, " +
		    "98: phi1, 99: phi2, 100: Delta_phi, 101: phih, " +
		    "102: phiR, 103: theta, " +
		    "104: DepA, 105: DepB, 106: DepC, 107: DepV, 108: DepW, " +
		    "109: Emiss2, 110: theta_gamma_gamma, " +
		    "111: pTmiss, 112: reconstructed.");

		println("Analyzing dvcs.");
		println("output text file is: $file");
	}

	writer.close();

	// End time
	long endTime = System.currentTimeMillis()
	// Calculate the elapsed time
	long elapsedTime = endTime - startTime
	// Print the elapsed time in milliseconds
	println("Elapsed time: ${elapsedTime} ms");
	println("The amount of variables stored in this tree is unecessarily large (many are not applicable to exclusive analyses).");
	println("If you're not explicitly using the .txt file it might be worth deleting it and just keeping the .root tree.");

}