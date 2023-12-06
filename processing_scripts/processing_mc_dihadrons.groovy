/*
 * author Timothy B. Hayward
 * 
 * SIDIS dihadron mc (generated and reconstructed variables saved)
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

	// Set the output file name based on the provided 4th argument or use the default name
	String output_file = args.length < 4 ? "hadron_dummy_out.txt" : args[3];
	if (args.length < 4) 
	    println("WARNING: Specify an output file name. Set to \"dihadron_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();
	BufferedWriter writer = new BufferedWriter(new FileWriter(file));

	// Set the number of files to process based on the provided 5th argument
	// use the size of the hipo_list if no argument provided
	int n_files = args.length < 5 || Integer.parseInt(args[4]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[4]);
	if (args.length < 5 || Integer.parseInt(args[4]) > hipo_list.size()) {
	    // Print warnings and information if the number of files is not specified or too large
	    println("WARNING: Number of files not specified or number too large.")
	    println("Setting # of files to be equal to number of files in the directory.");
	    println("There are $hipo_list.size files.");
	}

	int hadron_pair_counts = 0;
	GenericKinematicFitter research_fitter = new analysis_fitter(10.6041);
	// GenericKinematicFitter research_fitter=new proton_energy_loss_corrections_fitter(10.6041);
	GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6041);
	EventFilter filter = new EventFilter("11:"+p1_Str+":"+p2_Str+":X+:X-:Xn");

	// create a StringBuilder for accumulating lines
	StringBuilder batchLines = new StringBuilder();

	int num_events = 0;
	int current_file = 0;
	int max_lines = 1000;
	int lineCount = 0;
	while (current_file < n_files) {
		println(); println(); println("Opening file "+Integer.toString(current_file+1)
			+" out of "+n_files); println(); println();
		// limit to a certain number of files defined by n_files

		HipoDataSource reader = new HipoDataSource();

		reader.open(hipo_list[current_file]); // open next hipo file
		current_file++;
		HipoDataEvent event = reader.getNextEvent(); 

		while(reader.hasEvent()==true){
			num_events++; 
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
				HipoDataBank lundBank = (HipoDataBank) event.getBank("MC::Lund");
				HipoDataBank mcBank = (HipoDataBank) event.getBank("MC::Particle");

				int num_p1 = research_Event.countByPid(p1_Str.toInteger()); 
				int num_p2 = research_Event.countByPid(p2_Str.toInteger());

				for (int current_p1 = 0; current_p1 < num_p1; current_p1++) {
					for (int current_p2 = 0; current_p2 < num_p2; current_p2++) {

						Particle exp_e = research_Event.getParticleByPid(11,0);
						Particle exp_p1 = research_Event.getParticleByPid(p1_Str.toInteger(),current_p1);
						Particle exp_p2 = research_Event.getParticleByPid(p2_Str.toInteger(),current_p2);

						BeamEnergy Eb = new BeamEnergy(runnum, true);
						Dihadrons variables = new Dihadrons(event, research_Event, 
							p1_Str.toInteger(), current_p1, p2_Str.toInteger(), current_p2, 
							Eb.Eb());
						Dihadrons mc_variables = new Dihadrons(event, mc_Event, 
							p1_Str.toInteger(), current_p1, p2_Str.toInteger(), current_p2,
							Eb.Eb());

						if (variables.channel_test(variables)) {

							// lab kinematics data
							double e_p = variables.e_p();
							double e_theta = variables.e_theta();
							double e_phi = variables.e_phi();
							double p1_p = variables.p1_p();
							double p1_theta = variables.p1_theta();
							double p1_phi = variables.p1_phi();
							double p2_p = variables.p2_p();
							double p2_theta = variables.p2_theta();
							double p2_phi = variables.p2_phi();

							// lab kinematics MC
							double mc_e_p = mc_variables.e_p();
							double mc_e_theta = mc_variables.e_theta();
							double mc_e_phi = mc_variables.e_phi();
							double mc_p1_p = mc_variables.p1_p();
							double mc_p1_theta = mc_variables.p1_theta();
							double mc_p1_phi = mc_variables.p1_phi();
							double mc_p2_p = mc_variables.p2_p();
							double mc_p2_theta = mc_variables.p2_theta();
							double mc_p2_phi = mc_variables.p2_phi();

							// DIS variables data
							double Q2 = variables.Q2();
							double W = variables.W();
							double y = variables.y();
							double Mx = variables.Mx();
							double Mx1 = variables.Mx1();
							double Mx2 = variables.Mx2();

							// DIS variables MC
							double mc_Q2 = mc_variables.Q2();
							double mc_W = mc_variables.W(); 
							double mc_y = mc_variables.y();
							double mc_Mx = mc_variables.Mx();
							double mc_Mx1 = mc_variables.Mx1();
							double mc_Mx2 = mc_variables.Mx2();

							// SIDIS variables data
							double x = variables.x();
							double z = variables.z();
							double xF = variables.xF();
							double pT = variables.pT();
							double eta = variables.eta();

							// SIDIS variables MC
							double mc_x = mc_variables.x();
							double mc_z = mc_variables.z();
							double mc_xF = mc_variables.xF();
							double mc_pT = mc_variables.pT();
							double mc_eta = mc_variables.eta();

							// SIDIS dihadron variables data
							double z1 = variables.z1();
							double z2 = variables.z2();
							double xF1 = variables.xF1();
							double xF2 = variables.xF2();
							double zeta = variables.zeta();
							double zeta1 = variables.zeta1();
							double zeta2 = variables.zeta2();
							double Mh = variables.Mh();
							double pT1 = variables.pT1();
							double pT2 = variables.pT2();
							double pTpT = variables.pTpT();
							double eta1 = variables.eta1();
							double eta2 = variables.eta2();
							double Delta_eta = variables.Delta_eta();

							// SIDIS dihadron variables MC
							double mc_z1 = mc_variables.z1();
							double mc_z2 = mc_variables.z2();
							double mc_xF1 = mc_variables.xF1();
							double mc_xF2 = mc_variables.xF2();
							double mc_zeta = mc_variables.zeta(); 
							double mc_zeta1 = mc_variables.zeta1(); 
							double mc_zeta2 = mc_variables.zeta2(); 
							double mc_Mh = mc_variables.Mh();
							double mc_pT1 = mc_variables.pT1();
							double mc_pT2 = mc_variables.pT2();
							double mc_pTpT = mc_variables.pTpT();
							double mc_eta1 = mc_variables.eta1();
							double mc_eta2 = mc_variables.eta2();
							double mc_Delta_eta = mc_variables.Delta_eta();

							// angles data
							double phi1 = variables.phi1();
							double phi2 = variables.phi2();
							double Delta_phi = variables.Delta_phi();
							double phih = variables.phih();
							double phiR = variables.phiR();
							double theta = variables.theta();

							// angles MC
							double mc_phi1 = mc_variables.phi1();
							double mc_phi2 = mc_variables.phi2();
							double mc_Delta_phi = mc_variables.Delta_phi();
							double mc_phih = mc_variables.phih();
							double mc_phiR = mc_variables.phiR();
							double mc_theta = mc_variables.theta();

							// depolarization factors data
							double Depolarization_A = variables.Depolarization_A();
							double Depolarization_B = variables.Depolarization_B();
							double Depolarization_C = variables.Depolarization_C();
							double Depolarization_V = variables.Depolarization_V();
							double Depolarization_W = variables.Depolarization_W();

							// depolarization factors MC
							double mc_Depolarization_A = mc_variables.Depolarization_A();
							double mc_Depolarization_B = mc_variables.Depolarization_B();
							double mc_Depolarization_C = mc_variables.Depolarization_C();
							double mc_Depolarization_V = mc_variables.Depolarization_V();
							double mc_Depolarization_W = mc_variables.Depolarization_W();

							boolean matching_e = false;
							boolean matching_p1 = false;
							boolean matching_p2 = false;

							int matching_p1_pid = 0;
							int mc_p1_parent_index = 0;
							for (int current_part = 0; current_part < lundBank.rows(); current_part++) {
								int pid = lundBank.getInt("pid", current_part);
								if (matching_p1) { continue; }
								double mc_px = lundBank.getFloat("px", current_part);
								double mc_py = lundBank.getFloat("py", current_part);
								double mc_pz = lundBank.getFloat("pz", current_part);

								double mc_phi = phi_calculation(mc_px, mc_py);
								double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

								double exp_phi = phi_calculation(exp_p1.px(), exp_p1.py());
								double exp_theta = theta_calculation(exp_p1.px(), exp_p1.py(), 
									exp_p1.pz());

								matching_p1 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
									Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
								if (matching_p1) {
									matching_p1_pid = pid;
									mc_p1_parent_index = lundBank.getInt("parent", current_part)-1;
								}
							}
							int mc_p1_parent = lundBank.getInt("pid", mc_p1_parent_index);

							int matching_p2_pid = 0;
							int mc_p2_parent_index = 0;
							for (int current_part = 0; current_part < lundBank.rows(); current_part++) {
								int pid = lundBank.getInt("pid", current_part);
								if (matching_p2) { continue; }
								double mc_px = lundBank.getFloat("px", current_part);
								double mc_py = lundBank.getFloat("py", current_part);
								double mc_pz = lundBank.getFloat("pz", current_part);

								double mc_phi = phi_calculation(mc_px, mc_py);
								double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

								double exp_phi = phi_calculation(exp_p2.px(), exp_p2.py());
								double exp_theta = theta_calculation(exp_p2.px(), exp_p2.py(), 
									exp_p2.pz());

								matching_p2 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
									Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
								if (matching_p2) {
									matching_p2_pid = pid;
									mc_p2_parent_index = lundBank.getInt("parent", current_part)-1;
								} 
							}
							int mc_p2_parent = lundBank.getInt("pid", mc_p2_parent_index);

							matching_e = false;
							matching_p1 = false;
							matching_p2 = false;

							int matching_e_pid = 0;
							int mc_e_parent_index = 0;
							for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
								int pid = mcBank.getInt("pid", current_part);
								if (matching_e) { continue; }
								double mc_px = mcBank.getFloat("px", current_part);
								double mc_py = mcBank.getFloat("py", current_part);
								double mc_pz = mcBank.getFloat("pz", current_part);

								double mc_phi = phi_calculation(mc_px, mc_py);
								double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

								double exp_phi = phi_calculation(exp_e.px(), exp_e.py());
								double exp_theta = theta_calculation(exp_e.px(), exp_e.py(), exp_e.pz());

								matching_e = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
									Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
								if (matching_e) {
									matching_e_pid = pid;
								}
							}

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

							matching_p2_pid = 0;
							mc_p2_parent_index = 0;
							for (int current_part = 0; current_part < mcBank.rows(); current_part++) {
								int pid = mcBank.getInt("pid", current_part);
								if (matching_p2) { continue; }
								double mc_px = mcBank.getFloat("px", current_part);
								double mc_py = mcBank.getFloat("py", current_part);
								double mc_pz = mcBank.getFloat("pz", current_part);

								double mc_phi = phi_calculation(mc_px, mc_py);
								double mc_theta_lab = theta_calculation(mc_px, mc_py, mc_pz);

								double exp_phi = phi_calculation(exp_p2.px(), exp_p2.py());
								double exp_theta = theta_calculation(exp_p2.px(), exp_p2.py(), 
									exp_p2.pz());

								matching_p2 = Math.abs(exp_phi - mc_phi) < scale*3.0 && 
									Math.abs(exp_theta - mc_theta_lab) < scale*1.0;
								if (matching_p2) {
									matching_p2_pid = pid;
								} 
							}

							// Use a StringBuilder to append all data in a single call
			                StringBuilder line = new StringBuilder();
			                line.append(e_p).append(" ")			// 1
				                .append(mc_e_p).append(" ")			// 2
				                .append(e_theta).append(" ")		// 3
								.append(mc_e_theta).append(" ")		// 4
								.append(e_phi).append(" ")			// 5
								.append(mc_e_phi).append(" ")		// 6
								.append(p1_p).append(" ")			// 7
				                .append(mc_p1_p).append(" ")		// 8
				                .append(p1_theta).append(" ")		// 9 
								.append(mc_p1_theta).append(" ")	// 10
								.append(p1_phi).append(" ")			// 11
								.append(mc_p1_phi).append(" ")		// 12
								.append(p2_p).append(" ")			// 13
				                .append(mc_p2_p).append(" ")		// 14
				                .append(p2_theta).append(" ")		// 15
								.append(mc_p2_theta).append(" ")	// 16
								.append(p2_phi).append(" ")			// 17
								.append(mc_p2_phi).append(" ")		// 18
								.append(Q2).append(" ")				// 19
								.append(mc_Q2).append(" ")			// 20
								.append(W).append(" ")				// 21
								.append(mc_W).append(" ")			// 22
								.append(Mx).append(" ")				// 23
								.append(mc_Mx).append(" ")			// 24
								.append(Mx1).append(" ")			// 25
								.append(mc_Mx1).append(" ")			// 26
								.append(Mx2).append(" ")			// 27
								.append(mc_Mx2).append(" ")			// 28
								.append(x).append(" ")				// 29
								.append(mc_x).append(" ")			// 30
								.append(y).append(" ")				// 31
								.append(mc_y).append(" ")			// 32
								.append(z).append(" ")				// 33
								.append(mc_z).append(" ")			// 34
								.append(z1).append(" ")				// 35
								.append(mc_z1).append(" ")			// 36
								.append(z2).append(" ")				// 37
								.append(mc_z2).append(" ")			// 38
								.append(Mh).append(" ")				// 39
								.append(mc_Mh).append(" ")			// 40
								.append(xF).append(" ")				// 41
								.append(mc_xF).append(" ")			// 42
								.append(xF1).append(" ")			// 43
								.append(mc_xF1).append(" ")			// 44
								.append(xF2).append(" ")			// 45
								.append(mc_xF2).append(" ")			// 46
								.append(pT).append(" ")				// 47
								.append(mc_pT).append(" ")			// 48
								.append(pT1).append(" ")			// 49
								.append(mc_pT1).append(" ")			// 50
								.append(pT2).append(" ")			// 51
								.append(mc_pT2).append(" ")			// 52
								.append(pTpT).append(" ")			// 53
								.append(mc_pTpT).append(" ")		// 54
								.append(zeta).append(" ")			// 55
								.append(mc_zeta).append(" ")		// 56
								.append(zeta1).append(" ")			// 57
								.append(mc_zeta1).append(" ")		// 58
								.append(zeta2).append(" ")			// 59
								.append(mc_zeta2).append(" ")		// 60
								.append(eta).append(" ")			// 61
								.append(mc_eta).append(" ")			// 62
								.append(eta1).append(" ")			// 63
								.append(mc_eta1).append(" ")		// 64
								.append(eta2).append(" ")			// 65
								.append(mc_eta2).append(" ")		// 66
								.append(Delta_eta).append(" ")		// 67
								.append(mc_Delta_eta).append(" ")	// 68
								.append(phi1).append(" ")			// 69
								.append(mc_phi1).append(" ")		// 70
								.append(phi2).append(" ")			// 71
								.append(mc_phi2).append(" ")		// 72
								.append(Delta_phi).append(" ")		// 73
								.append(mc_Delta_phi).append(" ")	// 74
								.append(phih).append(" ")			// 75
								.append(mc_phih).append(" ")		// 76
								.append(phiR).append(" ")			// 77
								.append(mc_phiR).append(" ")		// 78
								.append(theta).append(" ")			// 79
								.append(mc_theta).append(" ")		// 80
								.append(Depolarization_A).append(" ")		// 81		
								.append(mc_Depolarization_A).append(" ")	// 82
								.append(Depolarization_B).append(" ")		// 83
								.append(mc_Depolarization_B).append(" ")	// 84
								.append(Depolarization_C).append(" ")		// 85
								.append(mc_Depolarization_C).append(" ")	// 86
								.append(Depolarization_V).append(" ")		// 87
								.append(mc_Depolarization_V).append(" ")	// 88
								.append(Depolarization_W).append(" ")		// 89
								.append(mc_Depolarization_W).append(" ")	// 90
								.append(matching_e_pid).append(" ")			// 91
								.append(matching_p1_pid).append(" ")		// 92
								.append(matching_p2_pid).append(" ")		// 93
								.append(mc_p1_parent).append(" ") 	// 94
								.append(mc_p2_parent).append("\n"); 	// 95

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
		println(); println();
		println("Set p1 PID = "+p1_Str+"\n");
		println("Set p2 PID = "+p2_Str+"\n");
		println("output file is: "+file);
		println();
	}
}
