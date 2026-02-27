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
	double phi = Math.toDegrees(Math.atan2(x,y));
	phi = phi - 90;
	if (phi < 0) phi = 360 + phi;
	phi = 360 - phi;
	return phi;	
}

public static double theta_calculation (double x, double y, double z) {
	double r = Math.sqrt(x*x + y*y + z*z);
	return (180/Math.PI)*Math.acos(z/r);
}

public static void main(String[] args) {

	long startTime = System.currentTimeMillis();

	if (!args) {
	    println("ERROR: Please enter a hipo file directory as the first argument");
	    System.exit(0);
	}

	def hipo_list = []
	(args[0] as File).eachFileRecurse(FileType.FILES) {
		if (it.name.endsWith('.hipo')) hipo_list << it
	}

	String p1_Str = args.length < 2 ? "211" : args[1];
	if (args.length < 2) println("WARNING: Specify a PDG PID for p1! Set to pi+ (211).");
	println("Set p1 PID = $p1_Str");
	int p1_int = p1_Str.toInteger();

	String output_file = args.length < 3 ? "hadron_dummy_out.txt" : args[2];
	if (args.length < 3)
	    println("WARNING: Specify an output file name. Set to \"hadron_dummy_out.txt\".");
	File file = new File(output_file);
	file.delete();

	int n_files = args.length < 4 || Integer.parseInt(args[3]) > hipo_list.size()
	    ? hipo_list.size() : Integer.parseInt(args[3]);

	double beam_energy = args.length < 5 ? 10.6 : Double.parseDouble(args[4]);

	Integer userProvidedRun = null
	if (args.length >= 6)
		userProvidedRun = Integer.parseInt(args[5]);

	GenericKinematicFitter research_fitter = new analysis_fitter(10.6041);
	GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6041);
	EventFilter filter = new EventFilter("11:"+p1_Str+":"+":X+:X-:Xn");

	StringBuilder batchLines = new StringBuilder();
	int num_events = 0;
	int max_lines = 1000;
	int lineCount = 0;

	for (current_file in 0..<n_files) {

		println("\n Opening file "+(current_file+1)+" out of "+n_files+".\n"); 
		
		HipoDataSource reader = new HipoDataSource();
		reader.open(hipo_list[current_file]);

		while(reader.hasEvent()){

			++num_events;
			event = reader.getNextEvent();

		    int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0);

		    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);
		    PhysicsEvent mc_Event = mc_fitter.getPhysicsEvent(event);

			if (!filter.isValid(research_Event)) continue;

			HipoDataBank recParticleBank = (HipoDataBank) event.getBank("REC::Particle");
			HipoDataBank recMatchBank    = event.hasBank("MC::RecMatch") ?
				(HipoDataBank) event.getBank("MC::RecMatch") : null;
			HipoDataBank mcBank          = event.hasBank("MC::Particle") ?
				(HipoDataBank) event.getBank("MC::Particle") : null;
			HipoDataBank lundBank        = event.hasBank("MC::Lund") ?
				(HipoDataBank) event.getBank("MC::Lund") : null;

			int num_p1 = research_Event.countByPid(p1_int);

			for (int current_p1 = 0; current_p1 < num_p1; current_p1++) {

				Particle exp_e  = research_Event.getParticleByPid(11,0);
				Particle exp_p1 = research_Event.getParticleByPid(p1_int,current_p1);

				BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false);
				double energy = (runnum == 11) ? beam_energy : Eb.Eb();

				TwoParticles variables = new TwoParticles(event, research_Event, 
					p1_int, current_p1, energy);
				TwoParticles mc_variables = new TwoParticles(event, mc_Event, 
					p1_int, current_p1, energy);

				if (!variables.channel_test(variables)) continue;

				// --- DATA & MC VARIABLES (UNCHANGED) ---
				double e_p = variables.e_p();
				double mc_e_p = mc_variables.e_p();
				double e_theta = variables.e_theta();
				double mc_e_theta = mc_variables.e_theta();
				double e_phi = variables.e_phi();
				double mc_e_phi = mc_variables.e_phi();
				double vz_e = variables.vz_e();
				double mc_vz_e = mc_variables.vz_e();

				double p_p = variables.p_p();
				double mc_p_p = mc_variables.p_p();
				double p_theta = variables.p_theta();
				double mc_p_theta = mc_variables.p_theta();
				double p_phi = variables.p_phi();
				double mc_p_phi = mc_variables.p_phi();
				double vz_p = variables.vz_p();
				double mc_vz_p = mc_variables.vz_p();

				double Q2 = variables.Q2();
				double mc_Q2 = mc_variables.Q2();
				double W = variables.W();
				double mc_W = mc_variables.W();
				double Mx = variables.Mx2();
				double mc_Mx = mc_variables.Mx();
				double Mx2 = variables.Mx2();
				double mc_Mx2 = mc_variables.Mx2();

				double x = variables.x();
				double mc_x = mc_variables.x();
				double y = variables.y();
				double mc_y = mc_variables.y();
				double z = variables.z();
				double mc_z = mc_variables.z();
				double xF = variables.xF();
				double mc_xF = mc_variables.xF();
				double pT = variables.pT();
				double mc_pT = mc_variables.pT();
				double eta = variables.eta();
				double mc_eta = mc_variables.eta();
				double xi = variables.xi();
				double mc_xi = mc_variables.xi();

				double trento_phi = variables.phi();
				double mc_trento_phi = mc_variables.phi();

				double Depolarization_A = variables.Depolarization_A();
				double mc_Depolarization_A = mc_variables.Depolarization_A();
				double Depolarization_B = variables.Depolarization_B();
				double mc_Depolarization_B = mc_variables.Depolarization_B();
				double Depolarization_C = variables.Depolarization_C();
				double mc_Depolarization_C = mc_variables.Depolarization_C();
				double Depolarization_V = variables.Depolarization_V();
				double mc_Depolarization_V = mc_variables.Depolarization_V();
				double Depolarization_W = variables.Depolarization_W();
				double mc_Depolarization_W = mc_variables.Depolarization_W();

				// =======================================
				// NEW TRUTH MATCHING USING MC::RecMatch
				// =======================================

				int matching_e_pid = 0;
				int matching_p1_pid = 0;
				int mc_p1_parent = 0;

				if (recMatchBank != null && mcBank != null && lundBank != null) {

					for (int row = 0; row < recMatchBank.rows(); row++) {

						int rec_pindex = recMatchBank.getInt("pindex", row);
						int mcindex    = recMatchBank.getInt("mcindex", row);

						if (mcindex < 0 || mcindex >= mcBank.rows()) continue;

						int rec_pid = recParticleBank.getInt("pid", rec_pindex);
						int truth_pid = mcBank.getInt("pid", mcindex);

						if (rec_pid == 11 && matching_e_pid == 0)
							matching_e_pid = truth_pid;

						if (rec_pid == p1_int && matching_p1_pid == 0) {

							matching_p1_pid = truth_pid;

							if (mcindex < lundBank.rows()) {
								int parent_index = lundBank.getInt("parent", mcindex);
								if (parent_index > 0) {
									int parent_row = parent_index - 1;
									if (parent_row >= 0 && parent_row < lundBank.rows())
										mc_p1_parent = lundBank.getInt("pid", parent_row);
								}
							}
						}
					}
				}

				// =======================================
				// END TRUTH MATCHING
				// =======================================

				StringBuilder line = new StringBuilder();
				line.append(e_p).append(" ")
					.append(mc_e_p).append(" ")
					.append(e_theta).append(" ")
					.append(mc_e_theta).append(" ")
					.append(e_phi).append(" ")
					.append(mc_e_phi).append(" ")
					.append(vz_e).append(" ")
					.append(mc_vz_e).append(" ")
					.append(p_p).append(" ")
					.append(mc_p_p).append(" ")
					.append(p_theta).append(" ")
					.append(mc_p_theta).append(" ")
					.append(p_phi).append(" ")
					.append(mc_p_phi).append(" ")
					.append(vz_p).append(" ")
					.append(mc_vz_p).append(" ")
					.append(Q2).append(" ")
					.append(mc_Q2).append(" ")
					.append(W).append(" ")
					.append(mc_W).append(" ")
					.append(Mx).append(" ")
					.append(mc_Mx).append(" ")
					.append(Mx2).append(" ")
					.append(mc_Mx2).append(" ")
					.append(x).append(" ")
					.append(mc_x).append(" ")
					.append(y).append(" ")
					.append(mc_y).append(" ")
					.append(z).append(" ")
					.append(mc_z).append(" ")
					.append(xF).append(" ")
					.append(mc_xF).append(" ")
					.append(pT).append(" ")
					.append(mc_pT).append(" ")
					.append(xi).append(" ")
					.append(mc_xi).append(" ")
					.append(eta).append(" ")
					.append(mc_eta).append(" ")
					.append(trento_phi).append(" ")
					.append(mc_trento_phi).append(" ")
					.append(Depolarization_A).append(" ")
					.append(mc_Depolarization_A).append(" ")
					.append(Depolarization_B).append(" ")
					.append(mc_Depolarization_B).append(" ")
					.append(Depolarization_C).append(" ")
					.append(mc_Depolarization_C).append(" ")
					.append(Depolarization_V).append(" ")
					.append(mc_Depolarization_V).append(" ")
					.append(Depolarization_W).append(" ")
					.append(mc_Depolarization_W).append(" ")
					.append(matching_e_pid).append(" ")
					.append(matching_p1_pid).append(" ")
					.append(mc_p1_parent).append("\n");

				batchLines.append(line.toString());
				lineCount++;

				if (lineCount >= max_lines) {
					file.append(batchLines.toString());
					batchLines.setLength(0);
					lineCount = 0;
				}
			}
		}
		reader.close();
	}

	if (batchLines.length() > 0)
		file.append(batchLines.toString());

	println("Set p1 PID = "+p1_Str);
	println("output file is: "+file);

	long endTime = System.currentTimeMillis();
	println("Elapsed time: "+(endTime-startTime)+" ms");
}