/*
 * author Timothy B. Hayward
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

public static double get_mc_event_weight(HipoDataEvent event) {
	double weight = 1.0;

	if (event.hasBank("MC::Event")) {
		def mc_event_bank = event.getBank("MC::Event");

		if (mc_event_bank.rows() > 0) {
			weight = mc_event_bank.getFloat("weight", 0);
		}
	}

	return weight;
}


// Generated and reconstructed quantities are evaluated with the SAME analyzer.
// Entries remain aligned by construction (one per HIPO event), not RUN::config numbers.
static List<String> pairedNames() { return ["fiducial_status","num_pos","num_neg","num_neutrals","runnum","evnum","helicity","detector1","detector2","e_p","e_theta","e_phi","vz_e","p1_p","p1_theta","p1_phi","vz_p1","p2_p","p2_theta","p2_phi","vz_p2","open_angle_ep","open_angle_ep1","open_angle_ep2","open_angle_p1p2","Q2","W","Mx2","Mx2_1","Mx2_2","x","t","t1","t2","tmin","y","z","z1","z2","Mh","xF","xF1","xF2","pT","pT1","pT2","pTpT","xi","xi1","xi2","eta","eta1","eta2","Delta_eta","eta1_gN","eta2_gN","phi1","phi2","Delta_phi","phih","phiR","theta","DepA","DepB","DepC","DepV","DepW","Mh_gammagamma","detector_gamma1","detector_gamma2","open_angle_egamma1","open_angle_egamma2","gamma_phi1","gamma_phi2","Emiss2","theta_pi0_pi0","pTmiss","weight"] }
static List<String> pairedMethods() { return ["get_fiducial_status","get_num_pos","get_num_neg","get_num_neutrals","__runnum","__evnum","get_helicity","get_detector1","get_detector2","e_p","e_theta","e_phi","vz_e","p1_p","p1_theta","p1_phi","vz_p1","p2_p","p2_theta","p2_phi","vz_p2","open_angle_ep","open_angle_ep1","open_angle_ep2","open_angle_p1p2","Q2","W","Mx2","Mx2_1","Mx2_2","x","t","t1","t2","tmin","y","z","z1","z2","Mh","xF","xF1","xF2","pT","pT1","pT2","pTpT","xi","xi1","xi2","eta","eta1","eta2","Delta_eta","eta1_gN","eta2_gN","phi1","phi2","Delta_phi","phih","phiR","theta","Depolarization_A","Depolarization_B","Depolarization_C","Depolarization_V","Depolarization_W","__pair_Mh","__pair_detector1","__pair_detector2","__pair_open1","__pair_open2","__pair_phi1","__pair_phi2","Emiss2","theta_gamma_gamma","pTmiss","__weight"] }
// Both kinematic sets are constructed through the same ThreeParticles analyzer.
// For the pi0-specific fields, pairAnalyzer is the two-photon analyzer.
static List<Double> pairedKinematics(def analyzer, def pairAnalyzer, int run, int ev, double w, boolean generated) {
    if (analyzer == null || pairAnalyzer == null) return pairedNames().collect { Double.NaN }
    List<Double> values = []
    for (int j = 0; j < pairedNames().size(); j++) {
        String name = pairedNames()[j]
        String method = pairedMethods()[j]
        if (name == 'runnum') { values.add((double)run); continue }
        if (name == 'evnum') { values.add((double)ev); continue }
        if (name == 'weight') { values.add(w); continue }
        if (method.startsWith('__pair_')) {
            if (method == '__pair_Mh') { values.add((double)pairAnalyzer.Mh()); continue }
            if (method == '__pair_detector1') { values.add(generated ? Double.NaN : (double)pairAnalyzer.get_detector1()); continue }
            if (method == '__pair_detector2') { values.add(generated ? Double.NaN : (double)pairAnalyzer.get_detector2()); continue }
            if (method == '__pair_open1') { values.add((double)pairAnalyzer.open_angle_ep1()); continue }
            if (method == '__pair_open2') { values.add((double)pairAnalyzer.open_angle_ep2()); continue }
            if (method == '__pair_phi1') { values.add((double)pairAnalyzer.phi1()); continue }
            if (method == '__pair_phi2') { values.add((double)pairAnalyzer.phi2()); continue }
        }
        if (generated && ['fiducial_status','detector1','detector2'].contains(name)) {
            values.add(Double.NaN); continue
        }
        def value = analyzer."${method}"()
        values.add(value instanceof Number ? ((Number)value).doubleValue() : Double.NaN)
    }
    return values
}

// Select a unique photon pair using the same gamma-gamma mass window as the
// legacy exclusive pi0 processing. The reconstructed pair additionally uses
// the original photon detector exclusion. This does not apply detector cuts
// to generator truth.
static Map choosePhotonPair(def event, def physicsEvent, double energy, boolean generated) {
    if (physicsEvent == null || physicsEvent.countByPid(11) < 1 ||
        physicsEvent.countByPid(2212) < 1 || physicsEvent.countByPid(22) < 2) return null
    Map best = null
    double bestDistance = Double.POSITIVE_INFINITY
    int n = physicsEvent.countByPid(22)
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            def pair = new ThreeParticles(event, physicsEvent, 22, i, 22, j, energy, generated)
            double mass = pair.Mh()
            if (!Double.isFinite(mass) || mass < 0.11d || mass > 0.16d) continue
            if (!generated && (pair.get_detector1() == 2 || pair.get_detector2() == 2)) continue
            double distance = Math.abs(mass - 0.1349768d)
            if (distance < bestDistance) {
                bestDistance = distance
                best = [first:i, second:j, analyzer:pair]
            }
        }
    }
    return best
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
	File sourceInput = new File(args[0]);
    if (sourceInput.isFile() && sourceInput.name.endsWith('.hipo')) hipo_list << sourceInput;
    else if (sourceInput.isDirectory()) sourceInput.eachFileRecurse(FileType.FILES) { if (it.name.endsWith('.hipo')) hipo_list << it };
    else throw new IllegalArgumentException('Input must be a HIPO file or directory: ' + args[0]);

    // Deterministic ordering matches the campaign progress manifest.
    hipo_list.sort { a, b -> a.absolutePath <=> b.absolutePath }
    if (hipo_list.size() != 1) throw new IllegalArgumentException("Paired pilot requires exactly one HIPO file");

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

	// Allow for QADB override (usually meaning you're processing MC)
	Integer userProvidedOverride = 0;
	if (args.length < 6) {
		println("No indication of QADB provided. Will use QADB.");
	} else {
		userProvidedOverride = Integer.parseInt(args[5]);
	}


    // Optional MC pilot controls: maxEvents (0=unlimited), correction period.
    long maxEvents = args.length > 6 ? Long.parseLong(args[6]) : 0L;
    Map correctionRuns = [fa18_inb:5000, fa18_out:5500, sp18_inb:3600,
                          sp18_out:3900, sp19_inb:6700];
    String correctionPeriod = args.length > 7 ? args[7].toLowerCase() : "none";
    if (correctionPeriod != "none" && !correctionRuns.containsKey(correctionPeriod)) {
        System.err.println("ERROR: invalid MC correction period: " + correctionPeriod);
        System.exit(2);
    }
    Integer correctionRun = correctionPeriod == "none" ? null : (Integer) correctionRuns[correctionPeriod];
    println("MC pi0 correction period: ${correctionPeriod}; correction run: ${correctionRun}; MC QADB bypass enabled.");
    // Event-aligned paired MC: stream one record per processed HIPO event.
    // No particle sidecar and no on-disk event TSV. Caller provides a FIFO path.
    BufferedWriter diagnosticEvents = new BufferedWriter(new FileWriter(output_file + ".events.tsv"), 1024*1024);
    List<String> eventColumns = ['event_sequence','weight','rec_kinematics_valid','gen_kinematics_valid'];
    eventColumns.addAll(pairedNames().collect { 'gen_' + it });
    eventColumns.addAll(pairedNames().collect { 'rec_' + it });
    diagnosticEvents.write(eventColumns.join(' ') + "\n");

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
	double weight;

	// load my kinematic fitter/PID
	// Initialize the fitter
	GenericKinematicFitter fitter;

	// Uncomment the desired fitter
	fitter = new analysis_fitter(beam_energy);
    monte_carlo_fitter truthFitter = new monte_carlo_fitter(beam_energy); 
	// fitter = new monte_carlo_fitter(10.6041);
	// fitter = new event_builder_fitter(10.6041);  

    if (correctionRun != null) ((analysis_fitter)fitter).setCorrectionRunOverride(correctionRun);

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
	
	// MC processing deliberately does not instantiate or query QADB.

	// create a StringBuilder for accumulating lines
	StringBuilder batchLines = new StringBuilder();

	int num_events = 0;
	int max_lines = 1000;
	int lineCount = 0;
	for (current_file in 0..<n_files) {
		// limit to a certain number of files defined by n_files
		println("[HIPO_FILE_START] " + hipo_list[current_file].absolutePath);
		println("\n Opening file "+Integer.toString(current_file+1)
			+" out of "+n_files+".\n"); 
		
		HipoDataSource reader = new HipoDataSource();
		reader.open(hipo_list[current_file]); // open next hipo file
		HipoDataEvent event = reader.getNextEvent(); 

		while (reader.hasEvent()) {
            if (maxEvents > 0 && num_events >= maxEvents) break;
		    ++num_events;
		    if (num_events % 1000000 == 0) { // not necessary, just updates output
		        print("processed: " + num_events + " events. ");
		    }

		    // get run and event numbers
		    event = reader.getNextEvent();
		    // collect info for QA
		    int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0);
		    int evnum = event.getBank("RUN::config").getInt('event', 0);
		    float torus = event.getBank("RUN::config").getFloat("torus", 0);
		    weight = get_mc_event_weight(event);

		    PhysicsEvent research_Event = fitter.getPhysicsEvent(event);
            PhysicsEvent generated_Event = truthFitter.getPhysicsEvent(event);
            Map generatedPair = choosePhotonPair(event, generated_Event, beam_energy, true);
            Map reconstructedPair = filter.isValid(research_Event) ?
                choosePhotonPair(event, research_Event, beam_energy, false) : null;
            // Use the composite pi0 (PID 111) created by the corresponding fitter.
            boolean truthValid = generatedPair != null && generated_Event.countByPid(111) > 0;
            boolean recoValid = reconstructedPair != null && research_Event.countByPid(111) > 0;
            def genAnalyzer = truthValid ? new ThreeParticles(event, generated_Event, 2212, 0, 111, 0, beam_energy, true) : null;
            def recAnalyzer = recoValid ? new ThreeParticles(event, research_Event, 2212, 0, 111, 0, beam_energy) : null;
            List<Double> genValues = pairedKinematics(genAnalyzer, generatedPair?.analyzer, runnum, evnum, weight, true);
            List<Double> recValues = pairedKinematics(recAnalyzer, reconstructedPair?.analyzer, runnum, evnum, weight, false);

            // The entry is created in the SAME HIPO event as both analyzers.
            // Missing reconstructed candidates are represented by NaN branches.
            diagnosticEvents.write(([num_events, weight, (recoValid ? 1 : 0), (truthValid ? 1 : 0)] + genValues + recValues).join(' ') + "\n");

	    	boolean process_event = filter.isValid(research_Event);
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
	                if (torus == -1 && Q2 < 1.456) { continue; }

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
	                    .append(pTmiss).append(" ")
	                    .append(weight).append("\n");

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
        println("[HIPO_FILE_DONE] " + hipo_list[current_file].absolutePath + " events_total=" + num_events);
        if (maxEvents > 0 && num_events >= maxEvents) break;
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
	    "75: Emiss2, 76: theta_pi0_pi0, 77: pTmiss, 78: weight.");

		println("output text file is: $file");
	}

    diagnosticEvents.close();

	writer.close();

	// End time
	long endTime = System.currentTimeMillis()
	// Calculate the elapsed time
	long elapsedTime = endTime - startTime
	// Print the elapsed time in milliseconds
	println("Elapsed time: ${elapsedTime} ms");

}