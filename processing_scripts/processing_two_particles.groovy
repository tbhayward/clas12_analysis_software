/*
 * author Timothy B. Hayward
 * 
 * SIDIS hadron 
 */

// import CLAS12 physics classes
import org.jlab.io.hipo.*
import org.jlab.io.base.DataEvent
import org.jlab.clas.physics.*
import org.jlab.clas12.physics.*

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.* 
import analyzers.*

import groovy.io.FileType

// dilks CLAS QA analysis
import clasqa.QADB 

public static double phi_calculation(double x, double y) {
    double phi = Math.toDegrees(Math.atan2(x, y))
    phi -= 90
    if (phi < 0) phi += 360
    phi = 360 - phi
    return phi
}

public static double theta_calculation(double x, double y, double z) {
    double r = Math.sqrt(x*x + y*y + z*z)
    return Math.toDegrees(Math.acos(z / r))
}

public static void main(String[] args) {

    long startTime = System.currentTimeMillis()

    if (!args) {
        println("ERROR: Please enter a hipo file directory as the first argument")
        System.exit(1)
    }

    // collect all .hipo files under input directory
    def hipo_list = []
    (args[0] as File).eachFileRecurse(FileType.FILES) {
        if (it.name.endsWith('.hipo')) hipo_list << it
    }
    if (hipo_list.isEmpty()) {
        println("ERROR: No .hipo files found in ${args[0]}")
        System.exit(1)
    }

    // parse PID list from second argument, default to "211"
    String pidArg = args.length < 2 ? "211" : args[1]
    if (args.length < 2) {
        println("WARNING: No PID(s) provided; defaulting to '$pidArg'.")
    }
    List<Integer> pids = pidArg.split(',').collect { it.trim().toInteger() }
    println("Processing PID(s): $pids")

    // build one EventFilter per PID
    Map<Integer,EventFilter> filterMap = pids.collectEntries { pid ->
        [ (pid): new EventFilter("11:${pid}:X+:X-:Xn") ]
    }

    // output file (third argument)
    String output_file = args.length < 3 ? "hadron_dummy_out.txt" : args[2]
    if (args.length < 3) {
        println("WARNING: No output file specified; defaulting to '$output_file'.")
    }
    File file = new File(output_file)
    file.delete()

    // number of files to process (fourth argument)
    int n_files = (args.length < 4 ||
                   Integer.parseInt(args[3]) == 0 ||
                   Integer.parseInt(args[3]) > hipo_list.size())
                  ? hipo_list.size()
                  : Integer.parseInt(args[3])
    if (args.length < 4 ||
        Integer.parseInt(args[3]) == 0 ||
        Integer.parseInt(args[3]) > hipo_list.size()) {
        println("WARNING: Invalid or missing file-count; processing all ${hipo_list.size()} files.")
    }

    // beam energy (fifth argument)
    double beam_energy = args.length < 5 ? 10.6 : Double.parseDouble(args[4])
    if (args.length < 5) {
        println("No beam energy provided; defaulting to 10.6 GeV.")
    }

    // optional run override (sixth argument)
    Integer userProvidedRun = null
    if (args.length >= 6) {
        userProvidedRun = Integer.parseInt(args[5])
    } else {
        println("Run number not provided; pulling from hipo files.")
    }

    // declare physics and RICH variables
    int helicity, detector
    int fiducial_status, num_pos, num_neg, num_neutrals
    double e_p, e_theta, e_phi, p_phi, p_p, p_theta, open_angle
    double vz_e, vz_p
    double Q2, W, x, y, Mx2, t, tmin, z, xF, pT, eta, xi, phi
    double Depolarization_A, Depolarization_B, Depolarization_C
    double Depolarization_V, Depolarization_W
    int emilay, emico, emqua, best_PID
    float RQ, ReQ, el_logl, pi_logl, k_logl, pr_logl, best_ch, best_c2, best_RL, best_ntot
    int EB_pid

    // kinematic fitter
    GenericKinematicFitter fitter = new analysis_fitter(10.6041)

    // QADB setup
    QADB qa = new QADB()
    ['TotalOutlier','TerminalOutlier','MarginalOutlier',
     'SectorLoss','LowLiveTime','Misc','ChargeHigh',
     'ChargeNegative','ChargeUnknown','PossiblyNoBeam']
    .each { qa.checkForDefect(it) }
    [5046,5047,5051,5128,5129,5130,5158,5159,5160,5163,5165,5166,5167,5168,
     5169,5180,5181,5182,5183,5400,5448,5495,5496,5505,5567,5610,5617,5621,
     5623,6736,6737,6738,6739,6740,6741,6742,6743,6744,6746,6747,6748,6749,
     6750,6751,6753,6754,6755,6756,6757,16194,16089,16185,16308,16184,16307,16309]
    .each { qa.allowMiscBit(it) }

    StringBuilder batchLines = new StringBuilder()
    int num_events = 0, max_lines = 1000, lineCount = 0

    for (current_file in 0..<n_files) {
        println("\nOpening file ${current_file+1} of $n_files")
        HipoDataSource reader = new HipoDataSource()
        reader.open(hipo_list[current_file])

        while (reader.hasEvent()) {
            ++num_events
            if (num_events % 500000 == 0) print("processed: $num_events events. ")

            HipoDataEvent event = reader.getNextEvent()
            int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0)
            if (runnum > 16600 && runnum < 16700) break
            int evnum = event.getBank("RUN::config").getInt('event', 0)

            PhysicsEvent research_Event = fitter.getPhysicsEvent(event)

            boolean passQA = (runnum == 11 || runnum < 5020 || runnum > 16772 || qa.pass(runnum, evnum))
            if (runnum > 17768) passQA = false
            if (!passQA) continue

            // loop over each requested PID
            for (int p1_int : pids) {
                EventFilter f = filterMap[p1_int]
                if (!f.isValid(research_Event)) continue

                int num_p1 = research_Event.countByPid(p1_int)
                for (int current_p1 = 0; current_p1 < num_p1; current_p1++) {
                    BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false)
                    double energy = (runnum == 11) ? beam_energy : Eb.Eb()
                    TwoParticles variables = new TwoParticles(event, research_Event, p1_int, current_p1, energy)
                    if (!variables.channel_test(variables)) continue

                    fiducial_status = variables.get_fiducial_status()
                    helicity        = variables.get_helicity()
                    detector        = variables.get_detector()
                    num_pos         = variables.get_num_pos()
                    num_neg         = variables.get_num_neg()
                    num_neutrals    = variables.get_num_neutrals()

                    e_p        = variables.e_p()
                    e_theta    = variables.e_theta()
                    e_phi      = variables.e_phi()
                    vz_e       = variables.vz_e()
                    open_angle = variables.open_angle()
                    p_p        = variables.p_p()
                    p_theta    = variables.p_theta()
                    p_phi      = variables.p_phi()
                    vz_p       = variables.vz_p()

                    Q2 = variables.Q2()
                    W  = variables.W()
                    Mx2 = variables.Mx2()
                    x   = variables.x()
                    t   = variables.t()
                    tmin= variables.tmin()
                    y   = variables.y()
                    z   = variables.z()
                    xF  = variables.xF()
                    pT  = variables.pT()
                    eta = variables.eta()
                    xi  = variables.xi()
                    phi = variables.phi()

                    Depolarization_A = variables.Depolarization_A()
                    Depolarization_B = variables.Depolarization_B()
                    Depolarization_C = variables.Depolarization_C()
                    Depolarization_V = variables.Depolarization_V()
                    Depolarization_W = variables.Depolarization_W()

                    emilay   = variables.emilay()
                    emico    = variables.emico()
                    emqua    = variables.emqua()
                    best_PID = variables.best_PID()
                    println(best_PID);
                    RQ       = variables.RQ()
                    ReQ      = variables.ReQ()
                    el_logl  = variables.el_logl()
                    pi_logl  = variables.pi_logl()
                    k_logl   = variables.k_logl()
                    pr_logl  = variables.pr_logl()
                    best_ch  = variables.best_ch()
                    best_c2  = variables.best_c2()
                    best_RL  = variables.best_RL()
                    best_ntot= variables.best_ntot()

                    EB_pid = p1_int

                    // build and batch the output line
                    StringBuilder line = new StringBuilder()
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
                        .append(Depolarization_W).append(" ")
                        .append(emilay).append(" ")
                        .append(emico).append(" ")
                        .append(emqua).append(" ")
                        .append(best_PID).append(" ")
                        .append(RQ).append(" ")
                        .append(ReQ).append(" ")
                        .append(el_logl).append(" ")
                        .append(pi_logl).append(" ")
                        .append(k_logl).append(" ")
                        .append(pr_logl).append(" ")
                        .append(best_ch).append(" ")
                        .append(best_c2).append(" ")
                        .append(best_RL).append(" ")
                        .append(best_ntot).append(" ")
                        .append(EB_pid).append("\n")

                    batchLines.append(line.toString())
                    if (++lineCount >= max_lines) {
                        file.append(batchLines.toString())
                        batchLines.setLength(0)
                        lineCount = 0
                    }
                }
            }
        }
        reader.close()
    }

    // flush remaining
    if (batchLines.length() > 0) {
        file.append(batchLines.toString())
    }

    // final header mapping
    println(
      "1: fiducial_status, 2: num_pos, 3: num_neg, 4: num_neutrals, " +
      "5: runnum, 6: evnum, 7: helicity, 8: detector, " +
      "9: e_p, 10: e_theta, 11: e_phi, 12: vz_e, 13: open_angle, " +
      "14: p_p, 15: p_theta, 16: p_phi, 17: vz_p, " +
      "18: Q2, 19: W, 20: Mx2, 21: x, 22: t, 23: tmin, " +
      "24: y, 25: z, 26: xF, 27: pT, 28: xi, 29: eta, 30: phi, " +
      "31: DepA, 32: DepB, 33: DepC, 34: DepV, 35: DepW, " +
      "36: emilay, 37: emico, 38: emqua, 39: best_PID, 40: RQ, " +
      "41: ReQ, 42: el_logl, 43: pi_logl, 44: k_logl, 45: pr_logl, " +
      "46: best_ch, 47: best_c2, 48: best_RL, 49: best_ntot, 50: EB_pid"
    )

    println("Output file: $file")
    println("Elapsed time: ${System.currentTimeMillis() - startTime} ms")
}