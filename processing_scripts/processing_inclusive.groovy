/*
 * author Timothy B. Hayward
 *
 * DIS (inclusive e-p)
 */

// import CLAS12 physics classes
import org.jlab.io.hipo.*
import org.jlab.io.base.DataEvent
import org.jlab.clas.physics.*
import org.jlab.clas12.physics.*

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*
import analyzers.*

// filetype for gathering files in directory
import groovy.io.FileType

// dilks CLAS QA analysis
import clasqa.QADB

public static void main(String[] args) {

    long startTime = System.currentTimeMillis()

    // ~~~~~~~~~~~~~~~~ set up input parameters ~~~~~~~~~~~~~~~~ //
    if (!args) {
        println("ERROR: Please enter a hipo file directory as the first argument")
        System.exit(0)
    }

    def hipo_list = []
    (args[0] as File).eachFileRecurse(FileType.FILES) {
        if (it.name.endsWith('.hipo')) hipo_list << it
    }

    String output_file = args.length < 2 ? "inclusive_dummy_out.txt" : args[1]
    if (args.length < 2) println('WARNING: Specify an output file name. Set to "inclusive_dummy_out.txt".')

    File file = new File(output_file)
    file.delete()
    BufferedWriter writer = new BufferedWriter(new FileWriter(file))

    int n_files = args.length < 3 || Integer.parseInt(args[2]) == 0 || Integer.parseInt(args[2]) > hipo_list.size()
        ? hipo_list.size() : Integer.parseInt(args[2])

    if (args.length < 3 || Integer.parseInt(args[2]) == 0 || Integer.parseInt(args[2]) > hipo_list.size()) {
        println("WARNING: Number of files not specified, set to 0, or number too large.")
        println("Setting # of files to be equal to number of files in the directory.")
        println("There are $hipo_list.size files.")
    }

    double beam_energy = args.length < 4 ? 10.6 : Double.parseDouble(args[3])
    if (args.length < 4) {
        println("No beam energy provided, defaulting to run number based beam energy.")
        println("All MC will use 10.604 GeV. You must manually enter a beam energy to change this.")
    }

    Integer userProvidedRun = null
    if (args.length < 5) {
        println("Run number not provided, will pull from hipo files.")
        println("Think carefully about this if you are processing MC.")
    } else {
        userProvidedRun = Integer.parseInt(args[4])
    }

    // Allow for QADB override (usually meaning you're processing MC)
    Integer userProvidedOverride = 0
    if (args.length < 6) {
        println("No indication of QADB provided. Will use QADB.")
    } else {
        userProvidedOverride = Integer.parseInt(args[5])
    }

    // ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //
    // event-level vars we will fill from Inclusive
    int fiducial_status, helicity
    int num_pos, num_neg, num_neutrals
    double e_p, e_theta, e_phi, vz_e
    double Q2, W, y, Mx2, x
    double Depolarization_A, Depolarization_B, Depolarization_C, Depolarization_V, Depolarization_W

    // GenericKinematicFitter fitter = new analysis_fitter(10.6041);
    GenericKinematicFitter fitter = new monte_carlo_fitter(10.6041);
    EventFilter filter = new EventFilter("11:X+:X-:Xn")

    // setup QA database
    QADB qa = new QADB("latest")
    qa.checkForDefect('TotalOutlier')
    qa.checkForDefect('TerminalOutlier')
    qa.checkForDefect('MarginalOutlier')
    qa.checkForDefect('SectorLoss')
    qa.checkForDefect('LowLiveTime')
    qa.checkForDefect('Misc')
    qa.checkForDefect('ChargeHigh')
    qa.checkForDefect('ChargeNegative')
    qa.checkForDefect('ChargeUnknown')

    [
        5418, 5419,
        5443,
        5444,
        6616,
        6736, 6737, 6738, 6739, 6740, 6741,
        6742, 6743, 6744, 6746, 6747, 6748,
        6749, 6750, 6751, 6753, 6754, 6755,
        6756, 6757
    ].each { run ->
        qa.allowMiscBit(run)
    }

    StringBuilder batchLines = new StringBuilder()

    int num_events = 0
    int max_lines = 1000
    int lineCount = 0

    for (current_file in 0..<n_files) {
        println("\n Opening file " + Integer.toString(current_file + 1) + " out of " + n_files + ".\n")

        HipoDataSource reader = new HipoDataSource()
        reader.open(hipo_list[current_file])

        while (reader.hasEvent()) {
            ++num_events
            if (num_events % 1000000 == 0) print("processed: " + num_events + " events. ")

            HipoDataEvent event = reader.getNextEvent()

            int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0)
            if (runnum > 16600 && runnum < 16700) break

            int evnum = event.getBank("RUN::config").getInt('event', 0)

            PhysicsEvent research_Event = fitter.getPhysicsEvent(event)

            boolean process_event = filter.isValid(research_Event) &&
                (runnum == 11 ||
                 userProvidedOverride == 1 ||
                 qa.pass(runnum, evnum))

            if (runnum == 5247) process_event = false
            if (runnum == 5345) process_event = false

            if (runnum == 5158 ||
                runnum == 5163 ||
                runnum == 5181 ||
                runnum == 5519 ||
                runnum == 5528 ||
                runnum == 5627 ||
                runnum == 3355 ||
                runnum == 3404 ||
                runnum == 3408 ||
                runnum == 3449 ||
                runnum == 3490 ||
                runnum == 3499 ||
                runnum == 3500 ||
                runnum == 3505 ||
                runnum == 3526 ||
                runnum == 3527 ||
                runnum == 3528 ||
                runnum == 3529 ||
                runnum == 3530 ||
                runnum == 3531 ||
                runnum == 3532 ||
                runnum == 3533 ||
                runnum == 3534 ||
                runnum == 3535 ||
                runnum == 3536 ||
                runnum == 3538 ||
                runnum == 3540 ||
                runnum == 3544 ||
                runnum == 3545 ||
                runnum == 3547 ||
                runnum == 3548 ||
                runnum == 3709 ||
                runnum == 3736 ||
                runnum == 3793 ||
                runnum == 3800 ||
                runnum == 3801 ||
                runnum == 3807 ||
                runnum == 3508 ||
                runnum == 3808 ||
                runnum == 3809 ||
                runnum == 3810 ||
                runnum == 3813 ||
                runnum == 3698 ||
                runnum == 3814 ||
                runnum == 3815 ||
                runnum == 3817 ||
                runnum == 4018 ||
                runnum == 4059 ||
                runnum == 4142 ||
                runnum == 4145 ||
                runnum == 4146 ||
                runnum == 4159 ||
                runnum == 4160 ||
                runnum == 4162 ||
                runnum == 4163 ||
                runnum == 4176 ||
                runnum == 4209 ||
                runnum == 4227 ||
                runnum == 4246 ||
                runnum == 4252 ||
                runnum == 4325 ||
                runnum == 3867 ||
                runnum == 3877 ||
                runnum == 3882 ||
                runnum == 3927 ||
                runnum == 3951 ||
                runnum == 3953 ||
                runnum == 3965 ||
                runnum == 3967 ||
                runnum == 3968 ||
                runnum == 3499 ||
                runnum == 3712 ||
                runnum == 3801 ||
                runnum == 3807 ||
                runnum == 3808 ||
                runnum == 3267 ||
                runnum == 3879 ||
                runnum == 3923 ||
                runnum == 3929 ||
                runnum == 3947) {
                process_event = false
            }

            // --- Toggle here ---
            // false -> baseline (no inverse-ISR)
            // true  -> enable inverse-ISR sampling (subtract R from q inside Inclusive)
            BeamEnergy Eb = new BeamEnergy(research_Event, runnum, false)

            double Ebeam = (runnum == 11) ? beam_energy : Eb.Eb()

            // If radiative sampling requested, prepare one-shot correction for Inclusive
            boolean applyInverseISR = Eb.isRadiativeApplied() && Eb.getEgammaGeV() > 0.0
            double Egamma = 0.0
            double isrTheta = 0.0
            double isrPhi = 0.0

            if (process_event) {
                if (applyInverseISR) {
                    Egamma = Eb.getEgammaGeV()
                    isrTheta = analyzers.ISRThetaKernel.sampleThetaRad(Egamma)
                    isrPhi = 2.0 * Math.PI * Math.random()
                    Inclusive.setNextInverseISRPhoton(Egamma, isrTheta, isrPhi)
                }

                Inclusive variables = new Inclusive(event, research_Event, Ebeam)

                if (Inclusive.channel_test(variables)) {
                    fiducial_status = variables.get_fiducial_status()
                    helicity = variables.get_helicity()
                    num_pos = variables.get_num_pos()
                    num_neg = variables.get_num_neg()
                    num_neutrals = variables.get_num_neutrals()

                    e_p = variables.e_p()
                    e_theta = variables.e_theta()
                    e_phi = variables.e_phi()
                    vz_e = variables.vz_e()

                    Q2 = variables.Q2()
                    W = variables.W()
                    x = variables.x()
                    y = variables.y()
                    Mx2 = variables.Mx2()

                    Depolarization_A = variables.Depolarization_A()
                    Depolarization_B = variables.Depolarization_B()
                    Depolarization_C = variables.Depolarization_C()
                    Depolarization_V = variables.Depolarization_V()
                    Depolarization_W = variables.Depolarization_W()

                    double isrTheta_deg = Math.toDegrees(isrTheta)
                    double isrPhi_deg = Math.toDegrees(isrPhi)

                    StringBuilder line = new StringBuilder()
                    line.append(fiducial_status).append(" ")
                        .append(num_pos).append(" ")
                        .append(num_neg).append(" ")
                        .append(num_neutrals).append(" ")
                        .append(runnum).append(" ")
                        .append(evnum).append(" ")
                        .append(helicity).append(" ")
                        .append(e_p).append(" ")
                        .append(e_theta).append(" ")
                        .append(e_phi).append(" ")
                        .append(vz_e).append(" ")
                        .append(applyInverseISR ? Egamma : 0.0).append(" ")
                        .append(applyInverseISR ? isrTheta_deg : 0.0).append(" ")
                        .append(applyInverseISR ? isrPhi_deg : 0.0).append(" ")
                        .append(Q2).append(" ")
                        .append(W).append(" ")
                        .append(Mx2).append(" ")
                        .append(x).append(" ")
                        .append(y).append(" ")
                        .append(Depolarization_A).append(" ")
                        .append(Depolarization_B).append(" ")
                        .append(Depolarization_C).append(" ")
                        .append(Depolarization_V).append(" ")
                        .append(Depolarization_W).append("\n")

                    batchLines.append(line.toString())
                    lineCount++

                    if (lineCount >= max_lines) {
                        writer.write(batchLines.toString())
                        batchLines.setLength(0)
                        lineCount = 0
                    }
                }
            }
        }

        reader.close()

        if (batchLines.length() > 0) {
            writer.write(batchLines.toString())
            batchLines.setLength(0)
            lineCount = 0
        }

        println(
            "1:  fiducial_status,  " +
            "2:  num_pos,          " +
            "3:  num_neg,          " +
            "4:  num_neutrals,     " +
            "5:  runnum,           " +
            "6:  evnum,            " +
            "7:  helicity,         " +
            "8:  e_p,              " +
            "9:  e_theta,          " +
            "10: e_phi,            " +
            "11: vz_e,             " +
            "12: Egamma (GeV),     " +
            "13: isrTheta (deg),   " +
            "14: isrPhi (deg),     " +
            "15: Q2,               " +
            "16: W,                " +
            "17: Mx2,              " +
            "18: x,                " +
            "19: y,                " +
            "20: DepA,             " +
            "21: DepB,             " +
            "22: DepC,             " +
            "23: DepV,             " +
            "24: DepW"
        )

        println("output text file is: $file")
    }

    writer.close()

    long endTime = System.currentTimeMillis()
    println("Elapsed time: ${endTime - startTime} ms")
}