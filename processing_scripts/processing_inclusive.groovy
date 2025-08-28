/*
 * author Timothy B. Hayward
 *
 * DIS (inclusive e−p)
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
    (args[0] as File).eachFileRecurse(FileType.FILES) { if (it.name.endsWith('.hipo')) hipo_list << it }

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

    // ~~~~~~~~~~~~~~~~ prepare physics analysis ~~~~~~~~~~~~~~~~ //
    // event-level vars we’ll fill from Inclusive
    int fiducial_status, helicity
    int num_pos, num_neg, num_neutrals
    double e_p, e_theta, e_phi, vz_e
    double Q2, W, y, Mx2, x
    double Depolarization_A, Depolarization_B, Depolarization_C, Depolarization_V, Depolarization_W

    GenericKinematicFitter fitter = new analysis_fitter(10.6041)
    EventFilter filter = new EventFilter("11:X+:X-:Xn")

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
    qa.checkForDefect('PossiblyNoBeam')
    [
        5046, 5047, 5051, 5128, 5129, 5130, 5158, 5159,
        5160, 5163, 5165, 5166, 5167, 5168, 5169, 5180,
        5181, 5182, 5183, 5400, 5448, 5495, 5496, 5505,
        5567, 5610, 5617, 5621, 5623, 6736, 6737, 6738,
        6739, 6740, 6741, 6742, 6743, 6744, 6746, 6747,
        6748, 6749, 6750, 6751, 6753, 6754, 6755, 6756,
        6757,
        16194, 16089, 16185, 16308, 16184, 16307, 16309,
        16872, 16975,
        17763, 17764, 17765, 17766, 17767, 17768,
        17179, 17180, 17181, 17182, 17183, 17188, 17189,
        17252
    ].each { run -> qa.allowMiscBit(run) }

    StringBuilder batchLines = new StringBuilder()

    int num_events = 0
    int max_lines = 1000
    int lineCount = 0

    for (current_file in 0..<n_files) {
        println("\n Opening file " + Integer.toString(current_file + 1) + " out of " + n_files + ".\n")

        HipoDataSource reader = new HipoDataSource()
        reader.open(hipo_list[current_file])

        HipoDataEvent event = reader.getNextEvent()

        while (reader.hasEvent()) {
            ++num_events
            if (num_events % 1000000 == 0) print("processed: " + num_events + " events. ")

            event = reader.getNextEvent()
            int runnum = userProvidedRun ?: event.getBank("RUN::config").getInt('run', 0)
            if (runnum > 16600 && runnum < 16700) break // Hall C bleedthrough
            int evnum = event.getBank("RUN::config").getInt('event', 0)

            PhysicsEvent research_Event = fitter.getPhysicsEvent(event)

            boolean process_event = filter.isValid(research_Event) &&
                    (runnum == 11 || runnum < 5020 || qa.pass(runnum, evnum))
            if (runnum > 17768) process_event = false
            if ([17331, 16987, 17079, 17190, 17639].contains(runnum)) process_event = false
            if ([16850, 16851, 16852, 16855, 16879].contains(runnum)) process_event = false

            // --- Toggle here ---
            // false -> baseline (no inverse-ISR)
            // true  -> enable inverse-ISR sampling (subtract R from q inside Inclusive)
            BeamEnergy Eb = new BeamEnergy(research_Event, runnum, /*isRadiative=*/false)

            double Ebeam = (runnum == 11) ? beam_energy : Eb.Eb()

            // If radiative sampling requested, prepare one-shot correction for Inclusive
            boolean applyInverseISR = Eb.isRadiativeApplied() && Eb.getEgammaGeV() > 0.0
            double Egamma = 0.0, isrTheta = 0.0, isrPhi = 0.0

            if (process_event) {
                if (applyInverseISR) {
                    Egamma   = Eb.getEgammaGeV()
                    isrTheta = analyzers.ISRThetaKernel.sampleThetaRad(Egamma) // radians
                    isrPhi   = 2.0 * Math.PI * Math.random()
                    Inclusive.setNextInverseISRPhoton(Egamma, isrTheta, isrPhi)
                }

                Inclusive variables = new Inclusive(event, research_Event, Ebeam)

                if (Inclusive.channel_test(variables)) {
                    fiducial_status = variables.get_fiducial_status()
                    helicity        = variables.get_helicity()
                    num_pos         = variables.get_num_pos()
                    num_neg         = variables.get_num_neg()
                    num_neutrals    = variables.get_num_neutrals()

                    // lab electron kinematics
                    e_p    = variables.e_p()
                    e_theta= variables.e_theta()
                    e_phi  = variables.e_phi()
                    vz_e   = variables.vz_e()

                    // DIS (already corrected if toggle on)
                    Q2 = variables.Q2()
                    W  = variables.W()
                    x  = variables.x()
                    y  = variables.y()
                    Mx2 = variables.Mx2()

                    Depolarization_A = variables.Depolarization_A()
                    Depolarization_B = variables.Depolarization_B()
                    Depolarization_C = variables.Depolarization_C()
                    Depolarization_V = variables.Depolarization_V()
                    Depolarization_W = variables.Depolarization_W()

                    double isrTheta_deg = Math.toDegrees(isrTheta)
                    double isrPhi_deg   = Math.toDegrees(isrPhi)

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
                        // NEW: simulated photon (R). Zeros in baseline mode.
                        .append(applyInverseISR ? Egamma       : 0.0).append(" ")
                        .append(applyInverseISR ? isrTheta_deg : 0.0).append(" ")
                        .append(applyInverseISR ? isrPhi_deg   : 0.0).append(" ")
                        // corrected (or baseline) DIS
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
                        file.append(batchLines.toString())
                        batchLines.setLength(0)
                        lineCount = 0
                    }
                } // channel test
            } // process_event
        } // while
        reader.close()

        if (batchLines.length() > 0) {
            file.append(batchLines.toString())
            batchLines.setLength(0)
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
    } // files

    writer.close()

    long endTime = System.currentTimeMillis()
    println("Elapsed time: ${endTime - startTime} ms")
}