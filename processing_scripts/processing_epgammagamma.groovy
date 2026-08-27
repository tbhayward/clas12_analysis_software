/*
 * Unified e'p'gamma gamma X processing for data and reconstructed MC.
 *
 * No pi0 mass cut is applied. Every UNIQUE photon pair is written once
 * (gamma2_index > gamma1_index), and Mh_gammagamma = Mh23 is saved explicitly.
 *
 * This deliberately retains both daughter photons separately. It is therefore
 * useful for:
 *   1) pi0/eta/continuum studies,
 *   2) tag-and-probe photon efficiency,
 *   3) BDT truth studies using each photon's MC::Lund parent.
 *
 * MC truth is event-by-event optional and is activated by MC::Particle.
 */

import org.jlab.io.hipo.*
import org.jlab.clas.physics.*
import org.jlab.clas12.physics.*
import extended_kinematic_fitters.*
import analyzers.*
import groovy.io.FileType
import clasqa.QADB


// -----------------------------------------------------------------------------
// Shared helpers
// -----------------------------------------------------------------------------

static int recIndexForPidOccurrence(HipoDataBank recBank, int pid, int occurrence) {
    int found = 0
    for (int row = 0; row < recBank.rows(); row++) {
        if (recBank.getInt("pid", row) != pid) continue
        if (found == occurrence) return row
        found++
    }
    return -1
}

static int occurrenceForMcRow(HipoDataBank mcBank, int pid, int mcRow) {
    if (mcBank == null || mcRow < 0 || mcRow >= mcBank.rows()) return -1
    int occurrence = 0
    for (int row = 0; row < mcBank.rows(); row++) {
        if (mcBank.getInt("pid", row) != pid) continue
        if (row == mcRow) return occurrence
        occurrence++
    }
    return -1
}

static Map truthMatch(
    HipoDataBank recMatchBank,
    HipoDataBank mcBank,
    HipoDataBank lundBank,
    int recPindex
) {
    Map result = [
        matching_pid: 0,
        mcindex: -1,
        mc_occurrence: -1,
        parent_index: 0,
        parent_pid: 0,
        grandparent_index: 0,
        grandparent_pid: 0
    ]

    if (recMatchBank == null || mcBank == null || recPindex < 0) return result

    int mcindex = -1
    for (int row = 0; row < recMatchBank.rows(); row++) {
        if (recMatchBank.getInt("pindex", row) == recPindex) {
            int candidate = recMatchBank.getInt("mcindex", row)
            if (candidate >= 0 && candidate < mcBank.rows()) {
                mcindex = candidate
                break
            }
        }
    }

    if (mcindex < 0) return result

    int truthPid = mcBank.getInt("pid", mcindex)
    result.matching_pid = truthPid
    result.mcindex = mcindex
    result.mc_occurrence = occurrenceForMcRow(mcBank, truthPid, mcindex)

    // MC::Lund parent is one-based in the existing processing convention.
    // This mirrors processing_mc_two_particles.groovy.
    if (lundBank != null && mcindex < lundBank.rows()) {
        int parentIndex = lundBank.getInt("parent", mcindex)
        result.parent_index = parentIndex

        if (parentIndex > 0) {
            int parentRow = parentIndex - 1
            if (parentRow >= 0 && parentRow < lundBank.rows()) {
                result.parent_pid = lundBank.getInt("pid", parentRow)

                int grandparentIndex = lundBank.getInt("parent", parentRow)
                result.grandparent_index = grandparentIndex

                if (grandparentIndex > 0) {
                    int grandparentRow = grandparentIndex - 1
                    if (grandparentRow >= 0 && grandparentRow < lundBank.rows()) {
                        result.grandparent_pid = lundBank.getInt("pid", grandparentRow)
                    }
                }
            }
        }
    }

    return result
}

static void appendRow(
    File file,
    StringBuilder batchLines,
    List values
) {
    batchLines.append(values.collect { value ->
        if (value == null) return "nan"
        if (value instanceof Double && Double.isNaN((Double)value)) return "nan"
        return value.toString()
    }.join(" "))
    batchLines.append("\n")
}

static String schemaLine(List<Map> schema) {
    return "#SCHEMA " + schema.collect { field ->
        field.name + ":" + field.type
    }.join(" ") + "\n"
}

static List<Double> nanList(int n) {
    return (0..<n).collect { Double.NaN }
}

static List<Double> threePhysicsValues(ThreeParticles v) {
    return [
        v.e_p(),
        v.e_theta(),
        v.e_phi(),
        v.vz_e(),
        v.p1_p(),
        v.p1_theta(),
        v.p1_phi(),
        v.vz_p1(),
        v.p2_p(),
        v.p2_theta(),
        v.p2_phi(),
        v.vz_p2(),
        v.open_angle_ep,
        v.open_angle_ep1,
        v.open_angle_ep2,
        v.open_angle_p1p2,
        v.Q2(),
        v.W(),
        v.Mx2(),
        v.Mx2_1(),
        v.Mx2_2(),
        v.x(),
        v.t(),
        v.t1(),
        v.t2(),
        v.tmin(),
        v.y(),
        v.z(),
        v.z1(),
        v.z2(),
        v.Mh(),
        v.xF(),
        v.xF1(),
        v.xF2(),
        v.pT(),
        v.pT1(),
        v.pT2(),
        v.pTpT(),
        v.xi(),
        v.xi1(),
        v.xi2(),
        v.eta(),
        v.eta1(),
        v.eta2(),
        v.Delta_eta(),
        v.eta1_gN(),
        v.eta2_gN(),
        v.phi1(),
        v.phi2(),
        v.Delta_phi(),
        v.phi(),
        v.phiR(),
        v.theta(),
        v.Depolarization_A(),
        v.Depolarization_B(),
        v.Depolarization_C(),
        v.Depolarization_V(),
        v.Depolarization_W(),
        v.Emiss2(),
        v.theta_gamma_gamma(),
        v.pTmiss()
    ]
}

static List<Double> fourPhysicsValues(FourParticles v) {
    return [
        v.e_p(),
        v.e_theta(),
        v.e_phi(),
        v.vz_e(),
        v.p1_p(),
        v.p1_theta(),
        v.p1_phi(),
        v.vz_p1(),
        v.p2_p(),
        v.p2_theta(),
        v.p2_phi(),
        v.vz_p2(),
        v.p3_p(),
        v.p3_theta(),
        v.p3_phi(),
        v.vz_p3(),
        v.open_angle_ep,
        v.open_angle_ep1,
        v.open_angle_ep2,
        v.open_angle_ep3,
        v.open_angle_p1p2,
        v.open_angle_p1p3,
        v.open_angle_p2p3,
        v.Q2(),
        v.W(),
        v.Mx2(),
        v.Mx2_1(),
        v.Mx2_2(),
        v.Mx2_3(),
        v.Mx2_12(),
        v.Mx2_13(),
        v.Mx2_23(),
        v.x(),
        v.t(),
        v.t1(),
        v.t2(),
        v.t3(),
        v.t12(),
        v.t13(),
        v.t23(),
        v.tmin(),
        v.y(),
        v.z(),
        v.z1(),
        v.z2(),
        v.z3(),
        v.z12(),
        v.z13(),
        v.z23(),
        v.Mh(),
        v.Mh12(),
        v.Mh13(),
        v.Mh23(),
        v.xF(),
        v.xF1(),
        v.xF2(),
        v.xF3(),
        v.xF12(),
        v.xF13(),
        v.xF23(),
        v.pT(),
        v.pT1(),
        v.pT2(),
        v.pT3(),
        v.pT12(),
        v.pT13(),
        v.pT23(),
        v.pTpT(),
        v.xi(),
        v.xi1(),
        v.xi2(),
        v.xi3(),
        v.xi12(),
        v.xi13(),
        v.xi23(),
        v.eta(),
        v.eta1(),
        v.eta2(),
        v.eta3(),
        v.eta12(),
        v.eta13(),
        v.eta23(),
        v.phi1(),
        v.phi2(),
        v.phi3(),
        v.phi12(),
        v.phi13(),
        v.phi23(),
        v.Delta_phi12(),
        v.Delta_phi13(),
        v.Delta_phi23(),
        v.phih(),
        v.phiR(),
        v.theta(),
        v.Depolarization_A(),
        v.Depolarization_B(),
        v.Depolarization_C(),
        v.Depolarization_V(),
        v.Depolarization_W()
    ]
}


static List<Double> singlePhotonTopologyValues(ThreeParticles v) {
    return [
        v.Mx2(), v.Mx2_1(), v.Mx2_2(),
        v.Emiss2(), v.theta_gamma_gamma(), v.pTmiss()
    ]
}

static boolean keepDataEvent(
    QADB qa,
    int runnum,
    int evnum,
    int qadbOverride
) {
    if (qadbOverride == 1) return true
    if (!qa.pass(runnum, evnum)) return false
    if (runnum == 5247) return false
    if (runnum == 5345) return false
    if (runnum > 17768 && runnum <= 17811) return false
    if ([17331,16987,17079,17190,17639].contains(runnum)) return false
    if ([16850,16851,16852,16855,16879].contains(runnum)) return false
    return true
}

static void main(String[] args) {
    long startTime = System.currentTimeMillis()

    if (args.length < 1) {
        println("Usage: processing_epgammagamma.groovy <hipo_dir> <output.txt> [n_files] [beam_energy] [run_override] [qadb_override] [max_output_rows]")
        System.exit(1)
    }

    File inputDirectory = args[0] as File

    if (!inputDirectory.exists() || !inputDirectory.isDirectory()) {
        println("ERROR: input path is not a directory: ${inputDirectory}")
        System.exit(1)
    }

    // Deliberately NON-RECURSIVE: process only .hipo files directly inside
    // the exact directory supplied by the user. Do not crawl subdirectories.
    List<File> hipoList = []
    inputDirectory.eachFile(FileType.FILES) {
        if (it.name.endsWith(".hipo")) hipoList << it
    }
    hipoList.sort { a, b -> a.absolutePath <=> b.absolutePath }

    if (hipoList.isEmpty()) {
        println("ERROR: no .hipo files found directly in: ${inputDirectory}")
        println("Subdirectories are intentionally not searched.")
        System.exit(1)
    }

    println("Found ${hipoList.size()} .hipo file(s) directly in ${inputDirectory}")

    String outputFile = args.length < 2 ? "epgammagamma_dummy_out.txt" : args[1]
    int requestedFiles = args.length < 3 ? 0 : Integer.parseInt(args[2])
    int nFiles = (requestedFiles <= 0 || requestedFiles > hipoList.size()) ? hipoList.size() : requestedFiles
    double beamEnergy = (args.length < 4 || args[3].trim().isEmpty()) ? 10.6 : Double.parseDouble(args[3])

    // A run override <= 0 means "use RUN::config", which makes processing.csh
    // convenient for both data and MC.
    Integer userProvidedRun = null
    if (args.length >= 5 && !args[4].trim().isEmpty()) {
        int parsedRun = Integer.parseInt(args[4])
        if (parsedRun > 0) userProvidedRun = parsedRun
    }

    int qadbOverride = (args.length < 6 || args[5].trim().isEmpty()) ? 0 : Integer.parseInt(args[5])
    long maxOutputRows = (args.length < 7 || args[6].trim().isEmpty()) ? 0L : Long.parseLong(args[6])

    List<Map> schema = [
        [name:"is_mc", type:"I"],
        [name:"fiducial_status", type:"I"],
        [name:"num_pos", type:"I"],
        [name:"num_neg", type:"I"],
        [name:"num_neutral", type:"I"],
        [name:"runnum", type:"I"],
        [name:"evnum", type:"I"],
        [name:"helicity", type:"I"],
        [name:"detector1", type:"I"],
        [name:"detector_gamma1", type:"I"],
        [name:"detector_gamma2", type:"I"],
        [name:"gamma1_index", type:"I"],
        [name:"gamma2_index", type:"I"],
        [name:"rec_e_pindex", type:"I"],
        [name:"rec_p_pindex", type:"I"],
        [name:"rec_gamma1_pindex", type:"I"],
        [name:"rec_gamma2_pindex", type:"I"],
        [name:"e_p", type:"D"],
        [name:"e_theta", type:"D"],
        [name:"e_phi", type:"D"],
        [name:"vz_e", type:"D"],
        [name:"p1_p", type:"D"],
        [name:"p1_theta", type:"D"],
        [name:"p1_phi", type:"D"],
        [name:"vz_p1", type:"D"],
        [name:"p2_p", type:"D"],
        [name:"p2_theta", type:"D"],
        [name:"p2_phi", type:"D"],
        [name:"vz_p2", type:"D"],
        [name:"p3_p", type:"D"],
        [name:"p3_theta", type:"D"],
        [name:"p3_phi", type:"D"],
        [name:"vz_p3", type:"D"],
        [name:"open_angle_ep", type:"D"],
        [name:"open_angle_ep1", type:"D"],
        [name:"open_angle_ep2", type:"D"],
        [name:"open_angle_ep3", type:"D"],
        [name:"open_angle_p1p2", type:"D"],
        [name:"open_angle_p1p3", type:"D"],
        [name:"open_angle_p2p3", type:"D"],
        [name:"Q2", type:"D"],
        [name:"W", type:"D"],
        [name:"Mx2", type:"D"],
        [name:"Mx2_1", type:"D"],
        [name:"Mx2_2", type:"D"],
        [name:"Mx2_3", type:"D"],
        [name:"Mx2_12", type:"D"],
        [name:"Mx2_13", type:"D"],
        [name:"Mx2_23", type:"D"],
        [name:"x", type:"D"],
        [name:"t", type:"D"],
        [name:"t1", type:"D"],
        [name:"t2", type:"D"],
        [name:"t3", type:"D"],
        [name:"t12", type:"D"],
        [name:"t13", type:"D"],
        [name:"t23", type:"D"],
        [name:"tmin", type:"D"],
        [name:"y", type:"D"],
        [name:"z", type:"D"],
        [name:"z1", type:"D"],
        [name:"z2", type:"D"],
        [name:"z3", type:"D"],
        [name:"z12", type:"D"],
        [name:"z13", type:"D"],
        [name:"z23", type:"D"],
        [name:"Mh", type:"D"],
        [name:"Mh12", type:"D"],
        [name:"Mh13", type:"D"],
        [name:"Mh23", type:"D"],
        [name:"xF", type:"D"],
        [name:"xF1", type:"D"],
        [name:"xF2", type:"D"],
        [name:"xF3", type:"D"],
        [name:"xF12", type:"D"],
        [name:"xF13", type:"D"],
        [name:"xF23", type:"D"],
        [name:"pT", type:"D"],
        [name:"pT1", type:"D"],
        [name:"pT2", type:"D"],
        [name:"pT3", type:"D"],
        [name:"pT12", type:"D"],
        [name:"pT13", type:"D"],
        [name:"pT23", type:"D"],
        [name:"pTpT", type:"D"],
        [name:"xi", type:"D"],
        [name:"xi1", type:"D"],
        [name:"xi2", type:"D"],
        [name:"xi3", type:"D"],
        [name:"xi12", type:"D"],
        [name:"xi13", type:"D"],
        [name:"xi23", type:"D"],
        [name:"eta", type:"D"],
        [name:"eta1", type:"D"],
        [name:"eta2", type:"D"],
        [name:"eta3", type:"D"],
        [name:"eta12", type:"D"],
        [name:"eta13", type:"D"],
        [name:"eta23", type:"D"],
        [name:"phi1", type:"D"],
        [name:"phi2", type:"D"],
        [name:"phi3", type:"D"],
        [name:"phi12", type:"D"],
        [name:"phi13", type:"D"],
        [name:"phi23", type:"D"],
        [name:"Delta_phi12", type:"D"],
        [name:"Delta_phi13", type:"D"],
        [name:"Delta_phi23", type:"D"],
        [name:"phih", type:"D"],
        [name:"phiR", type:"D"],
        [name:"theta", type:"D"],
        [name:"DepA", type:"D"],
        [name:"DepB", type:"D"],
        [name:"DepC", type:"D"],
        [name:"DepV", type:"D"],
        [name:"DepW", type:"D"],
        [name:"Mh_gammagamma", type:"D"],
        [name:"gamma1_epgamma_Mx2", type:"D"],
        [name:"gamma1_ep_Mx2", type:"D"],
        [name:"gamma1_egamma_Mx2", type:"D"],
        [name:"gamma1_Emiss2", type:"D"],
        [name:"gamma1_theta_gamma_gamma", type:"D"],
        [name:"gamma1_pTmiss", type:"D"],
        [name:"gamma2_epgamma_Mx2", type:"D"],
        [name:"gamma2_ep_Mx2", type:"D"],
        [name:"gamma2_egamma_Mx2", type:"D"],
        [name:"gamma2_Emiss2", type:"D"],
        [name:"gamma2_theta_gamma_gamma", type:"D"],
        [name:"gamma2_pTmiss", type:"D"],
        [name:"matching_e_pid", type:"I"],
        [name:"e_mcindex", type:"I"],
        [name:"matching_p_pid", type:"I"],
        [name:"p_mcindex", type:"I"],
        [name:"matching_gamma1_pid", type:"I"],
        [name:"gamma1_mcindex", type:"I"],
        [name:"gamma1_parent_index", type:"I"],
        [name:"gamma1_parent_pid", type:"I"],
        [name:"gamma1_grandparent_index", type:"I"],
        [name:"gamma1_grandparent_pid", type:"I"],
        [name:"matching_gamma2_pid", type:"I"],
        [name:"gamma2_mcindex", type:"I"],
        [name:"gamma2_parent_index", type:"I"],
        [name:"gamma2_parent_pid", type:"I"],
        [name:"gamma2_grandparent_index", type:"I"],
        [name:"gamma2_grandparent_pid", type:"I"],
        [name:"gen_valid", type:"I"],
        [name:"gen_e_p", type:"D"],
        [name:"gen_e_theta", type:"D"],
        [name:"gen_e_phi", type:"D"],
        [name:"gen_vz_e", type:"D"],
        [name:"gen_p1_p", type:"D"],
        [name:"gen_p1_theta", type:"D"],
        [name:"gen_p1_phi", type:"D"],
        [name:"gen_vz_p1", type:"D"],
        [name:"gen_p2_p", type:"D"],
        [name:"gen_p2_theta", type:"D"],
        [name:"gen_p2_phi", type:"D"],
        [name:"gen_vz_p2", type:"D"],
        [name:"gen_p3_p", type:"D"],
        [name:"gen_p3_theta", type:"D"],
        [name:"gen_p3_phi", type:"D"],
        [name:"gen_vz_p3", type:"D"],
        [name:"gen_open_angle_ep", type:"D"],
        [name:"gen_open_angle_ep1", type:"D"],
        [name:"gen_open_angle_ep2", type:"D"],
        [name:"gen_open_angle_ep3", type:"D"],
        [name:"gen_open_angle_p1p2", type:"D"],
        [name:"gen_open_angle_p1p3", type:"D"],
        [name:"gen_open_angle_p2p3", type:"D"],
        [name:"gen_Q2", type:"D"],
        [name:"gen_W", type:"D"],
        [name:"gen_Mx2", type:"D"],
        [name:"gen_Mx2_1", type:"D"],
        [name:"gen_Mx2_2", type:"D"],
        [name:"gen_Mx2_3", type:"D"],
        [name:"gen_Mx2_12", type:"D"],
        [name:"gen_Mx2_13", type:"D"],
        [name:"gen_Mx2_23", type:"D"],
        [name:"gen_x", type:"D"],
        [name:"gen_t", type:"D"],
        [name:"gen_t1", type:"D"],
        [name:"gen_t2", type:"D"],
        [name:"gen_t3", type:"D"],
        [name:"gen_t12", type:"D"],
        [name:"gen_t13", type:"D"],
        [name:"gen_t23", type:"D"],
        [name:"gen_tmin", type:"D"],
        [name:"gen_y", type:"D"],
        [name:"gen_z", type:"D"],
        [name:"gen_z1", type:"D"],
        [name:"gen_z2", type:"D"],
        [name:"gen_z3", type:"D"],
        [name:"gen_z12", type:"D"],
        [name:"gen_z13", type:"D"],
        [name:"gen_z23", type:"D"],
        [name:"gen_Mh", type:"D"],
        [name:"gen_Mh12", type:"D"],
        [name:"gen_Mh13", type:"D"],
        [name:"gen_Mh23", type:"D"],
        [name:"gen_xF", type:"D"],
        [name:"gen_xF1", type:"D"],
        [name:"gen_xF2", type:"D"],
        [name:"gen_xF3", type:"D"],
        [name:"gen_xF12", type:"D"],
        [name:"gen_xF13", type:"D"],
        [name:"gen_xF23", type:"D"],
        [name:"gen_pT", type:"D"],
        [name:"gen_pT1", type:"D"],
        [name:"gen_pT2", type:"D"],
        [name:"gen_pT3", type:"D"],
        [name:"gen_pT12", type:"D"],
        [name:"gen_pT13", type:"D"],
        [name:"gen_pT23", type:"D"],
        [name:"gen_pTpT", type:"D"],
        [name:"gen_xi", type:"D"],
        [name:"gen_xi1", type:"D"],
        [name:"gen_xi2", type:"D"],
        [name:"gen_xi3", type:"D"],
        [name:"gen_xi12", type:"D"],
        [name:"gen_xi13", type:"D"],
        [name:"gen_xi23", type:"D"],
        [name:"gen_eta", type:"D"],
        [name:"gen_eta1", type:"D"],
        [name:"gen_eta2", type:"D"],
        [name:"gen_eta3", type:"D"],
        [name:"gen_eta12", type:"D"],
        [name:"gen_eta13", type:"D"],
        [name:"gen_eta23", type:"D"],
        [name:"gen_phi1", type:"D"],
        [name:"gen_phi2", type:"D"],
        [name:"gen_phi3", type:"D"],
        [name:"gen_phi12", type:"D"],
        [name:"gen_phi13", type:"D"],
        [name:"gen_phi23", type:"D"],
        [name:"gen_Delta_phi12", type:"D"],
        [name:"gen_Delta_phi13", type:"D"],
        [name:"gen_Delta_phi23", type:"D"],
        [name:"gen_phih", type:"D"],
        [name:"gen_phiR", type:"D"],
        [name:"gen_theta", type:"D"],
        [name:"gen_DepA", type:"D"],
        [name:"gen_DepB", type:"D"],
        [name:"gen_DepC", type:"D"],
        [name:"gen_DepV", type:"D"],
        [name:"gen_DepW", type:"D"],
        [name:"gen_Mh_gammagamma", type:"D"],
        [name:"gen_gamma1_epgamma_Mx2", type:"D"],
        [name:"gen_gamma1_ep_Mx2", type:"D"],
        [name:"gen_gamma1_egamma_Mx2", type:"D"],
        [name:"gen_gamma1_Emiss2", type:"D"],
        [name:"gen_gamma1_theta_gamma_gamma", type:"D"],
        [name:"gen_gamma1_pTmiss", type:"D"],
        [name:"gen_gamma2_epgamma_Mx2", type:"D"],
        [name:"gen_gamma2_ep_Mx2", type:"D"],
        [name:"gen_gamma2_egamma_Mx2", type:"D"],
        [name:"gen_gamma2_Emiss2", type:"D"],
        [name:"gen_gamma2_theta_gamma_gamma", type:"D"],
        [name:"gen_gamma2_pTmiss", type:"D"]
    ]

    File file = new File(outputFile)
    file.delete()
    file.append(schemaLine(schema))

    GenericKinematicFitter recoFitter = new analysis_fitter(10.6041)
    GenericKinematicFitter mcFitter = new monte_carlo_fitter(10.6041)
    EventFilter recoFilter = new EventFilter("11:2212:22:22:Xn")

    QADB qa = new QADB("latest")
    qa.checkForDefect("TotalOutlier")
    qa.checkForDefect("TerminalOutlier")
    qa.checkForDefect("MarginalOutlier")
    qa.checkForDefect("SectorLoss")
    qa.checkForDefect("Misc")
    qa.checkForDefect("ChargeHigh")
    qa.checkForDefect("ChargeNegative")
    qa.checkForDefect("ChargeUnknown")
    qa.checkForDefect("PossiblyNoBeam")

    [6736,6737,6738,6739,6740,6741,6742,6743,6744,6746,6747,6748,6749,
     6750,6751,6753,6754,6755,6756,6757,16194,16089,16185,16308,16184,
     16307,16309,16872,16975,17763,17764,17765,17766,17767,17768,17179,
     17180,17181,17182,17183,17188,17189,17252].each { run -> qa.allowMiscBit(run) }

    StringBuilder batchLines = new StringBuilder()
    int bufferedRows = 0
    long inputEvents = 0L
    long outputRows = 0L
    boolean stopRequested = false

    for (int currentFile = 0; currentFile < nFiles && !stopRequested; currentFile++) {
        println("\nOpening file ${currentFile+1} out of ${nFiles}: ${hipoList[currentFile]}")

        HipoDataSource reader = new HipoDataSource()
        reader.open(hipoList[currentFile])

        while (reader.hasEvent()) {
            def event = reader.getNextEvent()
            inputEvents++

            if (inputEvents % 500000L == 0L) {
                println("processed HIPO events: ${inputEvents}, output epgammagamma rows: ${outputRows}")
            }

            if (!event.hasBank("REC::Particle") || !event.hasBank("RUN::config")) continue

            boolean isMC = event.hasBank("MC::Particle")
            HipoDataBank runBank = (HipoDataBank)event.getBank("RUN::config")
            int runnum = userProvidedRun != null ? userProvidedRun : runBank.getInt("run", 0)
            int evnum = runBank.getInt("event", 0)

            PhysicsEvent recoEvent = recoFitter.getPhysicsEvent(event)
            if (!recoFilter.isValid(recoEvent)) continue
            if (!isMC && !keepDataEvent(qa, runnum, evnum, qadbOverride)) continue

            double energy
            if (isMC || runnum == 11) {
                energy = beamEnergy
            } else {
                BeamEnergy Eb = new BeamEnergy(recoEvent, runnum, false)
                energy = Eb.Eb()
            }

            HipoDataBank recBank = (HipoDataBank)event.getBank("REC::Particle")
            HipoDataBank recMatchBank = event.hasBank("MC::RecMatch") ? (HipoDataBank)event.getBank("MC::RecMatch") : null
            HipoDataBank mcBank = isMC ? (HipoDataBank)event.getBank("MC::Particle") : null
            HipoDataBank lundBank = event.hasBank("MC::Lund") ? (HipoDataBank)event.getBank("MC::Lund") : null
            PhysicsEvent mcEvent = isMC ? mcFitter.getPhysicsEvent(event) : null

            int numGammas = recoEvent.countByPid(22)
            if (numGammas < 2) continue

            int recEIndex = recIndexForPidOccurrence(recBank, 11, 0)
            int recPIndex = recIndexForPidOccurrence(recBank, 2212, 0)
            Map eTruth = isMC ? truthMatch(recMatchBank, mcBank, lundBank, recEIndex) : truthMatch(null,null,null,-1)
            Map pTruth = isMC ? truthMatch(recMatchBank, mcBank, lundBank, recPIndex) : truthMatch(null,null,null,-1)

            for (int gamma1 = 0; gamma1 < numGammas-1; gamma1++) {
                for (int gamma2 = gamma1+1; gamma2 < numGammas; gamma2++) {
                    FourParticles recoVars
                    ThreeParticles recoSingle1
                    ThreeParticles recoSingle2

                    try {
                        recoVars = new FourParticles(event, recoEvent, 2212, 0, 22, gamma1, 22, gamma2, energy)
                        recoSingle1 = new ThreeParticles(event, recoEvent, 2212, 0, 22, gamma1, energy)
                        recoSingle2 = new ThreeParticles(event, recoEvent, 2212, 0, 22, gamma2, energy)
                    } catch (Exception ex) {
                        continue
                    }

                    if (!recoVars.channel_test(recoVars)) continue

                    int recGamma1Index = recIndexForPidOccurrence(recBank, 22, gamma1)
                    int recGamma2Index = recIndexForPidOccurrence(recBank, 22, gamma2)

                    Map gamma1Truth = isMC ? truthMatch(recMatchBank, mcBank, lundBank, recGamma1Index) : truthMatch(null,null,null,-1)
                    Map gamma2Truth = isMC ? truthMatch(recMatchBank, mcBank, lundBank, recGamma2Index) : truthMatch(null,null,null,-1)

                    int genValid = 0
                    List<Double> genFourValues = nanList(99)
                    List<Double> genSingle1Values = nanList(6)
                    List<Double> genSingle2Values = nanList(6)
                    double genMgg = Double.NaN

                    if (isMC && mcEvent != null &&
                        eTruth.matching_pid == 11 &&
                        pTruth.matching_pid == 2212 &&
                        gamma1Truth.matching_pid == 22 &&
                        gamma2Truth.matching_pid == 22 &&
                        pTruth.mc_occurrence >= 0 &&
                        gamma1Truth.mc_occurrence >= 0 &&
                        gamma2Truth.mc_occurrence >= 0 &&
                        gamma1Truth.mc_occurrence != gamma2Truth.mc_occurrence) {
                        try {
                            FourParticles genVars = new FourParticles(
                                event, mcEvent,
                                2212, (int)pTruth.mc_occurrence,
                                22, (int)gamma1Truth.mc_occurrence,
                                22, (int)gamma2Truth.mc_occurrence,
                                energy
                            )
                            ThreeParticles genSingle1 = new ThreeParticles(
                                event, mcEvent,
                                2212, (int)pTruth.mc_occurrence,
                                22, (int)gamma1Truth.mc_occurrence,
                                energy
                            )
                            ThreeParticles genSingle2 = new ThreeParticles(
                                event, mcEvent,
                                2212, (int)pTruth.mc_occurrence,
                                22, (int)gamma2Truth.mc_occurrence,
                                energy
                            )

                            genFourValues = fourPhysicsValues(genVars)
                            genSingle1Values = singlePhotonTopologyValues(genSingle1)
                            genSingle2Values = singlePhotonTopologyValues(genSingle2)
                            genMgg = genVars.Mh23()
                            genValid = 1
                        } catch (Exception ex) {
                            genValid = 0
                        }
                    }

                    List values = [
                        isMC ? 1 : 0,
                        recoVars.get_fiducial_status(),
                        recoVars.get_num_pos(),
                        recoVars.get_num_neg(),
                        recoVars.get_num_neutrals(),
                        runnum,
                        evnum,
                        recoVars.get_helicity(),
                        recoVars.get_detector1(),
                        recoVars.get_detector2(),
                        recoVars.get_detector3(),
                        gamma1,
                        gamma2,
                        recEIndex,
                        recPIndex,
                        recGamma1Index,
                        recGamma2Index
                    ]

                    values.addAll(fourPhysicsValues(recoVars))
                    values.add(recoVars.Mh23())
                    values.addAll(singlePhotonTopologyValues(recoSingle1))
                    values.addAll(singlePhotonTopologyValues(recoSingle2))

                    values.addAll([
                        eTruth.matching_pid, eTruth.mcindex,
                        pTruth.matching_pid, pTruth.mcindex,
                        gamma1Truth.matching_pid, gamma1Truth.mcindex,
                        gamma1Truth.parent_index, gamma1Truth.parent_pid,
                        gamma1Truth.grandparent_index, gamma1Truth.grandparent_pid,
                        gamma2Truth.matching_pid, gamma2Truth.mcindex,
                        gamma2Truth.parent_index, gamma2Truth.parent_pid,
                        gamma2Truth.grandparent_index, gamma2Truth.grandparent_pid,
                        genValid
                    ])

                    values.addAll(genFourValues)
                    values.add(genMgg)
                    values.addAll(genSingle1Values)
                    values.addAll(genSingle2Values)

                    appendRow(file, batchLines, values)
                    bufferedRows++
                    outputRows++

                    if (bufferedRows >= 1000) {
                        file.append(batchLines.toString())
                        batchLines.setLength(0)
                        bufferedRows = 0
                    }

                    if (maxOutputRows > 0L && outputRows >= maxOutputRows) {
                        stopRequested = true
                        break
                    }
                }
                if (stopRequested) break
            }
        }

        reader.close()
    }

    if (batchLines.length() > 0) file.append(batchLines.toString())

    println("output file: ${file}")
    println("input HIPO events processed: ${inputEvents}")
    println("output epgammagamma candidate rows: ${outputRows}")
    println("No M(gamma gamma) cut is applied; Mh_gammagamma is written for every unique pair.")
    println("Elapsed time: " + (System.currentTimeMillis()-startTime) + " ms")
}
