/*
 * Unified e'p'gammaX processing for data and reconstructed MC.
 *
 * - Automatically detects MC from MC::Particle.
 * - Uses reconstructed ThreeParticles kinematics for every photon candidate.
 * - For MC, uses MC::RecMatch + MC::Particle + MC::Lund to save the truth PID,
 *   immediate parent PID, and grandparent PID of the reconstructed photon.
 * - When the reconstructed e/p/gamma all have usable truth matches, a second
 *   ThreeParticles object is built from monte_carlo_fitter and the generated
 *   kinematics are saved with gen_* names.
 * - Optional final argument limits OUTPUT CANDIDATE ROWS. This is used by
 *   processing.csh to make the quick 1M tree before restarting for full stats.
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
        v.phih(),
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
        println("Usage: processing_epgamma.groovy <hipo_dir> <output.txt> [n_files] [beam_energy] [run_override] [qadb_override] [max_output_rows]")
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

    String outputFile = args.length < 2 ? "epgamma_dummy_out.txt" : args[1]
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

    println("Groovy epgamma resolved arguments:")
    println("  raw args count = ${args.length}")
    for (int i = 0; i < args.length; i++) {
        println("  args[${i}] = '${args[i]}'")
    }
    println("  requestedFiles = ${requestedFiles}")
    println("  beamEnergy = ${beamEnergy}")
    println("  runOverride = ${userProvidedRun}")
    println("  qadbOverride = ${qadbOverride}")
    println("  maxOutputRows = ${maxOutputRows}")

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
        [name:"detector2", type:"I"],
        [name:"gamma_index", type:"I"],
        [name:"rec_e_pindex", type:"I"],
        [name:"rec_p_pindex", type:"I"],
        [name:"rec_gamma_pindex", type:"I"],
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
        [name:"open_angle_ep", type:"D"],
        [name:"open_angle_ep1", type:"D"],
        [name:"open_angle_ep2", type:"D"],
        [name:"open_angle_p1p2", type:"D"],
        [name:"Q2", type:"D"],
        [name:"W", type:"D"],
        [name:"Mx2", type:"D"],
        [name:"Mx2_1", type:"D"],
        [name:"Mx2_2", type:"D"],
        [name:"x", type:"D"],
        [name:"t", type:"D"],
        [name:"t1", type:"D"],
        [name:"t2", type:"D"],
        [name:"tmin", type:"D"],
        [name:"y", type:"D"],
        [name:"z", type:"D"],
        [name:"z1", type:"D"],
        [name:"z2", type:"D"],
        [name:"Mh", type:"D"],
        [name:"xF", type:"D"],
        [name:"xF1", type:"D"],
        [name:"xF2", type:"D"],
        [name:"pT", type:"D"],
        [name:"pT1", type:"D"],
        [name:"pT2", type:"D"],
        [name:"pTpT", type:"D"],
        [name:"xi", type:"D"],
        [name:"xi1", type:"D"],
        [name:"xi2", type:"D"],
        [name:"eta", type:"D"],
        [name:"eta1", type:"D"],
        [name:"eta2", type:"D"],
        [name:"Delta_eta", type:"D"],
        [name:"eta1_gN", type:"D"],
        [name:"eta2_gN", type:"D"],
        [name:"phi1", type:"D"],
        [name:"phi2", type:"D"],
        [name:"Delta_phi", type:"D"],
        [name:"phih", type:"D"],
        [name:"phiR", type:"D"],
        [name:"theta", type:"D"],
        [name:"DepA", type:"D"],
        [name:"DepB", type:"D"],
        [name:"DepC", type:"D"],
        [name:"DepV", type:"D"],
        [name:"DepW", type:"D"],
        [name:"Emiss2", type:"D"],
        [name:"theta_gamma_gamma", type:"D"],
        [name:"pTmiss", type:"D"],
        [name:"matching_e_pid", type:"I"],
        [name:"e_mcindex", type:"I"],
        [name:"matching_p_pid", type:"I"],
        [name:"p_mcindex", type:"I"],
        [name:"matching_gamma_pid", type:"I"],
        [name:"gamma_mcindex", type:"I"],
        [name:"gamma_parent_index", type:"I"],
        [name:"gamma_parent_pid", type:"I"],
        [name:"gamma_grandparent_index", type:"I"],
        [name:"gamma_grandparent_pid", type:"I"],
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
        [name:"gen_open_angle_ep", type:"D"],
        [name:"gen_open_angle_ep1", type:"D"],
        [name:"gen_open_angle_ep2", type:"D"],
        [name:"gen_open_angle_p1p2", type:"D"],
        [name:"gen_Q2", type:"D"],
        [name:"gen_W", type:"D"],
        [name:"gen_Mx2", type:"D"],
        [name:"gen_Mx2_1", type:"D"],
        [name:"gen_Mx2_2", type:"D"],
        [name:"gen_x", type:"D"],
        [name:"gen_t", type:"D"],
        [name:"gen_t1", type:"D"],
        [name:"gen_t2", type:"D"],
        [name:"gen_tmin", type:"D"],
        [name:"gen_y", type:"D"],
        [name:"gen_z", type:"D"],
        [name:"gen_z1", type:"D"],
        [name:"gen_z2", type:"D"],
        [name:"gen_Mh", type:"D"],
        [name:"gen_xF", type:"D"],
        [name:"gen_xF1", type:"D"],
        [name:"gen_xF2", type:"D"],
        [name:"gen_pT", type:"D"],
        [name:"gen_pT1", type:"D"],
        [name:"gen_pT2", type:"D"],
        [name:"gen_pTpT", type:"D"],
        [name:"gen_xi", type:"D"],
        [name:"gen_xi1", type:"D"],
        [name:"gen_xi2", type:"D"],
        [name:"gen_eta", type:"D"],
        [name:"gen_eta1", type:"D"],
        [name:"gen_eta2", type:"D"],
        [name:"gen_Delta_eta", type:"D"],
        [name:"gen_eta1_gN", type:"D"],
        [name:"gen_eta2_gN", type:"D"],
        [name:"gen_phi1", type:"D"],
        [name:"gen_phi2", type:"D"],
        [name:"gen_Delta_phi", type:"D"],
        [name:"gen_phih", type:"D"],
        [name:"gen_phiR", type:"D"],
        [name:"gen_theta", type:"D"],
        [name:"gen_DepA", type:"D"],
        [name:"gen_DepB", type:"D"],
        [name:"gen_DepC", type:"D"],
        [name:"gen_DepV", type:"D"],
        [name:"gen_DepW", type:"D"],
        [name:"gen_Emiss2", type:"D"],
        [name:"gen_theta_gamma_gamma", type:"D"],
        [name:"gen_pTmiss", type:"D"]
    ]

    File file = new File(outputFile)
    file.delete()
    file.append(schemaLine(schema))

    // Use the standard reconstructed-event fitter available in
    // processing_classes.jar. The previous dvcs_fitter reference caused a
    // NoClassDefFoundError on the farm.
    GenericKinematicFitter recoFitter = new analysis_fitter(10.6041)
    GenericKinematicFitter mcFitter = new monte_carlo_fitter(10.6041)
    EventFilter recoFilter = new EventFilter("11:2212:22:Xn")

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

    for (
        int currentFile = 0;
        currentFile < nFiles
            && !stopRequested
            && (maxOutputRows <= 0L || outputRows < maxOutputRows);
        currentFile++
    ) {
        println("\nOpening file ${currentFile+1} out of ${nFiles}: ${hipoList[currentFile]}")

        HipoDataSource reader = new HipoDataSource()
        reader.open(hipoList[currentFile])

        while (
            reader.hasEvent()
            && !stopRequested
            && (maxOutputRows <= 0L || outputRows < maxOutputRows)
        ) {
            def event = reader.getNextEvent()
            inputEvents++

            if (inputEvents % 500000L == 0L) {
                println("processed HIPO events: ${inputEvents}, output epgamma rows: ${outputRows}")
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
            if (numGammas < 1) continue

            int recEIndex = recIndexForPidOccurrence(recBank, 11, 0)
            int recPIndex = recIndexForPidOccurrence(recBank, 2212, 0)
            Map eTruth = isMC ? truthMatch(recMatchBank, mcBank, lundBank, recEIndex) : truthMatch(null,null,null,-1)
            Map pTruth = isMC ? truthMatch(recMatchBank, mcBank, lundBank, recPIndex) : truthMatch(null,null,null,-1)

            for (int currentGamma = 0; currentGamma < numGammas; currentGamma++) {
                ThreeParticles recoVars = new ThreeParticles(event, recoEvent, 2212, 0, 22, currentGamma, energy)
                if (!recoVars.channel_test(recoVars)) continue

                int recGammaIndex = recIndexForPidOccurrence(recBank, 22, currentGamma)
                Map gammaTruth = isMC ? truthMatch(recMatchBank, mcBank, lundBank, recGammaIndex) : truthMatch(null,null,null,-1)

                int genValid = 0
                List<Double> genValues = nanList(61)

                if (isMC && mcEvent != null &&
                    eTruth.matching_pid == 11 &&
                    pTruth.matching_pid == 2212 &&
                    gammaTruth.matching_pid == 22 &&
                    pTruth.mc_occurrence >= 0 &&
                    gammaTruth.mc_occurrence >= 0) {
                    try {
                        ThreeParticles genVars = new ThreeParticles(
                            event, mcEvent,
                            2212, (int)pTruth.mc_occurrence,
                            22, (int)gammaTruth.mc_occurrence,
                            energy
                        )
                        genValues = threePhysicsValues(genVars)
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
                    currentGamma,
                    recEIndex,
                    recPIndex,
                    recGammaIndex
                ]

                values.addAll(threePhysicsValues(recoVars))
                values.addAll([
                    eTruth.matching_pid, eTruth.mcindex,
                    pTruth.matching_pid, pTruth.mcindex,
                    gammaTruth.matching_pid, gammaTruth.mcindex,
                    gammaTruth.parent_index, gammaTruth.parent_pid,
                    gammaTruth.grandparent_index, gammaTruth.grandparent_pid,
                    genValid
                ])
                values.addAll(genValues)

                appendRow(file, batchLines, values)
                bufferedRows++
                outputRows++

                if (bufferedRows >= 1000) {
                    file.append(batchLines.toString())
                    batchLines.setLength(0)
                    bufferedRows = 0
                }

                if (maxOutputRows > 0L && outputRows >= maxOutputRows) {
                    println(
                        "Reached requested output-row limit: ${outputRows} / ${maxOutputRows}. " +
                        "Stopping quick pass."
                    )
                    stopRequested = true
                    break
                }
            }
        }

        reader.close()
    }

    if (batchLines.length() > 0) file.append(batchLines.toString())

    println("output file: ${file}")
    println("input HIPO events processed: ${inputEvents}")
    println("output epgamma candidate rows: ${outputRows}")
    println("MC truth is filled only when MC::Particle exists; parent ancestry additionally requires MC::Lund.")
    println("Elapsed time: " + (System.currentTimeMillis()-startTime) + " ms")
}
