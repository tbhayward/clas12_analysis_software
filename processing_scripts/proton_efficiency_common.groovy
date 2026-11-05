/*
 * Shared event-aligned processor for the ep -> e pi+ pi- (p) proton
 * tag-and-probe study.  The data and CLASDIS entry points below call this
 * implementation with different modes but intentionally produce one schema.
 *
 * One tag row is written for every selected e pi+ pi- combination.  Every
 * Event Builder proton in that event is written to the proton sidecar.  No
 * proton matching window or proton-angle restriction is imposed here.
 */

import org.jlab.io.hipo.*
import org.jlab.clas.physics.*
import extended_kinematic_fitters.*
import analyzers.*
import groovy.io.FileType
import clasqa.QADB

class ProtonEfficiencyProcessor {
    static final double ELECTRON_MASS = 0.0005109989461d
    static final double PROTON_MASS = 0.9382720813d

    static double momentum(double x, double y, double z) {
        return Math.sqrt(x*x + y*y + z*z)
    }

    static double theta(double x, double y, double z) {
        double p = momentum(x, y, z)
        if (!(p > 0.0d)) return Double.NaN
        return Math.toDegrees(Math.acos(Math.max(-1.0d, Math.min(1.0d, z/p))))
    }

    static double phi(double x, double y) {
        double value = Math.toDegrees(Math.atan2(y, x))
        return value < 0.0d ? value + 360.0d : value
    }

    static double openingAngle(double ax, double ay, double az,
                               double bx, double by, double bz) {
        double am = momentum(ax, ay, az), bm = momentum(bx, by, bz)
        if (!(am > 0.0d) || !(bm > 0.0d)) return Double.NaN
        double cosine = (ax*bx + ay*by + az*bz)/(am*bm)
        return Math.toDegrees(Math.acos(Math.max(-1.0d, Math.min(1.0d, cosine))))
    }

    static int detectorFromStatus(int status) {
        int value = Math.abs(status)
        if (value >= 2000 && value < 4000) return 1 // FD
        if (value >= 4000 && value < 5000) return 2 // CD
        if (value >= 1000 && value < 2000) return 0 // FT
        return -1
    }

    static LorentzVector particleVector(Particle particle, double mass) {
        LorentzVector vector = new LorentzVector()
        vector.setPxPyPzM(particle.px(), particle.py(), particle.pz(), mass)
        return vector
    }

    static LorentzVector missingProton(PhysicsEvent physicsEvent, int pipIndex,
                                       int pimIndex, double beamEnergy) {
        Particle electron = physicsEvent.getParticleByPid(11, 0)
        Particle pip = physicsEvent.getParticleByPid(211, pipIndex)
        Particle pim = physicsEvent.getParticleByPid(-211, pimIndex)
        LorentzVector missing = new LorentzVector()
        double beamP = Math.sqrt(beamEnergy*beamEnergy - ELECTRON_MASS*ELECTRON_MASS)
        missing.setPxPyPzM(0.0d, 0.0d, beamP, ELECTRON_MASS)
        LorentzVector target = new LorentzVector()
        target.setPxPyPzM(0.0d, 0.0d, 0.0d, PROTON_MASS)
        missing.add(target)
        missing.sub(particleVector(electron, ELECTRON_MASS))
        missing.sub(particleVector(pip, 0.13957039d))
        missing.sub(particleVector(pim, 0.13957039d))
        return missing
    }

    static List<Map> protonCandidates(HipoDataBank recBank, HipoDataBank calBank,
                                      HipoDataBank trajBank, HipoDataBank runBank,
                                      analysis_fitter fitter, double electronVz,
                                      LorentzVector missing) {
        generic_tests tests = new generic_tests()
        fiducial_cuts fiducials = new fiducial_cuts()
        List<Map> candidates = []
        for (int row = 0; row < recBank.rows(); row++) {
            if (recBank.getInt('pid', row) != 2212) continue
            double px = recBank.getFloat('px', row)
            double py = recBank.getFloat('py', row)
            double pz = recBank.getFloat('pz', row)
            double p = momentum(px, py, pz)
            int status = recBank.getInt('status', row)
            boolean fd = tests.forward_detector_cut(row, recBank)
            boolean cd = tests.central_detector_cut(row, recBank)
            boolean vertexPass = tests.vertex_cut(row, recBank, runBank)
            boolean fiducialPass = (fd && fiducials.dc_fiducial_cut(row, recBank, trajBank, runBank)) ||
                                   (cd && fiducials.cvt_fiducial_cut(row, recBank, trajBank, 2))
            boolean momentumPass = p < 1.2d
            boolean selected = fitter.proton_test(row, 2212, recBank.getFloat('vz', row),
                                                  electronVz, recBank, calBank,
                                                  trajBank, runBank)
            candidates << [
                rec_index:row,
                selected:selected ? 1 : 0,
                vertex_pass:vertexPass ? 1 : 0,
                fiducial_pass:fiducialPass ? 1 : 0,
                momentum_pass:momentumPass ? 1 : 0,
                status:status,
                detector:detectorFromStatus(status),
                charge:recBank.getInt('charge', row),
                px:px, py:py, pz:pz, p:p,
                theta:theta(px, py, pz), phi:phi(px, py),
                vz:recBank.getFloat('vz', row),
                beta:recBank.getFloat('beta', row),
                chi2pid:recBank.getFloat('chi2pid', row),
                match_angle:openingAngle(missing.px(), missing.py(), missing.pz(), px, py, pz),
                delta_p:p - missing.p(),
                relative_delta_p:(missing.p() > 0.0d ? (p-missing.p())/missing.p() : Double.NaN)
            ]
        }
        return candidates
    }

    static Map closest(List<Map> candidates) {
        List<Map> finite = candidates.findAll { Double.isFinite((double)it.match_angle) }
        return finite ? finite.min { (double)it.match_angle } : null
    }

    static Map truthSummary(HipoDataEvent event, LorentzVector missing) {
        Map result = [n_e:0, n_pip:0, n_pim:0, n_proton:0, best:null]
        if (!event.hasBank('MC::Particle')) return result
        HipoDataBank bank = (HipoDataBank)event.getBank('MC::Particle')
        List<Map> protons = []
        for (int row = 0; row < bank.rows(); row++) {
            int pid = bank.getInt('pid', row)
            if (pid == 11) result.n_e++
            if (pid == 211) result.n_pip++
            if (pid == -211) result.n_pim++
            if (pid != 2212) continue
            result.n_proton++
            double px = bank.getFloat('px', row), py = bank.getFloat('py', row), pz = bank.getFloat('pz', row)
            protons << [index:row, px:px, py:py, pz:pz, p:momentum(px,py,pz),
                        theta:theta(px,py,pz), phi:phi(px,py),
                        match_angle:openingAngle(missing.px(), missing.py(), missing.pz(), px, py, pz)]
        }
        result.best = closest(protons)
        return result
    }

    static int matchedTruthIndex(HipoDataEvent event, int recIndex) {
        if (recIndex < 0 || !event.hasBank('MC::RecMatch')) return -1
        HipoDataBank matches = (HipoDataBank)event.getBank('MC::RecMatch')
        for (int row = 0; row < matches.rows(); row++) {
            if (matches.getInt('pindex', row) == recIndex) return matches.getInt('mcindex', row)
        }
        return -1
    }

    static Map matchedTruthTag(HipoDataEvent event, int eRec, int pipRec, int pimRec,
                               double beamEnergy, LorentzVector reconstructedMissing) {
        Map invalid = [valid:0, e_index:-1, pip_index:-1, pim_index:-1,
            e_p:Double.NaN,e_theta:Double.NaN,e_phi:Double.NaN,
            pip_p:Double.NaN,pip_theta:Double.NaN,pip_phi:Double.NaN,
            pim_p:Double.NaN,pim_theta:Double.NaN,pim_phi:Double.NaN,
            Q2:Double.NaN,W:Double.NaN,xB:Double.NaN,y:Double.NaN,
            missing:reconstructedMissing]
        if (!event.hasBank('MC::Particle')) return invalid
        HipoDataBank bank = (HipoDataBank)event.getBank('MC::Particle')
        int eIndex = matchedTruthIndex(event, eRec)
        int pipIndex = matchedTruthIndex(event, pipRec)
        int pimIndex = matchedTruthIndex(event, pimRec)
        if (eIndex < 0 || pipIndex < 0 || pimIndex < 0 ||
            eIndex >= bank.rows() || pipIndex >= bank.rows() || pimIndex >= bank.rows()) return invalid
        if (bank.getInt('pid',eIndex) != 11 || bank.getInt('pid',pipIndex) != 211 ||
            bank.getInt('pid',pimIndex) != -211) return invalid
        double ex=bank.getFloat('px',eIndex), ey=bank.getFloat('py',eIndex), ez=bank.getFloat('pz',eIndex)
        double ppx=bank.getFloat('px',pipIndex), ppy=bank.getFloat('py',pipIndex), ppz=bank.getFloat('pz',pipIndex)
        double pmx=bank.getFloat('px',pimIndex), pmy=bank.getFloat('py',pimIndex), pmz=bank.getFloat('pz',pimIndex)
        LorentzVector beam = new LorentzVector()
        beam.setPxPyPzM(0.0d,0.0d,Math.sqrt(beamEnergy*beamEnergy-ELECTRON_MASS*ELECTRON_MASS),ELECTRON_MASS)
        LorentzVector target = new LorentzVector()
        target.setPxPyPzM(0.0d,0.0d,0.0d,PROTON_MASS)
        LorentzVector electron = new LorentzVector(); electron.setPxPyPzM(ex,ey,ez,ELECTRON_MASS)
        LorentzVector pip = new LorentzVector(); pip.setPxPyPzM(ppx,ppy,ppz,0.13957039d)
        LorentzVector pim = new LorentzVector(); pim.setPxPyPzM(pmx,pmy,pmz,0.13957039d)
        LorentzVector q = new LorentzVector(beam); q.sub(electron)
        LorentzVector hadronic = new LorentzVector(q); hadronic.add(target)
        LorentzVector missing = new LorentzVector(hadronic); missing.sub(pip); missing.sub(pim)
        double Q2 = -q.mass2(), nu = q.e()
        return [valid:1,e_index:eIndex,pip_index:pipIndex,pim_index:pimIndex,
            e_p:momentum(ex,ey,ez),e_theta:theta(ex,ey,ez),e_phi:phi(ex,ey),
            pip_p:momentum(ppx,ppy,ppz),pip_theta:theta(ppx,ppy,ppz),pip_phi:phi(ppx,ppy),
            pim_p:momentum(pmx,pmy,pmz),pim_theta:theta(pmx,pmy,pmz),pim_phi:phi(pmx,pmy),
            Q2:Q2,W:hadronic.mass(),xB:(nu != 0.0d ? Q2/(2.0d*PROTON_MASS*nu) : Double.NaN),
            y:nu/beamEnergy,missing:missing]
    }

    static List truthTagValues(Map tag) {
        if ((int)tag.valid == 0) {
            return [0,tag.e_index,tag.pip_index,tag.pim_index] +
                   (0..<22).collect { Double.NaN }
        }
        LorentzVector missing = (LorentzVector)tag.missing
        return [tag.valid,tag.e_index,tag.pip_index,tag.pim_index,
            tag.e_p,tag.e_theta,tag.e_phi,(tag.e_p > 2.0d ? 1 : 0),tag.pip_p,tag.pip_theta,tag.pip_phi,
            tag.pim_p,tag.pim_theta,tag.pim_phi,tag.Q2,tag.W,tag.xB,tag.y,
            missing.px(),missing.py(),missing.pz(),missing.e(),missing.mass2(),missing.p(),
            theta(missing.px(),missing.py(),missing.pz()),phi(missing.px(),missing.py())]
    }

    static List bestValues(Map best) {
        if (best == null) return [-1, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                                  Double.NaN, Double.NaN, Double.NaN, Double.NaN]
        return [best.rec_index, best.px, best.py, best.pz, best.p, best.theta,
                best.phi, best.match_angle, best.relative_delta_p]
    }

    static List truthValues(Map truth) {
        Map best = (Map)truth.best
        if (best == null) return [truth.n_e, truth.n_pip, truth.n_pim, truth.n_proton,
                                  -1, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
                                  Double.NaN, Double.NaN, Double.NaN]
        return [truth.n_e, truth.n_pip, truth.n_pim, truth.n_proton, best.index,
                best.px, best.py, best.pz, best.p, best.theta, best.phi, best.match_angle]
    }

    static void writeRow(BufferedWriter writer, List values) {
        writer.write(values.collect { it == null ? 'NaN' : it.toString() }.join(' '))
        writer.newLine()
    }

    static QADB makeQa() {
        QADB qa = new QADB('latest')
        ['TotalOutlier','TerminalOutlier','MarginalOutlier','SectorLoss','LowLiveTime',
         'Misc','ChargeHigh','ChargeNegative','ChargeUnknown'].each { qa.checkForDefect(it) }
        [6736,6737,6738,6739,6740,6741,6742,6743,6744,6746,6748,6749,
         6750,6751,6753,6754,6755,6756,6757].each { qa.allowMiscBit(it) }
        return qa
    }

    static void run(String[] args, boolean mcMode) {
        if (args.length < 2) {
            System.err.println('Usage: processor input_hipo_path output_prefix [n_files=0] [beam_energy=10.6] [run_override=0] [qadb_override=0] [max_events=0]')
            System.exit(2)
        }
        File input = new File(args[0])
        String outputPrefix = args[1]
        int requestedFiles = args.length > 2 ? Integer.parseInt(args[2]) : 0
        double configuredBeam = args.length > 3 ? Double.parseDouble(args[3]) : 10.6d
        int runOverride = args.length > 4 ? Integer.parseInt(args[4]) : 0
        boolean qadbOverride = args.length > 5 && Integer.parseInt(args[5]) != 0
        long maxEvents = args.length > 6 ? Long.parseLong(args[6]) : 0L
        if (maxEvents < 0L) throw new IllegalArgumentException('max_events must be nonnegative')

        List<File> files = []
        if (input.isFile() && input.name.endsWith('.hipo')) files << input
        else input.eachFileRecurse(FileType.FILES) { if (it.name.endsWith('.hipo')) files << it }
        files.sort { a, b -> a.absolutePath <=> b.absolutePath }
        int nFiles = requestedFiles <= 0 ? files.size() : Math.min(requestedFiles, files.size())
        if (nFiles == 0) throw new IllegalArgumentException('No HIPO files found under ' + input)

        File tagFile = new File(outputPrefix + '.tags.tsv')
        File protonFile = new File(outputPrefix + '.protons.tsv')
        tagFile.parentFile?.mkdirs()
        protonFile.parentFile?.mkdirs()
        BufferedWriter tags = new BufferedWriter(new FileWriter(tagFile))
        BufferedWriter protons = new BufferedWriter(new FileWriter(protonFile))

        List tagHeader = ['file_index','event_sequence','tag_sequence','run','event','torus','helicity','beam_energy',
            'e_rec_index','pip_rec_index','pim_rec_index',
            'e_p','e_theta','e_phi','rec_e_p_gt_2','e_vz','pip_p','pip_theta','pip_phi','pip_vz','pim_p','pim_theta','pim_phi','pim_vz',
            'Q2','W','xB','y','missing_px','missing_py','missing_pz','missing_E','missing_M2','missing_p','missing_theta','missing_phi',
            'n_proton_eb','n_proton_selected',
            'best_eb_rec_index','best_eb_px','best_eb_py','best_eb_pz','best_eb_p','best_eb_theta','best_eb_phi','best_eb_match_angle','best_eb_relative_delta_p',
            'best_selected_rec_index','best_selected_px','best_selected_py','best_selected_pz','best_selected_p','best_selected_theta','best_selected_phi','best_selected_match_angle','best_selected_relative_delta_p',
            'truth_tag_match_valid','truth_e_index','truth_pip_index','truth_pim_index',
            'truth_e_p','truth_e_theta','truth_e_phi','truth_e_p_gt_2','truth_pip_p','truth_pip_theta','truth_pip_phi','truth_pim_p','truth_pim_theta','truth_pim_phi',
            'truth_Q2','truth_W','truth_xB','truth_y','truth_missing_px','truth_missing_py','truth_missing_pz','truth_missing_E','truth_missing_M2','truth_missing_p','truth_missing_theta','truth_missing_phi',
            'truth_n_e','truth_n_pip','truth_n_pim','truth_n_proton','truth_best_index','truth_best_px','truth_best_py','truth_best_pz','truth_best_p','truth_best_theta','truth_best_phi','truth_best_match_angle']
        List protonHeader = ['file_index','event_sequence','tag_sequence','run','event','rec_index','selected','vertex_pass','fiducial_pass','momentum_pass',
            'status','detector','charge','px','py','pz','p','theta','phi','vz','beta','chi2pid','match_angle','delta_p','relative_delta_p']
        tags.write(tagHeader.join(' ') + '\n')
        protons.write(protonHeader.join(' ') + '\n')

        analysis_fitter fitter = new analysis_fitter(configuredBeam)
        EventFilter tagFilter = new EventFilter('11:211:-211:X+:X-:Xn')
        QADB qa = mcMode ? null : makeQa()
        long eventSequence = 0L, tagSequence = 0L, selectedTags = 0L
        boolean stop = false

        for (int fileIndex = 0; fileIndex < nFiles && !stop; fileIndex++) {
            println("Opening ${fileIndex+1}/${nFiles}: ${files[fileIndex]}")
            HipoDataSource reader = new HipoDataSource()
            reader.open(files[fileIndex])
            // The first nSidis event is the intentional scaler/readout event.
            if (reader.hasEvent()) reader.getNextEvent()
            while (reader.hasEvent()) {
                if (maxEvents > 0L && eventSequence >= maxEvents) { stop = true; break }
                HipoDataEvent event = reader.getNextEvent()
                eventSequence++
                if (!event.hasBank('RUN::config') || !event.hasBank('REC::Particle')) continue
                HipoDataBank runBank = (HipoDataBank)event.getBank('RUN::config')
                int physicalRun = runBank.getInt('run', 0)
                int run = runOverride != 0 ? runOverride : physicalRun
                int eventNumber = runBank.getInt('event', 0)
                if (!mcMode && !qadbOverride && !qa.pass(run, eventNumber)) continue
                if (!mcMode && (run == 5247 || run == 5345)) continue

                PhysicsEvent physicsEvent = fitter.getPhysicsEvent(event)
                if (!tagFilter.isValid(physicsEvent)) continue
                HipoDataBank recBank = (HipoDataBank)event.getBank('REC::Particle')
                HipoDataBank calBank = (HipoDataBank)event.getBank('REC::Calorimeter')
                HipoDataBank trajBank = (HipoDataBank)event.getBank('REC::Traj')
                int nPip = physicsEvent.countByPid(211), nPim = physicsEvent.countByPid(-211)
                BeamEnergy beam = new BeamEnergy(physicsEvent, run, false)
                double beamEnergy = mcMode || physicalRun == 11 ? configuredBeam : beam.Eb()

                for (int pipIndex = 0; pipIndex < nPip; pipIndex++) {
                    for (int pimIndex = 0; pimIndex < nPim; pimIndex++) {
                        ThreeParticles k = new ThreeParticles(event, physicsEvent, 211, pipIndex, -211, pimIndex, beamEnergy)
                        if (!ThreeParticles.channel_test(k)) continue
                        tagSequence++
                        selectedTags++
                        LorentzVector missing = missingProton(physicsEvent, pipIndex, pimIndex, beamEnergy)
                        List<Map> candidates = protonCandidates(recBank, calBank, trajBank, runBank, fitter, k.vz_e(), missing)
                        Map bestEb = closest(candidates)
                        Map bestSelected = closest(candidates.findAll { it.selected == 1 })
                        int eRec = fitter.getLastRecIndexForPidOccurrence(11, 0)
                        int pipRec = fitter.getLastRecIndexForPidOccurrence(211, pipIndex)
                        int pimRec = fitter.getLastRecIndexForPidOccurrence(-211, pimIndex)
                        Map truthTag = mcMode ? matchedTruthTag(event,eRec,pipRec,pimRec,beamEnergy,missing) :
                            [valid:0,e_index:-1,pip_index:-1,pim_index:-1,e_p:Double.NaN,e_theta:Double.NaN,e_phi:Double.NaN,
                             pip_p:Double.NaN,pip_theta:Double.NaN,pip_phi:Double.NaN,pim_p:Double.NaN,pim_theta:Double.NaN,pim_phi:Double.NaN,
                             Q2:Double.NaN,W:Double.NaN,xB:Double.NaN,y:Double.NaN,missing:missing]
                        Map truth = mcMode ? truthSummary(event, (LorentzVector)truthTag.missing) : [n_e:0,n_pip:0,n_pim:0,n_proton:0,best:null]
                        int helicity = event.hasBank('REC::Event') ? event.getBank('REC::Event').getByte('helicity',0) : 0
                        List row = [fileIndex,eventSequence,tagSequence,run,eventNumber,runBank.getFloat('torus',0),helicity,beamEnergy,
                            eRec,pipRec,pimRec,k.e_p(),k.e_theta(),k.e_phi(),(k.e_p() > 2.0d ? 1 : 0),k.vz_e(),k.p1_p(),k.p1_theta(),k.p1_phi(),k.vz_p1(),
                            k.p2_p(),k.p2_theta(),k.p2_phi(),k.vz_p2(),k.Q2(),k.W(),k.x(),k.y(),missing.px(),missing.py(),missing.pz(),
                            missing.e(),missing.mass2(),missing.p(),theta(missing.px(),missing.py(),missing.pz()),phi(missing.px(),missing.py()),
                            candidates.size(),candidates.count { it.selected == 1 }]
                        row.addAll(bestValues(bestEb))
                        row.addAll(bestValues(bestSelected))
                        row.addAll(truthTagValues(truthTag))
                        row.addAll(truthValues(truth))
                        writeRow(tags, row)
                        for (Map candidate : candidates) {
                            writeRow(protons, [fileIndex,eventSequence,tagSequence,run,eventNumber,candidate.rec_index,candidate.selected,
                                candidate.vertex_pass,candidate.fiducial_pass,candidate.momentum_pass,candidate.status,candidate.detector,candidate.charge,
                                candidate.px,candidate.py,candidate.pz,candidate.p,candidate.theta,candidate.phi,candidate.vz,candidate.beta,candidate.chi2pid,
                                candidate.match_angle,candidate.delta_p,candidate.relative_delta_p])
                        }
                        if (selectedTags % 10000L == 0L) { tags.flush(); protons.flush() }
                    }
                }
                if (eventSequence % 100000L == 0L) println("processed=${eventSequence}, selected_tags=${selectedTags}")
            }
            reader.close()
        }
        tags.close()
        protons.close()
        println("Complete: events=${eventSequence}, tags=${selectedTags}, tag_file=${tagFile}, proton_file=${protonFile}")
    }
}

