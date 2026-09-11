/*
 * Dedicated photon-efficiency tag-and-probe skim.
 *
 * One output row is one e' p' gamma_tag hypothesis.  The expensive HIPO pass is
 * intentionally permissive: the final pi0 selection and probe matching are done
 * later from the ROOT trees.
 *
 * Arguments:
 *   0 input .hipo file (or directory)
 *   1 output text file
 *   2 MC beam energy (used for run 11; default 10.6041 GeV)
 *   3 run override (0 = use RUN::config)
 *   4 QADB override (0 = use QADB, 1 = skip QADB)
 *   5 is_mc (0/1)
 *   6 loose Mx2(ep) minimum (default -1.0 GeV^2)
 *   7 loose Mx2(ep) maximum (default  2.0 GeV^2)
 *
 * The tag photon is only required to be REC::Particle pid==22, p>=0.4 GeV,
 * and have FD/FT status.  beta/fiducial decisions are SAVED, not imposed.
 */

import org.jlab.io.hipo.*
import org.jlab.io.base.DataEvent
import org.jlab.clas.physics.*
import org.jlab.clas12.physics.*
import extended_kinematic_fitters.*
import analyzers.*
import groovy.io.FileType
import clasqa.QADB
import groovy.transform.Field

@Field static final double ME = 0.0005109989461
@Field static final double MP = 0.9382720813
@Field static final double DEG = 180.0 / Math.PI
@Field static final double SENT = -999.0
@Field static final int ISENT = -999
@Field static final int N_NEUTRAL_SAVE = 5
@Field static final int N_ANY_SAVE = 3

static double clamp(double x, double lo, double hi) { Math.max(lo, Math.min(hi, x)) }
static double p3(double x, double y, double z) { Math.sqrt(x*x + y*y + z*z) }
static double thetaDeg(double x, double y, double z) {
    double p = p3(x,y,z)
    return p > 0 ? Math.acos(clamp(z/p, -1.0, 1.0))*DEG : SENT
}
static double phiDeg(double x, double y) { Math.atan2(y,x)*DEG }
static double wrapPhi(double x) {
    while (x >= 180.0) x -= 360.0
    while (x < -180.0) x += 360.0
    return x
}
static double openingDeg(double ax, double ay, double az, double bx, double by, double bz) {
    double am = p3(ax,ay,az), bm = p3(bx,by,bz)
    if (am <= 0 || bm <= 0) return SENT
    return Math.acos(clamp((ax*bx + ay*by + az*bz)/(am*bm), -1.0, 1.0))*DEG
}
static int detectorFromStatus(int status) {
    int a = Math.abs(status)
    if (a >= 1000 && a < 2000) return 0 // FT
    if (a >= 2000 && a < 4000) return 1 // FD
    if (a >= 4000 && a < 5000) return 2 // CD
    return -1
}
static double nominalBeamEnergy(int runnum, double mcBeam) {
    if (runnum == 11) return mcBeam
    if (runnum >= 5032 && runnum <= 5666) return 10.6041
    if (runnum >= 3172 && runnum <= 3817) return 10.5940
    if (runnum >= 3863 && runnum <= 4326) return 10.5940
    if (runnum >= 6616 && runnum <= 6783) return 10.1998
    if (runnum >= 6120 && runnum <= 6399) return 10.5986
    if (runnum >= 6409 && runnum <= 6604) return 10.1998
    return mcBeam
}

static boolean looseProtonTest(int particleIndex, float vz, double electronVz,
                               HipoDataBank recBank, HipoDataBank calBank,
                               HipoDataBank trajBank, HipoDataBank runBank,
                               generic_tests gtests, fiducial_cuts fcuts) {
    double px=recBank.getFloat("px",particleIndex), py=recBank.getFloat("py",particleIndex), pz=recBank.getFloat("pz",particleIndex)
    double p=p3(px,py,pz)
    boolean fd=gtests.forward_detector_cut(particleIndex,recBank)
    boolean cd=gtests.central_detector_cut(particleIndex,recBank)
    double torus=runBank.getFloat("torus",0)
    return true
        && (cd ? p>0.3 : true)
        && (fd && torus>0 ? p>0.42 : true)
        && (fd && torus<0 ? p>0.50 : true)
        // Deliberately DO NOT impose analysis_fitter's p<1.2 GeV speed cut here.
        // That can be restored offline with p_pass_standard or p_corr_p.
        && thetaDeg(px,py,pz)<64.23
        && gtests.vertex_cut(particleIndex,recBank,runBank)
        && (fd ? fcuts.dc_fiducial_cut(particleIndex,recBank,trajBank,runBank) : true)
        && (cd ? fcuts.cvt_fiducial_cut(particleIndex,recBank,trajBank,2) : true)
}

static double mcWeight(HipoDataEvent event) {
    if (event.hasBank("MC::Event")) {
        HipoDataBank b = (HipoDataBank)event.getBank("MC::Event")
        if (b.rows() > 0) return b.getFloat("weight",0)
    }
    return 1.0
}

// Return detector response coordinates for a REC::Particle pindex.
// Prefer FT for FT-status objects, otherwise PCAL layer 1.  If no layer-1
// ECAL row exists, use the first associated calorimeter row.
static double[] responseXYZ(int pindex, int detector, HipoDataBank calBank, HipoDataBank ftBank) {
    if (detector == 0 && ftBank != null) {
        for (int r=0; r<ftBank.rows(); r++) {
            if (ftBank.getInt("pindex",r) == pindex) {
                return [ftBank.getFloat("x",r), ftBank.getFloat("y",r), ftBank.getFloat("z",r)] as double[]
            }
        }
    }
    if (calBank != null) {
        int fallback = -1
        for (int r=0; r<calBank.rows(); r++) {
            if (calBank.getInt("pindex",r) != pindex) continue
            if (fallback < 0) fallback = r
            if (calBank.getInt("layer",r) == 1) {
                return [calBank.getFloat("x",r), calBank.getFloat("y",r), calBank.getFloat("z",r)] as double[]
            }
        }
        if (fallback >= 0) {
            return [calBank.getFloat("x",fallback), calBank.getFloat("y",fallback), calBank.getFloat("z",fallback)] as double[]
        }
    }
    return [SENT,SENT,SENT] as double[]
}

static double[] candidateDirection(int i, HipoDataBank recBank, HipoDataBank calBank, HipoDataBank ftBank,
                                   double vxRef, double vyRef, double vzRef) {
    double px=recBank.getFloat("px",i), py=recBank.getFloat("py",i), pz=recBank.getFloat("pz",i)
    double pmag=p3(px,py,pz)
    if (pmag > 1.0e-6) return [px,py,pz] as double[]
    int det = detectorFromStatus(recBank.getInt("status",i))
    double[] xyz = responseXYZ(i,det,calBank,ftBank)
    if (xyz[0] < -900) return [0.0,0.0,0.0] as double[]
    return [xyz[0]-vxRef, xyz[1]-vyRef, xyz[2]-vzRef] as double[]
}

static List<Map> nearestCandidates(HipoDataBank recBank, HipoDataBank calBank, HipoDataBank ftBank,
                                   double predx, double predy, double predz,
                                   double vxRef, double vyRef, double vzRef,
                                   int electronIndex, int protonIndex, int tagIndex,
                                   boolean neutralOnly, int nsave) {
    List<Map> out=[]
    for (int i=0; i<recBank.rows(); i++) {
        if (i==electronIndex || i==protonIndex || i==tagIndex) continue
        int charge = recBank.getByte("charge",i)
        if (neutralOnly && charge != 0) continue
        double[] dir = candidateDirection(i,recBank,calBank,ftBank,vxRef,vyRef,vzRef)
        double da = openingDeg(predx,predy,predz,dir[0],dir[1],dir[2])
        if (da < -900) continue
        double pmag=p3(recBank.getFloat("px",i),recBank.getFloat("py",i),recBank.getFloat("pz",i))
        int status=recBank.getInt("status",i), det=detectorFromStatus(status)
        double[] xyz=responseXYZ(i,det,calBank,ftBank)
        out << [idx:i, pid:recBank.getInt("pid",i), charge:charge, status:status, det:det,
                p:pmag, theta:thetaDeg(dir[0],dir[1],dir[2]), phi:phiDeg(dir[0],dir[1]), da:da,
                x:xyz[0],y:xyz[1],z:xyz[2]]
    }
    out.sort { a,b -> a.da <=> b.da }
    if (out.size()>nsave) return out[0..<nsave]
    return out
}

static Map truthMatch(HipoDataEvent event, double dx, double dy, double dz, boolean requirePhoton, boolean requirePi0Parent) {
    Map best=[pid:ISENT,parent:ISENT,p:SENT,theta:SENT,phi:SENT,da:SENT,index:ISENT]
    if (!event.hasBank("MC::Lund")) return best
    HipoDataBank lund=(HipoDataBank)event.getBank("MC::Lund")
    double bestDa=1.0e9
    for (int i=0;i<lund.rows();i++) {
        int pid=lund.getInt("pid",i)
        if (requirePhoton && pid!=22) continue
        int parentPid=0
        int parent=lund.getInt("parent",i)
        if (parent>0 && parent-1<lund.rows()) parentPid=lund.getInt("pid",parent-1)
        if (requirePi0Parent && parentPid!=111) continue
        double px=lund.getFloat("px",i), py=lund.getFloat("py",i), pz=lund.getFloat("pz",i)
        double da=openingDeg(dx,dy,dz,px,py,pz)
        if (da>=0 && da<bestDa) {
            bestDa=da
            best=[pid:pid,parent:parentPid,p:p3(px,py,pz),theta:thetaDeg(px,py,pz),phi:phiDeg(px,py),da:da,index:i]
        }
    }
    return best
}

static String fmt(Object x) {
    if (x instanceof Integer || x instanceof Long || x instanceof Short || x instanceof Byte) return x.toString()
    return String.format(java.util.Locale.US,"%.10g",((Number)x).doubleValue())
}

static void appendCandidate(List vals, Map c) {
    if (c == null) {
        vals.add(ISENT); vals.add(ISENT); vals.add(ISENT); vals.add(ISENT); vals.add(ISENT)
        vals.add(SENT); vals.add(SENT); vals.add(SENT); vals.add(SENT); vals.add(SENT); vals.add(SENT); vals.add(SENT)
        return
    }
    vals.add(c.idx); vals.add(c.pid); vals.add(c.charge); vals.add(c.status); vals.add(c.det)
    vals.add(c.p); vals.add(c.theta); vals.add(c.phi); vals.add(c.da); vals.add(c.x); vals.add(c.y); vals.add(c.z)
}

static void main(String[] args) {
    if (args.length < 2) {
        println "Usage: process_photon_efficiency.groovy <input.hipo|dir> <output.txt> [mc_beam] [run_override] [qadb_override] [is_mc] [mx2_min] [mx2_max]"
        System.exit(1)
    }
    long start=System.currentTimeMillis()
    File input=new File(args[0])
    String output=args[1]
    double mcBeam=args.length>2 ? Double.parseDouble(args[2]) : 10.6041
    int runOverride=args.length>3 ? Integer.parseInt(args[3]) : 0
    int qaOverride=args.length>4 ? Integer.parseInt(args[4]) : 0
    int isMC=args.length>5 ? Integer.parseInt(args[5]) : 0
    double mx2Min=args.length>6 ? Double.parseDouble(args[6]) : -1.0
    double mx2Max=args.length>7 ? Double.parseDouble(args[7]) : 2.0

    List<File> hipos=[]
    if (input.isFile() && input.name.endsWith('.hipo')) hipos << input
    else if (input.isDirectory()) input.eachFileRecurse(FileType.FILES) { if (it.name.endsWith('.hipo')) hipos << it }
    hipos.sort { it.absolutePath }
    if (hipos.isEmpty()) { println "ERROR: no HIPO files found"; System.exit(2) }

    analysis_fitter fitter=new analysis_fitter(mcBeam)
    fiducial_cuts fcuts=new fiducial_cuts()
    pid_cuts pcuts=new pid_cuts()
    generic_tests gtests=new generic_tests()
    energy_loss_corrections eloss=new energy_loss_corrections()

    QADB qa=new QADB("latest")
    qa.checkForDefect('TotalOutlier'); qa.checkForDefect('TerminalOutlier'); qa.checkForDefect('MarginalOutlier')
    qa.checkForDefect('SectorLoss'); qa.checkForDefect('Misc'); qa.checkForDefect('ChargeHigh')
    qa.checkForDefect('ChargeNegative'); qa.checkForDefect('ChargeUnknown'); qa.checkForDefect('PossiblyNoBeam')
    [6736,6737,6738,6739,6740,6741,6742,6743,6744,6746,6747,6748,6749,6750,6751,6753,6754,6755,6756,6757].each { qa.allowMiscBit(it) }

    File outf=new File(output); outf.parentFile?.mkdirs(); outf.delete()
    BufferedWriter writer=new BufferedWriter(new FileWriter(outf))
    StringBuilder batch=new StringBuilder(); int lineCount=0
    long nevt=0, nrow=0, npair=0

    for (File hf : hipos) {
        println "Opening ${hf.absolutePath}"
        HipoDataSource reader=new HipoDataSource(); reader.open(hf)
        while (reader.hasEvent()) {
            HipoDataEvent event=(HipoDataEvent)reader.getNextEvent(); nevt++
            if (!event.hasBank("RUN::config") || !event.hasBank("REC::Particle") ||
                !event.hasBank("REC::Calorimeter") || !event.hasBank("REC::Traj") || !event.hasBank("REC::Cherenkov")) continue
            HipoDataBank run=(HipoDataBank)event.getBank("RUN::config")
            HipoDataBank rec=(HipoDataBank)event.getBank("REC::Particle")
            HipoDataBank cal=(HipoDataBank)event.getBank("REC::Calorimeter")
            HipoDataBank traj=(HipoDataBank)event.getBank("REC::Traj")
            HipoDataBank cc=(HipoDataBank)event.getBank("REC::Cherenkov")
            HipoDataBank ft=event.hasBank("REC::ForwardTagger") ? (HipoDataBank)event.getBank("REC::ForwardTagger") : null
            HipoDataBank evb=event.hasBank("REC::Event") ? (HipoDataBank)event.getBank("REC::Event") : null
            if (rec.rows()<1 || rec.getInt("pid",0)!=11) continue

            int runnum=runOverride!=0 ? runOverride : run.getInt("run",0)
            int evnum=run.getInt("event",0)
            if (!(runnum==11 || qaOverride==1 || qa.pass(runnum,evnum))) continue
            if (runnum==5247 || runnum==5345) continue

            double beamE=nominalBeamEnergy(runnum,mcBeam)
            double epx=rec.getFloat("px",0), epy=rec.getFloat("py",0), epz=rec.getFloat("pz",0)
            double ep=p3(epx,epy,epz), evz=rec.getFloat("vz",0)
            if (!fitter.electron_test(0,ep,rec,cal,traj,run,cc)) continue

            LorentzVector beam=new LorentzVector(); beam.setPxPyPzM(0,0,Math.sqrt(Math.max(0,beamE*beamE-ME*ME)),ME)
            LorentzVector target=new LorentzVector(); target.setPxPyPzM(0,0,0,MP)
            LorentzVector ele=new LorentzVector(); ele.setPxPyPzM(epx,epy,epz,ME)
            LorentzVector q=new LorentzVector(beam); q.sub(ele)
            double Q2=-q.mass2(), nu=q.e()
            double xB=nu>0 ? Q2/(2*MP*nu) : SENT
            double yy=beamE>0 ? nu/beamE : SENT
            LorentzVector hadsys=new LorentzVector(target); hadsys.add(q)
            double W2=hadsys.mass2(), W=W2>=0?Math.sqrt(W2):-Math.sqrt(-W2)

            int helicity=evb!=null && evb.rows()>0 ? evb.getByte("helicity",0) : 0
            float torus=run.getFloat("torus",0), solenoid=run.getFloat("solenoid",0)

            // broad multiplicity diagnostics from REC::Particle
            int nNeutral=0,nPid22=0,nPid2112=0,nPid0Neutral=0
            for (int ii=0;ii<rec.rows();ii++) {
                int charge=rec.getByte("charge",ii), pid=rec.getInt("pid",ii)
                if (charge==0) { nNeutral++; if(pid==22)nPid22++; if(pid==2112)nPid2112++; if(pid==0)nPid0Neutral++; }
            }

            for (int ip=0; ip<rec.rows(); ip++) {
                if (rec.getInt("pid",ip)!=2212) continue
                float pprawx=rec.getFloat("px",ip), pprawy=rec.getFloat("py",ip), pprawz=rec.getFloat("pz",ip)
                float pvz=rec.getFloat("vz",ip)
                boolean pPassStandard=fitter.proton_test(ip,2212,pvz,evz,rec,cal,traj,run)
                if (!looseProtonTest(ip,pvz,evz,rec,cal,traj,run,gtests,fcuts)) continue
                float[] pcorr=[pprawx,pprawy,pprawz] as float[]
                eloss.proton_energy_loss_corrections(ip,pcorr,rec,run)
                LorentzVector protRaw=new LorentzVector(); protRaw.setPxPyPzM(pprawx,pprawy,pprawz,MP)
                LorentzVector prot=new LorentzVector(); prot.setPxPyPzM(pcorr[0],pcorr[1],pcorr[2],MP)
                LorentzVector missEP=new LorentzVector(beam); missEP.add(target); missEP.sub(ele); missEP.sub(prot)
                double mx2ep=missEP.mass2()
                if (mx2ep<mx2Min || mx2ep>mx2Max) continue
                npair++

                LorentzVector tp=new LorentzVector(target); tp.sub(prot)
                double minusT=-tp.mass2()

                for (int ig=0; ig<rec.rows(); ig++) {
                    if (ig==0 || ig==ip) continue
                    if (rec.getInt("pid",ig)!=22) continue
                    double gpx=rec.getFloat("px",ig), gpy=rec.getFloat("py",ig), gpz=rec.getFloat("pz",ig), gp=p3(gpx,gpy,gpz)
                    if (gp<0.4) continue
                    int gstatus=rec.getInt("status",ig), gdet=detectorFromStatus(gstatus)
                    if (!(gdet==0 || gdet==1)) continue

                    boolean tagBeta=pcuts.beta_cut(ig,rec)
                    boolean tagFid=(gdet==1) ? fcuts.pcal_fiducial_cut(ig,2,run,rec,cal)
                                                   : (ft!=null && fcuts.forward_tagger_fiducial_cut(ig,rec,ft))
                    double[] tagXYZ=responseXYZ(ig,gdet,cal,ft)

                    LorentzVector tagRaw=new LorentzVector(); tagRaw.setPxPyPzM(gpx,gpy,gpz,0)
                    float[] gcorr=[(float)gpx,(float)gpy,(float)gpz] as float[]
                    eloss.sebastian_photon_energy_loss_corrections(ig,gcorr,rec,run)
                    LorentzVector tagCorr=new LorentzVector(); tagCorr.setPxPyPzM(gcorr[0],gcorr[1],gcorr[2],0)

                    LorentzVector probeRaw=new LorentzVector(beam); probeRaw.add(target); probeRaw.sub(ele); probeRaw.sub(prot); probeRaw.sub(tagRaw)
                    LorentzVector probeCorr=new LorentzVector(beam); probeCorr.add(target); probeCorr.sub(ele); probeCorr.sub(prot); probeCorr.sub(tagCorr)
                    double predx=probeRaw.px(), predy=probeRaw.py(), predz=probeRaw.pz()

                    List<Map> neutrals=nearestCandidates(rec,cal,ft,predx,predy,predz,0,0,evz,0,ip,ig,true,N_NEUTRAL_SAVE)
                    List<Map> any=nearestCandidates(rec,cal,ft,predx,predy,predz,0,0,evz,0,ip,ig,false,N_ANY_SAVE)

                    Map truthProbePi0=truthMatch(event,predx,predy,predz,true,true)
                    Map truthProbeAny=truthProbePi0.pid!=ISENT ? truthProbePi0 : truthMatch(event,predx,predy,predz,true,false)
                    Map truthTag=truthMatch(event,gpx,gpy,gpz,false,false)

                    List vals=[]
                    vals.add(runnum); vals.add(evnum); vals.add(helicity); vals.add(isMC)
                    vals.add(beamE); vals.add(torus); vals.add(solenoid); vals.add(mcWeight(event))
                    vals.add(rec.rows()); vals.add(nNeutral); vals.add(nPid22); vals.add(nPid2112); vals.add(nPid0Neutral)
                    vals.add(Q2); vals.add(W); vals.add(xB); vals.add(yy); vals.add(minusT)

                    vals.add(ep); vals.add(thetaDeg(epx,epy,epz)); vals.add(phiDeg(epx,epy)); vals.add(rec.getFloat("vx",0)); vals.add(rec.getFloat("vy",0)); vals.add(evz)

                    vals.add(ip); vals.add(detectorFromStatus(rec.getInt("status",ip))); vals.add(pPassStandard?1:0)
                    vals.add(p3(pprawx,pprawy,pprawz)); vals.add(thetaDeg(pprawx,pprawy,pprawz)); vals.add(phiDeg(pprawx,pprawy));
                    vals.add(p3(pcorr[0],pcorr[1],pcorr[2])); vals.add(thetaDeg(pcorr[0],pcorr[1],pcorr[2])); vals.add(phiDeg(pcorr[0],pcorr[1])); vals.add(pvz)

                    vals.add(ig); vals.add(gdet); vals.add(gstatus); vals.add(tagBeta?1:0); vals.add(tagFid?1:0)
                    vals.add(gp); vals.add(thetaDeg(gpx,gpy,gpz)); vals.add(phiDeg(gpx,gpy)); vals.add(rec.getFloat("beta",ig))
                    vals.add(p3(gcorr[0],gcorr[1],gcorr[2])); vals.add(thetaDeg(gcorr[0],gcorr[1],gcorr[2])); vals.add(phiDeg(gcorr[0],gcorr[1]));
                    vals.add(tagXYZ[0]); vals.add(tagXYZ[1]); vals.add(tagXYZ[2]); vals.add(tagXYZ[0]>-900?Math.hypot(tagXYZ[0],tagXYZ[1]):SENT)

                    vals.add(mx2ep); vals.add(mx2ep>=0?Math.sqrt(mx2ep):-Math.sqrt(-mx2ep))
                    vals.add(probeRaw.mass2()); vals.add(probeRaw.e()); vals.add(probeRaw.p()); vals.add(thetaDeg(predx,predy,predz)); vals.add(phiDeg(predx,predy))
                    vals.add(probeRaw.px()); vals.add(probeRaw.py()); vals.add(probeRaw.pz())
                    vals.add(probeCorr.mass2()); vals.add(probeCorr.e()); vals.add(probeCorr.p()); vals.add(thetaDeg(probeCorr.px(),probeCorr.py(),probeCorr.pz())); vals.add(phiDeg(probeCorr.px(),probeCorr.py()))
                    vals.add(probeCorr.px()); vals.add(probeCorr.py()); vals.add(probeCorr.pz())

                    for (int k=0;k<N_NEUTRAL_SAVE;k++) appendCandidate(vals,k<neutrals.size()?neutrals[k]:null)
                    for (int k=0;k<N_ANY_SAVE;k++) appendCandidate(vals,k<any.size()?any[k]:null)

                    vals.add(truthProbeAny.index); vals.add(truthProbeAny.pid); vals.add(truthProbeAny.parent); vals.add(truthProbeAny.p); vals.add(truthProbeAny.theta); vals.add(truthProbeAny.phi); vals.add(truthProbeAny.da)
                    vals.add(truthTag.index); vals.add(truthTag.pid); vals.add(truthTag.parent); vals.add(truthTag.p); vals.add(truthTag.theta); vals.add(truthTag.phi); vals.add(truthTag.da)

                    batch.append(vals.collect{fmt(it)}.join(' ')).append('\n'); lineCount++; nrow++
                    if (lineCount>=1000) { writer.write(batch.toString()); batch.setLength(0); lineCount=0 }
                } // tag
            } // proton
            if (nevt%1000000==0) println "processed ${nevt} events, wrote ${nrow} tag hypotheses"
        }
        reader.close()
    }
    if (batch.length()>0) writer.write(batch.toString())
    writer.close()
    double min=(System.currentTimeMillis()-start)/60000.0
    println String.format(java.util.Locale.US,"Done: events=%d accepted_ep_pairs=%d rows=%d time=%.2f min output=%s",nevt,npair,nrow,min,output)
}

// Groovy scripts executed with run-groovy do not automatically invoke a
// user-defined static main(String[]).  Explicitly enter the processor here.
main(args)
