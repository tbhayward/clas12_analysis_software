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
 *   8 loose Mx2(ep gamma_tag) minimum (default -0.25 GeV^2)
 *   9 loose Mx2(ep gamma_tag) maximum (default  0.25 GeV^2)
 *  10 sample kind: data|aaogen|clasdis|dvcsgen|auto (default auto)
 *
 * For sample kind clasdis, generator-exclusive e p pi0 events are intentionally
 * suppressed from the skim so the CLASDIS sample represents inclusive pi0 /
 * extra-particle topologies while AAOgen supplies the pure exclusive pi0 sample.
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
import java.util.zip.CRC32

@Field static final double ME = 0.0005109989461
@Field static final double MP = 0.9382720813
@Field static final double DEG = 180.0 / Math.PI
@Field static final double SENT = -999.0
@Field static final int ISENT = -999
@Field static final int N_NEUTRAL_SAVE = 5
@Field static final int N_ANY_SAVE = 3
@Field static final int N_GEN_GAMMA_SAVE = 12
@Field static final int SKIM_VERSION = 3

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

static void insertNearest(List<Map> out, Map cand, int nsave) {
    int pos=0
    while (pos<out.size() && out[pos].da <= cand.da) pos++
    out.add(pos,cand)
    if (out.size()>nsave) out.remove(out.size()-1)
}

// Find the nearest neutral candidates and nearest candidates of any charge in
// ONE REC::Particle scan.  Only the best N are retained as we go, so no large
// temporary candidate lists or full-list sorts are created for each tag.
static Map nearestCandidateSets(HipoDataBank recBank, HipoDataBank calBank, HipoDataBank ftBank,
                                double predx, double predy, double predz,
                                double vxRef, double vyRef, double vzRef,
                                int electronIndex, int protonIndex, int tagIndex) {
    List<Map> neutrals=[]
    List<Map> any=[]
    for (int i=0; i<recBank.rows(); i++) {
        if (i==electronIndex || i==protonIndex || i==tagIndex) continue
        int charge=recBank.getByte("charge",i)
        double[] dir=candidateDirection(i,recBank,calBank,ftBank,vxRef,vyRef,vzRef)
        double da=openingDeg(predx,predy,predz,dir[0],dir[1],dir[2])
        if (da < -900) continue
        double pmag=p3(recBank.getFloat("px",i),recBank.getFloat("py",i),recBank.getFloat("pz",i))
        int status=recBank.getInt("status",i), det=detectorFromStatus(status)
        double[] xyz=responseXYZ(i,det,calBank,ftBank)
        Map rs=recResponseSummary(i,recBank,calBank,ftBank)
        Map cand=[idx:i, pid:recBank.getInt("pid",i), charge:charge, status:status, det:det,
                  p:pmag, theta:thetaDeg(dir[0],dir[1],dir[2]), phi:phiDeg(dir[0],dir[1]), da:da,
                  x:xyz[0],y:xyz[1],z:xyz[2], beta:rs.beta,chi2pid:rs.chi2pid,responseE:rs.responseE,
                  pcalE:rs.pcalE,ecinE:rs.ecinE,ecoutE:rs.ecoutE,ecalE:rs.ecalE,
                  pcalSector:rs.pcalSector,pcalLu:rs.pcalLu,pcalLv:rs.pcalLv,pcalLw:rs.pcalLw,ftE:rs.ftE,ftR:rs.ftR]
        insertNearest(any,cand,N_ANY_SAVE)
        if (charge==0) insertNearest(neutrals,cand,N_NEUTRAL_SAVE)
    }
    return [neutrals:neutrals, any:any]
}

static Map emptyTruth() {
    return [pid:ISENT,parent:ISENT,p:SENT,theta:SENT,phi:SENT,da:SENT,index:ISENT]
}

// Determine both probe and tag truth matches in one MC::Lund scan.  The probe
// prefers a generated photon whose parent is a pi0; if none exists, it falls
// back to the nearest generated photon of any parent.  The tag diagnostic keeps
// the nearest generated particle, matching the previous behavior.
static Map truthMatches(HipoDataEvent event,
                        double predx, double predy, double predz,
                        double tagx, double tagy, double tagz) {
    Map bestProbePi0=emptyTruth(), bestProbeAny=emptyTruth(), bestTag=emptyTruth()
    if (!event.hasBank("MC::Lund")) return [probe:bestProbeAny, tag:bestTag]
    HipoDataBank lund=(HipoDataBank)event.getBank("MC::Lund")
    double bestProbePi0Da=1.0e9, bestProbeAnyDa=1.0e9, bestTagDa=1.0e9
    for (int i=0;i<lund.rows();i++) {
        int pid=lund.getInt("pid",i)
        int parentPid=0
        int parent=lund.getInt("parent",i)
        if (parent>0 && parent-1<lund.rows()) parentPid=lund.getInt("pid",parent-1)
        double px=lund.getFloat("px",i), py=lund.getFloat("py",i), pz=lund.getFloat("pz",i)
        double pp=p3(px,py,pz)
        double th=thetaDeg(px,py,pz), ph=phiDeg(px,py)

        double daTag=openingDeg(tagx,tagy,tagz,px,py,pz)
        if (daTag>=0 && daTag<bestTagDa) {
            bestTagDa=daTag
            bestTag=[pid:pid,parent:parentPid,p:pp,theta:th,phi:ph,da:daTag,index:i]
        }

        if (pid==22) {
            double daProbe=openingDeg(predx,predy,predz,px,py,pz)
            if (daProbe>=0 && daProbe<bestProbeAnyDa) {
                bestProbeAnyDa=daProbe
                bestProbeAny=[pid:pid,parent:parentPid,p:pp,theta:th,phi:ph,da:daProbe,index:i]
            }
            if (parentPid==111 && daProbe>=0 && daProbe<bestProbePi0Da) {
                bestProbePi0Da=daProbe
                bestProbePi0=[pid:pid,parent:parentPid,p:pp,theta:th,phi:ph,da:daProbe,index:i]
            }
        }
    }
    Map probe=(bestProbePi0.pid!=ISENT) ? bestProbePi0 : bestProbeAny
    return [probe:probe, tag:bestTag]
}


// Integer code stored in the ROOT tree.  Keep the string at the command line so
// production commands remain self-documenting.
static int sampleKindCode(String kind) {
    String k=(kind==null?"auto":kind.trim().toLowerCase())
    if (k=="data") return 0
    if (k=="aaogen") return 1
    if (k=="clasdis") return 2
    if (k=="dvcsgen") return 3
    return 9
}

static int lundParentPid(HipoDataBank lund, int lundIndex) {
    if (lund==null || lundIndex<0 || lundIndex>=lund.rows()) return ISENT
    int parent=lund.getInt("parent",lundIndex)
    if (parent<=0 || parent-1>=lund.rows()) return ISENT
    return lund.getInt("pid",parent-1)
}

// CLASDIS generator-topology classification is done entirely in MC::Lund.
// Do NOT assume that row i of MC::Particle corresponds to row i of MC::Lund:
// they are different banks with different row numbering.
//
// In the CLASDIS LUND record, type==1 denotes a stable final-state particle.
// The beam/target, virtual photon, resonances, and decayed pi0 therefore do not
// enter the final-state multiplicity count.  A pure exclusive e p pi0 event is
// identified after pi0 decay as exactly
//
//      e + p + gamma + gamma,
//
// where both final-state photons have a direct pi0 (pid 111) parent.
// Any additional stable hadron or photon makes the event part of the
// inclusive/extra-particle CLASDIS component.
//
// topology: -1 = not classified (non-CLASDIS sample),
//            0 = insufficient truth information,
//            1 = exclusive e p pi0 -> e p gamma gamma,
//            2 = inclusive / extra-particle final state.
static Map generatedTopology(HipoDataEvent event, boolean classifyClasdis) {
    Map out=[topology:(classifyClasdis?0:-1), ne:0, npositron:0, np:0, nantiproton:0,
             nneutron:0, ngamma:0, npi0gamma:0, npiplus:0, npiminus:0,
             nkplus:0, nkminus:0, npi0:0, nother:0,
             hasParticle:0, hasRecMatch:0, hasLund:0]

    if (event.hasBank("MC::Particle")) out.hasParticle=1
    if (event.hasBank("MC::RecMatch")) out.hasRecMatch=1
    if (event.hasBank("MC::Lund")) out.hasLund=1
    if (!event.hasBank("MC::Lund")) return out

    HipoDataBank lund=(HipoDataBank)event.getBank("MC::Lund")
    // Count pi0 ancestors whether or not they are stable final-state rows.
    for (int i=0;i<lund.rows();i++) if (lund.getInt("pid",i)==111) out.npi0++

    for (int i=0;i<lund.rows();i++) {
        int type=lund.getInt("type",i)
        if (type!=1) continue
        int pid=lund.getInt("pid",i)
        if (pid==11) out.ne++
        else if (pid==2212) out.np++
        else if (pid==22) { out.ngamma++; if (lundParentPid(lund,i)==111) out.npi0gamma++ }
        else {
            // Preserve the legacy meaning of gen_n_other: every stable final-state
            // particle other than e-, p, or gamma.  Detailed species counters are
            // additional diagnostics, not a redefinition of this branch.
            out.nother++
            if (pid==-11) out.npositron++
            else if (pid==-2212) out.nantiproton++
            else if (pid==2112) out.nneutron++
            else if (pid==211) out.npiplus++
            else if (pid==-211) out.npiminus++
            else if (pid==321) out.nkplus++
            else if (pid==-321) out.nkminus++
        }
    }

    // Only explicitly-declared CLASDIS is assigned exclusive/inclusive topology.
    if (classifyClasdis) {
        boolean exclusive=(out.ne==1 && out.np==1 && out.ngamma==2 && out.npi0gamma==2 && out.nother==0)
        out.topology=exclusive?1:2
    }
    return out
}

static Map emptyMCTruth() {
    return [index:ISENT,pid:ISENT,parent:ISENT,p:SENT,theta:SENT,phi:SENT,
            daPred:SENT,mcMatchCount:0,
            recIndex:ISENT,recPid:ISENT,recCharge:ISENT,recStatus:ISENT,recDetector:ISENT,
            recP:SENT,recTheta:SENT,recPhi:SENT,recDaTruth:SENT,
            nRecMatches:0,nRecPid22:0,nRecPid11:0,nRecPid2112:0,nRecNeutral:0,nRecFT:0,nRecFD:0,
            nRecPid22FT:0,nRecPid22FD:0,nRecPid11FT:0,nRecPid11FD:0,
            recBeta:SENT,recChi2pid:SENT,recX:SENT,recY:SENT,recZ:SENT,recResponseE:SENT,
            recPcalE:SENT,recEcinE:SENT,recEcoutE:SENT,recEcalE:SENT,recFtE:SENT,recFtR:SENT,
            mcDaRec:SENT]
}

// MC::Particle-based truth, available uniformly across the MC generators.
// A generated probe photon is chosen only for the tag-and-probe hypothesis; the
// companion event tree separately stores all generated photons for unbiased truth closure.
static Map mcParticleProbeMatch(HipoDataEvent event, HipoDataBank rec,
                                double predx,double predy,double predz, int tagRecIndex,
                                Map<Integer,Map> recMcCache, Map<Integer,Map> mcRecCache) {
    Map best=emptyMCTruth()
    if (!event.hasBank("MC::Particle")) return best
    HipoDataBank mc=(HipoDataBank)event.getBank("MC::Particle")
    Map tagTruth=recMcCache.containsKey(tagRecIndex) ? recMcCache[tagRecIndex] : emptyRecToMc()
    int tagMc=(int)tagTruth.index
    double bestDa=1.0e99
    int bestIndex=ISENT
    for (int i=0;i<mc.rows();i++) {
        if (i==tagMc || mc.getInt("pid",i)!=22) continue
        double px=mc.getFloat("px",i),py=mc.getFloat("py",i),pz=mc.getFloat("pz",i)
        double da=openingDeg(predx,predy,predz,px,py,pz)
        if (da>=0 && da<bestDa) { bestDa=da; bestIndex=i }
    }
    if (bestIndex==ISENT) return best
    double px=mc.getFloat("px",bestIndex),py=mc.getFloat("py",bestIndex),pz=mc.getFloat("pz",bestIndex)
    best.index=bestIndex; best.pid=mc.getInt("pid",bestIndex); best.parent=ISENT
    best.p=p3(px,py,pz); best.theta=thetaDeg(px,py,pz); best.phi=phiDeg(px,py); best.daPred=bestDa
    Map mr=mcRecCache.containsKey(bestIndex) ? mcRecCache[bestIndex] : emptyMcToRec()
    best.nRecMatches=mr.matchCount; best.nRecPid22=mr.nPid22; best.nRecPid11=mr.nPid11; best.nRecPid2112=mr.nPid2112
    best.nRecNeutral=mr.nNeutral; best.nRecFT=mr.nFT; best.nRecFD=mr.nFD
    best.nRecPid22FT=mr.nPid22FT; best.nRecPid22FD=mr.nPid22FD; best.nRecPid11FT=mr.nPid11FT; best.nRecPid11FD=mr.nPid11FD
    best.recIndex=mr.recIndex; best.recPid=mr.recPid; best.recCharge=mr.recCharge; best.recStatus=mr.recStatus; best.recDetector=mr.recDetector
    best.recP=mr.recP; best.recTheta=mr.recTheta; best.recPhi=mr.recPhi; best.recDaTruth=mr.recDaTruth
    best.recBeta=mr.recBeta; best.recChi2pid=mr.recChi2pid; best.recX=mr.recX; best.recY=mr.recY; best.recZ=mr.recZ; best.recResponseE=mr.recResponseE
    best.recPcalE=mr.recPcalE; best.recEcinE=mr.recEcinE; best.recEcoutE=mr.recEcoutE; best.recEcalE=mr.recEcalE; best.recFtE=mr.recFtE; best.recFtR=mr.recFtR
    return best
}

static Map mcParticleTagMatch(HipoDataEvent event, HipoDataBank rec, int tagRecIndex, Map<Integer,Map> recMcCache) {
    Map out=emptyMCTruth()
    Map m=recMcCache.containsKey(tagRecIndex) ? recMcCache[tagRecIndex] : emptyRecToMc()
    if ((int)m.index==ISENT) { out.mcMatchCount=m.matchCount; return out }
    out.index=m.index; out.pid=m.pid; out.p=m.p; out.theta=m.theta; out.phi=m.phi; out.mcMatchCount=m.matchCount; out.mcDaRec=m.da
    out.recIndex=tagRecIndex
    return out
}

static long sourceFileHash(String path) {
    CRC32 crc=new CRC32()
    byte[] bytes=path.getBytes("UTF-8")
    crc.update(bytes,0,bytes.length)
    return crc.getValue()
}

// Compact detector-response summary for one REC::Particle.  These quantities
// are intentionally stored now because they cannot be reconstructed later from
// the skim if we need threshold/fiducial/PID diagnostics.
static Map recResponseSummary(int pindex, HipoDataBank rec, HipoDataBank cal, HipoDataBank ft) {
    Map out=[beta:SENT,chi2pid:SENT,x:SENT,y:SENT,z:SENT,responseE:SENT,
             pcalE:0.0,ecinE:0.0,ecoutE:0.0,ecalE:0.0,pcalSector:ISENT,
             pcalLu:SENT,pcalLv:SENT,pcalLw:SENT,ftE:SENT,ftR:SENT]
    if (rec!=null && pindex>=0 && pindex<rec.rows()) {
        out.beta=rec.getFloat("beta",pindex)
        out.chi2pid=rec.getFloat("chi2pid",pindex)
        int det=detectorFromStatus(rec.getInt("status",pindex))
        double[] xyz=responseXYZ(pindex,det,cal,ft)
        out.x=xyz[0]; out.y=xyz[1]; out.z=xyz[2]
    }
    if (cal!=null) {
        for (int r=0;r<cal.rows();r++) {
            if (cal.getInt("pindex",r)!=pindex) continue
            int layer=cal.getInt("layer",r)
            double e=cal.getFloat("energy",r)
            out.ecalE=((double)out.ecalE)+e
            if (layer==1) {
                out.pcalE=((double)out.pcalE)+e
                out.pcalSector=cal.getInt("sector",r)
                out.pcalLu=cal.getFloat("lu",r); out.pcalLv=cal.getFloat("lv",r); out.pcalLw=cal.getFloat("lw",r)
            } else if (layer==4) {
                out.ecinE=((double)out.ecinE)+e
            } else if (layer==7) {
                out.ecoutE=((double)out.ecoutE)+e
            }
        }
    }
    if (ft!=null) {
        for (int r=0;r<ft.rows();r++) {
            if (ft.getInt("pindex",r)!=pindex) continue
            out.ftE=ft.getFloat("energy",r)
            out.ftR=ft.getFloat("radius",r)
            break
        }
    }
    int det=(rec!=null && pindex>=0 && pindex<rec.rows()) ? detectorFromStatus(rec.getInt("status",pindex)) : -1
    if (det==0 && ((double)out.ftE)>-900) out.responseE=out.ftE
    else if (det==1) out.responseE=out.ecalE
    return out
}

static List<Integer> uniqueMcMatchesForRec(HipoDataBank match, int recIndex) {
    List<Integer> out=[]
    if (match==null || recIndex<0) return out
    Set<Integer> seen=new LinkedHashSet<Integer>()
    for (int r=0;r<match.rows();r++) {
        if (match.getInt("pindex",r)!=recIndex) continue
        int mi=match.getInt("mcindex",r)
        if (mi>=0 && !seen.contains(mi)) { seen.add(mi); out.add(mi) }
    }
    return out
}

static List<Integer> uniqueRecMatchesForMc(HipoDataBank match, int mcIndex) {
    List<Integer> out=[]
    if (match==null || mcIndex<0) return out
    Set<Integer> seen=new LinkedHashSet<Integer>()
    for (int r=0;r<match.rows();r++) {
        if (match.getInt("mcindex",r)!=mcIndex) continue
        int ri=match.getInt("pindex",r)
        if (ri>=0 && !seen.contains(ri)) { seen.add(ri); out.add(ri) }
    }
    return out
}

static Map emptyRecToMc() {
    return [matchCount:0,index:ISENT,pid:ISENT,p:SENT,theta:SENT,phi:SENT,da:SENT]
}

// For one reconstructed object, select the angularly closest valid MC::Particle
// association while retaining the total number of distinct MC associations.
static Map bestMcForRec(HipoDataEvent event, HipoDataBank rec, HipoDataBank cal, HipoDataBank ft,
                        int recIndex, double vxRef,double vyRef,double vzRef) {
    Map out=emptyRecToMc()
    if (!event.hasBank("MC::Particle") || !event.hasBank("MC::RecMatch") ||
        recIndex<0 || recIndex>=rec.rows()) return out
    HipoDataBank mc=(HipoDataBank)event.getBank("MC::Particle")
    HipoDataBank match=(HipoDataBank)event.getBank("MC::RecMatch")
    List<Integer> mids=uniqueMcMatchesForRec(match,recIndex)
    out.matchCount=mids.size()
    double[] rdir=candidateDirection(recIndex,rec,cal,ft,vxRef,vyRef,vzRef)
    double bestDa=1.0e99
    for (int mi : mids) {
        if (mi<0 || mi>=mc.rows()) continue
        double px=mc.getFloat("px",mi),py=mc.getFloat("py",mi),pz=mc.getFloat("pz",mi)
        double da=openingDeg(rdir[0],rdir[1],rdir[2],px,py,pz)
        if (da<0) da=1.0e98
        if (out.index==ISENT || da<bestDa) {
            bestDa=da; out.index=mi; out.pid=mc.getInt("pid",mi)
            out.p=p3(px,py,pz); out.theta=thetaDeg(px,py,pz); out.phi=phiDeg(px,py)
            out.da=(da<1.0e90?da:SENT)
        }
    }
    return out
}

static Map emptyMcToRec() {
    return [matchCount:0,nPid22:0,nPid11:0,nPid2112:0,nNeutral:0,nFT:0,nFD:0,
            nPid22FT:0,nPid22FD:0,nPid11FT:0,nPid11FD:0,
            recIndex:ISENT,recPid:ISENT,recCharge:ISENT,recStatus:ISENT,recDetector:ISENT,
            recP:SENT,recTheta:SENT,recPhi:SENT,recDaTruth:SENT,
            recBeta:SENT,recChi2pid:SENT,recX:SENT,recY:SENT,recZ:SENT,recResponseE:SENT,
            recPcalE:SENT,recEcinE:SENT,recEcoutE:SENT,recEcalE:SENT,recFtE:SENT,recFtR:SENT]
}

// For one generated particle, summarize every distinct REC association and also
// keep the angularly closest reconstructed object.  This avoids the old
// arbitrary "first MC::RecMatch row" behavior while preserving enough counts to
// test alternative definitions offline.
static Map bestRecForMc(HipoDataEvent event, HipoDataBank rec, HipoDataBank cal, HipoDataBank ft,
                        int mcIndex, double vxRef,double vyRef,double vzRef) {
    Map out=emptyMcToRec()
    if (!event.hasBank("MC::Particle") || !event.hasBank("MC::RecMatch") || mcIndex<0) return out
    HipoDataBank mc=(HipoDataBank)event.getBank("MC::Particle")
    HipoDataBank match=(HipoDataBank)event.getBank("MC::RecMatch")
    if (mcIndex>=mc.rows()) return out
    double mpx=mc.getFloat("px",mcIndex),mpy=mc.getFloat("py",mcIndex),mpz=mc.getFloat("pz",mcIndex)
    List<Integer> rids=uniqueRecMatchesForMc(match,mcIndex)
    out.matchCount=rids.size()
    double bestDa=1.0e99
    for (int ri : rids) {
        if (ri<0 || ri>=rec.rows()) continue
        int pid=rec.getInt("pid",ri), charge=rec.getByte("charge",ri), status=rec.getInt("status",ri)
        int det=detectorFromStatus(status)
        if (pid==22) out.nPid22++
        if (pid==11) out.nPid11++
        if (pid==2112) out.nPid2112++
        if (charge==0) out.nNeutral++
        if (det==0) out.nFT++
        if (det==1) out.nFD++
        if (pid==22 && det==0) out.nPid22FT++
        if (pid==22 && det==1) out.nPid22FD++
        if (pid==11 && det==0) out.nPid11FT++
        if (pid==11 && det==1) out.nPid11FD++
        double[] rdir=candidateDirection(ri,rec,cal,ft,vxRef,vyRef,vzRef)
        double da=openingDeg(mpx,mpy,mpz,rdir[0],rdir[1],rdir[2])
        if (da<0) da=1.0e98
        if (out.recIndex==ISENT || da<bestDa) {
            bestDa=da
            double rpx=rec.getFloat("px",ri),rpy=rec.getFloat("py",ri),rpz=rec.getFloat("pz",ri)
            Map rs=recResponseSummary(ri,rec,cal,ft)
            out.recIndex=ri; out.recPid=pid; out.recCharge=charge; out.recStatus=status; out.recDetector=det
            out.recP=p3(rpx,rpy,rpz); out.recTheta=thetaDeg(rdir[0],rdir[1],rdir[2]); out.recPhi=phiDeg(rdir[0],rdir[1])
            out.recDaTruth=(da<1.0e90?da:SENT)
            out.recBeta=rs.beta; out.recChi2pid=rs.chi2pid; out.recX=rs.x; out.recY=rs.y; out.recZ=rs.z; out.recResponseE=rs.responseE
            out.recPcalE=rs.pcalE; out.recEcinE=rs.ecinE; out.recEcoutE=rs.ecoutE; out.recEcalE=rs.ecalE; out.recFtE=rs.ftE; out.recFtR=rs.ftR
        }
    }
    return out
}

// Event-level reverse-matching bookkeeping.  Counts reconstructed photon/electron
// candidates by the generated PID selected from MC::RecMatch.  These branches
// directly answer questions such as e_gen -> gamma_REC and gamma_gen -> e_REC.
static Map eventRecTruthCounts(HipoDataBank rec, Map<Integer,Map> recMcCache) {
    Map o=[pid22Total:0,pid22Unmatched:0,pid22Multi:0,pid22From22:0,pid22From11:0,pid22FromOther:0,
           ftPid22:0,ftPid22From22:0,ftPid22From11:0,ftPid22Unmatched:0,
           fdPid22:0,fdPid22From22:0,fdPid22From11:0,fdPid22Unmatched:0,
           ftPid11:0,ftPid11From22:0,fdPid11:0,fdPid11From22:0,
           fdPid2112:0,fdPid2112From22:0]
    for (int i=0;i<rec.rows();i++) {
        int pid=rec.getInt("pid",i), det=detectorFromStatus(rec.getInt("status",i))
        if (!(pid==22 || pid==11 || pid==2112)) continue
        Map m=recMcCache.containsKey(i) ? recMcCache[i] : emptyRecToMc()
        int mpid=(int)m.pid
        if (pid==22) {
            o.pid22Total++
            if ((int)m.matchCount==0) o.pid22Unmatched++
            if ((int)m.matchCount>1) o.pid22Multi++
            if (mpid==22) { o.pid22From22++ }
            else if (mpid==11) { o.pid22From11++ }
            else if (mpid!=ISENT) { o.pid22FromOther++ }
            if (det==0) {
                o.ftPid22++
                if (mpid==22) { o.ftPid22From22++ }
                else if (mpid==11) { o.ftPid22From11++ }
                else if ((int)m.matchCount==0) { o.ftPid22Unmatched++ }
            }
            if (det==1) {
                o.fdPid22++
                if (mpid==22) { o.fdPid22From22++ }
                else if (mpid==11) { o.fdPid22From11++ }
                else if ((int)m.matchCount==0) { o.fdPid22Unmatched++ }
            }
        }
        if (pid==11 && det==0) { o.ftPid11++; if(mpid==22)o.ftPid11From22++ }
        if (pid==11 && det==1) { o.fdPid11++; if(mpid==22)o.fdPid11From22++ }
        if (pid==2112 && det==1) { o.fdPid2112++; if(mpid==22)o.fdPid2112From22++ }
    }
    return o
}

// Store a compact list of all generated MC::Particle photons once per event in
// the companion PhotonEfficiencyEvents tree.  This makes a true generated-level
// efficiency closure possible without relying on the missing-vector probe choice.
static Map generatedPhotonRecords(HipoDataEvent event, HipoDataBank rec, HipoDataBank cal, HipoDataBank ft,
                                  double vxRef,double vyRef,double vzRef) {
    List<Map> out=[]
    Map<Integer,Map> byIndex=[:]
    int total=0
    if (!event.hasBank("MC::Particle")) return [total:0,overflow:0,records:out,byIndex:byIndex]
    HipoDataBank mc=(HipoDataBank)event.getBank("MC::Particle")
    for (int i=0;i<mc.rows();i++) {
        if (mc.getInt("pid",i)!=22) continue
        total++
        double px=mc.getFloat("px",i),py=mc.getFloat("py",i),pz=mc.getFloat("pz",i)
        Map mr=bestRecForMc(event,rec,cal,ft,i,vxRef,vyRef,vzRef)
        byIndex[i]=mr
        if (out.size()<N_GEN_GAMMA_SAVE) {
            Map g=[index:i,p:p3(px,py,pz),theta:thetaDeg(px,py,pz),phi:phiDeg(px,py),
                   vx:mc.getFloat("vx",i),vy:mc.getFloat("vy",i),vz:mc.getFloat("vz",i),rec:mr]
            out.add(g)
        }
    }
    return [total:total,overflow:Math.max(0,total-N_GEN_GAMMA_SAVE),records:out,byIndex:byIndex]
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


static Map correctedPhotonKinematics(Map c, HipoDataBank rec, HipoDataBank run, energy_loss_corrections eloss) {
    Map out=[p:SENT,theta:SENT,phi:SENT]
    if (c==null || (int)c.pid!=22) return out
    int idx=(int)c.idx
    if (idx<0 || idx>=rec.rows()) return out
    float[] v=[rec.getFloat("px",idx),rec.getFloat("py",idx),rec.getFloat("pz",idx)] as float[]
    eloss.sebastian_photon_energy_loss_corrections(idx,v,rec,run)
    out.p=p3(v[0],v[1],v[2]); out.theta=thetaDeg(v[0],v[1],v[2]); out.phi=phiDeg(v[0],v[1])
    return out
}

static void appendCandidateExtra(List vals, Map c, Map m, Map corr) {
    if (c==null) {
        for (int j=0;j<23;j++) vals.add(ISENT)
        return
    }
    vals.add(c.beta); vals.add(c.chi2pid); vals.add(c.responseE); vals.add(c.pcalE); vals.add(c.ecinE); vals.add(c.ecoutE); vals.add(c.ecalE)
    vals.add(c.pcalSector); vals.add(c.pcalLu); vals.add(c.pcalLv); vals.add(c.pcalLw); vals.add(c.ftE); vals.add(c.ftR)
    vals.add(corr!=null?corr.p:SENT); vals.add(corr!=null?corr.theta:SENT); vals.add(corr!=null?corr.phi:SENT)
    if (m==null) { vals.add(0); vals.add(ISENT); vals.add(ISENT); vals.add(SENT); vals.add(SENT); vals.add(SENT); vals.add(SENT) }
    else { vals.add(m.matchCount); vals.add(m.index); vals.add(m.pid); vals.add(m.p); vals.add(m.theta); vals.add(m.phi); vals.add(m.da) }
}

static void processPhotonEfficiency(String[] args) {
    if (args.length < 2) {
        println "Usage: process_photon_efficiency.groovy <input.hipo|dir> <output.txt> [mc_beam] [run_override] [qadb_override] [is_mc] [mx2_min] [mx2_max] [mx2_epg_min] [mx2_epg_max] [sample_kind]"
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
    double mx2EpgMin=args.length>8 ? Double.parseDouble(args[8]) : -0.25
    double mx2EpgMax=args.length>9 ? Double.parseDouble(args[9]) : 0.25
    String sampleKind=args.length>10 ? args[10].trim().toLowerCase() : "auto"
    if (!(sampleKind in ["auto","data","aaogen","clasdis","dvcsgen"])) {
        println "ERROR: sample_kind must be data|aaogen|clasdis|dvcsgen|auto"
        System.exit(3)
    }
    int sampleCode=sampleKindCode(sampleKind)

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
    File eventOutf=new File(output+".events"); eventOutf.delete()
    BufferedWriter eventWriter=new BufferedWriter(new FileWriter(eventOutf))
    StringBuilder batch=new StringBuilder(); int lineCount=0
    StringBuilder eventBatch=new StringBuilder(); int eventLineCount=0
    long nevt=0, nrow=0, neventRow=0, npair=0, ntag=0, nprobePositive=0, nprobeMx2=0
    long nClasdisExclusiveSkipped=0, nClasdisInclusiveKept=0, nClasdisUnknownSkipped=0

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

            // Explicit source selection is preferred.  `auto` never drops events:
            // it only records generic truth-bank availability.
            boolean classifyClasdis=(sampleKind=="clasdis")
            Map genTopo=generatedTopology(event,classifyClasdis)
            if (classifyClasdis) {
                if ((genTopo.topology as int)==1) { nClasdisExclusiveSkipped++; continue }
                if ((genTopo.topology as int)==2) nClasdisInclusiveKept++
                else { nClasdisUnknownSkipped++; continue }
            }

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

            // Companion event-level record.  This is written exactly once per accepted
            // electron event and therefore avoids the hypothesis duplication of the main tree.
            long sourceHash=sourceFileHash(hf.absolutePath)
            double evx=rec.getFloat("vx",0), evy=rec.getFloat("vy",0)
            Map<Integer,Map> recMcCache=[:]
            for (int ii=0;ii<rec.rows();ii++) recMcCache[ii]=bestMcForRec(event,rec,cal,ft,ii,evx,evy,evz)
            Map eTruthEvent=recMcCache.containsKey(0) ? recMcCache[0] : emptyRecToMc()
            Map recTruthCounts=eventRecTruthCounts(rec,recMcCache)
            Map genGammas=generatedPhotonRecords(event,rec,cal,ft,evx,evy,evz)
            Map<Integer,Map> mcRecCache=(Map<Integer,Map>)genGammas.byIndex
            List evVals=[]
            evVals.add(SKIM_VERSION); evVals.add(sourceHash); evVals.add(runnum); evVals.add(evnum); evVals.add(helicity); evVals.add(isMC); evVals.add(sampleCode)
            evVals.add(genTopo.topology); evVals.add(genTopo.hasParticle); evVals.add(genTopo.hasRecMatch); evVals.add(genTopo.hasLund)
            evVals.add(genTopo.ne); evVals.add(genTopo.npositron); evVals.add(genTopo.np); evVals.add(genTopo.nantiproton); evVals.add(genTopo.nneutron)
            evVals.add(genTopo.ngamma); evVals.add(genTopo.npi0gamma); evVals.add(genTopo.npiplus); evVals.add(genTopo.npiminus); evVals.add(genTopo.nkplus); evVals.add(genTopo.nkminus); evVals.add(genTopo.npi0); evVals.add(genTopo.nother)
            evVals.add(beamE); evVals.add(torus); evVals.add(solenoid); evVals.add(mcWeight(event))
            evVals.add(Q2); evVals.add(W); evVals.add(xB); evVals.add(yy)
            evVals.add(ep); evVals.add(thetaDeg(epx,epy,epz)); evVals.add(phiDeg(epx,epy)); evVals.add(evx); evVals.add(evy); evVals.add(evz)
            evVals.add(rec.rows()); evVals.add(nNeutral); evVals.add(nPid22); evVals.add(nPid2112); evVals.add(nPid0Neutral)
            evVals.add(eTruthEvent.matchCount); evVals.add(eTruthEvent.index); evVals.add(eTruthEvent.pid); evVals.add(eTruthEvent.p); evVals.add(eTruthEvent.theta); evVals.add(eTruthEvent.phi); evVals.add(eTruthEvent.da)
            evVals.add(recTruthCounts.pid22Total); evVals.add(recTruthCounts.pid22Unmatched); evVals.add(recTruthCounts.pid22Multi); evVals.add(recTruthCounts.pid22From22); evVals.add(recTruthCounts.pid22From11); evVals.add(recTruthCounts.pid22FromOther)
            evVals.add(recTruthCounts.ftPid22); evVals.add(recTruthCounts.ftPid22From22); evVals.add(recTruthCounts.ftPid22From11); evVals.add(recTruthCounts.ftPid22Unmatched)
            evVals.add(recTruthCounts.fdPid22); evVals.add(recTruthCounts.fdPid22From22); evVals.add(recTruthCounts.fdPid22From11); evVals.add(recTruthCounts.fdPid22Unmatched)
            evVals.add(recTruthCounts.ftPid11); evVals.add(recTruthCounts.ftPid11From22); evVals.add(recTruthCounts.fdPid11); evVals.add(recTruthCounts.fdPid11From22); evVals.add(recTruthCounts.fdPid2112); evVals.add(recTruthCounts.fdPid2112From22)
            evVals.add(genGammas.total); evVals.add(genGammas.overflow)
            List<Map> gg=(List<Map>)genGammas.records
            for (int k=0;k<N_GEN_GAMMA_SAVE;k++) {
                if (k<gg.size()) {
                    Map g=gg[k]; Map mr=(Map)g.rec
                    evVals.add(g.index); evVals.add(g.p); evVals.add(g.theta); evVals.add(g.phi); evVals.add(g.vx); evVals.add(g.vy); evVals.add(g.vz)
                    evVals.add(mr.matchCount); evVals.add(mr.nPid22); evVals.add(mr.nPid11); evVals.add(mr.nPid2112); evVals.add(mr.nNeutral); evVals.add(mr.nFT); evVals.add(mr.nFD); evVals.add(mr.nPid22FT); evVals.add(mr.nPid22FD); evVals.add(mr.nPid11FT); evVals.add(mr.nPid11FD)
                    evVals.add(mr.recIndex); evVals.add(mr.recPid); evVals.add(mr.recCharge); evVals.add(mr.recStatus); evVals.add(mr.recDetector); evVals.add(mr.recP); evVals.add(mr.recTheta); evVals.add(mr.recPhi); evVals.add(mr.recDaTruth); evVals.add(mr.recBeta); evVals.add(mr.recChi2pid); evVals.add(mr.recX); evVals.add(mr.recY); evVals.add(mr.recZ); evVals.add(mr.recResponseE); evVals.add(mr.recPcalE); evVals.add(mr.recEcinE); evVals.add(mr.recEcoutE); evVals.add(mr.recEcalE); evVals.add(mr.recFtE); evVals.add(mr.recFtR)
                } else {
                    // 39 fields per generated-photon slot: 7 truth + 11 association counts + 21 best-REC quantities.
                    for (int j=0;j<39;j++) evVals.add(ISENT)
                }
            }
            eventBatch.append(evVals.collect{fmt(it)}.join(' ')).append('\n'); eventLineCount++; neventRow++
            if (eventLineCount>=500) { eventWriter.write(eventBatch.toString()); eventBatch.setLength(0); eventLineCount=0 }

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
                    Map tagSummary=recResponseSummary(ig,rec,cal,ft)

                    LorentzVector tagRaw=new LorentzVector(); tagRaw.setPxPyPzM(gpx,gpy,gpz,0)
                    float[] gcorr=[(float)gpx,(float)gpy,(float)gpz] as float[]
                    eloss.sebastian_photon_energy_loss_corrections(ig,gcorr,rec,run)
                    LorentzVector tagCorr=new LorentzVector(); tagCorr.setPxPyPzM(gcorr[0],gcorr[1],gcorr[2],0)

                    ntag++
                    LorentzVector probeRaw=new LorentzVector(beam); probeRaw.add(target); probeRaw.sub(ele); probeRaw.sub(prot); probeRaw.sub(tagRaw)
                    // A missing photon hypothesis with non-positive energy is not a physical
                    // tag-and-probe denominator and is discarded before any matching work.
                    if (probeRaw.e() <= 0.0) continue
                    nprobePositive++

                    // probeRaw.mass2() is exactly Mx2(ep gamma_tag) for the raw tag-photon
                    // four-vector.  This broad configurable window suppresses generic SIDIS
                    // combinatorics while remaining vastly wider than the gamma missing-mass
                    // peak used in the final offline pi0 selection.
                    double mx2epg=probeRaw.mass2()
                    if (mx2epg<mx2EpgMin || mx2epg>mx2EpgMax) continue
                    nprobeMx2++

                    LorentzVector probeCorr=new LorentzVector(beam); probeCorr.add(target); probeCorr.sub(ele); probeCorr.sub(prot); probeCorr.sub(tagCorr)
                    double predx=probeRaw.px(), predy=probeRaw.py(), predz=probeRaw.pz()

                    Map candSets=nearestCandidateSets(rec,cal,ft,predx,predy,predz,0,0,evz,0,ip,ig)
                    List<Map> neutrals=(List<Map>)candSets.neutrals
                    List<Map> any=(List<Map>)candSets.any

                    Map truth=truthMatches(event,predx,predy,predz,gpx,gpy,gpz)
                    Map truthProbeAny=(Map)truth.probe
                    Map truthTag=(Map)truth.tag

                    // Uniform generator truth from MC::Particle plus the actual
                    // REC association from MC::RecMatch.  These are independent
                    // of the older MC::Lund-nearest-direction diagnostics above.
                    Map mcProbe=mcParticleProbeMatch(event,rec,predx,predy,predz,ig,recMcCache,mcRecCache)
                    Map mcTag=mcParticleTagMatch(event,rec,ig,recMcCache)
                    Map mcElectron=eTruthEvent
                    Map mcProton=recMcCache.containsKey(ip) ? recMcCache[ip] : emptyRecToMc()

                    List vals=[]
                    vals.add(runnum); vals.add(evnum); vals.add(helicity); vals.add(isMC)
                    vals.add(sampleCode); vals.add(genTopo.topology); vals.add(genTopo.hasParticle); vals.add(genTopo.hasRecMatch); vals.add(genTopo.hasLund)
                    vals.add(genTopo.ne); vals.add(genTopo.np); vals.add(genTopo.ngamma); vals.add(genTopo.npi0gamma); vals.add(genTopo.nother)
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
                    vals.add(mx2epg); vals.add(probeRaw.e()); vals.add(probeRaw.p()); vals.add(thetaDeg(predx,predy,predz)); vals.add(phiDeg(predx,predy))
                    vals.add(probeRaw.px()); vals.add(probeRaw.py()); vals.add(probeRaw.pz())
                    vals.add(probeCorr.mass2()); vals.add(probeCorr.e()); vals.add(probeCorr.p()); vals.add(thetaDeg(probeCorr.px(),probeCorr.py(),probeCorr.pz())); vals.add(phiDeg(probeCorr.px(),probeCorr.py()))
                    vals.add(probeCorr.px()); vals.add(probeCorr.py()); vals.add(probeCorr.pz())

                    for (int k=0;k<N_NEUTRAL_SAVE;k++) appendCandidate(vals,k<neutrals.size()?neutrals[k]:null)
                    for (int k=0;k<N_ANY_SAVE;k++) appendCandidate(vals,k<any.size()?any[k]:null)

                    vals.add(truthProbeAny.index); vals.add(truthProbeAny.pid); vals.add(truthProbeAny.parent); vals.add(truthProbeAny.p); vals.add(truthProbeAny.theta); vals.add(truthProbeAny.phi); vals.add(truthProbeAny.da)
                    vals.add(truthTag.index); vals.add(truthTag.pid); vals.add(truthTag.parent); vals.add(truthTag.p); vals.add(truthTag.theta); vals.add(truthTag.phi); vals.add(truthTag.da)

                    // MC::Particle / MC::RecMatch truth branches.
                    vals.add(mcProbe.index); vals.add(mcProbe.pid); vals.add(mcProbe.parent); vals.add(mcProbe.p); vals.add(mcProbe.theta); vals.add(mcProbe.phi); vals.add(mcProbe.daPred)
                    vals.add(mcProbe.recIndex); vals.add(mcProbe.recPid); vals.add(mcProbe.recCharge); vals.add(mcProbe.recStatus); vals.add(mcProbe.recDetector)
                    vals.add(mcProbe.recP); vals.add(mcProbe.recTheta); vals.add(mcProbe.recPhi); vals.add(mcProbe.recDaTruth)
                    vals.add(mcTag.index); vals.add(mcTag.pid); vals.add(mcTag.parent); vals.add(mcTag.p); vals.add(mcTag.theta); vals.add(mcTag.phi)

                    // Final-batch diagnostics appended after the legacy schema.
                    vals.add(SKIM_VERSION); vals.add(sourceHash)
                    vals.add(genTopo.npositron); vals.add(genTopo.nantiproton); vals.add(genTopo.nneutron); vals.add(genTopo.npiplus); vals.add(genTopo.npiminus); vals.add(genTopo.nkplus); vals.add(genTopo.nkminus); vals.add(genTopo.npi0)
                    vals.add(mcElectron.matchCount); vals.add(mcElectron.index); vals.add(mcElectron.pid); vals.add(mcElectron.p); vals.add(mcElectron.theta); vals.add(mcElectron.phi); vals.add(mcElectron.da)
                    vals.add(mcProton.matchCount); vals.add(mcProton.index); vals.add(mcProton.pid); vals.add(mcProton.p); vals.add(mcProton.theta); vals.add(mcProton.phi); vals.add(mcProton.da)
                    vals.add(mcTag.mcMatchCount); vals.add(mcTag.mcDaRec)
                    vals.add(mcProbe.nRecMatches); vals.add(mcProbe.nRecPid22); vals.add(mcProbe.nRecPid11); vals.add(mcProbe.nRecPid2112); vals.add(mcProbe.nRecNeutral); vals.add(mcProbe.nRecFT); vals.add(mcProbe.nRecFD); vals.add(mcProbe.nRecPid22FT); vals.add(mcProbe.nRecPid22FD); vals.add(mcProbe.nRecPid11FT); vals.add(mcProbe.nRecPid11FD)
                    vals.add(mcProbe.recBeta); vals.add(mcProbe.recChi2pid); vals.add(mcProbe.recX); vals.add(mcProbe.recY); vals.add(mcProbe.recZ); vals.add(mcProbe.recResponseE); vals.add(mcProbe.recPcalE); vals.add(mcProbe.recEcinE); vals.add(mcProbe.recEcoutE); vals.add(mcProbe.recEcalE); vals.add(mcProbe.recFtE); vals.add(mcProbe.recFtR)
                    vals.add(tagSummary.chi2pid); vals.add(tagSummary.responseE); vals.add(tagSummary.pcalE); vals.add(tagSummary.ecinE); vals.add(tagSummary.ecoutE); vals.add(tagSummary.ecalE); vals.add(tagSummary.pcalSector); vals.add(tagSummary.pcalLu); vals.add(tagSummary.pcalLv); vals.add(tagSummary.pcalLw); vals.add(tagSummary.ftE); vals.add(tagSummary.ftR)
                    for (int k=0;k<N_NEUTRAL_SAVE;k++) {
                        Map c=(k<neutrals.size()?neutrals[k]:null)
                        Map m=(c!=null && recMcCache.containsKey((int)c.idx) ? recMcCache[(int)c.idx] : null)
                        Map ck=correctedPhotonKinematics(c,rec,run,eloss)
                        appendCandidateExtra(vals,c,m,ck)
                    }
                    for (int k=0;k<N_ANY_SAVE;k++) {
                        Map c=(k<any.size()?any[k]:null)
                        Map m=(c!=null && recMcCache.containsKey((int)c.idx) ? recMcCache[(int)c.idx] : null)
                        Map ck=correctedPhotonKinematics(c,rec,run,eloss)
                        appendCandidateExtra(vals,c,m,ck)
                    }

                    batch.append(vals.collect{fmt(it)}.join(' ')).append('\n'); lineCount++; nrow++
                    if (lineCount>=1000) { writer.write(batch.toString()); batch.setLength(0); lineCount=0 }
                } // tag
            } // proton
            if (nevt%1000000==0) println "processed ${nevt} events: ep_pairs=${npair}, tags=${ntag}, positive_probe=${nprobePositive}, loose_Mx2_epg=${nprobeMx2}, rows=${nrow}"
        }
        reader.close()
    }
    if (batch.length()>0) writer.write(batch.toString())
    if (eventBatch.length()>0) eventWriter.write(eventBatch.toString())
    writer.close(); eventWriter.close()
    double min=(System.currentTimeMillis()-start)/60000.0
    println String.format(java.util.Locale.US,"Done: events=%d event_rows=%d accepted_ep_pairs=%d tag_candidates=%d positive_probe=%d loose_Mx2_epg=%d rows=%d time=%.2f min output=%s",nevt,neventRow,npair,ntag,nprobePositive,nprobeMx2,nrow,min,output)
    if (sampleKind=="clasdis") println String.format(java.util.Locale.US,"CLASDIS generator topology: exclusive_eppi0_skipped=%d inclusive_or_extra_kept=%d unknown_skipped=%d",nClasdisExclusiveSkipped,nClasdisInclusiveKept,nClasdisUnknownSkipped)
}

// Execute the processor from the Groovy script body.
// Do not name this method main(): Groovy generates its own script-class main(),
// and calling main(args) from the script body recursively re-enters run().
processPhotonEfficiency(args)
