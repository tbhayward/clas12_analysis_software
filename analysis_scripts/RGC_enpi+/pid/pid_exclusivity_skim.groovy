/*
 * RGC e pi+ n PID/exclusivity diagnostic skim
 *
 * Purpose:
 *   Follow every reconstructed positively charged FD hadron candidate through
 *   the e pi+ X channel selection WITHOUT requiring the REC-assigned hadron PID.
 *   All hadron-dependent kinematics are evaluated under the pi+ mass hypothesis.
 *
 * Output:
 *   whitespace-delimited worker file. One row = one positive FD hadron candidate.
 *   A companion Python script converts/merges workers into a ROOT TTree.
 *
 * Usage:
 *   run-groovy pid_exclusivity_skim.groovy input.hipo output.txt period
 *   period = Su22 | Fa22 | Sp23
 */

import org.jlab.io.hipo.*
import org.jlab.io.base.DataEvent
import org.jlab.clas.physics.*
import org.jlab.clas12.physics.*
import extended_kinematic_fitters.*
import analyzers.*
import clasqa.QADB

class PIDExclusivitySkim {

    static final double ME = 0.00051099895
    static final double MP = 0.9382720813
    static final double MN = 0.9395654133
    static final double MPI = 0.13957039

    // xB and -t' production bins: 4 x 6 = 24.
    static final double[] XB_EDGES = [0.10, 0.25, 0.35, 0.45, 0.60] as double[]
    static final double[] TP_EDGES = [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25] as double[]

    // Shared fitted Mx2 means, bins 1..24.
    static final double[] MX2_MU = [
        0.8646,0.8715,0.8630,0.8906,0.8493,0.8579,
        0.8781,0.8854,0.8806,0.8781,0.8708,0.8795,
        0.8861,0.8821,0.8816,0.8902,0.8832,0.8797,
        0.8863,0.8765,0.8833,0.8858,0.8858,0.8940
    ] as double[]

    static final Map<String,double[]> MX2_SIGMA = [
        'Su22': [0.0685,0.0756,0.0695,0.0883,0.0778,0.0329,
                 0.0756,0.0714,0.0807,0.0751,0.0556,0.0921,
                 0.0750,0.0704,0.1018,0.0649,0.0738,0.0347,
                 0.0731,0.0819,0.0614,0.0665,0.0524,0.0692] as double[],
        'Fa22': [0.0643,0.0669,0.0515,0.1144,0.0657,0.0642,
                 0.0632,0.0648,0.0565,0.0582,0.0682,0.0843,
                 0.0669,0.0654,0.0675,0.0744,0.0649,0.0690,
                 0.0651,0.0636,0.0665,0.0634,0.0756,0.0841] as double[],
        'Sp23': [0.0634,0.0556,0.0561,0.0959,0.0970,0.0616,
                 0.0635,0.0662,0.0683,0.0652,0.0679,0.0798,
                 0.0697,0.0678,0.0779,0.0714,0.0692,0.0561,
                 0.0528,0.0570,0.0593,0.0697,0.0505,0.0779] as double[]
    ]

    static double sq(double x) { return x*x }

    static double mass2(double E, double px, double py, double pz) {
        return E*E - px*px - py*py - pz*pz
    }

    static int analysisBin(double xb, double mtp) {
        int ix = -1, it = -1
        for (int i=0; i<XB_EDGES.length-1; ++i) {
            if (xb >= XB_EDGES[i] && xb < XB_EDGES[i+1]) { ix=i; break }
        }
        for (int i=0; i<TP_EDGES.length-1; ++i) {
            if (mtp >= TP_EDGES[i] && mtp < TP_EDGES[i+1]) { it=i; break }
        }
        return (ix >= 0 && it >= 0) ? ix*6 + it + 1 : -1
    }

    // Mesonic t_min for gamma* p -> pi+ n. This is the forward-pion value of t.
    static double mesonicTmin(double Q2, double W) {
        if (!(Q2 > 0.0) || !(W >= MPI + MN)) return Double.NaN
        double W2 = W*W
        double q0s = (W2 - MP*MP - Q2)/(2.0*W)
        double qs  = Math.sqrt(Math.max(q0s*q0s + Q2, 0.0))
        double epis = (W2 + MPI*MPI - MN*MN)/(2.0*W)
        double ppis = Math.sqrt(Math.max(epis*epis - MPI*MPI, 0.0))
        return -Q2 + MPI*MPI - 2.0*(q0s*epis - qs*ppis)
    }

    static boolean isFD(int status) {
        int s = Math.abs(status)
        return s >= 2000 && s < 4000
    }

    static QADB makeQA() {
        QADB qa = new QADB('latest')
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
        [16194,16089,16185,16308,16184,16307,16309,
         16872,16975,17763,17764,17765,17766,17767,17768,
         17179,17180,17181,17182,17183,17188,17189,17252].each { qa.allowMiscBit(it) }
        return qa
    }

    static boolean rejectRun(int run) {
        if (run > 16600 && run < 16700) return true       // Hall-C bleedthrough
        if (run > 17768 && run <= 17811) return true     // Sp23 outbending
        if ([17331,16987,17079,17190,17639].contains(run)) return true
        if ([16850,16851,16852,16855,16879].contains(run)) return true
        return false
    }

    static int electronIndex(HipoDataBank rec) {
        // analysis_fitter convention is led by the reconstructed electron. For this
        // diagnostic retain the highest-momentum REC pid=11 candidate if >1 exists.
        int best=-1; double bestp=-1
        for (int i=0; i<rec.rows(); ++i) {
            if (rec.getInt('pid',i) != 11) continue
            double px=rec.getFloat('px',i), py=rec.getFloat('py',i), pz=rec.getFloat('pz',i)
            double p=Math.sqrt(px*px+py*py+pz*pz)
            if (p>bestp) { bestp=p; best=i }
        }
        return best
    }

    static void main(String[] args) {
        if (args.length < 3) {
            System.err.println('Usage: run-groovy pid_exclusivity_skim.groovy input.hipo output.txt Su22|Fa22|Sp23')
            System.exit(2)
        }
        String input=args[0], output=args[1], period=args[2]
        if (!MX2_SIGMA.containsKey(period)) throw new IllegalArgumentException('Unknown period '+period)

        QADB qa=makeQA()
        GenericKinematicFitter fitter = new analysis_fitter(10.6041)
        BufferedWriter out = new BufferedWriter(new FileWriter(output))

        // Header is deliberately a comment; numpy/pandas can ignore it.
        out.write('# run event period assigned_pid charge status fd beta chi2pid p px py pz theta phi vz ' +
                  'e_p e_px e_py e_pz e_theta e_phi e_vz Ebeam Q2 W xB y t tmin minus_tprime Mx2 ' +
                  'analysis_bin mx2_mu mx2_sigma mx2_nsigma pass_qa pass_fd pass_p pass_chi2pid ' +
                  'pass_Q2 pass_W pass_y pass_phase_space pass_mx2 pass_final\n')

        HipoDataSource reader=new HipoDataSource()
        reader.open(input)
        long nev=0, ncand=0
        while (reader.hasEvent()) {
            DataEvent event=reader.getNextEvent(); ++nev
            if (!event.hasBank('RUN::config') || !event.hasBank('REC::Particle')) continue
            HipoDataBank runb=(HipoDataBank)event.getBank('RUN::config')
            HipoDataBank rec=(HipoDataBank)event.getBank('REC::Particle')
            int run=runb.getInt('run',0), ev=runb.getInt('event',0)
            if (rejectRun(run)) continue
            boolean passQA = qa.pass(run,ev)
            if (!passQA) continue

            // Use the same run-dependent beam-energy machinery as the production processor.
            PhysicsEvent research=fitter.getPhysicsEvent(event)
            if (research == null || research.countByPid(11) < 1) continue
            BeamEnergy ebObj = new BeamEnergy(research, run, false)
            double Ebeam=ebObj.Eb()

            int ie=electronIndex(rec)
            if (ie<0) continue
            double epx=rec.getFloat('px',ie), epy=rec.getFloat('py',ie), epz=rec.getFloat('pz',ie)
            double ep=Math.sqrt(epx*epx+epy*epy+epz*epz)
            double Ee=Math.sqrt(ep*ep+ME*ME)
            double etheta=Math.toDegrees(Math.acos(epz/ep))
            double ephi=Math.toDegrees(Math.atan2(epy,epx)); if (ephi<0) ephi+=360.0
            double evz=rec.getFloat('vz',ie)

            // q = k-k', target=(Mp,0). Beam electron mass is negligible here.
            double qE=Ebeam-Ee, qx=-epx, qy=-epy, qz=Ebeam-epz
            double Q2=-(qE*qE-qx*qx-qy*qy-qz*qz)
            double nu=qE
            double W2=MP*MP + 2.0*MP*nu - Q2
            double W=W2>0 ? Math.sqrt(W2) : Double.NaN
            double xB=(Q2>0 && nu>0) ? Q2/(2.0*MP*nu) : Double.NaN
            double y=nu/Ebeam
            double tmin=mesonicTmin(Q2,W)

            for (int ih=0; ih<rec.rows(); ++ih) {
                if (ih==ie) continue
                int charge=(int)rec.getByte('charge',ih)
                if (charge<=0) continue
                int pid=rec.getInt('pid',ih), status=rec.getInt('status',ih)
                boolean fd=isFD(status)

                double px=rec.getFloat('px',ih), py=rec.getFloat('py',ih), pz=rec.getFloat('pz',ih)
                double p=Math.sqrt(px*px+py*py+pz*pz)
                if (!(p>0)) continue
                double beta=rec.getFloat('beta',ih), chi2=rec.getFloat('chi2pid',ih), vz=rec.getFloat('vz',ih)
                double theta=Math.toDegrees(Math.acos(pz/p))
                double phi=Math.toDegrees(Math.atan2(py,px)); if (phi<0) phi+=360.0

                // CRITICAL: force every positive hadron candidate to the pion hypothesis.
                double Epi=Math.sqrt(p*p+MPI*MPI)

                // t = (q-p_pi)^2, i.e. the mesonic definition used for pi+ production.
                double dE=qE-Epi, dx=qx-px, dy=qy-py, dz=qz-pz
                double t=mass2(dE,dx,dy,dz)
                double minusTp=Double.isFinite(tmin) ? -(t-tmin) : Double.NaN

                // Mx2 = (P+q-p_pi)^2 under the SAME pion hypothesis.
                double mxE=MP+qE-Epi, mxpx=qx-px, mxpy=qy-py, mxpz=qz-pz
                double Mx2=mass2(mxE,mxpx,mxpy,mxpz)

                int abin=analysisBin(xB,minusTp)
                double mu=abin>0 ? MX2_MU[abin-1] : Double.NaN
                double sig=abin>0 ? MX2_SIGMA[period][abin-1] : Double.NaN
                double nsig=(abin>0 && sig>0) ? (Mx2-mu)/sig : Double.NaN

                boolean passP=(p>0.5 && p<5.0)
                boolean passChi=(Math.abs(chi2)<3.5)
                boolean passQ2=(Q2>1.0)
                boolean passW=(W>2.0)
                boolean passY=(y<0.80)
                boolean passPhase=(abin>0)
                boolean passMx=(abin>0 && Math.abs(nsig)<2.0)
                boolean passFinal=fd && passP && passChi && passQ2 && passW && passY && passPhase && passMx

                out.write(String.format(java.util.Locale.US,
                    '%d %d %s %d %d %d %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ' +
                    '%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g ' +
                    '%d %.8g %.8g %.8g %d %d %d %d %d %d %d %d %d %d%n',
                    run,ev,period,pid,charge,status,fd?1:0,beta,chi2,p,px,py,pz,theta,phi,vz,
                    ep,epx,epy,epz,etheta,ephi,evz,Ebeam,Q2,W,xB,y,t,tmin,minusTp,Mx2,
                    abin,mu,sig,nsig,1,fd?1:0,passP?1:0,passChi?1:0,passQ2?1:0,passW?1:0,
                    passY?1:0,passPhase?1:0,passMx?1:0,passFinal?1:0))
                ++ncand
            }
            if (nev%250000==0) println("${new File(input).name}: events=${nev}, candidates=${ncand}")
        }
        reader.close(); out.close()
        println("DONE ${input}: events=${nev}, positive-FD/all-positive candidates written=${ncand} -> ${output}")
    }
}
PIDExclusivitySkim.main(args)
