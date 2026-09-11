#include <TChain.h>
#include <TSystem.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Stats {
    long long n = 0;
    double sum = 0.0;
    double sum2 = 0.0;
    double minv = std::numeric_limits<double>::infinity();
    double maxv = -std::numeric_limits<double>::infinity();

    void fill(double x) {
        if (!std::isfinite(x) || x < -900.0) return;
        ++n;
        sum += x;
        sum2 += x*x;
        if (x < minv) minv = x;
        if (x > maxv) maxv = x;
    }

    double mean() const { return n ? sum/static_cast<double>(n) : 0.0; }
    double rms() const {
        if (!n) return 0.0;
        const double m = mean();
        return std::sqrt(std::max(0.0, sum2/static_cast<double>(n) - m*m));
    }
};

std::string pct(long long num, long long den) {
    std::ostringstream ss;
    if (den <= 0) {
        ss << "n/a";
    } else {
        ss << std::fixed << std::setprecision(3)
           << 100.0*static_cast<double>(num)/static_cast<double>(den) << "%";
    }
    return ss.str();
}

void print_stats(std::ostream& os, const char* label, const Stats& s, const char* units="") {
    os << std::left << std::setw(28) << label;
    if (!s.n) {
        os << " no valid entries\n";
        return;
    }
    os << " n=" << std::setw(10) << s.n
       << " mean=" << std::setw(12) << std::setprecision(6) << std::fixed << s.mean()
       << " RMS=" << std::setw(12) << s.rms()
       << " min=" << std::setw(12) << s.minv
       << " max=" << std::setw(12) << s.maxv;
    if (units && units[0] != '\0') os << " " << units;
    os << "\n";
}

bool has_branch(TChain& c, const char* name) {
    return c.GetBranch(name) != nullptr;
}

} // namespace

void photon_efficiency_diagnostics(const char* input,
                                   const char* output_txt="photon_efficiency_diagnostics.txt") {
    TString in(input);
    TString addPattern = in;

    // A directory may be supplied directly. Otherwise input can be one ROOT file
    // or a wildcard such as /path/nSidis_*_photon_efficiency.root.
    void* dirp = gSystem->OpenDirectory(in.Data());
    if (dirp) {
        gSystem->FreeDirectory(dirp);
        if (!in.EndsWith("/")) in += "/";
        addPattern = in + "*.root";
    }

    TChain chain("PhotonEfficiency");
    const int nFiles = chain.Add(addPattern.Data());
    const long long nEntries = chain.GetEntries();

    std::ofstream fout;
    if (output_txt && std::string(output_txt).size()) fout.open(output_txt);

    auto report = [&](const std::string& s) {
        std::cout << s;
        if (fout.is_open()) fout << s;
    };

    std::ostringstream head;
    head << "============================================================\n"
         << "Photon-efficiency skim diagnostic\n"
         << "Input: " << input << "\n"
         << "TChain pattern: " << addPattern.Data() << "\n"
         << "Files added: " << nFiles << "\n"
         << "Tree entries (tag hypotheses): " << nEntries << "\n"
         << "============================================================\n";
    report(head.str());

    if (nFiles <= 0 || nEntries <= 0) {
        report("No nonempty PhotonEfficiency entries found.\n");
        if (fout.is_open()) fout.close();
        return;
    }

    // Scalars.
    Int_t runnum=0, evnum=0, tag_detector=-1, tag_pass_beta=0, tag_pass_fiducial=0;
    Int_t n_neutral=0, n_pid22=0;
    Double_t Mx2_ep=0, Mx_ep=0, Mx2_epg_raw=0, Mx2_epg_corr=0;
    Double_t tag_corr_p=0, tag_corr_theta=0;
    Double_t probe_corr_E=0, probe_corr_p=0, probe_corr_theta=0, probe_corr_phi=0;

    // Nearest reconstructed candidates.
    Int_t neutral_idx[5]      = {-999,-999,-999,-999,-999};
    Int_t neutral_pid[5]      = {-999,-999,-999,-999,-999};
    Int_t neutral_charge[5]   = {-999,-999,-999,-999,-999};
    Int_t neutral_detector[5] = {-999,-999,-999,-999,-999};
    Double_t neutral_p[5]     = {-999,-999,-999,-999,-999};
    Double_t neutral_delta_alpha[5] = {-999,-999,-999,-999,-999};

    Int_t any_idx[3]      = {-999,-999,-999};
    Int_t any_pid[3]      = {-999,-999,-999};
    Int_t any_charge[3]   = {-999,-999,-999};
    Int_t any_detector[3] = {-999,-999,-999};
    Double_t any_delta_alpha[3] = {-999,-999,-999};

    // Check branches before doing anything so a converter/version mismatch is obvious.
    const std::vector<const char*> required = {
        "runnum","evnum","tag_detector","tag_pass_beta","tag_pass_fiducial",
        "n_neutral","n_pid22","Mx2_ep","Mx_ep","Mx2_epg_raw","Mx2_epg_corr",
        "tag_corr_p","tag_corr_theta","probe_corr_E","probe_corr_p",
        "probe_corr_theta","probe_corr_phi","neutral_idx","neutral_pid",
        "neutral_charge","neutral_detector","neutral_p","neutral_delta_alpha",
        "any_idx","any_pid","any_charge","any_detector","any_delta_alpha"
    };
    bool missing = false;
    for (const char* b : required) {
        if (!has_branch(chain,b)) {
            report(std::string("ERROR: missing branch: ") + b + "\n");
            missing = true;
        }
    }
    if (missing) {
        report("Stopping because the tree does not match this diagnostic macro.\n");
        if (fout.is_open()) fout.close();
        return;
    }

    chain.SetBranchAddress("runnum",&runnum);
    chain.SetBranchAddress("evnum",&evnum);
    chain.SetBranchAddress("tag_detector",&tag_detector);
    chain.SetBranchAddress("tag_pass_beta",&tag_pass_beta);
    chain.SetBranchAddress("tag_pass_fiducial",&tag_pass_fiducial);
    chain.SetBranchAddress("n_neutral",&n_neutral);
    chain.SetBranchAddress("n_pid22",&n_pid22);
    chain.SetBranchAddress("Mx2_ep",&Mx2_ep);
    chain.SetBranchAddress("Mx_ep",&Mx_ep);
    chain.SetBranchAddress("Mx2_epg_raw",&Mx2_epg_raw);
    chain.SetBranchAddress("Mx2_epg_corr",&Mx2_epg_corr);
    chain.SetBranchAddress("tag_corr_p",&tag_corr_p);
    chain.SetBranchAddress("tag_corr_theta",&tag_corr_theta);
    chain.SetBranchAddress("probe_corr_E",&probe_corr_E);
    chain.SetBranchAddress("probe_corr_p",&probe_corr_p);
    chain.SetBranchAddress("probe_corr_theta",&probe_corr_theta);
    chain.SetBranchAddress("probe_corr_phi",&probe_corr_phi);
    chain.SetBranchAddress("neutral_idx",neutral_idx);
    chain.SetBranchAddress("neutral_pid",neutral_pid);
    chain.SetBranchAddress("neutral_charge",neutral_charge);
    chain.SetBranchAddress("neutral_detector",neutral_detector);
    chain.SetBranchAddress("neutral_p",neutral_p);
    chain.SetBranchAddress("neutral_delta_alpha",neutral_delta_alpha);
    chain.SetBranchAddress("any_idx",any_idx);
    chain.SetBranchAddress("any_pid",any_pid);
    chain.SetBranchAddress("any_charge",any_charge);
    chain.SetBranchAddress("any_detector",any_detector);
    chain.SetBranchAddress("any_delta_alpha",any_delta_alpha);

    long long tagFT=0, tagFD=0, tagCD=0, tagUnknown=0;
    long long tagBeta=0, tagFid=0, tagBetaFid=0;
    long long probeFTangle=0, probeFDangle=0, probeOtherAngle=0;
    long long probeE04=0, probeE10=0, probeE20=0, probeE40=0, probeE60=0;
    long long ftTag_ftProbe=0, ftTag_fdProbe=0, fdTag_ftProbe=0, fdTag_fdProbe=0;

    long long nearestNeutralValid=0, nearestPid22=0, nearestPid22P04=0;
    long long nearestAnyCharged=0, nearestAnyElectron=0;

    const double radii[] = {0.5,1.0,2.0,3.0,5.0,10.0};
    const int NR = sizeof(radii)/sizeof(radii[0]);
    long long neutralWithin[NR] = {0};
    long long pid22Within[NR] = {0};
    long long pid22P04Within[NR] = {0};
    long long chargedWithin[NR] = {0};
    long long electronWithin[NR] = {0};

    // Broad diagnostic pi0 windows. These are NOT final analysis cuts.
    long long pi0Broad=0, pi0Core=0;
    long long pi0CoreFTProbe=0, pi0CoreFDProbe=0;
    long long pi0CorePid22Within2=0, pi0CorePid22Within3=0, pi0CorePid22P04Within3=0;
    long long pi0CoreFTPid22Within3=0, pi0CoreFDPid22Within3=0;

    Stats sMx2ep, sMxep, sMx2epgRaw, sMx2epgCorr;
    Stats sTagP, sTagTheta, sProbeE, sProbeP, sProbeTheta;
    Stats sNearNeutralDA, sNearPid22DA, sNearAnyDA;

    std::map<int,long long> runRows;
    std::set<unsigned long long> uniqueEvents;

    for (long long i=0; i<nEntries; ++i) {
        chain.GetEntry(i);

        runRows[runnum]++;
        const unsigned long long key = (static_cast<unsigned long long>(static_cast<unsigned int>(runnum)) << 32)
                                     | static_cast<unsigned int>(evnum);
        uniqueEvents.insert(key);

        if (tag_detector==0) ++tagFT;
        else if (tag_detector==1) ++tagFD;
        else if (tag_detector==2) ++tagCD;
        else ++tagUnknown;

        if (tag_pass_beta) ++tagBeta;
        if (tag_pass_fiducial) ++tagFid;
        if (tag_pass_beta && tag_pass_fiducial) ++tagBetaFid;

        const bool probeFT = (probe_corr_theta >= 0.0 && probe_corr_theta <= 5.5);
        const bool probeFD = (probe_corr_theta > 5.5 && probe_corr_theta <= 35.0);
        if (probeFT) ++probeFTangle;
        else if (probeFD) ++probeFDangle;
        else ++probeOtherAngle;

        if (tag_detector==0 && probeFT) ++ftTag_ftProbe;
        if (tag_detector==0 && probeFD) ++ftTag_fdProbe;
        if (tag_detector==1 && probeFT) ++fdTag_ftProbe;
        if (tag_detector==1 && probeFD) ++fdTag_fdProbe;

        if (probe_corr_E >= 0.4) ++probeE04;
        if (probe_corr_E >= 1.0) ++probeE10;
        if (probe_corr_E >= 2.0) ++probeE20;
        if (probe_corr_E >= 4.0) ++probeE40;
        if (probe_corr_E >= 6.0) ++probeE60;

        sMx2ep.fill(Mx2_ep);
        sMxep.fill(Mx_ep);
        sMx2epgRaw.fill(Mx2_epg_raw);
        sMx2epgCorr.fill(Mx2_epg_corr);
        sTagP.fill(tag_corr_p);
        sTagTheta.fill(tag_corr_theta);
        sProbeE.fill(probe_corr_E);
        sProbeP.fill(probe_corr_p);
        sProbeTheta.fill(probe_corr_theta);

        const bool nvalid = neutral_idx[0] >= 0 && neutral_delta_alpha[0] >= 0.0;
        const bool npid22 = nvalid && neutral_pid[0] == 22;
        const bool npid22p04 = npid22 && neutral_p[0] >= 0.4;
        if (nvalid) {
            ++nearestNeutralValid;
            sNearNeutralDA.fill(neutral_delta_alpha[0]);
        }
        if (npid22) {
            ++nearestPid22;
            sNearPid22DA.fill(neutral_delta_alpha[0]);
        }
        if (npid22p04) ++nearestPid22P04;

        const bool avalid = any_idx[0] >= 0 && any_delta_alpha[0] >= 0.0;
        if (avalid) sNearAnyDA.fill(any_delta_alpha[0]);
        if (avalid && any_charge[0] != 0) ++nearestAnyCharged;
        if (avalid && std::abs(any_pid[0]) == 11) ++nearestAnyElectron;

        for (int ir=0; ir<NR; ++ir) {
            if (nvalid && neutral_delta_alpha[0] < radii[ir]) ++neutralWithin[ir];
            if (npid22 && neutral_delta_alpha[0] < radii[ir]) ++pid22Within[ir];
            if (npid22p04 && neutral_delta_alpha[0] < radii[ir]) ++pid22P04Within[ir];
            if (avalid && any_charge[0] != 0 && any_delta_alpha[0] < radii[ir]) ++chargedWithin[ir];
            if (avalid && std::abs(any_pid[0]) == 11 && any_delta_alpha[0] < radii[ir]) ++electronWithin[ir];
        }

        const bool broad = (Mx_ep >= 0.08 && Mx_ep <= 0.20);
        const bool core  = (Mx_ep >= 0.10 && Mx_ep <= 0.17);
        if (broad) ++pi0Broad;
        if (core) {
            ++pi0Core;
            if (probeFT) ++pi0CoreFTProbe;
            if (probeFD) ++pi0CoreFDProbe;
            if (npid22 && neutral_delta_alpha[0] < 2.0) ++pi0CorePid22Within2;
            if (npid22 && neutral_delta_alpha[0] < 3.0) ++pi0CorePid22Within3;
            if (npid22p04 && neutral_delta_alpha[0] < 3.0) ++pi0CorePid22P04Within3;
            if (probeFT && npid22 && neutral_delta_alpha[0] < 3.0) ++pi0CoreFTPid22Within3;
            if (probeFD && npid22 && neutral_delta_alpha[0] < 3.0) ++pi0CoreFDPid22Within3;
        }
    }

    std::ostringstream os;
    os << std::fixed << std::setprecision(3);
    os << "\n--- Run / row accounting ------------------------------------\n";
    os << "Unique events represented: " << uniqueEvents.size() << "\n";
    os << "Mean tag hypotheses / represented event: "
       << (uniqueEvents.empty() ? 0.0 : static_cast<double>(nEntries)/uniqueEvents.size()) << "\n";
    for (const auto& kv : runRows) os << "  run " << kv.first << " : " << kv.second << " rows\n";

    os << "\n--- Actual reconstructed TAG detector -----------------------\n";
    os << "FT (detector=0): " << tagFT << "  (" << pct(tagFT,nEntries) << ")\n";
    os << "FD (detector=1): " << tagFD << "  (" << pct(tagFD,nEntries) << ")\n";
    os << "CD (detector=2): " << tagCD << "  (" << pct(tagCD,nEntries) << ")\n";
    os << "unknown:         " << tagUnknown << "  (" << pct(tagUnknown,nEntries) << ")\n";
    os << "tag beta pass:   " << tagBeta << "  (" << pct(tagBeta,nEntries) << ")\n";
    os << "tag fid pass:    " << tagFid << "  (" << pct(tagFid,nEntries) << ")\n";
    os << "both beta+fid:   " << tagBetaFid << "  (" << pct(tagBetaFid,nEntries) << ")\n";

    os << "\n--- Predicted PROBE angular region ---------------------------\n";
    os << "FT-angle theta <= 5.5 deg:       " << probeFTangle << "  (" << pct(probeFTangle,nEntries) << ")\n";
    os << "FD-angle 5.5 < theta <= 35 deg:  " << probeFDangle << "  (" << pct(probeFDangle,nEntries) << ")\n";
    os << "outside those angular regions:   " << probeOtherAngle << "  (" << pct(probeOtherAngle,nEntries) << ")\n";
    os << "NOTE: these probe labels use theta only; they are NOT probe fiducial cuts.\n";

    os << "\nTag/probe topology (probe defined by angle only):\n";
    os << "  FT tag + FT-angle probe : " << ftTag_ftProbe << "\n";
    os << "  FT tag + FD-angle probe : " << ftTag_fdProbe << "\n";
    os << "  FD tag + FT-angle probe : " << fdTag_ftProbe << "\n";
    os << "  FD tag + FD-angle probe : " << fdTag_fdProbe << "\n";

    os << "\n--- Predicted probe energy -----------------------------------\n";
    os << "E >= 0.4 GeV : " << probeE04 << "  (" << pct(probeE04,nEntries) << ")\n";
    os << "E >= 1.0 GeV : " << probeE10 << "  (" << pct(probeE10,nEntries) << ")\n";
    os << "E >= 2.0 GeV : " << probeE20 << "  (" << pct(probeE20,nEntries) << ")\n";
    os << "E >= 4.0 GeV : " << probeE40 << "  (" << pct(probeE40,nEntries) << ")\n";
    os << "E >= 6.0 GeV : " << probeE60 << "  (" << pct(probeE60,nEntries) << ")\n";

    os << "\n--- Nearest reconstructed candidate --------------------------\n";
    os << "Nearest neutral exists:       " << nearestNeutralValid << "  (" << pct(nearestNeutralValid,nEntries) << ")\n";
    os << "Nearest neutral is pid 22:    " << nearestPid22 << "  (" << pct(nearestPid22,nEntries) << ")\n";
    os << "Nearest pid22 and p>=0.4:     " << nearestPid22P04 << "  (" << pct(nearestPid22P04,nEntries) << ")\n";
    os << "Nearest ANY is charged:       " << nearestAnyCharged << "  (" << pct(nearestAnyCharged,nEntries) << ")\n";
    os << "Nearest ANY has |pid|=11:     " << nearestAnyElectron << "  (" << pct(nearestAnyElectron,nEntries) << ")\n";
    os << "\nOpening-angle matching scan:\n";
    os << "  dAlpha     neutral(any PID)      pid22             pid22,p>=0.4       charged(any)       |pid|=11\n";
    for (int ir=0; ir<NR; ++ir) {
        os << "  <" << std::setw(4) << radii[ir] << " deg  "
           << std::setw(10) << neutralWithin[ir] << " (" << std::setw(8) << pct(neutralWithin[ir],nEntries) << ")  "
           << std::setw(10) << pid22Within[ir] << " (" << std::setw(8) << pct(pid22Within[ir],nEntries) << ")  "
           << std::setw(10) << pid22P04Within[ir] << " (" << std::setw(8) << pct(pid22P04Within[ir],nEntries) << ")  "
           << std::setw(10) << chargedWithin[ir] << " (" << std::setw(8) << pct(chargedWithin[ir],nEntries) << ")  "
           << std::setw(10) << electronWithin[ir] << " (" << std::setw(8) << pct(electronWithin[ir],nEntries) << ")\n";
    }

    os << "\n--- Broad pi0-mass diagnostic (NOT final signal extraction) --\n";
    os << "0.08 <= Mx(ep) <= 0.20 GeV : " << pi0Broad << "  (" << pct(pi0Broad,nEntries) << ")\n";
    os << "0.10 <= Mx(ep) <= 0.17 GeV : " << pi0Core  << "  (" << pct(pi0Core,nEntries) << ")\n";
    os << "  core + FT-angle probe      : " << pi0CoreFTProbe << "  (" << pct(pi0CoreFTProbe,pi0Core) << " of core)\n";
    os << "  core + FD-angle probe      : " << pi0CoreFDProbe << "  (" << pct(pi0CoreFDProbe,pi0Core) << " of core)\n";
    os << "  core + pid22 dA<2 deg      : " << pi0CorePid22Within2 << "  (" << pct(pi0CorePid22Within2,pi0Core) << " of core)\n";
    os << "  core + pid22 dA<3 deg      : " << pi0CorePid22Within3 << "  (" << pct(pi0CorePid22Within3,pi0Core) << " of core)\n";
    os << "  core + pid22 p>=0.4,dA<3   : " << pi0CorePid22P04Within3 << "  (" << pct(pi0CorePid22P04Within3,pi0Core) << " of core)\n";
    os << "  FT-angle core recovered*   : " << pi0CoreFTPid22Within3 << " / " << pi0CoreFTProbe
       << "  (" << pct(pi0CoreFTPid22Within3,pi0CoreFTProbe) << ")\n";
    os << "  FD-angle core recovered*   : " << pi0CoreFDPid22Within3 << " / " << pi0CoreFDProbe
       << "  (" << pct(pi0CoreFDPid22Within3,pi0CoreFDProbe) << ")\n";
    os << "  *raw diagnostic only: nearest neutral pid22 within 3 deg; no background subtraction/final fiducial definition.\n";

    report(os.str());

    std::ostringstream statsOut;
    statsOut << "\n--- Kinematic summaries --------------------------------------\n";
    report(statsOut.str());

    // print_stats writes directly, so duplicate into the text output when requested.
    auto both_stats = [&](const char* label, const Stats& s, const char* units) {
        std::ostringstream tmp;
        print_stats(tmp,label,s,units);
        report(tmp.str());
    };
    both_stats("Mx2(ep)",sMx2ep,"GeV^2");
    both_stats("Mx(ep)",sMxep,"GeV");
    both_stats("Mx2(epg) raw",sMx2epgRaw,"GeV^2");
    both_stats("Mx2(epg) corrected",sMx2epgCorr,"GeV^2");
    both_stats("tag corrected p",sTagP,"GeV");
    both_stats("tag corrected theta",sTagTheta,"deg");
    both_stats("probe corrected E",sProbeE,"GeV");
    both_stats("probe corrected p",sProbeP,"GeV");
    both_stats("probe corrected theta",sProbeTheta,"deg");
    both_stats("nearest neutral dAlpha",sNearNeutralDA,"deg");
    both_stats("nearest pid22 dAlpha",sNearPid22DA,"deg");
    both_stats("nearest any dAlpha",sNearAnyDA,"deg");

    report("\nInterpretation note: this macro is intentionally diagnostic. The broad pi0 windows and angular FT/FD probe labels are not final efficiency selections.\n");

    if (fout.is_open()) {
        fout.close();
        std::cout << "\nSaved text report to: " << output_txt << "\n";
    }
}
