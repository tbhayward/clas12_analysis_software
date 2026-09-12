// photon_efficiency_valerii_reproduction.C
//
// Photon tag-and-probe analysis with two layers:
//   (1) stage-1 integrated diagnostics;
//   (2) a Valerii-style FD 7 x 3 x 6 efficiency/correction reproduction.
//
// The stage-1 layer remains useful QA; the Valerii layer is the physics path.
// NOTE: the June-2026 Delta-phi(p,gamma) exclusivity requirement is implemented
// as |delta phi_copl| < 5.7 deg using the zero-centered Trento-plane residual.
// Raw lab and unshifted Trento observables are retained as QA diagnostics.
//
//   * fully integrated FD extraction;
//   * fully integrated FT extraction;
//   * FT split into two broad energy bins: 0.4-2 GeV and >=2 GeV;
//   * 1-, 2-, and 3-sigma Delta-p raw recovery fractions;
//   * cut-flow accounting;
//   * pi0-parent-window stability scan;
//   * denominator p/theta/phi QA;
//   * automatic use of whatever AAOgen / CLASDIS / DVCSgen ROOT files exist;
//   * process-level parallel analysis of data / AAOgen / CLASDIS / DVCSgen;
//   * optional MC-truth closure diagnostics when truth branches are present.
//
// Every invocation deletes and rebuilds ./output.
//
// Run:
//   root -l -b -q 'photon_efficiency_valerii_reproduction.C()'
//
// The data/MC ROOT trees may keep accumulating between invocations; each new
// invocation chains every completed *.root file that exists at that moment.

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>
#include <TParameter.h>
#include <TNamed.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <numeric>

namespace pe {

void concise_publication_style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);
    gStyle->SetTitleSize(0.050,"XYZ");
    gStyle->SetLabelSize(0.042,"XYZ");
    gStyle->SetTitleOffset(1.05,"X");
    gStyle->SetTitleOffset(1.25,"Y");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(2);
    gStyle->SetEndErrorSize(4);
    gStyle->SetLegendBorderSize(0);
}

// -----------------------------------------------------------------------------
// Paths.
// -----------------------------------------------------------------------------
static const char* DATA_DIR =
    "/work/clas12/thayward/photon_efficiency/ROOT_trees/data/fa18_inb";
static const char* AAOGEN_DIR =
    "/work/clas12/thayward/photon_efficiency/ROOT_trees/aaogen/fa18_inb";
static const char* CLASDIS_DIR =
    "/work/clas12/thayward/photon_efficiency/ROOT_trees/clasdis/fa18_inb";
static const char* DVCSGEN_DIR =
    "/work/clas12/thayward/photon_efficiency/ROOT_trees/dvcsgen/fa18_inb";

// -----------------------------------------------------------------------------
// Stage-1 selection.
// -----------------------------------------------------------------------------

// Current nominal pi0-parent enrichment window.
// This remains explicitly a development choice until the exact production
// denominator selection is frozen.
static const double PI0_M_MIN = 0.08;
static const double PI0_M_MAX = 0.20;

static const double PROBE_P_MIN = 0.40;

// FD analysis range follows the Valerii analysis angular/momentum region.
static const double FD_THETA_MIN = 6.0;
static const double FD_THETA_MAX = 36.0;
static const double FD_P_MIN = 0.35;
static const double FD_P_MAX = 6.0;

// FT extension.  The main stage-1 physics question is energy dependence.
// Keep one integrated result plus two broad energy bins.
static const double FT_THETA_MIN = 2.0;
static const double FT_THETA_MAX = 5.5;
static const double FT_P_MIN = 0.40;
static const double FT_P_SPLIT = 2.0;
static const double FT_P_MAX = 8.5;

// FT nominal geometry from the existing CLAS12 fiducial implementation.
// The inferred probe is projected to the FT response plane and MUST land in
// the usable annulus and outside the excluded holes before it enters an FT
// denominator.  This is the FT analogue of applying the detector fiducial
// acceptance to the expected probe.
static const double FT_R_MIN = 8.5;   // cm
static const double FT_R_MAX = 15.5;  // cm

struct FTHole {
    double r, x, y;
};
static const FTHole FT_HOLES[] = {
    {1.60, -8.42,   9.89},
    {1.60, -9.89,  -5.33},
    {2.30, -6.15, -13.00},
    {2.00,  3.70,  -6.50}
};

// Residual histogram and fit.
static const int    DP_NBIN = 100;
static const double DP_HMIN = -1.0;
static const double DP_HMAX =  1.0;
static const double FIT_MIN = -0.60;
static const double FIT_MAX =  0.60;

static const int MIN_FIT_ENTRIES_FD = 100;
static const int MIN_FIT_ENTRIES_FT = 75;

static const double SIGMA_MIN = 0.010;
static const double SIGMA_MAX = 0.400;
static const double MAX_REL_SIGMA_ERR = 0.50;
static const double MAX_MEAN_ERR = 0.100;
static const double FIT_MEAN_MIN = -0.50;
static const double FIT_MEAN_MAX =  0.50;
static const double MEAN_BOUNDARY_TOL = 0.010;

// Parent-window scan.  Only a compact CSV is written.
struct MassWindow {
    const char* label;
    double lo, hi;
};
static const MassWindow MASS_WINDOWS[] = {
    {"tight_0p10_0p17", 0.10, 0.17},
    {"nominal_0p08_0p20", 0.08, 0.20},
    {"loose_0p05_0p25", 0.05, 0.25},
    {"very_loose_0p00_0p30", 0.00, 0.30}
};

// -----------------------------------------------------------------------------
// Utilities.
// -----------------------------------------------------------------------------
bool finite_good(double x) {
    return std::isfinite(x) && x > -900.0;
}

double deg2rad(double x) { return x*TMath::Pi()/180.0; }

double opening_angle_deg(double th1, double ph1, double th2, double ph2) {
    if (!finite_good(th1) || !finite_good(ph1) ||
        !finite_good(th2) || !finite_good(ph2)) return -999.0;

    const double t1=deg2rad(th1), p1=deg2rad(ph1);
    const double t2=deg2rad(th2), p2=deg2rad(ph2);

    double c=std::sin(t1)*std::sin(t2)*std::cos(p1-p2)
            +std::cos(t1)*std::cos(t2);
    c=std::max(-1.0,std::min(1.0,c));
    return std::acos(c)*180.0/TMath::Pi();
}

std::string make_pattern(const std::string& input) {
    TString in(input.c_str());
    if (in.Length()==0) return "";

    void* d=gSystem->OpenDirectory(in.Data());
    if (d) {
        gSystem->FreeDirectory(d);
        if (!in.EndsWith("/")) in += "/";
        in += "*.root";
    }
    return in.Data();
}

bool has_root_files(const std::string& dir) {
    void* dp=gSystem->OpenDirectory(dir.c_str());
    if (!dp) return false;

    bool found=false;
    const char* entry=nullptr;
    while ((entry=gSystem->GetDirEntry(dp))) {
        std::string n(entry);
        if (n.size()>=5 && n.substr(n.size()-5)==".root") {
            found=true;
            break;
        }
    } // endwhile
    gSystem->FreeDirectory(dp);
    return found;
}

void reset_output(const std::string& out) {
    if (gSystem->AccessPathName(out.c_str())==kFALSE)
        gSystem->Exec(("rm -rf "+out).c_str());
    gSystem->mkdir(out.c_str(),true);
}

double wrap_phi(double phi) {
    while (phi < -30.0) phi += 360.0;
    while (phi >= 330.0) phi -= 360.0;
    return phi;
}

// -----------------------------------------------------------------------------
// Tree branches.
// -----------------------------------------------------------------------------
struct Branches {
    Int_t runnum=0, evnum=0, is_mc=0;
    Int_t p_pass_standard=0;
    Int_t tag_detector=-1, tag_pass_beta=0, tag_pass_fiducial=0;
    Double_t mc_weight=1.0;

    Double_t Mx_ep=0, Mx2_ep=0;
    bool have_Mx2_ep=false;
    Double_t probe_corr_p=0, probe_corr_theta=0, probe_corr_phi=0;

    // Additional skim quantities used by the data-driven MC normalization.
    // These are optional for the stage-1 QA, but the normalization stage will
    // explicitly report which observables are unavailable if a legacy skim is used.
    Double_t beam_energy=0;
    Double_t e_p=0, e_theta=0, e_phi=0;
    Double_t p_corr_p=0, p_corr_theta=0, p_corr_phi=0;
    Double_t tag_corr_p=0, tag_corr_theta=0, tag_corr_phi=0;
    Double_t Mx2_epg_corr=0;
    Double_t probe_corr_E=0;
    bool have_beam_energy=false, have_e_kin=false, have_p_corr_kin=false;
    bool have_tag_corr_kin=false, have_Mx2_epg_corr=false, have_probe_corr_E=false;

    // Optional event vertex, if available in the current converter.
    Double_t e_vz=0;
    bool have_e_vz=false;

    Int_t neutral_idx[5];
    Int_t neutral_pid[5];
    Int_t neutral_charge[5];
    Int_t neutral_detector[5];
    Double_t neutral_p[5];
    Double_t neutral_theta[5];
    Double_t neutral_phi[5];
    Double_t neutral_delta_alpha[5];
    Double_t neutral_x[5];
    Double_t neutral_y[5];
    Double_t neutral_z[5];

    // Optional MC truth.
    Int_t truth_probe_index=-999;
    Int_t truth_probe_pid=-999;
    Int_t truth_probe_parent=-999;
    Double_t truth_probe_p=-999;
    Double_t truth_probe_theta=-999;
    Double_t truth_probe_phi=-999;
    Double_t truth_probe_delta_alpha=-999;
    bool have_truth=false;

    void reset_arrays() {
        for (int i=0;i<5;i++) {
            neutral_idx[i]=-999;
            neutral_pid[i]=-999;
            neutral_charge[i]=-999;
            neutral_detector[i]=-999;
            neutral_p[i]=-999;
            neutral_theta[i]=-999;
            neutral_phi[i]=-999;
            neutral_delta_alpha[i]=-999;
            neutral_x[i]=-999;
            neutral_y[i]=-999;
            neutral_z[i]=-999;
        } // endfor
    }
};

bool attach(TChain& c, Branches& b) {
    const char* required[] = {
        "runnum","evnum","is_mc","mc_weight","p_pass_standard",
        "tag_detector","tag_pass_beta","tag_pass_fiducial",
        "Mx_ep","probe_corr_p","probe_corr_theta","probe_corr_phi",
        "neutral_idx","neutral_pid","neutral_charge","neutral_detector",
        "neutral_p","neutral_theta","neutral_phi","neutral_delta_alpha",
        "neutral_x","neutral_y","neutral_z"
    };

    bool ok=true;
    for (const char* n : required) {
        if (!c.GetBranch(n)) {
            std::cerr << "ERROR: missing required branch " << n << "\n";
            ok=false;
        }
    } // endfor
    if (!ok) return false;

    // Read only branches used by this macro.  This materially reduces I/O for
    // multi-million-entry samples.
    c.SetBranchStatus("*",0);
    for (const char* n : required) c.SetBranchStatus(n,1);
    if (c.GetBranch("e_vz")) c.SetBranchStatus("e_vz",1);

    // Optional normalization branches.  They are present in the current
    // photon-efficiency converter, but keeping them optional preserves
    // compatibility with older skim files used for stage-1 diagnostics.
    const char* norm_optional[] = {
        "beam_energy","e_p","e_theta","e_phi",
        "p_corr_p","p_corr_theta","p_corr_phi",
        "tag_corr_p","tag_corr_theta","tag_corr_phi",
        "Mx2_ep","Mx2_epg_corr","probe_corr_E"
    };
    for (const char* n : norm_optional) if (c.GetBranch(n)) c.SetBranchStatus(n,1);

    c.SetBranchAddress("runnum",&b.runnum);
    c.SetBranchAddress("evnum",&b.evnum);
    c.SetBranchAddress("is_mc",&b.is_mc);
    c.SetBranchAddress("mc_weight",&b.mc_weight);
    c.SetBranchAddress("p_pass_standard",&b.p_pass_standard);
    c.SetBranchAddress("tag_detector",&b.tag_detector);
    c.SetBranchAddress("tag_pass_beta",&b.tag_pass_beta);
    c.SetBranchAddress("tag_pass_fiducial",&b.tag_pass_fiducial);
    c.SetBranchAddress("Mx_ep",&b.Mx_ep);
    c.SetBranchAddress("probe_corr_p",&b.probe_corr_p);
    c.SetBranchAddress("probe_corr_theta",&b.probe_corr_theta);
    c.SetBranchAddress("probe_corr_phi",&b.probe_corr_phi);

    if (c.GetBranch("Mx2_ep")) { c.SetBranchAddress("Mx2_ep",&b.Mx2_ep); b.have_Mx2_ep=true; }
    if (c.GetBranch("beam_energy")) { c.SetBranchAddress("beam_energy",&b.beam_energy); b.have_beam_energy=true; }
    if (c.GetBranch("e_p") && c.GetBranch("e_theta") && c.GetBranch("e_phi")) {
        c.SetBranchAddress("e_p",&b.e_p); c.SetBranchAddress("e_theta",&b.e_theta);
        c.SetBranchAddress("e_phi",&b.e_phi); b.have_e_kin=true;
    }
    if (c.GetBranch("p_corr_p") && c.GetBranch("p_corr_theta") && c.GetBranch("p_corr_phi")) {
        c.SetBranchAddress("p_corr_p",&b.p_corr_p); c.SetBranchAddress("p_corr_theta",&b.p_corr_theta);
        c.SetBranchAddress("p_corr_phi",&b.p_corr_phi); b.have_p_corr_kin=true;
    }
    if (c.GetBranch("tag_corr_p") && c.GetBranch("tag_corr_theta") && c.GetBranch("tag_corr_phi")) {
        c.SetBranchAddress("tag_corr_p",&b.tag_corr_p); c.SetBranchAddress("tag_corr_theta",&b.tag_corr_theta);
        c.SetBranchAddress("tag_corr_phi",&b.tag_corr_phi); b.have_tag_corr_kin=true;
    }
    if (c.GetBranch("Mx2_epg_corr")) { c.SetBranchAddress("Mx2_epg_corr",&b.Mx2_epg_corr); b.have_Mx2_epg_corr=true; }
    if (c.GetBranch("probe_corr_E")) { c.SetBranchAddress("probe_corr_E",&b.probe_corr_E); b.have_probe_corr_E=true; }

    c.SetBranchAddress("neutral_idx",b.neutral_idx);
    c.SetBranchAddress("neutral_pid",b.neutral_pid);
    c.SetBranchAddress("neutral_charge",b.neutral_charge);
    c.SetBranchAddress("neutral_detector",b.neutral_detector);
    c.SetBranchAddress("neutral_p",b.neutral_p);
    c.SetBranchAddress("neutral_theta",b.neutral_theta);
    c.SetBranchAddress("neutral_phi",b.neutral_phi);
    c.SetBranchAddress("neutral_delta_alpha",b.neutral_delta_alpha);
    c.SetBranchAddress("neutral_x",b.neutral_x);
    c.SetBranchAddress("neutral_y",b.neutral_y);
    c.SetBranchAddress("neutral_z",b.neutral_z);

    // Vertex branch name in the current converter.
    if (c.GetBranch("e_vz")) {
        c.SetBranchAddress("e_vz",&b.e_vz);
        b.have_e_vz=true;
    }

    // Optional truth branches: automatically active for MC skims.
    const char* truth_required[] = {
        "truth_probe_index","truth_probe_pid","truth_probe_parent",
        "truth_probe_p","truth_probe_theta","truth_probe_phi",
        "truth_probe_delta_alpha"
    };

    bool truth_ok=true;
    for (const char* n : truth_required)
        if (!c.GetBranch(n)) truth_ok=false;

    if (truth_ok) {
        for (const char* n : truth_required) c.SetBranchStatus(n,1);
        b.have_truth=true;
        c.SetBranchAddress("truth_probe_index",&b.truth_probe_index);
        c.SetBranchAddress("truth_probe_pid",&b.truth_probe_pid);
        c.SetBranchAddress("truth_probe_parent",&b.truth_probe_parent);
        c.SetBranchAddress("truth_probe_p",&b.truth_probe_p);
        c.SetBranchAddress("truth_probe_theta",&b.truth_probe_theta);
        c.SetBranchAddress("truth_probe_phi",&b.truth_probe_phi);
        c.SetBranchAddress("truth_probe_delta_alpha",&b.truth_probe_delta_alpha);
    }

    return true;
}

// -----------------------------------------------------------------------------
// Selection and candidate matching.
// -----------------------------------------------------------------------------
bool pass_before_mass(const Branches& b) {
    if (!b.p_pass_standard) return false;
    if (!b.tag_pass_beta) return false;
    if (!b.tag_pass_fiducial) return false;
    if (!(b.tag_detector==0 || b.tag_detector==1)) return false;
    if (!finite_good(b.probe_corr_p) || !finite_good(b.probe_corr_theta) ||
        !finite_good(b.probe_corr_phi)) return false;
    return true;
}

bool pass_nominal_mass(const Branches& b) {
    return finite_good(b.Mx_ep) &&
           b.Mx_ep>=PI0_M_MIN && b.Mx_ep<=PI0_M_MAX;
}

bool in_fd(const Branches& b) {
    return b.probe_corr_p>=FD_P_MIN && b.probe_corr_p<=FD_P_MAX &&
           b.probe_corr_theta>=FD_THETA_MIN &&
           b.probe_corr_theta<=FD_THETA_MAX;
}

bool in_ft(const Branches& b) {
    return b.probe_corr_p>=FT_P_MIN && b.probe_corr_p<=FT_P_MAX &&
           b.probe_corr_theta>=FT_THETA_MIN &&
           b.probe_corr_theta<=FT_THETA_MAX;
}

int best_probe_candidate(const Branches& b, int detector) {
    int best=-1;
    double best_da=1e9;

    for (int k=0;k<5;k++) {
        if (b.neutral_idx[k]<0) continue;
        if (b.neutral_charge[k]!=0) continue;
        if (b.neutral_pid[k]!=22) continue;
        if (b.neutral_detector[k]!=detector) continue;
        if (!finite_good(b.neutral_p[k]) || b.neutral_p[k]<PROBE_P_MIN) continue;
        if (!finite_good(b.neutral_delta_alpha[k])) continue;

        if (b.neutral_delta_alpha[k]<best_da) {
            best_da=b.neutral_delta_alpha[k];
            best=k;
        }
    } // endfor
    return best;
}


// -----------------------------------------------------------------------------
// FT-face projection preparation.
// -----------------------------------------------------------------------------
struct FTPlaneEstimate {
    bool valid=false;
    double z=0;
    long long n=0;
};

// Estimate the common FT response-plane z from reconstructed FT neutral/tag
// responses already present in the skim.  This avoids hard-coding a z-plane.
FTPlaneEstimate estimate_ft_plane(TChain& c, Branches& b) {
    std::vector<double> zs;
    zs.reserve(100000);

    for (Long64_t i=0;i<c.GetEntries();i++) {
        c.GetEntry(i);
        for (int k=0;k<5;k++) {
            if (b.neutral_idx[k]<0) continue;
            if (b.neutral_detector[k]!=0) continue;
            if (!finite_good(b.neutral_z[k])) continue;
            zs.push_back(b.neutral_z[k]);
            if (zs.size()>=100000) break;
        } // endfor
        if (zs.size()>=100000) break;
    } // endfor

    FTPlaneEstimate r;
    r.n=zs.size();
    if (zs.size()<10) return r;

    std::nth_element(zs.begin(),zs.begin()+zs.size()/2,zs.end());
    r.z=zs[zs.size()/2];
    r.valid=std::isfinite(r.z);
    return r;
}

struct FTProjection {
    bool valid=false;
    double x=-999,y=-999,r=-999;
    bool radial=false;
    bool hole_clear=false;
    bool fiducial=false;
};

FTProjection project_ft(const Branches& b, const FTPlaneEstimate& plane) {
    FTProjection q;
    if (!plane.valid) return q;
    if (!finite_good(b.probe_corr_theta) || !finite_good(b.probe_corr_phi))
        return q;

    const double vz=b.have_e_vz ? b.e_vz : 0.0;
    const double dz=plane.z-vz;
    if (!(dz>0)) return q;

    const double th=deg2rad(b.probe_corr_theta);
    const double ph=deg2rad(b.probe_corr_phi);
    const double rt=dz*std::tan(th);

    q.x=rt*std::cos(ph);
    q.y=rt*std::sin(ph);
    q.r=std::hypot(q.x,q.y);
    q.valid=std::isfinite(q.r);

    if (!q.valid) return q;

    q.radial=(q.r>FT_R_MIN && q.r<FT_R_MAX);
    q.hole_clear=true;
    for (const auto& h : FT_HOLES) {
        if (std::hypot(q.x-h.x,q.y-h.y)<h.r) {
            q.hole_clear=false;
            break;
        }
    } // endfor

    q.fiducial=q.radial && q.hole_clear;
    return q;
}

// -----------------------------------------------------------------------------
// Residual fitting.
// -----------------------------------------------------------------------------
struct FitResult {
    bool valid=false;
    int root_status=-999;
    long long candidates=0;

    // Gaussian signal parameters.
    double amplitude=0;
    double mean=0, mean_err=0;
    double sigma=0, sigma_err=0;

    // Linear-background parameters from gaus(0)+pol1(3).
    double bg0=0;
    double bg1=0;

    // Actual range used for the final fit.
    double fit_lo=0;
    double fit_hi=0;

    double chi2=0;
    int ndf=0;
    std::string reason="not fit";
};

FitResult fit_residual(TH1D* h, int min_entries) {
    FitResult r;
    if (!h) {
        r.reason="null histogram";
        return r;
    }

    r.candidates=static_cast<long long>(h->GetEntries());
    if (h->GetEntries()<min_entries) {
        r.reason="too few candidates";
        return r;
    }

    const double peak=h->GetBinCenter(h->GetMaximumBin());

    const double seed_lo=std::max(FIT_MIN,peak-0.20);
    const double seed_hi=std::min(FIT_MAX,peak+0.20);

    TF1 seed("seed_tmp","gaus",seed_lo,seed_hi);
    seed.SetParameters(h->GetMaximum(),peak,0.10);
    seed.SetParLimits(2,SIGMA_MIN,SIGMA_MAX);
    const int seed_status=h->Fit(&seed,"QNR");

    double mu=peak;
    double sg=0.10;
    if (seed_status==0) {
        mu=seed.GetParameter(1);
        sg=std::fabs(seed.GetParameter(2));
    }
    if (!(sg>SIGMA_MIN && sg<SIGMA_MAX)) sg=0.10;

    double lo=std::max(FIT_MIN,mu-4.0*sg);
    double hi=std::min(FIT_MAX,mu+4.0*sg);
    if (hi-lo<0.25) {
        lo=FIT_MIN;
        hi=FIT_MAX;
    }

    TF1 f("fit_tmp","gaus(0)+pol1(3)",lo,hi);
    f.SetParameters(h->GetMaximum(),mu,sg,
                    std::max(0.0,h->GetBinContent(1)),0.0);
    f.SetParLimits(1,FIT_MEAN_MIN,FIT_MEAN_MAX);
    f.SetParLimits(2,SIGMA_MIN,SIGMA_MAX);

    r.root_status=h->Fit(&f,"QNR");
    r.amplitude=f.GetParameter(0);
    r.mean=f.GetParameter(1);
    r.mean_err=f.GetParError(1);
    r.sigma=std::fabs(f.GetParameter(2));
    r.sigma_err=f.GetParError(2);
    r.bg0=f.GetParameter(3);
    r.bg1=f.GetParameter(4);
    r.fit_lo=lo;
    r.fit_hi=hi;
    r.chi2=f.GetChisquare();
    r.ndf=f.GetNDF();

    if (r.root_status!=0) {
        r.reason="ROOT fit status != 0";
        return r;
    }
    if (r.mean<=FIT_MEAN_MIN+MEAN_BOUNDARY_TOL ||
        r.mean>=FIT_MEAN_MAX-MEAN_BOUNDARY_TOL) {
        r.reason="mean at/near fit boundary";
        return r;
    }
    if (!(r.sigma>SIGMA_MIN*1.01 && r.sigma<SIGMA_MAX*0.99)) {
        r.reason="sigma at/near fit boundary";
        return r;
    }
    if (!(r.sigma_err>0) || r.sigma_err/r.sigma>MAX_REL_SIGMA_ERR) {
        r.reason="sigma uncertainty too large";
        return r;
    }
    if (!(r.mean_err>=0) || r.mean_err>MAX_MEAN_ERR) {
        r.reason="mean uncertainty too large";
        return r;
    }
    if (r.ndf<=0) {
        r.reason="non-positive NDF";
        return r;
    }

    r.valid=true;
    r.reason="valid";
    return r;
}

// -----------------------------------------------------------------------------
// Analysis objects.
// -----------------------------------------------------------------------------
struct RegionResult {
    std::string name;
    int detector=-1;
    double pmin=0,pmax=0;
    double thmin=0,thmax=0;

    std::unique_ptr<TH1D> residual;
    FitResult fit;

    long long denominator_rows=0;
    long long numerator_rows[3]={0,0,0};
    long long reconstructed_candidates=0;

    // Candidate residuals cached during the single tree pass. After the fit is
    // known, exact 1/2/3-sigma matched counts are obtained without rereading
    // the TChain. This roughly halves input I/O for large samples.
    std::vector<double> cached_dp;

    // MC truth closure.
    long long truth_pi0_probe=0;
    long long selected_reco_with_truth=0;
    long long reco_truth_da_lt1=0;
    long long reco_truth_da_lt2=0;
    long long reco_truth_da_lt3=0;
};

struct CutFlow {
    long long all=0;
    long long proton=0;
    long long tag_beta=0;
    long long tag_fid=0;
    long long finite_probe=0;
    long long nominal_mass=0;
    long long fd=0;
    long long ft=0;
    long long ft_projectable=0;
    long long ft_projected_fid=0;
};

struct WindowResult {
    MassWindow w;
    long long fd_denom=0, fd_reco=0;
    long long ft_denom=0, ft_reco=0;
    std::unique_ptr<TH1D> fd_h;
    std::unique_ptr<TH1D> ft_h;
    FitResult fd_fit, ft_fit;
};

struct SampleResult {
    std::string name,input;
    bool is_mc=false;
    long long entries=0;

    // Raw MC::Event.weight diagnostics only. These weights are NOT used in
    // stage-1 residuals or raw recovery fractions until their generator-specific
    // meaning and the full Valerii normalization prescription are established.
    long long weight_n=0, weight_nonfinite=0, weight_zero=0, weight_negative=0;
    double weight_sum=0, weight_sum2=0;
    double weight_min=0, weight_max=0;

    CutFlow cutflow;
    FTPlaneEstimate ft_plane;

    RegionResult fd;
    RegionResult ft_all;
    RegionResult ft_low;
    RegionResult ft_high;

    std::vector<WindowResult> windows;

    std::unique_ptr<TH1D> denom_p;
    std::unique_ptr<TH1D> denom_theta;
    std::unique_ptr<TH1D> denom_phi;
    std::unique_ptr<TH2D> ft_xy_projected;
};

void init_region(RegionResult& r,
                 const std::string& sample,
                 const std::string& name,
                 int detector,
                 double pmin,double pmax,
                 double thmin,double thmax) {
    r.name=name;
    r.detector=detector;
    r.pmin=pmin; r.pmax=pmax;
    r.thmin=thmin; r.thmax=thmax;

    r.residual.reset(new TH1D(
        Form("h_%s_%s_dp",sample.c_str(),name.c_str()),
        Form("%s %s;#Delta p_{#gamma2}=p_{rec}-p_{miss} (GeV);Candidates",
             sample.c_str(),name.c_str()),
        DP_NBIN,DP_HMIN,DP_HMAX));
    r.residual->Sumw2();
    r.residual->SetDirectory(nullptr);
    r.cached_dp.reserve(10000);
}

bool in_region(const Branches& b, const RegionResult& r) {
    return b.probe_corr_p>=r.pmin && b.probe_corr_p<=r.pmax &&
           b.probe_corr_theta>=r.thmin && b.probe_corr_theta<=r.thmax;
}

void init_sample(SampleResult& s) {
    init_region(s.fd,s.name,"FD",1,FD_P_MIN,FD_P_MAX,FD_THETA_MIN,FD_THETA_MAX);
    init_region(s.ft_all,s.name,"FT_all",0,FT_P_MIN,FT_P_MAX,FT_THETA_MIN,FT_THETA_MAX);
    init_region(s.ft_low,s.name,"FT_Elt2",0,FT_P_MIN,FT_P_SPLIT,FT_THETA_MIN,FT_THETA_MAX);
    init_region(s.ft_high,s.name,"FT_Ege2",0,FT_P_SPLIT,FT_P_MAX,FT_THETA_MIN,FT_THETA_MAX);

    s.denom_p.reset(new TH1D(Form("h_%s_denom_p",s.name.c_str()),
                             ";missing-#gamma_{2} p (GeV);Candidates",80,0,8));
    s.denom_theta.reset(new TH1D(Form("h_%s_denom_theta",s.name.c_str()),
                                 ";missing-#gamma_{2} #theta (deg);Candidates",80,0,40));
    s.denom_phi.reset(new TH1D(Form("h_%s_denom_phi",s.name.c_str()),
                               ";missing-#gamma_{2} wrapped #phi (deg);Candidates",72,-30,330));
    s.ft_xy_projected.reset(new TH2D(Form("h_%s_ft_xy",s.name.c_str()),
                                     ";projected FT x (cm);projected FT y (cm)",
                                     100,-20,20,100,-20,20));

    s.denom_p->Sumw2(); s.denom_theta->Sumw2(); s.denom_phi->Sumw2();
    s.denom_p->SetDirectory(nullptr);
    s.denom_theta->SetDirectory(nullptr);
    s.denom_phi->SetDirectory(nullptr);
    s.ft_xy_projected->SetDirectory(nullptr);

    for (const auto& mw : MASS_WINDOWS) {
        WindowResult wr;
        wr.w=mw;
        wr.fd_h.reset(new TH1D(Form("h_%s_%s_fd",s.name.c_str(),mw.label),
                               ";#Delta p (GeV);Candidates",DP_NBIN,DP_HMIN,DP_HMAX));
        wr.ft_h.reset(new TH1D(Form("h_%s_%s_ft",s.name.c_str(),mw.label),
                               ";#Delta p (GeV);Candidates",DP_NBIN,DP_HMIN,DP_HMAX));
        wr.fd_h->Sumw2(); wr.ft_h->Sumw2();
        wr.fd_h->SetDirectory(nullptr); wr.ft_h->SetDirectory(nullptr);
        s.windows.push_back(std::move(wr));
    } // endfor
}

double recovery_fraction(const RegionResult& r, int ns) {
    if (r.denominator_rows<=0) return -1;
    return static_cast<double>(r.numerator_rows[ns-1]) /
           static_cast<double>(r.denominator_rows);
}

double recovery_fraction_error(const RegionResult& r, int ns) {
    const double f=recovery_fraction(r,ns);
    if (!(f>=0 && f<=1) || r.denominator_rows<=0) return 0;
    return std::sqrt(std::max(0.0,f*(1.0-f)/
                             static_cast<double>(r.denominator_rows)));
}

double raw_weight_neff(const SampleResult& s) {
    if (s.weight_sum2<=0) return 0;
    return s.weight_sum*s.weight_sum/s.weight_sum2;
}

// -----------------------------------------------------------------------------
// Fill helpers.
// -----------------------------------------------------------------------------
void fill_region_residual(RegionResult& r, const Branches& b,
                          const FTProjection* ft_projection=nullptr) {
    if (!in_region(b,r)) return;

    // For FT regions, angular acceptance alone is insufficient because the FT
    // is annular and contains excluded holes.  Only inferred probes that
    // project into the usable FT fiducial area belong in the denominator.
    if (r.detector==0) {
        if (!ft_projection || !ft_projection->valid || !ft_projection->fiducial)
            return;
    }

    r.denominator_rows++;
    const int k=best_probe_candidate(b,r.detector);
    if (k<0) return;

    r.reconstructed_candidates++;
    const double dp=b.neutral_p[k]-b.probe_corr_p;
    if (std::isfinite(dp)) {
        r.residual->Fill(dp);
        r.cached_dp.push_back(dp);
    }

    if (b.have_truth &&
        b.truth_probe_pid==22 &&
        b.truth_probe_parent==111) {

        r.truth_pi0_probe++;
        const double da_truth=opening_angle_deg(
            b.neutral_theta[k],b.neutral_phi[k],
            b.truth_probe_theta,b.truth_probe_phi);

        if (finite_good(da_truth)) {
            r.selected_reco_with_truth++;
            if (da_truth<1.0) r.reco_truth_da_lt1++;
            if (da_truth<2.0) r.reco_truth_da_lt2++;
            if (da_truth<3.0) r.reco_truth_da_lt3++;
        }
    }
}

void finalize_region_counts(RegionResult& r) {
    if (!r.fit.valid || r.denominator_rows<=0) return;

    for (double dp : r.cached_dp) {
        const double d=std::fabs(dp-r.fit.mean);
        for (int ns=1;ns<=3;ns++) {
            if (d<ns*r.fit.sigma) r.numerator_rows[ns-1]++;
        } // endfor
    } // endfor
}

// -----------------------------------------------------------------------------
// Main per-sample analysis.
// -----------------------------------------------------------------------------
bool analyze_sample(SampleResult& s) {
    const std::string pattern=make_pattern(s.input);
    TChain c("PhotonEfficiency");
    const int nf=c.Add(pattern.c_str());
    s.entries=c.GetEntries();

    std::cout << "\n[" << s.name << "] files=" << nf
              << " entries=" << s.entries << "\n";

    if (nf<=0 || s.entries<=0) return false;

    Branches b;
    b.reset_arrays();
    if (!attach(c,b)) return false;
    c.SetCacheSize(64LL*1024LL*1024LL);
    c.AddBranchToCache("*",kTRUE);

    init_sample(s);
    s.ft_plane=estimate_ft_plane(c,b);

    std::cout << "[" << s.name << "] inferred FT response-plane z: ";
    if (s.ft_plane.valid)
        std::cout << s.ft_plane.z << " cm from " << s.ft_plane.n << " responses\n";
    else
        std::cout << "unavailable\n";

    // Pass 1: cut flow, residuals, kinematics, mass-window scan.
    for (Long64_t i=0;i<c.GetEntries();i++) {
        c.GetEntry(i);
        s.cutflow.all++;

        // Preserve MC::Event.weight only as a diagnostic.  It is deliberately
        // NOT applied to stage-1 histograms or recovery fractions.
        if (s.is_mc) {
            if (!std::isfinite(b.mc_weight)) {
                s.weight_nonfinite++;
            } else {
                if (s.weight_n==0) {
                    s.weight_min=b.mc_weight;
                    s.weight_max=b.mc_weight;
                } else {
                    s.weight_min=std::min(s.weight_min,b.mc_weight);
                    s.weight_max=std::max(s.weight_max,b.mc_weight);
                }
                s.weight_n++;
                if (b.mc_weight==0) s.weight_zero++;
                if (b.mc_weight<0) s.weight_negative++;
                s.weight_sum+=b.mc_weight;
                s.weight_sum2+=b.mc_weight*b.mc_weight;
            }
        }


        if (!b.p_pass_standard) continue;
        s.cutflow.proton++;

        if (!b.tag_pass_beta) continue;
        s.cutflow.tag_beta++;

        if (!b.tag_pass_fiducial) continue;
        s.cutflow.tag_fid++;

        if (!(b.tag_detector==0 || b.tag_detector==1)) continue;
        if (!finite_good(b.probe_corr_p) ||
            !finite_good(b.probe_corr_theta) ||
            !finite_good(b.probe_corr_phi)) continue;
        s.cutflow.finite_probe++;

        // Compute the projected FT intersection once per event.  It is used
        // both for the parent-window scan and for the nominal FT denominator.
        FTProjection ft_projection;
        const bool ft_angular=in_ft(b);
        if (ft_angular) ft_projection=project_ft(b,s.ft_plane);

        // Parent-window stability scan before imposing the nominal window.
        for (auto& wr : s.windows) {
            if (!finite_good(b.Mx_ep) ||
                b.Mx_ep<wr.w.lo || b.Mx_ep>wr.w.hi) continue;

            if (in_fd(b)) {
                wr.fd_denom++;
                const int k=best_probe_candidate(b,1);
                if (k>=0) {
                    wr.fd_reco++;
                    wr.fd_h->Fill(b.neutral_p[k]-b.probe_corr_p);
                }
            }

            // FT parent-window results now use the same projected-face
            // fiducial requirement as the nominal FT denominator.
            if (ft_angular && ft_projection.valid && ft_projection.fiducial) {
                wr.ft_denom++;
                const int k=best_probe_candidate(b,0);
                if (k>=0) {
                    wr.ft_reco++;
                    wr.ft_h->Fill(b.neutral_p[k]-b.probe_corr_p);
                }
            }
        } // endfor

        if (!pass_nominal_mass(b)) continue;
        s.cutflow.nominal_mass++;

        s.denom_p->Fill(b.probe_corr_p);
        s.denom_theta->Fill(b.probe_corr_theta);
        s.denom_phi->Fill(wrap_phi(b.probe_corr_phi));

        if (in_fd(b)) s.cutflow.fd++;
        if (ft_angular) {
            s.cutflow.ft++;
            if (ft_projection.valid) {
                s.cutflow.ft_projectable++;
                s.ft_xy_projected->Fill(ft_projection.x,ft_projection.y);
                if (ft_projection.fiducial) s.cutflow.ft_projected_fid++;
            }
        }

        fill_region_residual(s.fd,b);
        fill_region_residual(s.ft_all,b,&ft_projection);
        fill_region_residual(s.ft_low,b,&ft_projection);
        fill_region_residual(s.ft_high,b,&ft_projection);
    } // endfor

    s.fd.fit=fit_residual(s.fd.residual.get(),MIN_FIT_ENTRIES_FD);
    s.ft_all.fit=fit_residual(s.ft_all.residual.get(),MIN_FIT_ENTRIES_FT);
    s.ft_low.fit=fit_residual(s.ft_low.residual.get(),MIN_FIT_ENTRIES_FT);
    s.ft_high.fit=fit_residual(s.ft_high.residual.get(),MIN_FIT_ENTRIES_FT);

    for (auto& wr : s.windows) {
        wr.fd_fit=fit_residual(wr.fd_h.get(),MIN_FIT_ENTRIES_FD);
        wr.ft_fit=fit_residual(wr.ft_h.get(),MIN_FIT_ENTRIES_FT);
    } // endfor

    // Exact 1/2/3-sigma raw-recovery counts from residuals cached during
    // the single tree pass.  This avoids rereading every ROOT entry.
    finalize_region_counts(s.fd);
    finalize_region_counts(s.ft_all);
    finalize_region_counts(s.ft_low);
    finalize_region_counts(s.ft_high);

    return true;
}

// -----------------------------------------------------------------------------
// Output.
// -----------------------------------------------------------------------------
void draw_residual_pad(const SampleResult& s, const RegionResult& r) {
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.07);

    r.residual->SetLineWidth(2);
    r.residual->GetXaxis()->SetTitleSize(0.045);
    r.residual->GetYaxis()->SetTitleSize(0.045);
    r.residual->GetXaxis()->SetLabelSize(0.038);
    r.residual->GetYaxis()->SetLabelSize(0.038);
    r.residual->Draw("E");

    const double ymax=std::max(1.0,r.residual->GetMaximum());

    if (r.fit.valid) {
        TF1 fit_draw(
            Form("fit_draw_%s_%s",s.name.c_str(),r.name.c_str()),
            "gaus(0)+pol1(3)",
            r.fit.fit_lo,r.fit.fit_hi);
        fit_draw.SetParameters(
            r.fit.amplitude,
            r.fit.mean,
            r.fit.sigma,
            r.fit.bg0,
            r.fit.bg1);
        fit_draw.SetLineWidth(3);
        fit_draw.DrawCopy("SAME");

        for (int ns=1;ns<=3;ns++) {
            TLine l1(r.fit.mean-ns*r.fit.sigma,0,
                     r.fit.mean-ns*r.fit.sigma,0.92*ymax);
            TLine l2(r.fit.mean+ns*r.fit.sigma,0,
                     r.fit.mean+ns*r.fit.sigma,0.92*ymax);
            l1.SetLineStyle(ns);
            l2.SetLineStyle(ns);
            l1.DrawClone();
            l2.DrawClone();
        } // endfor
    }

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.034);
    // Deliberately lower than the previous version so the first line is
    // comfortably clear of the upper frame, even in the combined canvas.
    tx.DrawLatex(0.15,0.76,Form("%s %s",s.name.c_str(),r.name.c_str()));
    tx.DrawLatex(0.15,0.71,Form("N_{den}=%lld, N_{reco}=%lld",
                                r.denominator_rows,r.reconstructed_candidates));

    if (r.fit.valid) {
        tx.DrawLatex(0.15,0.66,
            Form("#mu=%.4f#pm%.4f GeV",r.fit.mean,r.fit.mean_err));
        tx.DrawLatex(0.15,0.61,
            Form("#sigma=%.4f#pm%.4f GeV",r.fit.sigma,r.fit.sigma_err));
        tx.DrawLatex(0.15,0.56,
            Form("#chi^{2}/ndf=%.1f/%d",r.fit.chi2,r.fit.ndf));
    } else {
        tx.DrawLatex(0.15,0.66,Form("FIT INVALID: %s",r.fit.reason.c_str()));
    }
}

void save_combined_residuals(
        const std::vector<std::unique_ptr<SampleResult>>& samples,
        const std::string& out) {
    if (samples.empty()) return;

    const int ncol=samples.size();
    const int nrow=4;
    TCanvas c("c_all_residuals","",520*ncol,420*nrow);
    c.Divide(ncol,nrow,0.001,0.001);

    for (int is=0;is<ncol;is++) {
        const SampleResult& s=*samples[is];
        const RegionResult* rr[] = {&s.fd,&s.ft_all,&s.ft_low,&s.ft_high};
        for (int ir=0;ir<nrow;ir++) {
            c.cd(ir*ncol+is+1);
            draw_residual_pad(s,*rr[ir]);
        } // endfor
    } // endfor

    c.SaveAs((out+"/delta_p_residuals.png").c_str());
}

void save_combined_kinematics(
        const std::vector<std::unique_ptr<SampleResult>>& samples,
        const std::string& out) {
    if (samples.empty()) return;

    const int nrow=samples.size();
    TCanvas c("c_all_kinematics","",1500,390*nrow);
    c.Divide(3,nrow,0.001,0.001);

    for (int is=0;is<nrow;is++) {
        const SampleResult& s=*samples[is];
        TH1D* hh[] = {s.denom_p.get(),s.denom_theta.get(),s.denom_phi.get()};
        for (int j=0;j<3;j++) {
            c.cd(is*3+j+1);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            hh[j]->SetTitle(Form("%s;%s;Candidates",
                                 s.name.c_str(),hh[j]->GetXaxis()->GetTitle()));
            hh[j]->Draw("E");
        } // endfor
    } // endfor

    c.SaveAs((out+"/denominator_kinematics.png").c_str());
}

void save_combined_ft_projection(
        const std::vector<std::unique_ptr<SampleResult>>& samples,
        const std::string& out) {
    int nvalid=0;
    for (const auto& sp : samples) if (sp->ft_plane.valid) nvalid++;
    if (nvalid==0) return;

    const int ncol=std::min(2,nvalid);
    const int nrow=(nvalid+ncol-1)/ncol;
    TCanvas c("c_all_ft_xy","",760*ncol,680*nrow);
    c.Divide(ncol,nrow,0.002,0.002);

    int ipad=0;
    for (const auto& sp : samples) {
        if (!sp->ft_plane.valid) continue;
        c.cd(++ipad);
        gPad->SetRightMargin(0.15);
        sp->ft_xy_projected->SetTitle(
            Form("%s;projected FT x (cm);projected FT y (cm)",sp->name.c_str()));
        sp->ft_xy_projected->Draw("COLZ");
    } // endfor

    c.SaveAs((out+"/FT_projected_xy.png").c_str());
}

void write_window_scan_all(
        const std::vector<std::unique_ptr<SampleResult>>& samples,
        const std::string& out) {
    std::ofstream f(out+"/parent_window_scan.csv");
    f << "sample,window,mmin,mmax,fd_denom,fd_reco,fd_fit_valid,fd_mean,fd_sigma,"
      << "ft_denom,ft_reco,ft_fit_valid,ft_mean,ft_sigma\n";

    for (const auto& sp : samples) {
        for (const auto& wr : sp->windows) {
            f << sp->name << "," << wr.w.label << "," << wr.w.lo << "," << wr.w.hi << ","
              << wr.fd_denom << "," << wr.fd_reco << ","
              << (wr.fd_fit.valid?1:0) << ","
              << wr.fd_fit.mean << "," << wr.fd_fit.sigma << ","
              << wr.ft_denom << "," << wr.ft_reco << ","
              << (wr.ft_fit.valid?1:0) << ","
              << wr.ft_fit.mean << "," << wr.ft_fit.sigma << "\n";
        } // endfor
    } // endfor
}

void write_region_csv_header(std::ofstream& f) {
    f << "sample,region,pmin,pmax,thetamin,thetamax,is_mc,"
      << "fit_valid,fit_reason,fit_candidates,mean,mean_err,sigma,sigma_err,"
      << "chi2,ndf,denominator_rows,reconstructed_rows,"
      << "matched_rows_1sigma,matched_rows_2sigma,matched_rows_3sigma,"
      << "raw_recovery_1sigma,raw_recovery_2sigma,raw_recovery_3sigma,"
      << "raw_recovery_err_1sigma,raw_recovery_err_2sigma,raw_recovery_err_3sigma,"
      << "truth_pi0_probe,selected_reco_with_truth,"
      << "reco_truth_da_lt1,reco_truth_da_lt2,reco_truth_da_lt3,"
      << "mc_weight_n,mc_weight_mean,mc_weight_min,mc_weight_max,mc_weight_neff,"
      << "mc_weight_zero,mc_weight_negative,mc_weight_nonfinite\n";
}

void append_region_csv(std::ofstream& f,
                       const SampleResult& s,
                       const RegionResult& r) {
    const double wmean=s.weight_n>0 ? s.weight_sum/static_cast<double>(s.weight_n) : 0.0;
    f << s.name << "," << r.name << ","
      << r.pmin << "," << r.pmax << ","
      << r.thmin << "," << r.thmax << ","
      << (s.is_mc?1:0) << ","
      << (r.fit.valid?1:0) << ",\"" << r.fit.reason << "\","
      << r.fit.candidates << ","
      << r.fit.mean << "," << r.fit.mean_err << ","
      << r.fit.sigma << "," << r.fit.sigma_err << ","
      << r.fit.chi2 << "," << r.fit.ndf << ","
      << r.denominator_rows << "," << r.reconstructed_candidates << ","
      << r.numerator_rows[0] << "," << r.numerator_rows[1] << "," << r.numerator_rows[2] << ","
      << recovery_fraction(r,1) << "," << recovery_fraction(r,2) << "," << recovery_fraction(r,3) << ","
      << recovery_fraction_error(r,1) << ","
      << recovery_fraction_error(r,2) << ","
      << recovery_fraction_error(r,3) << ","
      << r.truth_pi0_probe << "," << r.selected_reco_with_truth << ","
      << r.reco_truth_da_lt1 << "," << r.reco_truth_da_lt2 << ","
      << r.reco_truth_da_lt3 << ","
      << s.weight_n << "," << wmean << "," << s.weight_min << "," << s.weight_max << ","
      << raw_weight_neff(s) << "," << s.weight_zero << "," << s.weight_negative << ","
      << s.weight_nonfinite << "\n";
}

void write_all_summary(
        const std::vector<std::unique_ptr<SampleResult>>& samples,
        const std::string& out) {
    std::ofstream f(out+"/analysis_summary.txt");

    f << "Photon tag-and-probe stage-1 diagnostic summary\n"
      << "==============================================\n"
      << "Nominal parent window: " << PI0_M_MIN
      << " < Mx(ep) < " << PI0_M_MAX << " GeV\n\n"
      << "IMPORTANT INTERPRETATION\n"
      << "------------------------\n"
      << "The 1/2/3-sigma quantities below are RAW RECOVERY FRACTIONS, not\n"
      << "physical photon efficiencies.  In data, the ep-gamma denominator is\n"
      << "not yet separated into true ep-pi0 and non-pi0 components.  Therefore\n"
      << "no data/MC efficiency ratio or cross-section correction is reported at\n"
      << "this stage.  The physical efficiency will be formed only after the\n"
      << "Valerii-style normalized AAOgen + CLASDIS + DVCSgen composition /\n"
      << "background procedure is reproduced.\n\n"
      << "MC::Event.weight POLICY\n"
      << "-----------------------\n"
      << "MC::Event.weight is copied into the skim but its generator-specific\n"
      << "meaning has not yet been established.  It is NOT used in residual fits,\n"
      << "histograms, or raw recovery fractions.  Its distribution is printed only\n"
      << "as QA so that the later normalization prescription can be audited.\n\n";

    for (const auto& sp : samples) {
        const SampleResult& s=*sp;
        f << "============================================================\n"
          << "SAMPLE: " << s.name << "\n"
          << "Input: " << s.input << "\n"
          << "Tree entries: " << s.entries << "\n"
          << "MC sample: " << (s.is_mc?"yes":"no") << "\n\n";

        if (s.is_mc) {
            const double mean=s.weight_n>0 ? s.weight_sum/static_cast<double>(s.weight_n) : 0.0;
            f << "Raw MC::Event.weight QA\n"
              << "-----------------------\n"
              << "finite entries: " << s.weight_n << "\n"
              << "mean/min/max: " << mean << " / " << s.weight_min << " / " << s.weight_max << "\n"
              << "N_eff if these weights were used: " << raw_weight_neff(s) << "\n"
              << "zero / negative / nonfinite: " << s.weight_zero << " / "
              << s.weight_negative << " / " << s.weight_nonfinite << "\n"
              << "NOTE: these weights are not applied in stage 1.\n\n";
        }

        f << "Cut flow\n--------\n";
        auto row=[&](const char* label,long long n) {
            const double pct=s.cutflow.all>0 ? 100.0*n/s.cutflow.all : 0;
            f << std::left << std::setw(30) << label
              << std::right << std::setw(12) << n
              << "  " << std::fixed << std::setprecision(3) << pct
              << "% of skim rows\n";
        };
        row("all skim rows",s.cutflow.all);
        row("standard proton",s.cutflow.proton);
        row("tag beta",s.cutflow.tag_beta);
        row("tag fiducial",s.cutflow.tag_fid);
        row("finite probe kinematics",s.cutflow.finite_probe);
        row("nominal pi0 parent window",s.cutflow.nominal_mass);
        row("FD analysis region",s.cutflow.fd);
        row("FT angular region",s.cutflow.ft);
        row("FT projectable",s.cutflow.ft_projectable);
        row("FT projected fiducial",s.cutflow.ft_projected_fid);

        if (s.ft_plane.valid)
            f << "FT response-plane z: " << s.ft_plane.z
              << " cm from " << s.ft_plane.n << " reconstructed FT responses\n";
        f << "FT denominator requires projected FT fiducial acceptance.\n"
          << "The angular FT row above is retained as a geometry diagnostic; "
          << "the region denominators below use the projected-fiducial subset.\n\n";

        const RegionResult* rr[] = {&s.fd,&s.ft_all,&s.ft_low,&s.ft_high};
        for (const RegionResult* r : rr) {
            f << "--- " << r->name << " ---\n"
              << "denominator rows: " << r->denominator_rows << "\n"
              << "reconstructed candidates: " << r->reconstructed_candidates << "\n"
              << "fit valid: " << (r->fit.valid?1:0)
              << " (" << r->fit.reason << ")\n"
              << "mean: " << r->fit.mean << " +/- " << r->fit.mean_err << " GeV\n"
              << "sigma: " << r->fit.sigma << " +/- " << r->fit.sigma_err << " GeV\n";

            for (int ns=1;ns<=3;ns++) {
                f << ns << " sigma raw recovery fraction: "
                  << recovery_fraction(*r,ns) << " +/- "
                  << recovery_fraction_error(*r,ns) << "\n";
            } // endfor

            if (r->truth_pi0_probe>0) {
                f << "MC truth pi0-probe rows among reconstructed candidates: " << r->truth_pi0_probe << "\n"
                  << "selected reco candidates with truth comparison: "
                  << r->selected_reco_with_truth << "\n"
                  << "reco-to-truth dAlpha <1/<2/<3 deg: "
                  << r->reco_truth_da_lt1 << " / "
                  << r->reco_truth_da_lt2 << " / "
                  << r->reco_truth_da_lt3 << "\n";
            }
            f << "\n";
        } // endfor
    } // endfor

    f << "============================================================\n"
      << "No data/MC efficiency correction is produced in stage 1.\n";
}

// -----------------------------------------------------------------------------
// Worker-file serialization for process-level parallelism.
// Each child ROOT process image analyzes exactly one independent sample and
// writes one temporary ROOT file.  The parent waits, reloads those files, then
// creates the same compact combined output as before.  This uses POSIX fork(),
// not std::thread/std::async, avoiding the Cling __once_call TLS linker issue.
// -----------------------------------------------------------------------------
template <typename T>
void put_param(TDirectory* d, const std::string& name, T value) {
    d->cd();
    TParameter<T> p(name.c_str(),value);
    p.Write();
}

template <typename T>
bool get_param(TDirectory* d, const std::string& name, T& value) {
    auto* p=dynamic_cast<TParameter<T>*>(d->Get(name.c_str()));
    if (!p) return false;
    value=p->GetVal();
    return true;
}

void write_fit_state(TDirectory* d, const std::string& pre, const FitResult& r) {
    put_param<int>(d,pre+"valid",r.valid?1:0);
    put_param<int>(d,pre+"root_status",r.root_status);
    put_param<Long64_t>(d,pre+"candidates",r.candidates);
    put_param<double>(d,pre+"amplitude",r.amplitude);
    put_param<double>(d,pre+"mean",r.mean);
    put_param<double>(d,pre+"mean_err",r.mean_err);
    put_param<double>(d,pre+"sigma",r.sigma);
    put_param<double>(d,pre+"sigma_err",r.sigma_err);
    put_param<double>(d,pre+"bg0",r.bg0);
    put_param<double>(d,pre+"bg1",r.bg1);
    put_param<double>(d,pre+"fit_lo",r.fit_lo);
    put_param<double>(d,pre+"fit_hi",r.fit_hi);
    put_param<double>(d,pre+"chi2",r.chi2);
    put_param<int>(d,pre+"ndf",r.ndf);
    TNamed reason((pre+"reason").c_str(),r.reason.c_str());
    reason.Write();
}

void read_fit_state(TDirectory* d, const std::string& pre, FitResult& r) {
    int iv=0;
    get_param<int>(d,pre+"valid",iv); r.valid=(iv!=0);
    get_param<int>(d,pre+"root_status",r.root_status);
    Long64_t nc=0; get_param<Long64_t>(d,pre+"candidates",nc); r.candidates=nc;
    get_param<double>(d,pre+"amplitude",r.amplitude);
    get_param<double>(d,pre+"mean",r.mean);
    get_param<double>(d,pre+"mean_err",r.mean_err);
    get_param<double>(d,pre+"sigma",r.sigma);
    get_param<double>(d,pre+"sigma_err",r.sigma_err);
    get_param<double>(d,pre+"bg0",r.bg0);
    get_param<double>(d,pre+"bg1",r.bg1);
    get_param<double>(d,pre+"fit_lo",r.fit_lo);
    get_param<double>(d,pre+"fit_hi",r.fit_hi);
    get_param<double>(d,pre+"chi2",r.chi2);
    get_param<int>(d,pre+"ndf",r.ndf);
    auto* n=dynamic_cast<TNamed*>(d->Get((pre+"reason").c_str()));
    r.reason=n?n->GetTitle():"unavailable";
}

void write_region_state(TDirectory* parent, const RegionResult& r) {
    TDirectory* d=parent->mkdir(r.name.c_str());
    d->cd();
    r.residual->Write("residual");
    put_param<int>(d,"detector",r.detector);
    put_param<double>(d,"pmin",r.pmin); put_param<double>(d,"pmax",r.pmax);
    put_param<double>(d,"thmin",r.thmin); put_param<double>(d,"thmax",r.thmax);
    put_param<Long64_t>(d,"denominator_rows",r.denominator_rows);
    put_param<Long64_t>(d,"reconstructed_candidates",r.reconstructed_candidates);
    for (int i=0;i<3;i++) put_param<Long64_t>(d,Form("numerator_rows_%d",i),r.numerator_rows[i]);
    put_param<Long64_t>(d,"truth_pi0_probe",r.truth_pi0_probe);
    put_param<Long64_t>(d,"selected_reco_with_truth",r.selected_reco_with_truth);
    put_param<Long64_t>(d,"reco_truth_da_lt1",r.reco_truth_da_lt1);
    put_param<Long64_t>(d,"reco_truth_da_lt2",r.reco_truth_da_lt2);
    put_param<Long64_t>(d,"reco_truth_da_lt3",r.reco_truth_da_lt3);
    write_fit_state(d,"fit_",r.fit);
}

void read_region_state(TDirectory* parent, RegionResult& r) {
    TDirectory* d=dynamic_cast<TDirectory*>(parent->Get(r.name.c_str()));
    if (!d) return;
    if (auto* h=dynamic_cast<TH1D*>(d->Get("residual"))) {
        r.residual.reset(dynamic_cast<TH1D*>(h->Clone()));
        r.residual->SetDirectory(nullptr);
    }
    get_param<int>(d,"detector",r.detector);
    get_param<double>(d,"pmin",r.pmin); get_param<double>(d,"pmax",r.pmax);
    get_param<double>(d,"thmin",r.thmin); get_param<double>(d,"thmax",r.thmax);
    Long64_t ll=0;
    get_param<Long64_t>(d,"denominator_rows",ll); r.denominator_rows=ll;
    get_param<Long64_t>(d,"reconstructed_candidates",ll); r.reconstructed_candidates=ll;
    for (int i=0;i<3;i++) { get_param<Long64_t>(d,Form("numerator_rows_%d",i),ll); r.numerator_rows[i]=ll; }
    get_param<Long64_t>(d,"truth_pi0_probe",ll); r.truth_pi0_probe=ll;
    get_param<Long64_t>(d,"selected_reco_with_truth",ll); r.selected_reco_with_truth=ll;
    get_param<Long64_t>(d,"reco_truth_da_lt1",ll); r.reco_truth_da_lt1=ll;
    get_param<Long64_t>(d,"reco_truth_da_lt2",ll); r.reco_truth_da_lt2=ll;
    get_param<Long64_t>(d,"reco_truth_da_lt3",ll); r.reco_truth_da_lt3=ll;
    read_fit_state(d,"fit_",r.fit);
}


void write_root(const std::vector<std::unique_ptr<SampleResult>>& samples,
                const std::string& out) {
    TFile f((out+"/analysis_histograms.root").c_str(),"RECREATE");
    if (f.IsZombie()) {
        std::cerr << "WARNING: could not create " << out << "/analysis_histograms.root\n";
        return;
    }

    for (const auto& sp : samples) {
        if (!sp) continue;
        f.mkdir(sp->name.c_str());
        f.cd(sp->name.c_str());

        if (sp->fd.residual) sp->fd.residual->Write();
        if (sp->ft_all.residual) sp->ft_all.residual->Write();
        if (sp->ft_low.residual) sp->ft_low.residual->Write();
        if (sp->ft_high.residual) sp->ft_high.residual->Write();

        if (sp->denom_p) sp->denom_p->Write();
        if (sp->denom_theta) sp->denom_theta->Write();
        if (sp->denom_phi) sp->denom_phi->Write();
        if (sp->ft_xy_projected) sp->ft_xy_projected->Write();

        for (const auto& wr : sp->windows) {
            if (wr.fd_h) wr.fd_h->Write();
            if (wr.ft_h) wr.ft_h->Write();
        } // endfor
        f.cd();
    } // endfor

    f.Close();
}

bool save_worker_result(const SampleResult& s, const std::string& path) {
    TFile f(path.c_str(),"RECREATE");
    if (f.IsZombie()) {
        std::cerr << "ERROR: cannot create cache file " << path
                  << ". Check free space/quota on the cache filesystem.\n";
        return false;
    }
    TNamed nname("sample_name",s.name.c_str()); nname.Write();
    TNamed ninput("input",s.input.c_str()); ninput.Write();
    put_param<int>(&f,"is_mc",s.is_mc?1:0);
    put_param<Long64_t>(&f,"entries",s.entries);
    put_param<Long64_t>(&f,"weight_n",s.weight_n);
    put_param<Long64_t>(&f,"weight_nonfinite",s.weight_nonfinite);
    put_param<Long64_t>(&f,"weight_zero",s.weight_zero);
    put_param<Long64_t>(&f,"weight_negative",s.weight_negative);
    put_param<double>(&f,"weight_sum",s.weight_sum);
    put_param<double>(&f,"weight_sum2",s.weight_sum2);
    put_param<double>(&f,"weight_min",s.weight_min);
    put_param<double>(&f,"weight_max",s.weight_max);

    put_param<Long64_t>(&f,"cf_all",s.cutflow.all);
    put_param<Long64_t>(&f,"cf_proton",s.cutflow.proton);
    put_param<Long64_t>(&f,"cf_tag_beta",s.cutflow.tag_beta);
    put_param<Long64_t>(&f,"cf_tag_fid",s.cutflow.tag_fid);
    put_param<Long64_t>(&f,"cf_finite_probe",s.cutflow.finite_probe);
    put_param<Long64_t>(&f,"cf_nominal_mass",s.cutflow.nominal_mass);
    put_param<Long64_t>(&f,"cf_fd",s.cutflow.fd);
    put_param<Long64_t>(&f,"cf_ft",s.cutflow.ft);
    put_param<Long64_t>(&f,"cf_ft_projectable",s.cutflow.ft_projectable);
    put_param<Long64_t>(&f,"cf_ft_projected_fid",s.cutflow.ft_projected_fid);
    put_param<int>(&f,"ft_plane_valid",s.ft_plane.valid?1:0);
    put_param<double>(&f,"ft_plane_z",s.ft_plane.z);
    put_param<Long64_t>(&f,"ft_plane_n",s.ft_plane.n);

    f.mkdir("regions");
    TDirectory* rd=dynamic_cast<TDirectory*>(f.Get("regions"));
    write_region_state(rd,s.fd); write_region_state(rd,s.ft_all);
    write_region_state(rd,s.ft_low); write_region_state(rd,s.ft_high);

    f.mkdir("qa");
    TDirectory* qd=dynamic_cast<TDirectory*>(f.Get("qa")); qd->cd();
    s.denom_p->Write("denom_p"); s.denom_theta->Write("denom_theta");
    s.denom_phi->Write("denom_phi"); s.ft_xy_projected->Write("ft_xy_projected");

    f.mkdir("windows");
    TDirectory* wd=dynamic_cast<TDirectory*>(f.Get("windows"));
    for (size_t i=0;i<s.windows.size();i++) {
        const auto& wr=s.windows[i];
        TDirectory* d=wd->mkdir(Form("w%zu",i)); d->cd();
        TNamed label("label",wr.w.label); label.Write();
        put_param<double>(d,"lo",wr.w.lo); put_param<double>(d,"hi",wr.w.hi);
        put_param<Long64_t>(d,"fd_denom",wr.fd_denom); put_param<Long64_t>(d,"fd_reco",wr.fd_reco);
        put_param<Long64_t>(d,"ft_denom",wr.ft_denom); put_param<Long64_t>(d,"ft_reco",wr.ft_reco);
        wr.fd_h->Write("fd_h"); wr.ft_h->Write("ft_h");
        write_fit_state(d,"fd_fit_",wr.fd_fit); write_fit_state(d,"ft_fit_",wr.ft_fit);
    } // endfor
    f.Close();
    return true;
}

std::unique_ptr<SampleResult> load_worker_result(const std::string& path) {
    TFile f(path.c_str(),"READ");
    if (f.IsZombie()) return nullptr;
    auto* nn=dynamic_cast<TNamed*>(f.Get("sample_name"));
    auto* ni=dynamic_cast<TNamed*>(f.Get("input"));
    if (!nn || !ni) return nullptr;

    std::unique_ptr<SampleResult> s(new SampleResult);
    s->name=nn->GetTitle(); s->input=ni->GetTitle();
    int imc=0; get_param<int>(&f,"is_mc",imc); s->is_mc=(imc!=0);
    Long64_t ll=0;
    get_param<Long64_t>(&f,"entries",ll); s->entries=ll;
    get_param<Long64_t>(&f,"weight_n",ll); s->weight_n=ll;
    get_param<Long64_t>(&f,"weight_nonfinite",ll); s->weight_nonfinite=ll;
    get_param<Long64_t>(&f,"weight_zero",ll); s->weight_zero=ll;
    get_param<Long64_t>(&f,"weight_negative",ll); s->weight_negative=ll;
    get_param<double>(&f,"weight_sum",s->weight_sum); get_param<double>(&f,"weight_sum2",s->weight_sum2);
    get_param<double>(&f,"weight_min",s->weight_min); get_param<double>(&f,"weight_max",s->weight_max);

    get_param<Long64_t>(&f,"cf_all",ll); s->cutflow.all=ll;
    get_param<Long64_t>(&f,"cf_proton",ll); s->cutflow.proton=ll;
    get_param<Long64_t>(&f,"cf_tag_beta",ll); s->cutflow.tag_beta=ll;
    get_param<Long64_t>(&f,"cf_tag_fid",ll); s->cutflow.tag_fid=ll;
    get_param<Long64_t>(&f,"cf_finite_probe",ll); s->cutflow.finite_probe=ll;
    get_param<Long64_t>(&f,"cf_nominal_mass",ll); s->cutflow.nominal_mass=ll;
    get_param<Long64_t>(&f,"cf_fd",ll); s->cutflow.fd=ll;
    get_param<Long64_t>(&f,"cf_ft",ll); s->cutflow.ft=ll;
    get_param<Long64_t>(&f,"cf_ft_projectable",ll); s->cutflow.ft_projectable=ll;
    get_param<Long64_t>(&f,"cf_ft_projected_fid",ll); s->cutflow.ft_projected_fid=ll;
    int ipv=0; get_param<int>(&f,"ft_plane_valid",ipv); s->ft_plane.valid=(ipv!=0);
    get_param<double>(&f,"ft_plane_z",s->ft_plane.z); get_param<Long64_t>(&f,"ft_plane_n",ll); s->ft_plane.n=ll;

    init_sample(*s);
    TDirectory* rd=dynamic_cast<TDirectory*>(f.Get("regions"));
    if (rd) { read_region_state(rd,s->fd); read_region_state(rd,s->ft_all); read_region_state(rd,s->ft_low); read_region_state(rd,s->ft_high); }

    TDirectory* qd=dynamic_cast<TDirectory*>(f.Get("qa"));
    if (qd) {
        if (auto* h=dynamic_cast<TH1D*>(qd->Get("denom_p"))) { s->denom_p.reset(dynamic_cast<TH1D*>(h->Clone())); s->denom_p->SetDirectory(nullptr); }
        if (auto* h=dynamic_cast<TH1D*>(qd->Get("denom_theta"))) { s->denom_theta.reset(dynamic_cast<TH1D*>(h->Clone())); s->denom_theta->SetDirectory(nullptr); }
        if (auto* h=dynamic_cast<TH1D*>(qd->Get("denom_phi"))) { s->denom_phi.reset(dynamic_cast<TH1D*>(h->Clone())); s->denom_phi->SetDirectory(nullptr); }
        if (auto* h=dynamic_cast<TH2D*>(qd->Get("ft_xy_projected"))) { s->ft_xy_projected.reset(dynamic_cast<TH2D*>(h->Clone())); s->ft_xy_projected->SetDirectory(nullptr); }
    }

    TDirectory* wd=dynamic_cast<TDirectory*>(f.Get("windows"));
    if (wd) {
        for (size_t i=0;i<s->windows.size();i++) {
            TDirectory* d=dynamic_cast<TDirectory*>(wd->Get(Form("w%zu",i)));
            if (!d) continue;
            auto& wr=s->windows[i];
            get_param<Long64_t>(d,"fd_denom",ll); wr.fd_denom=ll;
            get_param<Long64_t>(d,"fd_reco",ll); wr.fd_reco=ll;
            get_param<Long64_t>(d,"ft_denom",ll); wr.ft_denom=ll;
            get_param<Long64_t>(d,"ft_reco",ll); wr.ft_reco=ll;
            if (auto* h=dynamic_cast<TH1D*>(d->Get("fd_h"))) { wr.fd_h.reset(dynamic_cast<TH1D*>(h->Clone())); wr.fd_h->SetDirectory(nullptr); }
            if (auto* h=dynamic_cast<TH1D*>(d->Get("ft_h"))) { wr.ft_h.reset(dynamic_cast<TH1D*>(h->Clone())); wr.ft_h->SetDirectory(nullptr); }
            read_fit_state(d,"fd_fit_",wr.fd_fit); read_fit_state(d,"ft_fit_",wr.ft_fit);
        } // endfor
    }
    f.Close();
    return s;
}

struct SampleSpec {
    std::string name,dir;
    bool is_mc=false;
};

static const char* CONCISE_CACHE_VERSION="20260912_concise_v9_ft_fitquality_scan";

std::uint64_t concise_hash(const std::string& s,std::uint64_t h=1469598103934665603ULL) {
    for (unsigned char c:s) {
        h ^= std::uint64_t(c);
        h *= 1099511628211ULL;
    } // endfor
    return h;
}

std::string concise_input_signature(const std::string& dir) {
    void* dp=gSystem->OpenDirectory(dir.c_str());
    if (!dp) return "missing";

    std::vector<std::string> files;
    const char* ent=nullptr;
    while ((ent=gSystem->GetDirEntry(dp))) {
        std::string n(ent);
        if (n.size()>=5 && n.substr(n.size()-5)==".root")
            files.push_back(n);
    } // endwhile
    gSystem->FreeDirectory(dp);
    std::sort(files.begin(),files.end());

    std::uint64_t h=concise_hash(CONCISE_CACHE_VERSION);
    for (const auto& n:files) {
        const std::string p=dir + (dir.empty() || dir.back()=='/' ? "" : "/") + n;
        Long_t id=0,size=0,flags=0,mtime=0;
        std::ostringstream ss;
        ss << n;
        if (gSystem->GetPathInfo(p.c_str(),&id,&size,&flags,&mtime)==0)
            ss << ":" << size << ":" << mtime;
        h=concise_hash(ss.str(),h);
    } // endfor

    std::ostringstream out;
    out << std::hex << h;
    return out.str();
}

std::string concise_cache_path(const SampleSpec& s) {
    // Large worker ROOT files belong on /work, not in the source/home filesystem.
    const std::string d="/work/clas12/thayward/photon_efficiency/cache/concise";
    gSystem->mkdir("/work/clas12/thayward/photon_efficiency/cache",true);
    gSystem->mkdir(d.c_str(),true);
    return d+"/"+s.name+"_"+concise_input_signature(s.dir)+".root";
}

bool concise_valid_root(const std::string& p) {
    if (gSystem->AccessPathName(p.c_str())!=kFALSE) return false;
    TFile f(p.c_str(),"READ");
    const bool ok=!f.IsZombie() && !f.TestBit(TFile::kRecovered);
    f.Close();
    return ok;
}

int analyze_worker(const SampleSpec& spec, const std::string& worker_dir) {
    SampleResult s;
    s.name=spec.name; s.input=spec.dir; s.is_mc=spec.is_mc;
    if (!analyze_sample(s)) return 2;
    const std::string path=worker_dir+"/"+spec.name+".root";
    return save_worker_result(s,path) ? 0 : 3;
}

std::vector<std::unique_ptr<SampleResult>> analyze_samples_parallel(const std::string& out) {
    const SampleSpec all_specs[] = {
        {"data",DATA_DIR,false},
        {"aaogen",AAOGEN_DIR,true},
        {"clasdis",CLASDIS_DIR,true},
        {"dvcsgen",DVCSGEN_DIR,true}
    };

    std::vector<SampleSpec> specs;
    for (const auto& spec : all_specs) {
        if (has_root_files(spec.dir)) {
            std::cout << "[AUTO] Found " << spec.name << " ROOT files in " << spec.dir << "\n";
            specs.push_back(spec);
        } else {
            std::cout << "[AUTO] No " << spec.name << " ROOT files yet; skipping.\n";
        }
    } // endfor

    std::vector<std::unique_ptr<SampleResult>> samples;
    if (specs.empty()) return samples;

    const std::string worker_dir=out+"/.workers";
    gSystem->mkdir(worker_dir.c_str(),true);
    std::cout.flush(); std::cerr.flush();

    struct Child { pid_t pid; SampleSpec spec; };
    std::vector<Child> children;

    for (const auto& spec : specs) {
        const pid_t pid=fork();
        if (pid==0) {
            const int rc=analyze_worker(spec,worker_dir);
            std::cout.flush(); std::cerr.flush();
            _exit(rc);
        } else if (pid>0) {
            children.push_back({pid,spec});
            std::cout << "[PARALLEL] launched " << spec.name << " worker pid=" << pid << "\n";
        } else {
            std::cerr << "WARNING: fork() failed for " << spec.name
                      << "; running that sample sequentially.\n";
            const int rc=analyze_worker(spec,worker_dir);
            if (rc!=0) std::cerr << "WARNING: sequential worker failed for " << spec.name << " rc=" << rc << "\n";
        }
    } // endfor

    for (const auto& child : children) {
        int status=0;
        waitpid(child.pid,&status,0);
        if (!WIFEXITED(status) || WEXITSTATUS(status)!=0) {
            std::cerr << "WARNING: worker " << child.spec.name << " failed";
            if (WIFEXITED(status)) std::cerr << " with rc=" << WEXITSTATUS(status);
            std::cerr << "\n";
        } else {
            std::cout << "[PARALLEL] finished " << child.spec.name << "\n";
        }
    } // endfor

    // Preserve deterministic column order: data, AAOgen, CLASDIS, DVCSgen.
    for (const auto& spec : specs) {
        const std::string path=worker_dir+"/"+spec.name+".root";
        if (gSystem->AccessPathName(path.c_str())!=kFALSE) continue;
        auto s=load_worker_result(path);
        if (s) samples.push_back(std::move(s));
        else std::cerr << "WARNING: could not reload worker result " << path << "\n";
    } // endfor

    return samples;
}


// -----------------------------------------------------------------------------
// Valerii-style FD reproduction.
//
// This is deliberately separated from the stage-1 integrated diagnostics above.
// It follows the workflow documented in Gamma_efficiency_Aug30_2026 as closely
// as the current skim permits:
//   * exact saved 7 x 3 x 6 missing-gamma2 p/theta/phi binning;
//   * FD/PCAL tag requirement and FD probe acceptance;
//   * data and weighted-total-MC Delta-p fits independently in every bin;
//   * nominal weighted MC = AAO + no-exclusivity CLASDIS + DVCS;
//   * unit-count MC component templates with normalization derived from this run;
//   * 1-, 2-, and 3-sigma matching using the fit belonging to data or MC;
//   * epsilon_data, epsilon_MC, and epsilon_data/epsilon_MC correction maps;
//   * data-driven normalization variations from observable-to-observable spread;
//   * normalization-only template shift/smearing nuisances (event-level MC is untouched).
//
// The existing stage-1 nominal Mx(ep) enrichment window is NOT imposed here.
// Valerii's ep-gamma-X denominator is intentionally treated as a mixed sample;
// the weighted MC composition handles the pi0/non-pi0 mixture.
// -----------------------------------------------------------------------------

static const int VAL_NP=7;
static const int VAL_NT=3;
static const int VAL_NPH=6;
static const int VAL_NBIN=VAL_NP*VAL_NT*VAL_NPH;
static const double VAL_P_EDGES[VAL_NP+1] = {0.35,0.50,1.10,1.70,2.30,2.90,3.70,6.00};
static const double VAL_T_EDGES[VAL_NT+1] = {6.0,20.0,27.0,36.0};
static const double VAL_PH_EDGES[VAL_NPH+1] = {-30.0,30.0,90.0,150.0,210.0,270.0,330.0};
static const int VAL_DP_NBIN=60;
static const double VAL_DP_MIN=-1.0;
static const double VAL_DP_MAX= 1.0;
static const int VAL_MIN_FIT_ENTRIES=30;

struct ValNormSet {
    std::string label;
    double aao=1.0;
    double clasdis=1.0;
    double dvcs=1.0;
};

// Historical June-2026 Valerii factors.  These are now comparison references
// only; the active factors are derived from this run's data and MC templates.
static const ValNormSet VAL_HISTORICAL_NOMINAL = {"Valerii_June_nominal",0.307,0.315,1.10};

static const int VAL_COUNT_NBIN=4000;
static const double VAL_COUNT_MIN=-4.0;
static const double VAL_COUNT_MAX= 4.0;

struct ValComponentBin {
    long long denom_rows=0;
    double denom_w=0;
    double denom_w2=0;
    std::unique_ptr<TH1D> residual_fit;
    std::unique_ptr<TH1D> residual_count;
};

enum CoarseRegionIndex {
    CR_FD_LOW=0,
    CR_FD_HIGH=1,
    CR_FT_LOW=2,
    CR_FT_HIGH=3,
    CR_N=4
};

static const char* CR_KEY[CR_N]={"FD_lowE","FD_highE","FT_lowE","FT_highE"};
static const char* CR_LABEL[CR_N]={
    "#splitline{FD}{E_{#gamma,probe}<2 GeV}",
    "#splitline{FD}{E_{#gamma,probe}#geq2 GeV}",
    "#splitline{FT}{E_{#gamma,probe}<2 GeV}",
    "#splitline{FT}{E_{#gamma,probe}#geq2 GeV}"
};

struct CoarseValRegion {
    long long denom_rows=0;
    long long truth_rows=0;
    long long truth_pi0_rows=0;
    std::unique_ptr<TH1D> residual;
};

struct ValComponent {
    std::string name;
    bool is_mc=false;
    long long entries=0;
    std::array<ValComponentBin,VAL_NBIN> bins;
    std::array<CoarseValRegion,CR_N> coarse;

    // Histograms used by the data-driven template-normalization stage.
    // Index order is defined by NormObs below.
    // Normalization QA stages.  "pre" is after the common ep-gamma-X/FD
    // baseline but before the June-2026 exclusivity cuts.  "nminus1" applies
    // all normalization cuts except the cut on the plotted quantity itself.
    // "full" applies all normalization cuts.  The low/high-E histograms use
    // the N-1 selection plus E_gamma<2 or >3 GeV and feed the template fits.
    std::vector<std::unique_ptr<TH1D>> norm_pre;
    // Sequential Step-1 diagnostic: distributions after the M_X^2(ep) cut
    // but before M_X^2(e gamma), coplanarity, and angle(gamma,X).
    std::vector<std::unique_ptr<TH1D>> norm_after_mx2ep;
    // Sequential Step-1D diagnostic: after M_X^2(ep) and M_X^2(e gamma),
    // before Trento coplanarity and angle(gamma,X).
    std::vector<std::unique_ptr<TH1D>> norm_after_mx2ep_mx2eg;
    // Sequential Step-1E diagnostic: after M_X^2(ep), M_X^2(e gamma),
    // and Trento coplanarity, before angle(gamma,X).
    std::vector<std::unique_ptr<TH1D>> norm_after_mx2ep_mx2eg_dphi;
    std::vector<std::unique_ptr<TH1D>> norm_nminus1;
    std::vector<std::unique_ptr<TH1D>> norm_full;
    std::vector<std::unique_ptr<TH1D>> norm_lowE;
    std::vector<std::unique_ptr<TH1D>> norm_highE;
    std::array<long long,6> norm_cutflow{{0,0,0,0,0,0}};

    // Independent FT-probe exclusivity/normalization sample.  The observed tag
    // photon remains FD, exactly as in the tag-and-probe workflow; the inferred
    // probe must project into the FT fiducial annulus.
    std::vector<std::unique_ptr<TH1D>> norm_ft_pre;
    std::vector<std::unique_ptr<TH1D>> norm_ft_after_mx2ep;
    std::vector<std::unique_ptr<TH1D>> norm_ft_after_mx2ep_mx2eg;
    std::vector<std::unique_ptr<TH1D>> norm_ft_after_mx2ep_mx2eg_dphi;
    std::vector<std::unique_ptr<TH1D>> norm_ft_nminus1;
    std::vector<std::unique_ptr<TH1D>> norm_ft_full;
    std::vector<std::unique_ptr<TH1D>> norm_ft_lowE;
    std::vector<std::unique_ptr<TH1D>> norm_ft_highE;
    std::array<long long,6> norm_ft_cutflow{{0,0,0,0,0,0}};

    // Baseline correlations with inferred probe energy, plus the sequential
    // probe-energy spectra used to check whether tighter cuts sculpt coverage.
    std::vector<std::unique_ptr<TH2D>> corr_fd;
    std::vector<std::unique_ptr<TH2D>> corr_ft;
    std::vector<std::unique_ptr<TH1D>> eprobe_stage_fd;
    std::vector<std::unique_ptr<TH1D>> eprobe_stage_ft;

    // FT N-1 distributions used for the exclusivity cut optimization.
    // x = inferred probe energy, y = scanned exclusivity variable.
    std::vector<std::unique_ptr<TH2D>> cutscan_ft;

    // FT coplanarity scan inputs for normalization-fit quality studies.
    // Each histogram stores the low-/high-tag-energy normalization observables
    // after all nominal FT cuts except that coplanarity is allowed to vary.
    std::vector<std::vector<std::unique_ptr<TH2D>>> fitqual_ft_low;
    std::vector<std::vector<std::unique_ptr<TH2D>>> fitqual_ft_high;

    bool normalization_branches_complete=false;
};

enum NormObs {
    NORM_MX2_EP=0,
    NORM_MX2_EPG=1,
    NORM_ANGLE_GX=2,
    NORM_ANGLE_EX=3,
    NORM_EGAMMA=4,
    NORM_MX2_EG=5,

    // Azimuthal QA plus the active zero-centered Trento coplanarity residual.
    // Raw lab phi and the individual/unshifted Trento quantities remain QA-only.
    NORM_DPHI_PG_RAW=6,
    NORM_PHI_P_TRENTO=7,
    NORM_PHI_G_TRENTO=8,
    NORM_DPHI_TRENTO=9,
    NORM_DPHI_TRENTO_SHIFT180=10,
    NORM_DELTA_T_PG=11,

    NORM_NOBS=12
};

struct NormObsDef {
    const char* key;
    const char* title;
    int nb; double lo,hi;
    bool use_low;
    bool use_high;
};

static const NormObsDef NORM_OBS[NORM_NOBS] = {
    {"Mx2_ep",     "M_{X}^{2}(ep);M_{X}^{2}(ep) (GeV^{2});Candidates",          90,-0.50,1.30,true, false},
    {"Mx2_epg",    "M_{X}^{2}(ep#gamma);M_{X}^{2}(ep#gamma) (GeV^{2});Candidates",80,-0.25,0.25,true, true},
    // Valerii slide 11 normalization observable.
    {"angle_gX",   "#angle(#gamma,X);#angle(#gamma,X) (deg);Candidates",       75,0.0,30.0,true, true},
    // Valerii slide 9 exclusivity cut.  This is intentionally NOT used as a
    // normalization-fit observable.
    {"angle_eX",   "#angle(e,X);#angle(e,X) (deg);Candidates",                 75,0.0,30.0,false,false},
    {"Egamma",     "Tag-photon energy;E_{#gamma} (GeV);Candidates",             95,0.4,8.0,true, false},
    {"Mx2_eg",     "M_{X}^{2}(e#gamma);M_{X}^{2}(e#gamma) (GeV^{2});Candidates",100,0.0,6.0,true, true},

    // Raw lab and individual/unshifted Trento quantities are QA-only.
    // The final zero-centered Trento coplanarity residual is also used as the
    // active exclusivity cut and as a high-E DVCS normalization diagnostic.
    {"dphi_pg_raw", "Raw lab azimuth difference;wrap(#phi_{p}^{lab}-#phi_{#gamma}^{lab}) (deg);Candidates",
                     180,-180.0,180.0,false,false},
    {"phi_p_trento", "Proton Trento azimuth;#phi_{p}^{Trento} (deg);Candidates",
                     180,-180.0,180.0,false,false},
    {"phi_g_trento", "Tag-photon Trento azimuth;#phi_{#gamma}^{Trento} (deg);Candidates",
                     180,-180.0,180.0,false,false},
    {"dphi_trento", "Trento azimuth difference;wrap(#phi_{p}^{Trento}-#phi_{#gamma}^{Trento}) (deg);Candidates",
                     180,-180.0,180.0,false,false},
    {"dphi_trento_shift180", "Trento coplanarity residual;#delta#phi_{copl} (deg);Candidates",
                     120,-30.0,30.0,false,true},
    {"delta_t_pg", "#Delta t(p,#gamma);t_{p}-t_{#gamma} (GeV^{2});Candidates",
                     120,-2.0,2.0,false,true}
};

// Exclusivity-vs-probe-energy diagnostics. These are filled at the detector-
// accepted ep-gamma-X baseline, BEFORE any active exclusivity cut, so cut
// optimization can be checked directly against the probe-energy coverage used
// by the efficiency map.
static const int CORR_NOBS=7;
static const int CORR_OBS[CORR_NOBS] = {
    NORM_MX2_EP,
    NORM_MX2_EPG,
    NORM_MX2_EG,
    NORM_DPHI_TRENTO_SHIFT180,
    NORM_ANGLE_GX,
    NORM_ANGLE_EX,
    NORM_DELTA_T_PG
};
static const int EPROBE_NBIN=40;
static const double EPROBE_MIN=0.40;
static const double EPROBE_MAX=8.40;
static const int EPROBE_NSTAGE=5;
static const char* EPROBE_STAGE_KEY[EPROBE_NSTAGE] = {
    "baseline","mx2_ep","mx2_eg","coplanarity","angle_gX"
};

// FT cut optimization.  Each scan is N-1: all of the other nominal
// exclusivity requirements are imposed while the listed variable is varied.
// Mx2(epgamma) is an additional tightening inside the loose skim window.
static const int CUTSCAN_NVAR=4;
static const int CUTSCAN_OBS[CUTSCAN_NVAR] = {
    NORM_MX2_EPG,
    NORM_DPHI_TRENTO_SHIFT180,
    NORM_MX2_EP,
    NORM_ANGLE_GX
};
static const char* CUTSCAN_KEY[CUTSCAN_NVAR] = {
    "Mx2_epg_halfwidth",
    "dphi_absmax",
    "Mx2_ep_halfwidth",
    "angle_gX_max"
};
static const char* CUTSCAN_XTITLE[CUTSCAN_NVAR] = {
    "|M_{X}^{2}(ep#gamma)| maximum (GeV^{2})",
    "|#Delta#phi_{copl}| maximum (deg)",
    "M_{X}^{2}(ep) half-width (GeV^{2})",
    "angle(#gamma,X) maximum (deg)"
};
static const int CUTSCAN_NPOINT=25;
static const double CUTSCAN_MX2EP_CENTER=0.039; // center of [-0.231,0.309] GeV^2

static const int FT_COPL_SCAN_N=9;
static const double FT_COPL_SCAN_VALUES[FT_COPL_SCAN_N] = {
    4.0,5.0,5.7,6.5,7.5,8.5,10.0,12.0,15.0
};


// June-2026 Valerii normalization/exclusivity selection.
// Slide 9 labels the final angular cut as Angle(e,X)<9.2 deg, but the literal
// scattered-electron/missing-probe opening angle removes essentially the whole
// reconstructed sample.  Until that shorthand is resolved, the active Step-1E
// selection is angle(gamma,X)<9.2 deg, which gives the expected smooth exclusive
// behavior and is also the angular observable used in Valerii's normalization plots.
// angle(e,X) remains available as QA only.
static const double NORM_MX2_EP_MIN=-0.231;
static const double NORM_MX2_EP_MAX= 0.309;
static const double NORM_MX2_EG_MIN= 1.4;
// The Trento/exclusive-plane QA shows the exclusive topology peaking at zero
// after shifting the proton-photon Trento-angle difference by 180 degrees.
// Use Valerii's quoted |Delta phi(p,gamma)| < 5.7 deg requirement on that
// zero-centered coplanarity residual.
static const double NORM_DPHI_TRENTO_MAX=5.7;
static const double NORM_ANGLE_GX_MAX=9.2;

// Template morphing deliberately follows the philosophy used in the DVCS
// exclusivity-selection suite: allow the reconstructed-MC template to shift and
// acquire extra Gaussian resolution when comparing it with data.  Here one
// COMMON morph is applied to all MC components for a given observable.  This
// avoids giving AAO/CLASDIS/DVCS independent shape freedom that could absorb
// the very component fractions the normalization fit is supposed to measure.
static const double NORM_MAX_SHIFT_BINS=6.0;
static const double NORM_MAX_SMEAR_BINS=10.0;
static const double NORM_MORPH_STEP_BINS=0.5;

struct NormFitPoint {
    std::string observable;
    bool valid=false;
    double aao=0,aao_err=0;
    double clasdis=0,clasdis_err=0;
    double dvcs=0,dvcs_err=0;
    double chi2=0; int ndf=0;
    // Common MC-template morph used only by the normalization fit.  The shift
    // and extra Gaussian width are in physical units of the plotted observable.
    // They are never applied to the event-level efficiency numerator/denominator.
    double morph_shift=0;
    double morph_sigma=0;
    double raw_chi2=0;
    int raw_ndf=0;
};

struct NormDerivation {
    bool valid=false;
    bool used_fallback=false;
    ValNormSet nominal;
    ValNormSet low;
    ValNormSet high;
    std::vector<NormFitPoint> low_points;
    std::vector<NormFitPoint> high_points;

    // Common low-E AAOgen/CLASDIS normalization across several observables.
    // Per-observable fits above remain QA only.
    NormFitPoint low_global;

    // FT-only sensitivity test: omit CLASDIS and renormalize AAOgen + DVCSgen.
    // This does not replace the nominal three-component extraction.
    bool no_clasdis_valid=false;
    ValNormSet no_clasdis;
};

struct ValEval {
    bool valid=false;
    double denom=0, denom_w2=0;
    double num=0, num_w2=0;
    double efficiency=-1, efficiency_err=0;
};

struct ValBinResult {
    FitResult data_fit;
    FitResult mc_fit;
    std::array<ValEval,3> data;
    std::array<ValEval,3> mc;
    std::array<double,3> correction{{-1,-1,-1}};
    std::array<double,3> correction_stat_err{{0,0,0}};
    double correction_norm_set1=-1;
    double correction_norm_set2=-1;
    double norm_systematic=0;
    double matching_systematic=0;
    double partial_total_unc=0;
    bool nominal_valid=false;
};

int val_axis_bin(double x, const double* e, int n) {
    for (int i=0;i<n;i++) {
        const bool last=(i==n-1);
        if (x>=e[i] && (x<e[i+1] || (last && x<=e[i+1]))) return i;
    } // endfor
    return -1;
}

int val_flat_bin(double p,double th,double ph) {
    const int ip=val_axis_bin(p,VAL_P_EDGES,VAL_NP);
    const int it=val_axis_bin(th,VAL_T_EDGES,VAL_NT);
    const double wph=wrap_phi(ph);
    const int iph=val_axis_bin(wph,VAL_PH_EDGES,VAL_NPH);
    if (ip<0 || it<0 || iph<0) return -1;
    return (ip*VAL_NT+it)*VAL_NPH+iph;
}

void val_unflatten(int ib,int& ip,int& it,int& iph) {
    iph=ib%VAL_NPH;
    const int q=ib/VAL_NPH;
    it=q%VAL_NT;
    ip=q/VAL_NT;
}

double val_component_scale(const std::string& name,const ValNormSet& n) {
    if (name=="aaogen") return n.aao;
    if (name=="clasdis") return n.clasdis;
    if (name=="dvcsgen") return n.dvcs;
    return 1.0;
}

FitResult fit_valerii_residual(TH1D* h) {
    FitResult r;
    if (!h) { r.reason="null histogram"; return r; }
    r.candidates=static_cast<long long>(h->GetEntries());
    if (h->GetEntries()<VAL_MIN_FIT_ENTRIES) {
        r.reason="too few candidates";
        return r;
    }

    const double peak=h->GetBinCenter(h->GetMaximumBin());
    const double seed_lo=std::max(VAL_DP_MIN,peak-0.20);
    const double seed_hi=std::min(VAL_DP_MAX,peak+0.20);
    TF1 seed("val_seed_tmp","gaus",seed_lo,seed_hi);
    seed.SetParameters(std::max(1.0,h->GetMaximum()),peak,0.10);
    seed.SetParLimits(2,0.010,0.500);
    const int seed_status=h->Fit(&seed,"QNR");

    double mu=peak;
    double sg=0.10;
    if (seed_status==0) {
        mu=seed.GetParameter(1);
        sg=std::fabs(seed.GetParameter(2));
    }
    if (!(sg>0.010 && sg<0.500)) sg=0.10;

    // The documented Valerii fit range is the full [-1,1] GeV residual range.
    // Keep the same Gaussian + linear-background form already validated in the
    // stage-1 diagnostic until the original production fitter source is imported.
    TF1 f("val_fit_tmp","gaus(0)+pol1(3)",VAL_DP_MIN,VAL_DP_MAX);
    // ROOT emits a noisy ParameterSettings warning when SetParameters is called
    // with a seed already outside a subsequently imposed bound.  Clamp the
    // seeds first, then install the limits, so every per-bin fit starts from a
    // legal point without changing the allowed fit region.
    mu=std::max(-0.749,std::min(0.749,mu));
    sg=std::max(0.0101,std::min(0.499,sg));
    f.SetParLimits(1,-0.75,0.75);
    f.SetParLimits(2,0.010,0.500);
    f.SetParameters(std::max(1.0,h->GetMaximum()),mu,sg,
                    std::max(0.0,h->GetBinContent(1)),0.0);
    r.root_status=h->Fit(&f,"QNR");
    r.amplitude=f.GetParameter(0);
    r.mean=f.GetParameter(1);
    r.mean_err=f.GetParError(1);
    r.sigma=std::fabs(f.GetParameter(2));
    r.sigma_err=f.GetParError(2);
    r.bg0=f.GetParameter(3);
    r.bg1=f.GetParameter(4);
    r.fit_lo=VAL_DP_MIN;
    r.fit_hi=VAL_DP_MAX;
    r.chi2=f.GetChisquare();
    r.ndf=f.GetNDF();

    if (r.root_status!=0) { r.reason="ROOT fit status != 0"; return r; }
    if (r.mean<=-0.74 || r.mean>=0.74) { r.reason="mean at/near fit boundary"; return r; }
    if (!(r.sigma>0.0101 && r.sigma<0.495)) { r.reason="sigma at/near fit boundary"; return r; }
    if (!(r.sigma_err>0) || r.sigma_err/r.sigma>0.50) { r.reason="sigma uncertainty too large"; return r; }
    if (!(r.mean_err>=0) || r.mean_err>0.15) { r.reason="mean uncertainty too large"; return r; }
    if (r.ndf<=0) { r.reason="non-positive NDF"; return r; }
    r.valid=true;
    r.reason="valid";
    return r;
}


std::unique_ptr<TH1D> make_norm_hist(int io,const char* region,const std::string& sample) {
    const auto& d=NORM_OBS[io];
    std::unique_ptr<TH1D> h(new TH1D(Form("norm_%s_%s_%s",region,d.key,sample.c_str()),d.title,d.nb,d.lo,d.hi));
    h->Sumw2(); h->SetDirectory(nullptr);
    return h;
}

std::unique_ptr<TH2D> make_corr_hist(int io,const char* detector,const std::string& sample) {
    const auto& d=NORM_OBS[io];
    const int ny=std::min(90,d.nb);
    std::unique_ptr<TH2D> h(new TH2D(
        Form("corr_%s_%s_%s",detector,d.key,sample.c_str()),
        Form(";E_{#gamma,probe} (GeV);%s",d.key),
        EPROBE_NBIN,EPROBE_MIN,EPROBE_MAX,ny,d.lo,d.hi));
    h->Sumw2();
    h->SetDirectory(nullptr);
    return h;
}

std::unique_ptr<TH1D> make_eprobe_stage_hist(const char* detector,int stage,
                                             const std::string& sample) {
    std::unique_ptr<TH1D> h(new TH1D(
        Form("eprobe_%s_%s_%s",detector,EPROBE_STAGE_KEY[stage],sample.c_str()),
        ";E_{#gamma,probe} (GeV);Candidates",
        EPROBE_NBIN,EPROBE_MIN,EPROBE_MAX));
    h->Sumw2();
    h->SetDirectory(nullptr);
    return h;
}

std::unique_ptr<TH2D> make_cutscan_hist(int iv,const std::string& sample) {
    const int io=CUTSCAN_OBS[iv];
    const auto& d=NORM_OBS[io];
    const int ny=std::min(120,d.nb);
    std::unique_ptr<TH2D> h(new TH2D(
        Form("cutscan_ft_%s_%s",CUTSCAN_KEY[iv],sample.c_str()),
        Form(";E_{#gamma,probe} (GeV);%s",d.key),
        EPROBE_NBIN,EPROBE_MIN,EPROBE_MAX,ny,d.lo,d.hi));
    h->Sumw2();
    h->SetDirectory(nullptr);
    return h;
}

double invariant_m2_from_epg(const Branches& b) {
    if (!b.have_beam_energy || !b.have_e_kin || !b.have_tag_corr_kin) return std::numeric_limits<double>::quiet_NaN();
    const double me=0.00051099895, mp=0.9382720813;
    const double et=b.e_p>0?std::sqrt(b.e_p*b.e_p+me*me):0.0;
    const double th=deg2rad(b.e_theta), ph=deg2rad(b.e_phi);
    const double ex=b.e_p*std::sin(th)*std::cos(ph), ey=b.e_p*std::sin(th)*std::sin(ph), ez=b.e_p*std::cos(th);
    const double gt=deg2rad(b.tag_corr_theta), gp=deg2rad(b.tag_corr_phi);
    const double gx=b.tag_corr_p*std::sin(gt)*std::cos(gp), gy=b.tag_corr_p*std::sin(gt)*std::sin(gp), gz=b.tag_corr_p*std::cos(gt);
    const double ebeam=b.beam_energy;
    const double pbeam=std::sqrt(std::max(0.0,ebeam*ebeam-me*me));
    const double px=-(ex+gx), py=-(ey+gy), pz=pbeam-ez-gz;
    const double E=ebeam+mp-et-b.tag_corr_p;
    return E*E-px*px-py*py-pz*pz;
}


double wrap180(double x) {
    while (x<=-180.0) x+=360.0;
    while (x> 180.0) x-=360.0;
    return x;
}

bool trento_phi_deg(const Branches& b,
                    double h_p,double h_theta_deg,double h_phi_deg,
                    double& phi_deg) {
    if (!b.have_beam_energy || !b.have_e_kin) return false;
    if (!(b.beam_energy>0) || !(b.e_p>0) || !(h_p>0)) return false;

    const double me=0.00051099895;
    const double pbeam=std::sqrt(std::max(0.0,b.beam_energy*b.beam_energy-me*me));

    // Beam and scattered-electron three-vectors in the CLAS/lab frame.
    TVector3 k(0.0,0.0,pbeam);
    const double eth=deg2rad(b.e_theta), eph=deg2rad(b.e_phi);
    TVector3 kp(b.e_p*std::sin(eth)*std::cos(eph),
                b.e_p*std::sin(eth)*std::sin(eph),
                b.e_p*std::cos(eth));

    TVector3 q=k-kp;
    if (q.Mag2()<=0) return false;
    TVector3 qhat=q.Unit();

    // Hadron/photon vector whose Trento azimuth is requested.
    const double hth=deg2rad(h_theta_deg), hph=deg2rad(h_phi_deg);
    TVector3 h(h_p*std::sin(hth)*std::cos(hph),
               h_p*std::sin(hth)*std::sin(hph),
               h_p*std::cos(hth));

    // Plane normals.  Only directions matter.
    TVector3 n_lep=q.Cross(k);
    TVector3 n_had=q.Cross(h);
    if (n_lep.Mag2()<=1e-18 || n_had.Mag2()<=1e-18) return false;
    n_lep=n_lep.Unit();
    n_had=n_had.Unit();

    // Signed angle from the lepton plane to the q-h plane about +q.
    // This is the standard Trento-style oriented-plane construction; the
    // overall sign convention is less important here than using the same
    // definition for proton and photon before taking their difference.
    const double c=std::max(-1.0,std::min(1.0,n_lep.Dot(n_had)));
    const double s=qhat.Dot(n_lep.Cross(n_had));
    phi_deg=wrap180(std::atan2(s,c)*180.0/M_PI);
    return std::isfinite(phi_deg);
}


double norm_delta_t_pg(const Branches& b) {
    if (!b.have_beam_energy || !b.have_e_kin || !b.have_p_corr_kin || !b.have_tag_corr_kin)
        return std::numeric_limits<double>::quiet_NaN();

    const double mp=0.9382720813;
    const double me=0.00051099895;

    const double Eb=b.beam_energy;
    const double pb=std::sqrt(std::max(0.0,Eb*Eb-me*me));
    TLorentzVector k(0.0,0.0,pb,Eb);
    TLorentzVector target(0.0,0.0,0.0,mp);

    const double eth=deg2rad(b.e_theta), eph=deg2rad(b.e_phi);
    TLorentzVector kp;
    kp.SetPxPyPzE(b.e_p*std::sin(eth)*std::cos(eph),
                  b.e_p*std::sin(eth)*std::sin(eph),
                  b.e_p*std::cos(eth),
                  std::sqrt(b.e_p*b.e_p+me*me));

    const double pth=deg2rad(b.p_corr_theta), pph=deg2rad(b.p_corr_phi);
    TLorentzVector pp;
    pp.SetPxPyPzE(b.p_corr_p*std::sin(pth)*std::cos(pph),
                  b.p_corr_p*std::sin(pth)*std::sin(pph),
                  b.p_corr_p*std::cos(pth),
                  std::sqrt(b.p_corr_p*b.p_corr_p+mp*mp));

    const double gth=deg2rad(b.tag_corr_theta), gph=deg2rad(b.tag_corr_phi);
    TLorentzVector g;
    g.SetPxPyPzE(b.tag_corr_p*std::sin(gth)*std::cos(gph),
                 b.tag_corr_p*std::sin(gth)*std::sin(gph),
                 b.tag_corr_p*std::cos(gth),
                 b.tag_corr_p);

    const TLorentzVector q=k-kp;

    const double tp=(target-pp).M2();
    const double tg=(q-g).M2();
    const double dt=tp-tg;

    return std::isfinite(dt) ? dt : std::numeric_limits<double>::quiet_NaN();
}

double norm_observable_value(const Branches& b,int io) {
    if (io==NORM_MX2_EP) return b.have_Mx2_ep?b.Mx2_ep:(b.Mx_ep>=0?b.Mx_ep*b.Mx_ep:-b.Mx_ep*b.Mx_ep);
    if (io==NORM_MX2_EPG) return b.have_Mx2_epg_corr?b.Mx2_epg_corr:std::numeric_limits<double>::quiet_NaN();
    if (io==NORM_ANGLE_GX) {
        if (!b.have_tag_corr_kin) return std::numeric_limits<double>::quiet_NaN();
        return opening_angle_deg(b.tag_corr_theta,b.tag_corr_phi,b.probe_corr_theta,b.probe_corr_phi);
    }
    if (io==NORM_ANGLE_EX) {
        if (!b.have_e_kin) return std::numeric_limits<double>::quiet_NaN();
        return opening_angle_deg(b.e_theta,b.e_phi,b.probe_corr_theta,b.probe_corr_phi);
    }
    if (io==NORM_EGAMMA) return b.have_tag_corr_kin?b.tag_corr_p:std::numeric_limits<double>::quiet_NaN();
    if (io==NORM_MX2_EG) return invariant_m2_from_epg(b);
    if (io==NORM_DPHI_PG_RAW) {
        if (!b.have_p_corr_kin || !b.have_tag_corr_kin)
            return std::numeric_limits<double>::quiet_NaN();
        return wrap180(b.p_corr_phi-b.tag_corr_phi);
    }

    if (io==NORM_DELTA_T_PG) {
        return norm_delta_t_pg(b);
    }

    if (io==NORM_PHI_P_TRENTO || io==NORM_PHI_G_TRENTO ||
        io==NORM_DPHI_TRENTO || io==NORM_DPHI_TRENTO_SHIFT180) {
        if (!b.have_p_corr_kin || !b.have_tag_corr_kin)
            return std::numeric_limits<double>::quiet_NaN();

        double phip=0.0, phig=0.0;
        if (!trento_phi_deg(b,b.p_corr_p,b.p_corr_theta,b.p_corr_phi,phip))
            return std::numeric_limits<double>::quiet_NaN();
        if (!trento_phi_deg(b,b.tag_corr_p,b.tag_corr_theta,b.tag_corr_phi,phig))
            return std::numeric_limits<double>::quiet_NaN();

        if (io==NORM_PHI_P_TRENTO) return phip;
        if (io==NORM_PHI_G_TRENTO) return phig;

        const double d=wrap180(phip-phig);
        if (io==NORM_DPHI_TRENTO) return d;

        // Two particles in the same exclusive hadronic plane can differ by
        // either ~0 or ~180 degrees depending on the particle/convention.
        // This residual explicitly tests the 180-degree interpretation.
        const double sign=(d>=0.0 ? 1.0 : -1.0);
        return sign*(180.0-std::fabs(d));
    }
    return std::numeric_limits<double>::quiet_NaN();
}


struct NormCutFlags {
    bool finite=false;
    bool mx2_ep=false;
    bool mx2_eg=false;
    bool dphi_trento=false;
    bool angle_gX=false;
    bool all=false;
};

NormCutFlags norm_cut_flags(const Branches& b) {
    NormCutFlags f;
    const double mx2ep=norm_observable_value(b,NORM_MX2_EP);
    const double mx2eg=norm_observable_value(b,NORM_MX2_EG);
    const double dphi =norm_observable_value(b,NORM_DPHI_TRENTO_SHIFT180);
    const double ang  =norm_observable_value(b,NORM_ANGLE_GX);

    f.finite=std::isfinite(mx2ep) && std::isfinite(mx2eg) &&
             std::isfinite(dphi) && std::isfinite(ang);
    if (!f.finite) return f;

    f.mx2_ep=(mx2ep>NORM_MX2_EP_MIN && mx2ep<NORM_MX2_EP_MAX);
    f.mx2_eg=(mx2eg>NORM_MX2_EG_MIN);
    f.dphi_trento=(std::fabs(dphi)<NORM_DPHI_TRENTO_MAX);
    f.angle_gX=(ang<NORM_ANGLE_GX_MAX);

    f.all=f.mx2_ep && f.mx2_eg && f.dphi_trento && f.angle_gX;
    return f;
}

bool norm_pass_nminus1(const NormCutFlags& f,int io) {
    if (!f.finite) return false;
    if (io!=NORM_MX2_EP && !f.mx2_ep) return false;
    if (io!=NORM_MX2_EG && !f.mx2_eg) return false;

    // When plotting/fitting the coplanarity residual itself, omit its own hard
    // cut.  All other observables see the active |delta phi_copl| < 5.7 deg cut.
    if (io!=NORM_DPHI_TRENTO_SHIFT180 && !f.dphi_trento) return false;

    if (io!=NORM_ANGLE_GX && !f.angle_gX) return false;
    return true;
}

std::unique_ptr<TH1D> morph_norm_hist(const TH1D* src,double shift,double sigma,const char* name) {
    if (!src) return nullptr;
    std::unique_ptr<TH1D> out((TH1D*)src->Clone(name));
    // Clone() carries the source histogram's Sumw2 array with it.  Calling
    // Sumw2() again here produces one ROOT warning for every morphed template.
    out->Reset("ICES");
    out->SetDirectory(nullptr);
    if (out->GetSumw2N()==0) out->Sumw2();
    const int nb=src->GetNbinsX();
    if (sigma<=1e-12) {
        for (int j=1;j<=nb;j++) {
            const double w=src->GetBinContent(j);
            const double e=src->GetBinError(j);
            if (w==0 && e==0) continue;
            const double x=src->GetBinCenter(j)+shift;
            const int i=out->FindBin(x);
            if (i<1 || i>nb) continue;
            out->SetBinContent(i,out->GetBinContent(i)+w);
            const double oe=out->GetBinError(i);
            out->SetBinError(i,std::sqrt(oe*oe+e*e));
        }
        const double sinteg=src->Integral(), ointeg=out->Integral();
        if (sinteg>0 && ointeg>0) out->Scale(sinteg/ointeg);
        return out;
    }
    const double inv=1.0/(std::sqrt(2.0)*sigma);
    for (int j=1;j<=nb;j++) {
        const double w=src->GetBinContent(j);
        const double e=src->GetBinError(j);
        if (w==0 && e==0) continue;
        const double mu=src->GetBinCenter(j)+shift;
        for (int i=1;i<=nb;i++) {
            const double lo=out->GetXaxis()->GetBinLowEdge(i);
            const double hi=out->GetXaxis()->GetBinUpEdge(i);
            const double prob=0.5*(std::erf((hi-mu)*inv)-std::erf((lo-mu)*inv));
            if (!(prob>0)) continue;
            out->SetBinContent(i,out->GetBinContent(i)+w*prob);
            const double oe=out->GetBinError(i);
            out->SetBinError(i,std::sqrt(oe*oe+e*e*prob*prob));
        }
    }
    const double sinteg=src->Integral(), ointeg=out->Integral();
    if (sinteg>0 && ointeg>0) out->Scale(sinteg/ointeg);
    return out;
}

NormFitPoint fit_two_templates_fixed_shapes(const TH1D* data,const TH1D* aao,const TH1D* cls,const std::string& key) {
    NormFitPoint r; r.observable=key;
    if (!data || !aao || !cls) return r;
    double saa=0,sbb=0,sab=0,sad=0,sbd=0; int used=0;
    for (int i=1;i<=data->GetNbinsX();i++) {
        const double d=data->GetBinContent(i), a=aao->GetBinContent(i), b=cls->GetBinContent(i);
        if (d<=0 && a<=0 && b<=0) continue;
        const double w=1.0/std::max(1.0,d);
        saa+=w*a*a; sbb+=w*b*b; sab+=w*a*b; sad+=w*a*d; sbd+=w*b*d; used++;
    }
    const double det=saa*sbb-sab*sab;
    if (used<5 || !(det>0)) return r;
    double A=(sad*sbb-sbd*sab)/det;
    double B=(sbd*saa-sad*sab)/det;
    if (A<0) { A=0; B=(sbb>0?sbd/sbb:0); }
    if (B<0) { B=0; A=(saa>0?sad/saa:0); }
    if (!(A>=0 && B>=0) || !std::isfinite(A) || !std::isfinite(B)) return r;
    double chi2=0; int n=0;
    for (int i=1;i<=data->GetNbinsX();i++) {
        const double d=data->GetBinContent(i), a=aao->GetBinContent(i), b=cls->GetBinContent(i);
        if (d<=0 && a<=0 && b<=0) continue;
        const double q=d-A*a-B*b; chi2+=q*q/std::max(1.0,d); n++;
    }
    r.aao=A; r.clasdis=B; r.aao_err=(det>0?std::sqrt(sbb/det):0); r.clasdis_err=(det>0?std::sqrt(saa/det):0);
    r.chi2=chi2; r.ndf=std::max(0,n-2); r.valid=(r.ndf>0);
    return r;
}

NormFitPoint fit_dvcs_fixed_shapes(const TH1D* data,const TH1D* aao,const TH1D* cls,const TH1D* dvc,
                                   double A,double B,const std::string& key) {
    NormFitPoint r; r.observable=key; r.aao=A; r.clasdis=B;
    if (!data || !aao || !cls || !dvc) return r;
    double sxx=0,sxy=0; int used=0;
    for (int i=1;i<=data->GetNbinsX();i++) {
        const double d=data->GetBinContent(i), x=dvc->GetBinContent(i);
        const double fixed=A*aao->GetBinContent(i)+B*cls->GetBinContent(i);
        if (d<=0 && x<=0 && fixed<=0) continue;
        const double w=1.0/std::max(1.0,d);
        sxx+=w*x*x; sxy+=w*x*(d-fixed); used++;
    }
    if (used<5 || !(sxx>0)) return r;
    double C=sxy/sxx; if (C<0) C=0;
    double chi2=0; int n=0;
    for (int i=1;i<=data->GetNbinsX();i++) {
        const double d=data->GetBinContent(i);
        const double m=A*aao->GetBinContent(i)+B*cls->GetBinContent(i)+C*dvc->GetBinContent(i);
        if (d<=0 && m<=0) continue;
        const double q=d-m; chi2+=q*q/std::max(1.0,d); n++;
    }
    r.dvcs=C; r.dvcs_err=std::sqrt(1.0/sxx); r.chi2=chi2; r.ndf=std::max(0,n-1); r.valid=(r.ndf>0 && std::isfinite(C));
    return r;
}

NormFitPoint fit_two_templates_morphed(const TH1D* data,const TH1D* aao,const TH1D* cls,const std::string& key) {
    NormFitPoint raw=fit_two_templates_fixed_shapes(data,aao,cls,key);
    NormFitPoint best; best.observable=key; best.raw_chi2=raw.chi2; best.raw_ndf=raw.ndf;
    double bestchi=std::numeric_limits<double>::infinity();
    if (!data || !aao || !cls) return best;
    const double bw=data->GetXaxis()->GetBinWidth(1);
    for (double sb=-NORM_MAX_SHIFT_BINS; sb<=NORM_MAX_SHIFT_BINS+1e-9; sb+=NORM_MORPH_STEP_BINS) {
        const double shift=sb*bw;
        for (double wb=0; wb<=NORM_MAX_SMEAR_BINS+1e-9; wb+=NORM_MORPH_STEP_BINS) {
            const double sigma=wb*bw;
            auto ma=morph_norm_hist(aao,shift,sigma,"tmp_morph_aao");
            auto mc=morph_norm_hist(cls,shift,sigma,"tmp_morph_cls");
            auto q=fit_two_templates_fixed_shapes(data,ma.get(),mc.get(),key);
            if (!q.valid || !(q.chi2<bestchi)) continue;
            bestchi=q.chi2; best=q; best.morph_shift=shift; best.morph_sigma=sigma;
            best.raw_chi2=raw.chi2; best.raw_ndf=raw.ndf;
        }
    }
    return best;
}

NormFitPoint fit_dvcs_template_morphed(const TH1D* data,const TH1D* aao,const TH1D* cls,const TH1D* dvc,
                                       double A,double B,const std::string& key) {
    NormFitPoint raw=fit_dvcs_fixed_shapes(data,aao,cls,dvc,A,B,key);
    NormFitPoint best; best.observable=key; best.aao=A; best.clasdis=B; best.raw_chi2=raw.chi2; best.raw_ndf=raw.ndf;
    double bestchi=std::numeric_limits<double>::infinity();
    if (!data || !aao || !cls || !dvc) return best;
    const double bw=data->GetXaxis()->GetBinWidth(1);
    for (double sb=-NORM_MAX_SHIFT_BINS; sb<=NORM_MAX_SHIFT_BINS+1e-9; sb+=NORM_MORPH_STEP_BINS) {
        const double shift=sb*bw;
        for (double wb=0; wb<=NORM_MAX_SMEAR_BINS+1e-9; wb+=NORM_MORPH_STEP_BINS) {
            const double sigma=wb*bw;
            auto ma=morph_norm_hist(aao,shift,sigma,"tmp_high_aao");
            auto mc=morph_norm_hist(cls,shift,sigma,"tmp_high_cls");
            auto mv=morph_norm_hist(dvc,shift,sigma,"tmp_high_dvc");
            auto q=fit_dvcs_fixed_shapes(data,ma.get(),mc.get(),mv.get(),A,B,key);
            if (!q.valid || !(q.chi2<bestchi)) continue;
            bestchi=q.chi2; best=q; best.morph_shift=shift; best.morph_sigma=sigma;
            best.raw_chi2=raw.chi2; best.raw_ndf=raw.ndf;
        }
    }
    return best;
}

// Normalization fits are defined above.  The active versions include a common
// shift + extra Gaussian smearing scan of the MC templates; the fixed-shape
// versions are retained internally for raw-versus-morphed diagnostics.

double mean_valid(const std::vector<NormFitPoint>& v,int which) {
    double s=0; int n=0;
    for (const auto& x:v) if (x.valid) { double q=(which==0?x.aao:(which==1?x.clasdis:x.dvcs)); if (q>=0 && std::isfinite(q)) {s+=q;n++;} }
    return n?s/n:0;
}

std::pair<double,double> range_valid(const std::vector<NormFitPoint>& v,int which,double fallback) {
    double lo=1e300,hi=-1e300; int n=0;
    for (const auto& x:v) if (x.valid) { double q=(which==0?x.aao:(which==1?x.clasdis:x.dvcs)); if (q>=0 && std::isfinite(q)) {lo=std::min(lo,q);hi=std::max(hi,q);n++;} }
    return n?std::make_pair(lo,hi):std::make_pair(fallback,fallback);
}

void style_norm_component(TH1D* h,int color,int width=2) {
    if (!h) return; h->SetLineColor(color); h->SetLineWidth(width); h->SetFillStyle(0); h->SetStats(0);
}

void draw_norm_fit(const TH1D* hd,const TH1D* ha,const TH1D* hc,const TH1D* hv,
                   double A,double B,double C,const std::string& title,const std::string& file,
                   double chi2=-1,int ndf=0,double morph_shift=0.0,double morph_sigma=0.0) {
    if (!hd || !ha || !hc || !hv) return;
    std::unique_ptr<TH1D> d((TH1D*)hd->Clone("norm_draw_data")); d->SetDirectory(nullptr);
    auto a=morph_norm_hist(ha,morph_shift,morph_sigma,"norm_draw_aao"); a->Scale(A);
    auto c=morph_norm_hist(hc,morph_shift,morph_sigma,"norm_draw_cls"); c->Scale(B);
    auto v=morph_norm_hist(hv,morph_shift,morph_sigma,"norm_draw_dvc"); v->Scale(C);
    std::unique_ptr<TH1D> tot((TH1D*)a->Clone("norm_draw_total")); tot->Add(c.get()); tot->Add(v.get()); tot->SetDirectory(nullptr);
    std::unique_ptr<TH1D> pull((TH1D*)d->Clone("norm_draw_pull")); pull->Reset(); pull->SetDirectory(nullptr);
    for (int i=1;i<=pull->GetNbinsX();i++) {
        const double den=std::sqrt(std::max(1.0,d->GetBinContent(i)));
        pull->SetBinContent(i,(d->GetBinContent(i)-tot->GetBinContent(i))/den);
    }
    d->SetMarkerStyle(20); d->SetMarkerSize(0.7); d->SetLineColor(kBlack);
    style_norm_component(a.get(),kRed+1); style_norm_component(c.get(),kOrange+7); style_norm_component(v.get(),kGreen+2); style_norm_component(tot.get(),kBlue+1,3);
    TCanvas can("c_norm_fit","",1100,850);
    TPad top("norm_top","",0,0.27,1,1); TPad bot("norm_bot","",0,0,1,0.27);
    top.SetBottomMargin(0.02); bot.SetTopMargin(0.03); bot.SetBottomMargin(0.33); top.Draw(); bot.Draw();
    top.cd();
    d->SetTitle(title.c_str()); d->GetXaxis()->SetLabelSize(0); d->SetMaximum(1.25*std::max(d->GetMaximum(),tot->GetMaximum()));
    d->Draw("E1"); a->Draw("HIST SAME"); c->Draw("HIST SAME"); v->Draw("HIST SAME"); tot->Draw("HIST SAME"); d->Draw("E1 SAME");
    TLegend leg(0.58,0.63,0.88,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(d.get(),"Data","lep"); leg.AddEntry(a.get(),Form("AAO #times %.4g",A),"l");
    leg.AddEntry(c.get(),Form("CLASDIS #times %.4g",B),"l"); leg.AddEntry(v.get(),Form("DVCSgen #times %.4g",C),"l"); leg.AddEntry(tot.get(),"Total MC","l"); leg.Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
    if (ndf>0) tx.DrawLatex(0.15,0.86,Form("#chi^{2}/ndf = %.1f/%d = %.2f",chi2,ndf,chi2/ndf));
    if (std::fabs(morph_shift)>1e-12 || morph_sigma>1e-12)
        tx.DrawLatex(0.15,0.81,Form("common MC morph: shift=%.4g, #sigma_{add}=%.4g",morph_shift,morph_sigma));
    bot.cd(); pull->SetTitle(""); pull->GetYaxis()->SetTitle("Pull"); pull->GetYaxis()->SetNdivisions(505); pull->GetYaxis()->SetTitleSize(0.10); pull->GetYaxis()->SetLabelSize(0.08); pull->GetYaxis()->SetTitleOffset(0.45);
    pull->GetXaxis()->SetTitle(hd->GetXaxis()->GetTitle()); pull->GetXaxis()->SetTitleSize(0.12); pull->GetXaxis()->SetLabelSize(0.09); pull->SetMinimum(-5); pull->SetMaximum(5); pull->Draw("HIST");
    TLine z(pull->GetXaxis()->GetXmin(),0,pull->GetXaxis()->GetXmax(),0); z.SetLineStyle(2); z.Draw();
    can.SaveAs(file.c_str());
}

void draw_norm_shape_overlay(const TH1D* hd,const TH1D* ha,const TH1D* hc,const TH1D* hv,
                             const std::string& title,const std::string& file) {
    if (!hd || !ha || !hc || !hv) return;
    std::unique_ptr<TH1D> d((TH1D*)hd->Clone("shape_data")); d->SetDirectory(nullptr);
    std::unique_ptr<TH1D> a((TH1D*)ha->Clone("shape_aao")); a->SetDirectory(nullptr);
    std::unique_ptr<TH1D> c((TH1D*)hc->Clone("shape_cls")); c->SetDirectory(nullptr);
    std::unique_ptr<TH1D> v((TH1D*)hv->Clone("shape_dvc")); v->SetDirectory(nullptr);
    auto unit=[](TH1D* h){ const double q=h?h->Integral():0; if (h && q>0) h->Scale(1.0/q); };
    unit(d.get()); unit(a.get()); unit(c.get()); unit(v.get());
    d->SetMarkerStyle(20); d->SetMarkerSize(0.65); d->SetLineColor(kBlack); d->SetStats(0);
    style_norm_component(a.get(),kRed+1); style_norm_component(c.get(),kOrange+7); style_norm_component(v.get(),kGreen+2);
    TCanvas can("c_norm_shape","",1050,760);
    d->SetTitle(title.c_str()); d->GetYaxis()->SetTitle("Unit-area candidates");
    d->SetMaximum(1.25*std::max({d->GetMaximum(),a->GetMaximum(),c->GetMaximum(),v->GetMaximum()}));
    d->Draw("E1"); a->Draw("HIST SAME"); c->Draw("HIST SAME"); v->Draw("HIST SAME"); d->Draw("E1 SAME");
    TLegend leg(0.62,0.68,0.88,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(d.get(),"Data","lep"); leg.AddEntry(a.get(),"AAO (unit area)","l");
    leg.AddEntry(c.get(),"CLASDIS (unit area)","l"); leg.AddEntry(v.get(),"DVCSgen (unit area)","l"); leg.Draw();
    can.SaveAs(file.c_str());
}

void draw_norm_cutflow(const std::vector<std::unique_ptr<ValComponent>>& vv,const std::string& file) {
    const char* labs[6]={"baseline","Mx2(ep)","Mx2(e#gamma)","|#delta#phi_{copl}|<5.7^{#circ}","angle(e,X)","all cuts"};
    TCanvas c("c_norm_cutflow","",1150,760);
    TLegend leg(0.68,0.68,0.90,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0);
    std::vector<std::unique_ptr<TH1D>> keep;
    bool first=true; int idx=0;
    const int colors[4]={kBlack,kRed+1,kOrange+7,kGreen+2};
    for (const auto& vp:vv) {
        if (!vp) continue;
        std::unique_ptr<TH1D> h(new TH1D(Form("cf_%d",idx),"Normalization/exclusivity cut flow;Selection stage;Survival fraction",6,0,6));
        h->SetDirectory(nullptr); h->SetStats(0); h->SetLineWidth(3); h->SetLineColor(colors[std::min(idx,3)]);
        h->SetMarkerColor(colors[std::min(idx,3)]); h->SetMarkerStyle(20+idx);
        const double n0=std::max(1LL,vp->norm_cutflow[0]);
        for (int i=0;i<6;i++) { h->GetXaxis()->SetBinLabel(i+1,labs[i]); h->SetBinContent(i+1,vp->norm_cutflow[i]/n0); }
        h->SetMinimum(0); h->SetMaximum(1.08); h->GetXaxis()->LabelsOption("v");
        h->Draw(first?"HIST P":"HIST P SAME"); first=false; leg.AddEntry(h.get(),vp->name.c_str(),"lp");
        keep.push_back(std::move(h)); idx++;
    }
    leg.Draw(); c.SetBottomMargin(0.25); c.SaveAs(file.c_str());
}

bool analyze_val_component_worker(const SampleSpec& spec,const std::string& path) {
    const std::string pattern=make_pattern(spec.dir);
    TChain c("PhotonEfficiency");
    const int nf=c.Add(pattern.c_str());
    if (nf<=0 || c.GetEntries()<=0) return false;

    Branches b;
    b.reset_arrays();
    if (!attach(c,b)) return false;
    c.SetCacheSize(64LL*1024LL*1024LL);
    c.AddBranchToCache("*",kTRUE);

    std::array<long long,VAL_NBIN> rows{};
    std::array<double,VAL_NBIN> sumw{};
    std::array<double,VAL_NBIN> sumw2{};

    std::array<long long,CR_N> coarse_rows{{0,0,0,0}};
    std::array<long long,CR_N> coarse_truth_rows{{0,0,0,0}};
    std::array<long long,CR_N> coarse_truth_pi0{{0,0,0,0}};
    std::array<std::unique_ptr<TH1D>,CR_N> coarse_h;
    for (int ir=0;ir<CR_N;ir++) {
        coarse_h[ir].reset(new TH1D(
            Form("coarse_%s",CR_KEY[ir]),
            Form(";%s;Candidates","#Delta p_{#gamma2}=p_{rec}-p_{miss} (GeV)"),
            VAL_COUNT_NBIN,VAL_COUNT_MIN,VAL_COUNT_MAX));
        coarse_h[ir]->Sumw2();
        coarse_h[ir]->SetDirectory(nullptr);
    } // endfor

    const FTPlaneEstimate ft_plane=estimate_ft_plane(c,b);

    // Do NOT write one TTree row per reconstructed candidate.  With millions of
    // skim rows that temporary representation can become multi-GB and exhaust
    // the ifarm home/output filesystem.  The analysis only needs the residual
    // distributions and their weighted sums, so persist fixed-size histograms.
    std::array<std::unique_ptr<TH1D>,VAL_NBIN> hfit;
    std::array<std::unique_ptr<TH1D>,VAL_NBIN> hcount;
    std::vector<std::unique_ptr<TH1D>> norm_pre, norm_after_mx2ep, norm_after_mx2ep_mx2eg,
        norm_after_mx2ep_mx2eg_dphi, norm_nminus1, norm_full, norm_low, norm_high;
    std::vector<std::unique_ptr<TH1D>> norm_ft_pre, norm_ft_after_mx2ep, norm_ft_after_mx2ep_mx2eg,
        norm_ft_after_mx2ep_mx2eg_dphi, norm_ft_nminus1, norm_ft_full, norm_ft_low, norm_ft_high;

    for (int io=0;io<NORM_NOBS;io++) {
        norm_pre.push_back(make_norm_hist(io,"fd_pre",spec.name));
        norm_after_mx2ep.push_back(make_norm_hist(io,"fd_after_mx2ep",spec.name));
        norm_after_mx2ep_mx2eg.push_back(make_norm_hist(io,"fd_after_mx2ep_mx2eg",spec.name));
        norm_after_mx2ep_mx2eg_dphi.push_back(make_norm_hist(io,"fd_after_mx2ep_mx2eg_dphi",spec.name));
        norm_nminus1.push_back(make_norm_hist(io,"fd_nminus1",spec.name));
        norm_full.push_back(make_norm_hist(io,"fd_full",spec.name));
        norm_low.push_back(make_norm_hist(io,"fd_lowE",spec.name));
        norm_high.push_back(make_norm_hist(io,"fd_highE",spec.name));

        norm_ft_pre.push_back(make_norm_hist(io,"ft_pre",spec.name));
        norm_ft_after_mx2ep.push_back(make_norm_hist(io,"ft_after_mx2ep",spec.name));
        norm_ft_after_mx2ep_mx2eg.push_back(make_norm_hist(io,"ft_after_mx2ep_mx2eg",spec.name));
        norm_ft_after_mx2ep_mx2eg_dphi.push_back(make_norm_hist(io,"ft_after_mx2ep_mx2eg_dphi",spec.name));
        norm_ft_nminus1.push_back(make_norm_hist(io,"ft_nminus1",spec.name));
        norm_ft_full.push_back(make_norm_hist(io,"ft_full",spec.name));
        norm_ft_low.push_back(make_norm_hist(io,"ft_lowE",spec.name));
        norm_ft_high.push_back(make_norm_hist(io,"ft_highE",spec.name));
    } // endfor
    std::vector<std::unique_ptr<TH2D>> corr_fd,corr_ft;
    for (int jc=0;jc<CORR_NOBS;jc++) {
        const int io=CORR_OBS[jc];
        corr_fd.push_back(make_corr_hist(io,"FD",spec.name));
        corr_ft.push_back(make_corr_hist(io,"FT",spec.name));
    } // endfor

    std::vector<std::unique_ptr<TH1D>> eprobe_stage_fd,eprobe_stage_ft;
    for (int is=0;is<EPROBE_NSTAGE;is++) {
        eprobe_stage_fd.push_back(make_eprobe_stage_hist("FD",is,spec.name));
        eprobe_stage_ft.push_back(make_eprobe_stage_hist("FT",is,spec.name));
    } // endfor

    std::vector<std::unique_ptr<TH2D>> cutscan_ft;
    for (int iv=0;iv<CUTSCAN_NVAR;iv++)
        cutscan_ft.push_back(make_cutscan_hist(iv,spec.name));

    std::vector<std::vector<std::unique_ptr<TH2D>>> fitqual_ft_low(FT_COPL_SCAN_N);
    std::vector<std::vector<std::unique_ptr<TH2D>>> fitqual_ft_high(FT_COPL_SCAN_N);
    for (int ic=0;ic<FT_COPL_SCAN_N;ic++) {
        for (int io=0;io<NORM_NOBS;io++) {
            const auto& od=NORM_OBS[io];
            fitqual_ft_low[ic].push_back(std::make_unique<TH2D>(
                Form("fitqual_low_c%02d_%s_%s",ic,od.key,spec.name.c_str()),
                Form(";E_{#gamma,probe} (GeV);%s",od.key),
                EPROBE_NBIN,EPROBE_MIN,EPROBE_MAX,
                std::min(120,od.nb),od.lo,od.hi));
            fitqual_ft_high[ic].push_back(std::make_unique<TH2D>(
                Form("fitqual_high_c%02d_%s_%s",ic,od.key,spec.name.c_str()),
                Form(";E_{#gamma,probe} (GeV);%s",od.key),
                EPROBE_NBIN,EPROBE_MIN,EPROBE_MAX,
                std::min(120,od.nb),od.lo,od.hi));
            fitqual_ft_low[ic].back()->Sumw2();
            fitqual_ft_high[ic].back()->Sumw2();
            fitqual_ft_low[ic].back()->SetDirectory(nullptr);
            fitqual_ft_high[ic].back()->SetDirectory(nullptr);
        } // endfor
    } // endfor

    std::array<long long,6> norm_cutflow{{0,0,0,0,0,0}};
    std::array<long long,6> norm_ft_cutflow{{0,0,0,0,0,0}};
    for (int ib=0;ib<VAL_NBIN;ib++) {
        hfit[ib].reset(new TH1D(Form("fit_b%03d",ib),";#Delta p_{#gamma2} (GeV);weighted candidates",
                                VAL_DP_NBIN,VAL_DP_MIN,VAL_DP_MAX));
        hcount[ib].reset(new TH1D(Form("count_b%03d",ib),";#Delta p_{#gamma2} (GeV);weighted candidates",
                                  VAL_COUNT_NBIN,VAL_COUNT_MIN,VAL_COUNT_MAX));
        hfit[ib]->Sumw2();
        hcount[ib]->Sumw2();
        hfit[ib]->SetDirectory(nullptr);
        hcount[ib]->SetDirectory(nullptr);
    } // endfor

    long long selected=0,reco=0,outside_count_range=0;
    for (Long64_t i=0;i<c.GetEntries();i++) {
        c.GetEntry(i);

        // Match the FD/PCAL workflow: the observed tag photon must itself be FD.
        // The mixed ep-gamma-X denominator is NOT restricted by the stage-1 pi0
        // parent window.
        if (!b.p_pass_standard) continue;
        if (!b.tag_pass_beta || !b.tag_pass_fiducial) continue;
        if (b.tag_detector!=1) continue;
        if (!finite_good(b.probe_corr_p) || !finite_good(b.probe_corr_theta) ||
            !finite_good(b.probe_corr_phi)) continue;

        // Compute the exclusivity selection BEFORE applying the FD-only
        // Valerii 7x3x6 binning.  This is essential for the coarse FT study:
        // the missing probe can lie in the FT even though the observed tag is FD.
        const double Eg=b.have_tag_corr_kin?b.tag_corr_p:std::numeric_limits<double>::quiet_NaN();
        const NormCutFlags ncf=norm_cut_flags(b);

        // Detector-specific inferred-probe acceptance.
        const bool probe_fd=(b.probe_corr_p>=PROBE_P_MIN &&
                             b.probe_corr_theta>=FD_THETA_MIN &&
                             b.probe_corr_theta<=FD_THETA_MAX);

        bool probe_ft=false;
        FTProjection ftp;
        if (b.probe_corr_p>=PROBE_P_MIN &&
            b.probe_corr_theta>=FT_THETA_MIN &&
            b.probe_corr_theta<=FT_THETA_MAX) {
            ftp=project_ft(b,ft_plane);
            probe_ft=ftp.valid && ftp.fiducial;
        } // endif

        auto fill_corr_sample = [&](bool accept,
                                    std::vector<std::unique_ptr<TH2D>>& corr,
                                    std::vector<std::unique_ptr<TH1D>>& stages) {
            if (!accept || !std::isfinite(b.probe_corr_p)) return;

            // Baseline correlations: no exclusivity requirement has yet been applied.
            for (int jc=0;jc<CORR_NOBS;jc++) {
                const double y=norm_observable_value(b,CORR_OBS[jc]);
                if (std::isfinite(y)) corr[jc]->Fill(b.probe_corr_p,y);
            } // endfor

            stages[0]->Fill(b.probe_corr_p);
            if (ncf.mx2_ep) {
                stages[1]->Fill(b.probe_corr_p);
                if (ncf.mx2_eg) {
                    stages[2]->Fill(b.probe_corr_p);
                    if (ncf.dphi_trento) {
                        stages[3]->Fill(b.probe_corr_p);
                        if (ncf.angle_gX) stages[4]->Fill(b.probe_corr_p);
                    } // endif
                } // endif
            } // endif
        };

        fill_corr_sample(probe_fd,corr_fd,eprobe_stage_fd);
        fill_corr_sample(probe_ft,corr_ft,eprobe_stage_ft);

        // FT N-1 cut-scan inputs.  These retain the Eprobe coordinate so every
        // prospective cut can be evaluated for both purity and energy sculpting.
        if (probe_ft) {
            for (int iv=0;iv<CUTSCAN_NVAR;iv++) {
                const int io=CUTSCAN_OBS[iv];
                bool pass_other=false;

                if (io==NORM_MX2_EPG) {
                    // Mx2(epgamma) is not one of the active June cuts; scan it
                    // as an additional tightening after the full nominal set.
                    pass_other=ncf.all;
                } else if (io==NORM_MX2_EP) {
                    pass_other=ncf.mx2_eg && ncf.dphi_trento && ncf.angle_gX;
                } else if (io==NORM_DPHI_TRENTO_SHIFT180) {
                    pass_other=ncf.mx2_ep && ncf.mx2_eg && ncf.angle_gX;
                } else if (io==NORM_ANGLE_GX) {
                    pass_other=ncf.mx2_ep && ncf.mx2_eg && ncf.dphi_trento;
                } // endif

                if (!pass_other) continue;
                const double y=norm_observable_value(b,io);
                if (std::isfinite(y))
                    cutscan_ft[iv]->Fill(b.probe_corr_p,y);
            } // endfor

            // Fit-quality scan: vary only coplanarity while retaining the
            // nominal Mx2(ep), Mx2(e#gamma), and angle(gamma,X) selections.
            // This directly tests the user's priority: maximize FT statistics
            // only where the AAO/DVCS normalization remains well constrained.
            if (ncf.mx2_ep && ncf.mx2_eg && ncf.angle_gX) {
                const double dphi=std::fabs(norm_observable_value(
                    b,NORM_DPHI_TRENTO_SHIFT180));
                for (int ic=0;ic<FT_COPL_SCAN_N;ic++) {
                    if (!(dphi<FT_COPL_SCAN_VALUES[ic])) continue;
                    const bool low_tag=(b.tag_corr_p<2.0);
                    const bool high_tag=(b.tag_corr_p>3.0);
                    if (!low_tag && !high_tag) continue;

                    for (int io=0;io<NORM_NOBS;io++) {
                        const double y=norm_observable_value(b,io);
                        if (!std::isfinite(y)) continue;
                        if (low_tag)
                            fitqual_ft_low[ic][io]->Fill(b.probe_corr_p,y);
                        if (high_tag)
                            fitqual_ft_high[ic][io]->Fill(b.probe_corr_p,y);
                    } // endfor
                } // endfor
            } // endif
        } // endif

        auto fill_norm_sample = [&](bool accept,
                                    std::vector<std::unique_ptr<TH1D>>& pre,
                                    std::vector<std::unique_ptr<TH1D>>& after_ep,
                                    std::vector<std::unique_ptr<TH1D>>& after_e_g,
                                    std::vector<std::unique_ptr<TH1D>>& after_dphi,
                                    std::vector<std::unique_ptr<TH1D>>& nminus1,
                                    std::vector<std::unique_ptr<TH1D>>& full,
                                    std::vector<std::unique_ptr<TH1D>>& low,
                                    std::vector<std::unique_ptr<TH1D>>& high,
                                    std::array<long long,6>& cutflow) {
            if (!accept) return;

            cutflow[0]++;
            if (ncf.mx2_ep) {
                cutflow[1]++;
                if (ncf.mx2_eg) {
                    cutflow[2]++;
                    if (ncf.dphi_trento) {
                        cutflow[3]++;
                        if (ncf.angle_gX) cutflow[4]++;
                    } // endif
                } // endif
            } // endif
            if (ncf.all) cutflow[5]++;

            for (int io=0;io<NORM_NOBS;io++) {
                const double x=norm_observable_value(b,io);
                if (!std::isfinite(x)) continue;

                pre[io]->Fill(x);
                if (ncf.mx2_ep) after_ep[io]->Fill(x);
                if (ncf.mx2_ep && ncf.mx2_eg) after_e_g[io]->Fill(x);
                if (ncf.mx2_ep && ncf.mx2_eg && ncf.dphi_trento)
                    after_dphi[io]->Fill(x);

                if (norm_pass_nminus1(ncf,io)) {
                    nminus1[io]->Fill(x);
                    // IMPORTANT: Valerii's normalization control regions use
                    // the OBSERVED TAG-photon energy, not the inferred probe.
                    if (std::isfinite(Eg) && Eg<2.0) low[io]->Fill(x);
                    if (std::isfinite(Eg) && Eg>3.0) high[io]->Fill(x);
                } // endif

                if (ncf.all) full[io]->Fill(x);
            } // endfor
        };

        // The same exclusivity sequence is evaluated independently for the
        // FD-probe and FT-probe populations.
        fill_norm_sample(probe_fd,
                         norm_pre,norm_after_mx2ep,norm_after_mx2ep_mx2eg,
                         norm_after_mx2ep_mx2eg_dphi,norm_nminus1,norm_full,
                         norm_low,norm_high,norm_cutflow);

        fill_norm_sample(probe_ft,
                         norm_ft_pre,norm_ft_after_mx2ep,norm_ft_after_mx2ep_mx2eg,
                         norm_ft_after_mx2ep_mx2eg_dphi,norm_ft_nminus1,norm_ft_full,
                         norm_ft_low,norm_ft_high,norm_ft_cutflow);

        // Coarse integrated efficiency regions use the corresponding detector
        // acceptance and the PROBE-energy split requested for the final result.
        if (ncf.all) {
            int cr=-1;
            int probe_detector=-1;
            if (probe_fd) {
                cr=(b.probe_corr_p<2.0 ? CR_FD_LOW : CR_FD_HIGH);
                probe_detector=1;
            } else if (probe_ft) {
                cr=(b.probe_corr_p<2.0 ? CR_FT_LOW : CR_FT_HIGH);
                probe_detector=0;
            } // endif

            if (cr>=0) {
                coarse_rows[cr]++;
                if (spec.is_mc && b.have_truth) {
                    coarse_truth_rows[cr]++;
                    if (b.truth_probe_pid==22 && b.truth_probe_parent==111)
                        coarse_truth_pi0[cr]++;
                } // endif

                const int kc=best_probe_candidate(b,probe_detector);
                if (kc>=0) {
                    const double dpc=b.neutral_p[kc]-b.probe_corr_p;
                    if (std::isfinite(dpc)) coarse_h[cr]->Fill(dpc);
                } // endif
            } // endif
        } // endif

        // Legacy detailed Valerii 7x3x6 output remains FD-only and is retained
        // internally for future granularity studies.  It no longer controls the
        // concise detector-integrated FD/FT normalization.
        const int ib=val_flat_bin(b.probe_corr_p,b.probe_corr_theta,b.probe_corr_phi);
        if (ib<0) continue;

        if (!ncf.all) continue;

        // Valerii FD reproduction: use UNIT event weights inside each MC sample.
        // The AAO/CLASDIS/DVCS relative normalizations are applied only when the
        // component histograms are combined below.  In particular, do not use
        // the skim's MC::Event.weight branch here: for CLASDIS that quantity is
        // not an event-statistics weight and would catastrophically distort the
        // mixture.
        const double tw=1.0;

        rows[ib]++;
        sumw[ib]+=tw;
        sumw2[ib]+=tw*tw;
        selected++;

        const int k=best_probe_candidate(b,1);
        if (k<0) continue;
        const double dp=b.neutral_p[k]-b.probe_corr_p;
        if (!std::isfinite(dp)) continue;
        hfit[ib]->Fill(dp,tw);
        hcount[ib]->Fill(dp,tw);
        if (dp<VAL_COUNT_MIN || dp>VAL_COUNT_MAX) outside_count_range++;
        reco++;
    } // endfor

    TFile f(path.c_str(),"RECREATE");
    if (f.IsZombie()) return false;
    TNamed nm("sample_name",spec.name.c_str()); nm.Write();
    put_param<int>(&f,"is_mc",spec.is_mc?1:0);
    put_param<Long64_t>(&f,"entries",c.GetEntries());
    put_param<Long64_t>(&f,"selected_denominator_rows",selected);
    put_param<Long64_t>(&f,"reconstructed_candidate_rows",reco);
    put_param<Long64_t>(&f,"residual_rows_outside_count_range",outside_count_range);

    TTree denominators("denominators","Valerii FD denominator sums by analysis bin");
    Int_t dbin=-1;
    Long64_t drows=0;
    Double_t dsumw=0,dsumw2=0;
    denominators.Branch("bin",&dbin,"bin/I");
    denominators.Branch("rows",&drows,"rows/L");
    denominators.Branch("sumw",&dsumw,"sumw/D");
    denominators.Branch("sumw2",&dsumw2,"sumw2/D");
    for (dbin=0;dbin<VAL_NBIN;dbin++) {
        drows=rows[dbin]; dsumw=sumw[dbin]; dsumw2=sumw2[dbin];
        denominators.Fill();
    } // endfor
    denominators.Write();

    TDirectory* cd=f.mkdir("coarse");
    if (cd) {
        cd->cd();
        for (int ir=0;ir<CR_N;ir++) {
            TDirectory* d=cd->mkdir(CR_KEY[ir]);
            d->cd();
            put_param<Long64_t>(d,"denom_rows",coarse_rows[ir]);
            put_param<Long64_t>(d,"truth_rows",coarse_truth_rows[ir]);
            put_param<Long64_t>(d,"truth_pi0_rows",coarse_truth_pi0[ir]);
            coarse_h[ir]->Write("residual");
            cd->cd();
        } // endfor
        f.cd();
    } // endif

    auto write_corr_dir = [&](const char* dirname,
                              const std::vector<std::unique_ptr<TH2D>>& corr,
                              const std::vector<std::unique_ptr<TH1D>>& stages) {
        TDirectory* d=f.mkdir(dirname);
        if (!d) return;
        d->cd();
        for (int jc=0;jc<CORR_NOBS;jc++)
            if (corr[jc]) corr[jc]->Write(Form("corr_%s",NORM_OBS[CORR_OBS[jc]].key));
        for (int is=0;is<EPROBE_NSTAGE;is++)
            if (stages[is]) stages[is]->Write(Form("eprobe_%s",EPROBE_STAGE_KEY[is]));
        f.cd();
    };

    write_corr_dir("correlations_FD",corr_fd,eprobe_stage_fd);
    write_corr_dir("correlations_FT",corr_ft,eprobe_stage_ft);

    {
        TDirectory* sd=f.mkdir("cutscan_FT");
        if (sd) {
            sd->cd();
            for (int iv=0;iv<CUTSCAN_NVAR;iv++)
                if (cutscan_ft[iv]) cutscan_ft[iv]->Write(CUTSCAN_KEY[iv]);
            f.cd();
        } // endif
    }

    {
        TDirectory* qd=f.mkdir("fitquality_FT");
        if (qd) {
            for (int ic=0;ic<FT_COPL_SCAN_N;ic++) {
                TDirectory* cd=qd->mkdir(Form("copl_%02d",ic));
                if (!cd) continue;
                TDirectory* ld=cd->mkdir("low");
                TDirectory* hd=cd->mkdir("high");
                if (ld) {
                    ld->cd();
                    for (int io=0;io<NORM_NOBS;io++)
                        fitqual_ft_low[ic][io]->Write(NORM_OBS[io].key);
                } // endif
                if (hd) {
                    hd->cd();
                    for (int io=0;io<NORM_NOBS;io++)
                        fitqual_ft_high[ic][io]->Write(NORM_OBS[io].key);
                } // endif
                f.cd();
            } // endfor
        } // endif
    }

    TDirectory* rd=f.mkdir("residuals");
    if (!rd) { f.Close(); return false; }
    rd->cd();
    for (int ib=0;ib<VAL_NBIN;ib++) {
        hfit[ib]->Write(Form("fit_b%03d",ib));
        hcount[ib]->Write(Form("count_b%03d",ib));
    } // endfor
    f.cd();
    auto write_norm_dir = [&](const char* dirname,
                              std::vector<std::unique_ptr<TH1D>>& pre,
                              std::vector<std::unique_ptr<TH1D>>& after_ep,
                              std::vector<std::unique_ptr<TH1D>>& after_e_g,
                              std::vector<std::unique_ptr<TH1D>>& after_dphi,
                              std::vector<std::unique_ptr<TH1D>>& nminus1,
                              std::vector<std::unique_ptr<TH1D>>& full,
                              std::vector<std::unique_ptr<TH1D>>& low,
                              std::vector<std::unique_ptr<TH1D>>& high) {
        TDirectory* nd=f.mkdir(dirname);
        if (!nd) return;
        nd->cd();
        for (int io=0;io<NORM_NOBS;io++) {
            pre[io]->Write(Form("pre_%s",NORM_OBS[io].key));
            after_ep[io]->Write(Form("after_mx2ep_%s",NORM_OBS[io].key));
            after_e_g[io]->Write(Form("after_mx2ep_mx2eg_%s",NORM_OBS[io].key));
            after_dphi[io]->Write(Form("after_mx2ep_mx2eg_dphi_%s",NORM_OBS[io].key));
            nminus1[io]->Write(Form("nminus1_%s",NORM_OBS[io].key));
            full[io]->Write(Form("full_%s",NORM_OBS[io].key));
            low[io]->Write(Form("lowE_%s",NORM_OBS[io].key));
            high[io]->Write(Form("highE_%s",NORM_OBS[io].key));
        } // endfor
        f.cd();
    };

    write_norm_dir("normalization_FD",
                   norm_pre,norm_after_mx2ep,norm_after_mx2ep_mx2eg,
                   norm_after_mx2ep_mx2eg_dphi,norm_nminus1,norm_full,
                   norm_low,norm_high);

    write_norm_dir("normalization_FT",
                   norm_ft_pre,norm_ft_after_mx2ep,norm_ft_after_mx2ep_mx2eg,
                   norm_ft_after_mx2ep_mx2eg_dphi,norm_ft_nminus1,norm_ft_full,
                   norm_ft_low,norm_ft_high);

    for (int ic=0;ic<6;ic++) {
        put_param<Long64_t>(&f,Form("norm_fd_cutflow_%d",ic),norm_cutflow[ic]);
        put_param<Long64_t>(&f,Form("norm_ft_cutflow_%d",ic),norm_ft_cutflow[ic]);
    } // endfor
    put_param<int>(&f,"normalization_branches_complete",
        (b.have_tag_corr_kin && b.have_Mx2_epg_corr && b.have_beam_energy && b.have_e_kin && b.have_p_corr_kin)?1:0);
    const Int_t nwrite=f.Write();
    const bool root_write_error=f.TestBit(TFile::kWriteError);
    f.Close();

    if (nwrite<=0 || root_write_error || gSystem->AccessPathName(path.c_str())) {
        std::cerr << "ERROR: failed to persist worker cache file " << path << "\n";
        gSystem->Unlink(path.c_str());
        return false;
    }

    // Reopen once before reporting success.  This catches truncated files from
    // quota/full-filesystem failures immediately instead of poisoning a rerun.
    {
        TFile check(path.c_str(),"READ");
        if (check.IsZombie() || check.TestBit(TFile::kRecovered)) {
            std::cerr << "ERROR: cache file is corrupt/incomplete: " << path << "\n";
            check.Close();
            gSystem->Unlink(path.c_str());
            return false;
        }
        check.Close();
    }

    return true;
}
std::unique_ptr<ValComponent> load_val_component(const std::string& path) {
    TFile f(path.c_str(),"READ");
    if (f.IsZombie()) return nullptr;
    auto* nn=dynamic_cast<TNamed*>(f.Get("sample_name"));
    if (!nn) return nullptr;

    std::unique_ptr<ValComponent> v(new ValComponent);
    v->name=nn->GetTitle();
    int imc=0; get_param<int>(&f,"is_mc",imc); v->is_mc=(imc!=0);
    Long64_t ent=0; get_param<Long64_t>(&f,"entries",ent); v->entries=ent;

    auto* den=dynamic_cast<TTree*>(f.Get("denominators"));
    if (!den) return nullptr;
    Int_t ib=-1; Long64_t rows=0; Double_t sw=0,sw2=0;
    den->SetBranchAddress("bin",&ib);
    den->SetBranchAddress("rows",&rows);
    den->SetBranchAddress("sumw",&sw);
    den->SetBranchAddress("sumw2",&sw2);
    for (Long64_t i=0;i<den->GetEntries();i++) {
        den->GetEntry(i);
        if (ib<0 || ib>=VAL_NBIN) continue;
        v->bins[ib].denom_rows=rows;
        v->bins[ib].denom_w=sw;
        v->bins[ib].denom_w2=sw2;
    } // endfor

    auto* cd=dynamic_cast<TDirectory*>(f.Get("coarse"));
    if (cd) {
        for (int ir=0;ir<CR_N;ir++) {
            TDirectory* d=dynamic_cast<TDirectory*>(cd->Get(CR_KEY[ir]));
            if (!d) continue;
            Long64_t q=0;
            get_param<Long64_t>(d,"denom_rows",q); v->coarse[ir].denom_rows=q;
            get_param<Long64_t>(d,"truth_rows",q); v->coarse[ir].truth_rows=q;
            get_param<Long64_t>(d,"truth_pi0_rows",q); v->coarse[ir].truth_pi0_rows=q;
            if (auto* h=dynamic_cast<TH1D*>(d->Get("residual"))) {
                v->coarse[ir].residual.reset(
                    dynamic_cast<TH1D*>(h->Clone(Form("%s_%s_residual",
                        v->name.c_str(),CR_KEY[ir]))));
                if (v->coarse[ir].residual) v->coarse[ir].residual->SetDirectory(nullptr);
            } // endif
        } // endfor
    } // endif

    auto* rd=dynamic_cast<TDirectory*>(f.Get("residuals"));
    if (!rd) return nullptr;
    for (ib=0;ib<VAL_NBIN;ib++) {
        auto* hf=dynamic_cast<TH1D*>(rd->Get(Form("fit_b%03d",ib)));
        auto* hc=dynamic_cast<TH1D*>(rd->Get(Form("count_b%03d",ib)));
        if (hf) {
            v->bins[ib].residual_fit.reset(dynamic_cast<TH1D*>(hf->Clone(Form("%s_fit_b%03d",v->name.c_str(),ib))));
            if (v->bins[ib].residual_fit) v->bins[ib].residual_fit->SetDirectory(nullptr);
        }
        if (hc) {
            v->bins[ib].residual_count.reset(dynamic_cast<TH1D*>(hc->Clone(Form("%s_count_b%03d",v->name.c_str(),ib))));
            if (v->bins[ib].residual_count) v->bins[ib].residual_count->SetDirectory(nullptr);
        }
    } // endfor
    int ncomplete=0;
    get_param<int>(&f,"normalization_branches_complete",ncomplete);
    v->normalization_branches_complete=(ncomplete!=0);

    auto load_norm_dir = [&](const char* dirname,
                             std::vector<std::unique_ptr<TH1D>>& pre,
                             std::vector<std::unique_ptr<TH1D>>& after_ep,
                             std::vector<std::unique_ptr<TH1D>>& after_e_g,
                             std::vector<std::unique_ptr<TH1D>>& after_dphi,
                             std::vector<std::unique_ptr<TH1D>>& nminus1,
                             std::vector<std::unique_ptr<TH1D>>& full,
                             std::vector<std::unique_ptr<TH1D>>& low,
                             std::vector<std::unique_ptr<TH1D>>& high) {
        auto* nd=dynamic_cast<TDirectory*>(f.Get(dirname));
        for (int io=0;io<NORM_NOBS;io++) {
            auto clone_one=[&](const char* region)->std::unique_ptr<TH1D> {
                if (!nd) return nullptr;
                auto* h=dynamic_cast<TH1D*>(nd->Get(Form("%s_%s",region,NORM_OBS[io].key)));
                if (!h) return nullptr;
                std::unique_ptr<TH1D> q(dynamic_cast<TH1D*>(h->Clone(
                    Form("%s_%s_%s_%s",v->name.c_str(),dirname,region,NORM_OBS[io].key))));
                if (q) q->SetDirectory(nullptr);
                return q;
            };

            pre.push_back(clone_one("pre"));
            after_ep.push_back(clone_one("after_mx2ep"));
            after_e_g.push_back(clone_one("after_mx2ep_mx2eg"));
            after_dphi.push_back(clone_one("after_mx2ep_mx2eg_dphi"));
            nminus1.push_back(clone_one("nminus1"));
            full.push_back(clone_one("full"));
            low.push_back(clone_one("lowE"));
            high.push_back(clone_one("highE"));
        } // endfor
    };

    load_norm_dir("normalization_FD",
                  v->norm_pre,v->norm_after_mx2ep,v->norm_after_mx2ep_mx2eg,
                  v->norm_after_mx2ep_mx2eg_dphi,v->norm_nminus1,v->norm_full,
                  v->norm_lowE,v->norm_highE);

    load_norm_dir("normalization_FT",
                  v->norm_ft_pre,v->norm_ft_after_mx2ep,v->norm_ft_after_mx2ep_mx2eg,
                  v->norm_ft_after_mx2ep_mx2eg_dphi,v->norm_ft_nminus1,v->norm_ft_full,
                  v->norm_ft_lowE,v->norm_ft_highE);

    for (int ic=0;ic<6;ic++) {
        Long64_t q=0;
        get_param<Long64_t>(&f,Form("norm_fd_cutflow_%d",ic),q);
        v->norm_cutflow[ic]=q;
        q=0;
        get_param<Long64_t>(&f,Form("norm_ft_cutflow_%d",ic),q);
        v->norm_ft_cutflow[ic]=q;
    } // endfor

    auto load_corr_dir = [&](const char* dirname,
                             std::vector<std::unique_ptr<TH2D>>& corr,
                             std::vector<std::unique_ptr<TH1D>>& stages) {
        auto* d=dynamic_cast<TDirectory*>(f.Get(dirname));
        for (int jc=0;jc<CORR_NOBS;jc++) {
            std::unique_ptr<TH2D> q;
            if (d) {
                auto* h=dynamic_cast<TH2D*>(d->Get(Form("corr_%s",NORM_OBS[CORR_OBS[jc]].key)));
                if (h) {
                    q.reset(dynamic_cast<TH2D*>(h->Clone(
                        Form("%s_%s_corr_%s",v->name.c_str(),dirname,NORM_OBS[CORR_OBS[jc]].key))));
                    if (q) q->SetDirectory(nullptr);
                } // endif
            } // endif
            corr.push_back(std::move(q));
        } // endfor
        for (int is=0;is<EPROBE_NSTAGE;is++) {
            std::unique_ptr<TH1D> q;
            if (d) {
                auto* h=dynamic_cast<TH1D*>(d->Get(Form("eprobe_%s",EPROBE_STAGE_KEY[is])));
                if (h) {
                    q.reset(dynamic_cast<TH1D*>(h->Clone(
                        Form("%s_%s_eprobe_%s",v->name.c_str(),dirname,EPROBE_STAGE_KEY[is]))));
                    if (q) q->SetDirectory(nullptr);
                } // endif
            } // endif
            stages.push_back(std::move(q));
        } // endfor
    };

    load_corr_dir("correlations_FD",v->corr_fd,v->eprobe_stage_fd);
    load_corr_dir("correlations_FT",v->corr_ft,v->eprobe_stage_ft);

    {
        auto* sd=dynamic_cast<TDirectory*>(f.Get("cutscan_FT"));
        for (int iv=0;iv<CUTSCAN_NVAR;iv++) {
            std::unique_ptr<TH2D> q;
            if (sd) {
                auto* h=dynamic_cast<TH2D*>(sd->Get(CUTSCAN_KEY[iv]));
                if (h) {
                    q.reset(dynamic_cast<TH2D*>(h->Clone(
                        Form("%s_cutscan_FT_%s",v->name.c_str(),CUTSCAN_KEY[iv]))));
                    if (q) q->SetDirectory(nullptr);
                } // endif
            } // endif
            v->cutscan_ft.push_back(std::move(q));
        } // endfor
    }

    {
        auto* qd=dynamic_cast<TDirectory*>(f.Get("fitquality_FT"));
        v->fitqual_ft_low.resize(FT_COPL_SCAN_N);
        v->fitqual_ft_high.resize(FT_COPL_SCAN_N);
        for (int ic=0;ic<FT_COPL_SCAN_N;ic++) {
            TDirectory* cd=qd?dynamic_cast<TDirectory*>(qd->Get(Form("copl_%02d",ic))):nullptr;
            TDirectory* ld=cd?dynamic_cast<TDirectory*>(cd->Get("low")):nullptr;
            TDirectory* hd=cd?dynamic_cast<TDirectory*>(cd->Get("high")):nullptr;
            for (int io=0;io<NORM_NOBS;io++) {
                std::unique_ptr<TH2D> lo,hi;
                if (ld) {
                    auto* h=dynamic_cast<TH2D*>(ld->Get(NORM_OBS[io].key));
                    if (h) {
                        lo.reset(dynamic_cast<TH2D*>(h->Clone(
                            Form("%s_fitqual_low_%02d_%s",v->name.c_str(),ic,NORM_OBS[io].key))));
                        if (lo) lo->SetDirectory(nullptr);
                    } // endif
                } // endif
                if (hd) {
                    auto* h=dynamic_cast<TH2D*>(hd->Get(NORM_OBS[io].key));
                    if (h) {
                        hi.reset(dynamic_cast<TH2D*>(h->Clone(
                            Form("%s_fitqual_high_%02d_%s",v->name.c_str(),ic,NORM_OBS[io].key))));
                        if (hi) hi->SetDirectory(nullptr);
                    } // endif
                } // endif
                v->fitqual_ft_low[ic].push_back(std::move(lo));
                v->fitqual_ft_high[ic].push_back(std::move(hi));
            } // endfor
        } // endfor
    }

    f.Close();
    return v;
}
std::vector<std::unique_ptr<ValComponent>> build_val_components_parallel(const std::string&) {
    const SampleSpec all_specs[] = {
        {"data",DATA_DIR,false},
        {"aaogen",AAOGEN_DIR,true},
        {"clasdis",CLASDIS_DIR,true},
        {"dvcsgen",DVCSGEN_DIR,true}
    };

    struct Child { pid_t pid; SampleSpec spec; std::string path; };
    std::vector<Child> children;
    std::vector<SampleSpec> specs;
    std::map<std::string,std::string> paths;

    for (const auto& spec:all_specs) {
        if (!has_root_files(spec.dir)) continue;
        specs.push_back(spec);
        const std::string path=concise_cache_path(spec);
        paths[spec.name]=path;

        if (concise_valid_root(path)) {
            std::cout << "[CACHE] " << spec.name << " -> " << path << "\n";
            continue;
        } // endif

        const pid_t pid=fork();
        if (pid==0) {
            const bool ok=analyze_val_component_worker(spec,path);
            std::cout.flush();
            std::cerr.flush();
            _exit(ok?0:2);
        } else if (pid>0) {
            children.push_back({pid,spec,path});
            std::cout << "[SCAN] launched " << spec.name << " pid=" << pid << "\n";
        } else {
            std::cerr << "WARNING: fork failed for " << spec.name << "; running sequentially\n";
            analyze_val_component_worker(spec,path);
        } // endif
    } // endfor

    for (const auto& ch:children) {
        int status=0;
        waitpid(ch.pid,&status,0);
        if (!WIFEXITED(status) || WEXITSTATUS(status)!=0)
            std::cerr << "WARNING: worker failed for " << ch.spec.name << "\n";
        else
            std::cout << "[SCAN] finished " << ch.spec.name << "\n";
    } // endfor

    std::vector<std::unique_ptr<ValComponent>> out;
    for (const auto& spec:specs) {
        const std::string path=paths[spec.name];
        if (!concise_valid_root(path)) continue;
        auto v=load_val_component(path);
        if (v) out.push_back(std::move(v));
    } // endfor
    return out;
}


const ValComponent* find_val_component(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                       const std::string& name) {
    for (const auto& v:vv) if (v && v->name==name) return v.get();
    return nullptr;
}

void val_integrate_window(const TH1D* h,double lo,double hi,double& sw,double& sw2) {
    sw=0; sw2=0;
    if (!h || !(hi>lo)) return;
    const int nb=h->GetNbinsX();
    for (int j=1;j<=nb;j++) {
        const double x=h->GetBinCenter(j);
        if (x<=lo || x>=hi) continue;
        sw += h->GetBinContent(j);
        const double e=h->GetBinError(j);
        sw2 += e*e;
    } // endfor
}

ValEval val_eval_data(const ValComponentBin& b,const FitResult& fit,int ns) {
    ValEval e;
    if (!fit.valid || b.denom_rows<=0 || !b.residual_count) return e;
    e.denom=static_cast<double>(b.denom_rows);
    e.denom_w2=e.denom;
    val_integrate_window(b.residual_count.get(),fit.mean-ns*fit.sigma,fit.mean+ns*fit.sigma,
                         e.num,e.num_w2);
    e.efficiency=e.num/e.denom;
    const double fail_w2=std::max(0.0,e.denom_w2-e.num_w2);
    const double var=((1.0-e.efficiency)*(1.0-e.efficiency)*e.num_w2 +
                      e.efficiency*e.efficiency*fail_w2)/(e.denom*e.denom);
    e.efficiency_err=std::sqrt(std::max(0.0,var));
    e.valid=std::isfinite(e.efficiency) && e.efficiency>=0;
    return e;
}

ValEval val_eval_mc(const std::vector<std::unique_ptr<ValComponent>>& vv,int ib,
                    const ValNormSet& norm,const FitResult& fit,int ns) {
    ValEval e;
    if (!fit.valid) return e;
    for (const auto& vp:vv) {
        if (!vp || !vp->is_mc) continue;
        const double scale=val_component_scale(vp->name,norm);
        const auto& b=vp->bins[ib];
        e.denom += scale*b.denom_w;
        e.denom_w2 += scale*scale*b.denom_w2;
        double n=0,n2=0;
        val_integrate_window(b.residual_count.get(),fit.mean-ns*fit.sigma,fit.mean+ns*fit.sigma,n,n2);
        e.num += scale*n;
        e.num_w2 += scale*scale*n2;
    } // endfor
    if (!(e.denom>0)) return e;
    e.efficiency=e.num/e.denom;
    const double fail_w2=std::max(0.0,e.denom_w2-e.num_w2);
    const double var=((1.0-e.efficiency)*(1.0-e.efficiency)*e.num_w2 +
                      e.efficiency*e.efficiency*fail_w2)/(e.denom*e.denom);
    e.efficiency_err=std::sqrt(std::max(0.0,var));
    e.valid=std::isfinite(e.efficiency) && e.efficiency>=0;
    return e;
}

std::unique_ptr<TH1D> val_make_data_hist(const ValComponent* data,int ib,const char* name) {
    if (!data || !data->bins[ib].residual_fit) return nullptr;
    std::unique_ptr<TH1D> h(dynamic_cast<TH1D*>(data->bins[ib].residual_fit->Clone(name)));
    if (h) h->SetDirectory(nullptr);
    return h;
}

std::unique_ptr<TH1D> val_make_mc_hist(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                       int ib,const ValNormSet& norm,const char* name) {
    std::unique_ptr<TH1D> h(new TH1D(name,";#Delta p_{#gamma2} (GeV);Weighted candidates",
                                     VAL_DP_NBIN,VAL_DP_MIN,VAL_DP_MAX));
    h->Sumw2(); h->SetDirectory(nullptr);
    for (const auto& vp:vv) {
        if (!vp || !vp->is_mc || !vp->bins[ib].residual_fit) continue;
        const double scale=val_component_scale(vp->name,norm);
        h->Add(vp->bins[ib].residual_fit.get(),scale);
    } // endfor
    return h;
}
std::array<ValBinResult,VAL_NBIN> evaluate_valerii_fd(
        const std::vector<std::unique_ptr<ValComponent>>& vv,
        const ValNormSet& norm,
        std::vector<std::unique_ptr<TH1D>>* data_hists=nullptr,
        std::vector<std::unique_ptr<TH1D>>* mc_hists=nullptr) {
    std::array<ValBinResult,VAL_NBIN> rr;
    const ValComponent* data=find_val_component(vv,"data");
    if (!data) return rr;

    for (int ib=0;ib<VAL_NBIN;ib++) {
        auto hd=val_make_data_hist(data,ib,Form("h_val_data_b%03d",ib));
        auto hm=val_make_mc_hist(vv,ib,norm,Form("h_val_mc_%s_b%03d",norm.label.c_str(),ib));
        rr[ib].data_fit=fit_valerii_residual(hd.get());
        rr[ib].mc_fit=fit_valerii_residual(hm.get());

        for (int ns=1;ns<=3;ns++) {
            rr[ib].data[ns-1]=val_eval_data(data->bins[ib],rr[ib].data_fit,ns);
            rr[ib].mc[ns-1]=val_eval_mc(vv,ib,norm,rr[ib].mc_fit,ns);
            const auto& ed=rr[ib].data[ns-1];
            const auto& em=rr[ib].mc[ns-1];
            if (ed.valid && em.valid && em.efficiency>0 && ed.efficiency>=0) {
                const double c=ed.efficiency/em.efficiency;
                rr[ib].correction[ns-1]=c;
                double rel2=0;
                if (ed.efficiency>0) rel2+=std::pow(ed.efficiency_err/ed.efficiency,2);
                if (em.efficiency>0) rel2+=std::pow(em.efficiency_err/em.efficiency,2);
                rr[ib].correction_stat_err[ns-1]=std::fabs(c)*std::sqrt(rel2);
            }
        } // endfor
        rr[ib].nominal_valid=rr[ib].data_fit.valid && rr[ib].mc_fit.valid &&
                             rr[ib].correction[1]>0 && std::isfinite(rr[ib].correction[1]);

        if (data_hists) data_hists->push_back(std::move(hd));
        if (mc_hists) mc_hists->push_back(std::move(hm));
    } // endfor
    return rr;
}

void val_draw_map(const std::array<ValBinResult,VAL_NBIN>& rr,
                  const std::string& what,const std::string& outfile) {
    // Match Valerii's note display convention: show the first six momentum bins;
    // retain the 3.7-6 GeV bin in CSV/ROOT products.
    //
    // IMPORTANT ROOT ownership detail: objects drawn on a TPad are referenced
    // by pointer.  A stack-local TH2D destroyed at the end of each loop
    // iteration leaves the pad holding a dangling pointer; the later SaveAs()
    // then produces a completely blank canvas.  Keep all six maps alive until
    // after the canvas has been painted and written.
    TCanvas c(Form("c_%s",what.c_str()),"",1500,920);
    c.Divide(3,2,0.003,0.003);
    std::vector<std::unique_ptr<TH2D>> maps;
    maps.reserve(6);
    for (int ip=0;ip<6;ip++) {
        maps.emplace_back(new TH2D(Form("hm_%s_%d",what.c_str(),ip),
               Form("%.2f < p < %.2f GeV;wrapped #phi (deg);#theta (deg)",
                    VAL_P_EDGES[ip],VAL_P_EDGES[ip+1]),
               VAL_NPH,VAL_PH_EDGES,VAL_NT,VAL_T_EDGES));
        TH2D* h=maps.back().get();
        h->SetDirectory(nullptr);
        for (int it=0;it<VAL_NT;it++) for (int iph=0;iph<VAL_NPH;iph++) {
            const int ib=(ip*VAL_NT+it)*VAL_NPH+iph;
            double z=0;
            if (what=="data_eff") z=rr[ib].data[1].valid?rr[ib].data[1].efficiency:0;
            else if (what=="mc_eff") z=rr[ib].mc[1].valid?rr[ib].mc[1].efficiency:0;
            else if (what=="correction") z=rr[ib].nominal_valid?rr[ib].correction[1]:0;
            else if (what=="data_sigma") z=rr[ib].data_fit.valid?rr[ib].data_fit.sigma:0;
            else if (what=="mc_sigma") z=rr[ib].mc_fit.valid?rr[ib].mc_fit.sigma:0;
            h->SetBinContent(iph+1,it+1,z);
        } // endfor
        h->SetStats(0);
        c.cd(ip+1);
        gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.13);
        h->Draw("COLZ TEXT");
        gPad->Modified();
        gPad->Update();
    } // endfor
    c.Modified();
    c.Update();
    c.SaveAs(outfile.c_str());
}



void draw_step1b_mx2_ep_overlay(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                const std::string& file) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;
    if (data->norm_pre.size()<=NORM_MX2_EP || aao->norm_pre.size()<=NORM_MX2_EP ||
        cls->norm_pre.size()<=NORM_MX2_EP || dvc->norm_pre.size()<=NORM_MX2_EP) return;

    const TH1D* hd0=data->norm_pre[NORM_MX2_EP].get();
    const TH1D* ha0=aao->norm_pre[NORM_MX2_EP].get();
    const TH1D* hc0=cls->norm_pre[NORM_MX2_EP].get();
    const TH1D* hv0=dvc->norm_pre[NORM_MX2_EP].get();
    if (!hd0 || !ha0 || !hc0 || !hv0) return;

    std::unique_ptr<TH1D> hd((TH1D*)hd0->Clone("step1b_data"));
    std::unique_ptr<TH1D> ha((TH1D*)ha0->Clone("step1b_aao"));
    std::unique_ptr<TH1D> hc((TH1D*)hc0->Clone("step1b_cls"));
    std::unique_ptr<TH1D> hv((TH1D*)hv0->Clone("step1b_dvc"));
    hd->SetDirectory(nullptr); ha->SetDirectory(nullptr);
    hc->SetDirectory(nullptr); hv->SetDirectory(nullptr);

    auto unit=[](TH1D* h) {
        const double q=h ? h->Integral() : 0.0;
        if (h && q>0.0) h->Scale(1.0/q);
    };
    unit(hd.get()); unit(ha.get()); unit(hc.get()); unit(hv.get());

    hd->SetStats(0); hd->SetMarkerStyle(20); hd->SetMarkerSize(0.65);
    hd->SetLineColor(kBlack); hd->SetMarkerColor(kBlack);
    style_norm_component(ha.get(),kRed+1);
    style_norm_component(hc.get(),kOrange+7);
    style_norm_component(hv.get(),kGreen+2);

    TCanvas can("c_step1b_mx2_overlay","",1100,780);
    hd->SetTitle("Step 1B: baseline M_{X}^{2}(ep) and low-mass exclusivity window");
    hd->GetXaxis()->SetTitle("M_{X}^{2}(ep) (GeV^{2})");
    hd->GetYaxis()->SetTitle("Unit-area candidates");
    hd->SetMaximum(1.30*std::max({hd->GetMaximum(),ha->GetMaximum(),hc->GetMaximum(),hv->GetMaximum()}));

    hd->Draw("E1");
    ha->Draw("HIST SAME");
    hc->Draw("HIST SAME");
    hv->Draw("HIST SAME");
    hd->Draw("E1 SAME");

    const double ymax=hd->GetMaximum()*1.24;
    TLine llo(NORM_MX2_EP_MIN,0.0,NORM_MX2_EP_MIN,ymax);
    TLine lhi(NORM_MX2_EP_MAX,0.0,NORM_MX2_EP_MAX,ymax);
    llo.SetLineStyle(2); lhi.SetLineStyle(2);
    llo.SetLineWidth(2); lhi.SetLineWidth(2);
    llo.Draw(); lhi.Draw();

    TLegend leg(0.60,0.64,0.89,0.89);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(hd.get(),"Data","lep");
    leg.AddEntry(ha.get(),"AAO (unit area)","l");
    leg.AddEntry(hc.get(),"CLASDIS (unit area)","l");
    leg.AddEntry(hv.get(),"DVCSgen (unit area)","l");
    leg.AddEntry(&llo,"M_{X}^{2}(ep) selection","l");
    leg.Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.036);
    tx.DrawLatex(0.14,0.86,Form("%.3f < M_{X}^{2}(ep) < %.3f GeV^{2}",
                               NORM_MX2_EP_MIN,NORM_MX2_EP_MAX));

    can.SaveAs(file.c_str());
}

void draw_step1b_mx2_ep_individual(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                   const std::string& file) {
    TCanvas can("c_step1b_mx2_individual","",1350,950);
    can.Divide(2,2);

    const char* names[4]={"data","aaogen","clasdis","dvcsgen"};
    const char* labels[4]={"Data","AAOgen","CLASDIS","DVCSgen"};
    const int colors[4]={kBlack,kRed+1,kOrange+7,kGreen+2};

    std::vector<std::unique_ptr<TH1D>> keep;
    for (int is=0;is<4;is++) {
        const ValComponent* v=find_val_component(vv,names[is]);
        if (!v || v->norm_pre.size()<=NORM_MX2_EP || !v->norm_pre[NORM_MX2_EP]) continue;

        std::unique_ptr<TH1D> h((TH1D*)v->norm_pre[NORM_MX2_EP]->Clone(Form("step1b_%s_counts",names[is])));
        h->SetDirectory(nullptr);
        h->SetStats(0);
        h->SetLineColor(colors[is]);
        h->SetMarkerColor(colors[is]);
        h->SetLineWidth(2);
        h->SetTitle(Form("%s baseline M_{X}^{2}(ep)",labels[is]));
        h->GetXaxis()->SetTitle("M_{X}^{2}(ep) (GeV^{2})");
        h->GetYaxis()->SetTitle("Candidates");

        can.cd(is+1);
        gPad->SetRightMargin(0.05);
        h->Draw("HIST");

        const double ymax=std::max(1.0,1.05*h->GetMaximum());
        TLine* llo=new TLine(NORM_MX2_EP_MIN,0.0,NORM_MX2_EP_MIN,ymax);
        TLine* lhi=new TLine(NORM_MX2_EP_MAX,0.0,NORM_MX2_EP_MAX,ymax);
        llo->SetLineStyle(2); lhi->SetLineStyle(2);
        llo->SetLineWidth(2); lhi->SetLineWidth(2);
        llo->Draw(); lhi->Draw();

        const long long n0=v->norm_cutflow[0];
        const long long n1=v->norm_cutflow[1];
        const double frac=(n0>0 ? double(n1)/double(n0) : 0.0);
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
        tx.DrawLatex(0.14,0.86,Form("pass = %lld / %lld = %.2f%%",n1,n0,100.0*frac));

        keep.push_back(std::move(h));
    }

    can.SaveAs(file.c_str());
}

void write_step1b_mx2_ep_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                 const std::string& dir) {
    gSystem->mkdir(dir.c_str(),true);

    std::ofstream csv(dir+"/mx2_ep_survival.csv");
    csv << "sample,baseline,pass_mx2_ep,fail_mx2_ep,survival_fraction\n";

    std::ofstream txt(dir+"/mx2_ep_summary.txt");
    txt << "Step 1B: M_X^2(ep) exclusivity requirement\n"
        << "==========================================\n"
        << "Definition: M_X^2(ep) = (k + p_target - k' - p')^2.\n"
        << "Active window: " << NORM_MX2_EP_MIN << " < M_X^2(ep) < "
        << NORM_MX2_EP_MAX << " GeV^2.\n"
        << "This is evaluated after the common reconstructed ep-gamma-X baseline and\n"
        << "before M_X^2(e-gamma), Trento coplanarity, and angle(gamma,X) cuts.\n\n";

    for (const auto& vp:vv) {
        if (!vp) continue;
        const long long n0=vp->norm_cutflow[0];
        const long long n1=vp->norm_cutflow[1];
        const long long nf=std::max(0LL,n0-n1);
        const double frac=(n0>0 ? double(n1)/double(n0) : 0.0);

        csv << vp->name << "," << n0 << "," << n1 << "," << nf << "," << frac << "\n";
        txt << vp->name << ": " << n1 << " / " << n0
            << " = " << 100.0*frac << "% survive\n";
    }

    csv.close();
    txt.close();

    draw_step1b_mx2_ep_overlay(vv,dir+"/mx2_ep_unit_area_with_cut.png");
    draw_step1b_mx2_ep_individual(vv,dir+"/mx2_ep_counts_with_cut.png");
}



void draw_step1c_mx2_eg_overlay(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                const std::string& file) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;

    auto geth=[](const ValComponent* v)->const TH1D* {
        if (!v || v->norm_after_mx2ep.size()<=NORM_MX2_EG) return nullptr;
        return v->norm_after_mx2ep[NORM_MX2_EG].get();
    };
    const TH1D *hd0=geth(data), *ha0=geth(aao), *hc0=geth(cls), *hv0=geth(dvc);
    if (!hd0 || !ha0 || !hc0 || !hv0) return;

    std::unique_ptr<TH1D> hd((TH1D*)hd0->Clone("step1c_data"));
    std::unique_ptr<TH1D> ha((TH1D*)ha0->Clone("step1c_aao"));
    std::unique_ptr<TH1D> hc((TH1D*)hc0->Clone("step1c_cls"));
    std::unique_ptr<TH1D> hv((TH1D*)hv0->Clone("step1c_dvc"));
    hd->SetDirectory(nullptr); ha->SetDirectory(nullptr);
    hc->SetDirectory(nullptr); hv->SetDirectory(nullptr);

    auto unit=[](TH1D* h) {
        const double q=h ? h->Integral() : 0.0;
        if (h && q>0.0) h->Scale(1.0/q);
    };
    unit(hd.get()); unit(ha.get()); unit(hc.get()); unit(hv.get());

    hd->SetStats(0); hd->SetMarkerStyle(20); hd->SetMarkerSize(0.65);
    hd->SetLineColor(kBlack); hd->SetMarkerColor(kBlack);
    style_norm_component(ha.get(),kRed+1);
    style_norm_component(hc.get(),kOrange+7);
    style_norm_component(hv.get(),kGreen+2);

    TCanvas can("c_step1c_mx2eg_overlay","",1100,780);
    hd->SetTitle("Step 1C: M_{X}^{2}(e#gamma) after the M_{X}^{2}(ep) requirement");
    hd->GetXaxis()->SetTitle("M_{X}^{2}(e#gamma) (GeV^{2})");
    hd->GetYaxis()->SetTitle("Unit-area candidates");
    hd->SetMaximum(1.28*std::max({hd->GetMaximum(),ha->GetMaximum(),hc->GetMaximum(),hv->GetMaximum()}));

    hd->Draw("E1");
    ha->Draw("HIST SAME");
    hc->Draw("HIST SAME");
    hv->Draw("HIST SAME");
    hd->Draw("E1 SAME");

    const double ymax=hd->GetMaximum()*1.22;
    TLine cut(NORM_MX2_EG_MIN,0.0,NORM_MX2_EG_MIN,ymax);
    cut.SetLineStyle(2); cut.SetLineWidth(2); cut.Draw();

    TLegend leg(0.59,0.64,0.89,0.89);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(hd.get(),"Data","lep");
    leg.AddEntry(ha.get(),"AAO (unit area)","l");
    leg.AddEntry(hc.get(),"CLASDIS (unit area)","l");
    leg.AddEntry(hv.get(),"DVCSgen (unit area)","l");
    leg.AddEntry(&cut,"M_{X}^{2}(e#gamma)>1.4 GeV^{2}","l");
    leg.Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
    tx.DrawLatex(0.14,0.86,"Input sample already passes the Step-1B low-mass exclusivity window");

    can.SaveAs(file.c_str());
}

void draw_step1c_mx2_eg_individual(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                   const std::string& file) {
    TCanvas can("c_step1c_mx2eg_individual","",1350,950);
    can.Divide(2,2);

    const char* names[4]={"data","aaogen","clasdis","dvcsgen"};
    const char* labels[4]={"Data","AAOgen","CLASDIS","DVCSgen"};
    const int colors[4]={kBlack,kRed+1,kOrange+7,kGreen+2};

    std::vector<std::unique_ptr<TH1D>> keep;
    for (int is=0;is<4;is++) {
        const ValComponent* v=find_val_component(vv,names[is]);
        if (!v || v->norm_after_mx2ep.size()<=NORM_MX2_EG ||
            !v->norm_after_mx2ep[NORM_MX2_EG]) continue;

        std::unique_ptr<TH1D> h((TH1D*)v->norm_after_mx2ep[NORM_MX2_EG]->Clone(Form("step1c_%s_counts",names[is])));
        h->SetDirectory(nullptr);
        h->SetStats(0); h->SetLineColor(colors[is]); h->SetMarkerColor(colors[is]);
        h->SetLineWidth(2);
        h->SetTitle(Form("%s: M_{X}^{2}(e#gamma) after M_{X}^{2}(ep)",labels[is]));
        h->GetXaxis()->SetTitle("M_{X}^{2}(e#gamma) (GeV^{2})");
        h->GetYaxis()->SetTitle("Candidates");

        can.cd(is+1);
        gPad->SetRightMargin(0.05);
        h->Draw("HIST");

        const double ymax=std::max(1.0,1.05*h->GetMaximum());
        TLine* cut=new TLine(NORM_MX2_EG_MIN,0.0,NORM_MX2_EG_MIN,ymax);
        cut->SetLineStyle(2); cut->SetLineWidth(2); cut->Draw();

        const long long nin=v->norm_cutflow[1];
        const long long nout=v->norm_cutflow[2];
        const double frac=(nin>0 ? double(nout)/double(nin) : 0.0);
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
        tx.DrawLatex(0.14,0.86,Form("pass = %lld / %lld = %.2f%%",nout,nin,100.0*frac));

        keep.push_back(std::move(h));
    }

    can.SaveAs(file.c_str());
}

void write_step1c_mx2_eg_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                 const std::string& dir) {
    gSystem->mkdir(dir.c_str(),true);

    std::ofstream csv(dir+"/mx2_eg_survival.csv");
    csv << "sample,input_after_mx2_ep,pass_mx2_eg,fail_mx2_eg,survival_fraction\n";

    std::ofstream txt(dir+"/mx2_eg_summary.txt");
    txt << "Step 1C: M_X^2(e gamma) exclusivity requirement\n"
        << "===============================================\n"
        << "Definition: M_X^2(e gamma) = (k + p_target - k' - gamma_tag)^2.\n"
        << "The missing system therefore contains the reconstructed proton plus any\n"
        << "additional undetected final-state particles, including the probe photon in\n"
        << "a true ep -> ep pi0 -> ep gamma gamma event.\n"
        << "Active requirement: M_X^2(e gamma) > " << NORM_MX2_EG_MIN << " GeV^2.\n"
        << "Input sample: events that already pass Step 1B's M_X^2(ep) window.\n"
        << "Later Trento-coplanarity and angle(gamma,X) cuts are NOT included in these\n"
        << "Step-1C diagnostic survival fractions.\n\n";

    for (const auto& vp:vv) {
        if (!vp) continue;
        const long long nin=vp->norm_cutflow[1];
        const long long nout=vp->norm_cutflow[2];
        const long long nf=std::max(0LL,nin-nout);
        const double frac=(nin>0 ? double(nout)/double(nin) : 0.0);

        csv << vp->name << "," << nin << "," << nout << "," << nf << "," << frac << "\n";
        txt << vp->name << ": " << nout << " / " << nin
            << " = " << 100.0*frac << "% survive\n";
    }

    csv.close();
    txt.close();

    draw_step1c_mx2_eg_overlay(vv,dir+"/mx2_eg_unit_area_after_mx2_ep.png");
    draw_step1c_mx2_eg_individual(vv,dir+"/mx2_eg_counts_after_mx2_ep.png");
}



void draw_step1d_coplanarity_overlay(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                     const std::string& file) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;

    auto geth=[](const ValComponent* v)->const TH1D* {
        if (!v || v->norm_after_mx2ep_mx2eg.size()<=NORM_DPHI_TRENTO_SHIFT180) return nullptr;
        return v->norm_after_mx2ep_mx2eg[NORM_DPHI_TRENTO_SHIFT180].get();
    };
    const TH1D *hd0=geth(data), *ha0=geth(aao), *hc0=geth(cls), *hv0=geth(dvc);
    if (!hd0 || !ha0 || !hc0 || !hv0) return;

    std::unique_ptr<TH1D> hd((TH1D*)hd0->Clone("step1d_data"));
    std::unique_ptr<TH1D> ha((TH1D*)ha0->Clone("step1d_aao"));
    std::unique_ptr<TH1D> hc((TH1D*)hc0->Clone("step1d_cls"));
    std::unique_ptr<TH1D> hv((TH1D*)hv0->Clone("step1d_dvc"));
    hd->SetDirectory(nullptr); ha->SetDirectory(nullptr);
    hc->SetDirectory(nullptr); hv->SetDirectory(nullptr);

    auto unit=[](TH1D* h) {
        const double q=h ? h->Integral() : 0.0;
        if (h && q>0.0) h->Scale(1.0/q);
    };
    unit(hd.get()); unit(ha.get()); unit(hc.get()); unit(hv.get());

    hd->SetStats(0); hd->SetMarkerStyle(20); hd->SetMarkerSize(0.65);
    hd->SetLineColor(kBlack); hd->SetMarkerColor(kBlack);
    style_norm_component(ha.get(),kRed+1);
    style_norm_component(hc.get(),kOrange+7);
    style_norm_component(hv.get(),kGreen+2);

    TCanvas can("c_step1d_copl_overlay","",1100,780);
    hd->SetTitle("Step 1D: Trento coplanarity after M_{X}^{2}(ep) and M_{X}^{2}(e#gamma)");
    hd->GetXaxis()->SetTitle("#Delta#phi_{copl} (deg)");
    hd->GetYaxis()->SetTitle("Unit-area candidates");
    hd->SetMaximum(1.30*std::max({hd->GetMaximum(),ha->GetMaximum(),hc->GetMaximum(),hv->GetMaximum()}));
    hd->Draw("E1"); ha->Draw("HIST SAME"); hc->Draw("HIST SAME"); hv->Draw("HIST SAME"); hd->Draw("E1 SAME");

    const double ymax=hd->GetMaximum()*1.24;
    TLine llo(-NORM_DPHI_TRENTO_MAX,0.0,-NORM_DPHI_TRENTO_MAX,ymax);
    TLine lhi(+NORM_DPHI_TRENTO_MAX,0.0,+NORM_DPHI_TRENTO_MAX,ymax);
    llo.SetLineStyle(2); lhi.SetLineStyle(2); llo.SetLineWidth(2); lhi.SetLineWidth(2);
    llo.Draw(); lhi.Draw();

    TLegend leg(0.60,0.64,0.89,0.89);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(hd.get(),"Data","lep");
    leg.AddEntry(ha.get(),"AAO (unit area)","l");
    leg.AddEntry(hc.get(),"CLASDIS (unit area)","l");
    leg.AddEntry(hv.get(),"DVCSgen (unit area)","l");
    leg.AddEntry(&llo,Form("|#Delta#phi_{copl}| < %.1f deg",NORM_DPHI_TRENTO_MAX),"l");
    leg.Draw();

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
    tx.DrawLatex(0.14,0.86,"Zero corresponds to back-to-back proton and tag-photon transverse directions");
    can.SaveAs(file.c_str());
}

void draw_step1d_coplanarity_individual(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                        const std::string& file) {
    TCanvas can("c_step1d_copl_individual","",1350,950);
    can.Divide(2,2);
    const char* names[4]={"data","aaogen","clasdis","dvcsgen"};
    const char* labels[4]={"Data","AAOgen","CLASDIS","DVCSgen"};
    const int colors[4]={kBlack,kRed+1,kOrange+7,kGreen+2};
    std::vector<std::unique_ptr<TH1D>> keep;

    for (int is=0;is<4;is++) {
        const ValComponent* v=find_val_component(vv,names[is]);
        if (!v || v->norm_after_mx2ep_mx2eg.size()<=NORM_DPHI_TRENTO_SHIFT180 ||
            !v->norm_after_mx2ep_mx2eg[NORM_DPHI_TRENTO_SHIFT180]) continue;
        std::unique_ptr<TH1D> h((TH1D*)v->norm_after_mx2ep_mx2eg[NORM_DPHI_TRENTO_SHIFT180]->Clone(Form("step1d_%s_counts",names[is])));
        h->SetDirectory(nullptr); h->SetStats(0); h->SetLineColor(colors[is]); h->SetLineWidth(2);
        h->SetTitle(Form("%s: Trento coplanarity after Steps 1B+1C",labels[is]));
        h->GetXaxis()->SetTitle("#Delta#phi_{copl} (deg)");
        h->GetYaxis()->SetTitle("Candidates");
        can.cd(is+1); h->Draw("HIST");
        const double ymax=std::max(1.0,1.05*h->GetMaximum());
        TLine* llo=new TLine(-NORM_DPHI_TRENTO_MAX,0.0,-NORM_DPHI_TRENTO_MAX,ymax);
        TLine* lhi=new TLine(+NORM_DPHI_TRENTO_MAX,0.0,+NORM_DPHI_TRENTO_MAX,ymax);
        llo->SetLineStyle(2); lhi->SetLineStyle(2); llo->SetLineWidth(2); lhi->SetLineWidth(2);
        llo->Draw(); lhi->Draw();

        const long long nin=v->norm_cutflow[2];
        const long long nout=v->norm_cutflow[3];
        const double frac=(nin>0 ? double(nout)/double(nin) : 0.0);
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
        tx.DrawLatex(0.14,0.86,Form("pass = %lld / %lld = %.2f%%",nout,nin,100.0*frac));
        keep.push_back(std::move(h));
    }
    can.SaveAs(file.c_str());
}

void write_step1d_coplanarity_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                      const std::string& dir) {
    gSystem->mkdir(dir.c_str(),true);
    std::ofstream csv(dir+"/coplanarity_survival.csv");
    csv << "sample,input_after_mx2_ep_and_mx2_eg,pass_coplanarity,fail_coplanarity,survival_fraction\n";
    std::ofstream txt(dir+"/coplanarity_summary.txt");
    txt << "Step 1D: Trento coplanarity requirement\n"
        << "=======================================\n"
        << "Input sample: events already passing M_X^2(ep) and M_X^2(e gamma).\n"
        << "Residual definition: wrapped Trento-style coplanarity residual with zero\n"
        << "at the back-to-back proton/tag-photon transverse configuration.\n"
        << "Active requirement: |Delta phi_copl| < " << NORM_DPHI_TRENTO_MAX << " deg.\n"
        << "The later angle(gamma,X) requirement is NOT included here.\n\n";
    for (const auto& vp:vv) {
        if (!vp) continue;
        const long long nin=vp->norm_cutflow[2], nout=vp->norm_cutflow[3];
        const long long nf=std::max(0LL,nin-nout);
        const double frac=(nin>0 ? double(nout)/double(nin) : 0.0);
        csv << vp->name << "," << nin << "," << nout << "," << nf << "," << frac << "\n";
        txt << vp->name << ": " << nout << " / " << nin << " = " << 100.0*frac << "% survive\n";
    }
    csv.close(); txt.close();
    draw_step1d_coplanarity_overlay(vv,dir+"/coplanarity_unit_area_after_mx2_cuts.png");
    draw_step1d_coplanarity_individual(vv,dir+"/coplanarity_counts_after_mx2_cuts.png");
}



void draw_step1e_angle_gX_overlay(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                  const std::string& file) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;

    auto geth=[](const ValComponent* v)->const TH1D* {
        if (!v || v->norm_after_mx2ep_mx2eg_dphi.size()<=NORM_ANGLE_GX) return nullptr;
        return v->norm_after_mx2ep_mx2eg_dphi[NORM_ANGLE_GX].get();
    };
    const TH1D *hd0=geth(data), *ha0=geth(aao), *hc0=geth(cls), *hv0=geth(dvc);
    if (!hd0 || !ha0 || !hc0 || !hv0) return;

    std::unique_ptr<TH1D> hd((TH1D*)hd0->Clone("step1e_data"));
    std::unique_ptr<TH1D> ha((TH1D*)ha0->Clone("step1e_aao"));
    std::unique_ptr<TH1D> hc((TH1D*)hc0->Clone("step1e_cls"));
    std::unique_ptr<TH1D> hv((TH1D*)hv0->Clone("step1e_dvc"));
    hd->SetDirectory(nullptr); ha->SetDirectory(nullptr);
    hc->SetDirectory(nullptr); hv->SetDirectory(nullptr);

    auto unit=[](TH1D* h) {
        const double q=h ? h->Integral() : 0.0;
        if (h && q>0.0) h->Scale(1.0/q);
    };
    unit(hd.get()); unit(ha.get()); unit(hc.get()); unit(hv.get());

    hd->SetStats(0);
    hd->SetMarkerStyle(20);
    hd->SetMarkerSize(0.65);
    hd->SetLineColor(kBlack);
    hd->SetMarkerColor(kBlack);
    style_norm_component(ha.get(),kRed+1);
    style_norm_component(hc.get(),kOrange+7);
    style_norm_component(hv.get(),kGreen+2);

    TCanvas can("c_step1e_anglegX_overlay","",1100,780);
    hd->SetTitle("Step 1E: angle(#gamma,X) after Steps 1B+1C+1D");
    hd->GetXaxis()->SetTitle("angle(#gamma,X) (deg)");
    hd->GetYaxis()->SetTitle("Unit-area candidates");
    hd->SetMaximum(1.30*std::max({hd->GetMaximum(),ha->GetMaximum(),hc->GetMaximum(),hv->GetMaximum()}));

    hd->Draw("E1");
    ha->Draw("HIST SAME");
    hc->Draw("HIST SAME");
    hv->Draw("HIST SAME");
    hd->Draw("E1 SAME");

    const double ymax=hd->GetMaximum()*1.24;
    TLine cut(NORM_ANGLE_GX_MAX,0.0,NORM_ANGLE_GX_MAX,ymax);
    cut.SetLineStyle(2);
    cut.SetLineWidth(2);
    cut.Draw();

    TLegend leg(0.60,0.64,0.89,0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hd.get(),"Data","lep");
    leg.AddEntry(ha.get(),"AAO (unit area)","l");
    leg.AddEntry(hc.get(),"CLASDIS (unit area)","l");
    leg.AddEntry(hv.get(),"DVCSgen (unit area)","l");
    leg.AddEntry(&cut,Form("angle(#gamma,X) < %.1f deg",NORM_ANGLE_GX_MAX),"l");
    leg.Draw();

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.035);
    tx.DrawLatex(0.14,0.86,"Input sample already passes M_{X}^{2}(ep), M_{X}^{2}(e#gamma), and Trento coplanarity");

    can.SaveAs(file.c_str());
}

void draw_step1e_angle_gX_individual(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                     const std::string& file) {
    TCanvas can("c_step1e_anglegX_individual","",1350,950);
    can.Divide(2,2);

    const char* names[4]={"data","aaogen","clasdis","dvcsgen"};
    const char* labels[4]={"Data","AAOgen","CLASDIS","DVCSgen"};
    const int colors[4]={kBlack,kRed+1,kOrange+7,kGreen+2};
    std::vector<std::unique_ptr<TH1D>> keep;

    for (int is=0;is<4;is++) {
        const ValComponent* v=find_val_component(vv,names[is]);
        if (!v || v->norm_after_mx2ep_mx2eg_dphi.size()<=NORM_ANGLE_GX ||
            !v->norm_after_mx2ep_mx2eg_dphi[NORM_ANGLE_GX]) continue;

        std::unique_ptr<TH1D> h(
            (TH1D*)v->norm_after_mx2ep_mx2eg_dphi[NORM_ANGLE_GX]->
                Clone(Form("step1e_%s_counts",names[is]))
        );
        h->SetDirectory(nullptr);
        h->SetStats(0);
        h->SetLineColor(colors[is]);
        h->SetMarkerColor(colors[is]);
        h->SetLineWidth(2);
        h->SetTitle(Form("%s: angle(#gamma,X) after Steps 1B+1C+1D",labels[is]));
        h->GetXaxis()->SetTitle("angle(#gamma,X) (deg)");
        h->GetYaxis()->SetTitle("Candidates");

        can.cd(is+1);
        h->Draw("HIST");

        const double ymax=std::max(1.0,1.05*h->GetMaximum());
        TLine* cut=new TLine(NORM_ANGLE_GX_MAX,0.0,NORM_ANGLE_GX_MAX,ymax);
        cut->SetLineStyle(2);
        cut->SetLineWidth(2);
        cut->Draw();

        const long long nin=v->norm_cutflow[3];
        const long long nout=v->norm_cutflow[4];
        const double frac=(nin>0 ? double(nout)/double(nin) : 0.0);

        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.042);
        tx.DrawLatex(0.14,0.86,Form("pass = %lld / %lld = %.2f%%",nout,nin,100.0*frac));

        keep.push_back(std::move(h));
    }

    can.SaveAs(file.c_str());
}

void write_step1e_angle_gX_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                   const std::string& dir) {
    gSystem->mkdir(dir.c_str(),true);

    std::ofstream csv(dir+"/angle_gX_survival.csv");
    csv << "sample,input_after_steps_1B_1C_1D,pass_angle_gX,fail_angle_gX,survival_fraction\n";

    std::ofstream txt(dir+"/angle_gX_summary.txt");
    txt << "Step 1E: angle(#gamma,X) requirement\n"
        << "===================================\n"
        << "Input sample: events already passing M_X^2(ep), M_X^2(e gamma), and\n"
        << "the Trento coplanarity requirement.\n"
        << "Definition used by the macro: opening angle between the reconstructed\n"
        << "tag-photon direction and the inferred missing-object direction X.\n"
        << "Active requirement: angle(#gamma,X) < " << NORM_ANGLE_GX_MAX << " deg.\n"
        << "This is the final sequential exclusivity requirement in Step 1.\n\n";

    for (const auto& vp:vv) {
        if (!vp) continue;

        const long long nin=vp->norm_cutflow[3];
        const long long nout=vp->norm_cutflow[4];
        const long long nf=std::max(0LL,nin-nout);
        const double frac=(nin>0 ? double(nout)/double(nin) : 0.0);

        csv << vp->name << "," << nin << "," << nout << "," << nf << "," << frac << "\n";
        txt << vp->name << ": " << nout << " / " << nin
            << " = " << 100.0*frac << "% survive\n";
    }

    csv.close();
    txt.close();

    draw_step1e_angle_gX_overlay(vv,dir+"/angle_gX_unit_area_after_steps_1B_1C_1D.png");
    draw_step1e_angle_gX_individual(vv,dir+"/angle_gX_counts_after_steps_1B_1C_1D.png");
}



void write_step2a_energy_region_diagnostics(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                            const std::string& dir) {
    gSystem->mkdir(dir.c_str(),true);

    const int io=NORM_EGAMMA;
    const char* names[4]={"data","aaogen","clasdis","dvcsgen"};
    const char* labels[4]={"Data","AAOgen","CLASDIS","DVCSgen"};
    const int colors[4]={kBlack,kRed+1,kOrange+7,kGreen+2};

    std::ofstream csv(dir+"/energy_region_counts.csv");
    csv << "sample,total_selected,Egamma_lt2,Egamma_2to3,Egamma_gt3,"
           "frac_lt2,frac_2to3,frac_gt3\n";

    std::vector<std::unique_ptr<TH1D>> hu;
    std::vector<std::unique_ptr<TH1D>> hc;
    hu.reserve(4); hc.reserve(4);

    for (int is=0;is<4;is++) {
        const ValComponent* v=find_val_component(vv,names[is]);
        if (!v || v->norm_full.size()<=io || !v->norm_full[io]) continue;

        std::unique_ptr<TH1D> counts((TH1D*)v->norm_full[io]->Clone(Form("step2a_counts_%s",names[is])));
        counts->SetDirectory(nullptr);
        counts->SetStats(0);

        const int b04=counts->GetXaxis()->FindBin(0.4+1e-9);
        const int b2m=counts->GetXaxis()->FindBin(2.0-1e-9);
        const int b2p=counts->GetXaxis()->FindBin(2.0+1e-9);
        const int b3m=counts->GetXaxis()->FindBin(3.0-1e-9);
        const int b3p=counts->GetXaxis()->FindBin(3.0+1e-9);
        const int b8 =counts->GetXaxis()->FindBin(8.0-1e-9);

        const double nlo=counts->Integral(b04,b2m);
        const double nmid=counts->Integral(b2p,b3m);
        const double nhi=counts->Integral(b3p,b8);
        const double nt=nlo+nmid+nhi;

        csv << names[is] << "," << nt << "," << nlo << "," << nmid << "," << nhi << ","
            << (nt>0?nlo/nt:0) << "," << (nt>0?nmid/nt:0) << "," << (nt>0?nhi/nt:0) << "\n";

        std::unique_ptr<TH1D> unit((TH1D*)counts->Clone(Form("step2a_unit_%s",names[is])));
        unit->SetDirectory(nullptr);
        const double q=unit->Integral();
        if (q>0) unit->Scale(1.0/q);

        if (is==0) {
            counts->SetMarkerStyle(20); counts->SetMarkerSize(0.65);
            counts->SetLineColor(kBlack); counts->SetMarkerColor(kBlack);
            unit->SetMarkerStyle(20); unit->SetMarkerSize(0.65);
            unit->SetLineColor(kBlack); unit->SetMarkerColor(kBlack);
        } else {
            style_norm_component(counts.get(),colors[is]);
            style_norm_component(unit.get(),colors[is]);
        }

        hc.push_back(std::move(counts));
        hu.push_back(std::move(unit));
    }
    csv.close();

    if (hu.size()==4) {
        TCanvas c("c_step2a_energy_unit","",1150,800);
        hu[0]->SetTitle("Step 2A: Valerii photon-energy normalization regions");
        hu[0]->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
        hu[0]->GetYaxis()->SetTitle("Unit-area selected candidates");
        hu[0]->SetMaximum(1.28*std::max({hu[0]->GetMaximum(),hu[1]->GetMaximum(),
                                         hu[2]->GetMaximum(),hu[3]->GetMaximum()}));
        hu[0]->Draw("E1");
        hu[1]->Draw("HIST SAME"); hu[2]->Draw("HIST SAME"); hu[3]->Draw("HIST SAME");
        hu[0]->Draw("E1 SAME");

        const double ymax=hu[0]->GetMaximum()*1.22;
        TLine l2(2.0,0,2.0,ymax), l3(3.0,0,3.0,ymax);
        l2.SetLineStyle(2); l3.SetLineStyle(2); l2.SetLineWidth(2); l3.SetLineWidth(2);
        l2.Draw(); l3.Draw();

        TLegend leg(0.60,0.62,0.89,0.89);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(hu[0].get(),"Data","lep");
        leg.AddEntry(hu[1].get(),"AAOgen","l");
        leg.AddEntry(hu[2].get(),"CLASDIS","l");
        leg.AddEntry(hu[3].get(),"DVCSgen","l");
        leg.AddEntry(&l2,"E_{#gamma}=2 GeV","l");
        leg.AddEntry(&l3,"E_{#gamma}=3 GeV","l");
        leg.Draw();

        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.033);
        tx.DrawLatex(0.14,0.86,"<2 GeV: fit AAO+CLASDIS only; 2-3 GeV: transition/closure only");
        tx.DrawLatex(0.14,0.81,">3 GeV: hold AAO+CLASDIS fixed and fit DVCS");
        c.SaveAs((dir+"/photon_energy_regions_unit_area.png").c_str());
    }

    if (hc.size()==4) {
        TCanvas c("c_step2a_energy_counts","",1350,950);
        c.Divide(2,2);
        for (int is=0;is<4;is++) {
            c.cd(is+1);
            hc[is]->SetTitle(Form("%s selected tag-photon energy",labels[is]));
            hc[is]->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
            hc[is]->GetYaxis()->SetTitle("Candidates");
            hc[is]->Draw(is==0?"E1":"HIST");
            const double ymax=std::max(1.0,1.05*hc[is]->GetMaximum());
            TLine* l2=new TLine(2.0,0,2.0,ymax);
            TLine* l3=new TLine(3.0,0,3.0,ymax);
            l2->SetLineStyle(2); l3->SetLineStyle(2);
            l2->SetLineWidth(2); l3->SetLineWidth(2);
            l2->Draw(); l3->Draw();
        }
        c.SaveAs((dir+"/photon_energy_regions_counts.png").c_str());
    }

    std::ofstream txt(dir+"/energy_region_definition.txt");
    txt << "Valerii sequential normalization regions\n"
        << "=======================================\n"
        << "The input is the fully selected ep-gamma-X sample after the Step-1 exclusivity cuts.\n\n"
        << "Low-energy region: E_gamma < 2 GeV\n"
        << "  * determine AAO and CLASDIS normalization factors;\n"
        << "  * DVCS is deliberately excluded from the fit, following Valerii's procedure.\n\n"
        << "Transition region: 2 <= E_gamma <= 3 GeV\n"
        << "  * do not derive any normalization factor here;\n"
        << "  * retain only as a closure/cross-check region.\n\n"
        << "High-energy region: E_gamma > 3 GeV\n"
        << "  * hold the low-E AAO and CLASDIS factors fixed;\n"
        << "  * determine the DVCS normalization factor;\n"
        << "  * include Delta t = t_p - t_gamma among the high-E diagnostics.\n";
    txt.close();
}



void write_lowE_pi0_balance_diagnostics(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                        const std::vector<NormFitPoint>& low_points,
                                        const std::string& dir) {
    gSystem->mkdir(dir.c_str(),true);

    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    if (!data || !aao || !cls) return;

    std::ofstream csv(dir+"/lowE_pi0_balance.csv");
    csv << "observable,"
           "aao_scale,clasdis_scale,"
           "data_entries,raw_aao_entries,raw_clasdis_entries,"
           "scaled_aao_yield,scaled_clasdis_yield,scaled_aao_plus_clasdis,"
           "model_over_data,aao_fraction_of_aao_plus_clasdis\n";

    struct Row {
        std::string obs;
        double model_over_data=0;
        double aao_fraction=0;
        double aao_scale=0;
        double cls_scale=0;
    };
    std::vector<Row> rows;

    auto obs_index=[](const std::string& key)->int {
        for (int io=0;io<NORM_NOBS;io++) {
            if (key==NORM_OBS[io].key) return io;
        }
        return -1;
    };

    for (const auto& q:low_points) {
        if (!q.valid) continue;
        const int io=obs_index(q.observable);
        if (io<0) continue;
        if (io>=(int)data->norm_lowE.size() ||
            io>=(int)aao->norm_lowE.size() ||
            io>=(int)cls->norm_lowE.size()) continue;
        const TH1D* hd=data->norm_lowE[io].get();
        const TH1D* ha=aao->norm_lowE[io].get();
        const TH1D* hc=cls->norm_lowE[io].get();
        if (!hd || !ha || !hc) continue;

        const double nd=hd->Integral(1,hd->GetNbinsX());
        const double na=ha->Integral(1,ha->GetNbinsX());
        const double nc=hc->Integral(1,hc->GetNbinsX());
        const double ya=q.aao*na;
        const double yc=q.clasdis*nc;
        const double ysum=ya+yc;
        const double ratio=(nd>0 ? ysum/nd : 0.0);
        const double fa=(ysum>0 ? ya/ysum : 0.0);

        csv << q.observable << ","
            << q.aao << "," << q.clasdis << ","
            << nd << "," << na << "," << nc << ","
            << ya << "," << yc << "," << ysum << ","
            << ratio << "," << fa << "\n";

        rows.push_back({q.observable,ratio,fa,q.aao,q.clasdis});
    }
    csv.close();

    std::ofstream txt(dir+"/lowE_pi0_balance_summary.txt");
    txt << "Low-E AAO/CLASDIS balance diagnostic\n"
        << "===================================\n"
        << "These are the per-observable Valerii-style E_gamma<2 GeV fits.\n"
        << "AAO and CLASDIS are kept as separate templates because their phase-space\n"
        << "coverage and shapes differ, and CLASDIS is not assumed to be a pure exclusive\n"
        << "pi0 sample.  For the photon-efficiency problem, however, both templates can\n"
        << "supply pi0-bearing events with a real partner photon.  Therefore the fitted\n"
        << "sum of their yields is monitored explicitly in addition to the individual\n"
        << "scale factors.\n\n";
    for (const auto& r:rows) {
        txt << r.obs
            << ": model/data=" << r.model_over_data
            << ", AAO fraction of (AAO+CLASDIS)=" << r.aao_fraction
            << ", scales=(" << r.aao_scale << "," << r.cls_scale << ")\n";
    }
    txt.close();

    if (rows.empty()) return;

    // Plot 1: does the total fitted AAO+CLASDIS yield remain stable even when
    // their individual coefficients trade against one another?
    {
        TCanvas c("c_lowE_pi0_total_closure","",1100,720);
        TH1D h("h_lowE_pi0_total_closure",
               "Low-E fitted AAO+CLASDIS total;Normalization observable;(AAO+CLASDIS fitted yield) / data yield",
               (int)rows.size(),0,(int)rows.size());
        h.SetStats(0);
        h.SetMarkerStyle(20);
        h.SetMarkerSize(1.1);
        for (int i=0;i<(int)rows.size();i++) {
            h.GetXaxis()->SetBinLabel(i+1,rows[i].obs.c_str());
            h.SetBinContent(i+1,rows[i].model_over_data);
        }
        h.SetMinimum(0.0);
        h.SetMaximum(std::max(1.5,1.20*h.GetMaximum()));
        h.Draw("P");
        TLine one(0,1.0,rows.size(),1.0);
        one.SetLineStyle(2);
        one.SetLineWidth(2);
        one.Draw();
        c.SetBottomMargin(0.22);
        h.GetXaxis()->LabelsOption("v");
        c.SaveAs((dir+"/lowE_pi0_total_closure.png").c_str());
    }

    // Plot 2: show the generator tradeoff directly.
    {
        TCanvas c("c_lowE_pi0_balance","",1100,720);
        TH1D h("h_lowE_pi0_balance",
               "Low-E AAO share of fitted AAO+CLASDIS yield;Normalization observable;AAO / (AAO+CLASDIS)",
               (int)rows.size(),0,(int)rows.size());
        h.SetStats(0);
        h.SetMarkerStyle(20);
        h.SetMarkerSize(1.1);
        for (int i=0;i<(int)rows.size();i++) {
            h.GetXaxis()->SetBinLabel(i+1,rows[i].obs.c_str());
            h.SetBinContent(i+1,rows[i].aao_fraction);
        }
        h.SetMinimum(0.0);
        h.SetMaximum(1.0);
        h.Draw("P");
        c.SetBottomMargin(0.22);
        h.GetXaxis()->LabelsOption("v");
        c.SaveAs((dir+"/lowE_aao_fraction_of_pi0_model.png").c_str());
    }
}


NormDerivation derive_normalization(const std::vector<std::unique_ptr<ValComponent>>& vv,const std::string& out) {
    NormDerivation R;
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    const std::string od=out+"/normalization"; gSystem->mkdir(od.c_str(),true);
    gSystem->mkdir((od+"/preselection_shapes").c_str(),true);
    gSystem->mkdir((od+"/nminus1_shapes").c_str(),true);
    gSystem->mkdir((od+"/lowE_fits").c_str(),true);
    gSystem->mkdir((od+"/highE_fits").c_str(),true);
    gSystem->mkdir((od+"/final_closure").c_str(),true);
    gSystem->mkdir((od+"/step1b_mx2_ep").c_str(),true);
    gSystem->mkdir((od+"/step1c_mx2_eg").c_str(),true);
    gSystem->mkdir((od+"/step1d_coplanarity").c_str(),true);
    gSystem->mkdir((od+"/step1e_angle_gX").c_str(),true);
    gSystem->mkdir((od+"/step2a_energy_regions").c_str(),true);
    gSystem->mkdir((od+"/step2b_lowE_pi0_balance").c_str(),true);

    if (!data || !aao || !cls || !dvc ||
        data->norm_lowE.size()!=NORM_NOBS || aao->norm_lowE.size()!=NORM_NOBS ||
        cls->norm_lowE.size()!=NORM_NOBS || dvc->norm_lowE.size()!=NORM_NOBS) {
        R.used_fallback=true; R.nominal=VAL_HISTORICAL_NOMINAL; R.low=R.nominal; R.high=R.nominal; return R;
    }

    // Shape-only QA before any exclusivity cut and with N-1 selections.  These
    // plots intentionally normalize every sample to unit area: they answer the
    // resolution/shape question without conflating it with MC sample statistics.
    for (int io=0;io<NORM_NOBS;io++) {
        draw_norm_shape_overlay(data->norm_pre[io].get(),aao->norm_pre[io].get(),cls->norm_pre[io].get(),dvc->norm_pre[io].get(),
            std::string("Preselection shape QA: ")+NORM_OBS[io].key,
            od+"/preselection_shapes/pre_"+NORM_OBS[io].key+".png");
        draw_norm_shape_overlay(data->norm_nminus1[io].get(),aao->norm_nminus1[io].get(),cls->norm_nminus1[io].get(),dvc->norm_nminus1[io].get(),
            std::string("N-1 exclusivity shape QA: ")+NORM_OBS[io].key,
            od+"/nminus1_shapes/nminus1_"+NORM_OBS[io].key+".png");
    }
    draw_norm_cutflow(vv,od+"/normalization_cutflow.png");

    // Step 1B diagnostics are deliberately written before any template
    // normalization is attempted.  They document exactly what the first
    // exclusivity requirement does to the baseline ep-gamma-X sample.
    write_step1b_mx2_ep_summary(vv,od+"/step1b_mx2_ep");

    // Step 1C diagnostics: inspect M_X^2(e gamma) sequentially after Step 1B
    // and before the later coplanarity / angle cuts.
    write_step1c_mx2_eg_summary(vv,od+"/step1c_mx2_eg");

    // Step 1D diagnostics: inspect the Trento coplanarity residual sequentially
    // after Steps 1B+1C and before angle(gamma,X).
    write_step1d_coplanarity_summary(vv,od+"/step1d_coplanarity");

    // Step 1E diagnostics: inspect angle(gamma,X) sequentially after all
    // preceding Step-1 exclusivity requirements and before its own cut.
    write_step1e_angle_gX_summary(vv,od+"/step1e_angle_gX");

    // Step 2A: reproduce Valerii's energy-region logic before fitting any
    // component normalization.  The 2-3 GeV transition region is explicitly
    // excluded from normalization derivation.
    write_step2a_energy_region_diagnostics(vv,od+"/step2a_energy_regions");

    // Stage A: low-E AAO + CLASDIS.  DVCS is intentionally absent from the fit,
    // exactly following Valerii's E_gamma<2 GeV normalization procedure.  For each observable scan one common
    // reconstructed-MC shift and extra Gaussian resolution.  The same morph is
    // applied to AAO and CLASDIS, then their non-negative scale factors are
    // solved analytically.  This mirrors the template-morph philosophy already
    // used in the DVCS exclusivity-selection suite without modifying event-level
    // MC kinematics or the later photon-efficiency numerator/denominator.
    for (int io=0;io<NORM_NOBS;io++) if (NORM_OBS[io].use_low) {
        auto raw=fit_two_templates_fixed_shapes(data->norm_lowE[io].get(),aao->norm_lowE[io].get(),cls->norm_lowE[io].get(),NORM_OBS[io].key);
        auto q=fit_two_templates_morphed(data->norm_lowE[io].get(),aao->norm_lowE[io].get(),cls->norm_lowE[io].get(),NORM_OBS[io].key);
        if (q.valid) R.low_points.push_back(q);
        draw_norm_fit(data->norm_lowE[io].get(),aao->norm_lowE[io].get(),cls->norm_lowE[io].get(),dvc->norm_lowE[io].get(),
                      raw.valid?raw.aao:0,raw.valid?raw.clasdis:0,0.0,
                      std::string("Low-E raw templates: ")+NORM_OBS[io].key+" (E_{#gamma}<2 GeV)",
                      od+"/lowE_fits/raw_"+NORM_OBS[io].key+".png",raw.chi2,raw.ndf);
        draw_norm_fit(data->norm_lowE[io].get(),aao->norm_lowE[io].get(),cls->norm_lowE[io].get(),dvc->norm_lowE[io].get(),
                      q.valid?q.aao:0,q.valid?q.clasdis:0,0.0,
                      std::string("Low-E morphed templates: ")+NORM_OBS[io].key+" (E_{#gamma}<2 GeV)",
                      od+"/lowE_fits/morphed_"+NORM_OBS[io].key+".png",q.chi2,q.ndf,q.morph_shift,q.morph_sigma);
    }
    // Before reducing the per-observable fits to Valerii's average
    // normalization factors, explicitly inspect whether AAO and CLASDIS are
    // trading against one another while their summed pi0-bearing yield remains
    // stable.  This is diagnostic only: the nominal factors below still follow
    // Valerii's per-observable-average prescription.
    write_lowE_pi0_balance_diagnostics(vv,R.low_points,od+"/step2b_lowE_pi0_balance");

    const double A=mean_valid(R.low_points,0), B=mean_valid(R.low_points,1);

    // Stage B: with A/B fixed to the low-E means, determine DVCS above 3 GeV.
    // The 2-3 GeV transition region is NOT used to derive any scale factor.
    // Delta-t is defined as t_p - t_gamma and is included here as one of
    // Valerii's high-E exclusivity/normalization diagnostics.
    // Again one common morph is scanned for all three reconstructed-MC pieces.
    for (int io=0;io<NORM_NOBS;io++) if (NORM_OBS[io].use_high) {
        auto raw=fit_dvcs_fixed_shapes(data->norm_highE[io].get(),aao->norm_highE[io].get(),cls->norm_highE[io].get(),dvc->norm_highE[io].get(),A,B,NORM_OBS[io].key);
        auto q=fit_dvcs_template_morphed(data->norm_highE[io].get(),aao->norm_highE[io].get(),cls->norm_highE[io].get(),dvc->norm_highE[io].get(),A,B,NORM_OBS[io].key);
        if (q.valid) R.high_points.push_back(q);
        draw_norm_fit(data->norm_highE[io].get(),aao->norm_highE[io].get(),cls->norm_highE[io].get(),dvc->norm_highE[io].get(),
                      A,B,raw.valid?raw.dvcs:0,
                      std::string("High-E raw templates: ")+NORM_OBS[io].key+" (E_{#gamma}>3 GeV)",
                      od+"/highE_fits/raw_"+NORM_OBS[io].key+".png",raw.chi2,raw.ndf);
        draw_norm_fit(data->norm_highE[io].get(),aao->norm_highE[io].get(),cls->norm_highE[io].get(),dvc->norm_highE[io].get(),
                      A,B,q.valid?q.dvcs:0,
                      std::string("High-E morphed templates: ")+NORM_OBS[io].key+" (E_{#gamma}>3 GeV)",
                      od+"/highE_fits/morphed_"+NORM_OBS[io].key+".png",q.chi2,q.ndf,q.morph_shift,q.morph_sigma);
    }
    const double C=mean_valid(R.high_points,2);
    R.valid=(A>0 && B>0 && C>=0 && R.low_points.size()>=2 && R.high_points.size()>=2);
    if (!R.valid) {
        R.used_fallback=true; R.nominal=VAL_HISTORICAL_NOMINAL; R.low=R.nominal; R.high=R.nominal;
    } else {
        R.nominal={"derived_nominal",A,B,C};
        auto ar=range_valid(R.low_points,0,A), br=range_valid(R.low_points,1,B), cr=range_valid(R.high_points,2,C);
        R.low={"derived_low",ar.first,br.first,cr.first};
        R.high={"derived_high",ar.second,br.second,cr.second};
    }

    auto find_low=[&](const std::string& k)->const NormFitPoint* { for (const auto& q:R.low_points) if (q.observable==k) return &q; return nullptr; };
    auto find_high=[&](const std::string& k)->const NormFitPoint* { for (const auto& q:R.high_points) if (q.observable==k) return &q; return nullptr; };

    // Full selected-sample closure.  Save raw and morphed versions side-by-side.
    // The morph is a normalization-fit nuisance only; the normalization factors,
    // not the morph, propagate to the photon-efficiency extraction.
    for (int io=0;io<NORM_NOBS;io++) {
        draw_norm_fit(data->norm_full[io].get(),aao->norm_full[io].get(),cls->norm_full[io].get(),dvc->norm_full[io].get(),
                      R.nominal.aao,R.nominal.clasdis,R.nominal.dvcs,
                      std::string("Final selected closure (raw MC): ")+NORM_OBS[io].key,
                      od+"/final_closure/raw_"+NORM_OBS[io].key+".png");
        const NormFitPoint* q=find_low(NORM_OBS[io].key);
        if (!q) q=find_high(NORM_OBS[io].key);
        const double sh=q?q->morph_shift:0.0, sg=q?q->morph_sigma:0.0;
        draw_norm_fit(data->norm_full[io].get(),aao->norm_full[io].get(),cls->norm_full[io].get(),dvc->norm_full[io].get(),
                      R.nominal.aao,R.nominal.clasdis,R.nominal.dvcs,
                      std::string("Final selected closure (morphed MC QA): ")+NORM_OBS[io].key,
                      od+"/final_closure/morphed_"+NORM_OBS[io].key+".png",-1,0,sh,sg);
    }

    // Dedicated energy-composition plot with the two normalization-region boundaries.
    {
        const int io=NORM_EGAMMA;
        std::unique_ptr<TH1D> d((TH1D*)data->norm_full[io]->Clone("energy_data")); d->SetDirectory(nullptr);
        std::unique_ptr<TH1D> a((TH1D*)aao->norm_full[io]->Clone("energy_aao")); a->Scale(R.nominal.aao); a->SetDirectory(nullptr);
        std::unique_ptr<TH1D> c((TH1D*)cls->norm_full[io]->Clone("energy_cls")); c->Scale(R.nominal.clasdis); c->SetDirectory(nullptr);
        std::unique_ptr<TH1D> v((TH1D*)dvc->norm_full[io]->Clone("energy_dvc")); v->Scale(R.nominal.dvcs); v->SetDirectory(nullptr);
        std::unique_ptr<TH1D> t((TH1D*)a->Clone("energy_total")); t->Add(c.get()); t->Add(v.get()); t->SetDirectory(nullptr);
        d->SetMarkerStyle(20); d->SetMarkerSize(0.65); d->SetLineColor(kBlack); d->SetStats(0);
        style_norm_component(a.get(),kRed+1); style_norm_component(c.get(),kOrange+7); style_norm_component(v.get(),kGreen+2); style_norm_component(t.get(),kBlue+1,3);
        TCanvas ce("c_energy_regions","",1100,760); d->SetTitle("Normalization regions and fitted MC composition;E_{#gamma} (GeV);Candidates");
        d->SetMaximum(1.25*std::max(d->GetMaximum(),t->GetMaximum())); d->Draw("E1"); a->Draw("HIST SAME"); c->Draw("HIST SAME"); v->Draw("HIST SAME"); t->Draw("HIST SAME"); d->Draw("E1 SAME");
        TLine l2(2.0,0,2.0,d->GetMaximum()); l2.SetLineStyle(2); l2.SetLineWidth(2); l2.Draw();
        TLine l3(3.0,0,3.0,d->GetMaximum()); l3.SetLineStyle(2); l3.SetLineWidth(2); l3.Draw();
        TLegend leg(0.62,0.66,0.89,0.89); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.AddEntry(d.get(),"Data","lep"); leg.AddEntry(a.get(),"AAO","l"); leg.AddEntry(c.get(),"CLASDIS","l"); leg.AddEntry(v.get(),"DVCSgen","l"); leg.AddEntry(t.get(),"Total MC","l"); leg.Draw();
        ce.SaveAs((od+"/normalization_energy_regions.png").c_str());
    }

    // Cut-flow CSV and component-fraction diagnostics.
    {
        std::ofstream cf(od+"/normalization_cutflow.csv");
        cf << "sample,baseline,mx2_ep,mx2_eg,dphi_pg,angle_eX,all_cuts,all_over_baseline\n";
        for (const auto& vp:vv) if (vp) {
            const double f=vp->norm_cutflow[0]>0?double(vp->norm_cutflow[5])/vp->norm_cutflow[0]:0;
            cf << vp->name; for (int i=0;i<6;i++) cf << ","<<vp->norm_cutflow[i]; cf << ","<<f<<"\n";
        }
    }
    {
        std::ofstream ff(od+"/component_fractions.csv");
        ff << "region,aao,clasdis,dvcs,total,aao_fraction,clasdis_fraction,dvcs_fraction\n";
        auto write_region=[&](const char* label,const TH1D* ha,const TH1D* hc,const TH1D* hv,double lo,double hi) {
            int ba=ha->GetXaxis()->FindBin(lo+1e-9), bb=ha->GetXaxis()->FindBin(hi-1e-9);
            const double ya=R.nominal.aao*ha->Integral(ba,bb), yc=R.nominal.clasdis*hc->Integral(ba,bb), yv=R.nominal.dvcs*hv->Integral(ba,bb);
            const double yt=ya+yc+yv; ff<<label<<","<<ya<<","<<yc<<","<<yv<<","<<yt<<","<<(yt?ya/yt:0)<<","<<(yt?yc/yt:0)<<","<<(yt?yv/yt:0)<<"\n";
        };
        const int io=NORM_EGAMMA;
        write_region("lowE",aao->norm_full[io].get(),cls->norm_full[io].get(),dvc->norm_full[io].get(),0.4,2.0);
        write_region("transition",aao->norm_full[io].get(),cls->norm_full[io].get(),dvc->norm_full[io].get(),2.0,3.0);
        write_region("highE",aao->norm_full[io].get(),cls->norm_full[io].get(),dvc->norm_full[io].get(),3.0,8.0);
        write_region("full",aao->norm_full[io].get(),cls->norm_full[io].get(),dvc->norm_full[io].get(),0.4,8.0);
    }

    // Summary CSVs.
    std::ofstream csv(od+"/normalization_fit_results.csv");
    csv << "stage,observable,valid,aao,aao_err,clasdis,clasdis_err,dvcs,dvcs_err,morph_shift,morph_sigma,raw_chi2,raw_ndf,raw_chi2_ndf,morphed_chi2,morphed_ndf,morphed_chi2_ndf\n";
    for (const auto& q:R.low_points) csv << "lowE,"<<q.observable<<","<<q.valid<<","<<q.aao<<","<<q.aao_err<<","<<q.clasdis<<","<<q.clasdis_err<<",0,0,"<<q.morph_shift<<","<<q.morph_sigma<<","<<q.raw_chi2<<","<<q.raw_ndf<<","<<(q.raw_ndf?q.raw_chi2/q.raw_ndf:0)<<","<<q.chi2<<","<<q.ndf<<","<<(q.ndf?q.chi2/q.ndf:0)<<"\n";
    for (const auto& q:R.high_points) csv << "highE,"<<q.observable<<","<<q.valid<<","<<q.aao<<",0,"<<q.clasdis<<",0,"<<q.dvcs<<","<<q.dvcs_err<<","<<q.morph_shift<<","<<q.morph_sigma<<","<<q.raw_chi2<<","<<q.raw_ndf<<","<<(q.raw_ndf?q.raw_chi2/q.raw_ndf:0)<<","<<q.chi2<<","<<q.ndf<<","<<(q.ndf?q.chi2/q.ndf:0)<<"\n";
    csv.close();

    std::ofstream fac(od+"/normalization_factors.csv");
    fac << "set,aao,clasdis,dvcs\n";
    fac << "derived_nominal,"<<R.nominal.aao<<","<<R.nominal.clasdis<<","<<R.nominal.dvcs<<"\n";
    fac << "derived_low,"<<R.low.aao<<","<<R.low.clasdis<<","<<R.low.dvcs<<"\n";
    fac << "derived_high,"<<R.high.aao<<","<<R.high.clasdis<<","<<R.high.dvcs<<"\n";
    fac << "Valerii_June_reference,0.307,0.315,1.10\n";
    fac.close();

    std::ofstream txt(od+"/normalization_summary.txt");
    txt << "Data-driven MC template normalization\n====================================\n"
        << "MC event weights: unit counts only; MC::Event.weight is not used.\n"
        << "Normalization selection: -0.231<Mx2(ep)<0.309 GeV2; Mx2(e-gamma)>1.4 GeV2; |DeltaPhi(p,gamma)|<5.7 deg; angle(gamma,X)<9.2 deg.\n"
        << "Fit histograms use an N-1 implementation of these cuts for the observable being fitted.\n"
        << "Low-E stage: AAO + CLASDIS for E_gamma<2 GeV. High-E stage: fix their means, then fit DVCS for E_gamma>3 GeV.\n"
        << "Normalization-only template morph: one common MC shift + extra Gaussian sigma per observable, scanned over +/-"<<NORM_MAX_SHIFT_BINS<<" and 0-"<<NORM_MAX_SMEAR_BINS<<" histogram bins.\n"
        << "The morph is NOT applied to event-level MC and does NOT alter the photon-efficiency numerator or denominator.\n"
        << "Nominal factors are unweighted means across valid observable fits; component-wise extrema define a conservative normalization variation.\n"
        << "The CLASDIS generated-level restriction mentioned in the June presentation is intentionally not reproduced here, per current analysis choice.\n\n"
        << "Derived nominal: AAO="<<R.nominal.aao<<" CLASDIS="<<R.nominal.clasdis<<" DVCS="<<R.nominal.dvcs<<"\n"
        << "Derived low:     AAO="<<R.low.aao<<" CLASDIS="<<R.low.clasdis<<" DVCS="<<R.low.dvcs<<"\n"
        << "Derived high:    AAO="<<R.high.aao<<" CLASDIS="<<R.high.clasdis<<" DVCS="<<R.high.dvcs<<"\n"
        << "Valerii June reference only: AAO=0.307 CLASDIS=0.315 DVCS=1.10\n"
        << "Low-E valid observable fits: "<<R.low_points.size()<<"; high-E valid observable fits: "<<R.high_points.size()<<"\n"        << "Delta-t definition: t_p - t_gamma, with t_p=(p_target-p')^2 and t_gamma=(q-gamma_tag)^2.\n"
        << "Fallback used: "<<(R.used_fallback?"YES":"NO")<<"\n";
    txt.close();

    // Summary factor plot.
    TCanvas cs("c_norm_summary","",1100,700); cs.SetGridy();
    TH1D frame("norm_summary_frame","Normalization factors by observable;Observable;Scale factor",(int)(R.low_points.size()+R.high_points.size()),0,(int)(R.low_points.size()+R.high_points.size()));
    frame.SetStats(0); frame.SetMinimum(0); double ymax=1.4*std::max({1.2,R.high.dvcs,R.high.aao,R.high.clasdis}); frame.SetMaximum(ymax);
    int ib=1; for (const auto& q:R.low_points) { frame.GetXaxis()->SetBinLabel(ib,("low "+q.observable).c_str()); frame.SetBinContent(ib,std::max(q.aao,q.clasdis)); ib++; }
    for (const auto& q:R.high_points) { frame.GetXaxis()->SetBinLabel(ib,("high "+q.observable).c_str()); frame.SetBinContent(ib,q.dvcs); ib++; }
    frame.LabelsOption("v","X"); frame.Draw("HIST TEXT"); cs.SetBottomMargin(0.30); cs.SaveAs((od+"/normalization_factor_summary.png").c_str());

    // Persist all normalization histograms for detailed offline inspection.
    TFile rf((od+"/normalization_histograms.root").c_str(),"RECREATE");
    if (!rf.IsZombie()) {
        for (const auto& vp:vv) if (vp) {
            TDirectory* d=rf.mkdir(vp->name.c_str()); if (!d) continue; d->cd();
            for (int io=0;io<NORM_NOBS;io++) {
                if (io<(int)vp->norm_pre.size() && vp->norm_pre[io]) vp->norm_pre[io]->Write(Form("pre_%s",NORM_OBS[io].key));
                if (io<(int)vp->norm_after_mx2ep.size() && vp->norm_after_mx2ep[io]) vp->norm_after_mx2ep[io]->Write(Form("after_mx2ep_%s",NORM_OBS[io].key));
                if (io<(int)vp->norm_after_mx2ep_mx2eg.size() && vp->norm_after_mx2ep_mx2eg[io])
                    vp->norm_after_mx2ep_mx2eg[io]->Write(Form("after_mx2ep_mx2eg_%s",NORM_OBS[io].key));
                if (io<(int)vp->norm_after_mx2ep_mx2eg_dphi.size() && vp->norm_after_mx2ep_mx2eg_dphi[io])
                    vp->norm_after_mx2ep_mx2eg_dphi[io]->Write(Form("after_mx2ep_mx2eg_dphi_%s",NORM_OBS[io].key));
                if (io<(int)vp->norm_nminus1.size() && vp->norm_nminus1[io])
                    vp->norm_nminus1[io]->Write(Form("nminus1_%s",NORM_OBS[io].key));
                if (io<(int)vp->norm_full.size() && vp->norm_full[io]) vp->norm_full[io]->Write(Form("full_%s",NORM_OBS[io].key));
                if (io<(int)vp->norm_lowE.size() && vp->norm_lowE[io]) vp->norm_lowE[io]->Write(Form("lowE_%s",NORM_OBS[io].key));
                if (io<(int)vp->norm_highE.size() && vp->norm_highE[io]) vp->norm_highE[io]->Write(Form("highE_%s",NORM_OBS[io].key));
            }
            rf.cd();
        }
        rf.Close();
    }
    return R;
}

void write_valerii_outputs(const std::vector<std::unique_ptr<ValComponent>>& vv,
                           const std::string& out) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) {
        std::ofstream s(out+"/valerii_fd_summary.txt");
        s << "Valerii FD reproduction not evaluated: data, AAOgen, CLASDIS, and DVCSgen are all required.\n";
        return;
    }

    // Derive the component normalizations from the actual files in this run
    // before constructing the weighted-total MC used for the efficiency maps.
    NormDerivation norm=derive_normalization(vv,out);
    std::vector<std::unique_ptr<TH1D>> hd_nom,hm_nom;
    auto nominal=evaluate_valerii_fd(vv,norm.nominal,&hd_nom,&hm_nom);
    auto set1=evaluate_valerii_fd(vv,norm.low);
    auto set2=evaluate_valerii_fd(vv,norm.high);

    for (int ib=0;ib<VAL_NBIN;ib++) {
        nominal[ib].correction_norm_set1=set1[ib].correction[1];
        nominal[ib].correction_norm_set2=set2[ib].correction[1];
        if (nominal[ib].nominal_valid) {
            const double c0=nominal[ib].correction[1];
            double sn=0;
            if (set1[ib].correction[1]>0) sn=std::max(sn,std::fabs(set1[ib].correction[1]-c0));
            if (set2[ib].correction[1]>0) sn=std::max(sn,std::fabs(set2[ib].correction[1]-c0));
            nominal[ib].norm_systematic=sn;
            double sm=0;
            if (nominal[ib].correction[0]>0) sm=std::max(sm,std::fabs(nominal[ib].correction[0]-c0));
            if (nominal[ib].correction[2]>0) sm=std::max(sm,std::fabs(nominal[ib].correction[2]-c0));
            nominal[ib].matching_systematic=sm;
            nominal[ib].partial_total_unc=std::sqrt(
                nominal[ib].correction_stat_err[1]*nominal[ib].correction_stat_err[1] + sn*sn + sm*sm);
        }
    } // endfor

    std::ofstream csv(out+"/valerii_fd_results.csv");
    csv << "bin,ip,it,iphi,p_lo,p_hi,theta_lo,theta_hi,phi_lo,phi_hi,"
        << "data_fit_valid,data_mu,data_mu_err,data_sigma,data_sigma_err,"
        << "mc_fit_valid,mc_mu,mc_mu_err,mc_sigma,mc_sigma_err,"
        << "data_denom,mc_denom,"
        << "data_eff_1s,data_eff_1s_err,mc_eff_1s,mc_eff_1s_err,corr_1s,corr_1s_staterr,"
        << "data_eff_2s,data_eff_2s_err,mc_eff_2s,mc_eff_2s_err,corr_2s,corr_2s_staterr,"
        << "data_eff_3s,data_eff_3s_err,mc_eff_3s,mc_eff_3s_err,corr_3s,corr_3s_staterr,"
        << "corr_2s_norm_set1,corr_2s_norm_set2,norm_syst,matching_syst,partial_total_unc,partial_rel_unc\n";
    for (int ib=0;ib<VAL_NBIN;ib++) {
        int ip,it,iph; val_unflatten(ib,ip,it,iph);
        const auto& r=nominal[ib];
        csv << ib << "," << ip << "," << it << "," << iph << ","
            << VAL_P_EDGES[ip] << "," << VAL_P_EDGES[ip+1] << ","
            << VAL_T_EDGES[it] << "," << VAL_T_EDGES[it+1] << ","
            << VAL_PH_EDGES[iph] << "," << VAL_PH_EDGES[iph+1] << ","
            << (r.data_fit.valid?1:0) << "," << r.data_fit.mean << "," << r.data_fit.mean_err << ","
            << r.data_fit.sigma << "," << r.data_fit.sigma_err << ","
            << (r.mc_fit.valid?1:0) << "," << r.mc_fit.mean << "," << r.mc_fit.mean_err << ","
            << r.mc_fit.sigma << "," << r.mc_fit.sigma_err << ","
            << r.data[1].denom << "," << r.mc[1].denom;
        for (int ns=0;ns<3;ns++) {
            csv << "," << r.data[ns].efficiency << "," << r.data[ns].efficiency_err
                << "," << r.mc[ns].efficiency << "," << r.mc[ns].efficiency_err
                << "," << r.correction[ns] << "," << r.correction_stat_err[ns];
        } // endfor
        const double rel=(r.nominal_valid && r.correction[1]!=0)?r.partial_total_unc/std::fabs(r.correction[1]):-1;
        csv << "," << r.correction_norm_set1 << "," << r.correction_norm_set2
            << "," << r.norm_systematic << "," << r.matching_systematic
            << "," << r.partial_total_unc << "," << rel << "\n";
    } // endfor
    csv.close();

    std::ofstream qa(out+"/valerii_fd_component_weight_qa.csv");
    qa << "component,is_mc,tree_entries,denom_rows,unit_sumw,unit_sumw2,mean_event_weight,nominal_component_scale,scaled_sumw,scaled_neff\n";
    for (const auto& vp:vv) {
        long long rows=0; double sw=0,sw2=0;
        for (const auto& b:vp->bins) { rows+=b.denom_rows; sw+=b.denom_w; sw2+=b.denom_w2; }
        const double sc=vp->is_mc?val_component_scale(vp->name,norm.nominal):1.0;
        const double ssw=sc*sw, ssw2=sc*sc*sw2;
        const double neff=ssw2>0?ssw*ssw/ssw2:0;
        const double meanw=rows>0?sw/static_cast<double>(rows):0.0;
        qa << vp->name << "," << (vp->is_mc?1:0) << "," << vp->entries << "," << rows << ","
           << sw << "," << sw2 << "," << meanw << "," << sc << "," << ssw << "," << neff << "\n";
    } // endfor
    qa.close();

    val_draw_map(nominal,"data_eff",out+"/valerii_fd_data_efficiency_2sigma.png");
    val_draw_map(nominal,"mc_eff",out+"/valerii_fd_weighted_mc_efficiency_2sigma.png");
    val_draw_map(nominal,"correction",out+"/valerii_fd_data_over_mc_correction_2sigma.png");
    val_draw_map(nominal,"data_sigma",out+"/valerii_fd_data_sigma.png");
    val_draw_map(nominal,"mc_sigma",out+"/valerii_fd_weighted_mc_sigma.png");

    TFile rf((out+"/valerii_fd_histograms.root").c_str(),"RECREATE");
    if (!rf.IsZombie()) {
        for (int ib=0;ib<VAL_NBIN;ib++) {
            int ip,it,iph; val_unflatten(ib,ip,it,iph);
            TDirectory* d=rf.mkdir(Form("bin_%03d_p%d_t%d_phi%d",ib,ip,it,iph));
            d->cd();
            if (ib<(int)hd_nom.size() && hd_nom[ib]) hd_nom[ib]->Write("data_dp");
            if (ib<(int)hm_nom.size() && hm_nom[ib]) hm_nom[ib]->Write("weighted_mc_dp");
            write_fit_state(d,"data_fit_",nominal[ib].data_fit);
            write_fit_state(d,"mc_fit_",nominal[ib].mc_fit);
        } // endfor
        rf.Close();
    }

    long long valid=0;
    double csum=0; long long cn=0;
    for (const auto& r:nominal) if (r.nominal_valid) { valid++; csum+=r.correction[1]; cn++; }
    std::ofstream summary(out+"/valerii_fd_summary.txt");
    summary << "Valerii-style FD photon-efficiency reproduction\n"
            << "==============================================\n"
            << "Binning: 7 p x 3 theta x 6 wrapped-phi = 126 bins\n"
            << "p edges (GeV): 0.35 0.50 1.10 1.70 2.30 2.90 3.70 6.00\n"
            << "theta edges (deg): 6 20 27 36\n"
            << "phi edges (deg): -30 30 90 150 210 270 330\n"
            << "Residual histograms: 60 bins on [-1,1] GeV\n"
            << "Active MC normalization is derived from this run's template fits.\n"
            << "Derived nominal: " << norm.nominal.aao << "*AAO + " << norm.nominal.clasdis
            << "*CLASDIS + " << norm.nominal.dvcs << "*DVCS\n"
            << "Normalization-variation low: (" << norm.low.aao << "," << norm.low.clasdis << "," << norm.low.dvcs << ")\n"
            << "Normalization-variation high: (" << norm.high.aao << "," << norm.high.clasdis << "," << norm.high.dvcs << ")\n"
            << "Valerii June reference only: (0.307,0.315,1.10)\n"
            << "MC events are unit weighted inside each component; only the component normalization constants are applied.\n"
            << "The skim MC::Event.weight branch is NOT used in the Valerii FD calculation.\n"
            << "Tag photon: FD/PCAL only; probe coordinates: missing gamma2.\n"
            << "June normalization/exclusivity cuts are applied to the Valerii denominator; the older stage-1 0.08<Mx(ep)<0.20 GeV development window is not used.\n"
            << "Normalization template morphing is QA/calibration only and does not smear event-level efficiency MC.\n"
            << "Numerator photon threshold remains the analysis skim threshold p>=0.4 GeV.\n\n"
            << "Nominal 2-sigma bins with valid data and weighted-MC fits: " << valid << " / " << VAL_NBIN << "\n";
    if (cn>0) summary << "Unweighted mean correction across valid bins (diagnostic only): " << csum/cn << "\n";
    summary << "\nIMPORTANT: partial_total_unc currently combines statistical, normalization-set, and\n"
            << "1/2/3-sigma matching-window variations only.  The PCAL/DC/SF configuration\n"
            << "variation from Valerii is not yet reproducible from the present skim, so the\n"
            << "final >30% reliability/neutral-correction rule is intentionally NOT applied yet.\n"
            << "The fit model remains Gaussian + linear background; replace it if the original\n"
            << "production fitter source establishes a different functional form.\n";
    summary.close();
}


const std::vector<std::unique_ptr<TH1D>>& norm_pre_for(const ValComponent* v,bool ft) {
    return ft ? v->norm_ft_pre : v->norm_pre;
}
const std::vector<std::unique_ptr<TH1D>>& norm_after_ep_for(const ValComponent* v,bool ft) {
    return ft ? v->norm_ft_after_mx2ep : v->norm_after_mx2ep;
}
const std::vector<std::unique_ptr<TH1D>>& norm_after_eg_for(const ValComponent* v,bool ft) {
    return ft ? v->norm_ft_after_mx2ep_mx2eg : v->norm_after_mx2ep_mx2eg;
}
const std::vector<std::unique_ptr<TH1D>>& norm_after_dphi_for(const ValComponent* v,bool ft) {
    return ft ? v->norm_ft_after_mx2ep_mx2eg_dphi : v->norm_after_mx2ep_mx2eg_dphi;
}
const std::vector<std::unique_ptr<TH1D>>& norm_low_for(const ValComponent* v,bool ft) {
    return ft ? v->norm_ft_lowE : v->norm_lowE;
}
const std::vector<std::unique_ptr<TH1D>>& norm_high_for(const ValComponent* v,bool ft) {
    return ft ? v->norm_ft_highE : v->norm_highE;
}
const std::vector<std::unique_ptr<TH1D>>& norm_full_for(const ValComponent* v,bool ft) {
    return ft ? v->norm_ft_full : v->norm_full;
}
const std::array<long long,6>& norm_cutflow_for(const ValComponent* v,bool ft) {
    return ft ? v->norm_ft_cutflow : v->norm_cutflow;
}

double rms_spread(const std::vector<double>& x,double mean) {
    if (x.size()<2) return 0.0;
    double s=0.0;
    for (double v:x) s+=(v-mean)*(v-mean);
    return std::sqrt(s/double(x.size()-1));
}


struct SimMorphState {
    int io=-1;
    double shift=0;
    double sigma=0;
};

NormFitPoint fit_two_templates_simultaneous_morphed(
        const std::vector<std::unique_ptr<TH1D>>& data,
        const std::vector<std::unique_ptr<TH1D>>& aao,
        const std::vector<std::unique_ptr<TH1D>>& cls,
        const std::vector<int>& obs) {
    NormFitPoint out;
    out.observable="simultaneous_lowE";
    if (obs.empty()) return out;

    std::vector<SimMorphState> states;
    states.reserve(obs.size());

    // Start each observable from its individually preferred common MC morph.
    for (int io:obs) {
        if (io<0 || io>=NORM_NOBS ||
            io>=int(data.size()) || io>=int(aao.size()) || io>=int(cls.size()) ||
            !data[io] || !aao[io] || !cls[io]) continue;
        auto q=fit_two_templates_morphed(data[io].get(),aao[io].get(),cls[io].get(),
                                         NORM_OBS[io].key);
        SimMorphState s;
        s.io=io;
        if (q.valid) {
            s.shift=q.morph_shift;
            s.sigma=q.morph_sigma;
        } // endif
        states.push_back(s);
    } // endfor
    if (states.size()<2) return out;

    auto solve_common = [&](double& A,double& B,double& eA,double& eB,
                            double& chi2,int& used)->bool {
        double saa=0,sbb=0,sab=0,sad=0,sbd=0;
        used=0;
        for (const auto& s:states) {
            auto ma=morph_norm_hist(aao[s.io].get(),s.shift,s.sigma,"sim_ma");
            auto mc=morph_norm_hist(cls[s.io].get(),s.shift,s.sigma,"sim_mc");
            const TH1D* hd=data[s.io].get();
            if (!ma || !mc || !hd) continue;
            for (int ib=1;ib<=hd->GetNbinsX();ib++) {
                const double d=hd->GetBinContent(ib);
                const double va=ma->GetBinContent(ib);
                const double vb=mc->GetBinContent(ib);
                if (d<=0 && va<=0 && vb<=0) continue;
                const double w=1.0/std::max(1.0,d);
                saa+=w*va*va;
                sbb+=w*vb*vb;
                sab+=w*va*vb;
                sad+=w*va*d;
                sbd+=w*vb*d;
                used++;
            } // endfor
        } // endfor

        const double det=saa*sbb-sab*sab;
        if (used<10 || !(det>0)) return false;
        A=(sad*sbb-sbd*sab)/det;
        B=(sbd*saa-sad*sab)/det;

        // Enforce physical non-negative component normalizations.
        if (A<0) {
            A=0;
            B=(sbb>0?sbd/sbb:0);
        } // endif
        if (B<0) {
            B=0;
            A=(saa>0?sad/saa:0);
        } // endif
        if (!(A>=0 && B>=0) || !std::isfinite(A) || !std::isfinite(B))
            return false;

        eA=(det>0?std::sqrt(std::max(0.0,sbb/det)):0);
        eB=(det>0?std::sqrt(std::max(0.0,saa/det)):0);

        chi2=0;
        for (const auto& s:states) {
            auto ma=morph_norm_hist(aao[s.io].get(),s.shift,s.sigma,"sim_ma_chi");
            auto mc=morph_norm_hist(cls[s.io].get(),s.shift,s.sigma,"sim_mc_chi");
            const TH1D* hd=data[s.io].get();
            if (!ma || !mc || !hd) continue;
            for (int ib=1;ib<=hd->GetNbinsX();ib++) {
                const double d=hd->GetBinContent(ib);
                const double m=A*ma->GetBinContent(ib)+B*mc->GetBinContent(ib);
                if (d<=0 && m<=0) continue;
                const double q=d-m;
                chi2+=q*q/std::max(1.0,d);
            } // endfor
        } // endfor
        return true;
    };

    double A=0,B=0,eA=0,eB=0,chi2=0;
    int used=0;
    if (!solve_common(A,B,eA,eB,chi2,used)) return out;

    // Alternate between one common normalization solution and an independent
    // detector-resolution morph for each observable.  The morph changes shape
    // only; A and B remain common to every observable.
    for (int iter=0;iter<3;iter++) {
        for (auto& s:states) {
            const TH1D* hd=data[s.io].get();
            const double bw=hd->GetXaxis()->GetBinWidth(1);
            double bestchi=std::numeric_limits<double>::infinity();
            double bestshift=s.shift,bestsigma=s.sigma;

            for (double sb=-NORM_MAX_SHIFT_BINS;
                 sb<=NORM_MAX_SHIFT_BINS+1e-9;
                 sb+=NORM_MORPH_STEP_BINS) {
                const double shift=sb*bw;
                for (double wb=0;
                     wb<=NORM_MAX_SMEAR_BINS+1e-9;
                     wb+=NORM_MORPH_STEP_BINS) {
                    const double sigma=wb*bw;
                    auto ma=morph_norm_hist(aao[s.io].get(),shift,sigma,"sim_scan_a");
                    auto mc=morph_norm_hist(cls[s.io].get(),shift,sigma,"sim_scan_c");
                    if (!ma || !mc) continue;
                    double c2=0;
                    for (int ib=1;ib<=hd->GetNbinsX();ib++) {
                        const double d=hd->GetBinContent(ib);
                        const double m=A*ma->GetBinContent(ib)+B*mc->GetBinContent(ib);
                        if (d<=0 && m<=0) continue;
                        const double q=d-m;
                        c2+=q*q/std::max(1.0,d);
                    } // endfor
                    if (c2<bestchi) {
                        bestchi=c2;
                        bestshift=shift;
                        bestsigma=sigma;
                    } // endif
                } // endfor
            } // endfor
            s.shift=bestshift;
            s.sigma=bestsigma;
        } // endfor

        if (!solve_common(A,B,eA,eB,chi2,used)) return out;
    } // endfor

    out.aao=A;
    out.clasdis=B;
    out.aao_err=eA;
    out.clasdis_err=eB;
    out.chi2=chi2;
    // Two common normalizations + two shape nuisances per observable.
    out.ndf=std::max(1,used-2-2*int(states.size()));
    out.valid=std::isfinite(A) && std::isfinite(B) && out.ndf>0;
    return out;
}

NormFitPoint fit_one_template_simultaneous_morphed(
        const std::vector<std::unique_ptr<TH1D>>& data,
        const std::vector<std::unique_ptr<TH1D>>& templ,
        const std::vector<int>& obs,
        const std::string& label) {
    NormFitPoint out;
    out.observable=label;
    if (obs.empty()) return out;

    struct OneState { int io=-1; double shift=0,sigma=0; };
    std::vector<OneState> states;
    for (int io:obs) {
        if (io<0 || io>=NORM_NOBS ||
            io>=int(data.size()) || io>=int(templ.size()) ||
            !data[io] || !templ[io]) continue;
        states.push_back({io,0,0});
    } // endfor
    if (states.size()<2) return out;

    auto solveA=[&](double& A,double& eA,double& chi2,int& used)->bool {
        double sxx=0,sxd=0;
        used=0;
        for (const auto& s:states) {
            auto mt=morph_norm_hist(templ[s.io].get(),s.shift,s.sigma,"one_mt");
            const TH1D* hd=data[s.io].get();
            if (!mt || !hd) continue;
            for (int ib=1;ib<=hd->GetNbinsX();ib++) {
                const double d=hd->GetBinContent(ib), x=mt->GetBinContent(ib);
                if (d<=0 && x<=0) continue;
                const double w=1.0/std::max(1.0,d);
                sxx+=w*x*x;
                sxd+=w*x*d;
                used++;
            } // endfor
        } // endfor
        if (!(sxx>0) || used<10) return false;
        A=std::max(0.0,sxd/sxx);
        eA=std::sqrt(1.0/sxx);
        chi2=0;
        for (const auto& s:states) {
            auto mt=morph_norm_hist(templ[s.io].get(),s.shift,s.sigma,"one_mt_chi");
            const TH1D* hd=data[s.io].get();
            for (int ib=1;ib<=hd->GetNbinsX();ib++) {
                const double d=hd->GetBinContent(ib), m=A*mt->GetBinContent(ib);
                if (d<=0 && m<=0) continue;
                const double q=d-m;
                chi2+=q*q/std::max(1.0,d);
            } // endfor
        } // endfor
        return true;
    };

    double A=0,eA=0,chi2=0;
    int used=0;
    if (!solveA(A,eA,chi2,used)) return out;

    for (int iter=0;iter<3;iter++) {
        for (auto& s:states) {
            const TH1D* hd=data[s.io].get();
            const double bw=hd->GetXaxis()->GetBinWidth(1);
            double bestchi=std::numeric_limits<double>::infinity();
            double bestshift=s.shift,bestsigma=s.sigma;
            for (double sb=-NORM_MAX_SHIFT_BINS;
                 sb<=NORM_MAX_SHIFT_BINS+1e-9;
                 sb+=NORM_MORPH_STEP_BINS) {
                const double shift=sb*bw;
                for (double wb=0;
                     wb<=NORM_MAX_SMEAR_BINS+1e-9;
                     wb+=NORM_MORPH_STEP_BINS) {
                    const double sigma=wb*bw;
                    auto mt=morph_norm_hist(templ[s.io].get(),shift,sigma,"one_scan");
                    if (!mt) continue;
                    double c2=0;
                    for (int ib=1;ib<=hd->GetNbinsX();ib++) {
                        const double d=hd->GetBinContent(ib);
                        const double m=A*mt->GetBinContent(ib);
                        if (d<=0 && m<=0) continue;
                        const double q=d-m;
                        c2+=q*q/std::max(1.0,d);
                    } // endfor
                    if (c2<bestchi) {
                        bestchi=c2;
                        bestshift=shift;
                        bestsigma=sigma;
                    } // endif
                } // endfor
            } // endfor
            s.shift=bestshift;
            s.sigma=bestsigma;
        } // endfor
        if (!solveA(A,eA,chi2,used)) return out;
    } // endfor

    out.aao=A;
    out.aao_err=eA;
    out.chi2=chi2;
    out.ndf=std::max(1,used-1-2*int(states.size()));
    out.valid=std::isfinite(A) && out.ndf>0;
    return out;
}

NormDerivation derive_normalization_concise(
        const std::vector<std::unique_ptr<ValComponent>>& vv,
        bool ft_probe) {
    NormDerivation R;
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");

    if (!data || !aao || !cls || !dvc) {
        R.used_fallback=true;
        R.nominal=VAL_HISTORICAL_NOMINAL;
        R.low=R.nominal;
        R.high=R.nominal;
        return R;
    } // endif

    const auto& dl=norm_low_for(data,ft_probe);
    const auto& al=norm_low_for(aao,ft_probe);
    const auto& cl=norm_low_for(cls,ft_probe);

    // Keep every Valerii-style one-observable result for QA.  These are no
    // longer averaged to determine the nominal AAOgen/CLASDIS normalization.
    for (int io=0;io<NORM_NOBS;io++) {
        if (!NORM_OBS[io].use_low) continue;
        auto q=fit_two_templates_morphed(
            dl[io].get(),al[io].get(),cl[io].get(),NORM_OBS[io].key);
        if (q.valid) R.low_points.push_back(q);
    } // endfor

    // Nominal common low-E normalization.  angle(gamma,X) is deliberately a
    // validation observable here: its baseline data shape contains a visibly
    // DVCS-like shoulder.  The low-E Valerii step assumes DVCS is negligible,
    // so forcing AAOgen/CLASDIS alone to absorb that shoulder biases B upward.
    // Once the high-E DVCS normalization is known, DVCSgen is overlaid on the
    // low-E panels as a closure check.
    const std::vector<int> low_global_obs = {
        NORM_MX2_EP,
        NORM_MX2_EPG,
        NORM_EGAMMA,
        NORM_MX2_EG
    };
    R.low_global=fit_two_templates_simultaneous_morphed(
        dl,al,cl,low_global_obs);

    if (!R.low_global.valid) {
        R.used_fallback=true;
        R.nominal=VAL_HISTORICAL_NOMINAL;
        R.low=R.nominal;
        R.high=R.nominal;
        return R;
    } // endif

    const double A=R.low_global.aao;
    const double B=R.low_global.clasdis;

    for (int io=0;io<NORM_NOBS;io++) {
        if (!NORM_OBS[io].use_high) continue;
        const auto& dh=norm_high_for(data,ft_probe);
        const auto& ah=norm_high_for(aao,ft_probe);
        const auto& ch=norm_high_for(cls,ft_probe);
        const auto& vh=norm_high_for(dvc,ft_probe);
        auto q=fit_dvcs_template_morphed(
            dh[io].get(),ah[io].get(),ch[io].get(),vh[io].get(),
            A,B,NORM_OBS[io].key);
        if (q.valid) R.high_points.push_back(q);
    } // endfor

    const double C=mean_valid(R.high_points,2);
    R.valid=(A>=0 && B>=0 && C>=0 &&
             R.low_global.valid && R.high_points.size()>=2);

    if (!R.valid) {
        R.used_fallback=true;
        R.nominal=VAL_HISTORICAL_NOMINAL;
        R.low=R.nominal;
        R.high=R.nominal;
        return R;
    } // endif

    R.nominal={"derived_simultaneous",A,B,C};

    // Normalization envelope: use the spread of the included low-E diagnostic
    // fits around the simultaneous solution.  angle(gamma,X) remains a closure
    // diagnostic and is not allowed to dominate this envelope.
    std::vector<double> av,bv,cv;
    for (const auto& q:R.low_points) {
        if (q.observable=="angle_gX") continue;
        av.push_back(q.aao);
        bv.push_back(q.clasdis);
    } // endfor
    for (const auto& q:R.high_points) cv.push_back(q.dvcs);

    const double sa=rms_spread(av,A);
    const double sb=rms_spread(bv,B);
    const double sc=rms_spread(cv,C);

    R.low={"simultaneous_minus_spread",
           std::max(0.0,A-sa),std::max(0.0,B-sb),std::max(0.0,C-sc)};
    R.high={"simultaneous_plus_spread",A+sa,B+sb,C+sc};

    // FT-only alternate model requested for the study: remove CLASDIS, refit
    // AAOgen simultaneously at low tag energy, then refit DVCSgen at high tag
    // energy with CLASDIS fixed exactly to zero.
    if (ft_probe) {
        const auto altA=fit_one_template_simultaneous_morphed(
            dl,al,low_global_obs,"simultaneous_lowE_AAO_only");
        if (altA.valid) {
            std::vector<NormFitPoint> alt_high;
            for (int io=0;io<NORM_NOBS;io++) {
                if (!NORM_OBS[io].use_high) continue;
                const auto& dh=norm_high_for(data,true);
                const auto& ah=norm_high_for(aao,true);
                const auto& ch=norm_high_for(cls,true);
                const auto& vh=norm_high_for(dvc,true);
                auto q=fit_dvcs_template_morphed(
                    dh[io].get(),ah[io].get(),ch[io].get(),vh[io].get(),
                    altA.aao,0.0,NORM_OBS[io].key);
                if (q.valid) alt_high.push_back(q);
            } // endfor
            const double Calt=mean_valid(alt_high,2);
            if (Calt>=0 && alt_high.size()>=2) {
                R.no_clasdis_valid=true;
                R.no_clasdis={"FT_AAO_plus_DVCS_only",altA.aao,0.0,Calt};
            } // endif
        } // endif
    } // endif

    return R;
}

const CoarseValRegion* coarse_region(const ValComponent* v,int ir) {
    if (!v || ir<0 || ir>=CR_N) return nullptr;
    return &v->coarse[ir];
}

double coarse_clasdis_pi0_fraction(const ValComponent* cls,int ir,double fallback=0.9877) {
    if (!cls || ir<0 || ir>=CR_N) return fallback;
    const auto& r=cls->coarse[ir];
    if (r.truth_rows>0)
        return double(r.truth_pi0_rows)/double(r.truth_rows);
    return fallback;
}

struct CoarseComposition {
    double f_pi0=-1;
    double ya=0,yc=0,yd=0;
    double f_clasdis_pi0=0;
};

CoarseComposition coarse_composition(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                     const ValNormSet& n,int ir) {
    CoarseComposition r;
    const ValComponent* a=find_val_component(vv,"aaogen");
    const ValComponent* c=find_val_component(vv,"clasdis");
    const ValComponent* d=find_val_component(vv,"dvcsgen");
    if (!a || !c || !d) return r;

    const double fc=coarse_clasdis_pi0_fraction(c,ir);
    r.f_clasdis_pi0=fc;
    r.ya=n.aao*double(a->coarse[ir].denom_rows);
    r.yc=n.clasdis*double(c->coarse[ir].denom_rows);
    r.yd=n.dvcs*double(d->coarse[ir].denom_rows);

    const double total=r.ya+r.yc+r.yd;
    const double pi0=r.ya+fc*r.yc;
    if (total>0) r.f_pi0=pi0/total;
    return r;
}

void concise_make_dirs(const std::string& out) {
    reset_output(out);
    gSystem->mkdir((out+"/1_exclusivity/FD").c_str(),true);
    gSystem->mkdir((out+"/1_exclusivity/FT").c_str(),true);
    gSystem->mkdir((out+"/2_normalization/FD").c_str(),true);
    gSystem->mkdir((out+"/2_normalization/FT").c_str(),true);
    gSystem->mkdir((out+"/3_pi0_fraction").c_str(),true);
    gSystem->mkdir((out+"/4_efficiency").c_str(),true);
}


std::string corr_axis_title(int io) {
    if (io==NORM_MX2_EP) return "M_{X}^{2}(ep) (GeV^{2})";
    if (io==NORM_MX2_EPG) return "M_{X}^{2}(ep#gamma) (GeV^{2})";
    if (io==NORM_MX2_EG) return "M_{X}^{2}(e#gamma) (GeV^{2})";
    if (io==NORM_DPHI_TRENTO_SHIFT180) return "#Delta#phi_{copl} (deg)";
    if (io==NORM_ANGLE_GX) return "angle(#gamma,X) (deg)";
    if (io==NORM_ANGLE_EX) return "angle(e,X) (deg)";
    if (io==NORM_DELTA_T_PG) return "#Delta t=t_{p}-t_{#gamma} (GeV^{2})";
    return NORM_OBS[io].key;
}

void draw_corr_cut_lines(int io,double xmin,double xmax) {
    auto hline=[&](double y) {
        TLine* l=new TLine(xmin,y,xmax,y);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw();
    };
    if (io==NORM_MX2_EP) {
        hline(NORM_MX2_EP_MIN);
        hline(NORM_MX2_EP_MAX);
    } else if (io==NORM_MX2_EG) {
        hline(NORM_MX2_EG_MIN);
    } else if (io==NORM_DPHI_TRENTO_SHIFT180) {
        hline(-NORM_DPHI_TRENTO_MAX);
        hline(+NORM_DPHI_TRENTO_MAX);
    } else if (io==NORM_ANGLE_GX) {
        hline(NORM_ANGLE_GX_MAX);
    } // endif
}

std::unique_ptr<TH2D> conditionalize_corr(const TH2D* src,const char* name) {
    if (!src) return nullptr;
    std::unique_ptr<TH2D> h(dynamic_cast<TH2D*>(src->Clone(name)));
    if (!h) return nullptr;
    h->SetDirectory(nullptr);
    h->SetStats(0);
    for (int ix=1;ix<=h->GetNbinsX();ix++) {
        double s=0;
        for (int iy=1;iy<=h->GetNbinsY();iy++)
            s+=h->GetBinContent(ix,iy);
        if (!(s>0)) continue;
        for (int iy=1;iy<=h->GetNbinsY();iy++) {
            h->SetBinContent(ix,iy,h->GetBinContent(ix,iy)/s);
            h->SetBinError(ix,iy,h->GetBinError(ix,iy)/s);
        } // endfor
    } // endfor
    return h;
}

void draw_exclusivity_correlations(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                   const std::string& dir,
                                   bool ft_probe) {
    const char* names[]={"data","aaogen","clasdis","dvcsgen"};
    const char* labels[]={"Data","AAOgen (#pi^{0})","CLASDIS","DVCSgen"};
    std::vector<std::unique_ptr<TH2D>> keep;

    TCanvas c(ft_probe?"c_corr_ft":"c_corr_fd","",2300,3250);
    c.Divide(4,CORR_NOBS,0.001,0.001);

    for (int jc=0;jc<CORR_NOBS;jc++) {
        const int io=CORR_OBS[jc];
        for (int is=0;is<4;is++) {
            const ValComponent* v=find_val_component(vv,names[is]);
            if (!v) continue;
            const auto& vec=ft_probe?v->corr_ft:v->corr_fd;
            if (jc>=int(vec.size()) || !vec[jc]) continue;

            auto h=conditionalize_corr(vec[jc].get(),
                Form("corr_draw_%s_%d_%d",ft_probe?"ft":"fd",jc,is));
            if (!h) continue;

            const int ipad=jc*4+is+1;
            c.cd(ipad);
            gPad->SetLeftMargin(is==0?0.20:0.13);
            gPad->SetRightMargin(0.035);
            gPad->SetBottomMargin(jc==CORR_NOBS-1?0.18:0.10);
            gPad->SetTopMargin(jc==0?0.18:0.07);
            gPad->SetTicks(1,1);

            h->SetTitle("");
            h->GetXaxis()->SetTitle(jc==CORR_NOBS-1?"E_{#gamma,probe} (GeV)":"");
            h->GetYaxis()->SetTitle(is==0?corr_axis_title(io).c_str():"");
            h->GetXaxis()->SetTitleSize(0.060);
            h->GetYaxis()->SetTitleSize(0.055);
            h->GetXaxis()->SetLabelSize(0.045);
            h->GetYaxis()->SetLabelSize(0.045);
            h->GetYaxis()->SetTitleOffset(1.55);

            // Compact display windows only; the stored correlations retain the
            // full histogram ranges.
            if (io==NORM_ANGLE_GX) h->GetYaxis()->SetRangeUser(0,20);
            if (io==NORM_ANGLE_EX) h->GetYaxis()->SetRangeUser(0,30);
            if (io==NORM_DPHI_TRENTO_SHIFT180) h->GetYaxis()->SetRangeUser(-30,30);
            if (io==NORM_DELTA_T_PG) h->GetYaxis()->SetRangeUser(-2,2);

            h->Draw("COL");
            draw_corr_cut_lines(io,EPROBE_MIN,EPROBE_MAX);

            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            if (jc==0) {
                tx.SetTextSize(0.062);
                tx.DrawLatex(0.17,0.92,labels[is]);
            } // endif
            if (is==3) {
                tx.SetTextSize(0.047);
                tx.DrawLatex(0.67,0.88,
                    Form("#rho = %.2f",vec[jc]->GetCorrelationFactor()));
            } // endif

            keep.push_back(std::move(h));
        } // endfor
    } // endfor

    c.SaveAs((dir+"/exclusivity_correlations_vs_Eprobe.png").c_str());
}

void draw_eprobe_cut_survival(const std::vector<std::unique_ptr<ValComponent>>& vv,
                              const std::string& dir,
                              bool ft_probe) {
    const char* names[]={"data","aaogen","clasdis","dvcsgen"};
    const char* labels[]={"Data","AAOgen (#pi^{0})","CLASDIS","DVCSgen"};
    const int cols[]={kBlack,kBlue+1,kMagenta+1,kRed+1,kGreen+2};
    const char* stage_labels[]={
        "baseline",
        "M_{X}^{2}(ep)",
        "M_{X}^{2}(e#gamma)",
        "coplanarity",
        "angle(#gamma,X)"
    };

    TCanvas c(ft_probe?"c_surv_ft":"c_surv_fd","",2200,1200);
    c.Divide(4,2,0.001,0.001);
    std::vector<std::unique_ptr<TH1D>> keep;

    auto make_ratio = [&](const TH1D* num,const TH1D* den,const char* name)->std::unique_ptr<TH1D> {
        if (!num || !den) return nullptr;
        std::unique_ptr<TH1D> h(dynamic_cast<TH1D*>(num->Clone(name)));
        if (!h) return nullptr;
        h->SetDirectory(nullptr);
        h->Reset("ICES");
        h->SetStats(0);

        for (int ib=1;ib<=num->GetNbinsX();ib++) {
            const double n=num->GetBinContent(ib);
            const double d=den->GetBinContent(ib);
            if (!(d>0)) continue;
            const double p=std::max(0.0,std::min(1.0,n/d));
            h->SetBinContent(ib,p);
            h->SetBinError(ib,std::sqrt(std::max(0.0,p*(1.0-p)/d)));
        } // endfor
        return h;
    };

    for (int isamp=0;isamp<4;isamp++) {
        const ValComponent* v=find_val_component(vv,names[isamp]);
        if (!v) continue;
        const auto& st=ft_probe?v->eprobe_stage_ft:v->eprobe_stage_fd;
        if (st.size()<EPROBE_NSTAGE || !st[0]) continue;

        // Top row: cumulative survival relative to detector-accepted baseline.
        c.cd(isamp+1);
        gPad->SetLeftMargin(isamp==0?0.16:0.11);
        gPad->SetRightMargin(0.03);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.14);
        gPad->SetTicks(1,1);

        TH1D frame(Form("surv_frame_%d",isamp),
                   ";E_{#gamma,probe} (GeV);Cumulative survival",
                   EPROBE_NBIN,EPROBE_MIN,EPROBE_MAX);
        frame.SetStats(0);
        frame.SetMinimum(0);
        frame.SetMaximum(1.05);
        frame.GetXaxis()->SetTitleSize(0.050);
        frame.GetYaxis()->SetTitleSize(0.050);
        frame.GetYaxis()->SetTitleOffset(isamp==0?1.35:1.0);
        frame.Draw("AXIS");

        std::vector<TH1D*> cumulative;
        for (int ist=0;ist<EPROBE_NSTAGE;ist++) {
            auto h=make_ratio(st[ist].get(),st[0].get(),
                              Form("surv_cum_%d_%d",isamp,ist));
            if (!h) continue;
            h->SetLineColor(cols[ist]);
            h->SetMarkerColor(cols[ist]);
            h->SetMarkerStyle(20+ist);
            h->SetMarkerSize(0.45);
            h->SetLineWidth(2);
            h->Draw("E1 SAME");
            cumulative.push_back(h.get());
            keep.push_back(std::move(h));
        } // endfor

        TLatex tx;
        tx.SetNDC();
        tx.SetTextFont(42);
        tx.SetTextSize(0.052);
        tx.DrawLatex(0.15,0.92,labels[isamp]);

        if (isamp==0) {
            TLegend* leg=new TLegend(0.40,0.50,0.96,0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.031);
            for (int ist=0;ist<int(cumulative.size());ist++)
                leg->AddEntry(cumulative[ist],stage_labels[ist],"lep");
            leg->Draw();
        } // endif

        // Bottom row: incremental efficiency of each cut relative to the stage
        // immediately before it.  This isolates which cut sculpts Eprobe.
        c.cd(4+isamp+1);
        gPad->SetLeftMargin(isamp==0?0.16:0.11);
        gPad->SetRightMargin(0.03);
        gPad->SetBottomMargin(0.17);
        gPad->SetTopMargin(0.10);
        gPad->SetTicks(1,1);

        TH1D iframe(Form("surv_iframe_%d",isamp),
                    ";E_{#gamma,probe} (GeV);Incremental cut efficiency",
                    EPROBE_NBIN,EPROBE_MIN,EPROBE_MAX);
        iframe.SetStats(0);
        iframe.SetMinimum(0);
        iframe.SetMaximum(1.05);
        iframe.GetXaxis()->SetTitleSize(0.050);
        iframe.GetYaxis()->SetTitleSize(0.050);
        iframe.GetYaxis()->SetTitleOffset(isamp==0?1.35:1.0);
        iframe.Draw("AXIS");

        std::vector<TH1D*> incremental;
        for (int ist=1;ist<EPROBE_NSTAGE;ist++) {
            auto h=make_ratio(st[ist].get(),st[ist-1].get(),
                              Form("surv_inc_%d_%d",isamp,ist));
            if (!h) continue;
            h->SetLineColor(cols[ist]);
            h->SetMarkerColor(cols[ist]);
            h->SetMarkerStyle(20+ist);
            h->SetMarkerSize(0.45);
            h->SetLineWidth(2);
            h->Draw("E1 SAME");
            incremental.push_back(h.get());
            keep.push_back(std::move(h));
        } // endfor

        if (isamp==0) {
            TLegend* leg=new TLegend(0.39,0.53,0.96,0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.031);
            for (int ist=1;ist<EPROBE_NSTAGE;ist++)
                if (ist-1<int(incremental.size()))
                    leg->AddEntry(incremental[ist-1],stage_labels[ist],"lep");
            leg->Draw();
        } // endif
    } // endfor

    c.SaveAs((dir+"/Eprobe_cut_survival.png").c_str());
}

struct FTCutScanPoint {
    int iv=-1;
    double cut=0;
    double pi0_purity=0;
    double two_gen_purity=0;
    double data_retention=0;
    double aao_retention=0;
    double aao_high_retention=0;
    double eprobe_distortion=0;
    double y_data=0;
    double y_aao=0;
    double y_clasdis=0;
    double y_dvcs=0;
};

double cutscan_value(int iv,int ip) {
    const double f=(CUTSCAN_NPOINT>1)?double(ip)/double(CUTSCAN_NPOINT-1):0.0;
    if (iv==0) return 0.03 + f*(0.25-0.03);  // |Mx2(epgamma)|
    if (iv==1) return 2.0  + f*(12.0-2.0);   // |Delta phi|
    if (iv==2) return 0.08 + f*(0.40-0.08);  // Mx2(ep) half-width
    if (iv==3) return 4.0  + f*(16.0-4.0);   // angle(gamma,X)
    return 0;
}

double cutscan_nominal_value(int iv) {
    if (iv==0) return 0.25; // loose skim limit; no tighter active cut yet
    if (iv==1) return NORM_DPHI_TRENTO_MAX;
    if (iv==2) return 0.5*(NORM_MX2_EP_MAX-NORM_MX2_EP_MIN);
    if (iv==3) return NORM_ANGLE_GX_MAX;
    return 0;
}

bool cutscan_pass(int iv,double y,double cut) {
    if (!std::isfinite(y)) return false;
    if (iv==0) return std::fabs(y)<cut;
    if (iv==1) return std::fabs(y)<cut;
    if (iv==2) return std::fabs(y-CUTSCAN_MX2EP_CENTER)<cut;
    if (iv==3) return y<cut;
    return false;
}

std::unique_ptr<TH1D> cutscan_project_eprobe(const TH2D* h2,int iv,double cut,const char* name) {
    if (!h2) return nullptr;
    std::unique_ptr<TH1D> h(new TH1D(name,";E_{#gamma,probe} (GeV);Candidates",
                                     EPROBE_NBIN,EPROBE_MIN,EPROBE_MAX));
    h->Sumw2();
    h->SetDirectory(nullptr);

    for (int ix=1;ix<=h2->GetNbinsX();ix++) {
        double s=0,s2=0;
        for (int iy=1;iy<=h2->GetNbinsY();iy++) {
            const double y=h2->GetYaxis()->GetBinCenter(iy);
            if (!cutscan_pass(iv,y,cut)) continue;
            s += h2->GetBinContent(ix,iy);
            const double e=h2->GetBinError(ix,iy);
            s2 += e*e;
        } // endfor
        h->SetBinContent(ix,s);
        h->SetBinError(ix,std::sqrt(std::max(0.0,s2)));
    } // endfor
    return h;
}

double hist_integral_range(const TH1D* h,double lo,double hi) {
    if (!h) return 0;
    double s=0;
    for (int ib=1;ib<=h->GetNbinsX();ib++) {
        const double x=h->GetBinCenter(ib);
        if (x>=lo && x<hi) s+=h->GetBinContent(ib);
    } // endfor
    return s;
}

// Shape-only difference between two Eprobe spectra.  Zero means the cut leaves
// the normalized energy spectrum unchanged; larger values mean more sculpting.
double eprobe_total_variation(const TH1D* selected,const TH1D* baseline) {
    if (!selected || !baseline) return 0;
    const double ns=selected->Integral();
    const double nb=baseline->Integral();
    if (!(ns>0) || !(nb>0)) return 0;

    double tv=0;
    for (int ib=1;ib<=selected->GetNbinsX();ib++) {
        const double ps=selected->GetBinContent(ib)/ns;
        const double pb=baseline->GetBinContent(ib)/nb;
        tv += std::fabs(ps-pb);
    } // endfor
    return 0.5*tv;
}

double selected_clasdis_pi0_fraction(const TH1D* selected,const ValComponent* cls) {
    if (!selected || !cls) return 1.0;
    const double nl=hist_integral_range(selected,EPROBE_MIN,2.0);
    const double nh=hist_integral_range(selected,2.0,EPROBE_MAX+1e-6);
    const double nt=nl+nh;
    if (!(nt>0)) return 1.0;
    const double fl=coarse_clasdis_pi0_fraction(cls,CR_FT_LOW);
    const double fh=coarse_clasdis_pi0_fraction(cls,CR_FT_HIGH);
    return (nl*fl+nh*fh)/nt;
}


std::unique_ptr<TH1D> project_fitqual_observable(const ValComponent* v,
                                                bool high_tag,
                                                int ic,int io,
                                                const char* name) {
    if (!v || ic<0 || ic>=FT_COPL_SCAN_N || io<0 || io>=NORM_NOBS)
        return nullptr;
    const auto& src=high_tag?v->fitqual_ft_high:v->fitqual_ft_low;
    if (src.size()<=size_t(ic) || src[ic].size()<=size_t(io) || !src[ic][io])
        return nullptr;

    TH2D* h2=src[ic][io].get();
    std::unique_ptr<TH1D> h(h2->ProjectionY(name,1,h2->GetNbinsX(),"e"));
    if (h) h->SetDirectory(nullptr);
    return h;
}

void draw_ft_cut_optimization(const std::vector<std::unique_ptr<ValComponent>>& vv,
                              const NormDerivation& R,
                              const std::string& dir) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;
    if (data->cutscan_ft.size()<CUTSCAN_NVAR ||
        aao->cutscan_ft.size()<CUTSCAN_NVAR ||
        cls->cutscan_ft.size()<CUTSCAN_NVAR ||
        dvc->cutscan_ft.size()<CUTSCAN_NVAR) return;
    if (aao->eprobe_stage_ft.empty() || data->eprobe_stage_ft.empty()) return;

    const TH1D* aao_baseline=aao->eprobe_stage_ft[0].get();
    const TH1D* data_baseline=data->eprobe_stage_ft[0].get();
    if (!aao_baseline || !data_baseline) return;

    const double aao_base=std::max(1.0,aao_baseline->Integral());
    const double data_base=std::max(1.0,data_baseline->Integral());
    const double aao_high_base=std::max(1.0,hist_integral_range(aao_baseline,2.0,EPROBE_MAX+1e-6));

    std::vector<FTCutScanPoint> pts;
    std::array<std::vector<FTCutScanPoint>,CUTSCAN_NVAR> byvar;

    for (int iv=0;iv<CUTSCAN_NVAR;iv++) {
        for (int ip=0;ip<CUTSCAN_NPOINT;ip++) {
            const double cut=cutscan_value(iv,ip);

            auto hd=cutscan_project_eprobe(data->cutscan_ft[iv].get(),iv,cut,
                                           Form("scan_data_%d_%d",iv,ip));
            auto ha=cutscan_project_eprobe(aao->cutscan_ft[iv].get(),iv,cut,
                                           Form("scan_aao_%d_%d",iv,ip));
            auto hc=cutscan_project_eprobe(cls->cutscan_ft[iv].get(),iv,cut,
                                           Form("scan_cls_%d_%d",iv,ip));
            auto hv=cutscan_project_eprobe(dvc->cutscan_ft[iv].get(),iv,cut,
                                           Form("scan_dvc_%d_%d",iv,ip));
            if (!hd || !ha || !hc || !hv) continue;

            FTCutScanPoint q;
            q.iv=iv;
            q.cut=cut;
            q.y_data=hd->Integral();
            q.y_aao=R.nominal.aao*ha->Integral();
            q.y_clasdis=R.nominal.clasdis*hc->Integral();
            q.y_dvcs=R.nominal.dvcs*hv->Integral();

            const double fc=selected_clasdis_pi0_fraction(hc.get(),cls);
            const double total=q.y_aao+q.y_clasdis+q.y_dvcs;
            if (total>0)
                q.pi0_purity=(q.y_aao+fc*q.y_clasdis)/total;

            const double two=q.y_aao+q.y_dvcs;
            if (two>0) q.two_gen_purity=q.y_aao/two;

            q.data_retention=q.y_data/data_base;
            q.aao_retention=ha->Integral()/aao_base;
            q.aao_high_retention=
                hist_integral_range(ha.get(),2.0,EPROBE_MAX+1e-6)/aao_high_base;
            q.eprobe_distortion=eprobe_total_variation(ha.get(),aao_baseline);

            pts.push_back(q);
            byvar[iv].push_back(q);
        } // endfor
    } // endfor

    const char* row_titles[CUTSCAN_NVAR]={
        "M_{X}^{2}(ep#gamma)",
        "Trento coplanarity",
        "M_{X}^{2}(ep)",
        "angle(#gamma,X)"
    };

    TCanvas c("c_ft_cut_optimization","",2100,2200);
    c.Divide(4,CUTSCAN_NVAR,0.001,0.001);

    for (int iv=0;iv<CUTSCAN_NVAR;iv++) {
        if (byvar[iv].empty()) continue;

        for (int metric=0;metric<4;metric++) {
            c.cd(iv*4+metric+1);
            gPad->SetLeftMargin(metric==0?0.18:0.13);
            gPad->SetRightMargin(0.035);
            gPad->SetBottomMargin(iv==CUTSCAN_NVAR-1?0.18:0.12);
            gPad->SetTopMargin(iv==0?0.18:0.08);
            gPad->SetTicks(1,1);

            std::unique_ptr<TGraph> g(new TGraph());
            int ipt=0;
            for (const auto& q:byvar[iv]) {
                double y=0;
                if (metric==0) y=q.two_gen_purity;
                if (metric==1) y=q.aao_retention;
                if (metric==2) y=q.aao_high_retention;
                if (metric==3) y=q.eprobe_distortion;
                g->SetPoint(ipt++,q.cut,y);
            } // endfor

            g->SetLineWidth(3);
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.65);
            g->SetTitle("");
            g->GetXaxis()->SetTitle(iv==CUTSCAN_NVAR-1?CUTSCAN_XTITLE[iv]:"");

            if (metric==0) {
                g->GetYaxis()->SetTitle("AAO/(AAO+DVCS)");
                g->SetMinimum(0); g->SetMaximum(1.02);
            } else if (metric==1) {
                g->GetYaxis()->SetTitle("AAO retention");
                g->SetMinimum(0); g->SetMaximum(1.02);
            } else if (metric==2) {
                g->GetYaxis()->SetTitle("AAO E_{probe}#geq2 retention");
                g->SetMinimum(0); g->SetMaximum(1.02);
            } else {
                g->GetYaxis()->SetTitle("E_{probe} shape distortion");
                g->SetMinimum(0);
                double ymax=0.05;
                for (const auto& q:byvar[iv]) ymax=std::max(ymax,q.eprobe_distortion);
                g->SetMaximum(std::min(1.0,1.20*ymax));
            } // endif

            g->GetXaxis()->SetTitleSize(0.050);
            g->GetYaxis()->SetTitleSize(0.047);
            g->GetXaxis()->SetLabelSize(0.040);
            g->GetYaxis()->SetLabelSize(0.040);
            g->GetYaxis()->SetTitleOffset(metric==0?1.65:1.35);
            g->Draw("ALP");

            const double xnom=cutscan_nominal_value(iv);
            const double ymin=g->GetYaxis()->GetXmin();
            const double ymax=g->GetYaxis()->GetXmax();
            TLine* l=new TLine(xnom,ymin,xnom,ymax);
            l->SetLineStyle(2);
            l->SetLineWidth(2);
            l->Draw();

            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            if (metric==0) {
                tx.SetTextSize(0.048);
                tx.DrawLatex(0.18,0.91,Form("(%c) %s",'a'+iv,row_titles[iv]));
            } // endif
            if (iv==0) {
                tx.SetTextSize(0.043);
                const char* hdr=
                    metric==0?"#pi^{0} purity proxy":
                    metric==1?"all-E AAO retention":
                    metric==2?"high-E AAO retention":
                              "energy sculpting";
                tx.DrawLatex(0.18,0.91,hdr);
            } // endif
        } // endfor
    } // endfor

    c.SaveAs((dir+"/FT_cut_optimization.png").c_str());

    // Append the full scan to the one existing FT exclusivity summary file.
    std::ofstream csv(dir+"/summary.csv",std::ios::app);
    csv << "\n# FT one-variable N-1 cut scan\n";
    csv << "scan_variable,cut_value,pi0_purity_nominal,two_generator_AAO_fraction,"
           "data_retention_from_detector_baseline,AAO_retention_from_detector_baseline,"
           "AAO_highE_retention_from_detector_baseline,Eprobe_shape_distortion,"
           "normalized_AAO_yield,normalized_CLASDIS_yield,normalized_DVCS_yield\n";
    csv << std::setprecision(10);
    for (const auto& q:pts) {
        csv << CUTSCAN_KEY[q.iv] << "," << q.cut << ","
            << q.pi0_purity << "," << q.two_gen_purity << ","
            << q.data_retention << "," << q.aao_retention << ","
            << q.aao_high_retention << "," << q.eprobe_distortion << ","
            << q.y_aao << "," << q.y_clasdis << "," << q.y_dvcs << "\n";
    } // endfor
    csv.close();
}

struct FTFitQualityPoint {
    double copl=0;
    double aao=0;
    double clasdis=0;
    double dvcs=0;
    double low_chi2ndf=0;
    double high_chi2ndf=0;
    double data_low=0;
    double data_high=0;
    double aao_highE_probe=0;
};

void draw_ft_fitquality_scan(const std::vector<std::unique_ptr<ValComponent>>& vv,
                             const std::string& dir) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;

    // Use the same low-E simultaneous observable set as the nominal fit:
    // Mx2(ep), Mx2(epgamma), tag E, Mx2(egamma).  angle(gamma,X) remains validation-only.
    const std::vector<int> low_obs={
        NORM_MX2_EP,
        NORM_MX2_EPG,
        NORM_EGAMMA,
        NORM_MX2_EG
    };
    const std::vector<int> high_obs={
        NORM_MX2_EPG,
        NORM_ANGLE_GX,
        NORM_MX2_EG,
        NORM_DPHI_TRENTO_SHIFT180,
        NORM_DELTA_T_PG
    };

    std::vector<FTFitQualityPoint> pts;

    for (int ic=0;ic<FT_COPL_SCAN_N;ic++) {
        FTFitQualityPoint q;
        q.copl=FT_COPL_SCAN_VALUES[ic];

        // Build temporary ValComponent-like normalization hist collections by
        // replacing the FT low/high 1D histograms with projections from this scan.
        std::vector<std::unique_ptr<TH1D>> dl,al,cl,vl;
        std::vector<std::unique_ptr<TH1D>> dh,ah,ch,vh;
        dl.resize(NORM_NOBS); al.resize(NORM_NOBS); cl.resize(NORM_NOBS); vl.resize(NORM_NOBS);
        dh.resize(NORM_NOBS); ah.resize(NORM_NOBS); ch.resize(NORM_NOBS); vh.resize(NORM_NOBS);

        for (int io=0;io<NORM_NOBS;io++) {
            dl[io]=project_fitqual_observable(data,false,ic,io,Form("dl_%d_%d",ic,io));
            al[io]=project_fitqual_observable(aao,false,ic,io,Form("al_%d_%d",ic,io));
            cl[io]=project_fitqual_observable(cls,false,ic,io,Form("cl_%d_%d",ic,io));
            vl[io]=project_fitqual_observable(dvc,false,ic,io,Form("vl_%d_%d",ic,io));

            dh[io]=project_fitqual_observable(data,true,ic,io,Form("dh_%d_%d",ic,io));
            ah[io]=project_fitqual_observable(aao,true,ic,io,Form("ah_%d_%d",ic,io));
            ch[io]=project_fitqual_observable(cls,true,ic,io,Form("ch_%d_%d",ic,io));
            vh[io]=project_fitqual_observable(dvc,true,ic,io,Form("vh_%d_%d",ic,io));
        } // endfor

        // Low-E simultaneous AAO+CLASDIS fit with common coefficients.
        // Start from a small grid and then coordinate-refine.  This mirrors the
        // nominal philosophy while remaining self-contained for the scan.
        double bestA=0,bestC=0,bestChi=1e300;
        int lowNdf=0;

        auto eval_low=[&](double A,double C)->std::pair<double,int> {
            double chi=0; int nbin=0;
            for (int io:low_obs) {
                if (!dl[io] || !al[io] || !cl[io]) continue;
                const int nb=dl[io]->GetNbinsX();
                for (int ib=1;ib<=nb;ib++) {
                    const double d=dl[io]->GetBinContent(ib);
                    const double ea=al[io]->GetBinError(ib);
                    const double ec=cl[io]->GetBinError(ib);
                    const double m=A*al[io]->GetBinContent(ib)+C*cl[io]->GetBinContent(ib);
                    const double ed=dl[io]->GetBinError(ib);
                    const double var=ed*ed+A*A*ea*ea+C*C*ec*ec;
                    if (!(var>0)) continue;
                    chi+=(d-m)*(d-m)/var;
                    nbin++;
                } // endfor
            } // endfor
            return {chi,nbin};
        };

        for (double A=0;A<=1.5;A+=0.05) {
            for (double C=0;C<=4.0;C+=0.10) {
                auto z=eval_low(A,C);
                if (z.first<bestChi) {
                    bestChi=z.first; bestA=A; bestC=C; lowNdf=z.second-2;
                } // endif
            } // endfor
        } // endfor
        for (double step: {0.02,0.01,0.005}) {
            double a0=bestA,c0=bestC;
            for (int ia=-5;ia<=5;ia++) {
                for (int icc=-5;icc<=5;icc++) {
                    double A=std::max(0.0,a0+ia*step);
                    double C=std::max(0.0,c0+icc*step*2.0);
                    auto z=eval_low(A,C);
                    if (z.first<bestChi) {
                        bestChi=z.first; bestA=A; bestC=C; lowNdf=z.second-2;
                    } // endif
                } // endfor
            } // endfor
        } // endfor

        q.aao=bestA;
        q.clasdis=bestC;
        q.low_chi2ndf=(lowNdf>0)?bestChi/lowNdf:0;

        // High-E DVCS-only scale with AAO/CLASDIS fixed to low-E values.
        double bestD=0,bestHChi=1e300;
        int highNdf=0;
        auto eval_high=[&](double D)->std::pair<double,int> {
            double chi=0; int nbin=0;
            for (int io:high_obs) {
                if (!dh[io] || !ah[io] || !ch[io] || !vh[io]) continue;
                const int nb=dh[io]->GetNbinsX();
                for (int ib=1;ib<=nb;ib++) {
                    const double d=dh[io]->GetBinContent(ib);
                    const double m=
                        bestA*ah[io]->GetBinContent(ib)+
                        bestC*ch[io]->GetBinContent(ib)+
                        D*vh[io]->GetBinContent(ib);
                    const double ed=dh[io]->GetBinError(ib);
                    const double ea=ah[io]->GetBinError(ib);
                    const double ec=ch[io]->GetBinError(ib);
                    const double ev=vh[io]->GetBinError(ib);
                    const double var=ed*ed+bestA*bestA*ea*ea+
                                     bestC*bestC*ec*ec+D*D*ev*ev;
                    if (!(var>0)) continue;
                    chi+=(d-m)*(d-m)/var;
                    nbin++;
                } // endfor
            } // endfor
            return {chi,nbin};
        };

        for (double D=0;D<=1.0;D+=0.01) {
            auto z=eval_high(D);
            if (z.first<bestHChi) {
                bestHChi=z.first; bestD=D; highNdf=z.second-1;
            } // endif
        } // endfor
        for (double step: {0.005,0.002}) {
            double d0=bestD;
            for (int id=-5;id<=5;id++) {
                double D=std::max(0.0,d0+id*step);
                auto z=eval_high(D);
                if (z.first<bestHChi) {
                    bestHChi=z.first; bestD=D; highNdf=z.second-1;
                } // endif
            } // endfor
        } // endfor

        q.dvcs=bestD;
        q.high_chi2ndf=(highNdf>0)?bestHChi/highNdf:0;

        if (dl[NORM_EGAMMA]) q.data_low=dl[NORM_EGAMMA]->Integral();
        if (dh[NORM_EGAMMA]) q.data_high=dh[NORM_EGAMMA]->Integral();

        // High-Eprobe AAO statistics after this coplanarity cut.
        auto ha=cutscan_project_eprobe(aao->cutscan_ft[1].get(),1,q.copl,
                                       Form("aao_probe_copl_%d",ic));
        if (ha) q.aao_highE_probe=hist_integral_range(ha.get(),2.0,EPROBE_MAX+1e-6);

        pts.push_back(q);
    } // endfor

    TCanvas c("c_ft_fitquality_scan","",1900,1050);
    c.Divide(3,2,0.001,0.001);

    auto make_graph=[&](int metric,const char* ytitle,double ymin,double ymax,int pad) {
        c.cd(pad);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.10);
        gPad->SetTicks(1,1);

        TGraph* g=new TGraph();
        int ip=0;
        for (const auto& q:pts) {
            double y=0;
            if (metric==0) y=q.low_chi2ndf;
            if (metric==1) y=q.high_chi2ndf;
            if (metric==2) y=q.aao;
            if (metric==3) y=q.clasdis;
            if (metric==4) y=q.dvcs;
            if (metric==5) y=q.aao_highE_probe;
            g->SetPoint(ip++,q.copl,y);
        } // endfor
        g->SetLineWidth(3);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.8);
        g->SetTitle(Form(";|#Delta#phi_{copl}| maximum (deg);%s",ytitle));
        g->GetXaxis()->SetTitleSize(0.050);
        g->GetYaxis()->SetTitleSize(0.050);
        g->GetYaxis()->SetTitleOffset(1.35);
        if (ymax>ymin) {
            g->SetMinimum(ymin);
            g->SetMaximum(ymax);
        } // endif
        g->Draw("ALP");

        const double ylo=g->GetYaxis()->GetXmin();
        const double yhi=g->GetYaxis()->GetXmax();
        TLine* l=new TLine(NORM_DPHI_TRENTO_MAX,ylo,NORM_DPHI_TRENTO_MAX,yhi);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw();
    };

    double maxLow=2.0,maxHigh=2.0,maxHighStat=1.0;
    for (const auto& q:pts) {
        maxLow=std::max(maxLow,q.low_chi2ndf*1.2);
        maxHigh=std::max(maxHigh,q.high_chi2ndf*1.2);
        maxHighStat=std::max(maxHighStat,q.aao_highE_probe*1.15);
    } // endfor

    make_graph(0,"Low-E simultaneous #chi^{2}/ndf",0,maxLow,1);
    make_graph(1,"High-E DVCS #chi^{2}/ndf",0,maxHigh,2);
    make_graph(2,"AAO normalization",0,0,3);
    make_graph(3,"CLASDIS normalization",0,0,4);
    make_graph(4,"DVCS normalization",0,0,5);
    make_graph(5,"AAO E_{probe}#geq2 raw events",0,maxHighStat,6);

    c.SaveAs((dir+"/FT_coplanarity_fit_quality_scan.png").c_str());

    std::ofstream csv(dir+"/summary.csv",std::ios::app);
    csv << "\n# FT coplanarity normalization-fit quality scan\n";
    csv << "coplanarity_absmax_deg,AAO_scale,CLASDIS_scale,DVCS_scale,"
           "lowE_chi2_ndf,highE_chi2_ndf,lowE_data_events,highE_data_events,"
           "AAO_highEprobe_raw_events\n";
    csv << std::setprecision(10);
    for (const auto& q:pts)
        csv << q.copl << "," << q.aao << "," << q.clasdis << "," << q.dvcs << ","
            << q.low_chi2ndf << "," << q.high_chi2ndf << ","
            << q.data_low << "," << q.data_high << "," << q.aao_highE_probe << "\n";
    csv.close();
}

void draw_exclusivity_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                              const std::string& dir,
                              bool ft_probe) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;

    struct P { int obs,stage; double c1,c2; bool two; const char* label; };
    const P pp[]={
        {NORM_MX2_EP,0,NORM_MX2_EP_MIN,NORM_MX2_EP_MAX,true,"(a) M_{X}^{2}(ep)"},
        {NORM_MX2_EG,1,NORM_MX2_EG_MIN,0,false,"(b) M_{X}^{2}(e#gamma)"},
        {NORM_DPHI_TRENTO_SHIFT180,2,-NORM_DPHI_TRENTO_MAX,NORM_DPHI_TRENTO_MAX,true,"(c) Trento coplanarity"},
        {NORM_ANGLE_GX,3,NORM_ANGLE_GX_MAX,0,false,"(d) angle(#gamma,X)"}
    };

    auto geth=[&](const ValComponent* v,int stage,int io)->const TH1D* {
        if (stage==0) return norm_pre_for(v,ft_probe)[io].get();
        if (stage==1) return norm_after_ep_for(v,ft_probe)[io].get();
        if (stage==2) return norm_after_eg_for(v,ft_probe)[io].get();
        return norm_after_dphi_for(v,ft_probe)[io].get();
    };

    const ValComponent* ss[]={data,aao,cls,dvc};
    const int col[]={kBlack,kRed+1,kOrange+7,kGreen+2};
    const char* lab[]={"Data","AAOgen","CLASDIS","DVCSgen"};
    std::vector<std::unique_ptr<TH1D>> keep;

    TCanvas c(ft_probe?"c_exclusivity_summary_ft":"c_exclusivity_summary_fd","",1500,1100);
    c.Divide(2,2);

    for (int ip=0;ip<4;ip++) {
        c.cd(ip+1);
        gPad->SetLeftMargin(0.13);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.14);

        std::vector<TH1D*> hh;
        double ymax=0;
        for (int is=0;is<4;is++) {
            const TH1D* h0=geth(ss[is],pp[ip].stage,pp[ip].obs);
            std::unique_ptr<TH1D> h((TH1D*)h0->Clone(Form("exc_%d_%d",ip,is)));
            h->SetDirectory(nullptr);
            const double q=h->Integral();
            if (q>0) h->Scale(1.0/q);
            h->SetStats(0);
            if (is==0) {
                h->SetMarkerStyle(20);
                h->SetMarkerSize(0.45);
                h->SetLineColor(kBlack);
            } else {
                h->SetLineColor(col[is]);
                h->SetLineWidth(2);
            } // endif
            ymax=std::max(ymax,h->GetMaximum());
            hh.push_back(h.get());
            keep.push_back(std::move(h));
        } // endfor

        hh[0]->SetTitle("");
        hh[0]->GetYaxis()->SetTitle("Unit-area candidates");

        // Display ranges only; event selection and normalization fits are unchanged.
        if (pp[ip].obs==NORM_ANGLE_GX) {
            for (auto* h:hh) if (h) h->GetXaxis()->SetRangeUser(0.0,20.0);
            hh[0]->GetXaxis()->SetTitle("angle(#gamma,X) (deg)");
        } else if (pp[ip].obs==NORM_DPHI_TRENTO_SHIFT180) {
            hh[0]->GetXaxis()->SetTitle("#Delta#phi_{copl} (deg)");
        }

        hh[0]->SetMaximum(1.27*ymax);
        hh[0]->Draw("E1");
        for (int is=1;is<4;is++) hh[is]->Draw("HIST SAME");
        hh[0]->Draw("E1 SAME");

        const double yy=1.18*ymax;
        TLine* l1=new TLine(pp[ip].c1,0,pp[ip].c1,yy);
        l1->SetLineStyle(2); l1->SetLineWidth(2); l1->Draw();
        if (pp[ip].two) {
            TLine* l2=new TLine(pp[ip].c2,0,pp[ip].c2,yy);
            l2->SetLineStyle(2); l2->SetLineWidth(2); l2->Draw();
        } // endif

        TLatex tx;
        tx.SetNDC();
        tx.SetTextFont(42);
        tx.SetTextSize(0.050);
        tx.DrawLatex(0.15,0.94,pp[ip].label);

        if (ip==0) {
            TLegend* leg=new TLegend(0.64,0.54,0.94,0.79);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            for (int is=0;is<4;is++)
                leg->AddEntry(hh[is],lab[is],is==0?"lep":"l");
            leg->Draw();
        } // endif
    } // endfor

    c.SaveAs((dir+"/exclusivity_summary.png").c_str());

    std::ofstream csv(dir+"/summary.csv");
    csv << "sample,baseline,mx2_ep,mx2_eg,coplanarity,angle_gX,all_cuts\n";
    for (const auto& v:vv) {
        csv << v->name;
        const auto& cf=norm_cutflow_for(v.get(),ft_probe);
        for (int i=0;i<6;i++) csv << "," << cf[i];
        csv << "\n";
    } // endfor
    csv.close();

    draw_exclusivity_correlations(vv,dir,ft_probe);
    draw_eprobe_cut_survival(vv,dir,ft_probe);
}


std::string pretty_norm_observable(const std::string& key) {
    if (key=="Mx2_ep") return "M_{X}^{2}(ep)";
    if (key=="Mx2_epg") return "M_{X}^{2}(ep#gamma)";
    if (key=="angle_gX") return "angle(#gamma,X)";
    if (key=="Egamma") return "E_{#gamma}";
    if (key=="Mx2_eg") return "M_{X}^{2}(e#gamma)";
    if (key=="dphi_trento_shift180") return "#Delta#phi_{copl}";
    if (key=="delta_t_pg") return "#Delta t(p,#gamma)";
    return key;
}

void draw_norm_panel(TH1D* frame,const TH1D* hd,const TH1D* ha,const TH1D* hc,const TH1D* hv,
                     double A,double B,double C,const NormFitPoint& q,const char* panel) {
    std::unique_ptr<TH1D> d((TH1D*)hd->Clone(Form("d_%s",panel))); d->SetDirectory(nullptr);
    auto a=morph_norm_hist(ha,q.morph_shift,q.morph_sigma,Form("a_%s",panel)); a->Scale(A);
    auto c=morph_norm_hist(hc,q.morph_shift,q.morph_sigma,Form("c_%s",panel)); c->Scale(B);
    auto v=morph_norm_hist(hv,q.morph_shift,q.morph_sigma,Form("v_%s",panel)); v->Scale(C);
    std::unique_ptr<TH1D> t((TH1D*)a->Clone(Form("t_%s",panel))); t->Add(c.get()); t->Add(v.get());

    d->SetStats(0);
    d->SetMarkerStyle(20);
    d->SetMarkerSize(0.45);
    d->SetLineColor(kBlack);
    style_norm_component(a.get(),kRed+1);
    style_norm_component(c.get(),kOrange+7);
    style_norm_component(v.get(),kGreen+2);
    style_norm_component(t.get(),kBlue+1,3);
    a->SetLineStyle(2);
    c->SetLineStyle(3);
    v->SetLineStyle(4);
    t->SetLineStyle(1);

    d->SetTitle("");
    d->GetYaxis()->SetTitleOffset(1.55);
    d->GetXaxis()->SetTitleOffset(1.12);

    // Compact display ranges only.  The template fits themselves use the
    // original full histogram ranges defined above.
    if (q.observable=="Mx2_ep") {
        d->GetXaxis()->SetRangeUser(-0.15,0.15);
        d->GetXaxis()->SetTitle("M_{X}^{2}(ep) (GeV^{2})");
    } else if (q.observable=="Mx2_epg") {
        d->GetXaxis()->SetRangeUser(-0.15,0.15);
        d->GetXaxis()->SetTitle("M_{X}^{2}(ep#gamma) (GeV^{2})");
    } else if (q.observable=="angle_gX") {
        d->GetXaxis()->SetRangeUser(0.0,20.0);
        d->GetXaxis()->SetTitle("angle(#gamma,X) (deg)");
    } else if (q.observable=="Egamma") {
        d->GetXaxis()->SetRangeUser(0.4,2.05);
        d->GetXaxis()->SetTitle("Tag-photon energy E_{#gamma} (GeV)");
    } else if (q.observable=="Mx2_eg") {
        d->GetXaxis()->SetRangeUser(1.0,5.5);
        d->GetXaxis()->SetTitle("M_{X}^{2}(e#gamma) (GeV^{2})");
    } else if (q.observable=="dphi_trento_shift180") {
        d->GetXaxis()->SetRangeUser(-20.0,20.0);
        d->GetXaxis()->SetTitle("#Delta#phi_{copl} (deg)");
    } else if (q.observable=="delta_t_pg") {
        d->GetXaxis()->SetRangeUser(-0.5,1.5);
        d->GetXaxis()->SetTitle("#Delta t=t_{p}-t_{#gamma} (GeV^{2})");
    }

    d->SetMaximum(1.30*std::max(d->GetMaximum(),t->GetMaximum()));
    TH1* dcopy=d->DrawCopy("E1");
    TH1* tcopy=t->DrawCopy("HIST SAME");
    TH1* acopy=a->DrawCopy("HIST SAME");
    TH1* ccopy=c->DrawCopy("HIST SAME");
    TH1* vcopy=v->DrawCopy("HIST SAME");
    d->DrawCopy("E1 SAME");

    if (std::string(panel).rfind("(a)",0)==0) {
        TLegend* leg=new TLegend(0.52,0.55,0.93,0.79);
        leg->SetNColumns(2);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.030);
        leg->AddEntry(dcopy,"Data","lep");
        leg->AddEntry(acopy,"AAOgen","l");
        leg->AddEntry(ccopy,"CLASDIS","l");
        leg->AddEntry(vcopy,"DVCSgen","l");
        leg->AddEntry(tcopy,"Total MC","l");
        leg->Draw();
    }

    TLatex tx;
    tx.SetNDC();
    tx.SetTextFont(42);
    tx.SetTextSize(0.045);
    tx.DrawLatex(0.15,0.92,panel);
    tx.SetTextSize(0.032);
    if (C==0.0)
        tx.DrawLatex(0.54,0.84,Form("AAO %.3f, CLASDIS %.3f",A,B));
    else
        tx.DrawLatex(0.64,0.84,Form("DVCS %.3f",C));
}

void draw_normalization_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                const NormDerivation& R,
                                const std::string& dir,
                                bool ft_probe) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;

    // Energy-region overview.
    {
        const int io=NORM_EGAMMA;
        const auto& df=norm_full_for(data,ft_probe);
        const auto& af=norm_full_for(aao,ft_probe);
        const auto& cf=norm_full_for(cls,ft_probe);
        const auto& vf=norm_full_for(dvc,ft_probe);
        std::unique_ptr<TH1D> d((TH1D*)df[io]->Clone("normE_d"));
        std::unique_ptr<TH1D> a((TH1D*)af[io]->Clone("normE_a"));
        std::unique_ptr<TH1D> c((TH1D*)cf[io]->Clone("normE_c"));
        std::unique_ptr<TH1D> v((TH1D*)vf[io]->Clone("normE_v"));
        d->SetDirectory(nullptr); a->SetDirectory(nullptr); c->SetDirectory(nullptr); v->SetDirectory(nullptr);
        a->Scale(R.nominal.aao); c->Scale(R.nominal.clasdis); v->Scale(R.nominal.dvcs);
        std::unique_ptr<TH1D> t((TH1D*)a->Clone("normE_t")); t->Add(c.get()); t->Add(v.get());
        d->SetStats(0); d->SetMarkerStyle(20); d->SetMarkerSize(0.55); d->SetLineColor(kBlack);
        style_norm_component(a.get(),kRed+1); style_norm_component(c.get(),kOrange+7);
        style_norm_component(v.get(),kGreen+2); style_norm_component(t.get(),kBlue+1,3);

        TCanvas ce("c_norm_energy","",1150,760);
        ce.SetLeftMargin(0.12); ce.SetRightMargin(0.04); ce.SetTopMargin(0.16);
        d->SetMaximum(1.28*std::max(d->GetMaximum(),t->GetMaximum()));
        d->GetXaxis()->SetTitle("Tag-photon energy E_{#gamma} (GeV)");
        d->GetYaxis()->SetTitle("Candidates");
        d->Draw("E1"); a->Draw("HIST SAME"); c->Draw("HIST SAME"); v->Draw("HIST SAME");
        t->Draw("HIST SAME"); d->Draw("E1 SAME");
        TLine l2(2.0,0,2.0,d->GetMaximum()); l2.SetLineStyle(2); l2.Draw();
        TLine l3(3.0,0,3.0,d->GetMaximum()); l3.SetLineStyle(2); l3.Draw();
        TLegend leg(0.64,0.57,0.94,0.82); leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(d.get(),"Data","lep"); leg.AddEntry(a.get(),"AAOgen","l");
        leg.AddEntry(c.get(),"CLASDIS","l"); leg.AddEntry(v.get(),"DVCSgen","l");
        leg.AddEntry(t.get(),"Total MC","l"); leg.Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.040);
        tx.DrawLatex(0.14,0.945,Form("%s-probe normalization: low tag-E fits AAOgen+CLASDIS; high tag-E fits DVCSgen",
                                      ft_probe?"FT":"FD"));
        ce.SaveAs((dir+"/energy_regions.png").c_str());
    }

    // The low- and high-energy canvases intentionally use different
    // observable sets because this follows Valerii's two-stage normalization:
    // low E constrains AAOgen+CLASDIS using Mx2(ep), Mx2(epgamma),
    // angle(gamma,X), tag-photon energy, and Mx2(egamma);
    // high E fixes those two normalizations and constrains DVCSgen using
    // Mx2(epgamma), angle(gamma,X), Mx2(egamma), Delta-phi, and Delta-t.
    // Common observables therefore appear in both stages; energy is low-E-only,
    // while Delta-phi and Delta-t are high-E-only.
    // Low-E fit canvas.
    {
        TCanvas c("c_low_norm","",1500,980);
        c.Divide(3,2);
        int pad=0;
        for (const auto& q:R.low_points) {
            int io=-1;
            for (int j=0;j<NORM_NOBS;j++) if (q.observable==NORM_OBS[j].key) io=j;
            if (io<0 || pad>=6) continue;
            c.cd(++pad);
            gPad->SetLeftMargin(0.20); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.21); gPad->SetTopMargin(0.12);
            draw_norm_panel(nullptr,norm_low_for(data,ft_probe)[io].get(),norm_low_for(aao,ft_probe)[io].get(),
                            norm_low_for(cls,ft_probe)[io].get(),norm_low_for(dvc,ft_probe)[io].get(),
                            R.nominal.aao,R.nominal.clasdis,R.nominal.dvcs,q,
                            Form("(%c) %s%s",'a'+pad-1,pretty_norm_observable(q.observable).c_str(),
                                 q.observable=="angle_gX"?" [validation]":""));
            if (pad==1) {
                TLegend* leg=new TLegend(0.57,0.58,0.94,0.86);
                leg->SetBorderSize(0); leg->SetFillStyle(0);
                leg->AddEntry((TObject*)nullptr,Form("common AAO = %.3f",R.nominal.aao),"");
                leg->AddEntry((TObject*)nullptr,Form("common CLASDIS = %.3f",R.nominal.clasdis),"");
                leg->Draw();
            } // endif
        } // endfor
        if (pad<6) {
            c.cd(6);
            gPad->SetLeftMargin(0.10);
            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            tx.SetTextSize(0.060);
            tx.DrawLatex(0.12,0.82,"Simultaneous low-E normalization");
            tx.SetTextSize(0.052);
            tx.DrawLatex(0.12,0.66,Form("AAOgen = %.3f",R.nominal.aao));
            tx.DrawLatex(0.12,0.54,Form("CLASDIS = %.3f",R.nominal.clasdis));
            if (R.no_clasdis_valid) {
                tx.SetTextSize(0.044);
                tx.DrawLatex(0.12,0.45,Form("FT AAO+DVCS-only: AAO = %.3f",R.no_clasdis.aao));
            } // endif
            tx.SetTextSize(0.040);
            tx.DrawLatex(0.12,0.38,"One common AAOgen/CLASDIS pair fits");
            tx.DrawLatex(0.12,0.32,"M_{X}^{2}(ep), M_{X}^{2}(ep#gamma), E_{#gamma}, M_{X}^{2}(e#gamma).");
            tx.DrawLatex(0.12,0.26,"angle(#gamma,X) is validation-only; DVCS overlay fixed from high E.");

            TLegend leg(0.10,0.05,0.92,0.25);
            leg.SetNColumns(2);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            TGraph gData,gAAO,gCLS,gDVCS,gTot;
            gData.SetMarkerStyle(20); gData.SetMarkerColor(kBlack);
            gAAO.SetLineColor(kRed+1); gAAO.SetLineWidth(2);
            gCLS.SetLineColor(kOrange+7); gCLS.SetLineWidth(2);
            gDVCS.SetLineColor(kGreen+2); gDVCS.SetLineWidth(2);
            gTot.SetLineColor(kBlue+1); gTot.SetLineWidth(3);
            leg.AddEntry(&gData,"Data","p");
            leg.AddEntry(&gAAO,"AAOgen","l");
            leg.AddEntry(&gCLS,"CLASDIS","l");
            leg.AddEntry(&gDVCS,"DVCSgen","l");
            leg.AddEntry(&gTot,"Total MC","l");
            leg.Draw();
        } // endif
        c.SaveAs((dir+"/lowE_fits.png").c_str());
    }

    // High-E fit canvas.
    {
        TCanvas c("c_high_norm","",1500,980);
        c.Divide(3,2);
        int pad=0;
        for (const auto& q:R.high_points) {
            int io=-1;
            for (int j=0;j<NORM_NOBS;j++) if (q.observable==NORM_OBS[j].key) io=j;
            if (io<0 || pad>=6) continue;
            c.cd(++pad);
            gPad->SetLeftMargin(0.20); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.21); gPad->SetTopMargin(0.12);
            draw_norm_panel(nullptr,norm_high_for(data,ft_probe)[io].get(),norm_high_for(aao,ft_probe)[io].get(),
                            norm_high_for(cls,ft_probe)[io].get(),norm_high_for(dvc,ft_probe)[io].get(),
                            R.nominal.aao,R.nominal.clasdis,q.dvcs,q,
                            Form("(%c) %s",'a'+pad-1,pretty_norm_observable(q.observable).c_str()));
        } // endfor
        if (pad<6) {
            c.cd(6);
            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            tx.SetTextSize(0.060);
            tx.DrawLatex(0.12,0.82,"Mean high-E normalization");
            tx.SetTextSize(0.052);
            tx.DrawLatex(0.12,0.64,Form("DVCSgen = %.3f",R.nominal.dvcs));
            tx.SetTextSize(0.040);
            tx.DrawLatex(0.12,0.46,"AAOgen + CLASDIS fixed to");
            tx.DrawLatex(0.12,0.40,"their low-E mean values.");

            TLegend leg(0.10,0.08,0.92,0.30);
            leg.SetNColumns(2);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            TGraph gData,gAAO,gCLS,gDVCS,gTot;
            gData.SetMarkerStyle(20); gData.SetMarkerColor(kBlack);
            gAAO.SetLineColor(kRed+1); gAAO.SetLineWidth(2);
            gCLS.SetLineColor(kOrange+7); gCLS.SetLineWidth(2);
            gDVCS.SetLineColor(kGreen+2); gDVCS.SetLineWidth(2);
            gTot.SetLineColor(kBlue+1); gTot.SetLineWidth(3);
            leg.AddEntry(&gData,"Data","p");
            leg.AddEntry(&gAAO,"AAOgen","l");
            leg.AddEntry(&gCLS,"CLASDIS","l");
            leg.AddEntry(&gDVCS,"DVCSgen","l");
            leg.AddEntry(&gTot,"Total MC","l");
            leg.Draw();
        } // endif
        c.SaveAs((dir+"/highE_fits.png").c_str());
    }

    std::ofstream csv(dir+"/summary.csv");
    csv << "stage,observable,aao,clasdis,dvcs,chi2_ndf,morph_shift,morph_sigma\n";

    csv << "lowE_global,simultaneous_without_angle_gX,"
        <<R.low_global.aao<<","<<R.low_global.clasdis<<",0,"
        <<(R.low_global.ndf?R.low_global.chi2/R.low_global.ndf:0)<<",0,0\n";

    for (const auto& q:R.low_points)
        csv << "lowE_diagnostic,"<<q.observable<<","<<q.aao<<","<<q.clasdis<<",0,"
            <<(q.ndf?q.chi2/q.ndf:0)<<","<<q.morph_shift<<","<<q.morph_sigma<<"\n";

    for (const auto& q:R.high_points)
        csv << "highE,"<<q.observable<<","<<R.nominal.aao<<","<<R.nominal.clasdis<<","<<q.dvcs<<","
            <<(q.ndf?q.chi2/q.ndf:0)<<","<<q.morph_shift<<","<<q.morph_sigma<<"\n";

    std::vector<double> av,bv,cv;
    for (const auto& q:R.low_points) {
        if (q.observable=="angle_gX") continue;
        av.push_back(q.aao);
        bv.push_back(q.clasdis);
    } // endfor
    for (const auto& q:R.high_points) cv.push_back(q.dvcs);

    csv << "nominal_simultaneous,ALL,"
        <<R.nominal.aao<<","<<R.nominal.clasdis<<","<<R.nominal.dvcs<<",0,0,0\n";
    csv << "diagnostic_stddev_without_angle_gX,ALL,"
        <<rms_spread(av,R.nominal.aao)<<","<<rms_spread(bv,R.nominal.clasdis)<<","
        <<rms_spread(cv,R.nominal.dvcs)<<",0,0,0\n";

    if (R.no_clasdis_valid)
        csv << "FT_alternate_no_CLASDIS,ALL,"
            <<R.no_clasdis.aao<<",0,"<<R.no_clasdis.dvcs<<",0,0,0\n";

}

void draw_pi0_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                      const NormDerivation& Rfd,
                      const NormDerivation& Rft,
                      const std::string& dir) {
    const ValComponent* cls=find_val_component(vv,"clasdis");

    double fc_global=0.9877;
    long long tr=0,tp=0;
    if (cls) {
        for (int ir=0;ir<CR_N;ir++) {
            tr+=cls->coarse[ir].truth_rows;
            tp+=cls->coarse[ir].truth_pi0_rows;
        } // endfor
    } // endif
    if (tr>0) fc_global=double(tp)/double(tr);

    TCanvas c("c_pi0_summary","",1050,720);
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.22);
    c.SetTopMargin(0.12);

    TH1D h("h_pi0_regions",
           ";Detector / probe-energy region;#pi^{0}-bearing fraction of selected ep#gammaX",
           CR_N,0,CR_N);
    h.SetStats(0);
    h.SetMinimum(0.0);
    h.SetMaximum(1.05);
    h.SetMarkerStyle(20);
    h.SetMarkerSize(1.45);

    for (int ir=0;ir<CR_N;ir++) {
        h.GetXaxis()->SetBinLabel(ir+1,CR_LABEL[ir]);
        const auto q=coarse_composition(vv,(ir<2?Rfd.nominal:Rft.nominal),ir);
        if (q.f_pi0>=0) h.SetBinContent(ir+1,q.f_pi0);
    } // endfor

    h.Draw("P");

    TGraphErrors gAlt;
    if (Rft.no_clasdis_valid) {
        int ip=0;
        for (int ir=CR_FT_LOW;ir<=CR_FT_HIGH;ir++) {
            const auto qa=coarse_composition(vv,Rft.no_clasdis,ir);
            if (qa.f_pi0<0) continue;
            gAlt.SetPoint(ip,ir+0.5,qa.f_pi0);
            gAlt.SetPointError(ip,0,0);
            ip++;
        } // endfor
        gAlt.SetMarkerStyle(24);
        gAlt.SetMarkerSize(1.5);
        gAlt.SetLineWidth(2);
        gAlt.Draw("P SAME");

        TLegend legAlt(0.55,0.77,0.94,0.89);
        legAlt.SetBorderSize(0);
        legAlt.SetFillStyle(0);
        legAlt.SetTextSize(0.032);
        legAlt.AddEntry(&h,"nominal AAO + CLASDIS + DVCS","p");
        legAlt.AddEntry(&gAlt,"FT sensitivity: AAO + DVCS only","p");
        legAlt.Draw();
    } // endif

    TLatex tx;
    tx.SetNDC();
    tx.SetTextFont(42);
    tx.SetTextSize(0.040);
    tx.DrawLatex(0.15,0.93,Form("AAOgen: 100%% #pi^{0};  CLASDIS truth: %.1f%% #pi^{0};  DVCSgen: 0%% #pi^{0}",
                                100.0*fc_global));

    c.SaveAs((dir+"/pi0_fraction.png").c_str());

    {
        TH1D hA("h_comp_aao",";Detector / probe-energy region;Normalized component fraction",CR_N,0,CR_N);
        TH1D hC("h_comp_cls","",CR_N,0,CR_N);
        TH1D hD("h_comp_dvcs","",CR_N,0,CR_N);
        TH1D hP("h_comp_pi0","",CR_N,0,CR_N);

        hA.SetStats(0); hA.SetMinimum(0.0); hA.SetMaximum(1.05);
        hA.SetLineColor(kRed+1); hA.SetLineStyle(2); hA.SetLineWidth(3);
        hC.SetLineColor(kOrange+7); hC.SetLineStyle(3); hC.SetLineWidth(3);
        hD.SetLineColor(kGreen+2); hD.SetLineStyle(4); hD.SetLineWidth(3);
        hP.SetMarkerColor(kBlack); hP.SetMarkerStyle(20); hP.SetMarkerSize(1.2);

        for (int ir=0;ir<CR_N;ir++) {
            hA.GetXaxis()->SetBinLabel(ir+1,CR_LABEL[ir]);
            const auto q=coarse_composition(vv,(ir<2?Rfd.nominal:Rft.nominal),ir);
            const double s=q.ya+q.yc+q.yd;
            if (s>0) {
                hA.SetBinContent(ir+1,q.ya/s);
                hC.SetBinContent(ir+1,q.yc/s);
                hD.SetBinContent(ir+1,q.yd/s);
            }
            if (q.f_pi0>=0) hP.SetBinContent(ir+1,q.f_pi0);
        }

        TCanvas cc("c_pi0_components","",1150,760);
        cc.SetLeftMargin(0.13); cc.SetRightMargin(0.04);
        cc.SetBottomMargin(0.25); cc.SetTopMargin(0.12);
        hA.Draw("HIST"); hC.Draw("HIST SAME"); hD.Draw("HIST SAME"); hP.Draw("P SAME");

        TGraphErrors gAltComp;
        if (Rft.no_clasdis_valid) {
            int ip=0;
            for (int ir=CR_FT_LOW;ir<=CR_FT_HIGH;ir++) {
                const auto qa=coarse_composition(vv,Rft.no_clasdis,ir);
                if (qa.f_pi0<0) continue;
                gAltComp.SetPoint(ip,ir+0.5,qa.f_pi0);
                gAltComp.SetPointError(ip,0,0);
                ip++;
            } // endfor
            gAltComp.SetMarkerStyle(24);
            gAltComp.SetMarkerSize(1.35);
            gAltComp.Draw("P SAME");
        } // endif

        TLegend leg(0.58,0.60,0.93,0.86);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(&hA,"AAOgen","l");
        leg.AddEntry(&hC,"CLASDIS","l");
        leg.AddEntry(&hD,"DVCSgen","l");
        leg.AddEntry(&hP,"#pi^{0}-bearing fraction","p");
        if (Rft.no_clasdis_valid)
            leg.AddEntry(&gAltComp,"FT #pi^{0} fraction, AAO + DVCS only","p");
        leg.Draw();

        TLatex txc; txc.SetNDC(); txc.SetTextFont(42); txc.SetTextSize(0.040);
        txc.DrawLatex(0.15,0.93,"Normalized MC composition behind the #pi^{0}-fraction estimate");
        cc.SaveAs((dir+"/composition_by_region.png").c_str());
    }

    std::ofstream csv(dir+"/summary.csv");
    csv << "region,aao_yield,clasdis_yield,dvcs_yield,clasdis_pi0_fraction,pi0_fraction,"
           "alt_no_clasdis_aao_yield,alt_no_clasdis_dvcs_yield,alt_no_clasdis_pi0_fraction\n";
    for (int ir=0;ir<CR_N;ir++) {
        const auto q=coarse_composition(vv,(ir<2?Rfd.nominal:Rft.nominal),ir);
        double aya=0,ayd=0,af=-1;
        if (ir>=CR_FT_LOW && Rft.no_clasdis_valid) {
            const auto qa=coarse_composition(vv,Rft.no_clasdis,ir);
            aya=qa.ya;
            ayd=qa.yd;
            af=qa.f_pi0;
        } // endif
        csv << CR_KEY[ir]<<","<<q.ya<<","<<q.yc<<","<<q.yd<<","
            <<q.f_clasdis_pi0<<","<<q.f_pi0<<","
            <<aya<<","<<ayd<<","<<af<<"\n";
    } // endfor
}


struct IntegratedEfficiencyResult {
    bool valid=false;
    bool data_fit_valid=false,mc_fit_valid=false;
    std::string data_fit_reason="not attempted";
    std::string mc_fit_reason="not attempted";
    double f_pi0=0;
    double data_denom=0,mc_pi0_denom=0;
    double data_num_raw=0,data_num_bg=0,data_num_pi0=0;
    double eff_data=0,eff_data_err=0;
    double eff_mc=0,eff_mc_err=0;
    double ratio=0,ratio_err=0;
    double mu_data=0,sigma_data=0,mu_mc=0,sigma_mc=0;
};

double hist_integral_window(const TH1D* h,double lo,double hi) {
    if (!h) return 0.0;
    double s=0;
    for (int i=1;i<=h->GetNbinsX();i++) {
        const double x=h->GetBinCenter(i);
        if (x>lo && x<hi) s+=h->GetBinContent(i);
    } // endfor
    return s;
}

IntegratedEfficiencyResult integrated_efficiency(
        const std::vector<std::unique_ptr<ValComponent>>& vv,
        const NormDerivation& R,int ir,
        std::unique_ptr<TH1D>* data_out=nullptr,
        std::unique_ptr<TH1D>* pi0mc_out=nullptr) {
    IntegratedEfficiencyResult out;
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* a=find_val_component(vv,"aaogen");
    const ValComponent* c=find_val_component(vv,"clasdis");
    const ValComponent* d=find_val_component(vv,"dvcsgen");
    if (!data || !a || !c || !d) return out;

    const auto comp=coarse_composition(vv,R.nominal,ir);
    if (!(comp.f_pi0>0)) return out;
    const double fc=comp.f_clasdis_pi0;
    out.f_pi0=comp.f_pi0;

    const auto& rd=data->coarse[ir];
    const auto& ra=a->coarse[ir];
    const auto& rc=c->coarse[ir];
    const auto& rv=d->coarse[ir];
    out.data_denom=rd.denom_rows;
    out.mc_pi0_denom=R.nominal.aao*ra.denom_rows +
                     R.nominal.clasdis*fc*rc.denom_rows;

    if (!rd.residual || !ra.residual || !rc.residual || !rv.residual || rd.denom_rows<=0)
        return out;

    std::unique_ptr<TH1D> hpi0((TH1D*)ra.residual->Clone(Form("pi0mc_%d",ir)));
    hpi0->SetDirectory(nullptr);
    hpi0->Scale(R.nominal.aao);
    std::unique_ptr<TH1D> hcpi((TH1D*)rc.residual->Clone(Form("cpi_%d",ir)));
    hcpi->SetDirectory(nullptr); hcpi->Scale(R.nominal.clasdis*fc);
    hpi0->Add(hcpi.get());

    std::unique_ptr<TH1D> hbg((TH1D*)rv.residual->Clone(Form("bgmc_%d",ir)));
    hbg->SetDirectory(nullptr); hbg->Scale(R.nominal.dvcs);
    std::unique_ptr<TH1D> hcback((TH1D*)rc.residual->Clone(Form("cbg_%d",ir)));
    hcback->SetDirectory(nullptr); hcback->Scale(R.nominal.clasdis*(1.0-fc));
    hbg->Add(hcback.get());

    std::unique_ptr<TH1D> hd((TH1D*)rd.residual->Clone(Form("dataeff_%d",ir)));
    hd->SetDirectory(nullptr);

    FitResult fd=fit_valerii_residual(hd.get());
    FitResult fm=fit_valerii_residual(hpi0.get());
    out.data_fit_valid=fd.valid;
    out.mc_fit_valid=fm.valid;
    out.data_fit_reason=fd.reason;
    out.mc_fit_reason=fm.reason;

    if (data_out) {
        data_out->reset((TH1D*)hd->Clone(Form("display_data_%d",ir)));
        (*data_out)->SetDirectory(nullptr);
    }
    if (pi0mc_out) {
        pi0mc_out->reset((TH1D*)hpi0->Clone(Form("display_pi0mc_%d",ir)));
        (*pi0mc_out)->SetDirectory(nullptr);
    }

    if (!fd.valid || !fm.valid) return out;

    const double data_num=hist_integral_window(hd.get(),fd.mean-3*fd.sigma,fd.mean+3*fd.sigma);
    const double bg_denom=R.nominal.dvcs*rv.denom_rows +
                          R.nominal.clasdis*(1.0-fc)*rc.denom_rows;
    const double bg_num=hist_integral_window(hbg.get(),fd.mean-3*fd.sigma,fd.mean+3*fd.sigma);
    const double bg_eff=(bg_denom>0 ? bg_num/bg_denom : 0.0);

    const double data_pi0_denom=comp.f_pi0*double(rd.denom_rows);
    const double predicted_bg_in_data=(1.0-comp.f_pi0)*double(rd.denom_rows)*bg_eff;
    const double data_pi0_num=data_num-predicted_bg_in_data;

    const double mc_pi0_denom=R.nominal.aao*ra.denom_rows +
                              R.nominal.clasdis*fc*rc.denom_rows;
    const double mc_pi0_num=hist_integral_window(hpi0.get(),fm.mean-3*fm.sigma,fm.mean+3*fm.sigma);

    if (!(data_pi0_denom>0 && mc_pi0_denom>0 && mc_pi0_num>0)) return out;

    out.valid=true;
    out.f_pi0=comp.f_pi0;
    out.data_denom=rd.denom_rows;
    out.data_num_raw=data_num;
    out.data_num_bg=predicted_bg_in_data;
    out.data_num_pi0=data_pi0_num;
    out.eff_data=data_pi0_num/data_pi0_denom;
    out.eff_mc=mc_pi0_num/mc_pi0_denom;
    out.ratio=(out.eff_mc>0 ? out.eff_data/out.eff_mc : 0);

    // Statistical uncertainties only.  Composition fractions and normalization
    // constants are held fixed here; their variation is treated later as a
    // systematic uncertainty.
    if (out.eff_data>=0.0 && out.eff_data<=1.0 && data_pi0_denom>0.0)
        out.eff_data_err=std::sqrt(out.eff_data*(1.0-out.eff_data)/data_pi0_denom);

    // Effective MC denominator for two constant-weight pi0-bearing components.
    const double wa=R.nominal.aao;
    const double wc=R.nominal.clasdis*fc;
    const double sw=wa*double(ra.denom_rows)+wc*double(rc.denom_rows);
    const double sw2=wa*wa*double(ra.denom_rows)+wc*wc*double(rc.denom_rows);
    const double neff=(sw2>0 ? sw*sw/sw2 : 0.0);
    if (out.eff_mc>=0.0 && out.eff_mc<=1.0 && neff>0.0)
        out.eff_mc_err=std::sqrt(out.eff_mc*(1.0-out.eff_mc)/neff);

    if (out.ratio>0 && out.eff_data>0 && out.eff_mc>0) {
        const double rd=out.eff_data_err/out.eff_data;
        const double rm=out.eff_mc_err/out.eff_mc;
        out.ratio_err=out.ratio*std::sqrt(rd*rd+rm*rm);
    } // endif

    out.mu_data=fd.mean; out.sigma_data=fd.sigma;
    out.mu_mc=fm.mean; out.sigma_mc=fm.sigma;

    return out;
}

void draw_efficiency_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                             const NormDerivation& Rfd,
                             const NormDerivation& Rft,
                             const std::string& dir) {
    std::array<IntegratedEfficiencyResult,CR_N> rr;
    std::array<std::unique_ptr<TH1D>,CR_N> hd,hm;

    for (int ir=0;ir<CR_N;ir++) {
        const NormDerivation& R=(ir<2 ? Rfd : Rft);
        rr[ir]=integrated_efficiency(vv,R,ir,&hd[ir],&hm[ir]);
    } // endfor

    // Residuals.
    {
        TCanvas c("c_eff_residuals","",1450,1050);
        c.Divide(2,2);
        for (int ir=0;ir<CR_N;ir++) {
            c.cd(ir+1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.18); gPad->SetTopMargin(0.13);
            if (!hd[ir] || !hm[ir]) {
                TLatex t0; t0.SetNDC(); t0.SetTextFont(42); t0.SetTextSize(0.050);
                t0.DrawLatex(0.18,0.55,"No residual histogram available");
                continue;
            }
            // Rebin display copies only; fitting/counting above uses the native binning.
            hd[ir]->Rebin(4);
            hm[ir]->Rebin(4);
            hd[ir]->SetStats(0);
            hd[ir]->SetMarkerStyle(20);
            hd[ir]->SetMarkerSize(0.45);
            hd[ir]->SetLineColor(kBlack);
            hm[ir]->SetLineColor(kBlue+1);
            hm[ir]->SetLineWidth(2);
            hd[ir]->SetMaximum(1.25*std::max(hd[ir]->GetMaximum(),hm[ir]->GetMaximum()));
            hd[ir]->GetXaxis()->SetTitle("#Delta p_{#gamma2} (GeV)");
            hd[ir]->GetYaxis()->SetTitle("Candidates / display bin");
            hd[ir]->Draw("E1");
            hm[ir]->Draw("HIST SAME");
            hd[ir]->Draw("E1 SAME");
            TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.048);
            tx.DrawLatex(0.16,0.93,Form("(%c) %s",'a'+ir,CR_LABEL[ir]));
            if (!rr[ir].valid) {
                tx.SetTextSize(0.030);
                tx.DrawLatex(0.16,0.84,Form("Data fit: %s",rr[ir].data_fit_reason.c_str()));
                tx.DrawLatex(0.16,0.78,Form("MC fit: %s",rr[ir].mc_fit_reason.c_str()));
                tx.DrawLatex(0.16,0.72,Form("Data denominator: %.0f",rr[ir].data_denom));
            }
            if (ir==0) {
                TLegend* leg=new TLegend(0.63,0.72,0.93,0.88);
                leg->SetBorderSize(0); leg->SetFillStyle(0);
                leg->AddEntry(hd[ir].get(),"Data","lep");
                leg->AddEntry(hm[ir].get(),"Normalized #pi^{0} MC","l");
                leg->Draw();
            } // endif
        } // endfor
        c.SaveAs((dir+"/residuals.png").c_str());
    }

    // Summary with statistical uncertainties only.
    {
        TCanvas c("c_eff_summary","",1500,700);
        c.Divide(2,1);

        c.cd(1);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.24);
        gPad->SetTopMargin(0.13);

        TGraphErrors gd,gm;
        gd.SetMarkerStyle(20); gd.SetMarkerSize(1.20);
        gm.SetMarkerStyle(24); gm.SetMarkerSize(1.20);

        for (int ir=0;ir<CR_N;ir++) if (rr[ir].valid) {
            int n=gd.GetN();
            gd.SetPoint(n,ir+0.5,rr[ir].eff_data);
            gd.SetPointError(n,0.0,rr[ir].eff_data_err);

            n=gm.GetN();
            gm.SetPoint(n,ir+0.5,rr[ir].eff_mc);
            gm.SetPointError(n,0.0,rr[ir].eff_mc_err);
        } // endfor

        TH1D axis("axis_eff",";Detector / probe-energy region;Photon efficiency",CR_N,0,CR_N);
        axis.SetStats(0);
        axis.SetMinimum(0.0);
        axis.SetMaximum(1.05);
        for (int ir=0;ir<CR_N;ir++) axis.GetXaxis()->SetBinLabel(ir+1,CR_LABEL[ir]);
        axis.GetXaxis()->SetLabelSize(0.040);
        axis.Draw();
        gd.Draw("P SAME");
        gm.Draw("P SAME");

        TLegend leg(0.50,0.70,0.91,0.84);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(&gd,"Data (#pi^{0}-corrected)","p");
        leg.AddEntry(&gm,"#pi^{0} MC","p");
        leg.Draw();

        TLatex tx;
        tx.SetNDC();
        tx.SetTextFont(42);
        tx.SetTextSize(0.050);
        tx.DrawLatex(0.16,0.93,"(a) Integrated 3#sigma photon efficiencies");

        c.cd(2);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.24);
        gPad->SetTopMargin(0.13);

        TH1D axisr("axis_ratio",";Detector / probe-energy region;#epsilon_{data}/#epsilon_{MC}",
                   CR_N,0,CR_N);
        axisr.SetStats(0);
        axisr.SetMinimum(0.45);
        axisr.SetMaximum(1.55);
        for (int ir=0;ir<CR_N;ir++) axisr.GetXaxis()->SetBinLabel(ir+1,CR_LABEL[ir]);
        axisr.GetXaxis()->SetLabelSize(0.040);
        axisr.Draw();

        TGraphErrors gr;
        gr.SetMarkerStyle(20);
        gr.SetMarkerSize(1.30);
        for (int ir=0;ir<CR_N;ir++) if (rr[ir].valid) {
            const int n=gr.GetN();
            gr.SetPoint(n,ir+0.5,rr[ir].ratio);
            gr.SetPointError(n,0.0,rr[ir].ratio_err);
        } // endfor
        gr.Draw("P SAME");

        TLine one(0,1.0,CR_N,1.0);
        one.SetLineStyle(2);
        one.Draw();

        TLatex tx2;
        tx2.SetNDC();
        tx2.SetTextFont(42);
        tx2.SetTextSize(0.050);
        tx2.DrawLatex(0.16,0.93,"(b) Data / MC photon-efficiency ratio");

        c.SaveAs((dir+"/efficiency_summary.png").c_str());
    }

    std::ofstream csv(dir+"/summary.csv");
    csv << "region,valid,data_fit_valid,data_fit_reason,mc_fit_valid,mc_fit_reason,"
           "pi0_fraction,data_denom,mc_pi0_denom,data_num_raw,predicted_background_num,"
           "data_pi0_num,eff_data,eff_data_stat,eff_mc,eff_mc_stat,"
           "data_over_mc,data_over_mc_stat,mu_data,sigma_data,mu_mc,sigma_mc\n";
    for (int ir=0;ir<CR_N;ir++) {
        const auto& q=rr[ir];
        csv << CR_KEY[ir]<<","<<q.valid<<","
            <<q.data_fit_valid<<",\""<<q.data_fit_reason<<"\","
            <<q.mc_fit_valid<<",\""<<q.mc_fit_reason<<"\","
            <<q.f_pi0<<","<<q.data_denom<<","<<q.mc_pi0_denom<<","
            <<q.data_num_raw<<","<<q.data_num_bg<<","<<q.data_num_pi0<<","
            <<q.eff_data<<","<<q.eff_data_err<<","
            <<q.eff_mc<<","<<q.eff_mc_err<<","
            <<q.ratio<<","<<q.ratio_err<<","
            <<q.mu_data<<","<<q.sigma_data<<","<<q.mu_mc<<","<<q.sigma_mc<<"\n";
    } // endfor
}

void run_concise_analysis(const std::string& out) {
    concise_make_dirs(out);

    auto vv=build_val_components_parallel(out);

    auto have_sample=[&](const std::string& name)->bool {
        return find_val_component(vv,name)!=nullptr;
    };

    const bool complete =
        have_sample("data") &&
        have_sample("aaogen") &&
        have_sample("clasdis") &&
        have_sample("dvcsgen");

    if (!complete) {
        std::cerr
            << "\nERROR: stopping analysis because the cache is incomplete.\n"
            << "All four samples are required: data, aaogen, clasdis, dvcsgen.\n"
            << "No normalization, pi0-fraction, or efficiency products will be "
            << "made from a partial MC mixture.\n";
        return;
    } // endif

    const NormDerivation norm_fd=derive_normalization_concise(vv,false);
    const NormDerivation norm_ft=derive_normalization_concise(vv,true);

    draw_exclusivity_summary(vv,out+"/1_exclusivity/FD",false);
    draw_exclusivity_summary(vv,out+"/1_exclusivity/FT",true);
    draw_ft_cut_optimization(vv,norm_ft,out+"/1_exclusivity/FT");
    draw_ft_fitquality_scan(vv,out+"/1_exclusivity/FT");

    draw_normalization_summary(vv,norm_fd,out+"/2_normalization/FD",false);
    draw_normalization_summary(vv,norm_ft,out+"/2_normalization/FT",true);

    draw_pi0_summary(vv,norm_fd,norm_ft,out+"/3_pi0_fraction");
    draw_efficiency_summary(vv,norm_fd,norm_ft,out+"/4_efficiency");

    std::cout << "\nConcise output written to:\n"
              << "  " << out << "/1_exclusivity/FD and FT/\n"
              << "  " << out << "/2_normalization/FD and FT/\n"
              << "  " << out << "/3_pi0_fraction/\n"
              << "  " << out << "/4_efficiency/\n";
}


void run_valerii_fd_reproduction(const std::string& out) {
    std::cout << "\n============================================================\n"
              << " Valerii-style FD 7x3x6 reproduction\n"
              << "============================================================\n"
              << "MC events are unit weighted inside each component.\n"
              << "AAO/CLASDIS/DVCS normalization factors are derived from this run's template fits.\n"
              << "June normalization cuts are applied to the ep-gamma-X denominator.\n"
              << "Template fits use N-1 cuts and a common MC shift + extra Gaussian smearing nuisance.\n"
              << "Low E_gamma<2 GeV determines AAO+CLASDIS with DVCS excluded from the fit.\n"
              << "The 2-3 GeV transition region is closure-only; high E_gamma>3 GeV determines DVCS with AAO/CLASDIS fixed.\n"
              << "Template morphing is normalization-only; event-level MC remains unsmeared.\n"
              << "The skim MC::Event.weight branch is NOT used in this path.\n"
              << "Data/MC fits are independent in every p/theta/phi bin.\n"
              << "Primary correction convention: epsilon_data / epsilon_MC.\n"
              << "Delta-phi exclusivity cut: |delta phi_copl| < 5.7 deg using the Trento-plane residual.\n"
              << "Raw lab dphi and the underlying Trento-angle quantities remain QA diagnostics.\n"
              << "High-E DVCS normalization diagnostics include Delta t = t_p - t_gamma.\n"
              << "============================================================\n";
    auto vv=build_val_components_parallel(out);
    write_valerii_outputs(vv,out);
    gSystem->Exec(("rm -rf "+out+"/.valerii_workers").c_str());
}

} // namespace pe

void photon_efficiency_valerii_reproduction() {
    using namespace pe;

    gROOT->SetBatch(kTRUE);
    concise_publication_style();

    std::cout
        << "\n============================================================\n"
        << " Concise photon-efficiency analysis\n"
        << "============================================================\n"
        << "1) independent FD-probe and FT-probe exclusivity selection\n"
        << "2) independent FD-probe and FT-probe Valerii-style normalization\n"
        << "3) pi0 fraction of selected ep-gamma-X events\n"
        << "4) integrated FD/FT photon efficiency, split at 2 GeV\n"
        << "One parallel tree scan per sample; persistent cache on reruns.\n"
        << "Cache: /work/clas12/thayward/photon_efficiency/cache/concise\n"
        << "============================================================\n";

    pe::run_concise_analysis("output");
}


