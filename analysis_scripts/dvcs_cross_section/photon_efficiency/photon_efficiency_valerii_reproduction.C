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

void set_publication_style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);
    gStyle->SetTitleSize(0.050,"XYZ");
    gStyle->SetLabelSize(0.043,"XYZ");
    gStyle->SetTitleOffset(1.05,"X");
    gStyle->SetTitleOffset(1.22,"Y");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(2);
    gStyle->SetEndErrorSize(4);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
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
    if (f.IsZombie()) return false;
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
    long long truth_rows=0;
    long long truth_pi0_rows=0;
    std::unique_ptr<TH1D> residual_fit;
    std::unique_ptr<TH1D> residual_count;
};

struct ValComponent {
    std::string name;
    bool is_mc=false;
    long long entries=0;
    std::array<ValComponentBin,VAL_NBIN> bins;

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
    leg.Draw(); c.SetBottomMargin(0.22); c.SaveAs(file.c_str());
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
    std::array<long long,VAL_NBIN> truth_rows{};
    std::array<long long,VAL_NBIN> truth_pi0_rows{};

    // Do NOT write one TTree row per reconstructed candidate.  With millions of
    // skim rows that temporary representation can become multi-GB and exhaust
    // the ifarm home/output filesystem.  The analysis only needs the residual
    // distributions and their weighted sums, so persist fixed-size histograms.
    std::array<std::unique_ptr<TH1D>,VAL_NBIN> hfit;
    std::array<std::unique_ptr<TH1D>,VAL_NBIN> hcount;
    std::vector<std::unique_ptr<TH1D>> norm_pre, norm_after_mx2ep, norm_after_mx2ep_mx2eg,
        norm_after_mx2ep_mx2eg_dphi, norm_nminus1, norm_full, norm_low, norm_high;
    for (int io=0;io<NORM_NOBS;io++) {
        norm_pre.push_back(make_norm_hist(io,"pre",spec.name));
        norm_after_mx2ep.push_back(make_norm_hist(io,"after_mx2ep",spec.name));
        norm_after_mx2ep_mx2eg.push_back(make_norm_hist(io,"after_mx2ep_mx2eg",spec.name));
        norm_after_mx2ep_mx2eg_dphi.push_back(make_norm_hist(io,"after_mx2ep_mx2eg_dphi",spec.name));
        norm_nminus1.push_back(make_norm_hist(io,"nminus1",spec.name));
        norm_full.push_back(make_norm_hist(io,"full",spec.name));
        norm_low.push_back(make_norm_hist(io,"lowE",spec.name));
        norm_high.push_back(make_norm_hist(io,"highE",spec.name));
    } // endfor
    std::array<long long,6> norm_cutflow{{0,0,0,0,0,0}};
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

        const int ib=val_flat_bin(b.probe_corr_p,b.probe_corr_theta,b.probe_corr_phi);
        if (ib<0) continue;

        // Data-driven normalization diagnostics.  First retain the baseline
        // candidate shapes, then apply the June-2026 exclusivity cuts.  For
        // each fit observable, the low/high-E template uses an N-1 selection:
        // every exclusivity cut is active except a cut directly on that plotted
        // observable.  This prevents a hard cut from manufacturing agreement in
        // the very distribution used to determine a normalization factor.
        const double Eg=b.have_tag_corr_kin?b.tag_corr_p:std::numeric_limits<double>::quiet_NaN();
        const NormCutFlags ncf=norm_cut_flags(b);
        norm_cutflow[0]++;
        if (ncf.mx2_ep) {
            norm_cutflow[1]++;
            if (ncf.mx2_eg) {
                norm_cutflow[2]++;
                if (ncf.dphi_trento) {
                    norm_cutflow[3]++;
                    if (ncf.angle_gX) norm_cutflow[4]++;
                }
            }
        }
        if (ncf.all) norm_cutflow[5]++;

        for (int io=0;io<NORM_NOBS;io++) {
            const double x=norm_observable_value(b,io);
            if (!std::isfinite(x)) continue;
            norm_pre[io]->Fill(x);
            if (ncf.mx2_ep) norm_after_mx2ep[io]->Fill(x);
            if (ncf.mx2_ep && ncf.mx2_eg) norm_after_mx2ep_mx2eg[io]->Fill(x);
            if (ncf.mx2_ep && ncf.mx2_eg && ncf.dphi_trento)
                norm_after_mx2ep_mx2eg_dphi[io]->Fill(x);
            if (norm_pass_nminus1(ncf,io)) {
                norm_nminus1[io]->Fill(x);
                if (std::isfinite(Eg) && Eg<2.0) norm_low[io]->Fill(x);
                if (std::isfinite(Eg) && Eg>3.0) norm_high[io]->Fill(x);
            }
            if (ncf.all) norm_full[io]->Fill(x);
        } // endfor

        // The efficiency denominator itself uses the same exclusivity selection
        // as the normalization stage.  This is the ep-gamma-X sample intended
        // to be enriched in ep-pi0-like events before asking whether the probe
        // photon was reconstructed.
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
        if (spec.is_mc && b.have_truth) {
            truth_rows[ib]++;
            if (b.truth_probe_pid==22 && b.truth_probe_parent==111)
                truth_pi0_rows[ib]++;
        }
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
    Long64_t drows=0,dtruth=0,dtruth_pi0=0;
    Double_t dsumw=0,dsumw2=0;
    denominators.Branch("bin",&dbin,"bin/I");
    denominators.Branch("rows",&drows,"rows/L");
    denominators.Branch("sumw",&dsumw,"sumw/D");
    denominators.Branch("sumw2",&dsumw2,"sumw2/D");
    denominators.Branch("truth_rows",&dtruth,"truth_rows/L");
    denominators.Branch("truth_pi0_rows",&dtruth_pi0,"truth_pi0_rows/L");
    for (dbin=0;dbin<VAL_NBIN;dbin++) {
        drows=rows[dbin]; dsumw=sumw[dbin]; dsumw2=sumw2[dbin];
        dtruth=truth_rows[dbin]; dtruth_pi0=truth_pi0_rows[dbin];
        denominators.Fill();
    } // endfor
    denominators.Write();

    TDirectory* rd=f.mkdir("residuals");
    if (!rd) { f.Close(); return false; }
    rd->cd();
    for (int ib=0;ib<VAL_NBIN;ib++) {
        hfit[ib]->Write(Form("fit_b%03d",ib));
        hcount[ib]->Write(Form("count_b%03d",ib));
    } // endfor
    f.cd();
    TDirectory* nd=f.mkdir("normalization");
    if (nd) {
        nd->cd();
        for (int io=0;io<NORM_NOBS;io++) {
            norm_pre[io]->Write(Form("pre_%s",NORM_OBS[io].key));
            norm_after_mx2ep[io]->Write(Form("after_mx2ep_%s",NORM_OBS[io].key));
            norm_after_mx2ep_mx2eg[io]->Write(Form("after_mx2ep_mx2eg_%s",NORM_OBS[io].key));
            norm_after_mx2ep_mx2eg_dphi[io]->Write(Form("after_mx2ep_mx2eg_dphi_%s",NORM_OBS[io].key));
            norm_nminus1[io]->Write(Form("nminus1_%s",NORM_OBS[io].key));
            norm_full[io]->Write(Form("full_%s",NORM_OBS[io].key));
            norm_low[io]->Write(Form("lowE_%s",NORM_OBS[io].key));
            norm_high[io]->Write(Form("highE_%s",NORM_OBS[io].key));
        } // endfor
        f.cd();
    }
    for (int ic=0;ic<6;ic++) put_param<Long64_t>(&f,Form("norm_cutflow_%d",ic),norm_cutflow[ic]);
    put_param<int>(&f,"normalization_branches_complete",
        (b.have_tag_corr_kin && b.have_Mx2_epg_corr && b.have_beam_energy && b.have_e_kin && b.have_p_corr_kin)?1:0);
    const Int_t nwrite=f.Write();
    f.Close();
    if (nwrite<=0 || gSystem->AccessPathName(path.c_str())) {
        std::cerr << "ERROR: failed to persist Valerii worker file " << path << "\n";
        return false;
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
    Int_t ib=-1;
    Long64_t rows=0,truthrows=0,truthpi0=0;
    Double_t sw=0,sw2=0;
    den->SetBranchAddress("bin",&ib);
    den->SetBranchAddress("rows",&rows);
    den->SetBranchAddress("sumw",&sw);
    den->SetBranchAddress("sumw2",&sw2);
    const bool have_truthrows=(den->GetBranch("truth_rows")!=nullptr);
    const bool have_truthpi0=(den->GetBranch("truth_pi0_rows")!=nullptr);
    if (have_truthrows) den->SetBranchAddress("truth_rows",&truthrows);
    if (have_truthpi0) den->SetBranchAddress("truth_pi0_rows",&truthpi0);
    for (Long64_t i=0;i<den->GetEntries();i++) {
        den->GetEntry(i);
        if (ib<0 || ib>=VAL_NBIN) continue;
        v->bins[ib].denom_rows=rows;
        v->bins[ib].denom_w=sw;
        v->bins[ib].denom_w2=sw2;
        v->bins[ib].truth_rows=have_truthrows?truthrows:0;
        v->bins[ib].truth_pi0_rows=have_truthpi0?truthpi0:0;
    } // endfor

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
    int ncomplete=0; get_param<int>(&f,"normalization_branches_complete",ncomplete); v->normalization_branches_complete=(ncomplete!=0);
    auto* nd=dynamic_cast<TDirectory*>(f.Get("normalization"));
    if (nd) {
        for (int io=0;io<NORM_NOBS;io++) {
            auto clone_one=[&](const char* region)->std::unique_ptr<TH1D> {
                auto* h=dynamic_cast<TH1D*>(nd->Get(Form("%s_%s",region,NORM_OBS[io].key)));
                if (!h) return nullptr;
                std::unique_ptr<TH1D> q(dynamic_cast<TH1D*>(h->Clone(Form("%s_%s_%s",v->name.c_str(),region,NORM_OBS[io].key))));
                if (q) q->SetDirectory(nullptr); return q;
            };
            v->norm_pre.push_back(clone_one("pre"));
            v->norm_after_mx2ep.push_back(clone_one("after_mx2ep"));
            v->norm_after_mx2ep_mx2eg.push_back(clone_one("after_mx2ep_mx2eg"));
            v->norm_after_mx2ep_mx2eg_dphi.push_back(clone_one("after_mx2ep_mx2eg_dphi"));
            v->norm_nminus1.push_back(clone_one("nminus1"));
            v->norm_full.push_back(clone_one("full"));
            v->norm_lowE.push_back(clone_one("lowE"));
            v->norm_highE.push_back(clone_one("highE"));
        } // endfor
    }
    for (int ic=0;ic<6;ic++) {
        Long64_t q=0; get_param<Long64_t>(&f,Form("norm_cutflow_%d",ic),q); v->norm_cutflow[ic]=q;
    } // endfor
    f.Close();
    return v;
}
std::vector<std::unique_ptr<ValComponent>> build_val_components_parallel(const std::string& out) {
    const SampleSpec all_specs[] = {
        {"data",DATA_DIR,false},
        {"aaogen",AAOGEN_DIR,true},
        {"clasdis",CLASDIS_DIR,true},
        {"dvcsgen",DVCSGEN_DIR,true}
    };
    const std::string wdir=out+"/.valerii_workers";
    // A previous interrupted/no-space run can leave corrupt worker files behind.
    // Never attempt to recover/reuse them.
    gSystem->Exec(("rm -rf "+wdir).c_str());
    gSystem->mkdir(wdir.c_str(),true);

    struct VChild { pid_t pid; SampleSpec spec; };
    std::vector<VChild> children;
    std::vector<SampleSpec> specs;
    std::cout.flush(); std::cerr.flush();

    for (const auto& spec:all_specs) {
        if (!has_root_files(spec.dir)) continue;
        specs.push_back(spec);
        const pid_t pid=fork();
        if (pid==0) {
            const std::string path=wdir+"/"+spec.name+".root";
            const bool ok=analyze_val_component_worker(spec,path);
            std::cout.flush(); std::cerr.flush();
            _exit(ok?0:2);
        } else if (pid>0) {
            children.push_back({pid,spec});
            std::cout << "[VALERII] launched " << spec.name << " worker pid=" << pid << "\n";
        } else {
            std::cerr << "WARNING: Valerii fork failed for " << spec.name << "; running sequentially.\n";
            analyze_val_component_worker(spec,wdir+"/"+spec.name+".root");
        }
    } // endfor

    for (const auto& ch:children) {
        int status=0; waitpid(ch.pid,&status,0);
        if (!WIFEXITED(status) || WEXITSTATUS(status)!=0)
            std::cerr << "WARNING: Valerii worker failed for " << ch.spec.name << "\n";
        else
            std::cout << "[VALERII] finished " << ch.spec.name << "\n";
    } // endfor

    std::vector<std::unique_ptr<ValComponent>> outv;
    for (const auto& spec:specs) {
        auto p=load_val_component(wdir+"/"+spec.name+".root");
        if (p) outv.push_back(std::move(p));
    } // endfor
    return outv;
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


struct Pi0TruthSummary {
    long long truth_rows=0;
    long long pi0_rows=0;
    double fraction=-1.0;
};

Pi0TruthSummary component_pi0_truth(const ValComponent* v) {
    Pi0TruthSummary r;
    if (!v || !v->is_mc) return r;
    for (int ib=0;ib<VAL_NBIN;ib++) {
        r.truth_rows += v->bins[ib].truth_rows;
        r.pi0_rows += v->bins[ib].truth_pi0_rows;
    }
    if (r.truth_rows>0) r.fraction=double(r.pi0_rows)/double(r.truth_rows);
    return r;
}

void write_pi0_truth_composition(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                 const NormDerivation& nd,
                                 const std::string& dir) {
    gSystem->mkdir(dir.c_str(),true);

    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");

    struct C { const char* name; const ValComponent* v; double scale; };
    const C cc[]={{"aaogen",aao,nd.nominal.aao},
                  {"clasdis",cls,nd.nominal.clasdis},
                  {"dvcsgen",dvc,nd.nominal.dvcs}};

    std::ofstream comp(dir+"/mc_component_pi0_truth.csv");
    comp << "component,truth_classified_rows,truth_pi0_rows,pi0_truth_fraction,nominal_scale\n";

    double model_total=0.0, model_pi0=0.0;
    for (const auto& c:cc) {
        const Pi0TruthSummary r=component_pi0_truth(c.v);
        comp << c.name << "," << r.truth_rows << "," << r.pi0_rows << ","
             << r.fraction << "," << c.scale << "\n";
        model_total += c.scale*double(r.truth_rows);
        model_pi0 += c.scale*double(r.pi0_rows);
    }
    comp.close();

    std::ofstream bins(dir+"/normalized_epgamma_pi0_fraction_by_bin.csv");
    bins << "bin,ip,itheta,iphi,p_lo,p_hi,theta_lo,theta_hi,phi_lo,phi_hi,"
            "weighted_total_truth,weighted_pi0_truth,pi0_fraction,"
            "aao_pi0_fraction,clasdis_pi0_fraction,dvcs_pi0_fraction\n";

    TH2D hmap("h_pi0_fraction_ptheta",
              ";Missing-probe p (GeV);Missing-probe #theta (deg)",
              VAL_NP,VAL_P_EDGES,VAL_NT,VAL_T_EDGES);
    hmap.SetStats(0);

    for (int ip=0;ip<VAL_NP;ip++) {
        for (int it=0;it<VAL_NT;it++) {
            double totpt=0,pi0pt=0;
            for (int iph=0;iph<VAL_NPH;iph++) {
                const int ib=(ip*VAL_NT+it)*VAL_NPH+iph;
                double tot=0,pi0=0,fa=-1,fc=-1,fd=-1;
                for (const auto& c:cc) {
                    if (!c.v) continue;
                    const auto& b=c.v->bins[ib];
                    if (b.truth_rows<=0) continue;
                    const double f=double(b.truth_pi0_rows)/double(b.truth_rows);
                    if (std::string(c.name)=="aaogen") fa=f;
                    else if (std::string(c.name)=="clasdis") fc=f;
                    else if (std::string(c.name)=="dvcsgen") fd=f;
                    tot += c.scale*double(b.truth_rows);
                    pi0 += c.scale*double(b.truth_pi0_rows);
                }
                const double frac=(tot>0?pi0/tot:-1.0);
                bins << ib << "," << ip << "," << it << "," << iph << ","
                     << VAL_P_EDGES[ip] << "," << VAL_P_EDGES[ip+1] << ","
                     << VAL_T_EDGES[it] << "," << VAL_T_EDGES[it+1] << ","
                     << VAL_PH_EDGES[iph] << "," << VAL_PH_EDGES[iph+1] << ","
                     << tot << "," << pi0 << "," << frac << ","
                     << fa << "," << fc << "," << fd << "\n";
                totpt += tot; pi0pt += pi0;
            }
            if (totpt>0) hmap.SetBinContent(ip+1,it+1,pi0pt/totpt);
        }
    }
    bins.close();

    TCanvas c("c_pi0_truth_summary","",1500,650);
    c.Divide(2,1);

    c.cd(1);
    gPad->SetLeftMargin(0.14);
    TH1D hcomp("h_component_pi0_truth",
               ";MC component;Truth #pi^{0}-photon fraction",3,0,3);
    hcomp.SetStats(0); hcomp.SetMinimum(-0.05); hcomp.SetMaximum(1.05);
    hcomp.SetMarkerStyle(20); hcomp.SetMarkerSize(1.4);
    for (int i=0;i<3;i++) {
        const Pi0TruthSummary r=component_pi0_truth(cc[i].v);
        hcomp.GetXaxis()->SetBinLabel(i+1,cc[i].name);
        if (r.fraction>=0) hcomp.SetBinContent(i+1,r.fraction);
    }
    hcomp.Draw("P");
    TLine one(0,1.0,3,1.0); one.SetLineStyle(3); one.Draw();
    TLatex t1; t1.SetNDC(); t1.SetTextFont(42); t1.SetTextSize(0.050);
    t1.DrawLatex(0.16,0.91,"(a) Truth composition after all cuts");

    c.cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.16);
    hmap.SetMinimum(0.0); hmap.SetMaximum(1.0);
    hmap.GetZaxis()->SetTitle("Normalized-model #pi^{0} fraction");
    hmap.Draw("COLZ TEXT");
    TLatex t2; t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.050);
    t2.DrawLatex(0.16,0.91,"(b) #pi^{0} fraction vs missing-probe kinematics");

    c.SaveAs((dir+"/epgamma_pi0_fraction_summary.png").c_str());

    std::ofstream txt(dir+"/epgamma_pi0_fraction_summary.txt");
    txt << "Truth-based pi0 composition of selected ep-gamma-X candidates\n"
        << "=============================================================\n"
        << "Classification: truth_probe_pid==22 and truth_probe_parent==111.\n"
        << "Nominal data-derived AAO/CLASDIS/DVCS factors weight the MC mixture.\n\n";
    for (const auto& c0:cc) {
        const Pi0TruthSummary r=component_pi0_truth(c0.v);
        txt << c0.name << ": " << r.pi0_rows << "/" << r.truth_rows
            << " = " << r.fraction << "; scale=" << c0.scale << "\n";
    }
    txt << "\nIntegrated normalized-model pi0 fraction = "
        << (model_total>0?model_pi0/model_total:-1.0) << "\n";
    txt.close();
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

    // Compact 1x2 summary.  Extend the right-panel range slightly beyond
    // [0,1] so solutions pinned exactly to 0 or 1 remain visible.
    {
        TCanvas c("c_lowE_pi0_balance_summary","",1500,650);
        c.Divide(2,1);

        c.cd(1);
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.25);
        TH1D htot("h_lowE_pi0_total_closure",
                  ";Normalization observable;(AAO+CLASDIS fitted yield) / data yield",
                  (int)rows.size(),0,(int)rows.size());
        htot.SetStats(0); htot.SetMarkerStyle(20); htot.SetMarkerSize(1.30);
        for (int i=0;i<(int)rows.size();i++) {
            htot.GetXaxis()->SetBinLabel(i+1,rows[i].obs.c_str());
            htot.SetBinContent(i+1,rows[i].model_over_data);
        }
        htot.SetMinimum(0.78); htot.SetMaximum(1.08);
        htot.GetXaxis()->LabelsOption("v");
        htot.Draw("P");
        TLine one(0,1.0,rows.size(),1.0);
        one.SetLineStyle(2); one.SetLineWidth(2); one.Draw();
        TLatex ta; ta.SetNDC(); ta.SetTextFont(42); ta.SetTextSize(0.050);
        ta.DrawLatex(0.17,0.91,"(a) Combined low-E normalization");

        c.cd(2);
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.25);
        TH1D hfrac("h_lowE_pi0_balance",
                   ";Normalization observable;AAO fraction of fitted AAO+CLASDIS yield",
                   (int)rows.size(),0,(int)rows.size());
        hfrac.SetStats(0); hfrac.SetMarkerStyle(20); hfrac.SetMarkerSize(1.30);
        for (int i=0;i<(int)rows.size();i++) {
            hfrac.GetXaxis()->SetBinLabel(i+1,rows[i].obs.c_str());
            hfrac.SetBinContent(i+1,rows[i].aao_fraction);
        }
        hfrac.SetMinimum(-0.08); hfrac.SetMaximum(1.08);
        hfrac.GetXaxis()->LabelsOption("v");
        hfrac.Draw("P");
        TLine zero(0,0.0,rows.size(),0.0), unity(0,1.0,rows.size(),1.0);
        zero.SetLineStyle(3); unity.SetLineStyle(3);
        zero.Draw(); unity.Draw();
        TLatex tb; tb.SetNDC(); tb.SetTextFont(42); tb.SetTextSize(0.050);
        tb.DrawLatex(0.17,0.91,"(b) AAO/CLASDIS tradeoff");

        c.SaveAs((dir+"/lowE_pi0_balance_summary.png").c_str());
    }
}



void draw_compact_step1_selection_summary(const std::vector<std::unique_ptr<ValComponent>>& vv,
                                          const std::string& file) {
    const ValComponent* data=find_val_component(vv,"data");
    const ValComponent* aao=find_val_component(vv,"aaogen");
    const ValComponent* cls=find_val_component(vv,"clasdis");
    const ValComponent* dvc=find_val_component(vv,"dvcsgen");
    if (!data || !aao || !cls || !dvc) return;

    struct P { int obs,stage; double c1,c2; bool two; const char* lab; };
    const P pp[]={
        {NORM_MX2_EP,0,NORM_MX2_EP_MIN,NORM_MX2_EP_MAX,true,"(a) M_{X}^{2}(ep)"},
        {NORM_MX2_EG,1,NORM_MX2_EG_MIN,0,false,"(b) M_{X}^{2}(e#gamma)"},
        {NORM_DPHI_TRENTO_SHIFT180,2,-NORM_DPHI_TRENTO_MAX,NORM_DPHI_TRENTO_MAX,true,"(c) Trento coplanarity"},
        {NORM_ANGLE_GX,3,NORM_ANGLE_GX_MAX,0,false,"(d) angle(#gamma,X)"}
    };
    auto geth=[](const ValComponent* v,int st,int io)->const TH1D* {
        if (!v) return nullptr;
        if (st==0 && io<(int)v->norm_pre.size()) return v->norm_pre[io].get();
        if (st==1 && io<(int)v->norm_after_mx2ep.size()) return v->norm_after_mx2ep[io].get();
        if (st==2 && io<(int)v->norm_after_mx2ep_mx2eg.size()) return v->norm_after_mx2ep_mx2eg[io].get();
        if (st==3 && io<(int)v->norm_after_mx2ep_mx2eg_dphi.size()) return v->norm_after_mx2ep_mx2eg_dphi[io].get();
        return nullptr;
    };

    TCanvas c("c_step1_compact","",1500,1100);
    c.Divide(2,2);
    std::vector<std::unique_ptr<TH1D>> keep;
    const ValComponent* vv4[]={data,aao,cls,dvc};
    const int cols[]={kBlack,kRed+1,kOrange+7,kGreen+2};
    const char* labs[]={"Data","AAOgen","CLASDIS","DVCSgen"};

    for (int ip=0;ip<4;ip++) {
        c.cd(ip+1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        std::vector<TH1D*> hs;
        double ymax=0;
        for (int j=0;j<4;j++) {
            const TH1D* h0=geth(vv4[j],pp[ip].stage,pp[ip].obs);
            if (!h0) { hs.push_back(nullptr); continue; }
            std::unique_ptr<TH1D> h((TH1D*)h0->Clone(Form("compact_%d_%d",ip,j)));
            h->SetDirectory(nullptr);
            const double q=h->Integral(); if (q>0) h->Scale(1.0/q);
            h->SetStats(0);
            if (j==0) {
                h->SetMarkerStyle(20); h->SetMarkerSize(0.55);
                h->SetLineColor(kBlack); h->SetMarkerColor(kBlack);
            } else { h->SetLineColor(cols[j]); h->SetLineWidth(2); }
            ymax=std::max(ymax,h->GetMaximum());
            hs.push_back(h.get()); keep.push_back(std::move(h));
        }
        if (!hs[0]) continue;
        hs[0]->SetTitle(""); hs[0]->GetYaxis()->SetTitle("Unit-area candidates");
        hs[0]->SetMaximum(1.28*ymax);
        hs[0]->Draw("E1");
        for (int j=1;j<4;j++) if (hs[j]) hs[j]->Draw("HIST SAME");
        hs[0]->Draw("E1 SAME");

        const double yy=1.20*ymax;
        TLine* l1=new TLine(pp[ip].c1,0,pp[ip].c1,yy);
        l1->SetLineStyle(2); l1->SetLineWidth(2); l1->Draw();
        if (pp[ip].two) {
            TLine* l2=new TLine(pp[ip].c2,0,pp[ip].c2,yy);
            l2->SetLineStyle(2); l2->SetLineWidth(2); l2->Draw();
        }
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.047);
        tx.DrawLatex(0.16,0.91,pp[ip].lab);
        if (ip==0) {
            TLegend* leg=new TLegend(0.58,0.64,0.88,0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            for (int j=0;j<4;j++) if (hs[j]) leg->AddEntry(hs[j],labs[j],j==0?"lep":"l");
            leg->Draw();
        }
    }
    c.SaveAs(file.c_str());
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
    draw_compact_step1_selection_summary(vv,od+"/step1_selection_summary.png");

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
    write_pi0_truth_composition(vv,norm,out+"/normalization/pi0_truth_composition");
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
    set_publication_style();

    const std::string out="output";
    reset_output(out);

    std::cout
        << "\n============================================================\n"
        << " Photon tag-and-probe stage-1 diagnostics\n"
        << "============================================================\n"
        << "FD: fully integrated\n"
        << "FT: projected-fiducial integrated + 0.4-2 GeV + >=2 GeV\n"
        << "Matching: Delta p_gamma2 = p_rec - p_miss\n"
        << "Primary numbers: unweighted RAW RECOVERY FRACTIONS (not efficiencies)\n"
        << "Stage-1 MC::Event.weight: QA only; Valerii FD path uses unit MC events x component normalization\n"
        << "Execution: independent samples analyzed in parallel child processes\n"
        << "Outputs are rebuilt in ./output each invocation.\n"
        << "============================================================\n";

    std::vector<std::unique_ptr<SampleResult>> samples=analyze_samples_parallel(out);

    // Run the second process-level pass before creating graphics in the parent.
    // Forking before canvases/files are opened avoids inheriting unnecessary ROOT
    // graphics state into the Valerii workers.
    run_valerii_fd_reproduction(out);

    std::ofstream csv(out+"/integrated_results.csv");
    write_region_csv_header(csv);
    for (const auto& sp : samples) {
        append_region_csv(csv,*sp,sp->fd);
        append_region_csv(csv,*sp,sp->ft_all);
        append_region_csv(csv,*sp,sp->ft_low);
        append_region_csv(csv,*sp,sp->ft_high);
    } // endfor
    csv.close();

    save_combined_residuals(samples,out);
    save_combined_kinematics(samples,out);
    save_combined_ft_projection(samples,out);
    write_window_scan_all(samples,out);
    write_all_summary(samples,out);
    write_root(samples,out);

    // Temporary worker files are implementation details; keep output compact.
    gSystem->Exec(("rm -rf "+out+"/.workers").c_str());

    std::cout << "\nFinished. Output products actually present on disk:\n";
    const char* expected[] = {
        "analysis_summary.txt",
        "integrated_results.csv",
        "parent_window_scan.csv",
        "delta_p_residuals.png",
        "denominator_kinematics.png",
        "FT_projected_xy.png",
        "analysis_histograms.root",
        "valerii_fd_summary.txt",
        "valerii_fd_results.csv",
        "valerii_fd_component_weight_qa.csv",
        "valerii_fd_data_efficiency_2sigma.png",
        "valerii_fd_weighted_mc_efficiency_2sigma.png",
        "valerii_fd_data_over_mc_correction_2sigma.png",
        "valerii_fd_data_sigma.png",
        "valerii_fd_weighted_mc_sigma.png",
        "valerii_fd_histograms.root",
        "normalization/normalization_summary.txt",
        "normalization/normalization_fit_results.csv",
        "normalization/normalization_factors.csv",
        "normalization/normalization_factor_summary.png",
        "normalization/normalization_cutflow.csv",
        "normalization/normalization_cutflow.png",
        "normalization/component_fractions.csv",
        "normalization/normalization_energy_regions.png",
        "normalization/normalization_histograms.root",
        "normalization/step1b_mx2_ep/mx2_ep_summary.txt",
        "normalization/step1b_mx2_ep/mx2_ep_survival.csv",
        "normalization/step1b_mx2_ep/mx2_ep_unit_area_with_cut.png",
        "normalization/step1b_mx2_ep/mx2_ep_counts_with_cut.png",
        "normalization/step1c_mx2_eg/mx2_eg_summary.txt",
        "normalization/step1c_mx2_eg/mx2_eg_survival.csv",
        "normalization/step1c_mx2_eg/mx2_eg_unit_area_after_mx2_ep.png",
        "normalization/step1c_mx2_eg/mx2_eg_counts_after_mx2_ep.png"
    };
    for (const char* fn:expected) {
        const std::string full=out+"/"+fn;
        if (!gSystem->AccessPathName(full.c_str()))
            std::cout << "  " << full << "\n";
    } // endfor
    if (gSystem->AccessPathName((out+"/valerii_fd_results.csv").c_str())) {
        std::cerr << "WARNING: Valerii FD products were not completed. Check worker errors above.\n";
    }
}

