// photon_efficiency_valerii_reproduction.C
//
// Photon tag-and-probe analysis with two layers:
//   (1) stage-1 integrated diagnostics;
//   (2) a Valerii-style FD 7 x 3 x 6 efficiency/correction reproduction.
//
// The stage-1 layer remains useful QA; the Valerii layer is the physics path.
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
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
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

namespace pe {

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

    Double_t Mx_ep=0;
    Double_t probe_corr_p=0, probe_corr_theta=0, probe_corr_phi=0;

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
        Form("%s %s;#Delta p_{#gamma2}=p_{rec}-p_{miss} [GeV];Candidates",
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
                             ";missing-#gamma_{2} p [GeV];Candidates",80,0,8));
    s.denom_theta.reset(new TH1D(Form("h_%s_denom_theta",s.name.c_str()),
                                 ";missing-#gamma_{2} #theta [deg];Candidates",80,0,40));
    s.denom_phi.reset(new TH1D(Form("h_%s_denom_phi",s.name.c_str()),
                               ";missing-#gamma_{2} wrapped #phi [deg];Candidates",72,-30,330));
    s.ft_xy_projected.reset(new TH2D(Form("h_%s_ft_xy",s.name.c_str()),
                                     ";projected FT x [cm];projected FT y [cm]",
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
                               ";#Delta p [GeV];Candidates",DP_NBIN,DP_HMIN,DP_HMAX));
        wr.ft_h.reset(new TH1D(Form("h_%s_%s_ft",s.name.c_str(),mw.label),
                               ";#Delta p [GeV];Candidates",DP_NBIN,DP_HMIN,DP_HMAX));
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
            Form("%s;projected FT x [cm];projected FT y [cm]",sp->name.c_str()));
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
//   * saved base/generator event weight multiplied by Valerii normalization;
//   * 1-, 2-, and 3-sigma matching using the fit belonging to data or MC;
//   * epsilon_data, epsilon_MC, and epsilon_data/epsilon_MC correction maps;
//   * two published alternative normalization sets as a normalization variation.
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
    const char* label;
    double aao;
    double clasdis;
    double dvcs;
};
static const ValNormSet VAL_NORMS[] = {
    {"nominal",0.307,0.315,1.10},
    {"set1",   0.266,0.303,1.10},
    {"set2",   0.330,0.320,1.10}
};

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

struct ValComponent {
    std::string name;
    bool is_mc=false;
    long long entries=0;
    std::array<ValComponentBin,VAL_NBIN> bins;
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

    // Do NOT write one TTree row per reconstructed candidate.  With millions of
    // skim rows that temporary representation can become multi-GB and exhaust
    // the ifarm home/output filesystem.  The analysis only needs the residual
    // distributions and their weighted sums, so persist fixed-size histograms.
    std::array<std::unique_ptr<TH1D>,VAL_NBIN> hfit;
    std::array<std::unique_ptr<TH1D>,VAL_NBIN> hcount;
    for (int ib=0;ib<VAL_NBIN;ib++) {
        hfit[ib].reset(new TH1D(Form("fit_b%03d",ib),";#Delta p_{#gamma2} [GeV];weighted candidates",
                                VAL_DP_NBIN,VAL_DP_MIN,VAL_DP_MAX));
        hcount[ib].reset(new TH1D(Form("count_b%03d",ib),";#Delta p_{#gamma2} [GeV];weighted candidates",
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

    TDirectory* rd=f.mkdir("residuals");
    if (!rd) { f.Close(); return false; }
    rd->cd();
    for (int ib=0;ib<VAL_NBIN;ib++) {
        hfit[ib]->Write(Form("fit_b%03d",ib));
        hcount[ib]->Write(Form("count_b%03d",ib));
    } // endfor
    f.cd();
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
    std::unique_ptr<TH1D> h(new TH1D(name,";#Delta p_{#gamma2} [GeV];Weighted candidates",
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
        auto hm=val_make_mc_hist(vv,ib,norm,Form("h_val_mc_%s_b%03d",norm.label,ib));
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
               Form("%.2f < p < %.2f GeV;wrapped #phi [deg];#theta [deg]",
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

    std::vector<std::unique_ptr<TH1D>> hd_nom,hm_nom;
    auto nominal=evaluate_valerii_fd(vv,VAL_NORMS[0],&hd_nom,&hm_nom);
    auto set1=evaluate_valerii_fd(vv,VAL_NORMS[1]);
    auto set2=evaluate_valerii_fd(vv,VAL_NORMS[2]);

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
        const double sc=vp->is_mc?val_component_scale(vp->name,VAL_NORMS[0]):1.0;
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
            << "p edges [GeV]: 0.35 0.50 1.10 1.70 2.30 2.90 3.70 6.00\n"
            << "theta edges [deg]: 6 20 27 36\n"
            << "phi edges [deg]: -30 30 90 150 210 270 330\n"
            << "Residual histograms: 60 bins on [-1,1] GeV\n"
            << "Nominal MC normalization: 0.307*AAO + 0.315*CLASDIS(no-exclusivity) + 1.10*DVCS\n"
            << "Alternative normalization sets: (0.266,0.303,1.10), (0.330,0.320,1.10)\n"
            << "MC events are unit weighted inside each component; only the component normalization constants are applied.\n"
            << "The skim MC::Event.weight branch is NOT used in the Valerii FD calculation.\n"
            << "Tag photon: FD/PCAL only; probe coordinates: missing gamma2.\n"
            << "Stage-1 0.08<Mx(ep)<0.20 GeV enrichment is NOT applied in this reproduction.\n"
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
              << "MC nominal weights: 0.307 AAO + 0.315 CLASDIS + 1.10 DVCS\n"
              << "MC events are unit weighted inside each component.\n"
              << "Only the AAO/CLASDIS/DVCS component normalization constants are applied.\n"
              << "The skim MC::Event.weight branch is NOT used in this path.\n"
              << "Data/MC fits are independent in every p/theta/phi bin.\n"
              << "Primary correction convention: epsilon_data / epsilon_MC.\n"
              << "============================================================\n";
    auto vv=build_val_components_parallel(out);
    write_valerii_outputs(vv,out);
    gSystem->Exec(("rm -rf "+out+"/.valerii_workers").c_str());
}

} // namespace pe

void photon_efficiency_valerii_reproduction() {
    using namespace pe;

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

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
        "valerii_fd_histograms.root"
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

