// photon_efficiency_valerii_reproduction.C
//
// Stage-1 photon tag-and-probe diagnostic analysis.
// Purposefully simple before moving to the full Valerii p x theta x phi maps:
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
#include <TParameter.h>
#include <TNamed.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
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
// We calculate a projected-face diagnostic but do NOT yet make it mandatory;
// that lets us validate the projection before changing the denominator.
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
    f.SetParLimits(1,-0.50,0.50);
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
void fill_region_residual(RegionResult& r, const Branches& b) {
    if (!in_region(b,r)) return;

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

            if (in_ft(b)) {
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
        if (in_ft(b)) {
            s.cutflow.ft++;

            const FTProjection q=project_ft(b,s.ft_plane);
            if (q.valid) {
                s.cutflow.ft_projectable++;
                s.ft_xy_projected->Fill(q.x,q.y);
                if (q.fiducial) s.cutflow.ft_projected_fid++;
            }
        }

        fill_region_residual(s.fd,b);
        fill_region_residual(s.ft_all,b);
        fill_region_residual(s.ft_low,b);
        fill_region_residual(s.ft_high,b);
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
        f << "FT projected fiducial remains diagnostic only at this stage.\n\n";

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
        << "FT: integrated + 0.4-2 GeV + >=2 GeV\n"
        << "Matching: Delta p_gamma2 = p_rec - p_miss\n"
        << "Primary numbers: unweighted RAW RECOVERY FRACTIONS (not efficiencies)\n"
        << "MC::Event.weight: QA only, not applied\n"
        << "Execution: independent samples analyzed in parallel child processes\n"
        << "Outputs are rebuilt in ./output each invocation.\n"
        << "============================================================\n";

    std::vector<std::unique_ptr<SampleResult>> samples=analyze_samples_parallel(out);

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

    std::cout
        << "\nFinished. Compact output products:\n"
        << "  output/analysis_summary.txt\n"
        << "  output/integrated_results.csv\n"
        << "  output/parent_window_scan.csv\n"
        << "  output/delta_p_residuals.png\n"
        << "  output/denominator_kinematics.png\n"
        << "  output/FT_projected_xy.png\n"
        << "  output/analysis_histograms.root\n";
}

