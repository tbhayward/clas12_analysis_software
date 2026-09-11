// photon_efficiency_valerii_reproduction.C
//
// Stage-1 photon-efficiency analysis.
// Purposefully simple before moving to the full Valerii p x theta x phi maps:
//
//   * fully integrated FD extraction;
//   * fully integrated FT extraction;
//   * FT split into two broad energy bins: 0.4-2 GeV and >=2 GeV;
//   * 1-, 2-, and 3-sigma Delta-p matching;
//   * cut-flow accounting;
//   * pi0-parent-window stability scan;
//   * denominator p/theta/phi QA;
//   * automatic use of whatever AAOgen / CLASDIS / DVCSgen ROOT files exist;
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
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

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
static const int MIN_FIT_ENTRIES_FT = 20;

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

double event_weight(const Branches& b, bool weighted) {
    if (!weighted) return 1.0;
    return std::isfinite(b.mc_weight) ? b.mc_weight : 0.0;
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

    double denominator=0, denominator_w2=0;
    double numerator[3]={0,0,0};
    double numerator_w2[3]={0,0,0};

    long long denominator_rows=0;
    long long reconstructed_candidates=0;

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
    bool weighted=false;
    long long entries=0;

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

double efficiency(const RegionResult& r, int ns) {
    if (r.denominator<=0) return -1;
    return r.numerator[ns-1]/r.denominator;
}

double efficiency_error(const RegionResult& r, int ns, bool weighted) {
    if (r.denominator<=0) return 0;

    const double e=efficiency(r,ns);
    if (!(e>=0 && e<=1)) return 0;

    double n=r.denominator;
    if (weighted) {
        if (r.denominator_w2<=0) return 0;
        n=r.denominator*r.denominator/r.denominator_w2;
    }
    return n>0 ? std::sqrt(std::max(0.0,e*(1.0-e)/n)) : 0;
}

// -----------------------------------------------------------------------------
// Fill helpers.
// -----------------------------------------------------------------------------
void fill_region_residual(RegionResult& r, const Branches& b, double w) {
    if (!in_region(b,r)) return;

    r.denominator_rows++;
    const int k=best_probe_candidate(b,r.detector);
    if (k<0) return;

    r.reconstructed_candidates++;
    const double dp=b.neutral_p[k]-b.probe_corr_p;
    if (std::isfinite(dp)) r.residual->Fill(dp,w);

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

void fill_region_counts(RegionResult& r, const Branches& b, double w) {
    if (!in_region(b,r)) return;

    r.denominator+=w;
    r.denominator_w2+=w*w;

    const int k=best_probe_candidate(b,r.detector);
    if (!r.fit.valid || k<0) return;

    const double dp=b.neutral_p[k]-b.probe_corr_p;
    const double d=std::fabs(dp-r.fit.mean);

    for (int ns=1;ns<=3;ns++) {
        if (d<ns*r.fit.sigma) {
            r.numerator[ns-1]+=w;
            r.numerator_w2[ns-1]+=w*w;
        }
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
        const double w=event_weight(b,s.weighted);
        for (auto& wr : s.windows) {
            if (!finite_good(b.Mx_ep) ||
                b.Mx_ep<wr.w.lo || b.Mx_ep>wr.w.hi) continue;

            if (in_fd(b)) {
                wr.fd_denom++;
                const int k=best_probe_candidate(b,1);
                if (k>=0) {
                    wr.fd_reco++;
                    wr.fd_h->Fill(b.neutral_p[k]-b.probe_corr_p,w);
                }
            }

            if (in_ft(b)) {
                wr.ft_denom++;
                const int k=best_probe_candidate(b,0);
                if (k>=0) {
                    wr.ft_reco++;
                    wr.ft_h->Fill(b.neutral_p[k]-b.probe_corr_p,w);
                }
            }
        } // endfor

        if (!pass_nominal_mass(b)) continue;
        s.cutflow.nominal_mass++;

        s.denom_p->Fill(b.probe_corr_p,w);
        s.denom_theta->Fill(b.probe_corr_theta,w);
        s.denom_phi->Fill(wrap_phi(b.probe_corr_phi),w);

        if (in_fd(b)) s.cutflow.fd++;
        if (in_ft(b)) {
            s.cutflow.ft++;

            const FTProjection q=project_ft(b,s.ft_plane);
            if (q.valid) {
                s.cutflow.ft_projectable++;
                s.ft_xy_projected->Fill(q.x,q.y,w);
                if (q.fiducial) s.cutflow.ft_projected_fid++;
            }
        }

        fill_region_residual(s.fd,b,w);
        fill_region_residual(s.ft_all,b,w);
        fill_region_residual(s.ft_low,b,w);
        fill_region_residual(s.ft_high,b,w);
    } // endfor

    s.fd.fit=fit_residual(s.fd.residual.get(),MIN_FIT_ENTRIES_FD);
    s.ft_all.fit=fit_residual(s.ft_all.residual.get(),MIN_FIT_ENTRIES_FT);
    s.ft_low.fit=fit_residual(s.ft_low.residual.get(),MIN_FIT_ENTRIES_FT);
    s.ft_high.fit=fit_residual(s.ft_high.residual.get(),MIN_FIT_ENTRIES_FT);

    for (auto& wr : s.windows) {
        wr.fd_fit=fit_residual(wr.fd_h.get(),MIN_FIT_ENTRIES_FD);
        wr.ft_fit=fit_residual(wr.ft_h.get(),MIN_FIT_ENTRIES_FT);
    } // endfor

    // Pass 2: 1/2/3 sigma counts.
    for (Long64_t i=0;i<c.GetEntries();i++) {
        c.GetEntry(i);
        if (!pass_before_mass(b)) continue;
        if (!pass_nominal_mass(b)) continue;

        const double w=event_weight(b,s.weighted);

        fill_region_counts(s.fd,b,w);
        fill_region_counts(s.ft_all,b,w);
        fill_region_counts(s.ft_low,b,w);
        fill_region_counts(s.ft_high,b,w);
    } // endfor

    return true;
}

// -----------------------------------------------------------------------------
// Output.
// -----------------------------------------------------------------------------
void save_residual_plot(const SampleResult& s,
                        const RegionResult& r,
                        const std::string& out) {
    TCanvas c(Form("c_%s_%s",s.name.c_str(),r.name.c_str()),"",950,720);
    r.residual->SetLineWidth(2);
    r.residual->Draw("E");

    const double ymax=std::max(1.0,r.residual->GetMaximum());

    std::unique_ptr<TF1> fit_draw;
    if (r.fit.valid) {
        // Reconstruct the exact fitted gaus(0)+pol1(3) function from the
        // parameters saved by fit_residual(), and overlay it on the data.
        fit_draw.reset(new TF1(
            Form("fit_draw_%s_%s",s.name.c_str(),r.name.c_str()),
            "gaus(0)+pol1(3)",
            r.fit.fit_lo,r.fit.fit_hi));
        fit_draw->SetParameters(
            r.fit.amplitude,
            r.fit.mean,
            r.fit.sigma,
            r.fit.bg0,
            r.fit.bg1);
        fit_draw->SetLineWidth(3);
        fit_draw->Draw("SAME");

        for (int ns=1;ns<=3;ns++) {
            TLine l1(r.fit.mean-ns*r.fit.sigma,0,
                     r.fit.mean-ns*r.fit.sigma,0.92*ymax);
            TLine l2(r.fit.mean+ns*r.fit.sigma,0,
                     r.fit.mean+ns*r.fit.sigma,0.92*ymax);
            l1.SetLineStyle(ns);
            l2.SetLineStyle(ns);
            l1.Draw();
            l2.Draw();
        } // endfor
    }

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.033);
    tx.DrawLatex(0.13,0.88,Form("%s %s",s.name.c_str(),r.name.c_str()));
    tx.DrawLatex(0.13,0.83,Form("denom = %lld, reconstructed = %lld",
                                r.denominator_rows,r.reconstructed_candidates));

    if (r.fit.valid) {
        tx.DrawLatex(0.13,0.78,
            Form("#mu = %.4f #pm %.4f GeV",r.fit.mean,r.fit.mean_err));
        tx.DrawLatex(0.13,0.73,
            Form("#sigma = %.4f #pm %.4f GeV",r.fit.sigma,r.fit.sigma_err));
        tx.DrawLatex(0.13,0.68,
            Form("#chi^{2}/ndf = %.1f/%d",r.fit.chi2,r.fit.ndf));
    } else {
        tx.DrawLatex(0.13,0.78,Form("FIT INVALID: %s",r.fit.reason.c_str()));
    }

    c.SaveAs((out+"/"+s.name+"_"+r.name+"_delta_p.png").c_str());
}

void save_kinematics(const SampleResult& s, const std::string& out) {
    TCanvas c(Form("c_%s_kin",s.name.c_str()),"",1500,450);
    c.Divide(3,1);
    c.cd(1); s.denom_p->Draw("E");
    c.cd(2); s.denom_theta->Draw("E");
    c.cd(3); s.denom_phi->Draw("E");
    c.SaveAs((out+"/"+s.name+"_denominator_kinematics.png").c_str());
}

void save_ft_projection(const SampleResult& s, const std::string& out) {
    if (!s.ft_plane.valid) return;

    TCanvas c(Form("c_%s_ft_xy",s.name.c_str()),"",800,750);
    s.ft_xy_projected->Draw("COLZ");
    c.SaveAs((out+"/"+s.name+"_FT_projected_xy.png").c_str());
}

void write_cutflow(const SampleResult& s, const std::string& out) {
    std::ofstream f(out+"/"+s.name+"_cutflow.txt");
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
        f << "\nInferred FT response-plane z = " << s.ft_plane.z
          << " cm from " << s.ft_plane.n << " reconstructed FT responses.\n";

    f << "\nNOTE: FT projected fiducial is diagnostic only in this stage;\n"
      << "the extraction denominator still uses the angular FT region.\n";
}

void write_window_scan(const SampleResult& s, const std::string& out) {
    std::ofstream f(out+"/"+s.name+"_parent_window_scan.csv");
    f << "window,mmin,mmax,fd_denom,fd_reco,fd_fit_valid,fd_mean,fd_sigma,"
      << "ft_denom,ft_reco,ft_fit_valid,ft_mean,ft_sigma\n";

    for (const auto& wr : s.windows) {
        f << wr.w.label << "," << wr.w.lo << "," << wr.w.hi << ","
          << wr.fd_denom << "," << wr.fd_reco << ","
          << (wr.fd_fit.valid?1:0) << ","
          << wr.fd_fit.mean << "," << wr.fd_fit.sigma << ","
          << wr.ft_denom << "," << wr.ft_reco << ","
          << (wr.ft_fit.valid?1:0) << ","
          << wr.ft_fit.mean << "," << wr.ft_fit.sigma << "\n";
    } // endfor
}

void write_region_csv_header(std::ofstream& f) {
    f << "sample,region,pmin,pmax,thetamin,thetamax,"
      << "fit_valid,fit_reason,fit_candidates,mean,mean_err,sigma,sigma_err,"
      << "chi2,ndf,denominator,reconstructed,"
      << "num_1sigma,num_2sigma,num_3sigma,"
      << "eff_1sigma,eff_2sigma,eff_3sigma,"
      << "err_1sigma,err_2sigma,err_3sigma,"
      << "truth_pi0_probe,selected_reco_with_truth,"
      << "reco_truth_da_lt1,reco_truth_da_lt2,reco_truth_da_lt3\n";
}

void append_region_csv(std::ofstream& f,
                       const SampleResult& s,
                       const RegionResult& r) {
    f << s.name << "," << r.name << ","
      << r.pmin << "," << r.pmax << ","
      << r.thmin << "," << r.thmax << ","
      << (r.fit.valid?1:0) << ",\"" << r.fit.reason << "\","
      << r.fit.candidates << ","
      << r.fit.mean << "," << r.fit.mean_err << ","
      << r.fit.sigma << "," << r.fit.sigma_err << ","
      << r.fit.chi2 << "," << r.fit.ndf << ","
      << r.denominator << "," << r.reconstructed_candidates << ","
      << r.numerator[0] << "," << r.numerator[1] << "," << r.numerator[2] << ","
      << efficiency(r,1) << "," << efficiency(r,2) << "," << efficiency(r,3) << ","
      << efficiency_error(r,1,s.weighted) << ","
      << efficiency_error(r,2,s.weighted) << ","
      << efficiency_error(r,3,s.weighted) << ","
      << r.truth_pi0_probe << "," << r.selected_reco_with_truth << ","
      << r.reco_truth_da_lt1 << "," << r.reco_truth_da_lt2 << ","
      << r.reco_truth_da_lt3 << "\n";
}

void write_sample_summary(const SampleResult& s, const std::string& out) {
    std::ofstream f(out+"/"+s.name+"_summary.txt");

    f << "Sample: " << s.name << "\n"
      << "Input: " << s.input << "\n"
      << "Tree entries: " << s.entries << "\n"
      << "Nominal parent window: " << PI0_M_MIN
      << " < Mx(ep) < " << PI0_M_MAX << " GeV\n\n";

    const RegionResult* rr[] = {&s.fd,&s.ft_all,&s.ft_low,&s.ft_high};
    for (const RegionResult* r : rr) {
        f << "--- " << r->name << " ---\n"
          << "denominator rows: " << r->denominator_rows << "\n"
          << "reconstructed candidates: " << r->reconstructed_candidates << "\n"
          << "fit valid: " << (r->fit.valid?1:0)
          << " (" << r->fit.reason << ")\n"
          << "mean: " << r->fit.mean << " +/- " << r->fit.mean_err << " GeV\n"
          << "sigma: " << r->fit.sigma << " +/- " << r->fit.sigma_err << " GeV\n";

        for (int ns=1;ns<=3;ns++)
            f << ns << " sigma matched fraction: "
              << efficiency(*r,ns) << " +/- "
              << efficiency_error(*r,ns,s.weighted) << "\n";

        if (r->truth_pi0_probe>0) {
            f << "MC truth pi0-probe rows: " << r->truth_pi0_probe << "\n"
              << "selected reco candidates with truth comparison: "
              << r->selected_reco_with_truth << "\n"
              << "reco-to-truth dAlpha <1/<2/<3 deg: "
              << r->reco_truth_da_lt1 << " / "
              << r->reco_truth_da_lt2 << " / "
              << r->reco_truth_da_lt3 << "\n";
        }
        f << "\n";
    } // endfor
}

void write_root(const std::vector<std::unique_ptr<SampleResult>>& samples,
                const std::string& out) {
    TFile f((out+"/analysis_histograms.root").c_str(),"RECREATE");

    for (const auto& sp : samples) {
        f.mkdir(sp->name.c_str());
        f.cd(sp->name.c_str());

        sp->fd.residual->Write();
        sp->ft_all.residual->Write();
        sp->ft_low.residual->Write();
        sp->ft_high.residual->Write();

        sp->denom_p->Write();
        sp->denom_theta->Write();
        sp->denom_phi->Write();
        sp->ft_xy_projected->Write();

        for (const auto& wr : sp->windows) {
            wr.fd_h->Write();
            wr.ft_h->Write();
        } // endfor
        f.cd();
    } // endfor

    f.Close();
}

void process_if_available(const std::string& name,
                          const std::string& dir,
                          bool weighted,
                          std::vector<std::unique_ptr<SampleResult>>& samples) {
    if (!has_root_files(dir)) {
        std::cout << "[AUTO] No " << name << " ROOT files yet; skipping.\n";
        return;
    }

    std::cout << "[AUTO] Found " << name << " ROOT files in " << dir << "\n";

    std::unique_ptr<SampleResult> s(new SampleResult);
    s->name=name;
    s->input=dir;
    s->weighted=weighted;

    if (!analyze_sample(*s)) {
        std::cerr << "WARNING: failed to analyze " << name << "\n";
        return;
    }
    samples.push_back(std::move(s));
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
        << " Photon-efficiency stage-1 analysis\n"
        << "============================================================\n"
        << "FD: fully integrated\n"
        << "FT: integrated + 0.4-2 GeV + >=2 GeV\n"
        << "Matching: Delta p_gamma2 = p_rec - p_miss\n"
        << "Outputs are rebuilt in ./output each invocation.\n"
        << "============================================================\n";

    std::vector<std::unique_ptr<SampleResult>> samples;

    process_if_available("data",DATA_DIR,false,samples);
    process_if_available("aaogen",AAOGEN_DIR,true,samples);
    process_if_available("clasdis",CLASDIS_DIR,true,samples);
    process_if_available("dvcsgen",DVCSGEN_DIR,true,samples);

    std::ofstream csv(out+"/integrated_results.csv");
    write_region_csv_header(csv);

    for (const auto& sp : samples) {
        save_residual_plot(*sp,sp->fd,out);
        save_residual_plot(*sp,sp->ft_all,out);
        save_residual_plot(*sp,sp->ft_low,out);
        save_residual_plot(*sp,sp->ft_high,out);

        save_kinematics(*sp,out);
        save_ft_projection(*sp,out);

        write_cutflow(*sp,out);
        write_window_scan(*sp,out);
        write_sample_summary(*sp,out);

        append_region_csv(csv,*sp,sp->fd);
        append_region_csv(csv,*sp,sp->ft_all);
        append_region_csv(csv,*sp,sp->ft_low);
        append_region_csv(csv,*sp,sp->ft_high);
    } // endfor

    csv.close();
    write_root(samples,out);

    std::cout
        << "\nFinished. Main products:\n"
        << "  output/integrated_results.csv\n"
        << "  output/data_cutflow.txt\n"
        << "  output/data_parent_window_scan.csv\n"
        << "  output/data_FD_delta_p.png\n"
        << "  output/data_FT_all_delta_p.png\n"
        << "  output/data_FT_Elt2_delta_p.png\n"
        << "  output/data_FT_Ege2_delta_p.png\n"
        << "  output/data_denominator_kinematics.png\n"
        << "  output/data_FT_projected_xy.png\n"
        << "  output/analysis_histograms.root\n";
}
