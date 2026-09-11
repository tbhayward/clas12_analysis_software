// photon_efficiency_valerii_reproduction.C
//
// First-pass reproduction of the Valerii missing-gamma2 photon-efficiency
// analysis using the PhotonEfficiency skim.
//
// FD:
//   * Uses Valerii's missing-gamma2 binning from the supplied note:
//       p [GeV]   : 0.35, 0.50, 1.10, 1.70, 2.30, 2.90, 3.90, 6.00
//       theta [deg]: 6, 20, 27, 36
//       wrapped phi: -30, 30, 90, 150, 210, 270, 330
//   * Fits Delta p_gamma2 = p_rec - p_miss separately in every p/theta/phi bin.
//   * Reuses fitted mean/sigma for 1-, 2-, and 3-sigma matching.
//   * Produces denominator/numerator/efficiency CSV products and QA plots.
//
// FT:
//   * Built alongside FD with deliberately coarser p/theta binning.
//   * Uses the same Delta-p matching logic.
//   * IMPORTANT: the current skim does not yet contain enough information to
//     impose the exact second-photon beta/fiducial selection, nor a projected
//     missing-probe FT face fiducial.  FT results are therefore preliminary.
//
// The current skim *does* contain the tag beta/fiducial flags, and those are
// imposed here.  For the reconstructed probe, the closest saved PID-22 neutral
// in the required detector with p>=0.4 GeV is used.
//
// Usage (data only):
// root -l -b -q 'photon_efficiency_valerii_reproduction.C(
//   "/work/.../data/fa18_inb",
//   "output/fa18_inb")'
//
// Usage once MC skims exist:
// root -l -b -q 'photon_efficiency_valerii_reproduction.C(
//   "/work/.../data/fa18_inb",
//   "output/fa18_inb",
//   "/work/.../aaogen/fa18_inb",
//   "/work/.../clasdis/fa18_inb",
//   "/work/.../dvcsgen/fa18_inb")'
//
// The first pass intentionally does NOT combine the three MC components with
// guessed external normalization constants.  It produces each component
// separately.  Once the exact Valerii normalization constants/configuration
// are carried over, the weighted-total MC can be formed without changing the
// event/matching machinery.

#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
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
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace per {

static const double PI = TMath::Pi();

// -----------------------------------------------------------------------------
// Analysis knobs.
// -----------------------------------------------------------------------------

// Broad pi0 parent window for this first reproduction pass.
// Keep configurable here rather than burying it in event logic.
// The exact production value should ultimately be copied verbatim from
// Valerii's selection source once that source/configuration is available.
static const double PI0_M_MIN = 0.08;
static const double PI0_M_MAX = 0.20;

// Reconstructed probe threshold used in the skim / CLAS12 photon selection.
static const double PROBE_P_MIN = 0.40;

// Delta-p histogram and fit range, matching the [-1,1] GeV range shown in
// Valerii's note.
static const double DP_HMIN = -1.0;
static const double DP_HMAX =  1.0;
static const int    DP_NBIN = 80;

// Fit only if enough reconstructed-probe candidates are available.
static const int MIN_FIT_ENTRIES_FD = 30;
static const int MIN_FIT_ENTRIES_FT = 20;

// Start the Gaussian+linear-background fit around the visible peak.
// The second pass can make this iterative if needed.
static const double FIT_MIN = -0.60;
static const double FIT_MAX =  0.60;

// Valerii FD binning.
static const std::vector<double> FD_P_EDGES =
    {0.35, 0.50, 1.10, 1.70, 2.30, 2.90, 3.90, 6.00};
static const std::vector<double> FD_THETA_EDGES =
    {6.0, 20.0, 27.0, 36.0};
static const std::vector<double> PHI_EDGES =
    {-30.0, 30.0, 90.0, 150.0, 210.0, 270.0, 330.0};

// Coarse FT bins for the first extension.  We intentionally do not explode
// 788-ish pi0-enriched candidates into a 7x3x6 grid.
static const std::vector<double> FT_P_EDGES =
    {0.40, 1.00, 2.00, 3.00, 4.50, 8.50};
static const std::vector<double> FT_THETA_EDGES =
    {2.0, 3.5, 5.5};

// -----------------------------------------------------------------------------
// Helpers.
// -----------------------------------------------------------------------------

double wrap_phi(double phi_deg) {
    while (phi_deg < -30.0) phi_deg += 360.0;
    while (phi_deg >= 330.0) phi_deg -= 360.0;
    return phi_deg;
}

int find_bin(double x, const std::vector<double>& edges) {
    if (edges.size() < 2) return -1;
    if (x < edges.front() || x > edges.back()) return -1;
    if (x == edges.back()) return static_cast<int>(edges.size()) - 2;
    for (size_t i=0; i+1<edges.size(); ++i) {
        if (x >= edges[i] && x < edges[i+1]) return static_cast<int>(i);
    }
    return -1;
}

std::string make_pattern(const char* input) {
    TString in(input ? input : "");
    if (in.Length()==0) return "";
    void* d = gSystem->OpenDirectory(in.Data());
    if (d) {
        gSystem->FreeDirectory(d);
        if (!in.EndsWith("/")) in += "/";
        in += "*.root";
    }
    return in.Data();
}

void ensure_dir(const std::string& d) {
    if (!d.empty()) gSystem->mkdir(d.c_str(), true);
}

bool finite_good(double x) {
    return std::isfinite(x) && x > -900.0;
}

std::string fmt_edge(double x) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(2) << x;
    return s.str();
}

struct BinDef {
    int ip=-1, it=-1, iph=-1;
    double p0=0,p1=0,t0=0,t1=0,ph0=-999,ph1=-999;

    std::string key(const std::string& detector) const {
        std::ostringstream s;
        s << detector << "_p" << ip << "_t" << it;
        if (iph >= 0) s << "_ph" << iph;
        return s.str();
    }

    std::string label() const {
        std::ostringstream s;
        s << std::fixed << std::setprecision(2)
          << p0 << "<p<" << p1 << " GeV, "
          << std::setprecision(1) << t0 << "<theta<" << t1 << " deg";
        if (iph >= 0) s << ", " << ph0 << "<phi<" << ph1 << " deg";
        return s.str();
    }
};

struct FitResult {
    bool valid=false;
    double mean=0;
    double mean_err=0;
    double sigma=0;
    double sigma_err=0;
    double chi2=0;
    int ndf=0;
    long long candidates=0;
};

struct CountResult {
    double denom=0;
    double denom_w2=0;
    double num[3] = {0,0,0};
    double num_w2[3] = {0,0,0};
};

struct Branches {
    Int_t runnum=0, evnum=0, is_mc=0;
    Int_t p_pass_standard=0;
    Int_t tag_detector=-1, tag_pass_beta=0, tag_pass_fiducial=0;
    Double_t mc_weight=1.0;
    Double_t Mx_ep=0, Mx2_epg_raw=0;
    Double_t probe_corr_p=0, probe_corr_theta=0, probe_corr_phi=0;
    Int_t neutral_idx[5];
    Int_t neutral_pid[5];
    Int_t neutral_charge[5];
    Int_t neutral_detector[5];
    Double_t neutral_p[5];
    Double_t neutral_theta[5];
    Double_t neutral_phi[5];
    Double_t neutral_delta_alpha[5];

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
        }
    }
};

bool attach(TChain& c, Branches& b) {
    const char* required[] = {
        "runnum","evnum","is_mc","mc_weight","p_pass_standard",
        "tag_detector","tag_pass_beta","tag_pass_fiducial",
        "Mx_ep","Mx2_epg_raw","probe_corr_p","probe_corr_theta","probe_corr_phi",
        "neutral_idx","neutral_pid","neutral_charge","neutral_detector",
        "neutral_p","neutral_theta","neutral_phi","neutral_delta_alpha"
    };
    bool ok=true;
    for (const char* n : required) {
        if (!c.GetBranch(n)) {
            std::cerr << "ERROR: missing branch " << n << "\n";
            ok=false;
        }
    }
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
    c.SetBranchAddress("Mx2_epg_raw",&b.Mx2_epg_raw);
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
    return true;
}

// Find the closest saved reconstructed PID-22 neutral in the requested detector.
// The skim's neutral array is already ordered by opening angle relative to the
// missing-probe direction, but scan all five so a non-photon neutral at index 0
// does not hide the nearest photon.
int best_probe_candidate(const Branches& b, int detector) {
    int best=-1;
    double best_da=1e9;
    for (int k=0;k<5;k++) {
        if (b.neutral_idx[k] < 0) continue;
        if (b.neutral_charge[k] != 0) continue;
        if (b.neutral_pid[k] != 22) continue;
        if (b.neutral_detector[k] != detector) continue;
        if (!finite_good(b.neutral_p[k]) || b.neutral_p[k] < PROBE_P_MIN) continue;
        if (!finite_good(b.neutral_delta_alpha[k])) continue;
        if (b.neutral_delta_alpha[k] < best_da) {
            best_da=b.neutral_delta_alpha[k];
            best=k;
        }
    }
    return best;
}

bool common_denominator_selection(const Branches& b) {
    if (!b.p_pass_standard) return false;
    if (!b.tag_pass_beta || !b.tag_pass_fiducial) return false;

    // A usable tag can be either FD or FT in the generalized study.
    if (!(b.tag_detector==0 || b.tag_detector==1)) return false;

    if (!finite_good(b.Mx_ep) || b.Mx_ep < PI0_M_MIN || b.Mx_ep > PI0_M_MAX)
        return false;

    if (!finite_good(b.probe_corr_p) || b.probe_corr_p < PROBE_P_MIN)
        return false;

    if (!finite_good(b.probe_corr_theta) || !finite_good(b.probe_corr_phi))
        return false;

    return true;
}

int fd_bin(const Branches& b, BinDef& out) {
    const int ip=find_bin(b.probe_corr_p,FD_P_EDGES);
    const int it=find_bin(b.probe_corr_theta,FD_THETA_EDGES);
    const double ph=wrap_phi(b.probe_corr_phi);
    const int iph=find_bin(ph,PHI_EDGES);
    if (ip<0 || it<0 || iph<0) return -1;

    out.ip=ip; out.it=it; out.iph=iph;
    out.p0=FD_P_EDGES[ip]; out.p1=FD_P_EDGES[ip+1];
    out.t0=FD_THETA_EDGES[it]; out.t1=FD_THETA_EDGES[it+1];
    out.ph0=PHI_EDGES[iph]; out.ph1=PHI_EDGES[iph+1];

    return (ip*3 + it)*6 + iph;
}

int ft_bin(const Branches& b, BinDef& out) {
    const int ip=find_bin(b.probe_corr_p,FT_P_EDGES);
    const int it=find_bin(b.probe_corr_theta,FT_THETA_EDGES);
    if (ip<0 || it<0) return -1;

    out.ip=ip; out.it=it; out.iph=-1;
    out.p0=FT_P_EDGES[ip]; out.p1=FT_P_EDGES[ip+1];
    out.t0=FT_THETA_EDGES[it]; out.t1=FT_THETA_EDGES[it+1];

    return ip*2 + it;
}

double event_weight(const Branches& b, bool weighted) {
    if (!weighted) return 1.0;
    if (!std::isfinite(b.mc_weight)) return 0.0;
    return b.mc_weight;
}

FitResult fit_residual(TH1D* h, int min_entries) {
    FitResult r;
    if (!h) return r;
    r.candidates = static_cast<long long>(h->GetEntries());
    if (h->GetEntries() < min_entries) return r;

    const int maxbin=h->GetMaximumBin();
    double peak=h->GetBinCenter(maxbin);

    // First narrow Gaussian seed around the observed peak.
    double seed_lo=std::max(FIT_MIN,peak-0.20);
    double seed_hi=std::min(FIT_MAX,peak+0.20);
    TF1 gseed("gseed_tmp","gaus",seed_lo,seed_hi);
    gseed.SetParameters(h->GetMaximum(),peak,0.08);
    int s0=h->Fit(&gseed,"QNR");
    double mu = (s0==0 ? gseed.GetParameter(1) : peak);
    double sg = (s0==0 ? std::fabs(gseed.GetParameter(2)) : 0.10);
    if (!(sg>0.005 && sg<0.5)) sg=0.10;

    double flo=std::max(FIT_MIN,mu-4.0*sg);
    double fhi=std::min(FIT_MAX,mu+4.0*sg);
    if (fhi-flo < 0.20) { flo=FIT_MIN; fhi=FIT_MAX; }

    // Gaussian signal + linear background.  This is deliberately simple for
    // the first pass and mirrors the signal+background structure visible in
    // Valerii's representative residual fits.
    TF1 f("fit_tmp","gaus(0)+pol1(3)",flo,fhi);
    f.SetParameters(h->GetMaximum(),mu,sg,
                    std::max(0.0,h->GetBinContent(1)),0.0);
    f.SetParLimits(2,0.005,0.50);
    int status=h->Fit(&f,"QNR");

    if (status!=0) return r;
    double sigma=std::fabs(f.GetParameter(2));
    if (!(sigma>0.005 && sigma<0.50)) return r;

    r.valid=true;
    r.mean=f.GetParameter(1);
    r.mean_err=f.GetParError(1);
    r.sigma=sigma;
    r.sigma_err=f.GetParError(2);
    r.chi2=f.GetChisquare();
    r.ndf=f.GetNDF();
    return r;
}

struct SampleResult {
    std::string name;
    std::string input;
    bool weighted=false;

    std::vector<std::unique_ptr<TH1D>> fd_h;
    std::vector<std::unique_ptr<TH1D>> ft_h;
    std::vector<FitResult> fd_fit;
    std::vector<FitResult> ft_fit;
    std::vector<CountResult> fd_count;
    std::vector<CountResult> ft_count;

    std::unique_ptr<TH1D> fd_global;
    std::unique_ptr<TH1D> ft_global;

    long long entries=0;
    long long denom_common=0;
    long long fd_denom_rows=0;
    long long ft_denom_rows=0;
    long long fd_probe_candidates=0;
    long long ft_probe_candidates=0;
};

void initialize_sample(SampleResult& r) {
    const int nfd=(FD_P_EDGES.size()-1)*(FD_THETA_EDGES.size()-1)*(PHI_EDGES.size()-1);
    const int nft=(FT_P_EDGES.size()-1)*(FT_THETA_EDGES.size()-1);

    r.fd_h.resize(nfd);
    r.ft_h.resize(nft);
    r.fd_fit.resize(nfd);
    r.ft_fit.resize(nft);
    r.fd_count.resize(nfd);
    r.ft_count.resize(nft);

    for (int i=0;i<nfd;i++) {
        r.fd_h[i].reset(new TH1D(Form("h_%s_fd_%03d",r.name.c_str(),i),
                                  ";#Delta p_{#gamma2}=p_{rec}-p_{miss} [GeV];Candidates",
                                  DP_NBIN,DP_HMIN,DP_HMAX));
        r.fd_h[i]->Sumw2();
        r.fd_h[i]->SetDirectory(nullptr);
    }

    for (int i=0;i<nft;i++) {
        r.ft_h[i].reset(new TH1D(Form("h_%s_ft_%03d",r.name.c_str(),i),
                                  ";#Delta p_{#gamma2}=p_{rec}-p_{miss} [GeV];Candidates",
                                  DP_NBIN,DP_HMIN,DP_HMAX));
        r.ft_h[i]->Sumw2();
        r.ft_h[i]->SetDirectory(nullptr);
    }

    r.fd_global.reset(new TH1D(Form("h_%s_fd_global",r.name.c_str()),
                               ";#Delta p_{#gamma2} [GeV];Candidates",
                               DP_NBIN,DP_HMIN,DP_HMAX));
    r.ft_global.reset(new TH1D(Form("h_%s_ft_global",r.name.c_str()),
                               ";#Delta p_{#gamma2} [GeV];Candidates",
                               DP_NBIN,DP_HMIN,DP_HMAX));
    r.fd_global->Sumw2();
    r.ft_global->Sumw2();
    r.fd_global->SetDirectory(nullptr);
    r.ft_global->SetDirectory(nullptr);
}

bool load_sample(SampleResult& r) {
    std::string pat=make_pattern(r.input.c_str());
    if (pat.empty()) return false;

    TChain c("PhotonEfficiency");
    int nf=c.Add(pat.c_str());
    r.entries=c.GetEntries();
    std::cout << "\n[" << r.name << "] files=" << nf
              << " entries=" << r.entries
              << " pattern=" << pat << "\n";
    if (nf<=0 || r.entries<=0) return false;

    Branches b;
    b.reset_arrays();
    if (!attach(c,b)) return false;

    initialize_sample(r);

    // ---------------------------------------------------------------------
    // Pass 1: build reconstructed-minus-missing momentum residuals.
    // ---------------------------------------------------------------------
    for (Long64_t i=0;i<c.GetEntries();i++) {
        c.GetEntry(i);
        if (!common_denominator_selection(b)) continue;
        r.denom_common++;

        const double w=event_weight(b,r.weighted);
        if (w==0.0) continue;

        BinDef bd;
        int ib=fd_bin(b,bd);
        if (ib>=0) {
            r.fd_denom_rows++;
            int k=best_probe_candidate(b,1);
            if (k>=0) {
                r.fd_probe_candidates++;
                const double dp=b.neutral_p[k]-b.probe_corr_p;
                if (std::isfinite(dp)) {
                    r.fd_h[ib]->Fill(dp,w);
                    r.fd_global->Fill(dp,w);
                }
            }
        }

        ib=ft_bin(b,bd);
        if (ib>=0) {
            r.ft_denom_rows++;
            int k=best_probe_candidate(b,0);
            if (k>=0) {
                r.ft_probe_candidates++;
                const double dp=b.neutral_p[k]-b.probe_corr_p;
                if (std::isfinite(dp)) {
                    r.ft_h[ib]->Fill(dp,w);
                    r.ft_global->Fill(dp,w);
                }
            }
        }
    } // endfor

    for (size_t i=0;i<r.fd_h.size();i++)
        r.fd_fit[i]=fit_residual(r.fd_h[i].get(),MIN_FIT_ENTRIES_FD);
    for (size_t i=0;i<r.ft_h.size();i++)
        r.ft_fit[i]=fit_residual(r.ft_h[i].get(),MIN_FIT_ENTRIES_FT);

    // ---------------------------------------------------------------------
    // Pass 2: denominator and matched numerator counts.
    //
    // Exactly one missing-probe hypothesis contributes to the denominator.
    // A numerator entry is counted if the reconstructed probe candidate has
    // Delta-p inside N*sigma of that sample/bin's fitted peak.
    // ---------------------------------------------------------------------
    for (Long64_t i=0;i<c.GetEntries();i++) {
        c.GetEntry(i);
        if (!common_denominator_selection(b)) continue;

        const double w=event_weight(b,r.weighted);
        if (w==0.0) continue;

        BinDef bd;
        int ib=fd_bin(b,bd);
        if (ib>=0) {
            CountResult& cr=r.fd_count[ib];
            cr.denom += w;
            cr.denom_w2 += w*w;

            const FitResult& fr=r.fd_fit[ib];
            int k=best_probe_candidate(b,1);
            if (fr.valid && k>=0) {
                const double dp=b.neutral_p[k]-b.probe_corr_p;
                const double d=std::fabs(dp-fr.mean);
                for (int ns=1;ns<=3;ns++) {
                    if (d < ns*fr.sigma) {
                        cr.num[ns-1]+=w;
                        cr.num_w2[ns-1]+=w*w;
                    }
                } // endfor
            }
        }

        ib=ft_bin(b,bd);
        if (ib>=0) {
            CountResult& cr=r.ft_count[ib];
            cr.denom += w;
            cr.denom_w2 += w*w;

            const FitResult& fr=r.ft_fit[ib];
            int k=best_probe_candidate(b,0);
            if (fr.valid && k>=0) {
                const double dp=b.neutral_p[k]-b.probe_corr_p;
                const double d=std::fabs(dp-fr.mean);
                for (int ns=1;ns<=3;ns++) {
                    if (d < ns*fr.sigma) {
                        cr.num[ns-1]+=w;
                        cr.num_w2[ns-1]+=w*w;
                    }
                } // endfor
            }
        }
    } // endfor

    return true;
}

double efficiency(const CountResult& c, int ns) {
    if (c.denom<=0) return -1;
    return c.num[ns-1]/c.denom;
}

double efficiency_error_unweighted(const CountResult& c, int ns) {
    if (c.denom<=0) return 0;
    double e=efficiency(c,ns);
    if (e<0 || e>1) return 0;
    return std::sqrt(std::max(0.0,e*(1.0-e)/c.denom));
}

// Approximate weighted efficiency uncertainty based on effective sample size.
// This is sufficient for QA; production weighted uncertainties should use the
// exact saved sumw/sumw2 treatment from Valerii's workflow.
double efficiency_error_weighted(const CountResult& c, int ns) {
    if (c.denom<=0 || c.denom_w2<=0) return 0;
    double e=efficiency(c,ns);
    if (e<0 || e>1) return 0;
    double neff=c.denom*c.denom/c.denom_w2;
    if (neff<=0) return 0;
    return std::sqrt(std::max(0.0,e*(1.0-e)/neff));
}

void write_fd_csv(const SampleResult& r, const std::string& out) {
    std::ofstream f(out);
    f << "sample,pbin,thetabin,phibin,pmin,pmax,thetamin,thetamax,phimin,phimax,"
      << "fit_valid,fit_candidates,mean,mean_err,sigma,sigma_err,chi2,ndf,"
      << "denominator,num_1sigma,num_2sigma,num_3sigma,"
      << "eff_1sigma,eff_2sigma,eff_3sigma,"
      << "efferr_1sigma,efferr_2sigma,efferr_3sigma\n";

    int idx=0;
    for (int ip=0;ip<(int)FD_P_EDGES.size()-1;ip++) {
        for (int it=0;it<(int)FD_THETA_EDGES.size()-1;it++) {
            for (int iph=0;iph<(int)PHI_EDGES.size()-1;iph++,idx++) {
                const FitResult& fr=r.fd_fit[idx];
                const CountResult& cr=r.fd_count[idx];
                f << r.name << "," << ip << "," << it << "," << iph << ","
                  << FD_P_EDGES[ip] << "," << FD_P_EDGES[ip+1] << ","
                  << FD_THETA_EDGES[it] << "," << FD_THETA_EDGES[it+1] << ","
                  << PHI_EDGES[iph] << "," << PHI_EDGES[iph+1] << ","
                  << (fr.valid?1:0) << "," << fr.candidates << ","
                  << fr.mean << "," << fr.mean_err << ","
                  << fr.sigma << "," << fr.sigma_err << ","
                  << fr.chi2 << "," << fr.ndf << ","
                  << cr.denom << "," << cr.num[0] << "," << cr.num[1] << "," << cr.num[2];

                for (int ns=1;ns<=3;ns++) f << "," << efficiency(cr,ns);
                for (int ns=1;ns<=3;ns++)
                    f << "," << (r.weighted ? efficiency_error_weighted(cr,ns)
                                             : efficiency_error_unweighted(cr,ns));
                f << "\n";
            } // endfor
        } // endfor
    } // endfor
}

void write_ft_csv(const SampleResult& r, const std::string& out) {
    std::ofstream f(out);
    f << "sample,pbin,thetabin,pmin,pmax,thetamin,thetamax,"
      << "fit_valid,fit_candidates,mean,mean_err,sigma,sigma_err,chi2,ndf,"
      << "denominator,num_1sigma,num_2sigma,num_3sigma,"
      << "eff_1sigma,eff_2sigma,eff_3sigma,"
      << "efferr_1sigma,efferr_2sigma,efferr_3sigma\n";

    int idx=0;
    for (int ip=0;ip<(int)FT_P_EDGES.size()-1;ip++) {
        for (int it=0;it<(int)FT_THETA_EDGES.size()-1;it++,idx++) {
            const FitResult& fr=r.ft_fit[idx];
            const CountResult& cr=r.ft_count[idx];
            f << r.name << "," << ip << "," << it << ","
              << FT_P_EDGES[ip] << "," << FT_P_EDGES[ip+1] << ","
              << FT_THETA_EDGES[it] << "," << FT_THETA_EDGES[it+1] << ","
              << (fr.valid?1:0) << "," << fr.candidates << ","
              << fr.mean << "," << fr.mean_err << ","
              << fr.sigma << "," << fr.sigma_err << ","
              << fr.chi2 << "," << fr.ndf << ","
              << cr.denom << "," << cr.num[0] << "," << cr.num[1] << "," << cr.num[2];

            for (int ns=1;ns<=3;ns++) f << "," << efficiency(cr,ns);
            for (int ns=1;ns<=3;ns++)
                f << "," << (r.weighted ? efficiency_error_weighted(cr,ns)
                                         : efficiency_error_unweighted(cr,ns));
            f << "\n";
        } // endfor
    } // endfor
}

void save_global_residual(const SampleResult& r,
                          const std::string& outdir,
                          bool ft) {
    TH1D* h = ft ? r.ft_global.get() : r.fd_global.get();
    if (!h) return;

    TCanvas c(Form("c_%s_%s_global",r.name.c_str(),ft?"ft":"fd"),
              "",900,700);
    h->SetLineWidth(2);
    h->Draw("E");

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.035);
    tx.DrawLatex(0.15,0.86,Form("%s, %s",r.name.c_str(),ft?"FT":"FD"));
    tx.DrawLatex(0.15,0.81,Form("entries = %.0f",h->GetEntries()));

    const std::string path=outdir+"/"+r.name+"_"+(ft?"ft":"fd")+"_global_delta_p.png";
    c.SaveAs(path.c_str());
}

void save_fd_sigma_map(const SampleResult& r, const std::string& outdir) {
    // Phi-integrated event-weighted mean fitted sigma for each p/theta bin.
    TH2D h(Form("h_%s_fd_sigma_map",r.name.c_str()),
           ";p_{miss} [GeV];#theta_{miss} [deg]",
           FD_P_EDGES.size()-1,FD_P_EDGES.data(),
           FD_THETA_EDGES.size()-1,FD_THETA_EDGES.data());

    for (int ip=0;ip<(int)FD_P_EDGES.size()-1;ip++) {
        for (int it=0;it<(int)FD_THETA_EDGES.size()-1;it++) {
            double sw=0, ss=0;
            for (int iph=0;iph<(int)PHI_EDGES.size()-1;iph++) {
                int idx=(ip*3+it)*6+iph;
                const FitResult& fr=r.fd_fit[idx];
                if (!fr.valid) continue;
                double w=std::max(1LL,fr.candidates);
                sw+=w;
                ss+=w*fr.sigma;
            } // endfor
            if (sw>0) h.SetBinContent(ip+1,it+1,ss/sw);
        } // endfor
    } // endfor

    TCanvas c(Form("c_%s_fd_sigma",r.name.c_str()),"",1000,650);
    gStyle->SetPaintTextFormat(".3f");
    h.Draw("COLZ TEXT");
    c.SaveAs((outdir+"/"+r.name+"_fd_sigma_map.png").c_str());
}

void save_ft_efficiency(const SampleResult& r, const std::string& outdir) {
    // Theta-integrated efficiency versus missing-probe p for 1/2/3 sigma.
    const int np=FT_P_EDGES.size()-1;
    std::vector<double> x(np), ex(np), y1(np),y2(np),y3(np),e1(np),e2(np),e3(np);

    for (int ip=0;ip<np;ip++) {
        CountResult sum;
        for (int it=0;it<(int)FT_THETA_EDGES.size()-1;it++) {
            int idx=ip*2+it;
            const CountResult& c=r.ft_count[idx];
            sum.denom+=c.denom; sum.denom_w2+=c.denom_w2;
            for (int k=0;k<3;k++) {
                sum.num[k]+=c.num[k]; sum.num_w2[k]+=c.num_w2[k];
            } // endfor
        } // endfor

        x[ip]=0.5*(FT_P_EDGES[ip]+FT_P_EDGES[ip+1]);
        ex[ip]=0.5*(FT_P_EDGES[ip+1]-FT_P_EDGES[ip]);

        double* yy[3]={&y1[ip],&y2[ip],&y3[ip]};
        double* ee[3]={&e1[ip],&e2[ip],&e3[ip]};
        for (int ns=1;ns<=3;ns++) {
            *yy[ns-1]=efficiency(sum,ns);
            *ee[ns-1]=r.weighted ? efficiency_error_weighted(sum,ns)
                                  : efficiency_error_unweighted(sum,ns);
        } // endfor
    } // endfor

    TGraphErrors g1(np,x.data(),y1.data(),ex.data(),e1.data());
    TGraphErrors g2(np,x.data(),y2.data(),ex.data(),e2.data());
    TGraphErrors g3(np,x.data(),y3.data(),ex.data(),e3.data());

    TCanvas c(Form("c_%s_ft_eff",r.name.c_str()),"",900,700);
    TH1D frame("frame_ft",";p_{miss} [GeV];raw matched fraction",
               100,FT_P_EDGES.front(),FT_P_EDGES.back());
    frame.SetMinimum(0);
    frame.SetMaximum(1.05);
    frame.Draw();

    g1.SetMarkerStyle(20); g1.Draw("P SAME");
    g2.SetMarkerStyle(21); g2.Draw("P SAME");
    g3.SetMarkerStyle(22); g3.Draw("P SAME");

    TLegend leg(0.68,0.18,0.88,0.35);
    leg.AddEntry(&g1,"1#sigma","p");
    leg.AddEntry(&g2,"2#sigma","p");
    leg.AddEntry(&g3,"3#sigma","p");
    leg.Draw();

    c.SaveAs((outdir+"/"+r.name+"_ft_raw_efficiency_vs_p.png").c_str());
}

void save_summary(const SampleResult& r, const std::string& outdir) {
    std::ofstream f(outdir+"/"+r.name+"_summary.txt");
    auto both=[&](const std::string& s) {
        std::cout << s;
        f << s;
    };

    std::ostringstream s;
    s << "\n============================================================\n"
      << "Sample: " << r.name << "\n"
      << "Input: " << r.input << "\n"
      << "Tree entries: " << r.entries << "\n"
      << "Common pi0/tag-quality denominator candidates: " << r.denom_common << "\n"
      << "FD denominator rows in Valerii p/theta/phi range: " << r.fd_denom_rows << "\n"
      << "FD rows with reconstructed PID22 candidate: " << r.fd_probe_candidates << "\n"
      << "FT denominator rows in coarse p/theta range: " << r.ft_denom_rows << "\n"
      << "FT rows with reconstructed PID22 candidate: " << r.ft_probe_candidates << "\n"
      << "FD valid residual fits: "
      << std::count_if(r.fd_fit.begin(),r.fd_fit.end(),[](const FitResult& q){return q.valid;})
      << " / " << r.fd_fit.size() << "\n"
      << "FT valid residual fits: "
      << std::count_if(r.ft_fit.begin(),r.ft_fit.end(),[](const FitResult& q){return q.valid;})
      << " / " << r.ft_fit.size() << "\n"
      << "============================================================\n";
    both(s.str());
}

void process_and_write(const std::string& name,
                       const std::string& input,
                       bool weighted,
                       const std::string& outdir,
                       std::vector<std::unique_ptr<SampleResult>>& results) {
    if (input.empty()) return;

    std::unique_ptr<SampleResult> r(new SampleResult);
    r->name=name;
    r->input=input;
    r->weighted=weighted;

    if (!load_sample(*r)) {
        std::cerr << "WARNING: sample " << name << " could not be loaded; skipping.\n";
        return;
    }

    write_fd_csv(*r,outdir+"/"+name+"_fd_bins.csv");
    write_ft_csv(*r,outdir+"/"+name+"_ft_bins.csv");
    save_global_residual(*r,outdir,false);
    save_global_residual(*r,outdir,true);
    save_fd_sigma_map(*r,outdir);
    save_ft_efficiency(*r,outdir);
    save_summary(*r,outdir);

    results.push_back(std::move(r));
}

} // namespace per

void photon_efficiency_valerii_reproduction(
        const char* data_input,
        const char* output_dir="valerii_reproduction_output",
        const char* aaogen_input="",
        const char* clasdis_input="",
        const char* dvcsgen_input="") {

    using namespace per;

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    std::string out=output_dir ? output_dir : "valerii_reproduction_output";
    ensure_dir(out);

    std::cout
        << "\n"
        << "============================================================\n"
        << " Valerii-style photon-efficiency reproduction: first pass\n"
        << "============================================================\n"
        << "FD binning:\n"
        << "  p     = [0.35,0.50,1.10,1.70,2.30,2.90,3.90,6.00] GeV\n"
        << "  theta = [6,20,27,36] deg\n"
        << "  phi   = [-30,30,90,150,210,270,330] deg\n"
        << "Matching observable: Delta p_gamma2 = p_rec - p_miss\n"
        << "Matching windows: 1, 2, 3 fitted sigma\n"
        << "Parent pi0 window (FIRST-PASS knob): "
        << PI0_M_MIN << " < Mx(ep) < " << PI0_M_MAX << " GeV\n"
        << "\nIMPORTANT LIMITATION:\n"
        << "  The current skim does not save beta/fiducial flags for the\n"
        << "  reconstructed probe candidate, and does not yet project the\n"
        << "  missing FT probe to the FT face.  Therefore this run is a\n"
        << "  reproduction/QA pass, not the final correction map.\n"
        << "============================================================\n";

    std::vector<std::unique_ptr<SampleResult>> results;

    process_and_write("data",data_input?data_input:"",false,out,results);
    process_and_write("aaogen",aaogen_input?aaogen_input:"",true,out,results);
    process_and_write("clasdis",clasdis_input?clasdis_input:"",true,out,results);
    process_and_write("dvcsgen",dvcsgen_input?dvcsgen_input:"",true,out,results);

    // Write all residual histograms into one ROOT file for subsequent inspection.
    TFile fout((out+"/valerii_reproduction_histograms.root").c_str(),"RECREATE");
    for (const auto& rp : results) {
        fout.mkdir(rp->name.c_str());
        fout.cd(rp->name.c_str());

        if (rp->fd_global) rp->fd_global->Write();
        if (rp->ft_global) rp->ft_global->Write();
        for (const auto& h : rp->fd_h) if (h) h->Write();
        for (const auto& h : rp->ft_h) if (h) h->Write();

        fout.cd();
    } // endfor
    fout.Close();

    std::cout
        << "\nFinished. Products written to: " << out << "\n"
        << "First things to inspect:\n"
        << "  data_summary.txt\n"
        << "  data_fd_global_delta_p.png\n"
        << "  data_fd_sigma_map.png\n"
        << "  data_fd_bins.csv\n"
        << "  data_ft_global_delta_p.png\n"
        << "  data_ft_raw_efficiency_vs_p.png\n"
        << "  data_ft_bins.csv\n";
}
