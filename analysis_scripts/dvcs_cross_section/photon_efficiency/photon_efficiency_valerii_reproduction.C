// photon_efficiency_valerii_reproduction.C
//
// First-stage Valerii-style photon-efficiency extraction.
// Deliberately FULLY INTEGRATED in FD and FT.
//
// We are starting simple:
//   1. one integrated FD residual fit;
//   2. one integrated FT residual fit;
//   3. 1-, 2-, and 3-sigma matched fractions from those fits;
//   4. identical machinery for data / AAOgen / CLASDIS / DVCSgen.
//
// The residual is Valerii's momentum matching variable:
//      Delta p_gamma2 = p_rec - p_miss
//
// IMPORTANT:
//   * The current skim does not yet save candidate-level probe beta/fiducial
//     flags, so the reconstructed probe requirement is PID 22, correct detector,
//     p >= 0.4 GeV.
//   * The FT denominator is still theta-defined (2 <= theta <= 5.5 deg), not a
//     projected FT-face fiducial.  FT is therefore preliminary.
//   * The parent pi0 window remains a first-pass enrichment cut.
//
// Usage, data only:
// root -l -b -q 'photon_efficiency_valerii_reproduction.C(
//   "/work/clas12/thayward/photon_efficiency/ROOT_trees/data/fa18_inb",
//   "valerii_reproduction_fa18_inb")'
//
// Usage with MC:
// root -l -b -q 'photon_efficiency_valerii_reproduction.C(
//   "/work/.../data/fa18_inb",
//   "valerii_reproduction_fa18_inb",
//   "/work/.../aaogen/fa18_inb",
//   "/work/.../clasdis/fa18_inb",
//   "/work/.../dvcsgen/fa18_inb")'

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
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
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace pe {

// -----------------------------------------------------------------------------
// First-pass selection knobs.
// -----------------------------------------------------------------------------
static const double PI0_M_MIN = 0.08;
static const double PI0_M_MAX = 0.20;

static const double PROBE_P_MIN = 0.40;

// FD region follows the Valerii analysis angular range.
static const double FD_THETA_MIN = 6.0;
static const double FD_THETA_MAX = 36.0;
static const double FD_P_MIN = 0.35;
static const double FD_P_MAX = 6.0;

// FT first-pass angular region.
static const double FT_THETA_MIN = 2.0;
static const double FT_THETA_MAX = 5.5;
static const double FT_P_MIN = 0.40;
static const double FT_P_MAX = 8.5;

// Residual histogram / fit range.
static const int    DP_NBIN = 100;
static const double DP_HMIN = -1.0;
static const double DP_HMAX =  1.0;
static const double FIT_MIN = -0.60;
static const double FIT_MAX =  0.60;

static const int MIN_FIT_ENTRIES_FD = 100;
static const int MIN_FIT_ENTRIES_FT = 20;

// Strict fit-quality requirements.  ROOT fit status alone is not enough.
static const double SIGMA_MIN = 0.010;
static const double SIGMA_MAX = 0.400;
static const double MAX_REL_SIGMA_ERR = 0.50;
static const double MAX_MEAN_ERR = 0.100;

// -----------------------------------------------------------------------------
// Helpers.
// -----------------------------------------------------------------------------
bool finite_good(double x) {
    return std::isfinite(x) && x > -900.0;
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

struct Branches {
    Int_t runnum=0, evnum=0, is_mc=0;
    Int_t p_pass_standard=0;
    Int_t tag_detector=-1, tag_pass_beta=0, tag_pass_fiducial=0;
    Double_t mc_weight=1.0;

    Double_t Mx_ep=0;
    Double_t probe_corr_p=0, probe_corr_theta=0, probe_corr_phi=0;

    Int_t neutral_idx[5];
    Int_t neutral_pid[5];
    Int_t neutral_charge[5];
    Int_t neutral_detector[5];
    Double_t neutral_p[5];
    Double_t neutral_delta_alpha[5];

    void reset_arrays() {
        for (int i=0;i<5;i++) {
            neutral_idx[i]=-999;
            neutral_pid[i]=-999;
            neutral_charge[i]=-999;
            neutral_detector[i]=-999;
            neutral_p[i]=-999;
            neutral_delta_alpha[i]=-999;
        } // endfor
    }
};

bool attach(TChain& c, Branches& b) {
    const char* required[] = {
        "runnum","evnum","is_mc","mc_weight","p_pass_standard",
        "tag_detector","tag_pass_beta","tag_pass_fiducial",
        "Mx_ep","probe_corr_p","probe_corr_theta","probe_corr_phi",
        "neutral_idx","neutral_pid","neutral_charge","neutral_detector",
        "neutral_p","neutral_delta_alpha"
    };

    bool ok=true;
    for (const char* n : required) {
        if (!c.GetBranch(n)) {
            std::cerr << "ERROR: missing branch " << n << "\n";
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
    c.SetBranchAddress("neutral_delta_alpha",b.neutral_delta_alpha);

    return true;
}

bool common_selection(const Branches& b) {
    if (!b.p_pass_standard) return false;
    if (!b.tag_pass_beta || !b.tag_pass_fiducial) return false;
    if (!(b.tag_detector==0 || b.tag_detector==1)) return false;

    if (!finite_good(b.Mx_ep)) return false;
    if (b.Mx_ep < PI0_M_MIN || b.Mx_ep > PI0_M_MAX) return false;

    if (!finite_good(b.probe_corr_p) || !finite_good(b.probe_corr_theta))
        return false;

    return true;
}

bool in_fd(const Branches& b) {
    return b.probe_corr_p >= FD_P_MIN &&
           b.probe_corr_p <= FD_P_MAX &&
           b.probe_corr_theta >= FD_THETA_MIN &&
           b.probe_corr_theta <= FD_THETA_MAX;
}

bool in_ft(const Branches& b) {
    return b.probe_corr_p >= FT_P_MIN &&
           b.probe_corr_p <= FT_P_MAX &&
           b.probe_corr_theta >= FT_THETA_MIN &&
           b.probe_corr_theta <= FT_THETA_MAX;
}

// Closest saved PID-22 neutral in the requested detector.
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
    } // endfor

    return best;
}

double event_weight(const Branches& b, bool weighted) {
    if (!weighted) return 1.0;
    if (!std::isfinite(b.mc_weight)) return 0.0;
    return b.mc_weight;
}

struct FitResult {
    bool valid=false;
    int root_status=-999;
    long long candidates=0;

    double mean=0;
    double mean_err=0;
    double sigma=0;
    double sigma_err=0;
    double chi2=0;
    int ndf=0;

    std::string reason="not fit";
};

FitResult fit_integrated_residual(TH1D* h, int min_entries) {
    FitResult r;
    if (!h) {
        r.reason="null histogram";
        return r;
    }

    r.candidates=static_cast<long long>(h->GetEntries());

    if (h->GetEntries() < min_entries) {
        r.reason="too few candidates";
        return r;
    }

    // Use the visible maximum as the seed.
    const int maxbin=h->GetMaximumBin();
    const double peak=h->GetBinCenter(maxbin);

    // Seed Gaussian in a narrow window around the peak.
    const double seed_lo=std::max(FIT_MIN,peak-0.20);
    const double seed_hi=std::min(FIT_MAX,peak+0.20);

    TF1 seed("seed_integrated","gaus",seed_lo,seed_hi);
    seed.SetParameters(h->GetMaximum(),peak,0.10);
    seed.SetParLimits(2,SIGMA_MIN,SIGMA_MAX);

    const int seed_status=h->Fit(&seed,"QNR");

    double mu=peak;
    double sg=0.10;
    if (seed_status==0) {
        mu=seed.GetParameter(1);
        sg=std::fabs(seed.GetParameter(2));
    }

    if (!(sg>SIGMA_MIN && sg<SIGMA_MAX))
        sg=0.10;

    // Final fit: Gaussian signal + linear background.
    double lo=std::max(FIT_MIN,mu-4.0*sg);
    double hi=std::min(FIT_MAX,mu+4.0*sg);
    if (hi-lo < 0.25) {
        lo=FIT_MIN;
        hi=FIT_MAX;
    }

    TF1 f("fit_integrated","gaus(0)+pol1(3)",lo,hi);
    f.SetParameters(h->GetMaximum(),mu,sg,
                    std::max(0.0,h->GetBinContent(1)),0.0);
    f.SetParLimits(1,-0.50,0.50);
    f.SetParLimits(2,SIGMA_MIN,SIGMA_MAX);

    r.root_status=h->Fit(&f,"QNR");

    r.mean=f.GetParameter(1);
    r.mean_err=f.GetParError(1);
    r.sigma=std::fabs(f.GetParameter(2));
    r.sigma_err=f.GetParError(2);
    r.chi2=f.GetChisquare();
    r.ndf=f.GetNDF();

    // Strict validity checks.
    if (r.root_status!=0) {
        r.reason="ROOT fit status != 0";
        return r;
    }

    if (!(r.sigma > SIGMA_MIN*1.01 && r.sigma < SIGMA_MAX*0.99)) {
        r.reason="sigma at/near fit boundary";
        return r;
    }

    if (!(r.sigma_err>0) || r.sigma_err/r.sigma > MAX_REL_SIGMA_ERR) {
        r.reason="sigma uncertainty too large";
        return r;
    }

    if (!(r.mean_err>=0) || r.mean_err > MAX_MEAN_ERR) {
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

struct DetectorResult {
    std::string name;
    int detector=-1;

    std::unique_ptr<TH1D> residual;
    FitResult fit;

    double denominator=0;
    double denominator_w2=0;
    double numerator[3]={0,0,0};
    double numerator_w2[3]={0,0,0};

    long long denominator_rows=0;
    long long reconstructed_candidates=0;
};

struct SampleResult {
    std::string name;
    std::string input;
    bool weighted=false;

    long long entries=0;
    long long common_rows=0;

    DetectorResult fd;
    DetectorResult ft;
};

void initialize_detector(DetectorResult& d,
                         const std::string& sample,
                         const std::string& det,
                         int detector) {
    d.name=det;
    d.detector=detector;
    d.residual.reset(new TH1D(
        Form("h_%s_%s_delta_p",sample.c_str(),det.c_str()),
        Form("%s %s;#Delta p_{#gamma2}=p_{rec}-p_{miss} [GeV];Candidates",
             sample.c_str(),det.c_str()),
        DP_NBIN,DP_HMIN,DP_HMAX));
    d.residual->Sumw2();
    d.residual->SetDirectory(nullptr);
}

bool load_sample(SampleResult& r) {
    const std::string pattern=make_pattern(r.input.c_str());
    if (pattern.empty()) return false;

    TChain c("PhotonEfficiency");
    const int nf=c.Add(pattern.c_str());
    r.entries=c.GetEntries();

    std::cout << "\n[" << r.name << "] files=" << nf
              << " entries=" << r.entries
              << " pattern=" << pattern << "\n";

    if (nf<=0 || r.entries<=0) return false;

    Branches b;
    b.reset_arrays();
    if (!attach(c,b)) return false;

    initialize_detector(r.fd,r.name,"FD",1);
    initialize_detector(r.ft,r.name,"FT",0);

    // ---------------------------------------------------------------------
    // Pass 1: build the fully integrated residual distributions.
    // ---------------------------------------------------------------------
    for (Long64_t i=0;i<c.GetEntries();i++) {
        c.GetEntry(i);

        if (!common_selection(b)) continue;
        r.common_rows++;

        const double w=event_weight(b,r.weighted);
        if (w==0.0) continue;

        if (in_fd(b)) {
            r.fd.denominator_rows++;

            const int k=best_probe_candidate(b,1);
            if (k>=0) {
                r.fd.reconstructed_candidates++;
                const double dp=b.neutral_p[k]-b.probe_corr_p;
                if (std::isfinite(dp))
                    r.fd.residual->Fill(dp,w);
            }
        }

        if (in_ft(b)) {
            r.ft.denominator_rows++;

            const int k=best_probe_candidate(b,0);
            if (k>=0) {
                r.ft.reconstructed_candidates++;
                const double dp=b.neutral_p[k]-b.probe_corr_p;
                if (std::isfinite(dp))
                    r.ft.residual->Fill(dp,w);
            }
        }
    } // endfor

    r.fd.fit=fit_integrated_residual(r.fd.residual.get(),MIN_FIT_ENTRIES_FD);
    r.ft.fit=fit_integrated_residual(r.ft.residual.get(),MIN_FIT_ENTRIES_FT);

    // ---------------------------------------------------------------------
    // Pass 2: fully integrated denominator and 1/2/3-sigma numerators.
    // ---------------------------------------------------------------------
    for (Long64_t i=0;i<c.GetEntries();i++) {
        c.GetEntry(i);

        if (!common_selection(b)) continue;

        const double w=event_weight(b,r.weighted);
        if (w==0.0) continue;

        if (in_fd(b)) {
            r.fd.denominator+=w;
            r.fd.denominator_w2+=w*w;

            const int k=best_probe_candidate(b,1);
            if (r.fd.fit.valid && k>=0) {
                const double dp=b.neutral_p[k]-b.probe_corr_p;
                const double d=std::fabs(dp-r.fd.fit.mean);

                for (int ns=1;ns<=3;ns++) {
                    if (d < ns*r.fd.fit.sigma) {
                        r.fd.numerator[ns-1]+=w;
                        r.fd.numerator_w2[ns-1]+=w*w;
                    }
                } // endfor
            }
        }

        if (in_ft(b)) {
            r.ft.denominator+=w;
            r.ft.denominator_w2+=w*w;

            const int k=best_probe_candidate(b,0);
            if (r.ft.fit.valid && k>=0) {
                const double dp=b.neutral_p[k]-b.probe_corr_p;
                const double d=std::fabs(dp-r.ft.fit.mean);

                for (int ns=1;ns<=3;ns++) {
                    if (d < ns*r.ft.fit.sigma) {
                        r.ft.numerator[ns-1]+=w;
                        r.ft.numerator_w2[ns-1]+=w*w;
                    }
                } // endfor
            }
        }
    } // endfor

    return true;
}

double eff(const DetectorResult& d, int ns) {
    if (d.denominator<=0) return -1;
    return d.numerator[ns-1]/d.denominator;
}

double eff_error(const DetectorResult& d, int ns, bool weighted) {
    if (d.denominator<=0) return 0;

    const double e=eff(d,ns);
    if (!(e>=0 && e<=1)) return 0;

    double n=d.denominator;
    if (weighted) {
        if (d.denominator_w2<=0) return 0;
        n=d.denominator*d.denominator/d.denominator_w2;
    }

    if (n<=0) return 0;
    return std::sqrt(std::max(0.0,e*(1.0-e)/n));
}

void save_residual_plot(const SampleResult& r,
                        const DetectorResult& d,
                        const std::string& outdir) {
    if (!d.residual) return;

    TCanvas c(Form("c_%s_%s",r.name.c_str(),d.name.c_str()),"",950,720);

    d.residual->SetLineWidth(2);
    d.residual->Draw("E");

    double ymax=d.residual->GetMaximum();

    if (d.fit.valid) {
        for (int ns=1;ns<=3;ns++) {
            const double lo=d.fit.mean-ns*d.fit.sigma;
            const double hi=d.fit.mean+ns*d.fit.sigma;

            TLine l1(lo,0,lo,ymax*0.92);
            TLine l2(hi,0,hi,ymax*0.92);
            l1.SetLineStyle(ns);
            l2.SetLineStyle(ns);
            l1.Draw();
            l2.Draw();
        } // endfor
    }

    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.033);

    tx.DrawLatex(0.13,0.88,Form("%s %s",r.name.c_str(),d.name.c_str()));
    tx.DrawLatex(0.13,0.83,Form("reconstructed candidates = %lld",
                                d.reconstructed_candidates));

    if (d.fit.valid) {
        tx.DrawLatex(0.13,0.78,
            Form("#mu = %.4f #pm %.4f GeV",d.fit.mean,d.fit.mean_err));
        tx.DrawLatex(0.13,0.73,
            Form("#sigma = %.4f #pm %.4f GeV",d.fit.sigma,d.fit.sigma_err));
        tx.DrawLatex(0.13,0.68,
            Form("#chi^{2}/ndf = %.1f/%d",d.fit.chi2,d.fit.ndf));
    } else {
        tx.DrawLatex(0.13,0.78,Form("FIT INVALID: %s",d.fit.reason.c_str()));
    }

    c.SaveAs((outdir+"/"+r.name+"_"+d.name+"_integrated_delta_p.png").c_str());
}

void write_summary(const SampleResult& r, const std::string& outdir) {
    std::ofstream f(outdir+"/"+r.name+"_integrated_summary.txt");

    auto print_detector=[&](const DetectorResult& d) {
        f << "\n--- " << d.name << " ----------------------------------------\n";
        f << "denominator rows: " << d.denominator_rows << "\n";
        f << "reconstructed PID22 candidates: " << d.reconstructed_candidates << "\n";
        f << "fit valid: " << (d.fit.valid?1:0) << "\n";
        f << "fit reason: " << d.fit.reason << "\n";
        f << "fit ROOT status: " << d.fit.root_status << "\n";
        f << "mean: " << d.fit.mean << " +/- " << d.fit.mean_err << " GeV\n";
        f << "sigma: " << d.fit.sigma << " +/- " << d.fit.sigma_err << " GeV\n";
        f << "chi2/ndf: " << d.fit.chi2 << "/" << d.fit.ndf << "\n";
        f << "weighted denominator: " << d.denominator << "\n";

        for (int ns=1;ns<=3;ns++) {
            f << ns << " sigma numerator: " << d.numerator[ns-1] << "\n";
            f << ns << " sigma matched fraction: "
              << eff(d,ns) << " +/- " << eff_error(d,ns,r.weighted) << "\n";
        } // endfor
    };

    f << "============================================================\n";
    f << "Fully integrated Valerii-style first-stage extraction\n";
    f << "sample: " << r.name << "\n";
    f << "input: " << r.input << "\n";
    f << "tree entries: " << r.entries << "\n";
    f << "common selected rows: " << r.common_rows << "\n";
    f << "pi0 enrichment window: " << PI0_M_MIN << " < Mx(ep) < "
      << PI0_M_MAX << " GeV\n";

    print_detector(r.fd);
    print_detector(r.ft);

    f << "\nIMPORTANT: matched fractions are preliminary.  Probe-level beta/fiducial\n"
      << "flags and a projected FT-face denominator are not yet available in the skim.\n";
}

void write_csv(const std::vector<std::unique_ptr<SampleResult>>& results,
               const std::string& outdir) {
    std::ofstream f(outdir+"/integrated_results.csv");
    f << "sample,detector,fit_valid,fit_reason,fit_candidates,"
      << "mean,mean_err,sigma,sigma_err,chi2,ndf,"
      << "denominator,reconstructed_candidates,"
      << "num_1sigma,num_2sigma,num_3sigma,"
      << "eff_1sigma,eff_2sigma,eff_3sigma,"
      << "err_1sigma,err_2sigma,err_3sigma\n";

    for (const auto& rp : results) {
        const DetectorResult* ds[2]={&rp->fd,&rp->ft};

        for (const DetectorResult* d : ds) {
            f << rp->name << "," << d->name << ","
              << (d->fit.valid?1:0) << ","
              << "\"" << d->fit.reason << "\"" << ","
              << d->fit.candidates << ","
              << d->fit.mean << "," << d->fit.mean_err << ","
              << d->fit.sigma << "," << d->fit.sigma_err << ","
              << d->fit.chi2 << "," << d->fit.ndf << ","
              << d->denominator << ","
              << d->reconstructed_candidates << ","
              << d->numerator[0] << "," << d->numerator[1] << "," << d->numerator[2] << ","
              << eff(*d,1) << "," << eff(*d,2) << "," << eff(*d,3) << ","
              << eff_error(*d,1,rp->weighted) << ","
              << eff_error(*d,2,rp->weighted) << ","
              << eff_error(*d,3,rp->weighted) << "\n";
        } // endfor
    } // endfor
}

void process_sample(const std::string& name,
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
        std::cerr << "WARNING: could not load " << name << "\n";
        return;
    }

    save_residual_plot(*r,r->fd,outdir);
    save_residual_plot(*r,r->ft,outdir);
    write_summary(*r,outdir);

    results.push_back(std::move(r));
}

} // namespace pe

void photon_efficiency_valerii_reproduction(
        const char* data_input,
        const char* aaogen_input="",
        const char* clasdis_input="",
        const char* dvcsgen_input="") {

    using namespace pe;

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // Keep one clean analysis-output directory.  Every invocation replaces
    // the previous contents so iterative development does not accumulate
    // versioned output directories.
    std::string out="output";
    if (gSystem->AccessPathName(out.c_str()) == kFALSE) {
        gSystem->Exec(("rm -rf " + out).c_str());
    }
    ensure_dir(out);

    std::cout
        << "\n============================================================\n"
        << " Fully integrated Valerii-style photon-efficiency extraction\n"
        << "============================================================\n"
        << "FD: " << FD_P_MIN << " <= p <= " << FD_P_MAX
        << " GeV, " << FD_THETA_MIN << " <= theta <= "
        << FD_THETA_MAX << " deg\n"
        << "FT: " << FT_P_MIN << " <= p <= " << FT_P_MAX
        << " GeV, " << FT_THETA_MIN << " <= theta <= "
        << FT_THETA_MAX << " deg\n"
        << "Matching: Delta p_gamma2 = p_rec - p_miss\n"
        << "Windows: 1, 2, 3 fitted sigma\n"
        << "============================================================\n";

    std::vector<std::unique_ptr<SampleResult>> results;

    process_sample("data",data_input?data_input:"",false,out,results);
    process_sample("aaogen",aaogen_input?aaogen_input:"",true,out,results);
    process_sample("clasdis",clasdis_input?clasdis_input:"",true,out,results);
    process_sample("dvcsgen",dvcsgen_input?dvcsgen_input:"",true,out,results);

    write_csv(results,out);

    TFile fout((out+"/integrated_residuals.root").c_str(),"RECREATE");
    for (const auto& rp : results) {
        if (rp->fd.residual) rp->fd.residual->Write();
        if (rp->ft.residual) rp->ft.residual->Write();
    } // endfor
    fout.Close();

    std::cout
        << "\nFinished. Main outputs:\n"
        << "  " << out << "/integrated_results.csv\n"
        << "  " << out << "/data_integrated_summary.txt\n"
        << "  " << out << "/data_FD_integrated_delta_p.png\n"
        << "  " << out << "/data_FT_integrated_delta_p.png\n"
        << "  " << out << "/integrated_residuals.root\n";
}
