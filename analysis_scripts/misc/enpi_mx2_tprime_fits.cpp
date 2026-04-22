#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TPad.h>
#include <TLine.h>
#include <TROOT.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

struct XBSlice {
    std::string name;
    double xmin;
    double xmax;
};

struct TPrimeBin {
    double tmin;
    double tmax;
};

struct FitResult {
    std::string slice_name;
    double x_min = 0.0;
    double x_max = 0.0;
    double t_min = 0.0;
    double t_max = 0.0;
    double mu = std::numeric_limits<double>::quiet_NaN();
    double mu_err = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double sigma_err = std::numeric_limits<double>::quiet_NaN();
    double amp = std::numeric_limits<double>::quiet_NaN();
    double p0 = std::numeric_limits<double>::quiet_NaN();
    double p1 = std::numeric_limits<double>::quiet_NaN();
    double p2 = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    int fit_status = -1;
    int cov_status = -1;
    bool fit_success = false;
    long long n_entries = 0;
};

struct FitConfig {
    double mu_min;
    double mu_max;
    double sigma_min;
    double sigma_max;
    double peak_search_min;
    double peak_search_max;
    double prefit_halfwidth;
    double final_fit_min;
    double final_fit_max;
    double sigma_seed;
};

struct SummaryPoint {
    double x_min;
    double x_max;
    double t_min;
    double t_max;
    double mu;
    double mu_err;
    double sigma;
    double sigma_err;
};

static const std::string kInputFile = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_combined_inb_NH3_epi+_2.root";
static const std::string kTreeName = "PhysicsEvents";
static const std::string kOutputDir = "output/enpi+_Mx2_fits";

static constexpr double kMx2Min = 0.70;
static constexpr double kMx2Max = 1.30;
static constexpr int kNMx2Bins = 15;

static constexpr double kMuMin = 0.85;
static constexpr double kMuMax = 0.95;
static constexpr double kSigmaMin = 0.055;
static constexpr double kSigmaMax = 0.12;

static constexpr double kPeakSearchMin = 0.84;
static constexpr double kPeakSearchMax = 0.99;
static constexpr double kPrefitHalfWidth = 0.055;
static constexpr double kFinalFitMin = 0.70;
static constexpr double kFinalFitMax = 1.30;
static constexpr double kSigmaSeed = 0.06;

/*
 * Special carve-out 1:
 * lowest xB, lowest -tprime  -> (ix,it) = (0,0)
 * Problem: Gaussian center was too far left.
 */
static constexpr double kSpecial1MuMin = 0.90;
static constexpr double kSpecial1MuMax = 1.02;
static constexpr double kSpecial1SigmaMin = 0.060;
static constexpr double kSpecial1SigmaMax = 0.16;
static constexpr double kSpecial1PeakSearchMin = 0.90;
static constexpr double kSpecial1PeakSearchMax = 1.05;
static constexpr double kSpecial1PrefitHalfWidth = 0.080;
static constexpr double kSpecial1FinalFitMin = 0.78;
static constexpr double kSpecial1FinalFitMax = 1.20;
static constexpr double kSpecial1SigmaSeed = 0.085;

/*
 * Special carve-out 2:
 * second row, third column -> (ix,it) = (1,2)
 * Problem: Gaussian center is too far right.
 */
static constexpr double kSpecial2MuMin = 0.88;
static constexpr double kSpecial2MuMax = 0.94;
static constexpr double kSpecial2SigmaMin = 0.060;
static constexpr double kSpecial2SigmaMax = 0.12;
static constexpr double kSpecial2PeakSearchMin = 0.86;
static constexpr double kSpecial2PeakSearchMax = 0.97;
static constexpr double kSpecial2PrefitHalfWidth = 0.045;
static constexpr double kSpecial2FinalFitMin = 0.70;
static constexpr double kSpecial2FinalFitMax = 1.18;
static constexpr double kSpecial2SigmaSeed = 0.075;

static const std::vector<XBSlice> kXBSlices = {
    {"Low",     0.10, 0.25},
    {"MidLow",  0.25, 0.35},
    {"MidHigh", 0.35, 0.45},
    {"High",    0.45, 0.60}
};

static const std::vector<TPrimeBin> kTPrimeBins = {
    {0.05, 0.25},
    {0.25, 0.45},
    {0.45, 0.65},
    {0.65, 0.85},
    {0.85, 1.05},
    {1.05, 1.25}
};

bool is_special_bin_1(int ix, int it) {
    return (ix == 0 && it == 0);
}

bool is_special_bin_2(int ix, int it) {
    return (ix == 1 && it == 2);
}

FitConfig get_fit_config(int ix, int it) {
    if (is_special_bin_1(ix, it)) {
        return {
            kSpecial1MuMin,
            kSpecial1MuMax,
            kSpecial1SigmaMin,
            kSpecial1SigmaMax,
            kSpecial1PeakSearchMin,
            kSpecial1PeakSearchMax,
            kSpecial1PrefitHalfWidth,
            kSpecial1FinalFitMin,
            kSpecial1FinalFitMax,
            kSpecial1SigmaSeed
        };
    }

    if (is_special_bin_2(ix, it)) {
        return {
            kSpecial2MuMin,
            kSpecial2MuMax,
            kSpecial2SigmaMin,
            kSpecial2SigmaMax,
            kSpecial2PeakSearchMin,
            kSpecial2PeakSearchMax,
            kSpecial2PrefitHalfWidth,
            kSpecial2FinalFitMin,
            kSpecial2FinalFitMax,
            kSpecial2SigmaSeed
        };
    }

    return {
        kMuMin,
        kMuMax,
        kSigmaMin,
        kSigmaMax,
        kPeakSearchMin,
        kPeakSearchMax,
        kPrefitHalfWidth,
        kFinalFitMin,
        kFinalFitMax,
        kSigmaSeed
    };
}

int find_xbin(double x) {
    for (int i = 0; i < static_cast<int>(kXBSlices.size()); i++) {
        if (x > kXBSlices[i].xmin && x < kXBSlices[i].xmax) {
            return i;
        }
    }
    return -1;
}

int find_tbin(double minus_tprime) {
    for (int i = 0; i < static_cast<int>(kTPrimeBins.size()); i++) {
        if (minus_tprime >= kTPrimeBins[i].tmin && minus_tprime < kTPrimeBins[i].tmax) {
            return i;
        }
    }
    return -1;
}

std::string make_hist_name(int ix, int it) {
    std::ostringstream oss;
    oss << "h_mx2_x" << ix << "_t" << it;
    return oss.str();
}

std::string make_fit_name(int ix, int it) {
    std::ostringstream oss;
    oss << "f_x" << ix << "_t" << it;
    return oss.str();
}

std::string make_prefit_name(int ix, int it) {
    std::ostringstream oss;
    oss << "gprefit_x" << ix << "_t" << it;
    return oss.str();
}

std::string format_double(double x, int prec = 4) {
    if (!std::isfinite(x)) return "--";
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(prec) << x;
    return oss.str();
}

void setup_global_style() {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetLabelSize(0.040, "XY");
    gStyle->SetTitleSize(0.045, "XY");
    gStyle->SetTitleOffset(1.05, "X");
    gStyle->SetTitleOffset(1.15, "Y");
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
}

void ensure_output_dir(const std::string& dir) {
    gSystem->mkdir(dir.c_str(), true);
}

int find_peak_bin_in_window(TH1D* hist, double xmin, double xmax) {
    int bin_min = hist->FindBin(xmin);
    int bin_max = hist->FindBin(xmax);
    bin_min = std::max(bin_min, 1);
    bin_max = std::min(bin_max, hist->GetNbinsX());

    int best_bin = bin_min;
    double best_val = -1.0;

    for (int ibin = bin_min; ibin <= bin_max; ibin++) {
        double val = hist->GetBinContent(ibin);
        if (val > best_val) {
            best_val = val;
            best_bin = ibin;
        }
    }

    return best_bin;
}

double average_sideband_level(TH1D* hist) {
    double sum = 0.0;
    int n = 0;

    for (int ibin = 1; ibin <= hist->GetNbinsX(); ibin++) {
        double x = hist->GetBinCenter(ibin);
        if ((x >= 0.70 && x <= 0.80) || (x >= 1.02 && x <= 1.10)) {
            sum += hist->GetBinContent(ibin);
            n++;
        }
    }

    if (n == 0) return 0.0;
    return sum / static_cast<double>(n);
}

FitResult fit_histogram(TH1D* hist, int ix, int it) {
    FitResult result;
    result.slice_name = kXBSlices[ix].name;
    result.x_min = kXBSlices[ix].xmin;
    result.x_max = kXBSlices[ix].xmax;
    result.t_min = kTPrimeBins[it].tmin;
    result.t_max = kTPrimeBins[it].tmax;
    result.n_entries = static_cast<long long>(hist->GetEntries());

    if (result.n_entries < 50) {
        return result;
    }

    FitConfig cfg = get_fit_config(ix, it);

    int peak_bin = find_peak_bin_in_window(hist, cfg.peak_search_min, cfg.peak_search_max);
    double peak_x = hist->GetBinCenter(peak_bin);
    double peak_y = hist->GetBinContent(peak_bin);

    double bg_guess = average_sideband_level(hist);
    double amp_guess = std::max(peak_y - bg_guess, 1.0);

    double prefit_min = std::max(kMx2Min, peak_x - cfg.prefit_halfwidth);
    double prefit_max = std::min(kMx2Max, peak_x + cfg.prefit_halfwidth);

    std::string gprefit_name = make_prefit_name(ix, it);
    TF1 gprefit(gprefit_name.c_str(), "gaus", prefit_min, prefit_max);
    gprefit.SetParameter(0, amp_guess);
    gprefit.SetParameter(1, std::max(cfg.mu_min + 0.002, std::min(peak_x, cfg.mu_max - 0.002)));
    gprefit.SetParameter(2, cfg.sigma_seed);
    gprefit.SetParLimits(1, cfg.mu_min, cfg.mu_max);
    gprefit.SetParLimits(2, cfg.sigma_min, cfg.sigma_max);

    hist->Fit(&gprefit, "Q0R");

    double mu_guess = gprefit.GetParameter(1);
    double sigma_guess = gprefit.GetParameter(2);

    if (!std::isfinite(mu_guess) || mu_guess < cfg.mu_min || mu_guess > cfg.mu_max) {
        mu_guess = std::max(cfg.mu_min + 0.002, std::min(peak_x, cfg.mu_max - 0.002));
    }

    if (!std::isfinite(sigma_guess) || sigma_guess < cfg.sigma_min || sigma_guess > cfg.sigma_max) {
        sigma_guess = cfg.sigma_seed;
    }

    std::string fname = make_fit_name(ix, it);
    TF1 fit_fn(fname.c_str(), "gaus(0) + pol2(3)", cfg.final_fit_min, cfg.final_fit_max);
    fit_fn.SetLineColor(kRed + 1);
    fit_fn.SetLineStyle(2);
    fit_fn.SetLineWidth(2);

    fit_fn.SetParameter(0, std::max(gprefit.GetParameter(0), 1.0));
    fit_fn.SetParameter(1, mu_guess);
    fit_fn.SetParameter(2, sigma_guess);
    fit_fn.SetParameter(3, bg_guess);
    fit_fn.SetParameter(4, 0.0);
    fit_fn.SetParameter(5, 0.0);

    fit_fn.SetParLimits(0, 0.0, 1.0e12);
    fit_fn.SetParLimits(1, cfg.mu_min, cfg.mu_max);
    fit_fn.SetParLimits(2, cfg.sigma_min, cfg.sigma_max);

    TFitResultPtr fit_ptr = hist->Fit(&fit_fn, "QRSL");
    int fit_status = static_cast<int>(fit_ptr);
    int cov_status = -1;
    if (fit_ptr.Get()) {
        cov_status = fit_ptr->CovMatrixStatus();
    }

    result.fit_status = fit_status;
    result.cov_status = cov_status;

    result.amp = fit_fn.GetParameter(0);
    result.mu = fit_fn.GetParameter(1);
    result.mu_err = fit_fn.GetParError(1);
    result.sigma = fit_fn.GetParameter(2);
    result.sigma_err = fit_fn.GetParError(2);
    result.p0 = fit_fn.GetParameter(3);
    result.p1 = fit_fn.GetParameter(4);
    result.p2 = fit_fn.GetParameter(5);
    result.chi2 = fit_fn.GetChisquare();
    result.ndf = fit_fn.GetNDF();

    bool finite_ok =
        std::isfinite(result.mu) &&
        std::isfinite(result.mu_err) &&
        std::isfinite(result.sigma) &&
        std::isfinite(result.sigma_err);

    bool in_range_ok =
        (result.mu >= cfg.mu_min) &&
        (result.mu <= cfg.mu_max) &&
        (result.sigma >= cfg.sigma_min) &&
        (result.sigma <= cfg.sigma_max);

    result.fit_success = (finite_ok && in_range_ok);

    return result;
}

std::vector<SummaryPoint> get_hardcoded_summary_points() {
    std::vector<SummaryPoint> pts = {
        {0.1000, 0.2500, 0.0500, 0.2500, 0.9583, 0.0016, 0.0723, 0.0023},
        {0.1000, 0.2500, 0.2500, 0.4500, 0.9439, 0.0020, 0.0812, 0.0030},
        {0.1000, 0.2500, 0.4500, 0.6500, 0.9348, 0.0027, 0.0732, 0.0039},
        {0.1000, 0.2500, 0.6500, 0.8500, 0.9296, 0.0037, 0.0807, 0.0055},
        {0.1000, 0.2500, 0.8500, 1.0500, 0.9292, 0.0049, 0.0805, 0.0102},
        {0.1000, 0.2500, 1.0500, 1.2500, 0.9221, 0.0051, 0.0738, 0.0069},

        {0.2500, 0.3500, 0.0500, 0.2500, 0.9500, 0.0000, 0.0723, 0.0008},
        {0.2500, 0.3500, 0.2500, 0.4500, 0.9340, 0.0013, 0.0803, 0.0020},
        {0.2500, 0.3500, 0.4500, 0.6500, 0.9263, 0.0018, 0.0815, 0.0000},
        {0.2500, 0.3500, 0.6500, 0.8500, 0.9151, 0.0025, 0.0761, 0.0038},
        {0.2500, 0.3500, 0.8500, 1.0500, 0.9052, 0.0033, 0.0734, 0.0047},
        {0.2500, 0.3500, 1.0500, 1.2500, 0.9054, 0.0051, 0.0769, 0.0076},

        {0.3500, 0.4500, 0.0500, 0.2500, 0.9377, 0.0006, 0.0787, 0.0009},
        {0.3500, 0.4500, 0.2500, 0.4500, 0.9159, 0.0014, 0.0703, 0.0021},
        {0.3500, 0.4500, 0.4500, 0.6500, 0.9139, 0.0018, 0.0790, 0.0026},
        {0.3500, 0.4500, 0.6500, 0.8500, 0.9066, 0.0024, 0.0821, 0.0038},
        {0.3500, 0.4500, 0.8500, 1.0500, 0.8983, 0.0029, 0.0740, 0.0045},
        {0.3500, 0.4500, 1.0500, 1.2500, 0.9079, 0.0042, 0.0847, 0.0066},

        {0.4500, 0.6000, 0.0500, 0.2500, 0.9166, 0.0007, 0.0849, 0.0010},
        {0.4500, 0.6000, 0.2500, 0.4500, 0.9067, 0.0026, 0.0800, 0.0002},
        {0.4500, 0.6000, 0.4500, 0.6500, 0.9114, 0.0022, 0.0790, 0.0033},
        {0.4500, 0.6000, 0.6500, 0.8500, 0.9023, 0.0026, 0.0763, 0.0037},
        {0.4500, 0.6000, 0.8500, 1.0500, 0.8977, 0.0028, 0.0740, 0.0041},
        {0.4500, 0.6000, 1.0500, 1.2500, 0.8979, 0.0032, 0.0667, 0.0047}
    };
    return pts;
}

void write_csv(const std::vector<FitResult>& results, const std::string& path) {
    std::ofstream out(path.c_str());
    out << "slice_name,x_min,x_max,t_min,t_max,mu,mu_err,sigma,sigma_err,n_entries,fit_success,fit_status,cov_status,chi2,ndf\n";
    for (const auto& r : results) {
        out << r.slice_name << ","
            << std::fixed << std::setprecision(4) << r.x_min << ","
            << std::fixed << std::setprecision(4) << r.x_max << ","
            << std::fixed << std::setprecision(4) << r.t_min << ","
            << std::fixed << std::setprecision(4) << r.t_max << ",";

        if (std::isfinite(r.mu)) {
            out << std::fixed << std::setprecision(6) << r.mu << ","
                << std::fixed << std::setprecision(6) << r.mu_err << ","
                << std::fixed << std::setprecision(6) << r.sigma << ","
                << std::fixed << std::setprecision(6) << r.sigma_err << ",";
        } else {
            out << ",,,,";
        }

        out << r.n_entries << ","
            << (r.fit_success ? 1 : 0) << ","
            << r.fit_status << ","
            << r.cov_status << ",";

        if (std::isfinite(r.chi2)) {
            out << std::fixed << std::setprecision(6) << r.chi2 << ",";
        } else {
            out << ",";
        }

        out << r.ndf << "\n";
    }
}

void write_latex_table(const std::vector<FitResult>& results, const std::string& path) {
    std::ofstream out(path.c_str());

    out << "\\begin{table}[htbp]\n";
    out << "  \\centering\n";
    out << "  \\caption{Fit results in $(x_B,-t^\\prime)$ bins. Shown are the fitted mean $\\mu$, its uncertainty $\\delta\\mu$, the fitted width $\\sigma$, and its uncertainty $\\delta\\sigma$.}\n";
    out << "  \\label{tab:mx2_tprime_fit_results}\n";
    out << "  \\begin{tabular}{cccccccc}\n";
    out << "    \\hline\n";
    out << "    $x_{B,\\min}$ & $x_{B,\\max}$ & $(-t^\\prime)_{\\min}$ & $(-t^\\prime)_{\\max}$ & $\\mu$ & $\\delta\\mu$ & $\\sigma$ & $\\delta\\sigma$ \\\\\n";
    out << "    \\hline\n";

    std::string previous_slice = "";
    for (const auto& r : results) {
        if (!previous_slice.empty() && r.slice_name != previous_slice) {
            out << "    \\hline\n";
        }

        out << "    "
            << std::fixed << std::setprecision(4) << r.x_min << " & "
            << std::fixed << std::setprecision(4) << r.x_max << " & "
            << std::fixed << std::setprecision(4) << r.t_min << " & "
            << std::fixed << std::setprecision(4) << r.t_max << " & "
            << format_double(r.mu, 4) << " & "
            << format_double(r.mu_err, 4) << " & "
            << format_double(r.sigma, 4) << " & "
            << format_double(r.sigma_err, 4) << " \\\\\n";

        previous_slice = r.slice_name;
    }

    out << "    \\hline\n";
    out << "  \\end{tabular}\n";
    out << "\\end{table}\n";
}

void draw_fit_grid(const std::vector<std::vector<TH1D*> >& hists, const std::vector<FitResult>& results, const std::string& png_path, const std::string& pdf_path) {
    TCanvas canvas("c_grid", "Mx2 grid", 2400, 1400);
    canvas.Divide(6, 4, 0.001, 0.001);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.040);

    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            int ipad = ix * static_cast<int>(kTPrimeBins.size()) + it + 1;
            canvas.cd(ipad);

            gPad->SetGrid();
            gPad->SetLeftMargin(0.13);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin(0.10);
            gPad->SetBottomMargin(0.12);

            TH1D* h = hists[ix][it];
            const FitResult& r = results[ix * static_cast<int>(kTPrimeBins.size()) + it];

            h->SetLineColor(kBlack);
            h->SetMarkerColor(kBlack);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.7);
            h->SetLineWidth(1);

            h->GetXaxis()->SetTitle("M_{X}^{2} (GeV^{2})");
            h->GetYaxis()->SetTitle("Counts");
            h->GetXaxis()->CenterTitle();
            h->GetYaxis()->CenterTitle();
            h->GetXaxis()->SetTitleSize(0.050);
            h->GetYaxis()->SetTitleSize(0.050);
            h->GetXaxis()->SetLabelSize(0.040);
            h->GetYaxis()->SetLabelSize(0.040);

            h->Draw("E1");

            std::ostringstream title;
            title << "x_{B} #in [" << std::fixed << std::setprecision(2)
                  << r.x_min << ", " << r.x_max << "], -t' #in ["
                  << r.t_min << ", " << r.t_max << "]";
            latex.DrawLatex(0.12, 0.92, title.str().c_str());

            TLegend leg(0.18, 0.13, 0.82, 0.28);
            leg.SetBorderSize(1);
            leg.SetFillStyle(1001);
            leg.SetFillColor(kWhite);
            leg.SetTextSize(0.038);

            TF1* f = h->GetFunction(make_fit_name(ix, it).c_str());
            if (f) {
                f->SetLineColor(kRed + 1);
                f->SetLineStyle(2);
                f->SetLineWidth(2);
                f->Draw("same");
            }

            if (std::isfinite(r.mu) && std::isfinite(r.sigma)) {
                std::ostringstream fit_label;
                fit_label << "#mu = " << std::fixed << std::setprecision(4) << r.mu << " (GeV^{2})"
                          << ", #sigma = " << std::fixed << std::setprecision(4) << r.sigma << " (GeV^{2})";
                leg.AddEntry(h, "Data", "lep");
                leg.AddEntry((TObject*)0, fit_label.str().c_str(), "");
            } else {
                leg.AddEntry(h, "Data", "lep");
                leg.AddEntry((TObject*)0, "No extracted values", "");
            }

            leg.Draw();
        }
    }

    canvas.SaveAs(png_path.c_str());
    canvas.SaveAs(pdf_path.c_str());
}

void draw_summary_plot(const std::vector<FitResult>& results, const std::string& png_path, const std::string& pdf_path) {
    TCanvas canvas("c_summary", "mu sigma summary", 1600, 1000);
    canvas.Divide(1, 2, 0.001, 0.001);

    int n_total = static_cast<int>(results.size());
    int colors[4] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2};

    canvas.cd(1);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.12);

    TH1D frame_mu("frame_mu", "", n_total, 0.5, n_total + 0.5);
    frame_mu.SetMinimum(0.84);
    frame_mu.SetMaximum(1.03);
    frame_mu.GetXaxis()->SetTitle("Bin");
    frame_mu.GetYaxis()->SetTitle("#mu (GeV^{2})");
    frame_mu.GetXaxis()->CenterTitle();
    frame_mu.GetYaxis()->CenterTitle();
    frame_mu.GetYaxis()->SetTitleSize(0.050);
    frame_mu.GetYaxis()->SetLabelSize(0.040);
    frame_mu.GetYaxis()->SetTitleOffset(0.95);
    frame_mu.Draw();

    TLegend leg_mu(0.12, 0.77, 0.52, 0.93);
    leg_mu.SetNColumns(2);
    leg_mu.SetBorderSize(1);
    leg_mu.SetFillStyle(1001);
    leg_mu.SetFillColor(kWhite);
    leg_mu.SetTextSize(0.035);

    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        std::vector<double> xx, yy, ex, ey;
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            int idx = ix * static_cast<int>(kTPrimeBins.size()) + it;
            if (!std::isfinite(results[idx].mu) || !std::isfinite(results[idx].mu_err)) continue;
            xx.push_back(idx + 1);
            ex.push_back(0.0);
            yy.push_back(results[idx].mu);
            ey.push_back(results[idx].mu_err);
        }

        if (xx.empty()) continue;

        TGraphErrors* g = new TGraphErrors(static_cast<int>(xx.size()), xx.data(), yy.data(), ex.data(), ey.data());
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);
        g->SetMarkerColor(colors[ix]);
        g->SetLineColor(colors[ix]);
        g->SetLineWidth(2);
        g->Draw("P SAME");

        std::ostringstream label;
        label << "x_{B} #in [" << std::fixed << std::setprecision(2)
              << kXBSlices[ix].xmin << ", " << kXBSlices[ix].xmax << "]";
        leg_mu.AddEntry(g, label.str().c_str(), "lp");
    }

    for (int xline = 6; xline <= 18; xline += 6) {
        TLine* line = new TLine(xline + 0.5, frame_mu.GetMinimum(), xline + 0.5, frame_mu.GetMaximum());
        line->SetLineStyle(3);
        line->SetLineColor(kGray + 1);
        line->Draw("same");
    }

    leg_mu.Draw();

    canvas.cd(2);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);

    TH1D frame_sigma("frame_sigma", "", n_total, 0.5, n_total + 0.5);
    frame_sigma.SetMinimum(0.04);
    frame_sigma.SetMaximum(0.16);
    frame_sigma.GetXaxis()->SetTitle("Bin");
    frame_sigma.GetYaxis()->SetTitle("#sigma (GeV^{2})");
    frame_sigma.GetXaxis()->CenterTitle();
    frame_sigma.GetYaxis()->CenterTitle();
    frame_sigma.GetYaxis()->SetTitleSize(0.050);
    frame_sigma.GetYaxis()->SetLabelSize(0.040);
    frame_sigma.GetYaxis()->SetTitleOffset(0.95);
    frame_sigma.Draw();

    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        std::vector<double> xx, yy, ex, ey;
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            int idx = ix * static_cast<int>(kTPrimeBins.size()) + it;
            if (!std::isfinite(results[idx].sigma) || !std::isfinite(results[idx].sigma_err)) continue;
            xx.push_back(idx + 1);
            ex.push_back(0.0);
            yy.push_back(results[idx].sigma);
            ey.push_back(results[idx].sigma_err);
        }

        if (xx.empty()) continue;

        TGraphErrors* g = new TGraphErrors(static_cast<int>(xx.size()), xx.data(), yy.data(), ex.data(), ey.data());
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);
        g->SetMarkerColor(colors[ix]);
        g->SetLineColor(colors[ix]);
        g->SetLineWidth(2);
        g->Draw("P SAME");
    }

    for (int xline = 6; xline <= 18; xline += 6) {
        TLine* line = new TLine(xline + 0.5, frame_sigma.GetMinimum(), xline + 0.5, frame_sigma.GetMaximum());
        line->SetLineStyle(3);
        line->SetLineColor(kGray + 1);
        line->Draw("same");
    }

    canvas.SaveAs(png_path.c_str());
    canvas.SaveAs(pdf_path.c_str());
}

void draw_summary_plot_hardcoded(const std::vector<SummaryPoint>& points, const std::string& png_path, const std::string& pdf_path) {
    TCanvas canvas("c_summary_hardcoded", "mu sigma summary hardcoded", 1600, 1000);
    canvas.Divide(1, 2, 0.001, 0.001);

    int n_total = static_cast<int>(points.size());
    int colors[4] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2};

    double sigma_max_data = 0.0;
    for (size_t i = 0; i < points.size(); i++) {
        sigma_max_data = std::max(sigma_max_data, points[i].sigma + points[i].sigma_err);
    }
    double sigma_plot_max = std::max(0.16, 1.05 * sigma_max_data);

    canvas.cd(1);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.12);

    TH1D frame_mu("frame_mu_hardcoded", "", n_total, 0.5, n_total + 0.5);
    frame_mu.SetMinimum(0.84);
    frame_mu.SetMaximum(1.03);
    frame_mu.GetXaxis()->SetTitle("Bin");
    frame_mu.GetYaxis()->SetTitle("#mu");
    frame_mu.GetXaxis()->CenterTitle();
    frame_mu.GetYaxis()->CenterTitle();
    frame_mu.GetYaxis()->SetTitleSize(0.050);
    frame_mu.GetYaxis()->SetLabelSize(0.040);
    frame_mu.GetYaxis()->SetTitleOffset(0.95);
    frame_mu.Draw();

    TLegend leg_mu(0.12, 0.77, 0.52, 0.93);
    leg_mu.SetNColumns(2);
    leg_mu.SetBorderSize(1);
    leg_mu.SetFillStyle(1001);
    leg_mu.SetFillColor(kWhite);
    leg_mu.SetTextSize(0.035);

    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        std::vector<double> xx, yy, ex, ey;
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            int idx = ix * static_cast<int>(kTPrimeBins.size()) + it;
            xx.push_back(idx + 1);
            ex.push_back(0.0);
            yy.push_back(points[idx].mu);
            ey.push_back(points[idx].mu_err);
        }

        TGraphErrors* g = new TGraphErrors(static_cast<int>(xx.size()), xx.data(), yy.data(), ex.data(), ey.data());
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);
        g->SetMarkerColor(colors[ix]);
        g->SetLineColor(colors[ix]);
        g->SetLineWidth(2);
        g->Draw("P SAME");

        std::ostringstream label;
        label << "x_{B} #in [" << std::fixed << std::setprecision(2)
              << kXBSlices[ix].xmin << ", " << kXBSlices[ix].xmax << "]";
        leg_mu.AddEntry(g, label.str().c_str(), "lp");
    }

    for (int xline = 6; xline <= 18; xline += 6) {
        TLine* line = new TLine(xline + 0.5, frame_mu.GetMinimum(), xline + 0.5, frame_mu.GetMaximum());
        line->SetLineStyle(3);
        line->SetLineColor(kGray + 1);
        line->Draw("same");
    }

    leg_mu.Draw();

    canvas.cd(2);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);

    TH1D frame_sigma("frame_sigma_hardcoded", "", n_total, 0.5, n_total + 0.5);
    frame_sigma.SetMinimum(0.04);
    frame_sigma.SetMaximum(sigma_plot_max);
    frame_sigma.GetXaxis()->SetTitle("Bin");
    frame_sigma.GetYaxis()->SetTitle("#sigma");
    frame_sigma.GetXaxis()->CenterTitle();
    frame_sigma.GetYaxis()->CenterTitle();
    frame_sigma.GetYaxis()->SetTitleSize(0.050);
    frame_sigma.GetYaxis()->SetLabelSize(0.040);
    frame_sigma.GetYaxis()->SetTitleOffset(0.95);
    frame_sigma.Draw();

    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        std::vector<double> xx, yy, ex, ey;
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            int idx = ix * static_cast<int>(kTPrimeBins.size()) + it;
            xx.push_back(idx + 1);
            ex.push_back(0.0);
            yy.push_back(points[idx].sigma);
            ey.push_back(points[idx].sigma_err);
        }

        TGraphErrors* g = new TGraphErrors(static_cast<int>(xx.size()), xx.data(), yy.data(), ex.data(), ey.data());
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);
        g->SetMarkerColor(colors[ix]);
        g->SetLineColor(colors[ix]);
        g->SetLineWidth(2);
        g->Draw("P SAME");
    }

    for (int xline = 6; xline <= 18; xline += 6) {
        TLine* line = new TLine(xline + 0.5, frame_sigma.GetMinimum(), xline + 0.5, frame_sigma.GetMaximum());
        line->SetLineStyle(3);
        line->SetLineColor(kGray + 1);
        line->Draw("same");
    }

    canvas.SaveAs(png_path.c_str());
    canvas.SaveAs(pdf_path.c_str());
}

int main() {
    setup_global_style();
    ensure_output_dir(kOutputDir);

    std::cout << "Opening file: " << kInputFile << std::endl;
    TFile* file = TFile::Open(kInputFile.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: could not open input file." << std::endl;
        return 1;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get(kTreeName.c_str()));
    if (!tree) {
        std::cerr << "ERROR: could not find tree " << kTreeName << std::endl;
        file->Close();
        return 1;
    }

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("x", 1);
    tree->SetBranchStatus("tprime", 1);
    tree->SetBranchStatus("Mx2", 1);

    std::vector<std::vector<TH1D*> > hists(
        kXBSlices.size(),
        std::vector<TH1D*>(kTPrimeBins.size(), nullptr)
    );

    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            std::string hname = make_hist_name(ix, it);
            std::ostringstream htitle;
            htitle << ";M_{X}^{2} (GeV^{2});Counts";
            hists[ix][it] = new TH1D(hname.c_str(), htitle.str().c_str(), kNMx2Bins, kMx2Min, kMx2Max);
            hists[ix][it]->Sumw2();
        }
    }

    std::cout << "Beginning single event loop..." << std::endl;

    TTreeReader reader(tree);
    TTreeReaderValue<double> x(reader, "x");
    TTreeReaderValue<double> tprime(reader, "tprime");
    TTreeReaderValue<double> mx2(reader, "Mx2");

    long long n_processed = 0;
    long long n_filled = 0;

    while (reader.Next()) {
        n_processed++;

        if (!std::isfinite(*x) || !std::isfinite(*tprime) || !std::isfinite(*mx2)) {
            continue;
        }

        double minus_tprime = -(*tprime);

        int ix = find_xbin(*x);
        if (ix < 0) continue;

        int it = find_tbin(minus_tprime);
        if (it < 0) continue;

        hists[ix][it]->Fill(*mx2);
        n_filled++;
    }

    std::cout << "Processed events: " << n_processed << std::endl;
    std::cout << "Filled events:    " << n_filled << std::endl;

    std::vector<FitResult> results;
    results.reserve(kXBSlices.size() * kTPrimeBins.size());

    std::cout << "Fitting histograms..." << std::endl;
    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            FitResult r = fit_histogram(hists[ix][it], ix, it);
            results.push_back(r);

            std::cout
                << std::setw(8) << r.slice_name
                << "  x=[" << std::fixed << std::setprecision(2) << r.x_min << "," << r.x_max << "]"
                << "  -tprime=[" << std::fixed << std::setprecision(2) << r.t_min << "," << r.t_max << "]"
                << "  N=" << r.n_entries
                << "  fit_success=" << (r.fit_success ? 1 : 0)
                << "  status=" << r.fit_status
                << "  cov=" << r.cov_status
                << "  mu=" << format_double(r.mu, 6)
                << "  mu_err=" << format_double(r.mu_err, 6)
                << "  sigma=" << format_double(r.sigma, 6)
                << "  sigma_err=" << format_double(r.sigma_err, 6)
                << std::endl;
        }
    }

    std::vector<SummaryPoint> hardcoded_points = get_hardcoded_summary_points();

    std::string csv_path = kOutputDir + "/Mx2_fit_params_tprime.csv";
    std::string tex_path = kOutputDir + "/Mx2_fit_params_tprime_table.tex";
    std::string grid_png = kOutputDir + "/Mx2_tprime_fit_grid.png";
    std::string grid_pdf = kOutputDir + "/Mx2_tprime_fit_grid.pdf";
    std::string summary_png = kOutputDir + "/Mx2_tprime_mu_sigma_summary.png";
    std::string summary_pdf = kOutputDir + "/Mx2_tprime_mu_sigma_summary.pdf";
    std::string summary_hardcoded_png = kOutputDir + "/Mx2_tprime_mu_sigma_summary_hardcoded.png";
    std::string summary_hardcoded_pdf = kOutputDir + "/Mx2_tprime_mu_sigma_summary_hardcoded.pdf";

    std::cout << "Writing CSV..." << std::endl;
    write_csv(results, csv_path);

    std::cout << "Writing LaTeX table..." << std::endl;
    write_latex_table(results, tex_path);

    std::cout << "Drawing fit grid..." << std::endl;
    draw_fit_grid(hists, results, grid_png, grid_pdf);

    std::cout << "Drawing fit-derived summary plot..." << std::endl;
    draw_summary_plot(results, summary_png, summary_pdf);

    std::cout << "Drawing hard-coded summary plot..." << std::endl;
    draw_summary_plot_hardcoded(hardcoded_points, summary_hardcoded_png, summary_hardcoded_pdf);

    std::cout << "Done." << std::endl;
    std::cout << "CSV:                   " << csv_path << std::endl;
    std::cout << "LaTeX table:           " << tex_path << std::endl;
    std::cout << "Grid PNG:              " << grid_png << std::endl;
    std::cout << "Grid PDF:              " << grid_pdf << std::endl;
    std::cout << "Summary PNG:           " << summary_png << std::endl;
    std::cout << "Summary PDF:           " << summary_pdf << std::endl;
    std::cout << "Hard-coded Summary PNG:" << summary_hardcoded_png << std::endl;
    std::cout << "Hard-coded Summary PDF:" << summary_hardcoded_pdf << std::endl;

    file->Close();
    return 0;
}