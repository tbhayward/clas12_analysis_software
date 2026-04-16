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

static const std::string kInputFile = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_combined_inb_NH3_epi+_2.root";
static const std::string kTreeName = "PhysicsEvents";
static const std::string kOutputDir = "output/enpi+_Mx2_fits";

static constexpr double kMx2Min = 0.70;
static constexpr double kMx2Max = 1.10;
static constexpr int kNMx2Bins = 60;

static constexpr double kMuMin = 0.85;
static constexpr double kMuMax = 0.95;
static constexpr double kSigmaMin = 0.03;
static constexpr double kSigmaMax = 0.08;

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

FitResult fit_histogram(TH1D* hist, int ix, int it) {
    FitResult result;
    result.slice_name = kXBSlices[ix].name;
    result.x_min = kXBSlices[ix].xmin;
    result.x_max = kXBSlices[ix].xmax;
    result.t_min = kTPrimeBins[it].tmin;
    result.t_max = kTPrimeBins[it].tmax;
    result.n_entries = static_cast<long long>(hist->GetEntries());

    int peak_bin = hist->GetMaximumBin();
    double peak_x = hist->GetBinCenter(peak_bin);

    double edge_background = 0.0;
    int n_edge = 0;
    for (int i = 1; i <= 5; i++) {
        edge_background += hist->GetBinContent(i);
        n_edge++;
    }
    for (int i = kNMx2Bins - 4; i <= kNMx2Bins; i++) {
        edge_background += hist->GetBinContent(i);
        n_edge++;
    }
    if (n_edge > 0) {
        edge_background /= static_cast<double>(n_edge);
    }

    double amplitude_guess = std::max(hist->GetMaximum() - edge_background, 1.0);
    double mu_guess = std::max(kMuMin + 0.002, std::min(peak_x, kMuMax - 0.002));
    double sigma_guess = 0.07;

    std::string fname = make_fit_name(ix, it);
    TF1 fit_fn(fname.c_str(), "gaus(0) + pol2(3)", kMx2Min, kMx2Max);
    fit_fn.SetLineColor(kRed + 1);
    fit_fn.SetLineStyle(2);
    fit_fn.SetLineWidth(2);

    fit_fn.SetParameter(0, amplitude_guess);
    fit_fn.SetParameter(1, mu_guess);
    fit_fn.SetParameter(2, sigma_guess);
    fit_fn.SetParameter(3, edge_background);
    fit_fn.SetParameter(4, 0.0);
    fit_fn.SetParameter(5, 0.0);

    fit_fn.SetParLimits(0, 0.0, 1.0e12);
    fit_fn.SetParLimits(1, kMuMin, kMuMax);
    fit_fn.SetParLimits(2, kSigmaMin, kSigmaMax);

    TFitResultPtr fit_ptr = hist->Fit(&fit_fn, "QRS");
    int fit_status = static_cast<int>(fit_ptr);
    int cov_status = -1;
    if (fit_ptr.Get()) {
        cov_status = fit_ptr->CovMatrixStatus();
    }

    result.fit_status = fit_status;
    result.cov_status = cov_status;
    result.fit_success = (fit_status == 0);

    if (result.fit_success) {
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
    }

    return result;
}

void write_csv(const std::vector<FitResult>& results, const std::string& path) {
    std::ofstream out(path.c_str());
    out << "slice_name,x_min,x_max,t_min,t_max,mu,mu_err,sigma,sigma_err,n_entries,fit_success\n";
    for (const auto& r : results) {
        out << r.slice_name << ","
            << std::fixed << std::setprecision(4) << r.x_min << ","
            << std::fixed << std::setprecision(4) << r.x_max << ","
            << std::fixed << std::setprecision(4) << r.t_min << ","
            << std::fixed << std::setprecision(4) << r.t_max << ",";
        if (r.fit_success) {
            out << std::fixed << std::setprecision(6) << r.mu << ","
                << std::fixed << std::setprecision(6) << r.mu_err << ","
                << std::fixed << std::setprecision(6) << r.sigma << ","
                << std::fixed << std::setprecision(6) << r.sigma_err << ",";
        } else {
            out << ",,,,";
        }
        out << r.n_entries << ","
            << (r.fit_success ? 1 : 0) << "\n";
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

void draw_fit_grid(const std::vector<std::vector<TH1D*>>& hists, const std::vector<FitResult>& results, const std::string& png_path, const std::string& pdf_path) {
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

            if (r.fit_success) {
                TF1* f = h->GetFunction(make_fit_name(ix, it).c_str());
                if (f) {
                    f->SetLineColor(kRed + 1);
                    f->SetLineStyle(2);
                    f->SetLineWidth(2);
                    f->Draw("same");
                }

                std::ostringstream fit_label;
                fit_label << "#mu = " << std::fixed << std::setprecision(4) << r.mu
                          << ", #sigma = " << std::fixed << std::setprecision(4) << r.sigma;
                leg.AddEntry(h, "Data", "lep");
                leg.AddEntry((TObject*)0, fit_label.str().c_str(), "");
            } else {
                leg.AddEntry(h, "Data", "lep");
                leg.AddEntry((TObject*)0, "Fit failed", "");
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
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.12);

    TH1D frame_mu("frame_mu", "", n_total, 0.5, n_total + 0.5);
    frame_mu.SetMinimum(0.84);
    frame_mu.SetMaximum(0.97);
    frame_mu.GetXaxis()->SetTitle("Bin");
    frame_mu.GetYaxis()->SetTitle("#mu");
    frame_mu.GetXaxis()->CenterTitle();
    frame_mu.GetYaxis()->CenterTitle();
    frame_mu.Draw();

    TLegend leg_mu(0.12, 0.77, 0.52, 0.93);
    leg_mu.SetNColumns(2);
    leg_mu.SetBorderSize(1);
    leg_mu.SetFillStyle(1001);
    leg_mu.SetFillColor(kWhite);
    leg_mu.SetTextSize(0.035);

    std::vector<TGraphErrors*> mu_graphs;
    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        std::vector<double> xx, yy, ex, ey;
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            int idx = ix * static_cast<int>(kTPrimeBins.size()) + it;
            xx.push_back(idx + 1);
            ex.push_back(0.0);
            yy.push_back(results[idx].mu);
            ey.push_back(results[idx].mu_err);
        }

        TGraphErrors* g = new TGraphErrors(static_cast<int>(xx.size()), xx.data(), yy.data(), ex.data(), ey.data());
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);
        g->SetMarkerColor(colors[ix]);
        g->SetLineColor(colors[ix]);
        g->SetLineWidth(2);
        g->Draw("P SAME");
        mu_graphs.push_back(g);

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
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);

    TH1D frame_sigma("frame_sigma", "", n_total, 0.5, n_total + 0.5);
    frame_sigma.SetMinimum(0.04);
    frame_sigma.SetMaximum(0.09);
    frame_sigma.GetXaxis()->SetTitle("Bin");
    frame_sigma.GetYaxis()->SetTitle("#sigma");
    frame_sigma.GetXaxis()->CenterTitle();
    frame_sigma.GetYaxis()->CenterTitle();
    frame_sigma.Draw();

    std::vector<TGraphErrors*> sigma_graphs;
    for (int ix = 0; ix < static_cast<int>(kXBSlices.size()); ix++) {
        std::vector<double> xx, yy, ex, ey;
        for (int it = 0; it < static_cast<int>(kTPrimeBins.size()); it++) {
            int idx = ix * static_cast<int>(kTPrimeBins.size()) + it;
            xx.push_back(idx + 1);
            ex.push_back(0.0);
            yy.push_back(results[idx].sigma);
            ey.push_back(results[idx].sigma_err);
        }

        TGraphErrors* g = new TGraphErrors(static_cast<int>(xx.size()), xx.data(), yy.data(), ex.data(), ey.data());
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.1);
        g->SetMarkerColor(colors[ix]);
        g->SetLineColor(colors[ix]);
        g->SetLineWidth(2);
        g->Draw("P SAME");
        sigma_graphs.push_back(g);
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
                << "  mu=" << format_double(r.mu, 6)
                << "  mu_err=" << format_double(r.mu_err, 6)
                << "  sigma=" << format_double(r.sigma, 6)
                << "  sigma_err=" << format_double(r.sigma_err, 6)
                << std::endl;
        }
    }

    std::string csv_path = kOutputDir + "/Mx2_fit_params_tprime.csv";
    std::string tex_path = kOutputDir + "/Mx2_fit_params_tprime_table.tex";
    std::string grid_png = kOutputDir + "/Mx2_tprime_fit_grid.png";
    std::string grid_pdf = kOutputDir + "/Mx2_tprime_fit_grid.pdf";
    std::string summary_png = kOutputDir + "/Mx2_tprime_mu_sigma_summary.png";
    std::string summary_pdf = kOutputDir + "/Mx2_tprime_mu_sigma_summary.pdf";

    std::cout << "Writing CSV..." << std::endl;
    write_csv(results, csv_path);

    std::cout << "Writing LaTeX table..." << std::endl;
    write_latex_table(results, tex_path);

    std::cout << "Drawing fit grid..." << std::endl;
    draw_fit_grid(hists, results, grid_png, grid_pdf);

    std::cout << "Drawing summary plot..." << std::endl;
    draw_summary_plot(results, summary_png, summary_pdf);

    std::cout << "Done." << std::endl;
    std::cout << "CSV:         " << csv_path << std::endl;
    std::cout << "LaTeX table: " << tex_path << std::endl;
    std::cout << "Grid PNG:    " << grid_png << std::endl;
    std::cout << "Grid PDF:    " << grid_pdf << std::endl;
    std::cout << "Summary PNG: " << summary_png << std::endl;
    std::cout << "Summary PDF: " << summary_pdf << std::endl;

    file->Close();
    return 0;
}