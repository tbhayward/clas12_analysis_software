// dvcsgen_ycol_diagnostic.cpp
// -----------------------------------------------------------------------------
// Determine dvcsgen --ycol semantics by checking which variable has a hard edge
// at 0.005 on a sample generated with --ycol 0.005.
//
// This version processes at most the first 5,000,000 events (or fewer if the
// tree has fewer entries).
//
// Produces a 2-pad canvas:
//   left:  P1_neg = - (k - qgamma)^{2} / Q^{2}   (dimensionless)
//   right: ycol_minus_y = ycol - y              (dimensionless)
//
// Saves:
//   /u/home/thayward/dvcsgen_ycol_diagnostic.png
//   /u/home/thayward/dvcsgen_ycol_diagnostic.pdf
//
// Branches required (angles in radians):
//   e_p, e_theta, e_phi
//   p2_p, p2_theta, p2_phi     (p2 = photon)
//   x, Q2, t1, y
//
// Tree assumed: PhysicsEvents
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>

static void fatal(const std::string &msg) {
    std::cerr << "FATAL: " << msg << std::endl;
    std::exit(1);
}

static void require_branch(TTree *t, const char *name) {
    if (!t) fatal("require_branch: null tree");
    if (!t->GetBranch(name)) {
        fatal(std::string("missing required branch: ") + name);
    }
}

static TLorentzVector make_beam_electron(double Ebeam_GeV) {
    const double me = 0.00051099895; // (GeV)
    if (Ebeam_GeV <= me) {
        fatal("beam energy <= electron mass (GeV)");
    }

    const double pz = std::sqrt(Ebeam_GeV * Ebeam_GeV - me * me);
    return TLorentzVector(0.0, 0.0, pz, Ebeam_GeV);
}

static TLorentzVector make_p4_from_p_theta_phi(double p,
                                               double theta_rad,
                                               double phi_rad,
                                               double mass_GeV) {
    if (!(p >= 0.0)) {
        fatal("make_p4_from_p_theta_phi: p < 0");
    }

    const double px = p * std::sin(theta_rad) * std::cos(phi_rad);
    const double py = p * std::sin(theta_rad) * std::sin(phi_rad);
    const double pz = p * std::cos(theta_rad);
    const double E  = std::sqrt(p * p + mass_GeV * mass_GeV);

    return TLorentzVector(px, py, pz, E);
}

// ycol = (Q2 - |t|) / (Q2 - x*|t|) where x is Bjorken x (0<x<1), |t| = |t1|
static double compute_ycol(double x, double Q2, double t1) {
    const double tabs = std::fabs(t1);

    const double denom = Q2 - x * tabs;
    if (!(denom > 0.0)) {
        // Inconsistent kinematics for DVCS definitions; treat as failing.
        return -1.0;
    }

    return (Q2 - tabs) / denom;
}

int main(int argc, char **argv) {
    const std::string in_file = (argc > 1)
        ? std::string(argv[1])
        : std::string("/volatile/clas12/thayward/dvcs_event_builder_test/gen_dvcsgen_sp18_out.root");

    const double Ebeam = (argc > 2)
        ? std::atof(argv[2])
        : 10.594;

    const std::string out_dir = (argc > 3)
        ? std::string(argv[3])
        : std::string("/u/home/thayward");

    const std::string tree_name = (argc > 4)
        ? std::string(argv[4])
        : std::string("PhysicsEvents");

    const double ycol_cut = 0.005;

    // Hard cap for speed
    const Long64_t max_events = 5000000;

    TFile *f = TFile::Open(in_file.c_str(), "READ");
    if (!f || f->IsZombie()) {
        fatal("cannot open input ROOT file: " + in_file);
    }

    TTree *t = dynamic_cast<TTree *>(f->Get(tree_name.c_str()));
    if (!t) {
        fatal("cannot find TTree '" + tree_name + "' in file: " + in_file);
    }

    // Required branches (explicit; fail fast if names mismatch)
    require_branch(t, "e_p");
    require_branch(t, "e_theta");
    require_branch(t, "e_phi");
    require_branch(t, "p2_p");
    require_branch(t, "p2_theta");
    require_branch(t, "p2_phi");
    require_branch(t, "x");
    require_branch(t, "Q2");
    require_branch(t, "t1");
    require_branch(t, "y");

    double e_p = 0.0, e_theta = 0.0, e_phi = 0.0;
    double p2_p = 0.0, p2_theta = 0.0, p2_phi = 0.0;
    double x = 0.0, Q2 = 0.0, t1 = 0.0, y = 0.0;

    t->SetBranchAddress("e_p", &e_p);
    t->SetBranchAddress("e_theta", &e_theta);
    t->SetBranchAddress("e_phi", &e_phi);

    t->SetBranchAddress("p2_p", &p2_p);
    t->SetBranchAddress("p2_theta", &p2_theta);
    t->SetBranchAddress("p2_phi", &p2_phi);

    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("Q2", &Q2);
    t->SetBranchAddress("t1", &t1);
    t->SetBranchAddress("y", &y);

    // Histograms
    TH1D *h_P1neg = new TH1D("h_P1neg",
                            "Candidate A: P1_neg = - (k - qgamma)^{2} / Q^{2};P1_neg;Counts",
                            400, 0.0, 0.08);

    TH1D *h_ycolmy = new TH1D("h_ycolmy",
                             "Candidate B: ycol - y; ycol - y;Counts",
                             400, -0.02, 0.08);

    TH1D *h_P1abs  = new TH1D("h_P1abs",
                             "P1_abs = |(k - qgamma)^{2}/Q^{2}|;P1_abs;Counts",
                             400, 0.0, 0.08);

    const TLorentzVector k = make_beam_electron(Ebeam);
    const double me = 0.00051099895; // (GeV)

    const Long64_t nentries_total = t->GetEntries();
    const Long64_t nentries = (nentries_total < max_events) ? nentries_total : max_events;

    Long64_t n_bad_Q2 = 0;
    Long64_t n_bad_ycol = 0;

    Long64_t n_P1neg_below = 0;
    Long64_t n_ycolmy_below = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        if (!(Q2 > 0.0) || !std::isfinite(Q2)) {
            n_bad_Q2++;
            continue;
        }

        const TLorentzVector kprime = make_p4_from_p_theta_phi(e_p,  e_theta,  e_phi,  me);
        const TLorentzVector qgamma = make_p4_from_p_theta_phi(p2_p, p2_theta, p2_phi, 0.0);

        const double P1_scaled = (k - qgamma).M2() / Q2;
        const double P1_neg    = -P1_scaled;
        const double P1_abs    = std::fabs(P1_scaled);

        const double ycol = compute_ycol(x, Q2, t1);
        if (!(ycol > -0.5)) {
            n_bad_ycol++;
            continue;
        }
        const double ycol_minus_y = ycol - y;

        h_P1neg->Fill(P1_neg);
        h_P1abs->Fill(P1_abs);
        h_ycolmy->Fill(ycol_minus_y);

        if (P1_neg < ycol_cut) n_P1neg_below++;
        if (ycol_minus_y < ycol_cut) n_ycolmy_below++;
    } //endfor

    std::cout << "Input file: " << in_file << std::endl;
    std::cout << "Tree:       " << tree_name << std::endl;
    std::cout << "Entries (total):     " << nentries_total << std::endl;
    std::cout << "Entries (processed): " << nentries << std::endl;
    std::cout << "Ebeam:      " << Ebeam << " (GeV)" << std::endl;
    std::cout << "Cut value:  " << ycol_cut << std::endl;
    std::cout << "Skipped (bad Q2):    " << n_bad_Q2 << std::endl;
    std::cout << "Skipped (bad ycol):  " << n_bad_ycol << std::endl;
    std::cout << "Count P1_neg < 0.005:      " << n_P1neg_below << std::endl;
    std::cout << "Count (ycol - y) < 0.005:  " << n_ycolmy_below << std::endl;

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1110);

    TCanvas *c = new TCanvas("c", "dvcsgen ycol diagnostic", 1400, 600);
    c->Divide(2, 1);

    c->cd(1);
    gPad->SetLogy();
    h_P1neg->Draw("HIST");
    TLine *l1 = new TLine(ycol_cut, 0.0, ycol_cut, h_P1neg->GetMaximum());
    l1->SetLineWidth(2);
    l1->Draw("SAME");
    TLatex tx1;
    tx1.SetNDC();
    tx1.SetTextSize(0.040);
    tx1.DrawLatex(0.12, 0.86, "Vertical line at 0.005");
    tx1.DrawLatex(0.12, 0.80, "Candidate A: P1_neg = - (k - qgamma)^{2} / Q^{2}");

    c->cd(2);
    gPad->SetLogy();
    h_ycolmy->Draw("HIST");
    TLine *l2 = new TLine(ycol_cut, 0.0, ycol_cut, h_ycolmy->GetMaximum());
    l2->SetLineWidth(2);
    l2->Draw("SAME");
    TLatex tx2;
    tx2.SetNDC();
    tx2.SetTextSize(0.040);
    tx2.DrawLatex(0.12, 0.86, "Vertical line at 0.005");
    tx2.DrawLatex(0.12, 0.80, "Candidate B: ycol - y");

    const std::string out_png = out_dir + "/dvcsgen_ycol_diagnostic.png";
    const std::string out_pdf = out_dir + "/dvcsgen_ycol_diagnostic.pdf";

    c->SaveAs(out_png.c_str());
    c->SaveAs(out_pdf.c_str());

    std::cout << "Wrote: " << out_png << std::endl;
    std::cout << "Wrote: " << out_pdf << std::endl;

    delete c;
    delete h_P1neg;
    delete h_P1abs;
    delete h_ycolmy;
    f->Close();
    delete f;

    return 0;
}