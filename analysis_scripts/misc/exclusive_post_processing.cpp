// postprocess_enpi_skim.cpp
// Tight skim for enpi+ PhysicsEvents:
//   Cuts: 0.65 < Mx2 < 1.125, fiducial_status == 111, 0.09 < x < 0.61
//   Drops from OUTPUT: fiducial_status, num_pos, num_neg, num_neutral,
//                     evnum, detector, xi, eta,
//                     any existing t/tmin/tprime/sinthetagamma (and mc_* versions)
//   Recomputes and ADDS: t, tmin (exact), tprime, sinthetagamma
//
// Compile (tcsh):
//   g++ -O2 -std=c++17 postprocess_enpi_skim.cpp `root-config --cflags --libs` -o postprocess_enpi_skim
//
// Run:
//   ./postprocess_enpi_skim /path/to/input.root
//
// Output filename is input with "_2" inserted before ".root" (or appended "_2.root" if no .root).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TSystem.h>

// -----------------------------------------------------------------------------
// constants
// -----------------------------------------------------------------------------
static constexpr double M_E  = 0.000511;        // (GeV)
static constexpr double M_PI = 0.139570;        // (GeV)
static constexpr double M_N  = 0.9382720813;    // (GeV) proton mass

// -----------------------------------------------------------------------------
// beam energy map (run -> Eb)
// -----------------------------------------------------------------------------
static double beamEnergy(int runnum) {
    if (runnum >= 6616  && runnum <= 6783)   return 10.1998;
    if (runnum >= 16042 && runnum <= 17065)  return 10.5473;
    if (runnum >= 17067 && runnum <= 17724)  return 10.5563;
    if (runnum >= 17725 && runnum <= 17811)  return 10.5593;
    return 10.5563;
}

// -----------------------------------------------------------------------------
// kinematics helpers
// -----------------------------------------------------------------------------
static double compute_t_scalar(int runnum,
                               double e_p, double e_theta, double e_phi,
                               double p_p, double p_theta, double p_phi) {
    double Eb = beamEnergy(runnum);
    if (Eb <= 0.0) return 1e9;

    // electron 3-momentum
    double E_e   = std::sqrt(e_p*e_p + M_E*M_E);
    double se    = std::sin(e_theta), ce = std::cos(e_theta);
    double ex    = e_p * se * std::cos(e_phi);
    double ey    = e_p * se * std::sin(e_phi);
    double ez    = e_p * ce;

    // pion 3-momentum
    double E_pi  = std::sqrt(p_p*p_p + M_PI*M_PI);
    double sp    = std::sin(p_theta), cp = std::cos(p_theta);
    double px    = p_p * sp * std::cos(p_phi);
    double py    = p_p * sp * std::sin(p_phi);
    double pz    = p_p * cp;

    // q = k - k' = (Eb - E_e, -ex, -ey, Eb - ez)
    double dE = (Eb - E_e) - E_pi;
    double dx = -ex - px;
    double dy = -ey - py;
    double dz = (Eb - ez) - pz;

    // t = (q - p_pi)^2 already folded into the dE,dx,dy,dz above
    return dE*dE - (dx*dx + dy*dy + dz*dz);
}

static double compute_sin_theta_gamma(double y, double xB, double Q2) {
    if (Q2 <= 0.0) return 0.0;
    const double Q     = std::sqrt(Q2);
    const double gamma = (Q > 0.0) ? (2.0 * xB * M_N / Q) : 0.0;
    const double num   = 1.0 - y - 0.25 * y * y * gamma * gamma;
    const double den   = 1.0 + gamma * gamma;
    if (den <= 0.0) return 0.0;
    const double ratio = std::max(0.0, num / den); // protect sqrt
    return gamma * std::sqrt(ratio);
}

// Exact tmin for DVCS/DVMP:
static double compute_tmin_exact(double xB, double Q2) {
    const bool xb_ok = (xB > 0.0 && xB < 1.0);
    if (Q2 <= 0.0 || !xb_ok) {
        if (xb_ok) {
            const double denom = (1.0 - xB);
            if (denom > 0.0) return - (M_N*xB)*(M_N*xB) / denom; // high-energy approx
        }
        return 0.0;
    }
    const double eps2 = 4.0 * M_N * M_N * xB * xB / Q2;
    const double root = std::sqrt(1.0 + eps2);
    const double num  = Q2 * ( 2.0*(1.0 - xB)*(1.0 - root) + eps2 );
    const double den  = 4.0*xB*(1.0 - xB) + eps2;
    if (den == 0.0) return 0.0;
    return - num / den; // negative
}

// -----------------------------------------------------------------------------
// tiny helpers
// -----------------------------------------------------------------------------
static bool has_branch(TTree* t, const char* name) {
    return t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name);
}

static bool has_all_branches(TTree* t, const std::vector<const char*>& names) {
    for (auto* n : names) {
        if (!has_branch(t, n)) return false;
    }
    return true;
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input.root>\n";
        return 1;
    }

    std::string infile(argv[1]);
    std::string outfile = infile;
    auto pos = outfile.rfind(".root");
    if (pos != std::string::npos) outfile.insert(pos, "_2");
    else                          outfile += "_2.root";

    // open input
    TFile* fin = TFile::Open(infile.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: could not open input file " << infile << "\n";
        return 1;
    }

    TTree* tin = static_cast<TTree*>(fin->Get("PhysicsEvents"));
    if (!tin) {
        std::cerr << "Error: tree PhysicsEvents not found in " << infile << "\n";
        fin->Close();
        return 1;
    }

    // -------------------------------------------------------------------------
    // Build output schema by cloning only the branches we want to KEEP.
    // Temporarily DISABLE branches we do not want in the output, clone an empty
    // structure, then re-enable for reading.
    // -------------------------------------------------------------------------
    std::vector<const char*> drop_from_output = {
        "fiducial_status","num_pos","num_neg","num_neutral",
        "evnum","detector","xi","eta",
        // Also drop any preexisting computed fields if present:
        "t","tmin","tprime","sinthetagamma",
        "mc_t","mc_tmin","mc_tprime","mc_sinthetagamma"
    };

    std::vector<const char*> to_disable;
    to_disable.reserve(drop_from_output.size());
    for (auto* nm : drop_from_output) {
        if (has_branch(tin, nm)) to_disable.push_back(nm);
    }

    tin->SetBranchStatus("*", 1);
    for (auto* nm : to_disable) tin->SetBranchStatus(nm, 0);

    // Create output file (use strong compression to save space)
    TFile* fout = TFile::Open(outfile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: could not create output file " << outfile << "\n";
        fin->Close();
        return 1;
    }

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
    fout->SetCompressionSettings(ROOT::RCompressionSetting::EDefaults::kUseZSTD);
#else
    fout->SetCompressionAlgorithm(ROOT::kZLIB);
    fout->SetCompressionLevel(9);
#endif

    // Clone structure (0 entries) with only enabled branches copied
    TTree* tout = tin->CloneTree(0, "fast");

    // Re-enable everything for reading from the input now
    tin->SetBranchStatus("*", 1);

    // -------------------------------------------------------------------------
    // Set up input addresses we need for cuts and computations
    // -------------------------------------------------------------------------
    Int_t    runnum = 0;
    Int_t    fiducial_status = 0;
    Double_t e_p=0, e_theta=0, e_phi=0;
    Double_t p_p=0, p_theta=0, p_phi=0;
    Double_t Mx2=0;
    Double_t xB=0, Q2=0, y=0;

    tin->SetBranchAddress("runnum",          &runnum);
    tin->SetBranchAddress("fiducial_status", &fiducial_status);
    tin->SetBranchAddress("e_p",             &e_p);
    tin->SetBranchAddress("e_theta",         &e_theta);
    tin->SetBranchAddress("e_phi",           &e_phi);
    tin->SetBranchAddress("p_p",             &p_p);
    tin->SetBranchAddress("p_theta",         &p_theta);
    tin->SetBranchAddress("p_phi",           &p_phi);
    tin->SetBranchAddress("Mx2",             &Mx2);
    tin->SetBranchAddress("x",               &xB);
    tin->SetBranchAddress("Q2",              &Q2);
    tin->SetBranchAddress("y",               &y);

    // Detect MC kinematics availability
    const bool mc_min = has_branch(tin, "mc_e_p");
    bool mc_present = false;
    Double_t mc_e_p=0, mc_e_theta=0, mc_e_phi=0;
    Double_t mc_p_p=0, mc_p_theta=0, mc_p_phi=0;
    Double_t mc_xB=0, mc_Q2=0, mc_y=0;

    if (mc_min) {
        mc_present = has_all_branches(tin, {
            "mc_e_p","mc_e_theta","mc_e_phi",
            "mc_p_p","mc_p_theta","mc_p_phi",
            "mc_x","mc_Q2","mc_y"
        });
        if (mc_present) {
            tin->SetBranchAddress("mc_e_p",     &mc_e_p);
            tin->SetBranchAddress("mc_e_theta", &mc_e_theta);
            tin->SetBranchAddress("mc_e_phi",   &mc_e_phi);
            tin->SetBranchAddress("mc_p_p",     &mc_p_p);
            tin->SetBranchAddress("mc_p_theta", &mc_p_theta);
            tin->SetBranchAddress("mc_p_phi",   &mc_p_phi);
            tin->SetBranchAddress("mc_x",       &mc_xB);
            tin->SetBranchAddress("mc_Q2",      &mc_Q2);
            tin->SetBranchAddress("mc_y",       &mc_y);
        } else {
            std::cerr << "[WARN] Found some mc_* but not all required; MC recalculation will be skipped.\n";
        }
    }

    // -------------------------------------------------------------------------
    // Add new computed branches to the OUTPUT tree
    // -------------------------------------------------------------------------
    Double_t t_val = 0.0;
    Double_t tmin_val = 0.0;
    Double_t tprime_val = 0.0;
    Double_t sinthetagamma_val = 0.0;

    tout->Branch("t",             &t_val,              "t/D");
    tout->Branch("tmin",          &tmin_val,           "tmin/D");
    tout->Branch("tprime",        &tprime_val,         "tprime/D");
    tout->Branch("sinthetagamma", &sinthetagamma_val,  "sinthetagamma/D");

    Double_t mc_t_val=0.0, mc_tmin_val=0.0, mc_tprime_val=0.0, mc_sinthetagamma_val=0.0;
    if (mc_present) {
        tout->Branch("mc_t",             &mc_t_val,             "mc_t/D");
        tout->Branch("mc_tmin",          &mc_tmin_val,          "mc_tmin/D");
        tout->Branch("mc_tprime",        &mc_tprime_val,        "mc_tprime/D");
        tout->Branch("mc_sinthetagamma", &mc_sinthetagamma_val, "mc_sinthetagamma/D");
    }

    // Optional I/O tuning for size/performance
    tout->SetAutoSave(0);
    tout->SetAutoFlush(100000);

    // -------------------------------------------------------------------------
    // Event loop with requested cuts
    //   Cuts:
    //     0.65 < Mx2 < 1.125
    //     fiducial_status == 111
    //     0.09 < x < 0.61
    // -------------------------------------------------------------------------
    const double x_lo = 0.090;
    const double x_hi = 0.610;

    const Long64_t n = tin->GetEntries();
    Long64_t n_kept = 0;
    for (Long64_t i = 0; i < n; ++i) {
        tin->GetEntry(i);

        if (!(Mx2 > 0.65 && Mx2 < 1.125)) continue;
        if (fiducial_status != 111) continue;
        if (!(xB > x_lo && xB < x_hi)) continue;

        // Compute data values
        t_val             = compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi);
        tmin_val          = compute_tmin_exact(xB, Q2);
        tprime_val        = t_val - tmin_val;
        sinthetagamma_val = compute_sin_theta_gamma(y, xB, Q2);

        // MC values if available
        if (mc_present) {
            mc_t_val             = compute_t_scalar(runnum, mc_e_p, mc_e_theta, mc_e_phi,
                                                    mc_p_p, mc_p_theta, mc_p_phi);
            mc_tmin_val          = compute_tmin_exact(mc_xB, mc_Q2);
            mc_tprime_val        = mc_t_val - mc_tmin_val;
            mc_sinthetagamma_val = compute_sin_theta_gamma(mc_y, mc_xB, mc_Q2);
        }

        tout->Fill();
        ++n_kept;
    }

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Input:  " << infile  << "\n";
    std::cout << "Output: " << outfile << "\n";
    std::cout << "Events kept: " << n_kept << " / " << n << "\n";
    std::cout << "Applied cuts: 0.65 < Mx2 < 1.125, fiducial_status == 111, 0.09 < x < 0.61\n";
    std::cout << "Dropped from OUTPUT: fiducial_status, num_pos, num_neg, num_neutral, evnum, detector, xi, eta\n";
    std::cout << "Added branches: t, tmin(exact), tprime, sinthetagamma"
              << (mc_present ? ", mc_t, mc_tmin(exact), mc_tprime, mc_sinthetagamma\n" : "\n");
    std::cout << "Reminder: for plots usually use -tprime.\n";
    return 0;
}