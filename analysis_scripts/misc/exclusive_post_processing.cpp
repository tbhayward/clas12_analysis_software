// exclusive_post_processing.cpp
//
// Skimmer for CLAS12 enpi+ PhysicsEvents.
//
// Selection (applied via CopyTree):
//   0.65 < Mx2 < 1.125
//   fiducial_status == 111
//   0.09 < x < 0.61
//
// Output policy:
//   - Keep all branches EXCEPT the following eight:
//       fiducial_status, num_pos, num_neg, num_neutral, evnum, detector, xi, eta
//   - Preserve existing t, tmin, tprime, sinthetagamma if they are already in the input.
//   - If any of those four are missing, compute and add ONLY the missing ones.
//
// Sanity printouts:
//   - Input entries
//   - Expected entries from the selection expression
//   - Actual entries written to output
//   - Mx2 min/mean in the OUTPUT tree (should satisfy the lower cut)
//
// Build (tcsh/bash):
//   g++ -O2 -std=c++17 exclusive_post_processing.cpp `root-config --cflags --libs` -o exclusive_post_processing
//
// Run:
//   ./exclusive_post_processing /path/to/input.root
//
// Output file name:
//   Inserts "_2" before ".root" (or appends "_2.root" if no .root suffix).

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>
#include <RVersion.h>

// ----------------------------- constants -------------------------------------
// (Use names that won't collide with <cmath> macros like M_E, M_PI)
static constexpr double MASS_E  = 0.000511;        // GeV (electron)
static constexpr double MASS_PI = 0.139570;        // GeV (charged pion)
static constexpr double MASS_N  = 0.9382720813;    // GeV (proton)

// ------------------------------ helpers --------------------------------------
static double beamEnergy(int runnum) {
    if (runnum >= 6616  && runnum <= 6783)   return 10.1998;
    if (runnum >= 16042 && runnum <= 17065)  return 10.5473;
    if (runnum >= 17067 && runnum <= 17724)  return 10.5563;
    if (runnum >= 17725 && runnum <= 17811)  return 10.5593;
    return 10.5563;
}

static bool has_branch(TTree* t, const char* name) {
    return t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name);
}

static double compute_t_scalar(int runnum,
                               double e_p, double e_theta, double e_phi,
                               double p_p, double p_theta, double p_phi) {
    const double Eb = beamEnergy(runnum);
    if (Eb <= 0.0) return 1e9;

    const double E_e = std::sqrt(e_p*e_p + MASS_E*MASS_E);
    const double se  = std::sin(e_theta), ce = std::cos(e_theta);
    const double ex  = e_p * se * std::cos(e_phi);
    const double ey  = e_p * se * std::sin(e_phi);
    const double ez  = e_p * ce;

    const double E_pi = std::sqrt(p_p*p_p + MASS_PI*MASS_PI);
    const double sp   = std::sin(p_theta), cp = std::cos(p_theta);
    const double px   = p_p * sp * std::cos(p_phi);
    const double py   = p_p * sp * std::sin(p_phi);
    const double pz   = p_p * cp;

    // q = k - k' = (Eb - E_e, -ex, -ey, Eb - ez)
    const double dE = (Eb - E_e) - E_pi;
    const double dx = -ex - px;
    const double dy = -ey - py;
    const double dz = (Eb - ez) - pz;

    return dE*dE - (dx*dx + dy*dy + dz*dz);
}

static double compute_sin_theta_gamma(double y, double xB, double Q2) {
    if (Q2 <= 0.0) return 0.0;
    const double Q     = std::sqrt(Q2);
    const double gamma = (Q > 0.0) ? (2.0 * xB * MASS_N / Q) : 0.0;
    const double num   = 1.0 - y - 0.25 * y * y * gamma * gamma;
    const double den   = 1.0 + gamma * gamma;
    if (den <= 0.0) return 0.0;
    const double ratio = std::max(0.0, num / den);
    return gamma * std::sqrt(ratio);
}

static double compute_tmin_exact(double xB, double Q2) {
    const bool xb_ok = (xB > 0.0 && xB < 1.0);
    if (Q2 <= 0.0 || !xb_ok) {
        if (xb_ok) {
            const double denom = (1.0 - xB);
            if (denom > 0.0) return - (MASS_N*xB)*(MASS_N*xB) / denom; // high-E approx
        }
        return 0.0;
    }
    const double eps2 = 4.0*MASS_N*MASS_N*xB*xB / Q2;
    const double root = std::sqrt(1.0 + eps2);
    const double num  = Q2 * ( 2.0*(1.0 - xB)*(1.0 - root) + eps2 );
    const double den  = 4.0*xB*(1.0 - xB) + eps2;
    if (den == 0.0) return 0.0;
    return - num / den;
}

// -------------------------------- main ---------------------------------------
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

    // Open input
    TFile* fin = TFile::Open(infile.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: could not open input " << infile << "\n";
        return 1;
    }
    TTree* tin = static_cast<TTree*>(fin->Get("PhysicsEvents"));
    if (!tin) {
        std::cerr << "Error: PhysicsEvents not found in " << infile << "\n";
        fin->Close();
        return 1;
    }

    // Selection (used by CopyTree)
    const char* CUTS = "Mx2>0.65 && Mx2<1.125 && fiducial_status==111 && x>0.09 && x<0.61";

    // Counts for sanity
    const Long64_t n_in     = tin->GetEntries();
    const Long64_t n_expect = tin->GetEntries(CUTS);

    // Perform the skim (this enforces the cuts)
    TTree* tskim = tin->CopyTree(CUTS);
    if (!tskim) {
        std::cerr << "CopyTree returned null.\n";
        fin->Close();
        return 1;
    }

    // Determine which computed branches already exist in the skim
    const bool have_t       = has_branch(tskim, "t");
    const bool have_tmin    = has_branch(tskim, "tmin");
    const bool have_tprime  = has_branch(tskim, "tprime");
    const bool have_stg     = has_branch(tskim, "sinthetagamma");

    // Disable the 8 columns to drop in the OUTPUT, then clone kept structure
    const char* DROP[] = {"fiducial_status","num_pos","num_neg","num_neutral",
                          "evnum","detector","xi","eta"};
    tskim->SetBranchStatus("*", 1);
    for (const char* nm : DROP) if (has_branch(tskim, nm)) tskim->SetBranchStatus(nm, 0);

    // Create output file with portable compression (zlib level 9 works everywhere)
    TFile* fout = TFile::Open(outfile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: could not create output " << outfile << "\n";
        fin->Close();
        return 1;
    }
    fout->SetCompressionAlgorithm(ROOT::kZLIB);
    fout->SetCompressionLevel(9);

    // Clone structure of active (kept) branches only; 0 entries to start
    TTree* tout = tskim->CloneTree(0, "fast");

    // If any of (t, tmin, tprime, sinthetagamma) are missing, add them to OUTPUT
    Double_t t_val=0, tmin_val=0, tprime_val=0, sinthetagamma_val=0;
    TBranch *b_t=nullptr, *b_tmin=nullptr, *b_tprime=nullptr, *b_stg=nullptr;

    if (!have_t)      b_t      = tout->Branch("t",             &t_val,              "t/D");
    if (!have_tmin)   b_tmin   = tout->Branch("tmin",          &tmin_val,           "tmin/D");
    if (!have_tprime) b_b_tprime = nullptr; // placeholder to keep naming consistent
    if (!have_tprime) b_tprime = tout->Branch("tprime",        &tprime_val,         "tprime/D");
    if (!have_stg)    b_stg    = tout->Branch("sinthetagamma", &sinthetagamma_val,  "sinthetagamma/D");

    // Set up addresses for inputs needed to compute missing outputs
    Int_t    runnum=0;
    Double_t e_p=0,e_theta=0,e_phi=0;
    Double_t p_p=0,p_theta=0,p_phi=0;
    Double_t xB=0,Q2=0,yv=0;

    bool need_compute = (!have_t) || (!have_tmin) || (!have_tprime) || (!have_stg);
    if (need_compute) {
        if (!tskim->GetBranch("runnum") ||
            !tskim->GetBranch("e_p")    || !tskim->GetBranch("e_theta") || !tskim->GetBranch("e_phi") ||
            !tskim->GetBranch("p_p")    || !tskim->GetBranch("p_theta") || !tskim->GetBranch("p_phi") ||
            !tskim->GetBranch("x")      || !tskim->GetBranch("Q2")      || !tskim->GetBranch("y")) {
            std::cerr << "Missing inputs for computing derived branches in skim.\n";
            fout->Close(); fin->Close();
            return 1;
        }
        tskim->SetBranchStatus("runnum", 1);
        tskim->SetBranchStatus("e_p",    1);
        tskim->SetBranchStatus("e_theta",1);
        tskim->SetBranchStatus("e_phi",  1);
        tskim->SetBranchStatus("p_p",    1);
        tskim->SetBranchStatus("p_theta",1);
        tskim->SetBranchStatus("p_phi",  1);
        tskim->SetBranchStatus("x",      1);
        tskim->SetBranchStatus("Q2",     1);
        tskim->SetBranchStatus("y",      1);

        tskim->SetBranchAddress("runnum",  &runnum);
        tskim->SetBranchAddress("e_p",     &e_p);
        tskim->SetBranchAddress("e_theta", &e_theta);
        tskim->SetBranchAddress("e_phi",   &e_phi);
        tskim->SetBranchAddress("p_p",     &p_p);
        tskim->SetBranchAddress("p_theta", &p_theta);
        tskim->SetBranchAddress("p_phi",   &p_phi);
        tskim->SetBranchAddress("x",       &xB);
        tskim->SetBranchAddress("Q2",      &Q2);
        tskim->SetBranchAddress("y",       &yv);
    }

    // I/O tuning (optional)
    tout->SetAutoSave(0);
    tout->SetAutoFlush(100000);

    // Loop over skimmed entries: fill kept branches; compute and fill only missing outputs
    const Long64_t nsel = tskim->GetEntries();
    double mx2_min = 1e300, mx2_sum = 0.0;
    Double_t mx2_tmp = 0.0;
    tskim->SetBranchStatus("Mx2", 1);
    tskim->SetBranchAddress("Mx2", &mx2_tmp);

    for (Long64_t i=0; i<nsel; ++i) {
        tskim->GetEntry(i);

        if (need_compute) {
            if (!have_t || !have_tmin || !have_tprime) {
                const double tt   = compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi);
                const double tmin = compute_tmin_exact(xB, Q2);
                if (!have_t)      t_val      = tt;
                if (!have_tmin)   tmin_val   = tmin;
                if (!have_tprime) tprime_val = tt - tmin;
            }
            if (!have_stg) {
                sinthetagamma_val = compute_sin_theta_gamma(yv, xB, Q2);
            }
        }

        if (mx2_tmp < mx2_min) mx2_min = mx2_tmp;
        mx2_sum += mx2_tmp;

        if (b_t)      b_t->Fill();
        if (b_tmin)   b_tmin->Fill();
        if (b_tprime) b_tprime->Fill();
        if (b_stg)    b_stg->Fill();

        tout->Fill();
    }

    // Write output
    fout->cd();
    tout->Write("PhysicsEvents");  // keep same name
    fout->Close();
    fin->Close();

    // Final sanity summary
    const double mx2_mean = (nsel>0 ? mx2_sum / nsel : 0.0);
    std::cout << "Input entries : " << n_in << "\n";
    std::cout << "Selected (exp): " << n_expect << " (by CopyTree expression)\n";
    std::cout << "Selected (out): " << nsel << " (actual written)\n";
    std::cout << "Mx2 in output : min=" << mx2_min << "  mean=" << mx2_mean << "\n";
    std::cout << "Dropped cols  : fiducial_status, num_pos, num_neg, num_neutral, evnum, detector, xi, eta\n";
    std::cout << "Preserved     : existing t/tmin/tprime/sinthetagamma; added only if missing.\n";

    return 0;
}