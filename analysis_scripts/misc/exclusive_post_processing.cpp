#include <iostream>
#include <string>
#include <cmath>
#include <TFile.h>
#include <TTree.h>

// constants
static constexpr double m_e  = 0.000511;   // GeV
static constexpr double m_pi = 0.139570;   // GeV

// map run -> beam energy
static double beamEnergy(int runnum) {
    if (runnum >= 6616  && runnum <= 6783)   return 10.1998;
    if (runnum >= 16042 && runnum <= 17065)  return 10.5473;
    if (runnum >= 17067 && runnum <= 17724)  return 10.5563;
    if (runnum >= 17725 && runnum <= 17811)  return 10.5593;
    return 0.0;
}

// compute t = (q - p_pi)^2
static double compute_t_scalar(int runnum,
                               double e_p, double e_theta, double e_phi,
                               double p_p, double p_theta, double p_phi)
{
    double Eb = beamEnergy(runnum);
    if (Eb <= 0.0) return 1e6;

    // electron 4-vector
    double E_e   = std::sqrt(e_p*e_p + m_e*m_e);
    double sin_e = std::sin(e_theta), cos_e = std::cos(e_theta);
    double ex = e_p*sin_e*std::cos(e_phi);
    double ey = e_p*sin_e*std::sin(e_phi);
    double ez = e_p*cos_e;

    // pion 4-vector
    double E_pi  = std::sqrt(p_p*p_p + m_pi*m_pi);
    double sin_p = std::sin(p_theta), cos_p = std::cos(p_theta);
    double px = p_p*sin_p*std::cos(p_phi);
    double py = p_p*sin_p*std::sin(p_phi);
    double pz = p_p*cos_p;

    // q = (Eb - E_e, -ex, -ey, Eb - ez)
    double dE = (Eb - E_e) - E_pi;
    double dx = -ex - px;
    double dy = -ey - py;
    double dz = (Eb - ez) - pz;

    return dE*dE - (dx*dx + dy*dy + dz*dz);
}

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

    TFile *fin = TFile::Open(infile.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: could not open input file " << infile << "\n";
        return 1;
    }

    TTree *tin = static_cast<TTree*>(fin->Get("PhysicsEvents"));
    if (!tin) {
        std::cerr << "Error: tree PhysicsEvents not found in " << infile << "\n";
        return 1;
    }

    // we will recompute t, so disable input branch t if present
    tin->SetBranchStatus("t", 0);

    // set addresses
    Int_t    runnum;
    Double_t e_p, e_theta, e_phi;
    Double_t p_p, p_theta, p_phi;
    Double_t Mx2;  // needed for the cut Mx2 < 1.4

    tin->SetBranchAddress("runnum", &runnum);
    tin->SetBranchAddress("e_p", &e_p);
    tin->SetBranchAddress("e_theta", &e_theta);
    tin->SetBranchAddress("e_phi", &e_phi);
    tin->SetBranchAddress("p_p", &p_p);
    tin->SetBranchAddress("p_theta", &p_theta);
    tin->SetBranchAddress("p_phi", &p_phi);
    tin->SetBranchAddress("Mx2", &Mx2);

    TFile *fout = TFile::Open(outfile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: could not create output file " << outfile << "\n";
        return 1;
    }

    // clone structure (with 0 entries) and add our computed t
    TTree *tout = tin->CloneTree(0);
    Double_t t_val;
    tout->Branch("t", &t_val, "t/D");

    Long64_t n = tin->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tin->GetEntry(i);

        // apply the requested cut: only write events with Mx2 < 1.4
        if (Mx2 >= 1.4) continue;

        t_val = compute_t_scalar(runnum, e_p, e_theta, e_phi,
                                 p_p, p_theta, p_phi);
        tout->Fill();
    } // end for

    fout->Write();
    fout->Close();
    fin->Close();
    std::cout << "Wrote: " << outfile << "\n";
    return 0;
}