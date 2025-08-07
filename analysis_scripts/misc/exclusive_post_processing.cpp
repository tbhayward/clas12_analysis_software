#include <iostream>
#include <string>
#include <cmath>
#include <TFile.h>
#include <TTree.h>

// constants
static constexpr double m_e  = 0.000511;  // GeV
static constexpr double m_pi = 0.139570;  // GeV

//==============================================================================
// beamEnergy(run):
//============================================================================== 
static double beamEnergy(int runnum) {
    if (runnum >= 6616  && runnum <= 6783)   return 10.1998;
    if (runnum >= 16042 && runnum <= 17065)  return 10.5473;
    if (runnum >= 17067 && runnum <= 17724)  return 10.5563;
    if (runnum >= 17725 && runnum <= 17811)  return 10.5593;
    return 0.0;
}

//==============================================================================
// compute_t_scalar(...):
//============================================================================== 
static double compute_t_scalar(int runnum,
                               double e_p,     double e_theta, double e_phi,
                               double p_p,     double p_theta, double p_phi)
{
    double Eb = beamEnergy(runnum);
    if (Eb <= 0.0) return 1e6; // invalid run → flag

    // scattered electron 4-vector
    double E_e  = std::sqrt(e_p*e_p + m_e*m_e);
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

    // virtual photon q = p_beam - p_e'
    double E_q = Eb - E_e;
    double qx  = -ex;
    double qy  = -ey;
    double qz  = Eb - ez;

    // Δ = q - p_pi
    double dE = E_q - E_pi;
    double dx = qx  - px;
    double dy = qy  - py;
    double dz = qz  - pz;

    // t = Δ^2 = (dE)^2 - |Δp|^2
    return dE*dE - (dx*dx + dy*dy + dz*dz);
}

int main(int argc, char** argv) {
    if (argc!=2) {
        std::cerr<<"Usage: "<< argv[0] <<" <input.root>\n";
        return 1;
    }
    std::string infile = argv[1];
    // build output name by inserting "_2" before ".root"
    std::string outfile = infile;
    auto pos = outfile.rfind(".root");
    if (pos!=std::string::npos) outfile.insert(pos, "_2");
    else                    outfile += "_2.root";

    // open input
    TFile *fin = TFile::Open(infile.c_str(),"READ");
    if (!fin||fin->IsZombie()) {
        std::cerr<<"Error: cannot open "<< infile <<"\n";
        return 2;
    }
    TTree *tin = (TTree*)fin->Get("PhysicsEvents");
    if (!tin) {
        std::cerr<<"Error: no PhysicsEvents tree in "<< infile <<"\n";
        return 3;
    }

    // disable old t branch
    tin->SetBranchStatus("t",0);

    // set up branch addresses for inputs
    Int_t   runnum;
    Double_t e_p,   e_theta,   e_phi;
    Double_t p_p,   p_theta,   p_phi;
    tin->SetBranchAddress("runnum",  &runnum);
    tin->SetBranchAddress("e_p",      &e_p);
    tin->SetBranchAddress("e_theta",  &e_theta);
    tin->SetBranchAddress("e_phi",    &e_phi);
    tin->SetBranchAddress("p_p",      &p_p);
    tin->SetBranchAddress("p_theta",  &p_theta);
    tin->SetBranchAddress("p_phi",    &p_phi);

    // clone tree structure (no entries) into new file
    TFile *fout = TFile::Open(outfile.c_str(),"RECREATE");
    TTree *tout = tin->CloneTree(0);

    // create new t branch
    Double_t t_val;
    tout->Branch("t",&t_val,"t/D");

    // loop and fill
    Long64_t n = tin->GetEntries();
    for (Long64_t i=0; i<n; ++i) {
        tin->GetEntry(i);
        t_val = compute_t_scalar(runnum,
                                 e_p, e_theta, e_phi,
                                 p_p, p_theta, p_phi);
        tout->Fill();
    }

    // write and close
    fout->Write();
    fout->Close();
    fin ->Close();

    std::cout<<"Wrote updated file: "<< outfile <<"\n";
    return 0;
}