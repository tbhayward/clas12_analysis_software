#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>

// -----------------------------------------------------------------------------
// constants
// -----------------------------------------------------------------------------
static constexpr double m_e  = 0.000511;        // GeV
static constexpr double m_pi = 0.139570;        // GeV
static constexpr double M_N  = 0.9382720813;    // GeV (nucleon mass, proton)

// map run -> beam energy
static double beamEnergy(int runnum) {
    if (runnum >= 6616  && runnum <= 6783)   return 10.1998;
    if (runnum >= 16042 && runnum <= 17065)  return 10.5473;
    if (runnum >= 17067 && runnum <= 17724)  return 10.5563;
    if (runnum >= 17725 && runnum <= 17811)  return 10.5593;
    return 10.5563;
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

// compute sin(theta_gamma):
//   sin θγ = γ * sqrt( (1 - y - (1/4) y^2 γ^2) / (1 + γ^2) ),  γ = 2 xB M_N / Q
static double compute_sin_theta_gamma(double y, double xB, double Q2)
{
    if (Q2 <= 0.0) return 0.0;
    const double Q      = std::sqrt(Q2);
    const double gamma  = (Q > 0.0) ? (2.0 * xB * M_N / Q) : 0.0;
    const double num    = 1.0 - y - 0.25 * y * y * gamma * gamma;
    const double den    = 1.0 + gamma * gamma;
    if (den <= 0.0) return 0.0;
    const double ratio  = std::max(0.0, num / den);   // protect sqrt
    return gamma * std::sqrt(ratio);
}

// Q^2-dependent tmin (exact DVCS/DVMP form):
//   eps^2 = 4 M^2 xB^2 / Q^2
//   tmin  = - Q^2 [ 2(1-xB)(1 - sqrt(1 + eps^2)) + eps^2 ] / [ 4 xB(1-xB) + eps^2 ]
// Fallback to high-energy approx if Q2<=0 or xB out of (0,1):
//   tmin ≈ - (M xB)^2 / (1 - xB)
static double compute_tmin_exact(double xB, double Q2)
{
    const bool xb_ok = (xB > 0.0 && xB < 1.0);
    if (Q2 <= 0.0 || !xb_ok) {
        if (xb_ok) {
            const double denom = (1.0 - xB);
            if (denom > 0.0) return -(M_N*xB)*(M_N*xB) / denom;
        }
        return 0.0;
    }

    const double eps2  = 4.0 * M_N * M_N * xB * xB / Q2;
    const double root  = std::sqrt(1.0 + eps2);
    const double num   = Q2 * ( 2.0*(1.0 - xB)*(1.0 - root) + eps2 );
    const double den   = 4.0*xB*(1.0 - xB) + eps2;
    if (den == 0.0) return 0.0;
    return - num / den;  // tmin is negative
}

static bool has_branch(TTree* t, const char* name) {
    return t && t->GetListOfBranches() && t->GetListOfBranches()->FindObject(name);
}

// Convenience: check a list of branch names are all present
static bool has_all_branches(TTree* t, const std::vector<const char*>& names) {
    for (auto* n : names) {
        if (!has_branch(t, n)) return false;
    }
    return true;
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

    // We will recompute t, tmin, tprime; do not read (or clone) any existing versions.
    if (has_branch(tin, "t"))             tin->SetBranchStatus("t", 0);
    if (has_branch(tin, "tmin"))          tin->SetBranchStatus("tmin", 0);
    if (has_branch(tin, "sinthetagamma")) tin->SetBranchStatus("sinthetagamma", 0);

    // Detect if MC kinematics exist; if so, also recompute MC versions
    const bool mc_present_min = has_branch(tin, "mc_e_p");
    bool mc_present = false;
    if (mc_present_min) {
        if (has_branch(tin, "mc_t"))             tin->SetBranchStatus("mc_t", 0);
        if (has_branch(tin, "mc_tmin"))          tin->SetBranchStatus("mc_tmin", 0);
        if (has_branch(tin, "mc_sinthetagamma")) tin->SetBranchStatus("mc_sinthetagamma", 0);

        mc_present = has_all_branches(tin, {
            "mc_e_p","mc_e_theta","mc_e_phi",
            "mc_p_p","mc_p_theta","mc_p_phi",
            "mc_x","mc_Q2","mc_y"
        });
        if (!mc_present) {
            std::cerr << "[WARN] Found 'mc_e_p' but not all required mc_* branches; "
                         "MC recalculation will be skipped.\n";
        }
    }

    // set addresses (need x, Q2, y for sin(theta_gamma) and tmin)
    Int_t    runnum;
    Double_t e_p, e_theta, e_phi;
    Double_t p_p, p_theta, p_phi;
    Double_t Mx2;
    Double_t xB, Q2, y;

    tin->SetBranchAddress("runnum",   &runnum);
    tin->SetBranchAddress("e_p",      &e_p);
    tin->SetBranchAddress("e_theta",  &e_theta);
    tin->SetBranchAddress("e_phi",    &e_phi);
    tin->SetBranchAddress("p_p",      &p_p);
    tin->SetBranchAddress("p_theta",  &p_theta);
    tin->SetBranchAddress("p_phi",    &p_phi);
    tin->SetBranchAddress("Mx2",      &Mx2);
    tin->SetBranchAddress("x",        &xB);
    tin->SetBranchAddress("Q2",       &Q2);
    tin->SetBranchAddress("y",        &y);

    // MC addresses (if available)
    Double_t mc_e_p=0, mc_e_theta=0, mc_e_phi=0;
    Double_t mc_p_p=0, mc_p_theta=0, mc_p_phi=0;
    Double_t mc_xB=0, mc_Q2=0, mc_y=0;

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
    }

    TFile *fout = TFile::Open(outfile.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: could not create output file " << outfile << "\n";
        return 1;
    }

    // clone structure (with 0 entries) and add our computed variables
    TTree *tout = tin->CloneTree(0);

    // Data outputs
    Double_t t_val          = 0.0;
    Double_t tmin_val       = 0.0;   // recalculated exact tmin(x,Q2)
    Double_t tprime_val     = 0.0;   // t - tmin (≤ 0 typically)
    Double_t sinthetagamma  = 0.0;

    tout->Branch("t",             &t_val,         "t/D");
    tout->Branch("tmin",          &tmin_val,      "tmin/D");              // overwrite with exact form
    tout->Branch("tprime",        &tprime_val,    "tprime/D");
    tout->Branch("sinthetagamma", &sinthetagamma, "sinthetagamma/D");

    // MC outputs (if inputs exist)
    Double_t mc_t_val = 0.0;
    Double_t mc_tmin_val = 0.0;
    Double_t mc_tprime_val = 0.0;
    Double_t mc_sinthetagamma = 0.0;

    if (mc_present) {
        tout->Branch("mc_t",             &mc_t_val,         "mc_t/D");
        tout->Branch("mc_tmin",          &mc_tmin_val,      "mc_tmin/D");
        tout->Branch("mc_tprime",        &mc_tprime_val,    "mc_tprime/D");
        tout->Branch("mc_sinthetagamma", &mc_sinthetagamma, "mc_sinthetagamma/D");
    }

    const Long64_t n = tin->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tin->GetEntry(i);

        // write only events with Mx2 < 1.3 (original post-process filter)
        if (Mx2 >= 1.3) continue;

        // DATA: compute new values
        t_val          = compute_t_scalar(runnum, e_p, e_theta, e_phi, p_p, p_theta, p_phi);
        tmin_val       = compute_tmin_exact(xB, Q2);
        tprime_val     = t_val - tmin_val;                 // t' = t - tmin  (usually plot -t')
        sinthetagamma  = compute_sin_theta_gamma(y, xB, Q2);

        // MC: compute if available
        if (mc_present) {
            mc_t_val         = compute_t_scalar(runnum, mc_e_p, mc_e_theta, mc_e_phi,
                                                mc_p_p, mc_p_theta, mc_p_phi);
            mc_tmin_val      = compute_tmin_exact(mc_xB, mc_Q2);
            mc_tprime_val    = mc_t_val - mc_tmin_val;
            mc_sinthetagamma = compute_sin_theta_gamma(mc_y, mc_xB, mc_Q2);
        }

        tout->Fill();
    }

    fout->Write();
    fout->Close();
    fin->Close();
    std::cout << "Wrote: " << outfile << "\n";
    std::cout << "Added/overwrote branches: t, tmin(exact), tprime, sinthetagamma"
              << (mc_present ? ", mc_t, mc_tmin(exact), mc_tprime, mc_sinthetagamma\n" : "\n");
    std::cout << "Reminder: for plots use  -tprime  (i.e. -(t - tmin)).\n";
    return 0;
}