// Include necessary ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>

using namespace std;

// function to get tmin 
double gettmin(double x) {
    double mp = 0.938272; // proton mass in GeV
    return -pow((mp*x),2)/(1-x);
}

// function get t
double gett(double p, double theta) {
    double mp = 0.938272; // proton mass in GeV
    double E = mp; // target proton energy (written as E to help checking calculation but at rest)
    return 2*mp*(E - p) - 2*sqrt(mp*mp + E*E)*sqrt(mp*mp + p*p) +
          2*sqrt(mp*mp + E*E)*sqrt(mp*mp + p*p)*cos(theta);
}

// function to get the polarization value
double getPol(int runnum) {
  double pol = 0.86; 
    if (runnum == 11 ) { pol = 0.86; } // runnum == 11 indicates Monte Carlo in CLAS12
    else if (runnum >= 5032 && runnum < 5333) { pol = 0.8592; } 
    else if (runnum >= 5333 && runnum <= 5666) { pol = 0.8922; }
    else if (runnum >= 6616 && runnum <= 6783) { pol = 0.8453; }
    else if (runnum >= 6142 && runnum <= 6149) { pol = 0.81132; }
    else if (runnum >= 6150 && runnum <= 6188) { pol = 0.82137; }
    else if (runnum >= 6189 && runnum <= 6260) { pol = 0.83598; }
    else if (runnum >= 6261 && runnum <= 6339) { pol = 0.80770; }
    else if (runnum >= 6340 && runnum <= 6342) { pol = 0.85536; }
    else if (runnum >= 6344 && runnum <= 6399) { pol = 0.87038; }
    else if (runnum >= 6420 && runnum <= 6476) { pol = 0.88214; }
    else if (runnum >= 6479 && runnum <= 6532) { pol = 0.86580; }
    else if (runnum >= 6533 && runnum <= 6603) { pol = 0.87887; }
    else if (runnum >= 11013 && runnum <= 11309) { pol = 0.84983; }
    else if (runnum >= 11323 && runnum <= 11334) { pol = 0.87135; }
    else if (runnum >= 11335 && runnum <= 11387) { pol = 0.85048; }
    else if (runnum >= 11389 && runnum <= 11571) { pol = 0.84262; }
    else if (runnum >= 16000) { pol = 0.83534; } // RGC +/- 0.01440
  return pol;
}

int main(int argc, char *argv[]) {
    // Check for correct number of command line arguments
    if (argc != 4) {
        cout << "Usage: " << argv[0];
        cout << " <input_text_file> <output_root_file> <hadron_count>" << endl;
        return 1;
    }
    
    // Open input text file for reading
    ifstream infile(argv[1]);
    if (!infile.is_open()) {
        cout << "Error opening input text file." << endl;
        return 1;
    }

    // Create a new ROOT file and TTree for output
    TFile *outfile = new TFile(argv[2], "RECREATE");
    TTree *tree = new TTree("PhysicsEvents", "Physics Events Tree");

    // Determine the hadron count from the command line argument
    int hadron_count = atoi(argv[3]);

    // Declare common variables
    int runnum, evnum, helicity;
    double beam_pol, e_p, e_theta, e_phi, vz_e, Q2, W, Mx, Mx2, x, y, t, tmin;
    double z, xF, pT, zeta, eta, phi, DepA, DepB, DepC, DepV, DepW;
    double p_p, p_theta, p_phi, vz_p;
    // Additional variables for two hadrons
    double p1_p, p1_theta, p1_phi, vz_p1, p2_p, p2_theta, p2_phi, vz_p2;
    double z1, z2, Mh, xF1, xF2, pT1, pT2, pTpT, zeta1, zeta2;
    double eta1, eta2, Delta_eta, eta1_gN, eta2_gN;
    double phi1, phi2, Delta_phi, phih, phiR, theta;
    double DepA, DepB, DepC, DepV, DepW;

    // Case for zero hadrons (inclusive)
    if (hadron_count == 0) {
        // Link TTree branches to variables for zero hadrons
        tree->Branch("runnum", &runnum, "runnum/I");
        tree->Branch("evnum", &evnum, "evnum/I");
        tree->Branch("helicity", &helicity, "helicity/I");
        tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("tmin", &tmin, "tmin/D");
    }

    // Case for one hadron
    else if (hadron_count == 1) {
        // Link TTree branches to variables for one hadron
        tree->Branch("runnum", &runnum, "runnum/I");
        tree->Branch("evnum", &evnum, "evnum/I");
        tree->Branch("helicity", &helicity, "helicity/I");
        tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("p_p", &p_p, "p_p/D");
        tree->Branch("p_theta", &p_theta, "p_theta/D");
        tree->Branch("p_phi", &p_phi, "p_phi/D");
        tree->Branch("vz_p", &vz_p, "vz_p/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("z", &z, "z/D");
        tree->Branch("xF", &xF, "xF/D");
        tree->Branch("pT", &pT, "pT/D");
        tree->Branch("zeta", &zeta, "zeta/D");
        tree->Branch("eta", &eta, "eta/D");
        tree->Branch("phi", &phi, "phi/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
    }

    // Case for two hadrons (dihadrons)
    else if (hadron_count == 2) {

        // Link TTree branches to variables for two hadrons
        tree->Branch("runnum", &runnum, "runnum/I");
        tree->Branch("evnum", &evnum, "evnum/I");
        tree->Branch("helicity", &helicity, "helicity/I");
        tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
        tree->Branch("e_p", &e_p, "e_p/D");
        tree->Branch("e_theta", &e_theta, "e_theta/D");
        tree->Branch("e_phi", &e_phi, "e_phi/D");
        tree->Branch("vz_e", &vz_e, "vz_e/D");
        tree->Branch("p1_p", &p1_p, "p1_p/D");
        tree->Branch("p1_theta", &p1_theta, "p1_theta/D");
        tree->Branch("p1_phi", &p1_phi, "p1_phi/D");
        tree->Branch("vz_p1", &vz_p1, "vz_p1/D");
        tree->Branch("p2_p", &p2_p, "p2_p/D");
        tree->Branch("p2_theta", &p2_theta, "p2_theta/D");
        tree->Branch("p2_phi", &p2_phi, "p2_phi/D");
        tree->Branch("vz_p2", &vz_p2, "vz_p2/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W", &W, "W/D");
        tree->Branch("Mx", &Mx, "Mx/D");
        tree->Branch("Mx2", &Mx2, "Mx2/D");
        tree->Branch("x", &x, "x/D");
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("tmin", &tmin, "tmin/D");
        tree->Branch("z1", &z1, "z1/D");
        tree->Branch("z2", &z2, "z2/D");
        tree->Branch("Mh", &Mh, "Mh/D");
        tree->Branch("xF1", &xF1, "xF1/D");
        tree->Branch("xF2", &xF2, "xF2/D");
        tree->Branch("pT1", &pT1, "pT1/D");
        tree->Branch("pT2", &pT2, "pT2/D");
        tree->Branch("pTpT", &pTpT, "pTpT/D");
        tree->Branch("zeta1", &zeta1, "zeta1/D");
        tree->Branch("zeta2", &zeta2, "zeta2/D");
        tree->Branch("eta1", &eta1, "eta1/D");
        tree->Branch("eta2", &eta2, "eta2/D");
        tree->Branch("Delta_eta", &Delta_eta, "Delta_eta/D");
        tree->Branch("eta1_gN", &eta1_gN, "eta1_gN/D");
        tree->Branch("eta2_gN", &eta2_gN, "eta2_gN/D");
        tree->Branch("phi1", &phi1, "phi1/D");
        tree->Branch("phi2", &phi2, "phi2/D");
        tree->Branch("Delta_phi", &Delta_phi, "Delta_phi/D");
        tree->Branch("phih", &phih, "phih/D");
        tree->Branch("phiR", &phiR, "phiR/D");
        tree->Branch("theta", &theta, "theta/D");
        tree->Branch("DepA", &DepA, "DepA/D");
        tree->Branch("DepB", &DepB, "DepB/D");
        tree->Branch("DepC", &DepC, "DepC/D");
        tree->Branch("DepV", &DepV, "DepV/D");
        tree->Branch("DepW", &DepW, "DepW/D");
    }

    // Loop to read each line from the text file and fill the TTree based on hadron_count
    if (hadron_count == 0) {
        while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
            Q2 >> W >> Mx >> Mx2 >> x >> y) {
            beam_pol = getPol(runnum);
            t = gett(e_p, e_theta); // for inclusive we calculate t with electron kinematics
            tmin = gettmin(x);  
            tree->Fill(); // Fill the tree with the read data
        }
    } 
    else if (hadron_count == 1) {
        while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
            p_p >> p_theta >> p_phi >> vz_p >> Q2 >> W >> Mx >> Mx2 >> x >> y >> z >> xF >> 
            pT >> zeta >> eta >> phi >> DepA >> DepB >> DepC >> DepV >> DepW) {

            beam_pol = getPol(runnum);
            t = gett(p_p, p_theta); // for SIDIS we calculate t with proton kinematics
            tmin = gettmin(x); 
            tree->Fill(); // Fill the tree with the read data
        }
    } 
    else if (hadron_count == 2) {
        while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
            p1_p >> p1_theta >> p1_phi >> vz_p1 >> p2_p >> p2_theta >> p2_phi >> vz_p2 >> 
            Q2 >> W >> Mx >> Mx2 >> x >> y >> z1 >> z2 >> Mh >> xF1 >> xF2 >> pT1 >> pT2 >> 
            pTpT >> zeta1 >> zeta2 >> eta1 >> eta2 >> Delta_eta >> eta1_gN >> eta2_gN >> 
            phi1 >> phi2 >> Delta_phi >> phih >> phiR >> theta >> 
            DepA >> DepB >> DepC >> DepV >> DepW) {
            beam_pol = getPol(runnum);
            t = gett(p2_p, p2_theta); // for SIDIS we calculate t with proton kinematics
            tmin = gettmin(x); 
            tree->Fill(); // Fill the tree with the read data
        }
    }

    // Write the TTree to the ROOT file and close it
    tree->Write();
    outfile->Close();

    // Close the input text file
    infile.close();

    return 0;
}
