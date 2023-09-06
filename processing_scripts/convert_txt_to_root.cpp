// Include necessary ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]) {
    // Check for correct number of command line arguments
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <input_text_file> <output_root_file>" << endl;
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

    // Declare variables to store data for each event
    int runnum, evnum, helicity;
    double e_p, e_theta, e_phi, vz_e, p_p, p_theta, p_phi, vz_p, Q2, W, Mx, Mx2;
    double x, y, z, xF, pT, zeta, eta, phi, DepA, DepB, DepC, DepV, DepW;

    // Link TTree branches to variables
    tree->Branch("runnum", &runnum, "runnum/I");
    tree->Branch("evnum", &evnum, "evnum/I");
    tree->Branch("helicity", &helicity, "helicity/I");
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

    // Loop to read each line from the text file and fill the TTree
    while (infile >> runnum >> evnum >> helicity >> e_p >> e_theta >> e_phi >> vz_e >> 
        p_p >> p_theta >> p_phi >> vz_p >> Q2 >> W >> Mx >> Mx2 >> x >> y >> z >> xF >> pT >> 
        zeta >> eta >> phi >> DepA >> DepB >> DepC >> DepV >> DepW) {
        tree->Fill(); // Fill the tree with the read data
    }

    // Write the TTree to the ROOT file and close it
    tree->Write();
    outfile->Close();

    // Close the input text file
    infile.close();

    return 0;
}
