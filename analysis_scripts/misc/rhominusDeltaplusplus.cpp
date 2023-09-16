// Created 9/6/23
// Studies of an epi+pi-pX final state to look for rho- Delta++ events

#include <TCanvas.h>
#include <TH1F.h>
#include <TObjArray.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <iostream>
#include <string>
#include <map>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

// function get t
double gett(double p, double theta) {
    double mp = 0.938272; // proton mass in GeV
    double E = mp; // target proton energy (written as E to help checking calculation but at rest)
    return 2*mp*(E - p) - 2*sqrt(mp*mp + E*E)*sqrt(mp*mp + p*p) +
          2*sqrt(mp*mp + E*E)*sqrt(mp*mp + p*p)*cos(theta);
}

void createHistograms(TTree* tree, const char* outDir) {
	// Declare kinematic variables to read from the tree
	double e_p, e_theta, e_phi, p1_p, p1_theta, p1_phi;
	double p2_p, p2_theta, p2_phi, p3_p, p3_theta, p3_phi;
	// Declare derived variables to read from the tree
	double Mh, Mh12, Mh13, Mh23;
	double Mx, Mx1, Mx2, Mx3, Mx12, Mx13, Mx23;
	// Declare variables for asymmetry calculations
	int helicity;
	double beam_pol, phi1, phi2, phi3, phi12, phi13, phi23, DepW, DepA;
	double x, z13, z2x, t13, t2x;
	// Set branch addresses
	tree->SetBranchAddress("helicity", &helicity);
	tree->SetBranchAddress("beam_pol", &beam_pol);
    tree->SetBranchAddress("e_p", &e_p);
    tree->SetBranchAddress("e_theta", &e_theta);
    tree->SetBranchAddress("e_phi", &e_phi);
    tree->SetBranchAddress("p1_p", &p1_p);
    tree->SetBranchAddress("p1_theta", &p1_theta);
    tree->SetBranchAddress("p1_phi", &p1_phi);
    tree->SetBranchAddress("p2_p", &p2_p);
    tree->SetBranchAddress("p2_theta", &p2_theta);
    tree->SetBranchAddress("p2_phi", &p2_phi);
    tree->SetBranchAddress("p3_p", &p3_p);
    tree->SetBranchAddress("p3_theta", &p3_theta);
    tree->SetBranchAddress("p3_phi", &p3_phi);
    tree->SetBranchAddress("phi1", &phi1);
    tree->SetBranchAddress("phi2", &phi2);
    tree->SetBranchAddress("phi3", &phi3);
    tree->SetBranchAddress("phi12", &phi12);
    tree->SetBranchAddress("phi13", &phi13);
    tree->SetBranchAddress("phi23", &phi23);
    tree->SetBranchAddress("DepW", &DepW);
    tree->SetBranchAddress("DepA", &DepA);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("z13", &z13);

    // Declare new variables to store missing particle information
	float px_p, px_theta, px_phi, Mx1x, Mx2x, Mx3x, Mh1x, Mh2x, Mh3x;

	// Define initial state 4-momentum (10.1998 GeV electron beam and stationary proton)
    TLorentzVector p_initial(0, 0, 10.1998, 10.1998 + 0.938); // (px, py, pz, E)
	for (Long64_t entry = 0; entry < tree->GetEntries(); entry++) {
		tree->GetEntry(entry);

    	// Create 4-momentum vectors for final state particles
        TLorentzVector p_e, p1, p2, p3;
        p_e.SetXYZM(e_p*sin(e_theta)*cos(e_phi), e_p*sin(e_theta)*sin(e_phi), 
        	e_p*cos(e_theta), 0.511e-3);
        p1.SetXYZM(p1_p*sin(p1_theta)*cos(p1_phi), p1_p*sin(p1_theta)*sin(p1_phi), 
        	p1_p*cos(p1_theta), 0.139570);
        p2.SetXYZM(p2_p*sin(p2_theta)*cos(p2_phi), p2_p*sin(p2_theta)*sin(p2_phi), 
        	p2_p*cos(p2_theta), 0.139570);
        p3.SetXYZM(p3_p*sin(p3_theta)*cos(p3_phi), p3_p*sin(p3_theta)*sin(p3_phi), 
        	p3_p*cos(p3_theta), 0.938272);

        // Calculate 4-momentum of missing particle
        TLorentzVector p_x = p_initial - (p_e + p1 + p2 + p3);

        // Populate missing particle variables
        px_p = p_x.P();
        px_theta = p_x.Theta();
        px_phi = p_x.Phi();

        // Calculate missing mass variables
        Mx1x = (p_initial - (p_e + p1)).M();
        Mx2x = (p_initial - (p_e + p2)).M();
        Mx3x = (p_initial - (p_e + p3)).M();

        // Calculate invariant mass variables
        Mh1x = (p1 + p_x).M();
        Mh2x = (p2 + p_x).M();
        Mh3x = (p3 + p_x).M();

        cout << Mh2x << endl;

	}
}

void rhominusDeltaplusplus(std::string root_file_path) {
    gStyle->SetCanvasColor(0);

    TFile* file = new TFile(root_file_path.c_str(), "READ");

    if (!file->IsOpen()) {
        cout << "Error opening ROOT file (is the location correct?). Exiting." << endl;
    }

    TTree* tree = (TTree*)file->Get("PhysicsEvents");

    if (!tree) {
        cout << "Error getting trees from ROOT files." << endl;
    }

    createHistograms(tree, "output");

    file->Close(); delete file;
}