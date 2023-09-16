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

struct HistConfig {
    int bins;
    double min;
    double max;
};

std::map<std::string, HistConfig> histConfigs = {
    {"Mh12", {200, 0.00, 3.00}},
    {"Mh13", {200, 0.00, 1.50}}
};

void createHistograms(TTreeReader &dataReader, const char* outDir) {
	// Declare derived variables to read from the tree
	double z2x, t13, t2x;
	// Declare reader locations
	TTreeReaderValue<int> helicity(dataReader, "helicity");
	TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
    TTreeReaderValue<double> e_p(dataReader, "e_p");
    TTreeReaderValue<double> e_theta(dataReader, "e_theta");
    TTreeReaderValue<double> e_phi(dataReader, "e_phi");
    TTreeReaderValue<double> p1_p(dataReader, "p1_p");
    TTreeReaderValue<double> p1_theta(dataReader, "p1_theta");
    TTreeReaderValue<double> p1_phi(dataReader, "p1_phi");
    TTreeReaderValue<double> p2_p(dataReader, "p2_p");
    TTreeReaderValue<double> p2_theta(dataReader, "p2_theta");
    TTreeReaderValue<double> p2_phi(dataReader, "p2_phi");
    TTreeReaderValue<double> p3_p(dataReader, "p3_p");
    TTreeReaderValue<double> p3_theta(dataReader, "p3_theta");
    TTreeReaderValue<double> p3_phi(dataReader, "p3_phi");
    TTreeReaderValue<double> phi1(dataReader, "phi1");
    TTreeReaderValue<double> phi2(dataReader, "phi2");
    TTreeReaderValue<double> phi3(dataReader, "phi3");
    TTreeReaderValue<double> phi12(dataReader, "phi12");
    TTreeReaderValue<double> phi13(dataReader, "phi13");
    TTreeReaderValue<double> phi23(dataReader, "phi23");
    TTreeReaderValue<double> DepW(dataReader, "DepW");
    TTreeReaderValue<double> DepA(dataReader, "DepA");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> z13(dataReader, "z");

    // Declare new variables to store missing particle information
	float px_p, px_theta, px_phi, Mx1x, Mx2x, Mx3x, Mh1x, Mh2x, Mh3x;

	// Define initial state 4-momentum (10.1998 GeV electron beam and stationary proton)
    TLorentzVector p_initial(0, 0, 10.1998, 10.1998 + 0.938); // (px, py, pz, E)

    // histograms
    HistConfig config = histConfigs[branchName];
    TH1F histMh12("Mh12", "", config.bins, config.min, config.max);
    TH1F histMh13("Mh13", "", config.bins, config.min, config.max);
	KinematicCuts kinematicCuts(dataReader);  // Create an instance of the KinematicCuts class
	int counter = 0;
	while (dataReader.Next()) {
		counter++;

    	// Create 4-momentum vectors for final state particles
        TLorentzVector p_e, p1, p2, p3;
        p_e.SetXYZM(*e_p*sin(*e_theta)*cos(*e_phi), *e_p*sin(*e_theta)*sin(*e_phi), 
        	*e_p*cos(*e_theta), 0.511e-3);
        p1.SetXYZM(*p1_p*sin(*p1_theta)*cos(*p1_phi), *p1_p*sin(*p1_theta)*sin(*p1_phi), 
        	*p1_p*cos(*p1_theta), 0.139570);
        p2.SetXYZM(*p2_p*sin(*p2_theta)*cos(*p2_phi), *p2_p*sin(*p2_theta)*sin(*p2_phi), 
        	*p2_p*cos(*p2_theta), 0.139570);
        p3.SetXYZM(*p3_p*sin(*p3_theta)*cos(*p3_phi), *p3_p*sin(*p3_theta)*sin(*p3_phi), 
        	*p3_p*cos(*p3_theta), 0.938272);

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

        cout << counter << " " << Mh2x << endl;

	}
	dataReader.Restart();  // Reset the TTreeReader at the end of the function
}

void rhominusDeltaplusplus(std::string root_file_path) {
	// Start the timer
	auto start_time = std::chrono::high_resolution_clock::now();
    gStyle->SetCanvasColor(0);

    TFile* file = new TFile(root_file_path.c_str(), "READ");

    if (!file->IsOpen()) {
        cout << "Error opening ROOT file (is the location correct?). Exiting." << endl;
    }

    TTree* tree = (TTree*)file->Get("PhysicsEvents");

    if (!tree) {
        cout << "Error getting trees from ROOT files." << endl;
    }

    TTreeReader dataReader(tree); // Create a TTreeReader for the data tree

    createHistograms(dataReader, "output");

    file->Close(); delete file;

    // Stop the timer
	auto end_time = std::chrono::high_resolution_clock::now();
	// Calculate the elapsed time in seconds and microseconds
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - 
	start_time).count();
	double seconds = duration / 1e6;
	// Convert to hours, minutes, and seconds
	int hours = static_cast<int>(seconds) / 3600;
	int remaining_time = static_cast<int>(seconds) % 3600;
	int mins = remaining_time / 60;
	int remaining_seconds = remaining_time % 60;
	// Print the elapsed time
	cout << "Time elapsed: ";
	cout << hours << " hours, " << mins << " mins, " << remaining_seconds << " seconds." << endl;
	return 0;
}