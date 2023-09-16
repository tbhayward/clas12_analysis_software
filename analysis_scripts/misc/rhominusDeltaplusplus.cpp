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
    {"Mh13", {200, 0.00, 1.50}},
    {"Mh23", {200, 0.00, 3.50}},
    {"Mh1x", {200, 0.00, 1.50}},
    {"Mh2x", {200, 0.00, 1.50}},
    {"Mh3x", {200, 0.00, 1.50}}
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
    TTreeReaderValue<double> Mx(dataReader, "Mx");
    TTreeReaderValue<double> Mx12(dataReader, "Mx12");
    TTreeReaderValue<double> Mx13(dataReader, "Mx13");
    TTreeReaderValue<double> Mx23(dataReader, "Mx23");
    TTreeReaderValue<double> Mh12(dataReader, "Mh12");
    TTreeReaderValue<double> Mh13(dataReader, "Mh13");
    TTreeReaderValue<double> Mh23(dataReader, "Mh23");

    // Declare new variables to store missing particle information
	float px_p, px_theta, px_phi, Mx1x, Mx2x, Mx3x, Mh1x, Mh2x, Mh3x;

	// Define initial state 4-momentum (10.1998 GeV electron beam and stationary proton)
    TLorentzVector p_initial(0, 0, 10.1998, 10.1998 + 0.938); // (px, py, pz, E)

    TCanvas canvas("Invariant Masses", "Canvas", 1600, 600);  
    // Width doubled for side-by-side panels
    TPad *pad1 = new TPad("pad1", "The pad with the function",0.0,0.0,0.33,1.0,21);
    pad1->SetLeftMargin(0.2); pad1->SetBottomMargin(0.2);
    pad1->SetFillColor(0);  // Set the fill color to white for pad1
    pad1->Draw();
    pad1->cd();  // Set current pad to pad1

    // histograms
    HistConfig configMh12 = histConfigs["Mh12"];
    TH1F histMh12("Mh12", "", configMh12.bins, configMh12.min, configMh12.max);
    HistConfig configMh13 = histConfigs["Mh13"];
    TH1F histMh13("Mh13", "", configMh13.bins, configMh13.min, configMh13.max);
    HistConfig configMh23 = histConfigs["Mh23"];
    TH1F histMh23("Mh23", "", configMh23.bins, configMh23.min, configMh23.max);
    HistConfig configMh1x = histConfigs["Mh1x"];
    TH1F histMh1x("Mh1x", "", configMh1x.bins, configMh1x.min, configMh1x.max);
    HistConfig configMh2x = histConfigs["Mh2x"];
    TH1F histMh2x("Mh2x", "", configMh2x.bins, configMh2x.min, configMh2x.max);
    HistConfig configMh3x = histConfigs["Mh3x"];
    TH1F histMh3x("Mh3x", "", configMh3x.bins, configMh3x.min, configMh3x.max);

	int counter = 0;
	while (dataReader.Next()) {
		counter++;
		if (*Mx < 0 || *Mx12 < 0 || *Mx13 < 0 || *Mx23 < 0) { continue; }

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

        histMh12.Fill(*Mh12); histMh13.Fill(*Mh13); histMh23.Fill(*Mh23);
        histMh1x.Fill(Mh1x); histMh2x.Fill(Mh2x); histMh3x.Fill(Mh3x); 

        cout << counter << " " << Mh2x << endl;
	}
	dataReader.Restart();  // Reset the TTreeReader at the end of the function

	histMh12.SetLineColor(kBlack); histMh13.SetLineColor(kBlack); 
    histMh23.SetLineColor(kBlack); histMh1x.SetLineColor(kBlack);
    histMh2x.SetLineColor(kBlack); histMh3x.SetLineColor(kBlack);
    histMh12.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
    histMh12.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
    histMh13.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
    histMh13.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
    histMh23.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
    histMh23.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
    histMh1x.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
    histMh1x.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
    histMh2x.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
    histMh2x.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
    histMh3x.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
    histMh3x.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size

    histMh12.Draw(""); histMh12.SetStats(0); 
    histMh12.GetXaxis()->SetTitle("{M}_{h}");
    histMh12.GetYaxis()->SetTitle("Counts");

	// Save the canvas
    canvas.SaveAs("output.png", outDir, branchName);
	delete pad1;
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