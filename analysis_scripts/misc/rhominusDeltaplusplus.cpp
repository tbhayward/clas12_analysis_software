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

// Function to set the style of a histogram
void setHistStyle(TH1* hist, const char* xTitle, const char* yTitle) {
	hist->SetLineColor(kBlack);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleSize(0.07);
    hist->GetYaxis()->SetTitleSize(0.07);
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);
    hist->SetStats(0);
}

// Function to set the style of a canvas and its pads
void setCanvasStyle(TCanvas &canvas, int nCols, int nRows) {
    canvas.Divide(nCols, nRows);
    for (int i = 1; i <= nCols * nRows; ++i) {
        canvas.cd(i);
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.185);
        gPad->SetRightMargin(0.175);
    }
}

struct HistConfig {
    int bins;
    double min;
    double max;
};

std::map<std::string, HistConfig> histConfigs = {
    {"z1", {500, 0.00, 1.00}},
    {"p1_p", {500, 0.00, 3.50}},
    {"p3_p", {500, 0.00, 3.50}},
    {"p13_theta", {500, 0.00, 45.00}},
    {"Mh12", {500, 0.00, 2.00}},
    {"Mh1x", {500, 0.00, 2.00}},
    {"Mh2x", {500, 0.00, 2.00}},
    {"Mh13", {500, 1.00, 2.00}},
    {"Mh23", {500, 1.00, 3.00}},
    {"Mh1x", {500, 0.00, 2.50}},
    {"Mh2x", {500, 0.00, 2.50}},
    {"Mh3x", {500, 0.00, 2.50}},
    {"Mx", {500, 0.00, 2.00}},
    {"Mx13", {500, 0.00, 2.50}},
    {"Mx2x", {500, 1.00, 3.50}},
    {"xF13", {500, -1.00, 1.00}},
    {"Delta_E2EX", {500, -3.00, 3.00}}
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
    TTreeReaderValue<double> z1(dataReader, "z1");
    TTreeReaderValue<double> z13(dataReader, "z13");
    TTreeReaderValue<double> Mx(dataReader, "Mx");
    TTreeReaderValue<double> Mx12(dataReader, "Mx12");
    TTreeReaderValue<double> Mx13(dataReader, "Mx13");
    TTreeReaderValue<double> Mx23(dataReader, "Mx23");
    TTreeReaderValue<double> Mh12(dataReader, "Mh12");
    TTreeReaderValue<double> Mh13(dataReader, "Mh13");
    TTreeReaderValue<double> Mh23(dataReader, "Mh23");
    TTreeReaderValue<double> xF13(dataReader, "xF13");
    TTreeReaderValue<double> Delta_phi13(dataReader, "Delta_phi13");

    // Declare new variables to store missing particle information
	float px_p, px_theta, px_phi, Mx1x, Mx2x, Mx3x, Mh1x, Mh2x, Mh3x;

	// Define initial state 4-momentum (10.1998 GeV electron beam and stationary proton)
    TLorentzVector p_initial(0, 0, 10.1998, 10.1998 + 0.938); // (px, py, pz, E)

    // Create canvas and set its style
    TCanvas canvas("Invariant Masses", "Canvas", 1600, 1000);
    setCanvasStyle(canvas, 3, 2);

	// Create test canvas and set its style
    TCanvas test_canvas("Cuts Test", "Canvas", 1600, 1000);
    setCanvasStyle(test_canvas, 3, 2);

    // 1D histograms
    HistConfig configMh13 = histConfigs["Mh13"];
    TH1F histMh13("Mh13", "", configMh13.bins, configMh13.min, configMh13.max);
    TH1F histMh13_cuts("Mh13, cuts", "", configMh13.bins, configMh13.min, configMh13.max);
    //
    HistConfig configMh2x = histConfigs["Mh2x"];
    TH1F histMh2x("Mh2x", "", configMh2x.bins, configMh2x.min, configMh2x.max);
    TH1F histMh2x_cuts("Mh2x, cuts", "", configMh2x.bins, configMh2x.min, configMh2x.max);
    //
    HistConfig configMx = histConfigs["Mx"];
    TH1F histMx("Mx", "", configMx.bins, configMx.min, configMx.max);

    // Add a 2D histogram for Mh13 vs Mh2x
    TH2F histMh13vsMh2x("Mh13vsMh2x", "", configMh2x.bins/10, configMh2x.min, configMh2x.max,
    	configMh13.bins/10, configMh13.min, configMh13.max); 

    ////////////////////////////////////////
    // test histograms
    HistConfig configMh12 = histConfigs["Mh12"];
    TH1F histMh12("Mh12", "", configMh12.bins, configMh12.min, configMh12.max);

    HistConfig configMh1x = histConfigs["Mh1x"];
    TH1F histMh1x("Mh1x", "", configMh1x.bins, configMh1x.min, configMh1x.max);

    HistConfig configp13_theta = histConfigs["p13_theta"];
    TH2F histMh13vsp13_theta("Mh13vsp13_theta", "", configp13_theta.bins/5, configp13_theta.min, 
    	configp13_theta.max,
    	configMh13.bins/5, configMh13.min, configMh13.max);

    TH2F histMh13vsMh12("Mh13vsMh12", "", configMh12.bins/5, configMh12.min, configMh12.max,
    	configMh13.bins/5, configMh13.min, configMh13.max);

    // HistConfig configMh1x = histConfigs["Mh1x"];
    TH2F histMh13vsMh1x("Mh13vsMh1x", "", configMh1x.bins/5, configMh1x.min, configMh1x.max,
    	configMh13.bins/5, configMh13.min, configMh13.max);

    HistConfig configDelta_E2EX = histConfigs["Delta_E2EX"];
    TH2F histMh13vsDelta_E2EX("Mh13vsE2EX", "", configDelta_E2EX.bins/5, configDelta_E2EX.min, 
    	configDelta_E2EX.max,
    	configMh13.bins/5, configMh13.min, configMh13.max);


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

        // Fill histograms without cuts
        histMh13.Fill(*Mh13); histMh2x.Fill(Mh2x); histMx.Fill(*Mx); 


        if (*Mx < 0.30 && (*Delta_phi13 < 0.1 || Delta_phi13 > 0.5)) {
        	histMh13_cuts.Fill(*Mh13); histMh2x_cuts.Fill(Mh2x);

        	histMh13vsMh2x.Fill(Mh2x, *Mh13);

        	histMh13vsp13_theta.Fill(TMath::RadToDeg()*(p1+p3).Theta(), *Mh13);

        	histMh13vsMh12.Fill(*Mh12, *Mh13);

        	histMh13vsMh1x.Fill(Mh1x, *Mh13);

        	histMh12.Fill(*Mh12);

        	histMh1x.Fill(Mh1x);

        	histMh13vsDelta_E2EX.Fill(p1.E()-p2.E(), *Mh13);
        }

	}
	dataReader.Restart();  // Reset the TTreeReader at the end of the function


	// Draw histograms on the canvas sub-pads
    canvas.cd(1);
    setHistStyle(&histMh13, "#it{M}_{h(#pi^{+}p)} (GeV)", "Counts");
    histMh13.Draw(); // Draw Mh13 in first pad
    //
    canvas.cd(2);
    setHistStyle(&histMh2x, "#it{M}_{h(#pi^{-}X)} (GeV)", "Counts");
    histMh2x.Draw(); // Draw Mh2x in second pad
    //
    canvas.cd(3);
    setHistStyle(&histMx, "#it{M}_{X(ep -> e'#pi^{+}#pi^{-}p[X])} (GeV)", "Counts");
    histMx.Draw(); // Draw Mx in third pad
    //
    canvas.cd(4);
    histMh13_cuts.SetTitle("#it{M}_{X} < 0.3 GeV");
    setHistStyle(&histMh13_cuts, "#it{M}_{h(#pi^{+}p)} (GeV)", "Counts");
    histMh13_cuts.Draw(); // Draw Mh13_cuts in fourth pad
    //
    canvas.cd(5);
    setHistStyle(&histMh2x_cuts, "#it{M}_{h(#pi^{-}X)} (GeV)", "Counts");
    histMh2x_cuts.SetTitle("#it{M}_{X} < 0.3 GeV");
    histMh2x_cuts.Draw(); // Draw Mh2x_cuts in fifth pad
    //
    // Draw the 2D histogram in the sixth panel
    canvas.cd(6);
    setHistStyle(&histMh13vsMh2x, "#it{M}_{h(#pi^{-}X)} (GeV)", "#it{M}_{h(#pi^{+}p)} (GeV)");
    histMh13vsMh2x.SetTitle("#it{M}_{X} < 0.3 GeV");
    histMh13vsMh2x.Draw("colz");  // Draw using color to represent the bin content
    histMh13vsMh2x.SetStats(0);
	
	// Save the canvas
    canvas.SaveAs(Form("%s/%s.png", outDir, "output"));





    // test histograms
	// Draw histograms on the canvas sub-pads
    test_canvas.cd(1);
    setHistStyle(&histMh12, "#it{M}_{h(#pi^{+}#pi^{-})} (GeV)", "Counts");
    histMh12.Draw(); 
    //
    test_canvas.cd(2);
    setHistStyle(&histMh1x, "#it{M}_{h(#pi^{+}X)} (GeV)", "Counts");
    histMh1x.Draw(); 
    //
    test_canvas.cd(3);
    setHistStyle(&histMh13vsp13_theta, "#it{#theta}_{#pi^{+}p}", "#it{M}_{h(#pi^{+}p)} (GeV)");
    histMh13vsp13_theta.Draw("colz"); 
    //
    test_canvas.cd(4);
    setHistStyle(&histMh13vsMh12, "#it{M}_{h(#pi^{+}#pi^{-})} (GeV)","#it{M}_{h(#pi^{+}p)} (GeV)");
    histMh13vsMh12.Draw("colz"); 
    //
    test_canvas.cd(5);
    setHistStyle(&histMh13vsMh1x, "#it{M}_{h(#pi^{+}X)} (GeV)","#it{M}_{h(#pi^{+}p)} (GeV)");
    histMh13vsMh1x.Draw("colz"); 
    //
    test_canvas.cd(6);
    setHistStyle(&histMh13vsDelta_E2EX, "#it{E}_{#pi^{-}} - #it{E}_{X} (GeV)", 
    	"#it{M}_{h(#pi^{+}p)} (GeV)");
    histMh13vsDelta_E2EX.Draw("colz"); 
    //

    // Save the canvas
    test_canvas.SaveAs(Form("%s/%s.png", outDir, "output_test"));
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