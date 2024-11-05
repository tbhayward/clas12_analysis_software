void radiative_mx_distribution(const char* filename) {
    // Include necessary headers
    #include <TFile.h>
    #include <TTree.h>
    #include <TH1F.h>
    #include <TCanvas.h>
    #include <TLegend.h>
    #include <TSystem.h>
    #include <TStyle.h>
    #include <iostream>

    using namespace std;

    // Open the file
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        cout << "Error opening file " << filename << endl;
        return;
    }

    // Get the tree
    TTree* tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        cout << "Error: PhysicsEvents tree not found in file " << filename << endl;
        return;
    }

    // Set branch addresses
    Float_t e_p, e_theta, e_phi;
    Float_t p_p, p_theta, p_phi;

    tree->SetBranchAddress("e_p", &e_p);
    tree->SetBranchAddress("e_theta", &e_theta);
    tree->SetBranchAddress("e_phi", &e_phi);

    tree->SetBranchAddress("p_p", &p_p);
    tree->SetBranchAddress("p_theta", &p_theta);
    tree->SetBranchAddress("p_phi", &p_phi);

    // Constants
    Double_t m_e = 0.000511; // Electron mass in GeV/c^2
    Double_t M_p = 0.93827;  // Proton mass in GeV/c^2

    // Beam energies
    Double_t beam_energies[] = {10.6, 10.3, 10.0, 8.0, 6.0};
    const int n_beam = 5;

    // Histograms
    TH1F* hMx2[n_beam];

    for (int i = 0; i < n_beam; ++i) {
        TString hname = Form("hMx2_%d", i);
        TString htitle = Form("Mx2 for Beam Energy %.1f GeV", beam_energies[i]);
        hMx2[i] = new TH1F(hname, htitle, 200, -4, 8); // Adjust bins and range as needed
    }

    // Colors
    Color_t colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta};

    // Loop over events
    Long64_t nentries = tree->GetEntries();

    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
        tree->GetEntry(ientry);

        // Get outgoing electron and proton four-momenta
        // Angles are in radians

        for (int j = 0; j < n_beam; ++j) {
            Double_t E_beam = beam_energies[j];

            // Initial electron
            Double_t E_electron_initial = sqrt(E_beam * E_beam + m_e * m_e);
            Double_t p_electron_initial_x = 0;
            Double_t p_electron_initial_y = 0;
            Double_t p_electron_initial_z = sqrt(E_beam * E_beam - m_e * m_e);

            // Initial proton
            Double_t E_proton_initial = M_p;
            Double_t p_proton_initial_x = 0;
            Double_t p_proton_initial_y = 0;
            Double_t p_proton_initial_z = 0;

            // Total initial four-momentum
            Double_t E_initial = E_electron_initial + E_proton_initial;
            Double_t p_initial_x = p_electron_initial_x + p_proton_initial_x;
            Double_t p_initial_y = p_electron_initial_y + p_proton_initial_y;
            Double_t p_initial_z = p_electron_initial_z + p_proton_initial_z;

            // Scattered electron
            Double_t E_electron_scattered = sqrt(e_p * e_p + m_e * m_e);
            Double_t p_electron_scattered_x = e_p * sin(e_theta) * cos(e_phi);
            Double_t p_electron_scattered_y = e_p * sin(e_theta) * sin(e_phi);
            Double_t p_electron_scattered_z = e_p * cos(e_theta);

            // Scattered proton
            Double_t E_proton_scattered = sqrt(p_p * p_p + M_p * M_p);
            Double_t p_proton_scattered_x = p_p * sin(p_theta) * cos(p_phi);
            Double_t p_proton_scattered_y = p_p * sin(p_theta) * sin(p_phi);
            Double_t p_proton_scattered_z = p_p * cos(p_theta);

            // Total final four-momentum
            Double_t E_final = E_electron_scattered + E_proton_scattered;
            Double_t p_final_x = p_electron_scattered_x + p_proton_scattered_x;
            Double_t p_final_y = p_electron_scattered_y + p_proton_scattered_y;
            Double_t p_final_z = p_electron_scattered_z + p_proton_scattered_z;

            // Missing four-momentum
            Double_t E_missing = E_initial - E_final;
            Double_t p_missing_x = p_initial_x - p_final_x;
            Double_t p_missing_y = p_initial_y - p_final_y;
            Double_t p_missing_z = p_initial_z - p_final_z;

            // Compute Mx2
            Double_t Mx2 = E_missing * E_missing - (p_missing_x * p_missing_x + p_missing_y * p_missing_y + p_missing_z * p_missing_z);

            // Fill histogram
            hMx2[j]->Fill(Mx2);
        }
    }

    // Create canvas
    TCanvas* c1 = new TCanvas("c1", "Radiative Mx2 Distribution", 1200, 800);
    c1->Divide(1,2);

    // First pad
    c1->cd(1);
    hMx2[0]->GetXaxis()->SetRangeUser(-1, 1);
    hMx2[0]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    hMx2[0]->GetYaxis()->SetTitle("Counts");
    hMx2[0]->SetLineColor(colors[0]);
    hMx2[0]->SetStats(0);
    hMx2[0]->Draw();

    for (int i = 1; i < n_beam; ++i) {
        hMx2[i]->SetLineColor(colors[i]);
        hMx2[i]->SetStats(0);
        hMx2[i]->Draw("same");
    }

    // Create legend
    TLegend* legend1 = new TLegend(0.7,0.7,0.9,0.9);
    for (int i = 0; i < n_beam; ++i) {
        TString entry = Form("%.1f GeV", beam_energies[i]);
        legend1->AddEntry(hMx2[i], entry, "l");
    }
    legend1->Draw();

    // Second pad
    c1->cd(2);
    hMx2[0]->GetXaxis()->SetRangeUser(-3, 6);
    hMx2[0]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    hMx2[0]->GetYaxis()->SetTitle("Counts");
    hMx2[0]->SetLineColor(colors[0]);
    hMx2[0]->SetStats(0);
    hMx2[0]->Draw();

    for (int i = 1; i < n_beam; ++i) {
        hMx2[i]->SetLineColor(colors[i]);
        hMx2[i]->SetStats(0);
        hMx2[i]->Draw("same");
    }

    // Create legend for second pad
    TLegend* legend2 = new TLegend(0.7,0.7,0.9,0.9);
    for (int i = 0; i < n_beam; ++i) {
        TString entry = Form("%.1f GeV", beam_energies[i]);
        legend2->AddEntry(hMx2[i], entry, "l");
    }
    legend2->Draw();

    // Save the canvas
    // Check if output directory exists
    if (gSystem->AccessPathName("output")) {
        gSystem->Exec("mkdir output");
    }
    c1->SaveAs("output/rad_distribution.png");

    // Clean up
    file->Close();
}