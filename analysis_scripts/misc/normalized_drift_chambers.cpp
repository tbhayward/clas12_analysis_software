/*****************************************************************************
 * File: normalized_drift_chambers.cpp
 *
 * Compile with:
 *   g++ normalized_drift_chambers.cpp -o normalized_drift_chambers \
 *       $(root-config --cflags --libs)
 *
 * Run with:
 *   ./normalized_drift_chambers [NeventsSIDISDVCS]
 *   (If NeventsSIDISDVCS is omitted, *all* SIDIS-DVCS events will be processed.)
 *
 * Explanation:
 *   - This program merges three SIDIS-DVCS ROOT files into one TChain and
 *     reads a separate CLASDIS ROOT file into another TChain.
 *   - It creates 2D histograms for each of the three drift chamber regions
 *     (Region 1, 2, 3) in both datasets, then takes the ratio of
 *     SIDISDVCS/CLASDIS and plots them side by side.
 *   - The Z-axis is set to log scale in each pad, each pad has extra right 
 *     padding for color scale, and we apply the cuts:
 *         (particle_pid == 11) &&
 *         (edge_1 > 5) && (edge_2 > 5) && (edge_3 > 10)
 *     before filling histograms.
 *****************************************************************************/

#include <iostream>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"

int main(int argc, char** argv)
{
    //--------------------------------------------------------------------------
    // 1. Parse optional input argument for number of SIDIS-DVCS events
    //--------------------------------------------------------------------------
    Long64_t maxEventsSIDIS = -1; // -1 means process all events
    if (argc > 1) {
        maxEventsSIDIS = std::stoll(argv[1]);
        std::cout << "Will process up to " << maxEventsSIDIS 
                  << " events for the SIDIS-DVCS data set." << std::endl;
    } else {
        std::cout << "No event limit specified; will process *all* SIDIS-DVCS events." 
                  << std::endl;
    }

    //--------------------------------------------------------------------------
    // 2. Create TChains and add the input files
    //--------------------------------------------------------------------------
    // SIDIS-DVCS (these three files will be merged)
    TChain chain_sidisdvcs("PhysicsEvents");
    chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration.root");
    chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration_1.root");
    chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration_2.root");

    // CLASDIS (separate file)
    TChain chain_clasdis("PhysicsEvents");
    chain_clasdis.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/clasdis_rgc_su22_inb_neutron_calibration.root");

    //--------------------------------------------------------------------------
    // 3. Set up branch addresses for the variables we need
    //--------------------------------------------------------------------------
    Int_t    particle_pid;
    Double_t traj_x_6,  traj_y_6;
    Double_t traj_x_18, traj_y_18;
    Double_t traj_x_36, traj_y_36;

    // Additional variables for the new cuts
    Double_t traj_edge_6, traj_edge_18, traj_edge_36;

    // --- SIDIS-DVCS chain
    chain_sidisdvcs.SetBranchAddress("particle_pid", &particle_pid);
    chain_sidisdvcs.SetBranchAddress("traj_x_6",    &traj_x_6);
    chain_sidisdvcs.SetBranchAddress("traj_y_6",    &traj_y_6);
    chain_sidisdvcs.SetBranchAddress("traj_x_18",   &traj_x_18);
    chain_sidisdvcs.SetBranchAddress("traj_y_18",   &traj_y_18);
    chain_sidisdvcs.SetBranchAddress("traj_x_36",   &traj_x_36);
    chain_sidisdvcs.SetBranchAddress("traj_y_36",   &traj_y_36);
    chain_sidisdvcs.SetBranchAddress("traj_edge_6",      &traj_edge_6);
    chain_sidisdvcs.SetBranchAddress("traj_edge_18",      &traj_edge_18);
    chain_sidisdvcs.SetBranchAddress("traj_edge_36",      &traj_edge_36);

    // --- CLASDIS chain
    chain_clasdis.SetBranchAddress("particle_pid",  &particle_pid);
    chain_clasdis.SetBranchAddress("traj_x_6",      &traj_x_6);
    chain_clasdis.SetBranchAddress("traj_y_6",      &traj_y_6);
    chain_clasdis.SetBranchAddress("traj_x_18",     &traj_x_18);
    chain_clasdis.SetBranchAddress("traj_y_18",     &traj_y_18);
    chain_clasdis.SetBranchAddress("traj_x_36",     &traj_x_36);
    chain_clasdis.SetBranchAddress("traj_y_36",     &traj_y_36);
    chain_clasdis.SetBranchAddress("traj_edge_6",        &traj_edge_6);
    chain_clasdis.SetBranchAddress("traj_edge_18",        &traj_edge_18);
    chain_clasdis.SetBranchAddress("traj_edge_36",        &traj_edge_36);

    //--------------------------------------------------------------------------
    // 4. Create 2D histograms for each region, for each dataset
    //--------------------------------------------------------------------------
    // Region 1: -300 to 300 for x and y
    // Region 2: -400 to 400 for x and y
    // Region 3: -600 to 600 for x and y
    int nbins = 250;

    TH2D* h2_r1_sidisdvcs = new TH2D("h2_r1_sidisdvcs", "Region 1; x; y", 
                                     nbins, -180.0, 300.0, nbins, -180.0, 300.0);
    TH2D* h2_r2_sidisdvcs = new TH2D("h2_r2_sidisdvcs", "Region 2; x; y", 
                                     nbins, -280.0, 400.0, nbins, -280.0, 400.0);
    TH2D* h2_r3_sidisdvcs = new TH2D("h2_r3_sidisdvcs", "Region 3; x; y", 
                                     nbins, -320.0, 600.0, nbins, -320.0, 600.0);

    TH2D* h2_r1_clasdis = new TH2D("h2_r1_clasdis", "Region 1; x; y", 
                                   nbins, -180.0, 300.0, nbins, -180.0, 300.0);
    TH2D* h2_r2_clasdis = new TH2D("h2_r2_clasdis", "Region 2; x; y", 
                                   nbins, -280.0, 400.0, nbins, -280.0, 400.0);
    TH2D* h2_r3_clasdis = new TH2D("h2_r3_clasdis", "Region 3; x; y", 
                                   nbins, -320.0, 600.0, nbins, -320.0, 600.0);

    //--------------------------------------------------------------------------
    // 5. Fill histograms from the SIDIS-DVCS chain
    //--------------------------------------------------------------------------
    Long64_t nEntriesSIDIS = chain_sidisdvcs.GetEntries();
    if (maxEventsSIDIS > 0 && maxEventsSIDIS < nEntriesSIDIS) {
        nEntriesSIDIS = maxEventsSIDIS; 
    }

    for (Long64_t i = 0; i < nEntriesSIDIS; i++) {
        chain_sidisdvcs.GetEntry(i);

        // Filter 1: must be an electron
        if (particle_pid != 11) continue; // #endif

        // Filter 2: must satisfy edge_1 > 5, edge_2 > 5, edge_3 > 10
        if (traj_edge_6 <= 5 || traj_edge_18 <= 5 || traj_edge_36 <= 10) continue; // #endif

        // Fill Region 1 if traj_x_6 != -9999
        if (traj_x_6 != -9999) {
            h2_r1_sidisdvcs->Fill(traj_x_6, traj_y_6);
        } //endfor

        // Fill Region 2 if traj_x_18 != -9999
        if (traj_x_18 != -9999) {
            h2_r2_sidisdvcs->Fill(traj_x_18, traj_y_18);
        } //endfor

        // Fill Region 3 if traj_x_36 != -9999
        if (traj_x_36 != -9999) {
            h2_r3_sidisdvcs->Fill(traj_x_36, traj_y_36);
        } //endfor
    } //endfor

    //--------------------------------------------------------------------------
    // 6. Fill histograms from the CLASDIS chain (no event limit specified)
    //--------------------------------------------------------------------------
    Long64_t nEntriesCLASDIS = chain_clasdis.GetEntries();
    for (Long64_t i = 0; i < nEntriesCLASDIS; i++) {
        chain_clasdis.GetEntry(i);

        // Filter 1: must be an electron
        if (particle_pid != 11) continue; // #endif

        // Filter 2: must satisfy edge_1 > 5, edge_2 > 5, edge_3 > 10
        if (traj_edge_6 <= 5 || traj_edge_18 <= 5 || traj_edge_36 <= 10) continue; // #endif

        // Fill Region 1 if traj_x_6 != -9999
        if (traj_x_6 != -9999) {
            h2_r1_clasdis->Fill(traj_x_6, traj_y_6);
        } //endfor

        // Fill Region 2 if traj_x_18 != -9999
        if (traj_x_18 != -9999) {
            h2_r2_clasdis->Fill(traj_x_18, traj_y_18);
        } //endfor

        // Fill Region 3 if traj_x_36 != -9999
        if (traj_x_36 != -9999) {
            h2_r3_clasdis->Fill(traj_x_36, traj_y_36);
        } //endfor
    } //endfor

    //--------------------------------------------------------------------------
    // 7. Compute ratio histograms = SIDISDVCS / CLASDIS for each region
    //--------------------------------------------------------------------------
    TH2D* h2_r1_ratio = (TH2D*) h2_r1_sidisdvcs->Clone("h2_r1_ratio");
    h2_r1_ratio->Divide(h2_r1_clasdis);
    h2_r1_ratio->SetTitle("Region 1");

    TH2D* h2_r2_ratio = (TH2D*) h2_r2_sidisdvcs->Clone("h2_r2_ratio");
    h2_r2_ratio->Divide(h2_r2_clasdis);
    h2_r2_ratio->SetTitle("Region 2");

    TH2D* h2_r3_ratio = (TH2D*) h2_r3_sidisdvcs->Clone("h2_r3_ratio");
    h2_r3_ratio->Divide(h2_r3_clasdis);
    h2_r3_ratio->SetTitle("Region 3");

    // Remove stat boxes
    gStyle->SetOptStat(0);

    //--------------------------------------------------------------------------
    // 8. Draw the three ratio histograms on a single canvas, 1 row × 3 columns
    //--------------------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "Normalized Drift Chambers", 1800, 600);
    c->Divide(3,1);

    // Region 1
    c->cd(1);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);   // Add padding on the right so the color scale doesn't clip
    // gPad->SetLogz();              // Log scale for the Z axis
    h2_r1_ratio->Draw("COLZ");

    // Region 2
    c->cd(2);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    // gPad->SetLogz();
    h2_r2_ratio->Draw("COLZ");

    // Region 3
    c->cd(3);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    // gPad->SetLogz();
    h2_r3_ratio->Draw("COLZ");

    //--------------------------------------------------------------------------
    // 9. Save to file
    //--------------------------------------------------------------------------
    c->SaveAs("output/normalized_drift_chambers.png");

    // Cleanup
    delete c;
    delete h2_r1_sidisdvcs;
    delete h2_r2_sidisdvcs;
    delete h2_r3_sidisdvcs;
    delete h2_r1_clasdis;
    delete h2_r2_clasdis;
    delete h2_r3_clasdis;

    return 0;
}