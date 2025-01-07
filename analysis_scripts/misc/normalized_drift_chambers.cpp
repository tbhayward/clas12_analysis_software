/*****************************************************************************
 * File: normalized_drift_chambers.cpp
 *
 * Compile with:
 *   g++ normalized_drift_chambers.cpp -o normalized_drift_chambers \
 *       $(root-config --cflags --libs)
 *
 * Run with:
 *   ./normalized_drift_chambers [NeventsSIDISDVCS] [dataFile] [mcFile]
 *
 * Notes on the arguments:
 *   1) NeventsSIDISDVCS: 
 *      - If 0, processes *all* events 
 *      - If omitted, also processes all events by default.
 *      - Otherwise, processes up to NeventsSIDISDVCS events.
 *   2) dataFile: 
 *      - If provided, uses this file for the SIDIS-DVCS chain (no merging).
 *      - If omitted, merges the three hard-coded SIDIS-DVCS files.
 *   3) mcFile:
 *      - If provided, uses this file for the CLASDIS chain (skips the hard-coded file).
 *      - If omitted, uses the single hard-coded CLASDIS file.
 *
 * Explanation:
 *   - Creates TChains and reads the user-specified (or default/hard-coded)
 *     data (SIDIS-DVCS) and mc (CLASDIS) files into "PhysicsEvents" TTrees.
 *   - Applies filters on PID, trajectory edges, and -9999 track checks.
 *   - Creates 2D histograms per region, accumulates them for each dataset,
 *     then takes the ratio SIDIS-DVCS/CLASDIS in each region.
 *   - Draws a 1×3 canvas of ratio plots with a log Z-axis and saves it.
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
    // 1. Parse input arguments
    //--------------------------------------------------------------------------
    Long64_t maxEventsSIDIS = -1; // Default: process all events
    
    if (argc > 1) {
        // If the user specified a first argument, convert it to a number
        maxEventsSIDIS = std::stoll(argv[1]);
        // If that number is 0, treat it as "process all" (same as -1)
        if (maxEventsSIDIS == 0) {
            maxEventsSIDIS = -1;
        }
        std::cout << "First argument (NeventsSIDISDVCS) = " 
                  << argv[1] 
                  << "; interpreted as maxEventsSIDIS = " 
                  << maxEventsSIDIS 
                  << std::endl;
    } else {
        std::cout << "No event limit argument specified; will process *all* events." 
                  << std::endl;
    }

    // If second argument is provided, it's the user-specified SIDIS-DVCS data file
    bool useUserDataFile = false;
    std::string userDataFile;
    if (argc > 2) {
        useUserDataFile = true;
        userDataFile = argv[2];
        std::cout << "Second argument provided (dataFile): " << userDataFile << std::endl;
    }

    // If third argument is provided, it's the user-specified CLASDIS MC file
    bool useUserMCFile = false;
    std::string userMCFile;
    if (argc > 3) {
        useUserMCFile = true;
        userMCFile = argv[3];
        std::cout << "Third argument provided (mcFile): " << userMCFile << std::endl;
    }

    //--------------------------------------------------------------------------
    // 2. Create TChains and add the input files
    //--------------------------------------------------------------------------
    TChain chain_sidisdvcs("PhysicsEvents");
    TChain chain_clasdis("PhysicsEvents");

    if (useUserDataFile) {
        // If a second argument was provided, use that single file for data
        chain_sidisdvcs.Add(userDataFile.c_str());
    } else {
        // Otherwise, revert to the original 3-file merge
        chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration.root");
        chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration_1.root");
        chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration_2.root");
    }

    if (useUserMCFile) {
        // If a third argument was provided, use that single file for MC
        chain_clasdis.Add(userMCFile.c_str());
    } else {
        // Otherwise, revert to the original hard-coded single file
        chain_clasdis.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/clasdis_rgc_su22_inb_neutron_calibration.root");
    }

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
    chain_sidisdvcs.SetBranchAddress("traj_edge_6", &traj_edge_6);
    chain_sidisdvcs.SetBranchAddress("traj_edge_18",&traj_edge_18);
    chain_sidisdvcs.SetBranchAddress("traj_edge_36",&traj_edge_36);

    // --- CLASDIS chain
    chain_clasdis.SetBranchAddress("particle_pid",  &particle_pid);
    chain_clasdis.SetBranchAddress("traj_x_6",      &traj_x_6);
    chain_clasdis.SetBranchAddress("traj_y_6",      &traj_y_6);
    chain_clasdis.SetBranchAddress("traj_x_18",     &traj_x_18);
    chain_clasdis.SetBranchAddress("traj_y_18",     &traj_y_18);
    chain_clasdis.SetBranchAddress("traj_x_36",     &traj_x_36);
    chain_clasdis.SetBranchAddress("traj_y_36",     &traj_y_36);
    chain_clasdis.SetBranchAddress("traj_edge_6",   &traj_edge_6);
    chain_clasdis.SetBranchAddress("traj_edge_18",  &traj_edge_18);
    chain_clasdis.SetBranchAddress("traj_edge_36",  &traj_edge_36);

    //--------------------------------------------------------------------------
    // 4. Create 2D histograms for each region, for each dataset
    //--------------------------------------------------------------------------
    // (User has changed the x,y ranges for each region here. 
    //  Below is the version in your snippet, but you can adjust as needed.)
    int nbins = 250;
    TH2D* h2_r1_sidisdvcs = new TH2D("h2_r1_sidisdvcs", "Region 1; x; y", 
                                     nbins, -180.0, 180, nbins, -180.0, 180);
    TH2D* h2_r2_sidisdvcs = new TH2D("h2_r2_sidisdvcs", "Region 2; x; y", 
                                     nbins, -280.0, 280, nbins, -280.0, 280);
    TH2D* h2_r3_sidisdvcs = new TH2D("h2_r3_sidisdvcs", "Region 3; x; y", 
                                     nbins, -320.0, 320, nbins, -320.0, 320);

    TH2D* h2_r1_clasdis = new TH2D("h2_r1_clasdis", "Region 1; x; y", 
                                   nbins, -180.0, 180, nbins, -180.0, 180.0);
    TH2D* h2_r2_clasdis = new TH2D("h2_r2_clasdis", "Region 2; x; y", 
                                   nbins, -280.0, 280, nbins, -280.0, 280);
    TH2D* h2_r3_clasdis = new TH2D("h2_r3_clasdis", "Region 3; x; y", 
                                   nbins, -320.0, 320, nbins, -320.0, 320);

    //--------------------------------------------------------------------------
    // 5. Fill histograms from the SIDIS-DVCS chain
    //--------------------------------------------------------------------------
    Long64_t nEntriesSIDIS = chain_sidisdvcs.GetEntries();
    // If user gave a specific limit that is >0, use it, else process all
    if (maxEventsSIDIS > 0 && maxEventsSIDIS < nEntriesSIDIS) {
        nEntriesSIDIS = maxEventsSIDIS; 
    }

    for (Long64_t i = 0; i < nEntriesSIDIS; i++) {
        chain_sidisdvcs.GetEntry(i);

        // Filter 1: must be an electron
        if (particle_pid != 11) continue; // #endif
        // (Or commented line for protons if you want to revert.)

        // Filter 2: must satisfy edge_6 > 5, edge_18 > 5, edge_36 > 10
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

        // // Filter 1: must be an electron
        // if (particle_pid != 11) continue; // #endif
        // Filter 1: must be a proton
        if (particle_pid != 2212) continue; // #endif

        // Filter 2: must satisfy edge_6 > 5, edge_18 > 5, edge_36 > 10
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
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();  // Log scale for the Z axis
    h2_r1_ratio->Draw("COLZ");

    // Region 2
    c->cd(2);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r2_ratio->Draw("COLZ");

    // Region 3
    c->cd(3);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r3_ratio->Draw("COLZ");

    //--------------------------------------------------------------------------
    // 9. Save to file
    //--------------------------------------------------------------------------
    c->SaveAs("output/normalized_drift_chambers_test.png");

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