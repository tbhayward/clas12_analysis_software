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
 *     then normalizes them individually before taking the ratio in 2D:
 *        ratio = (normalized SIDIS-DVCS) / (normalized CLASDIS).
 *   - Also creates 1D histograms of the edge variables themselves
 *     (traj_edge_6, traj_edge_18, traj_edge_36), normalizes them, 
 *     and takes their ratio. Each region's edge variable is shown
 *     in a 1×3 layout on a second canvas.
 *****************************************************************************/

#include <iostream>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
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
    // The x,y ranges here reflect your snippet:
    int nbins2D = 250;
    TH2D* h2_r1_sidisdvcs = new TH2D("h2_r1_sidisdvcs", "Region 1; x; y", 
                                     nbins2D, -180.0, 180, nbins2D, -180.0, 180);
    TH2D* h2_r2_sidisdvcs = new TH2D("h2_r2_sidisdvcs", "Region 2; x; y", 
                                     nbins2D, -280.0, 280, nbins2D, -280.0, 280);
    TH2D* h2_r3_sidisdvcs = new TH2D("h2_r3_sidisdvcs", "Region 3; x; y", 
                                     nbins2D, -320.0, 320, nbins2D, -320.0, 320);

    TH2D* h2_r1_clasdis = new TH2D("h2_r1_clasdis", "Region 1; x; y", 
                                   nbins2D, -180.0, 180, nbins2D, -180.0, 180.0);
    TH2D* h2_r2_clasdis = new TH2D("h2_r2_clasdis", "Region 2; x; y", 
                                   nbins2D, -280.0, 280, nbins2D, -280.0, 280);
    TH2D* h2_r3_clasdis = new TH2D("h2_r3_clasdis", "Region 3; x; y", 
                                   nbins2D, -320.0, 320, nbins2D, -320.0, 320);

    //--------------------------------------------------------------------------
    // 4a. Create 1D histograms for the edge variables 
    //     (region 1 -> edge_6, region 2 -> edge_18, region 3 -> edge_36).
    //     Adjust binning/range if needed (example: 100 bins from 0 to 50 cm).
    //--------------------------------------------------------------------------
    int nbins1D = 100;
    double edgeMin = 0.0;
    double edgeMax = 100.0;

    // SIDIS-DVCS 1D
    TH1D* h1_edge_r1_sidis = new TH1D("h1_edge_r1_sidis", "Region 1; edge (cm); counts", 
                                      nbins1D, edgeMin, edgeMax);
    TH1D* h1_edge_r2_sidis = new TH1D("h1_edge_r2_sidis", "Region 2; edge (cm); counts", 
                                      nbins1D, edgeMin, edgeMax);
    TH1D* h1_edge_r3_sidis = new TH1D("h1_edge_r3_sidis", "Region 3; edge (cm); counts", 
                                      nbins1D, edgeMin, edgeMax);

    // CLASDIS 1D
    TH1D* h1_edge_r1_clas = new TH1D("h1_edge_r1_clas", "Region 1; edge (cm); counts", 
                                     nbins1D, edgeMin, edgeMax);
    TH1D* h1_edge_r2_clas = new TH1D("h1_edge_r2_clas", "Region 2; edge (cm); counts", 
                                     nbins1D, edgeMin, edgeMax);
    TH1D* h1_edge_r3_clas = new TH1D("h1_edge_r3_clas", "Region 3; edge (cm); counts", 
                                     nbins1D, edgeMin, edgeMax);

    //--------------------------------------------------------------------------
    // 5. Fill histograms from the SIDIS-DVCS chain
    //--------------------------------------------------------------------------
    Long64_t nEntriesSIDIS = chain_sidisdvcs.GetEntries();
    // If user gave a specific limit that is >0, use it; else process all
    if (maxEventsSIDIS > 0 && maxEventsSIDIS < nEntriesSIDIS) {
        nEntriesSIDIS = maxEventsSIDIS; 
    }

    for (Long64_t i = 0; i < nEntriesSIDIS; i++) {
        chain_sidisdvcs.GetEntry(i);

        // Filter 1: must be an electron
        if (particle_pid != 11) continue; // #endif

        // Filter 2: must satisfy edge_6 > 5, edge_18 > 5, edge_36 > 10
        if (traj_edge_6 <= 5 || traj_edge_18 <= 5 || traj_edge_36 <= 10) continue; // #endif

        // Fill 2D Region 1 if traj_x_6 != -9999
        if (traj_x_6 != -9999) {
            h2_r1_sidisdvcs->Fill(traj_x_6, traj_y_6);
        } //endfor

        // Fill 2D Region 2 if traj_x_18 != -9999
        if (traj_x_18 != -9999) {
            h2_r2_sidisdvcs->Fill(traj_x_18, traj_y_18);
        } //endfor

        // Fill 2D Region 3 if traj_x_36 != -9999
        if (traj_x_36 != -9999) {
            h2_r3_sidisdvcs->Fill(traj_x_36, traj_y_36);
        } //endfor

        // Fill 1D edge histograms
        // Region 1: traj_edge_6
        h1_edge_r1_sidis->Fill(traj_edge_6);

        // Region 2: traj_edge_18
        h1_edge_r2_sidis->Fill(traj_edge_18);

        // Region 3: traj_edge_36
        h1_edge_r3_sidis->Fill(traj_edge_36);
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

        // Fill 2D Region 1 if traj_x_6 != -9999
        if (traj_x_6 != -9999) {
            h2_r1_clasdis->Fill(traj_x_6, traj_y_6);
        } //endfor

        // Fill 2D Region 2 if traj_x_18 != -9999
        if (traj_x_18 != -9999) {
            h2_r2_clasdis->Fill(traj_x_18, traj_y_18);
        } //endfor

        // Fill 2D Region 3 if traj_x_36 != -9999
        if (traj_x_36 != -9999) {
            h2_r3_clasdis->Fill(traj_x_36, traj_y_36);
        } //endfor

        // Fill 1D edge histograms
        // Region 1: traj_edge_6
        h1_edge_r1_clas->Fill(traj_edge_6);

        // Region 2: traj_edge_18
        h1_edge_r2_clas->Fill(traj_edge_18);

        // Region 3: traj_edge_36
        h1_edge_r3_clas->Fill(traj_edge_36);
    } //endfor

    //--------------------------------------------------------------------------
    // 7. Normalize each 2D histogram individually, then compute ratio
    //--------------------------------------------------------------------------
    //--- Normalize the SIDIS-DVCS histograms
    double int_r1_sidis = h2_r1_sidisdvcs->Integral();
    if (int_r1_sidis > 0) h2_r1_sidisdvcs->Scale(1.0 / int_r1_sidis);

    double int_r2_sidis = h2_r2_sidisdvcs->Integral();
    if (int_r2_sidis > 0) h2_r2_sidisdvcs->Scale(1.0 / int_r2_sidis);

    double int_r3_sidis = h2_r3_sidisdvcs->Integral();
    if (int_r3_sidis > 0) h2_r3_sidisdvcs->Scale(1.0 / int_r3_sidis);

    //--- Normalize the CLASDIS histograms
    double int_r1_clas = h2_r1_clasdis->Integral();
    if (int_r1_clas > 0) h2_r1_clasdis->Scale(1.0 / int_r1_clas);

    double int_r2_clas = h2_r2_clasdis->Integral();
    if (int_r2_clas > 0) h2_r2_clasdis->Scale(1.0 / int_r2_clas);

    double int_r3_clas = h2_r3_clasdis->Integral();
    if (int_r3_clas > 0) h2_r3_clasdis->Scale(1.0 / int_r3_clas);

    // Now compute ratio = (normalized SIDIS) / (normalized CLASDIS)
    TH2D* h2_r1_ratio = (TH2D*) h2_r1_sidisdvcs->Clone("h2_r1_ratio");
    h2_r1_ratio->SetTitle("Region 1");
    h2_r1_ratio->Divide(h2_r1_clasdis);

    TH2D* h2_r2_ratio = (TH2D*) h2_r2_sidisdvcs->Clone("h2_r2_ratio");
    h2_r2_ratio->SetTitle("Region 2");
    h2_r2_ratio->Divide(h2_r2_clasdis);

    TH2D* h2_r3_ratio = (TH2D*) h2_r3_sidisdvcs->Clone("h2_r3_ratio");
    h2_r3_ratio->SetTitle("Region 3");
    h2_r3_ratio->Divide(h2_r3_clasdis);

    //--------------------------------------------------------------------------
    // 8. Draw the three ratio histograms on a single canvas, 1 row × 3 columns
    //--------------------------------------------------------------------------
    gStyle->SetOptStat(0);  // remove stat boxes

    TCanvas *c2D = new TCanvas("c2D", "Normalized 2D Drift Chambers", 1800, 600);
    c2D->Divide(3,1);

    // Region 1
    c2D->cd(1);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();  // Log scale for the Z axis
    h2_r1_ratio->Draw("COLZ");

    // Region 2
    c2D->cd(2);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r2_ratio->Draw("COLZ");

    // Region 3
    c2D->cd(3);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r3_ratio->Draw("COLZ");

    c2D->SaveAs("output/normalized_drift_chambers_test.png");


    //--------------------------------------------------------------------------
    // 9. Normalize each 1D edge histogram individually, then compute ratio
    //--------------------------------------------------------------------------
    //--- SIDIS
    double e1_sidis = h1_edge_r1_sidis->Integral();
    if (e1_sidis > 0) h1_edge_r1_sidis->Scale(1.0 / e1_sidis);
    double e2_sidis = h1_edge_r2_sidis->Integral();
    if (e2_sidis > 0) h1_edge_r2_sidis->Scale(1.0 / e2_sidis);
    double e3_sidis = h1_edge_r3_sidis->Integral();
    if (e3_sidis > 0) h1_edge_r3_sidis->Scale(1.0 / e3_sidis);

    //--- CLAS
    double e1_clas = h1_edge_r1_clas->Integral();
    if (e1_clas > 0) h1_edge_r1_clas->Scale(1.0 / e1_clas);
    double e2_clas = h1_edge_r2_clas->Integral();
    if (e2_clas > 0) h1_edge_r2_clas->Scale(1.0 / e2_clas);
    double e3_clas = h1_edge_r3_clas->Integral();
    if (e3_clas > 0) h1_edge_r3_clas->Scale(1.0 / e3_clas);

    // Ratios
    TH1D* h1_edge_r1_ratio = (TH1D*) h1_edge_r1_sidis->Clone("h1_edge_r1_ratio");
    h1_edge_r1_ratio->SetTitle("Region 1");
    h1_edge_r1_ratio->Divide(h1_edge_r1_clas);

    TH1D* h1_edge_r2_ratio = (TH1D*) h1_edge_r2_sidis->Clone("h1_edge_r2_ratio");
    h1_edge_r2_ratio->SetTitle("Region 2");
    h1_edge_r2_ratio->Divide(h1_edge_r2_clas);

    TH1D* h1_edge_r3_ratio = (TH1D*) h1_edge_r3_sidis->Clone("h1_edge_r3_ratio");
    h1_edge_r3_ratio->SetTitle("Region 3");
    h1_edge_r3_ratio->Divide(h1_edge_r3_clas);


    //--------------------------------------------------------------------------
    // 10. Draw the 1×3 ratio plots for the edge variables
    //--------------------------------------------------------------------------
    TCanvas *c1D = new TCanvas("c1D", "Edge Ratios", 1800, 600);
    c1D->Divide(3,1);

    // Region 1 edge
    c1D->cd(1);
    // Setting axis titles
    h1_edge_r1_ratio->GetXaxis()->SetTitle("edge (cm)");
    h1_edge_r1_ratio->GetYaxis()->SetTitle("normalized density ratio");
    h1_edge_r1_ratio->Draw("HIST");  // 1D ratio plot

    // Region 2 edge
    c1D->cd(2);
    h1_edge_r2_ratio->GetXaxis()->SetTitle("edge (cm)");
    h1_edge_r2_ratio->GetYaxis()->SetTitle("normalized density ratio");
    h1_edge_r2_ratio->Draw("HIST");

    // Region 3 edge
    c1D->cd(3);
    h1_edge_r3_ratio->GetXaxis()->SetTitle("edge (cm)");
    h1_edge_r3_ratio->GetYaxis()->SetTitle("normalized density ratio");
    h1_edge_r3_ratio->Draw("HIST");

    // Save the new edge-ratio canvas
    c1D->SaveAs("output/normalized_drift_chambers_edges.png");

    //--------------------------------------------------------------------------
    // Cleanup
    //--------------------------------------------------------------------------
    delete c2D;
    delete c1D;

    // 2D hists
    delete h2_r1_sidisdvcs;
    delete h2_r2_sidisdvcs;
    delete h2_r3_sidisdvcs;
    delete h2_r1_clasdis;
    delete h2_r2_clasdis;
    delete h2_r3_clasdis;
    delete h2_r1_ratio;
    delete h2_r2_ratio;
    delete h2_r3_ratio;

    // 1D hists
    delete h1_edge_r1_sidis;
    delete h1_edge_r2_sidis;
    delete h1_edge_r3_sidis;
    delete h1_edge_r1_clas;
    delete h1_edge_r2_clas;
    delete h1_edge_r3_clas;
    delete h1_edge_r1_ratio;
    delete h1_edge_r2_ratio;
    delete h1_edge_r3_ratio;

    return 0;
}