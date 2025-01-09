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
 *      - If provided, uses this file for the CLASDIS chain (skips the
 *        hard-coded file).
 *      - If omitted, uses the single hard-coded CLASDIS file.
 *
 * Explanation:
 *   - Creates TChains and reads the user-specified (or default/hard-coded)
 *     data (SIDIS-DVCS) and mc (CLASDIS) files into "PhysicsEvents" TTrees.
 *   - Applies filters on PID, trajectory edges, and -9999 track checks.
 *   - Creates 2D histograms per region (both uncut & cut) for each dataset,
 *     then normalizes them individually before taking the ratio:
 *        ratio = (normalized SIDIS-DVCS) / (normalized CLASDIS).
 *   - We produce three 2D canvases:
 *       (1) Uncut (no circles)
 *       (2) Same data, but overlay red circles on each region
 *       (3) A "cut" version, only data inside the circle radii is used
 *   - Also creates 1D histograms of the edge variables themselves
 *     (traj_edge_6, traj_edge_18, traj_edge_36), normalizes them,
 *     and takes their ratio. This 1D edge ratio plot is not changed by
 *     the circular cuts or circles overlay; it's the same as before.
 *
 *   - Finally, we produce a 4th "rad-phi" canvas with 1×3 subplots. Each region
 *     shows the ratio of normalized SIDIS vs CLASDIS, but now binned in
 *     radius (rho) and phi, where phi=0 means pointing left (x<0,y=0).
 *
 * Output images are written into:
 *   output/normalization/
 *
 *   - normalized_drift_chambers_uncut.png
 *   - normalized_drift_chambers_circle.png
 *   - normalized_drift_chambers_cut.png
 *   - normalized_drift_chambers_edges.png
 *   - normalized_drift_chambers_radphi.png
 *****************************************************************************/

#include <iostream>
#include <string>
#include <cmath>          // for sqrt, atan2, etc.
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "TEllipse.h"     // Use TEllipse for drawing circles

// Helper function to unify the color scale for three TH2D histograms
void SetSame2DScale(TH2D* h1, TH2D* h2, TH2D* h3)
{
    double minVal = +1e9;
    double maxVal = -1e9;

    // Find global min & max across the three
    for (auto* hh : {h1, h2, h3}) {
        double locMin = hh->GetMinimum();
        double locMax = hh->GetMaximum();
        if (locMin < minVal) minVal = locMin;
        if (locMax > maxVal) maxVal = locMax;
    }

    // Apply them
    for (auto* hh : {h1, h2, h3}) {
        hh->SetMinimum(minVal);
        hh->SetMaximum(maxVal);
    }
}

// Helper function to unify the y-scale for three TH1D histograms
void SetSame1DScale(TH1D* h1, TH1D* h2, TH1D* h3)
{
    double minVal = +1e9;
    double maxVal = -1e9;

    for (auto* hh : {h1, h2, h3}) {
        double locMin = hh->GetMinimum();
        double locMax = hh->GetMaximum();
        if (locMin < minVal) minVal = locMin;
        if (locMax > maxVal) maxVal = locMax;
    }

    for (auto* hh : {h1, h2, h3}) {
        hh->SetMinimum(minVal);
        hh->SetMaximum(maxVal);
    }
}

int main(int argc, char** argv)
{
    //--------------------------------------------------------------------------
    // 1. Parse input arguments
    //--------------------------------------------------------------------------
    Long64_t maxEventsSIDIS = -1; // Default: process all events
    
    if (argc > 1) {
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
        chain_sidisdvcs.Add(userDataFile.c_str());
    } else {
        chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration.root");
        chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration_1.root");
        chain_sidisdvcs.Add("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/calibration/sidisdvcs_rgc_su22_inb_calibration_2.root");
    }

    if (useUserMCFile) {
        chain_clasdis.Add(userMCFile.c_str());
    } else {
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

    // // Filter line (commented out so we can revert):
    // // if (particle_pid != 11) continue; // electron
    // // if (particle_pid != 2212) continue; // proton

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
    // 4. Create 2D histograms for each region, for each dataset (uncut & cut)
    //--------------------------------------------------------------------------
    int nbins2D = 200;
    // -- Uncut
    TH2D* h2_r1_sidisdvcs = new TH2D("h2_r1_sidisdvcs", "Region 1; x; y", 
                                     nbins2D, -180.0, 180, nbins2D, -180.0, 180);
    TH2D* h2_r2_sidisdvcs = new TH2D("h2_r2_sidisdvcs", "Region 2; x; y", 
                                     nbins2D, -280.0, 280, nbins2D, -280.0, 280);
    TH2D* h2_r3_sidisdvcs = new TH2D("h2_r3_sidisdvcs", "Region 3; x; y", 
                                     nbins2D, -450.0, 450, nbins2D, -450.0, 450);

    TH2D* h2_r1_clasdis = new TH2D("h2_r1_clasdis", "Region 1; x; y", 
                                   nbins2D, -180.0, 180, nbins2D, -180.0, 180.0);
    TH2D* h2_r2_clasdis = new TH2D("h2_r2_clasdis", "Region 2; x; y", 
                                   nbins2D, -280.0, 280, nbins2D, -280.0, 280);
    TH2D* h2_r3_clasdis = new TH2D("h2_r3_clasdis", "Region 3; x; y", 
                                   nbins2D, -450.0, 450, nbins2D, -450.0, 450);

    // -- Cut
    TH2D* h2_r1_sidisdvcs_cut = new TH2D("h2_r1_sidisdvcs_cut", "Region 1 (cut); x; y", 
                                         nbins2D, -180.0, 180, nbins2D, -180.0, 180);
    TH2D* h2_r2_sidisdvcs_cut = new TH2D("h2_r2_sidisdvcs_cut", "Region 2 (cut); x; y", 
                                         nbins2D, -280.0, 280, nbins2D, -280.0, 280);
    TH2D* h2_r3_sidisdvcs_cut = new TH2D("h2_r3_sidisdvcs_cut", "Region 3 (cut); x; y", 
                                         nbins2D, -450.0, 450, nbins2D, -450.0, 450);

    TH2D* h2_r1_clasdis_cut = new TH2D("h2_r1_clasdis_cut", "Region 1 (cut); x; y", 
                                       nbins2D, -180.0, 180, nbins2D, -180.0, 180.0);
    TH2D* h2_r2_clasdis_cut = new TH2D("h2_r2_clasdis_cut", "Region 2 (cut); x; y", 
                                       nbins2D, -280.0, 280, nbins2D, -280.0, 280);
    TH2D* h2_r3_clasdis_cut = new TH2D("h2_r3_clasdis_cut", "Region 3 (cut); x; y", 
                                       nbins2D, -450.0, 450, nbins2D, -450.0, 450);

    //--------------------------------------------------------------------------
    // 4a. Create 1D histograms for the edge variables
    //--------------------------------------------------------------------------
    int nbins1D = 125;
    double edgeMin = 0.0;
    double edgeMax = 125.0;

    TH1D* h1_edge_r1_sidis = new TH1D("h1_edge_r1_sidis", "Region 1; edge (cm); counts", 
                                      nbins1D, edgeMin, edgeMax);
    TH1D* h1_edge_r2_sidis = new TH1D("h1_edge_r2_sidis", "Region 2; edge (cm); counts", 
                                      nbins1D, edgeMin, edgeMax);
    TH1D* h1_edge_r3_sidis = new TH1D("h1_edge_r3_sidis", "Region 3; edge (cm); counts", 
                                      nbins1D, edgeMin, edgeMax);

    TH1D* h1_edge_r1_clas = new TH1D("h1_edge_r1_clas", "Region 1; edge (cm); counts", 
                                     nbins1D, edgeMin, edgeMax);
    TH1D* h1_edge_r2_clas = new TH1D("h1_edge_r2_clas", "Region 2; edge (cm); counts", 
                                     nbins1D, edgeMin, edgeMax);
    TH1D* h1_edge_r3_clas = new TH1D("h1_edge_r3_clas", "Region 3; edge (cm); counts", 
                                     nbins1D, edgeMin, edgeMax);

    // --------------------------------------------------------------------------
    // 4b. Create new 2D histograms for (radius, phi) for each region
    //     Each region has its own radial binning, from 40 up to the region's max.
    //     Number of radial bins = (rMax - 40). We'll use 180 bins for phi in [0,360).
    // --------------------------------------------------------------------------

    // ---- Region 1: r in [40, 170]
    int nbinsR1    = 170 - 20;  // = 130
    double rMin1   = 20.0;
    double rMax1   = 170.0;

    // ---- Region 2: r in [40, 250]
    int nbinsR2    = 250 - 20;  // = 210
    double rMin2   = 20.0;
    double rMax2   = 250.0;

    // ---- Region 3: r in [40, 420]
    int nbinsR3    = 420 - 20;  // = 380
    double rMin3   = 20.0;
    double rMax3   = 420.0;

    // Phi binning for all regions
    int    nbinsPhi = 180;
    double phiMin   = 0.0;
    double phiMax   = 360.0;

    // ------------------- SIDIS rad-phi histograms -------------------
    TH2D* h2_r1_sidis_radphi = new TH2D(
        "h2_r1_sidis_radphi",
        "Region 1 rad-phi; radius (cm); phi (deg)",
        nbinsR1, rMin1, rMax1,
        nbinsPhi, phiMin, phiMax
    );

    TH2D* h2_r2_sidis_radphi = new TH2D(
        "h2_r2_sidis_radphi",
        "Region 2 rad-phi; radius (cm); phi (deg)",
        nbinsR2, rMin2, rMax2,
        nbinsPhi, phiMin, phiMax
    );

    TH2D* h2_r3_sidis_radphi = new TH2D(
        "h2_r3_sidis_radphi",
        "Region 3 rad-phi; radius (cm); phi (deg)",
        nbinsR3, rMin3, rMax3,
        nbinsPhi, phiMin, phiMax
    );

    // ------------------- CLAS rad-phi histograms -------------------
    TH2D* h2_r1_clas_radphi = new TH2D(
        "h2_r1_clas_radphi",
        "Region 1 rad-phi; radius (cm); phi (deg)",
        nbinsR1, rMin1, rMax1,
        nbinsPhi, phiMin, phiMax
    );

    TH2D* h2_r2_clas_radphi = new TH2D(
        "h2_r2_clas_radphi",
        "Region 2 rad-phi; radius (cm); phi (deg)",
        nbinsR2, rMin2, rMax2,
        nbinsPhi, phiMin, phiMax
    );

    TH2D* h2_r3_clas_radphi = new TH2D(
        "h2_r3_clas_radphi",
        "Region 3 rad-phi; radius (cm); phi (deg)",
        nbinsR3, rMin3, rMax3,
        nbinsPhi, phiMin, phiMax
    );

    //--------------------------------------------------------------------------
    // Radii for the 3 circles
    //--------------------------------------------------------------------------
    const double r1 = 140.0;  // Region 1 radius
    const double r2 = 215.0;  // Region 2 radius
    const double r3 = 290.0;  // Region 3 radius
    // const double r3 = 1000.0;  // Region 3 radius

    //--------------------------------------------------------------------------
    // 5. Fill histograms from the SIDIS-DVCS chain
    //--------------------------------------------------------------------------
    Long64_t nEntriesSIDIS = chain_sidisdvcs.GetEntries();
    if (maxEventsSIDIS > 0 && maxEventsSIDIS < nEntriesSIDIS) {
        nEntriesSIDIS = maxEventsSIDIS; 
    }

    for (Long64_t i = 0; i < nEntriesSIDIS; i++) {
        chain_sidisdvcs.GetEntry(i);

        // // Filter line (commented out so we can revert):
        // if (particle_pid != 11) continue; // electron
        if (particle_pid != 2212) continue; // proton

        // if (traj_edge_6 <= 5 || traj_edge_18 <= 5 || traj_edge_36 <= 10) continue; // #endif
        if (traj_edge_6 <= 3 || traj_edge_18 <= 3 || traj_edge_36 <= 7) continue; // #endif

        // Region 1 (x_6,y_6)
        if (traj_x_6 != -9999) {
            h2_r1_sidisdvcs->Fill(traj_x_6, traj_y_6);
            double dist1 = std::sqrt(traj_x_6*traj_x_6 + traj_y_6*traj_y_6);
            if (dist1 < r1) {
                h2_r1_sidisdvcs_cut->Fill(traj_x_6, traj_y_6);
            }

            // Fill rad-phi
            double rawAngle = std::atan2(traj_y_6, traj_x_6); 
            double phi = rawAngle - M_PI;  // shift so left is 0 deg
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180./M_PI;
            h2_r1_sidis_radphi->Fill(dist1, phiDeg);
        }

        // Region 2 (x_18,y_18)
        if (traj_x_18 != -9999) {
            h2_r2_sidisdvcs->Fill(traj_x_18, traj_y_18);
            double dist2 = std::sqrt(traj_x_18*traj_x_18 + traj_y_18*traj_y_18);
            if (dist2 < r2) {
                h2_r2_sidisdvcs_cut->Fill(traj_x_18, traj_y_18);
            }

            double rawAngle = std::atan2(traj_y_18, traj_x_18); 
            double phi = rawAngle - M_PI;
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180./M_PI;
            h2_r2_sidis_radphi->Fill(dist2, phiDeg);
        }

        // Region 3 (x_36,y_36)
        if (traj_x_36 != -9999) {
            h2_r3_sidisdvcs->Fill(traj_x_36, traj_y_36);
            double dist3 = std::sqrt(traj_x_36*traj_x_36 + traj_y_36*traj_y_36);
            if (dist3 < r3) {
                h2_r3_sidisdvcs_cut->Fill(traj_x_36, traj_y_36);
            }

            double rawAngle = std::atan2(traj_y_36, traj_x_36); 
            double phi = rawAngle - M_PI;
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180./M_PI;
            h2_r3_sidis_radphi->Fill(dist3, phiDeg);
        }

        // Fill 1D edge histograms
        h1_edge_r1_sidis->Fill(traj_edge_6);
        h1_edge_r2_sidis->Fill(traj_edge_18);
        h1_edge_r3_sidis->Fill(traj_edge_36);
    }

    //--------------------------------------------------------------------------
    // 6. Fill histograms from the CLASDIS chain (no event limit)
    //--------------------------------------------------------------------------
    Long64_t nEntriesCLASDIS = chain_clasdis.GetEntries();
    for (Long64_t i = 0; i < nEntriesCLASDIS; i++) {
        chain_clasdis.GetEntry(i);

        // // Filter line (commented out so we can revert):
        // if (particle_pid != 11) continue; // electron
        if (particle_pid != 2212) continue; // proton

        // if (traj_edge_6 <= 5 || traj_edge_18 <= 5 || traj_edge_36 <= 10) continue; // #endif
        if (traj_edge_6 <= 3 || traj_edge_18 <= 3 || traj_edge_36 <= 7) continue; // #endif

        // Region 1
        if (traj_x_6 != -9999) {
            h2_r1_clasdis->Fill(traj_x_6, traj_y_6);
            double dist1 = std::sqrt(traj_x_6*traj_x_6 + traj_y_6*traj_y_6);
            if (dist1 < r1) {
                h2_r1_clasdis_cut->Fill(traj_x_6, traj_y_6);
            }

            double rawAngle = std::atan2(traj_y_6, traj_x_6); 
            double phi = rawAngle - M_PI;
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180./M_PI;
            h2_r1_clas_radphi->Fill(dist1, phiDeg);
        }

        // Region 2
        if (traj_x_18 != -9999) {
            h2_r2_clasdis->Fill(traj_x_18, traj_y_18);
            double dist2 = std::sqrt(traj_x_18*traj_x_18 + traj_y_18*traj_y_18);
            if (dist2 < r2) {
                h2_r2_clasdis_cut->Fill(traj_x_18, traj_y_18);
            }

            double rawAngle = std::atan2(traj_y_18, traj_x_18); 
            double phi = rawAngle - M_PI;
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180./M_PI;
            h2_r2_clas_radphi->Fill(dist2, phiDeg);
        }

        // Region 3
        if (traj_x_36 != -9999) {
            h2_r3_clasdis->Fill(traj_x_36, traj_y_36);
            double dist3 = std::sqrt(traj_x_36*traj_x_36 + traj_y_36*traj_y_36);
            if (dist3 < r3) {
                h2_r3_clasdis_cut->Fill(traj_x_36, traj_y_36);
            }

            double rawAngle = std::atan2(traj_y_36, traj_x_36); 
            double phi = rawAngle - M_PI;
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180./M_PI;
            h2_r3_clas_radphi->Fill(dist3, phiDeg);
        }

        // Fill 1D edges
        h1_edge_r1_clas->Fill(traj_edge_6);
        h1_edge_r2_clas->Fill(traj_edge_18);
        h1_edge_r3_clas->Fill(traj_edge_36);
    }

    //--------------------------------------------------------------------------
    // 7. Normalize each 2D histogram (uncut), then compute ratio
    //--------------------------------------------------------------------------
    double int_r1_sidis = h2_r1_sidisdvcs->Integral();
    if (int_r1_sidis > 0) h2_r1_sidisdvcs->Scale(1.0 / int_r1_sidis);
    double int_r2_sidis = h2_r2_sidisdvcs->Integral();
    if (int_r2_sidis > 0) h2_r2_sidisdvcs->Scale(1.0 / int_r2_sidis);
    double int_r3_sidis = h2_r3_sidisdvcs->Integral();
    if (int_r3_sidis > 0) h2_r3_sidisdvcs->Scale(1.0 / int_r3_sidis);

    double int_r1_clas = h2_r1_clasdis->Integral();
    if (int_r1_clas > 0) h2_r1_clasdis->Scale(1.0 / int_r1_clas);
    double int_r2_clas = h2_r2_clasdis->Integral();
    if (int_r2_clas > 0) h2_r2_clasdis->Scale(1.0 / int_r2_clas);
    double int_r3_clas = h2_r3_clasdis->Integral();
    if (int_r3_clas > 0) h2_r3_clasdis->Scale(1.0 / int_r3_clas);

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
    // 7a. Unify scale for the three uncut ratio histograms
    //--------------------------------------------------------------------------
    SetSame2DScale(h2_r1_ratio, h2_r2_ratio, h2_r3_ratio);

    //--------------------------------------------------------------------------
    // 8. Draw the uncut ratio histograms (Canvas 1)
    //--------------------------------------------------------------------------
    gStyle->SetOptStat(0);

    TCanvas *c2D_uncut = new TCanvas("c2D_uncut", "Normalized 2D Drift Chambers (Uncut)", 1800, 600);
    c2D_uncut->Divide(3,1);

    c2D_uncut->cd(1);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r1_ratio->Draw("COLZ");

    c2D_uncut->cd(2);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r2_ratio->Draw("COLZ");

    c2D_uncut->cd(3);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r3_ratio->Draw("COLZ");

    // c2D_uncut->SaveAs("output/normalization/electron_normalized_drift_chambers_uncut.png");
    c2D_uncut->SaveAs("output/normalization/proton_normalized_drift_chambers_uncut.png");

    //--------------------------------------------------------------------------
    // 9. Draw the same ratio histograms with circles (Canvas 2)
    //--------------------------------------------------------------------------
    TCanvas *c2D_circle = new TCanvas("c2D_circle", "Normalized 2D Drift Chambers (Circle Overlay)", 1800, 600);
    c2D_circle->Divide(3,1);

    int circleColor    = kRed;
    int circleLineWidth = 2;
    int circleFillStyle = 0;

    // Region 1
    c2D_circle->cd(1);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r1_ratio->Draw("COLZ");
    {
      TEllipse* circ1 = new TEllipse(0.0, 0.0, r1, r1);
      circ1->SetLineColor(circleColor);
      circ1->SetLineWidth(circleLineWidth);
      circ1->SetFillStyle(circleFillStyle);
      circ1->Draw("same");
    }

    // Region 2
    c2D_circle->cd(2);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r2_ratio->Draw("COLZ");
    {
      TEllipse* circ2 = new TEllipse(0.0, 0.0, r2, r2);
      circ2->SetLineColor(circleColor);
      circ2->SetLineWidth(circleLineWidth);
      circ2->SetFillStyle(circleFillStyle);
      circ2->Draw("same");
    }

    // Region 3
    c2D_circle->cd(3);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r3_ratio->Draw("COLZ");
    {
      TEllipse* circ3 = new TEllipse(0.0, 0.0, r3, r3);
      circ3->SetLineColor(circleColor);
      circ3->SetLineWidth(circleLineWidth);
      circ3->SetFillStyle(circleFillStyle);
      circ3->Draw("same");
    }

    // c2D_circle->SaveAs("output/normalization/electron_normalized_drift_chambers_circle.png");
    c2D_circle->SaveAs("output/normalization/proton_normalized_drift_chambers_circle.png");

    //--------------------------------------------------------------------------
    // 10. Normalize the "cut" histograms, then ratio
    //--------------------------------------------------------------------------
    double int_r1_sidis_cut = h2_r1_sidisdvcs_cut->Integral();
    if (int_r1_sidis_cut > 0) h2_r1_sidisdvcs_cut->Scale(1.0 / int_r1_sidis_cut);
    double int_r2_sidis_cut = h2_r2_sidisdvcs_cut->Integral();
    if (int_r2_sidis_cut > 0) h2_r2_sidisdvcs_cut->Scale(1.0 / int_r2_sidis_cut);
    double int_r3_sidis_cut = h2_r3_sidisdvcs_cut->Integral();
    if (int_r3_sidis_cut > 0) h2_r3_sidisdvcs_cut->Scale(1.0 / int_r3_sidis_cut);

    double int_r1_clas_cut = h2_r1_clasdis_cut->Integral();
    if (int_r1_clas_cut > 0) h2_r1_clasdis_cut->Scale(1.0 / int_r1_clas_cut);
    double int_r2_clas_cut = h2_r2_clasdis_cut->Integral();
    if (int_r2_clas_cut > 0) h2_r2_clasdis_cut->Scale(1.0 / int_r2_clas_cut);
    double int_r3_clas_cut = h2_r3_clasdis_cut->Integral();
    if (int_r3_clas_cut > 0) h2_r3_clasdis_cut->Scale(1.0 / int_r3_clas_cut);

    TH2D* h2_r1_ratio_cut = (TH2D*) h2_r1_sidisdvcs_cut->Clone("h2_r1_ratio_cut");
    h2_r1_ratio_cut->SetTitle("Region 1 (cut)");
    h2_r1_ratio_cut->Divide(h2_r1_clasdis_cut);

    TH2D* h2_r2_ratio_cut = (TH2D*) h2_r2_sidisdvcs_cut->Clone("h2_r2_ratio_cut");
    h2_r2_ratio_cut->SetTitle("Region 2 (cut)");
    h2_r2_ratio_cut->Divide(h2_r2_clasdis_cut);

    TH2D* h2_r3_ratio_cut = (TH2D*) h2_r3_sidisdvcs_cut->Clone("h2_r3_ratio_cut");
    h2_r3_ratio_cut->SetTitle("Region 3 (cut)");
    h2_r3_ratio_cut->Divide(h2_r3_clasdis_cut);

    //--------------------------------------------------------------------------
    // 10a. Unify scale for the three "cut" ratio histograms
    //--------------------------------------------------------------------------
    SetSame2DScale(h2_r1_ratio_cut, h2_r2_ratio_cut, h2_r3_ratio_cut);

    //--------------------------------------------------------------------------
    // 11. Draw the "cut" ratio histograms (Canvas 3)
    //--------------------------------------------------------------------------
    TCanvas *c2D_cut = new TCanvas("c2D_cut", "Normalized 2D Drift Chambers (Cut)", 1800, 600);
    c2D_cut->Divide(3,1);

    c2D_cut->cd(1);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r1_ratio_cut->Draw("COLZ");

    c2D_cut->cd(2);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r2_ratio_cut->Draw("COLZ");

    c2D_cut->cd(3);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r3_ratio_cut->Draw("COLZ");

    // c2D_cut->SaveAs("output/normalization/electron_normalized_drift_chambers_cut.png");
    c2D_cut->SaveAs("output/normalization/proton_normalized_drift_chambers_cut.png");

    //--------------------------------------------------------------------------
    // 12. Normalize each 1D edge histogram individually, then compute ratio
    //--------------------------------------------------------------------------
    double e1_sidis = h1_edge_r1_sidis->Integral();
    if (e1_sidis > 0) h1_edge_r1_sidis->Scale(1.0 / e1_sidis);
    double e2_sidis = h1_edge_r2_sidis->Integral();
    if (e2_sidis > 0) h1_edge_r2_sidis->Scale(1.0 / e2_sidis);
    double e3_sidis = h1_edge_r3_sidis->Integral();
    if (e3_sidis > 0) h1_edge_r3_sidis->Scale(1.0 / e3_sidis);

    double e1_clas = h1_edge_r1_clas->Integral();
    if (e1_clas > 0) h1_edge_r1_clas->Scale(1.0 / e1_clas);
    double e2_clas = h1_edge_r2_clas->Integral();
    if (e2_clas > 0) h1_edge_r2_clas->Scale(1.0 / e2_clas);
    double e3_clas = h1_edge_r3_clas->Integral();
    if (e3_clas > 0) h1_edge_r3_clas->Scale(1.0 / e3_clas);

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
    // 12a. Unify scale for the three 1D edge ratio histograms
    //--------------------------------------------------------------------------
    SetSame1DScale(h1_edge_r1_ratio, h1_edge_r2_ratio, h1_edge_r3_ratio);

    //--------------------------------------------------------------------------
    // 13. Draw the 1×3 ratio plots for the edge variables
    //--------------------------------------------------------------------------
    TCanvas *c1D = new TCanvas("c1D", "Edge Ratios", 1800, 600);
    c1D->Divide(3,1);

    // Region 1 edge
    c1D->cd(1);
    h1_edge_r1_ratio->GetXaxis()->SetTitle("edge (cm)");
    h1_edge_r1_ratio->GetYaxis()->SetTitle("normalized density ratio");
    h1_edge_r1_ratio->Draw("HIST");

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

    // c1D->SaveAs("output/normalization/electron_normalized_drift_chambers_edges.png");
    c1D->SaveAs("output/normalization/proton_normalized_drift_chambers_edges.png");

    //--------------------------------------------------------------------------
    // 14. Now produce the final "rad-phi" ratio for each region (uncut).
    //--------------------------------------------------------------------------
    // Normalize SIDIS rad-phi
    double int_r1_sidis_radphi = h2_r1_sidis_radphi->Integral();
    if (int_r1_sidis_radphi > 0) h2_r1_sidis_radphi->Scale(1.0/int_r1_sidis_radphi);

    double int_r2_sidis_radphi = h2_r2_sidis_radphi->Integral();
    if (int_r2_sidis_radphi > 0) h2_r2_sidis_radphi->Scale(1.0/int_r2_sidis_radphi);

    double int_r3_sidis_radphi = h2_r3_sidis_radphi->Integral();
    if (int_r3_sidis_radphi > 0) h2_r3_sidis_radphi->Scale(1.0/int_r3_sidis_radphi);

    // Normalize CLAS rad-phi
    double int_r1_clas_radphi = h2_r1_clas_radphi->Integral();
    if (int_r1_clas_radphi > 0) h2_r1_clas_radphi->Scale(1.0/int_r1_clas_radphi);

    double int_r2_clas_radphi = h2_r2_clas_radphi->Integral();
    if (int_r2_clas_radphi > 0) h2_r2_clas_radphi->Scale(1.0/int_r2_clas_radphi);

    double int_r3_clas_radphi = h2_r3_clas_radphi->Integral();
    if (int_r3_clas_radphi > 0) h2_r3_clas_radphi->Scale(1.0/int_r3_clas_radphi);

    // Ratio
    TH2D* h2_r1_ratio_radphi = (TH2D*) h2_r1_sidis_radphi->Clone("h2_r1_ratio_radphi");
    h2_r1_ratio_radphi->SetTitle("Region 1 (radius vs phi)");
    h2_r1_ratio_radphi->Divide(h2_r1_clas_radphi);

    TH2D* h2_r2_ratio_radphi = (TH2D*) h2_r2_sidis_radphi->Clone("h2_r2_ratio_radphi");
    h2_r2_ratio_radphi->SetTitle("Region 2 (radius vs phi)");
    h2_r2_ratio_radphi->Divide(h2_r2_clas_radphi);

    TH2D* h2_r3_ratio_radphi = (TH2D*) h2_r3_sidis_radphi->Clone("h2_r3_ratio_radphi");
    h2_r3_ratio_radphi->SetTitle("Region 3 (radius vs phi)");
    h2_r3_ratio_radphi->Divide(h2_r3_clas_radphi);

    // Set same scale
    SetSame2DScale(h2_r1_ratio_radphi, h2_r2_ratio_radphi, h2_r3_ratio_radphi);

    // Draw final rad-phi ratio canvas (Canvas 4)
    TCanvas* cRadPhi = new TCanvas("cRadPhi","Ratio of radius vs phi",1800,600);
    cRadPhi->Divide(3,1);

    cRadPhi->cd(1);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r1_ratio_radphi->Draw("COLZ");

    cRadPhi->cd(2);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r2_ratio_radphi->Draw("COLZ");

    cRadPhi->cd(3);
    gPad->SetRightMargin(0.20);
    gPad->SetLeftMargin(0.20);
    gPad->SetLogz();
    h2_r3_ratio_radphi->Draw("COLZ");

    // cRadPhi->SaveAs("output/normalization/electron_normalized_drift_chambers_radphi.png");
    cRadPhi->SaveAs("output/normalization/proton_normalized_drift_chambers_radphi.png");

    //--------------------------------------------------------------------------
    // Cleanup
    //--------------------------------------------------------------------------
    delete c2D_uncut;
    delete c2D_circle;
    delete c2D_cut;
    delete c1D;
    delete cRadPhi;

    // 2D hists (uncut)
    delete h2_r1_sidisdvcs;
    delete h2_r2_sidisdvcs;
    delete h2_r3_sidisdvcs;
    delete h2_r1_clasdis;
    delete h2_r2_clasdis;
    delete h2_r3_clasdis;
    delete h2_r1_ratio;
    delete h2_r2_ratio;
    delete h2_r3_ratio;

    // 2D hists (cut)
    delete h2_r1_sidisdvcs_cut;
    delete h2_r2_sidisdvcs_cut;
    delete h2_r3_sidisdvcs_cut;
    delete h2_r1_clasdis_cut;
    delete h2_r2_clasdis_cut;
    delete h2_r3_clasdis_cut;
    delete h2_r1_ratio_cut;
    delete h2_r2_ratio_cut;
    delete h2_r3_ratio_cut;

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

    // Rad-phi hists
    delete h2_r1_sidis_radphi;
    delete h2_r2_sidis_radphi;
    delete h2_r3_sidis_radphi;
    delete h2_r1_clas_radphi;
    delete h2_r2_clas_radphi;
    delete h2_r3_clas_radphi;
    delete h2_r1_ratio_radphi;
    delete h2_r2_ratio_radphi;
    delete h2_r3_ratio_radphi;

    return 0;
}