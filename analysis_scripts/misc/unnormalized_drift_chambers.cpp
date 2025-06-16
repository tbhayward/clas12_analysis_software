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
 *      - If provided, uses this file for the CLASDIS chain.
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
#include "TEllipse.h"

// ----------------------------------------------------------------------------
// Helper function to unify the color scale for three TH2D histograms
// ----------------------------------------------------------------------------
void SetSame2DScale(TH2D* h1, TH2D* h2, TH2D* h3) {
    double minVal = +1e9;
    double maxVal = -1e9;
    for (auto* hh : {h1, h2, h3}) {
        minVal = std::min(minVal, hh->GetMinimum());
        maxVal = std::max(maxVal, hh->GetMaximum());
    }
    for (auto* hh : {h1, h2, h3}) {
        hh->SetMinimum(minVal);
        hh->SetMaximum(maxVal);
    }
}

// ----------------------------------------------------------------------------
// Helper function to unify the y-scale for three TH1D histograms
// ----------------------------------------------------------------------------
void SetSame1DScale(TH1D* h1, TH1D* h2, TH1D* h3) {
    double minVal = +1e9;
    double maxVal = -1e9;
    for (auto* hh : {h1, h2, h3}) {
        minVal = std::min(minVal, hh->GetMinimum());
        maxVal = std::max(maxVal, hh->GetMaximum());
    }
    for (auto* hh : {h1, h2, h3}) {
        hh->SetMinimum(minVal);
        hh->SetMaximum(maxVal);
    }
}

int main(int argc, char** argv) {
    //--------------------------------------------------------------------------
    // 1. Parse input arguments
    //--------------------------------------------------------------------------
    Long64_t maxEventsSIDIS = -1; // default: all events
    if (argc > 1) {
        maxEventsSIDIS = std::stoll(argv[1]);
        if (maxEventsSIDIS == 0) maxEventsSIDIS = -1;
    }

    bool useUserDataFile = (argc > 2);
    std::string userDataFile = useUserDataFile ? argv[2] : "";

    bool useUserMCFile = (argc > 3);
    std::string userMCFile = useUserMCFile ? argv[3] : "";

    //--------------------------------------------------------------------------
    // 2. Create TChains and add the input files
    //--------------------------------------------------------------------------
    TChain chain_sidisdvcs("PhysicsEvents");
    TChain chain_clasdis  ("PhysicsEvents");

    if (useUserDataFile) {
        chain_sidisdvcs.Add(userDataFile.c_str());
    } else {
        chain_sidisdvcs.Add("/work/clas12/thayward/.../sidisdvcs_rgc_su22_inb_calibration.root");
        chain_sidisdvcs.Add("/work/clas12/thayward/.../sidisdvcs_rgc_su22_inb_calibration_1.root");
        chain_sidisdvcs.Add("/work/clas12/thayward/.../sidisdvcs_rgc_su22_inb_calibration_2.root");
    }

    if (useUserMCFile) {
        chain_clasdis.Add(userMCFile.c_str());
    } else {
        chain_clasdis.Add("/work/clas12/thayward/.../clasdis_rgc_su22_inb_neutron_calibration.root");
    }

    //--------------------------------------------------------------------------
    // 3. Set up branch addresses for the variables we need
    //--------------------------------------------------------------------------
    Int_t    particle_pid;
    Double_t traj_x_6,  traj_y_6;
    Double_t traj_x_18, traj_y_18;
    Double_t traj_x_36, traj_y_36;
    Double_t traj_edge_6, traj_edge_18, traj_edge_36;

    chain_sidisdvcs.SetBranchAddress("particle_pid", &particle_pid);
    chain_sidisdvcs.SetBranchAddress("traj_x_6",     &traj_x_6);
    chain_sidisdvcs.SetBranchAddress("traj_y_6",     &traj_y_6);
    chain_sidisdvcs.SetBranchAddress("traj_x_18",    &traj_x_18);
    chain_sidisdvcs.SetBranchAddress("traj_y_18",    &traj_y_18);
    chain_sidisdvcs.SetBranchAddress("traj_x_36",    &traj_x_36);
    chain_sidisdvcs.SetBranchAddress("traj_y_36",    &traj_y_36);
    chain_sidisdvcs.SetBranchAddress("traj_edge_6",  &traj_edge_6);
    chain_sidisdvcs.SetBranchAddress("traj_edge_18", &traj_edge_18);
    chain_sidisdvcs.SetBranchAddress("traj_edge_36", &traj_edge_36);

    chain_clasdis.SetBranchAddress("particle_pid", &particle_pid);
    chain_clasdis.SetBranchAddress("traj_x_6",     &traj_x_6);
    chain_clasdis.SetBranchAddress("traj_y_6",     &traj_y_6);
    chain_clasdis.SetBranchAddress("traj_x_18",    &traj_x_18);
    chain_clasdis.SetBranchAddress("traj_y_18",    &traj_y_18);
    chain_clasdis.SetBranchAddress("traj_x_36",    &traj_x_36);
    chain_clasdis.SetBranchAddress("traj_y_36",    &traj_y_36);
    chain_clasdis.SetBranchAddress("traj_edge_6",  &traj_edge_6);
    chain_clasdis.SetBranchAddress("traj_edge_18", &traj_edge_18);
    chain_clasdis.SetBranchAddress("traj_edge_36", &traj_edge_36);

    //--------------------------------------------------------------------------
    // 4. Book all histograms
    //--------------------------------------------------------------------------
    // 4a) 2D uncut
    int nbins2D = 200;
    TH2D* h2_r1_sidisdvcs = new TH2D("h2_r1_sidisdvcs", "Region 1; x; y",
                                     nbins2D, -180.0, 180.0, nbins2D, -180.0, 180.0);
    TH2D* h2_r2_sidisdvcs = new TH2D("h2_r2_sidisdvcs", "Region 2; x; y",
                                     nbins2D, -280.0, 280.0, nbins2D, -280.0, 280.0);
    TH2D* h2_r3_sidisdvcs = new TH2D("h2_r3_sidisdvcs", "Region 3; x; y",
                                     nbins2D, -450.0, 450.0, nbins2D, -450.0, 450.0);

    TH2D* h2_r1_clasdis = new TH2D("h2_r1_clasdis", "Region 1; x; y",
                                   nbins2D, -180.0, 180.0, nbins2D, -180.0, 180.0);
    TH2D* h2_r2_clasdis = new TH2D("h2_r2_clasdis", "Region 2; x; y",
                                   nbins2D, -280.0, 280.0, nbins2D, -280.0, 280.0);
    TH2D* h2_r3_clasdis = new TH2D("h2_r3_clasdis", "Region 3; x; y",
                                   nbins2D, -450.0, 450.0, nbins2D, -450.0, 450.0);

    // 4b) 2D cut
    TH2D* h2_r1_sidisdvcs_cut = new TH2D("h2_r1_sidisdvcs_cut", "Region 1 (cut); x; y",
                                         nbins2D, -180.0, 180.0, nbins2D, -180.0, 180.0);
    TH2D* h2_r2_sidisdvcs_cut = new TH2D("h2_r2_sidisdvcs_cut", "Region 2 (cut); x; y",
                                         nbins2D, -280.0, 280.0, nbins2D, -280.0, 280.0);
    TH2D* h2_r3_sidisdvcs_cut = new TH2D("h2_r3_sidisdvcs_cut", "Region 3 (cut); x; y",
                                         nbins2D, -450.0, 450.0, nbins2D, -450.0, 450.0);

    TH2D* h2_r1_clasdis_cut = new TH2D("h2_r1_clasdis_cut", "Region 1 (cut); x; y",
                                       nbins2D, -180.0, 180.0, nbins2D, -180.0, 180.0);
    TH2D* h2_r2_clasdis_cut = new TH2D("h2_r2_clasdis_cut", "Region 2 (cut); x; y",
                                       nbins2D, -280.0, 280.0, nbins2D, -280.0, 280.0);
    TH2D* h2_r3_clasdis_cut = new TH2D("h2_r3_clasdis_cut", "Region 3 (cut); x; y",
                                       nbins2D, -450.0, 450.0, nbins2D, -450.0, 450.0);

    // 4c) 1D edge histograms
    int nbins1D = 125;
    double edgeMin = 0.0, edgeMax = 125.0;
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

    // 4d) 2D rad-phi histograms
    int nbinsR1 = 150, nbinsR2 = 230, nbinsR3 = 400, nbinsPhi = 180;
    double rMin1 = 20.0, rMax1 = 170.0;
    double rMin2 = 20.0, rMax2 = 250.0;
    double rMin3 = 20.0, rMax3 = 420.0;
    double phiMin = 0.0, phiMax = 360.0;

    TH2D* h2_r1_sidis_radphi = new TH2D("h2_r1_sidis_radphi",
        "Region 1 rad-phi; radius (cm); phi (deg)",
        nbinsR1, rMin1, rMax1, nbinsPhi, phiMin, phiMax);
    TH2D* h2_r2_sidis_radphi = new TH2D("h2_r2_sidis_radphi",
        "Region 2 rad-phi; radius (cm); phi (deg)",
        nbinsR2, rMin2, rMax2, nbinsPhi, phiMin, phiMax);
    TH2D* h2_r3_sidis_radphi = new TH2D("h2_r3_sidis_radphi",
        "Region 3 rad-phi; radius (cm); phi (deg)",
        nbinsR3, rMin3, rMax3, nbinsPhi, phiMin, phiMax);

    TH2D* h2_r1_clas_radphi = new TH2D("h2_r1_clas_radphi",
        "Region 1 rad-phi; radius (cm); phi (deg)",
        nbinsR1, rMin1, rMax1, nbinsPhi, phiMin, phiMax);
    TH2D* h2_r2_clas_radphi = new TH2D("h2_r2_clas_radphi",
        "Region 2 rad-phi; radius (cm); phi (deg)",
        nbinsR2, rMin2, rMax2, nbinsPhi, phiMin, phiMax);
    TH2D* h2_r3_clas_radphi = new TH2D("h2_r3_clas_radphi",
        "Region 3 rad-phi; radius (cm); phi (deg)",
        nbinsR3, rMin3, rMax3, nbinsPhi, phiMin, phiMax);

    // Radii for the three circular cuts
    const double r1 = 140.0;
    const double r2 = 215.0;
    const double r3 = 290.0;

    gStyle->SetOptStat(0);

    //--------------------------------------------------------------------------
    // 5. Determine common event limit for both chains
    //--------------------------------------------------------------------------
    Long64_t nEntriesSIDIS = chain_sidisdvcs.GetEntries();
    Long64_t nEntriesCLAS  = chain_clasdis.GetEntries();
    Long64_t nLimit = nEntriesSIDIS;
    if (maxEventsSIDIS > 0 && maxEventsSIDIS < nLimit) {
        nLimit = maxEventsSIDIS;
    } // #endif
    if (nEntriesCLAS < nLimit) {
        nLimit = nEntriesCLAS;
    } // #endif

    std::cout << "Processing " << nLimit
              << " events from each chain (SIDIS-DVCS & CLASDIS)\n";

    //--------------------------------------------------------------------------
    // 6. Fill SIDIS-DVCS histograms
    //--------------------------------------------------------------------------
    for (Long64_t i = 0; i < nLimit; ++i) {
        chain_sidisdvcs.GetEntry(i);

        // // Filter line (commented out so we can revert):
        if (particle_pid != 11) continue; // electron
        // if (particle_pid != 2212) continue; // proton

        if (traj_edge_6 <= 3 || traj_edge_18 <= 3 || traj_edge_36 <= 10) continue; // #endif
        // if (traj_edge_6 <= 3 || traj_edge_18 <= 3 || traj_edge_36 <= 10) continue; // #endif

        // Region 1
        if (traj_x_6 != -9999) {
            h2_r1_sidisdvcs->Fill(traj_x_6, traj_y_6);
            double dist1 = std::sqrt(traj_x_6*traj_x_6 + traj_y_6*traj_y_6);
            if (dist1 < r1) {
                h2_r1_sidisdvcs_cut->Fill(traj_x_6, traj_y_6);
            }
            double rawAngle = std::atan2(traj_y_6, traj_x_6);
            double phi = rawAngle - M_PI;
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180.0 / M_PI;
            h2_r1_sidis_radphi->Fill(dist1, phiDeg);
        }

        // Region 2
        if (traj_x_18 != -9999) {
            h2_r2_sidisdvcs->Fill(traj_x_18, traj_y_18);
            double dist2 = std::sqrt(traj_x_18*traj_x_18 + traj_y_18*traj_y_18);
            if (dist2 < r2) {
                h2_r2_sidisdvcs_cut->Fill(traj_x_18, traj_y_18);
            }
            double rawAngle = std::atan2(traj_y_18, traj_x_18);
            double phi = rawAngle - M_PI;
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180.0 / M_PI;
            h2_r2_sidis_radphi->Fill(dist2, phiDeg);
        }

        // Region 3
        if (traj_x_36 != -9999) {
            h2_r3_sidisdvcs->Fill(traj_x_36, traj_y_36);
            double dist3 = std::sqrt(traj_x_36*traj_x_36 + traj_y_36*traj_y_36);
            if (dist3 < r3) {
                h2_r3_sidisdvcs_cut->Fill(traj_x_36, traj_y_36);
            }
            double rawAngle = std::atan2(traj_y_36, traj_x_36);
            double phi = rawAngle - M_PI;
            if (phi < 0) phi += 2.*M_PI;
            double phiDeg = phi * 180.0 / M_PI;
            h2_r3_sidis_radphi->Fill(dist3, phiDeg);
        }

        // Fill 1D edge histograms
        h1_edge_r1_sidis->Fill(traj_edge_6);
        h1_edge_r2_sidis->Fill(traj_edge_18);
        h1_edge_r3_sidis->Fill(traj_edge_36);
    }

    //--------------------------------------------------------------------------
    // 7. Fill CLASDIS histograms
    //--------------------------------------------------------------------------
    for (Long64_t i = 0; i < nLimit; ++i) {
        chain_clasdis.GetEntry(i);

        // // Filter line (commented out so we can revert):
        if (particle_pid != 11) continue; // electron
        // if (particle_pid != 2212) continue; // proton

        if (traj_edge_6 <= 3 || traj_edge_18 <= 3 || traj_edge_36 <= 10) continue; // #endif
        // if (traj_edge_6 <= 3 || traj_edge_18 <= 3 || traj_edge_36 <= 10) continue; // #endif

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
            double phiDeg = phi * 180.0 / M_PI;
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
            double phiDeg = phi * 180.0 / M_PI;
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
            double phiDeg = phi * 180.0 / M_PI;
            h2_r3_clas_radphi->Fill(dist3, phiDeg);
        }

        // Fill 1D edge histograms
        h1_edge_r1_clas->Fill(traj_edge_6);
        h1_edge_r2_clas->Fill(traj_edge_18);
        h1_edge_r3_clas->Fill(traj_edge_36);
    }

    //--------------------------------------------------------------------------
    // 8. Normalize, compute ratios, and draw canvases
    //--------------------------------------------------------------------------
    // 8a) Normalize uncut 2D
    double I_r1_sidis = h2_r1_sidisdvcs->Integral();
    if (I_r1_sidis > 0) h2_r1_sidisdvcs->Scale(1.0/I_r1_sidis);
    double I_r2_sidis = h2_r2_sidisdvcs->Integral();
    if (I_r2_sidis > 0) h2_r2_sidisdvcs->Scale(1.0/I_r2_sidis);
    double I_r3_sidis = h2_r3_sidisdvcs->Integral();
    if (I_r3_sidis > 0) h2_r3_sidisdvcs->Scale(1.0/I_r3_sidis);

    double I_r1_clas = h2_r1_clasdis->Integral();
    if (I_r1_clas > 0) h2_r1_clasdis->Scale(1.0/I_r1_clas);
    double I_r2_clas = h2_r2_clasdis->Integral();
    if (I_r2_clas > 0) h2_r2_clasdis->Scale(1.0/I_r2_clas);
    double I_r3_clas = h2_r3_clasdis->Integral();
    if (I_r3_clas > 0) h2_r3_clasdis->Scale(1.0/I_r3_clas);

    // 8b) Clone & ratio uncut
    TH2D* h2_r1_ratio = (TH2D*)h2_r1_sidisdvcs->Clone("h2_r1_ratio");
    h2_r1_ratio->SetTitle("Region 1");
    h2_r1_ratio->Divide(h2_r1_clasdis);

    TH2D* h2_r2_ratio = (TH2D*)h2_r2_sidisdvcs->Clone("h2_r2_ratio");
    h2_r2_ratio->SetTitle("Region 2");
    h2_r2_ratio->Divide(h2_r2_clasdis);

    TH2D* h2_r3_ratio = (TH2D*)h2_r3_sidisdvcs->Clone("h2_r3_ratio");
    h2_r3_ratio->SetTitle("Region 3");
    h2_r3_ratio->Divide(h2_r3_clasdis);

    // 8c) Unify scale & draw uncut (Canvas 1)
    SetSame2DScale(h2_r1_ratio, h2_r2_ratio, h2_r3_ratio);

    TCanvas* c2D_uncut = new TCanvas("c2D_uncut", "Normalized 2D Drift Chambers (Uncut)", 1800, 600);
    c2D_uncut->Divide(3,1);
    for (int i = 0; i < 3; ++i) {
        c2D_uncut->cd(i+1);
        gPad->SetLeftMargin(0.20);
        gPad->SetRightMargin(0.20);
        gPad->SetLogz();
        (i==0? h2_r1_ratio : i==1? h2_r2_ratio : h2_r3_ratio)->Draw("COLZ");
    }
    c2D_uncut->SaveAs("output/normalization/normalized_drift_chambers_uncut.png");

    // 8d) Draw with circles (Canvas 2)
    TCanvas* c2D_circle = new TCanvas("c2D_circle", "Normalized 2D Drift Chambers (Circle Overlay)", 1800, 600);
    c2D_circle->Divide(3,1);
    int circleColor = kRed, circleLW = 2, circleFS = 0;
    for (int i = 0; i < 3; ++i) {
        c2D_circle->cd(i+1);
        gPad->SetLeftMargin(0.20);
        gPad->SetRightMargin(0.20);
        gPad->SetLogz();
        TH2D* hist = (i==0? h2_r1_ratio : i==1? h2_r2_ratio : h2_r3_ratio);
        double R = (i==0? r1 : i==1? r2 : r3);
        hist->Draw("COLZ");
        TEllipse* cir = new TEllipse(0,0,R,R);
        cir->SetLineColor(circleColor);
        cir->SetLineWidth(circleLW);
        cir->SetFillStyle(circleFS);
        cir->Draw("same");
    }
    c2D_circle->SaveAs("output/normalization/normalized_drift_chambers_circle.png");

    //--------------------------------------------------------------------------
    // 9. Normalize & ratio the "cut" histograms
    //--------------------------------------------------------------------------
    double I_r1_sidis_cut = h2_r1_sidisdvcs_cut->Integral();
    if (I_r1_sidis_cut > 0) h2_r1_sidisdvcs_cut->Scale(1.0/I_r1_sidis_cut);
    double I_r2_sidis_cut = h2_r2_sidisdvcs_cut->Integral();
    if (I_r2_sidis_cut > 0) h2_r2_sidisdvcs_cut->Scale(1.0/I_r2_sidis_cut);
    double I_r3_sidis_cut = h2_r3_sidisdvcs_cut->Integral();
    if (I_r3_sidis_cut > 0) h2_r3_sidisdvcs_cut->Scale(1.0/I_r3_sidis_cut);

    double I_r1_clas_cut = h2_r1_clasdis_cut->Integral();
    if (I_r1_clas_cut > 0) h2_r1_clasdis_cut->Scale(1.0/I_r1_clas_cut);
    double I_r2_clas_cut = h2_r2_clasdis_cut->Integral();
    if (I_r2_clas_cut > 0) h2_r2_clasdis_cut->Scale(1.0/I_r2_clas_cut);
    double I_r3_clas_cut = h2_r3_clasdis_cut->Integral();
    if (I_r3_clas_cut > 0) h2_r3_clasdis_cut->Scale(1.0/I_r3_clas_cut);

    TH2D* h2_r1_ratio_cut = (TH2D*)h2_r1_sidisdvcs_cut->Clone("h2_r1_ratio_cut");
    h2_r1_ratio_cut->SetTitle("Region 1 (cut)");
    h2_r1_ratio_cut->Divide(h2_r1_clasdis_cut);

    TH2D* h2_r2_ratio_cut = (TH2D*)h2_r2_sidisdvcs_cut->Clone("h2_r2_ratio_cut");
    h2_r2_ratio_cut->SetTitle("Region 2 (cut)");
    h2_r2_ratio_cut->Divide(h2_r2_clasdis_cut);

    TH2D* h2_r3_ratio_cut = (TH2D*)h2_r3_sidisdvcs_cut->Clone("h2_r3_ratio_cut");
    h2_r3_ratio_cut->SetTitle("Region 3 (cut)");
    h2_r3_ratio_cut->Divide(h2_r3_clasdis_cut);

    SetSame2DScale(h2_r1_ratio_cut, h2_r2_ratio_cut, h2_r3_ratio_cut);

    TCanvas* c2D_cut = new TCanvas("c2D_cut", "Normalized 2D Drift Chambers (Cut)", 1800, 600);
    c2D_cut->Divide(3,1);
    for (int i = 0; i < 3; ++i) {
        c2D_cut->cd(i+1);
        gPad->SetLeftMargin(0.20);
        gPad->SetRightMargin(0.20);
        gPad->SetLogz();
        (i==0? h2_r1_ratio_cut : i==1? h2_r2_ratio_cut : h2_r3_ratio_cut)->Draw("COLZ");
    }
    c2D_cut->SaveAs("output/normalization/normalized_drift_chambers_cut.png");

    //--------------------------------------------------------------------------
    // 10. 1D edge histograms: normalize & ratio
    //--------------------------------------------------------------------------
    double E1_sidis = h1_edge_r1_sidis->Integral();
    if (E1_sidis > 0) h1_edge_r1_sidis->Scale(1.0/E1_sidis);
    double E2_sidis = h1_edge_r2_sidis->Integral();
    if (E2_sidis > 0) h1_edge_r2_sidis->Scale(1.0/E2_sidis);
    double E3_sidis = h1_edge_r3_sidis->Integral();
    if (E3_sidis > 0) h1_edge_r3_sidis->Scale(1.0/E3_sidis);

    double E1_clas = h1_edge_r1_clas->Integral();
    if (E1_clas > 0) h1_edge_r1_clas->Scale(1.0/E1_clas);
    double E2_clas = h1_edge_r2_clas->Integral();
    if (E2_clas > 0) h1_edge_r2_clas->Scale(1.0/E2_clas);
    double E3_clas = h1_edge_r3_clas->Integral();
    if (E3_clas > 0) h1_edge_r3_clas->Scale(1.0/E3_clas);

    TH1D* h1_edge_r1_ratio = (TH1D*)h1_edge_r1_sidis->Clone("h1_edge_r1_ratio");
    h1_edge_r1_ratio->SetTitle("Region 1");
    h1_edge_r1_ratio->Divide(h1_edge_r1_clas);

    TH1D* h1_edge_r2_ratio = (TH1D*)h1_edge_r2_sidis->Clone("h1_edge_r2_ratio");
    h1_edge_r2_ratio->SetTitle("Region 2");
    h1_edge_r2_ratio->Divide(h1_edge_r2_clas);

    TH1D* h1_edge_r3_ratio = (TH1D*)h1_edge_r3_sidis->Clone("h1_edge_r3_ratio");
    h1_edge_r3_ratio->SetTitle("Region 3");
    h1_edge_r3_ratio->Divide(h1_edge_r3_clas);

    SetSame1DScale(h1_edge_r1_ratio, h1_edge_r2_ratio, h1_edge_r3_ratio);

    TCanvas* c1D = new TCanvas("c1D", "Edge Ratios", 1800, 600);
    c1D->Divide(3,1);
    // Region 1
    c1D->cd(1);
    h1_edge_r1_ratio->GetXaxis()->SetTitle("edge (cm)");
    h1_edge_r1_ratio->GetYaxis()->SetTitle("normalized density ratio");
    h1_edge_r1_ratio->Draw("HIST");
    // Region 2
    c1D->cd(2);
    h1_edge_r2_ratio->GetXaxis()->SetTitle("edge (cm)");
    h1_edge_r2_ratio->GetYaxis()->SetTitle("normalized density ratio");
    h1_edge_r2_ratio->Draw("HIST");
    // Region 3
    c1D->cd(3);
    h1_edge_r3_ratio->GetXaxis()->SetTitle("edge (cm)");
    h1_edge_r3_ratio->GetYaxis()->SetTitle("normalized density ratio");
    h1_edge_r3_ratio->Draw("HIST");

    c1D->SaveAs("output/normalization/normalized_drift_chambers_edges.png");

    //--------------------------------------------------------------------------
    // 11. Final "rad-phi" ratio canvases
    //--------------------------------------------------------------------------
    double RP1_sidis = h2_r1_sidis_radphi->Integral();
    if (RP1_sidis > 0) h2_r1_sidis_radphi->Scale(1.0/RP1_sidis);
    double RP2_sidis = h2_r2_sidis_radphi->Integral();
    if (RP2_sidis > 0) h2_r2_sidis_radphi->Scale(1.0/RP2_sidis);
    double RP3_sidis = h2_r3_sidis_radphi->Integral();
    if (RP3_sidis > 0) h2_r3_sidis_radphi->Scale(1.0/RP3_sidis);

    double RP1_clas = h2_r1_clas_radphi->Integral();
    if (RP1_clas > 0) h2_r1_clas_radphi->Scale(1.0/RP1_clas);
    double RP2_clas = h2_r2_clas_radphi->Integral();
    if (RP2_clas > 0) h2_r2_clas_radphi->Scale(1.0/RP2_clas);
    double RP3_clas = h2_r3_clas_radphi->Integral();
    if (RP3_clas > 0) h2_r3_clas_radphi->Scale(1.0/RP3_clas);

    TH2D* h2_r1_ratio_radphi = (TH2D*)h2_r1_sidis_radphi->Clone("h2_r1_ratio_radphi");
    h2_r1_ratio_radphi->SetTitle("Region 1 (radius vs phi)");
    h2_r1_ratio_radphi->Divide(h2_r1_clas_radphi);

    TH2D* h2_r2_ratio_radphi = (TH2D*)h2_r2_sidis_radphi->Clone("h2_r2_ratio_radphi");
    h2_r2_ratio_radphi->SetTitle("Region 2 (radius vs phi)");
    h2_r2_ratio_radphi->Divide(h2_r2_clas_radphi);

    TH2D* h2_r3_ratio_radphi = (TH2D*)h2_r3_sidis_radphi->Clone("h2_r3_ratio_radphi");
    h2_r3_ratio_radphi->SetTitle("Region 3 (radius vs phi)");
    h2_r3_ratio_radphi->Divide(h2_r3_clas_radphi);

    SetSame2DScale(h2_r1_ratio_radphi, h2_r2_ratio_radphi, h2_r3_ratio_radphi);

    TCanvas* cRadPhi = new TCanvas("cRadPhi", "Ratio of radius vs phi", 1800, 600);
    cRadPhi->Divide(3,1);
    for (int i = 0; i < 3; ++i) {
        cRadPhi->cd(i+1);
        gPad->SetLeftMargin(0.20);
        gPad->SetRightMargin(0.20);
        gPad->SetLogz();
        (i==0? h2_r1_ratio_radphi : i==1? h2_r2_ratio_radphi : h2_r3_ratio_radphi)->Draw("COLZ");
    }
    cRadPhi->SaveAs("output/normalization/normalized_drift_chambers_radphi.png");

    //--------------------------------------------------------------------------
    // Cleanup
    //--------------------------------------------------------------------------
    delete c2D_uncut;
    delete c2D_circle;
    delete c2D_cut;
    delete c1D;
    delete cRadPhi;

    // 2D hists (uncut)
    delete h2_r1_sidisdvcs; delete h2_r2_sidisdvcs; delete h2_r3_sidisdvcs;
    delete h2_r1_clasdis;   delete h2_r2_clasdis;   delete h2_r3_clasdis;
    delete h2_r1_ratio;     delete h2_r2_ratio;     delete h2_r3_ratio;

    // 2D hists (cut)
    delete h2_r1_sidisdvcs_cut; delete h2_r2_sidisdvcs_cut; delete h2_r3_sidisdvcs_cut;
    delete h2_r1_clasdis_cut;   delete h2_r2_clasdis_cut;   delete h2_r3_clasdis_cut;
    delete h2_r1_ratio_cut;     delete h2_r2_ratio_cut;     delete h2_r3_ratio_cut;

    // 1D hists
    delete h1_edge_r1_sidis; delete h1_edge_r2_sidis; delete h1_edge_r3_sidis;
    delete h1_edge_r1_clas;  delete h1_edge_r2_clas;  delete h1_edge_r3_clas;
    delete h1_edge_r1_ratio; delete h1_edge_r2_ratio; delete h1_edge_r3_ratio;

    // Rad-phi hists
    delete h2_r1_sidis_radphi; delete h2_r2_sidis_radphi; delete h2_r3_sidis_radphi;
    delete h2_r1_clas_radphi;  delete h2_r2_clas_radphi;  delete h2_r3_clas_radphi;
    delete h2_r1_ratio_radphi; delete h2_r2_ratio_radphi; delete h2_r3_ratio_radphi;

    return 0;
}