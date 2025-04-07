// example.cpp
// This script demonstrates basic ROOT functionality by loading a ROOT file
// containing a TTree named "PhysicsEvents", printing its branches,
// filling histograms from branch data, and drawing the results on a canvas.
//
// To compile this script, use the following command (assuming ROOT is installed):
//   g++ example.cpp `root-config --cflags --libs` -o example
//
// To run the compiled executable (provide your ROOT file as argument):
//   ./example path/to/your/rootfile.root

// Include ROOT headers for file handling, tree processing, histogramming, and drawing.
#include <TFile.h>       // For opening ROOT files
#include <TTree.h>       // For handling TTree objects
#include <TBranch.h>     // For handling branches of a TTree
#include <TCanvas.h>     // For drawing canvases
#include <TH1F.h>        // For 1D histograms
#include <TH2F.h>        // For 2D histograms
#include <TApplication.h>// For interactive application display (if needed)
#include <TMath.h>       // For mathematical functions (e.g., Pi)
#include <iostream>      // For standard I/O operations
#include <TROOT.h>       // load gROOT to run in batch mode

using namespace std;

int main(int argc, char **argv) {

    gROOT->SetBatch(kTRUE); // Run ROOT in batch mode (no visual displays etc.)

    //----------------------------------------------------------------------------------
    // Check command-line arguments: the user must supply a ROOT file name.
    //----------------------------------------------------------------------------------
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <ROOT file>" << endl;
        return 1;
    }

    //----------------------------------------------------------------------------------
    // Open the ROOT file provided by the user.
    // TFile::Open returns a pointer to the file. If the file cannot be opened or is
    // invalid (IsZombie), we exit the program with an error message.
    //----------------------------------------------------------------------------------
    TFile *file = TFile::Open(argv[1]);
    if (!file || file->IsZombie()) {
        cout << "Error: Cannot open file " << argv[1] << endl;
        return 1;
    }

    //----------------------------------------------------------------------------------
    // Retrieve the TTree named "PhysicsEvents" from the file.
    // If the tree isn't found, print an error and exit.
    //----------------------------------------------------------------------------------
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        cout << "Error: TTree 'PhysicsEvents' not found in file " << argv[1] << endl;
        file->Close();
        return 1;
    }

    //----------------------------------------------------------------------------------
    // Print all branch names in the TTree.
    // This loops over all branches in the tree and outputs their names.
    //----------------------------------------------------------------------------------
    cout << "Branches in the TTree:" << endl;
    TObjArray *branchList = tree->GetListOfBranches();
    for (int i = 0; i < branchList->GetEntries(); i++) {
        TBranch *branch = (TBranch*)branchList->At(i);
        cout << "  " << branch->GetName() << endl;
    }

    //----------------------------------------------------------------------------------
    // Define variables to hold the data for the branches we are interested in.
    //----------------------------------------------------------------------------------
    double e_p = 0.0;      // Variable to hold energy (GeV)
    double e_theta = 0.0;  // Variable to hold theta (radians)

    //----------------------------------------------------------------------------------
    // Set branch addresses so that data from the TTree is loaded into our variables.
    //----------------------------------------------------------------------------------
    tree->SetBranchAddress("e_p", &e_p);
    tree->SetBranchAddress("e_theta", &e_theta);

    //----------------------------------------------------------------------------------
    // Create a 1D histogram for the "e_p" branch.
    // The histogram has 100 bins ranging from 0 to 10.
    // The x-axis label is set with LaTeX formatting: "#it{#e_{p}} (GeV)"
    //----------------------------------------------------------------------------------
    TH1F *hist_e_p = new TH1F("hist_e_p", "Histogram of e_{p};#it{e_{p}} (GeV);Entries", 100, 0, 10);

    //----------------------------------------------------------------------------------
    // Create a 2D histogram to plot "e_theta" (in degrees) versus "e_p".
    // e_theta is converted from radians to degrees.
    // The x-axis is for e_p (from 0 to 10 GeV), and the y-axis is for e_theta in degrees (0 to 180).
    //----------------------------------------------------------------------------------
    TH2F *hist2D = new TH2F("hist2D", "e_{theta} (deg) vs e_{p};#it{e_{p}} (GeV);#it{#theta} (degrees)", 100, 0, 10, 100, 0, 30);

    //----------------------------------------------------------------------------------
    // Loop over all entries in the TTree to fill the histograms.
    // For each entry, we:
    //  - Load the current data into our variables.
    //  - Fill the 1D histogram with e_p.
    //  - Convert e_theta from radians to degrees and fill the 2D histogram.
    //----------------------------------------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);                   // Load entry 'i'
        hist_e_p->Fill(e_p);                 // Fill the 1D histogram with e_p

        // Convert e_theta from radians to degrees.
        double e_theta_deg = e_theta * 180.0 / TMath::Pi();
        // Fill the 2D histogram: e_p on the x-axis, e_theta (in degrees) on the y-axis.
        hist2D->Fill(e_p, e_theta_deg);
    }

    //----------------------------------------------------------------------------------
    // Create a canvas to display the plots.
    // The canvas is divided into 1 row and 2 columns.
    //----------------------------------------------------------------------------------
    TCanvas *canvas = new TCanvas("canvas", "ROOT Example", 1200, 600);
    canvas->Divide(2, 1);

    //----------------------------------------------------------------------------------
    // Draw the 1D histogram (e_p) on the left pad.
    //----------------------------------------------------------------------------------
    canvas->cd(1);                        // Switch to pad 1
    hist_e_p->SetLineColor(kBlue);        // Set line color to blue
    hist_e_p->Draw();                     // Draw the histogram

    //----------------------------------------------------------------------------------
    // Draw the 2D histogram (e_theta vs. e_p) on the right pad.
    // The "COLZ" option displays the plot with a color palette.
    //----------------------------------------------------------------------------------
    canvas->cd(2);                        // Switch to pad 2
    hist2D->Draw("COLZ");                 // Draw the 2D histogram with color mapping

    //----------------------------------------------------------------------------------
    // Save the canvas as a PNG file in the current working directory.
    //----------------------------------------------------------------------------------
    canvas->SaveAs("example.png");

    //----------------------------------------------------------------------------------
    // Optional: If you wish to interact with the canvas (zoom, move, etc.), you can run
    // the TApplication event loop. If you prefer batch mode (i.e., not displaying a window),
    // you can comment out the following lines.
    //----------------------------------------------------------------------------------
    // TApplication app("app", &argc, argv);
    // canvas->Update();  // Update the canvas to ensure it is drawn
    // app.Run();         // Start the application event loop

    //----------------------------------------------------------------------------------
    // Clean up by closing the file
    //----------------------------------------------------------------------------------
    file->Close();

    // End of program.
    return 0;
}