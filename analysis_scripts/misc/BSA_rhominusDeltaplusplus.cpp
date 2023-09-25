// Created 9/25/23
// get BSA of rho- Delta++ 

#include <TTree.h>
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


void createBSAPlot(TTreeReader &dataReader, const char* outDir) {
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
}

void BSA_rhominusDeltaplusplus(std::string root_file_path) {
    // Start the timer
    auto start_time = std::chrono::high_resolution_clock::now();
    gStyle->SetCanvasColor(0);

    TFile* file = new TFile(root_file_path.c_str(), "READ");

    if (!file->IsOpen()) {
        cout << "Error opening ROOT file (is the location correct?). Exiting." << endl;
    }

    TTree* tree = (TTree*)file->Get("PhysicsEvents");

    if (!tree) {
        cout << "Error getting trees from ROOT file." << endl;
    }

    TTreeReader dataReader(tree); // Create a TTreeReader for the data tree

    createBSAPlot(dataReader, "output");

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