// Created 9/6/23
// Quick script to add branches to trihadron tree for the missing particles

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>


void addMissingParticleVariables(std::string root_file_path) {
    // Open ROOT file and get the tree
    TFile* file = new TFile(root_file_path.c_str(), "UPDATE");
    TTree* tree = (TTree*)file->Get("PhysicsEvents");

    // Declare variables to read from the tree
    double e_p, e_theta, e_phi, p1_p, p1_theta, p1_phi, p2_p, p2_theta, p2_phi, p3_p, p3_theta; 
    double p3_phi;

    // Set branch addresses
    tree->SetBranchAddress("e_p", &e_p);
    tree->SetBranchAddress("e_theta", &e_theta);
    tree->SetBranchAddress("e_phi", &e_phi);
    tree->SetBranchAddress("p1_p", &p1_p);
    tree->SetBranchAddress("p1_theta", &p1_theta);
    tree->SetBranchAddress("p1_phi", &p1_phi);
    tree->SetBranchAddress("p2_p", &p2_p);
    tree->SetBranchAddress("p2_theta", &p2_theta);
    tree->SetBranchAddress("p2_phi", &p2_phi);
    tree->SetBranchAddress("p3_p", &p3_p);
    tree->SetBranchAddress("p3_theta", &p3_theta);
    tree->SetBranchAddress("p3_phi", &p3_phi);

    // Declare new variables to store missing particle information
    float px_p, px_theta, px_phi, Mx1x, Mx2x, Mx3x, Mh1x, Mh2x, Mh3x;

    // Create new branches
    TTree* newTree = tree->CloneTree(0); // clone only structure, not events
    newTree->Branch("px_p", &px_p, "px_p/D");
    newTree->Branch("px_theta", &px_theta, "px_theta/D");
    newTree->Branch("px_phi", &px_phi, "px_phi/D");
    newTree->Branch("Mx1x", &Mx1x, "Mx1x/D");
    newTree->Branch("Mx2x", &Mx2x, "Mx2x/D");
    newTree->Branch("Mx3x", &Mx3x, "Mx3x/D");
    newTree->Branch("Mh1x", &Mh1x, "Mh1x/D");
    newTree->Branch("Mh2x", &Mh2x, "Mh2x/D");
    newTree->Branch("Mh3x", &Mh3x, "Mh3x/D");

    // Define initial state 4-momentum (10.1998 GeV electron beam and stationary proton)
    TLorentzVector p_initial(0, 0, 10.1998, 10.1998 + 0.938); // (px, py, pz, E)

    cout << "passed" << endl;
    // Loop over all events in the tree
    Long64_t nEntries = tree->GetEntries();
    for(Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Create 4-momentum vectors for final state particles
        TLorentzVector p_e, p1, p2, p3;
        p_e.SetXYZM(e_p * sin(e_theta) * cos(e_phi), e_p * sin(e_theta) * sin(e_phi), e_p * cos(e_theta), 0.511e-3);
        p1.SetXYZM(p1_p * sin(p1_theta) * cos(p1_phi), p1_p * sin(p1_theta) * sin(p1_phi), p1_p * cos(p1_theta), 0.13957);
        p2.SetXYZM(p2_p * sin(p2_theta) * cos(p2_phi), p2_p * sin(p2_theta) * sin(p2_phi), p2_p * cos(p2_theta), 0.13957);
        p3.SetXYZM(p3_p * sin(p3_theta) * cos(p3_phi), p3_p * sin(p3_theta) * sin(p3_phi), p3_p * cos(p3_theta), 0.938);

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

        // Fill new tree
        newTree->Fill();
    }

    // Extract the file name, remove the ".root" extension and append "_MissingParticles.root"
    std::string new_file_name = 
        root_file_path.substr(0, root_file_path.find(".root")) + "_MissingParticles.root";

    // Create a new ROOT file to store the new tree
    TFile* newFile = new TFile(new_file_name.c_str(), "RECREATE");

    // Write new tree to the new file
    newTree->Write();

    // Close all the files
    newFile->Close();
    file->Close();

    // Memory cleanup
    delete newTree;
    delete tree;
    delete newFile;
    delete file;
}
