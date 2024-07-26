#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iostream>
#include <string>

int DetermineQ2yBin(float Q2, float y) {
    // Q2-y Bins 1-4
    if (Q2 > 1.000 && Q2 <= 2.000) {
        if (y > 0.650 && y <= 0.750) return 1;
        if (y > 0.550 && y <= 0.650) return 2;
        if (y > 0.450 && y <= 0.550) return 3;
        if (y > 0.300 && y <= 0.450) return 4;
    }
    // Q2-y Bins 5-8
    else if (Q2 > 2.000 && Q2 <= 3.000) {
        if (y > 0.650 && y <= 0.750) return 5;
        if (y > 0.550 && y <= 0.650) return 6;
        if (y > 0.450 && y <= 0.550) return 7;
        if (y > 0.300 && y <= 0.450) return 8;
    }
    // Q2-y Bins 9-12
    else if (Q2 > 3.000 && Q2 <= 4.000) {
        if (y > 0.650 && y <= 0.750) return 9;
        if (y > 0.550 && y <= 0.650) return 10;
        if (y > 0.450 && y <= 0.550) return 11;
        if (y > 0.300 && y <= 0.450) return 12;
    }
    // Q2-y Bins 13-15
    else if (Q2 > 4.000 && Q2 <= 5.000) {
        if (y > 0.650 && y <= 0.750) return 13;
        if (y > 0.550 && y <= 0.650) return 14;
        if (y > 0.450 && y <= 0.550) return 15;
    }
    // Q2-y Bins 16-17
    else if (Q2 > 5.000 && Q2 <= 7.000) {
        if (y > 0.650 && y <= 0.750) return 16;
        if (y > 0.550 && y <= 0.650) return 17;
    }
    return 0;
}

void process_file(const char* input_filename) {
    TFile *input_file = TFile::Open(input_filename);
    if (!input_file || input_file->IsZombie()) {
        std::cerr << "Error opening file: " << input_filename << std::endl;
        return;
    }

    TTree *input_tree;
    input_file->GetObject("PhysicsEvents", input_tree);
    if (!input_tree) {
        std::cerr << "Error: PhysicsEvents tree not found!" << std::endl;
        input_file->Close();
        return;
    }

    // Define the output filenames
    std::string base_name(input_filename);
    size_t dot_pos = base_name.find_last_of('.');
    if (dot_pos != std::string::npos) {
        base_name = base_name.substr(0, dot_pos);
    }

    // Create output files and trees for each Q2-y bin
    std::vector<TFile*> output_files(18);
    std::vector<TTree*> output_trees(18);
    for (int i = 0; i < 18; ++i) {
        std::string output_filename = base_name + "_" + std::to_string(i) + ".root";
        output_files[i] = TFile::Open(output_filename.c_str(), "RECREATE");
        output_trees[i] = input_tree->CloneTree(0); // Clone structure, no entries yet
    }

    // Create a TTreeReader to read the input tree
    TTreeReader reader(input_tree);
    tree->Branch("runnum", &runnum, "runnum/I");
    tree->Branch("evnum", &evnum, "evnum/I");
    tree->Branch("helicity", &helicity, "helicity/I");
    tree->Branch("beam_pol", &beam_pol, "beam_pol/D");
    tree->Branch("target_pol", &target_pol, "target_pol/D");
    tree->Branch("e_p", &e_p, "e_p/D");
    tree->Branch("e_theta", &e_theta, "e_theta/D");
    tree->Branch("e_phi", &e_phi, "e_phi/D");
    tree->Branch("vz_e", &vz_e, "vz_e/D");
    tree->Branch("p_p", &p_p, "p_p/D");
    tree->Branch("p_theta", &p_theta, "p_theta/D");
    tree->Branch("p_phi", &p_phi, "p_phi/D");
    tree->Branch("vz_p", &vz_p, "vz_p/D");
    tree->Branch("Q2", &Q2, "Q2/D");
    tree->Branch("W", &W, "W/D");
    tree->Branch("Mx", &Mx, "Mx/D");
    tree->Branch("Mx2", &Mx2, "Mx2/D");
    tree->Branch("x", &x, "x/D");
    tree->Branch("y", &y, "y/D");
    tree->Branch("t", &t, "t/D");
    tree->Branch("tmin", &tmin, "tmin/D");
    tree->Branch("z", &z, "z/D");
    tree->Branch("xF", &xF, "xF/D");
    tree->Branch("pT", &pT, "pT/D");
    tree->Branch("zeta", &zeta, "zeta/D");
    tree->Branch("eta", &eta, "eta/D");
    tree->Branch("phi", &phi, "phi/D");
    tree->Branch("DepA", &DepA, "DepA/D");
    tree->Branch("DepB", &DepB, "DepB/D");
    tree->Branch("DepC", &DepC, "DepC/D");
    tree->Branch("DepV", &DepV, "DepV/D");
    tree->Branch("DepW", &DepW, "DepW/D");

    // Loop over entries and fill the corresponding output trees
    while (reader.Next()) {
        if (*Mx > 1.4) {
            // Determine the Q2-y bin and fill the corresponding tree
            int bin = DetermineQ2yBin(*Q2, *y);
            if (bin > 0 && bin < 18) {
                output_trees[bin]->Fill();
            }
            // Also fill the general tree for events passing Mx cut
            output_trees[0]->Fill();
        }
    }

    // Write and close the output files
    for (int i = 0; i < 18; ++i) {
        output_files[i]->cd();
        output_trees[i]->Write();
        output_files[i]->Close();
        delete output_files[i];
    }

    input_file->Close();
    delete input_file;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_root_file>" << std::endl;
        return 1;
    }

    const char* input_filename = argv[1];
    process_file(input_filename);

    return 0;
}