#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iostream>
#include <string>


int DetermineQ2xBin(float Q2, float x) {
    // Bin 1: x < 0.1
    if (x < 0.1) {
        return 1;
    }

    // Bins 2-4: 0.1 < x < 0.14
    if (x > 0.1 && x < 0.14) {
        if (Q2 < 1.50) return 2;
        if (Q2 >= 1.50 && Q2 < 1.70) return 3;
        if (Q2 >= 1.70) return 4;
    }

    // Bins 5-8: 0.14 < x < 0.21
    if (x > 0.14 && x < 0.21) {
        if (Q2 < 1.50) return 5;
        if (Q2 >= 1.50 && Q2 < 1.70) return 6;
        if (Q2 >= 1.70 && Q2 < 2.00) return 7;
        if (Q2 >= 2.00) return 8;
    }

    // Bins 9-12: 0.21 < x < 0.30
    if (x > 0.21 && x < 0.30) {
        if (Q2 < 2.20) return 9;
        if (Q2 >= 2.20 && Q2 < 2.60) return 10;
        if (Q2 >= 2.60) return 11;
    }

    // Bins 13-14: 0.30 < x < 0.42
    if (x > 0.30 && x < 0.42) {
        if (Q2 < 3.20) return 12;
        if (Q2 >= 3.20) return 13;
    }

    // Bin 15: x > 0.42
    if (x >= 0.42) {
        return 14;
    }

    // If no conditions are met, return 0
    return -2;
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
    for (int i = 0; i < 15; ++i) {
        std::string output_filename = base_name + "_" + std::to_string(i) + ".root";
        output_files[i] = TFile::Open(output_filename.c_str(), "RECREATE");
        output_trees[i] = input_tree->CloneTree(0); // Clone the structure of the input tree
    }

    // Create a TTreeReader to read the input tree
    TTreeReader reader(input_tree);
    TTreeReaderValue<int> fiducial_status(reader, "fiducial_status");
    TTreeReaderValue<int> num_pos(reader, "num_pos");
    TTreeReaderValue<int> num_neg(reader, "num_neg");
    TTreeReaderValue<int> num_neutrals(reader, "num_neutrals");
    TTreeReaderValue<int> runnum(reader, "runnum");
    TTreeReaderValue<int> evnum(reader, "evnum");
    TTreeReaderValue<int> helicity(reader, "helicity");
    TTreeReaderValue<int> detector(reader, "detector");
    TTreeReaderValue<double> beam_pol(reader, "beam_pol");
    TTreeReaderValue<double> target_pol(reader, "target_pol");
    TTreeReaderValue<double> e_p(reader, "e_p");
    TTreeReaderValue<double> e_theta(reader, "e_theta");
    TTreeReaderValue<double> e_phi(reader, "e_phi");
    TTreeReaderValue<double> vz_e(reader, "vz_e");
    TTreeReaderValue<double> p_p(reader, "p_p");
    TTreeReaderValue<double> p_theta(reader, "p_theta");
    TTreeReaderValue<double> p_phi(reader, "p_phi");
    TTreeReaderValue<double> vz_p(reader, "vz_p");
    TTreeReaderValue<double> open_angle(reader, "open_angle");
    TTreeReaderValue<double> Q2(reader, "Q2");
    TTreeReaderValue<double> W(reader, "W");
    TTreeReaderValue<double> Mx2(reader, "Mx2");
    TTreeReaderValue<double> x(reader, "x");
    TTreeReaderValue<double> y(reader, "y");
    TTreeReaderValue<double> t(reader, "t");
    TTreeReaderValue<double> tmin(reader, "tmin");
    TTreeReaderValue<double> z(reader, "z");
    TTreeReaderValue<double> xF(reader, "xF");
    TTreeReaderValue<double> pT(reader, "pT");
    TTreeReaderValue<double> zeta(reader, "zeta");
    TTreeReaderValue<double> eta(reader, "eta");
    TTreeReaderValue<double> phi(reader, "phi");
    TTreeReaderValue<double> DepA(reader, "DepA");
    TTreeReaderValue<double> DepB(reader, "DepB");
    TTreeReaderValue<double> DepC(reader, "DepC");
    TTreeReaderValue<double> DepV(reader, "DepV");
    TTreeReaderValue<double> DepW(reader, "DepW");

    // Loop over entries and fill the corresponding output trees
    while (reader.Next()) {
        int random_int = *fiducial_status + *num_pos + *num_neg + *num_neutrals + *runnum + 
            *evnum + *helicity + *detector;
        double random = *beam_pol + *target_pol + *e_p + *e_theta + *e_phi + 
            *vz_e + *p_p + *p_theta + *p_phi + *vz_p + *open_angle + *Q2 + *W  + 
            *Mx2 + *x + *y + *t + *tmin + *z + *xF + *pT + *zeta + *eta + *phi + *DepA +
            *DepB + *DepC + *DepV + *DepW;
        if (*Mx2 < 0.00) {
            continue;
        }
        if (*Mx2 > 1.7) {
            // Determine the Q2-y bin and fill the corresponding tree
            int bin = DetermineQ2xBin(*Q2, *x);
            if (bin > 0 && bin <= 14) {
                output_trees[bin]->Fill();
            }
        }t
        // Also fill the general tree for one dimensional dependencies
        output_trees[0]->Fill();
    }

    // Write and close the output files
    for (int i = 0; i < 15; ++i) {
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