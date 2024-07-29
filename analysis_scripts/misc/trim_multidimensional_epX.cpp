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
        output_trees[i] = input_tree->CloneTree(0); // Clone the structure of the input tree
    }

    // Create a TTreeReader to read the input tree
    TTreeReader reader(input_tree);
    TTreeReaderValue<int> runnum(reader, "runnum");
    TTreeReaderValue<int> evnum(reader, "evnum");
    TTreeReaderValue<int> helicity(reader, "helicity");
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
    TTreeReaderValue<double> Q2(reader, "Q2");
    TTreeReaderValue<double> W(reader, "W");
    TTreeReaderValue<double> Mx(reader, "Mx");
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
        int random_int = *runnum + *evnum + *helicity;
        double random = *beam_pol + *target_pol + *e_p + *e_theta + *e_phi + 
            *vz_e + *p_p + *p_theta + *p_phi + *vz_p + *Q2 + *W + *Mx + 
            *Mx2 + *x + *y + *t + *tmin + *z + *xF + *pT + *zeta + *eta + *phi + *DepA +
            *DepB + *DepC + *DepV + *DepW;
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