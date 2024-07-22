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

    // Create a TChain to use TTreeReader
    TChain chain("PhysicsEvents");
    chain.Add(input_filename);

    TTreeReader reader(&chain);
    TTreeReaderValue<double> Mx(reader, "Mx");
    TTreeReaderValue<double> Q2(reader, "Q2");
    TTreeReaderValue<double> y(reader, "y");

    for (int i = 0; i < 18; ++i) {
        std::string output_filename = base_name + "_" + std::to_string(i) + ".root";
        TFile *output_file = TFile::Open(output_filename.c_str(), "RECREATE");
        TTree *output_tree = input_tree->CloneTree(0);

        while (reader.Next()) {
            if (*Mx > 1.4) {
                if (i == 0 || DetermineQ2yBin(*Q2, *y) == i) {
                    output_tree->Fill();
                }
            }
        }

        output_file->Write();
        output_file->Close();
    }

    input_file->Close();
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