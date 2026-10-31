#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TObject.h>
#include <TTree.h>

#include "reconstructed_data_schema.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool validate_written_file(const char* path, const std::vector<std::string>& names,
                           unsigned long long expected_rows) {
    std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot reopen output ROOT file\n";
        return false;
    }
    TTree* tree = nullptr;
    file->GetObject("PhysicsEvents", tree);
    if (!tree) {
        std::cerr << "Output is missing PhysicsEvents\n";
        return false;
    }
    if (tree->GetEntries() != static_cast<Long64_t>(expected_rows)) {
        std::cerr << "Output row count changed after write\n";
        return false;
    }
    auto* branches = tree->GetListOfBranches();
    if (!branches || branches->GetEntries() != static_cast<int>(names.size())) {
        std::cerr << "Output branch count does not match the selected schema\n";
        return false;
    }
    for (std::size_t i = 0; i < names.size(); ++i) {
        auto* branch = dynamic_cast<TBranch*>(branches->At(static_cast<int>(i)));
        if (!branch || names[i] != branch->GetName()) {
            std::cerr << "Output branch mismatch at position " << (i + 1) << "\n";
            return false;
        }
        TLeaf* leaf = branch->GetLeaf(names[i].c_str());
        if (!leaf || std::string(leaf->GetTypeName()) != "Double_t") {
            std::cerr << "Output branch " << names[i] << " is not a Double_t scalar\n";
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: convert_data_stream FIFO OUTPUT.root "
                     "dvcs_data|pi0_data|pi0_as_dvcs\n";
        return 2;
    }

    const std::string mode = argv[3];
    const auto* names = reconstructed_data_schema::for_mode(mode);
    if (!names) {
        std::cerr << "Unknown mode: " << mode << "\n";
        return 2;
    }
    if (!reconstructed_data_schema::definition_is_valid(mode, *names)) {
        std::cerr << "Internal schema definition is invalid for " << mode << "\n";
        return 5;
    }

    std::ifstream input(argv[1]);
    if (!input) {
        std::cerr << "Cannot open input stream\n";
        return 2;
    }

    bool bad = false;
    unsigned long long rows = 0;
    {
        TFile output(argv[2], "RECREATE");
        if (output.IsZombie()) return 3;

        TTree tree("PhysicsEvents", "Reprocessed CLAS12 physics events");
        std::vector<double> values(names->size());
        for (std::size_t i = 0; i < names->size(); ++i) {
            const std::string leaf = (*names)[i] + "/D";
            if (!tree.Branch((*names)[i].c_str(), &values[i], leaf.c_str())) {
                std::cerr << "Cannot create branch " << (*names)[i] << "\n";
                bad = true;
                break;
            }
        }
        tree.SetAutoFlush(100000);

        std::string line;
        while (!bad && std::getline(input, line)) {
            if (line.empty()) continue;
            std::istringstream fields(line);
            std::string token;
            for (std::size_t i = 0; i < names->size(); ++i) {
                if (!(fields >> token)) {
                    bad = true;
                    break;
                }
                char* end = nullptr;
                errno = 0;
                values[i] = std::strtod(token.c_str(), &end);
                if (end == token.c_str() || *end != '\0') {
                    bad = true;
                    break;
                }
            }
            if (!bad && (fields >> token)) bad = true;
            if (bad) {
                std::cerr << "Invalid column count or numeric token at selected row "
                          << (rows + 1) << "; expected " << names->size() << " columns\n";
                break;
            }
            tree.Fill();
            ++rows;
        }
        if (input.bad()) bad = true;
        if (!bad && tree.Write("", TObject::kOverwrite) <= 0) {
            std::cerr << "Failed to write PhysicsEvents\n";
            bad = true;
        }
        output.Close();
    }

    if (!bad && validate_written_file(argv[2], *names, rows)) {
        std::cout << "[DATA_ROOT_ROWS] " << rows << "\n";
        std::cout << "[DATA_ROOT_SCHEMA] mode=" << mode << " tree=PhysicsEvents branches="
                  << names->size() << " status=OK\n";
        return 0;
    }

    std::remove(argv[2]);
    return 4;
}
