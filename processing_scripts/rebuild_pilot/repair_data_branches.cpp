// Restore the original PhysicsEvents schema without reprocessing HIPO files.
// Usage: repair_data_branches ROOT_FILE SCHEMA_TEXT_FILE
#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TObject.h>
#include <TTree.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: repair_data_branches file.root schema.txt\n";
        return 2;
    }

    std::ifstream schema_file(argv[2]);
    if (!schema_file) {
        std::cerr << "Cannot open schema file: " << argv[2] << "\n";
        return 2;
    }

    std::vector<std::string> names;
    std::string name;
    while (schema_file >> name) {
        names.push_back(name);
    }
    if (names.empty()) {
        std::cerr << "Schema is empty\n";
        return 2;
    }

    TFile file(argv[1], "UPDATE");
    if (file.IsZombie()) {
        std::cerr << "Cannot open ROOT file: " << argv[1] << "\n";
        return 3;
    }

    TTree* tree = nullptr;
    file.GetObject("DataEvents", tree);
    if (!tree) {
        file.GetObject("PhysicsEvents", tree);
    }
    if (!tree) {
        std::cerr << "Neither DataEvents nor PhysicsEvents was found\n";
        return 4;
    }

    const auto nbranches = static_cast<size_t>(tree->GetListOfBranches()->GetEntries());
    if (nbranches != names.size()) {
        std::cerr << "Branch count mismatch: found " << nbranches
                  << ", expected " << names.size() << "\n";
        return 5;
    }

    bool generic_names = true;
    bool descriptive_names = true;
    for (size_t i = 0; i < names.size(); ++i) {
        auto* branch = static_cast<TBranch*>(tree->GetListOfBranches()->At(i));
        char generic[32];
        std::snprintf(generic, sizeof(generic), "col%02zu", i + 1);
        const std::string actual(branch->GetName());
        generic_names = generic_names && (actual == generic);
        descriptive_names = descriptive_names && (actual == names[i]);
        if (branch->GetListOfLeaves()->GetEntries() != 1) {
            std::cerr << "Unexpected number of leaves in branch " << actual << "\n";
            return 7;
        }
    }

    if (descriptive_names && std::string(tree->GetName()) == "PhysicsEvents") {
        std::cout << "[ALREADY_COMPATIBLE] " << argv[1] << "\n";
        return 0;
    }
    if (!generic_names && !descriptive_names) {
        std::cerr << "Mixed or unexpected branch names; refusing modification\n";
        return 6;
    }

    const Long64_t entries_before = tree->GetEntries();
    if (generic_names) {
        for (size_t i = 0; i < names.size(); ++i) {
            auto* branch = static_cast<TBranch*>(tree->GetListOfBranches()->At(i));
            auto* leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
            branch->SetName(names[i].c_str());
            leaf->SetName(names[i].c_str());
        }
    }

    tree->SetName("PhysicsEvents");
    file.cd();
    if (tree->Write("", TObject::kOverwrite) <= 0) {
        std::cerr << "Failed to write repaired tree metadata\n";
        return 8;
    }
    file.Close();

    TFile verify(argv[1], "READ");
    TTree* checked = nullptr;
    verify.GetObject("PhysicsEvents", checked);
    if (!checked || checked->GetEntries() != entries_before) {
        std::cerr << "Verification failed: tree missing or entry count changed\n";
        return 9;
    }
    for (const auto& expected : names) {
        if (!checked->GetBranch(expected.c_str())) {
            std::cerr << "Verification failed: missing branch " << expected << "\n";
            return 10;
        }
    }
    std::cout << "[REPAIRED] " << argv[1]
              << " entries=" << entries_before << "\n";
    return 0;
}
