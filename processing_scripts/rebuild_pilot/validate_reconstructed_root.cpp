#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include "reconstructed_data_schema.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

std::string positional_name(std::size_t index) {
    std::ostringstream stream;
    stream << "col" << std::setw(2) << std::setfill('0') << (index + 1);
    return stream.str();
}

bool validate_branches(TTree* tree, const std::vector<std::string>& expected,
                       bool allow_legacy, std::string& layout) {
    auto* branches = tree->GetListOfBranches();
    if (!branches || branches->GetEntries() != static_cast<int>(expected.size())) {
        std::cerr << "Branch count mismatch: got "
                  << (branches ? branches->GetEntries() : -1) << ", expected "
                  << expected.size() << "\n";
        return false;
    }

    bool descriptive = true;
    bool positional = allow_legacy;
    std::unordered_set<std::string> seen;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        auto* branch = dynamic_cast<TBranch*>(branches->At(static_cast<int>(i)));
        if (!branch) return false;
        const std::string actual = branch->GetName();
        if (!seen.insert(actual).second) {
            std::cerr << "Duplicate branch: " << actual << "\n";
            return false;
        }
        descriptive = descriptive && actual == expected[i];
        positional = positional && actual == positional_name(i);
    }
    if (!descriptive && !positional) {
        std::cerr << "Branch names or order do not match the " << expected.size()
                  << "-column schema\n";
        return false;
    }

    layout = descriptive ? "descriptive" : "legacy-positional";
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const std::string actual = descriptive ? expected[i] : positional_name(i);
        auto* branch = tree->GetBranch(actual.c_str());
        auto* leaf = branch ? branch->GetLeaf(actual.c_str()) : nullptr;
        if (!leaf || std::string(leaf->GetTypeName()) != "Double_t" ||
            leaf->GetLeafCount() != nullptr || leaf->GetLenStatic() != 1) {
            std::cerr << "Branch " << actual << " is not a scalar Double_t\n";
            return false;
        }
        if (!descriptive) {
            const char* alias = tree->GetAlias(expected[i].c_str());
            if (alias && std::string(alias) != actual) {
                std::cerr << "Incorrect alias for " << expected[i] << "\n";
                return false;
            }
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4 || argc > 5) {
        std::cerr << "usage: validate_reconstructed_root ROOT_FILE "
                     "dvcs_data|pi0_data|pi0_as_dvcs EXPECTED_ROWS [--allow-legacy]\n";
        return 2;
    }
    const std::string path = argv[1];
    const std::string mode = argv[2];
    const auto* expected = reconstructed_data_schema::for_mode(mode);
    if (!expected) {
        std::cerr << "Unknown mode: " << mode << "\n";
        return 2;
    }
    if (!reconstructed_data_schema::definition_is_valid(mode, *expected)) {
        std::cerr << "Internal schema definition is invalid for " << mode << "\n";
        return 2;
    }
    char* end = nullptr;
    const long long expected_rows = std::strtoll(argv[3], &end, 10);
    if (!end || *end != '\0' || expected_rows < -1) {
        std::cerr << "EXPECTED_ROWS must be -1 or a nonnegative integer\n";
        return 2;
    }
    const bool allow_legacy = argc == 5 && std::string(argv[4]) == "--allow-legacy";
    if (argc == 5 && !allow_legacy) {
        std::cerr << "Unknown option: " << argv[4] << "\n";
        return 2;
    }

    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open ROOT file: " << path << "\n";
        return 1;
    }

    TTree* tree = nullptr;
    file->GetObject("PhysicsEvents", tree);
    std::string tree_name = "PhysicsEvents";
    if (!tree && allow_legacy) {
        file->GetObject("DataEvents", tree);
        tree_name = "DataEvents";
    }
    if (!tree) {
        std::cerr << "Missing PhysicsEvents"
                  << (allow_legacy ? " (and legacy DataEvents)" : "") << " tree\n";
        return 1;
    }

    const long long rows = tree->GetEntries();
    if (expected_rows >= 0 && rows != expected_rows) {
        std::cerr << "Row count mismatch: got " << rows << ", expected " << expected_rows << "\n";
        return 1;
    }
    std::string layout;
    if (!validate_branches(tree, *expected, allow_legacy, layout)) return 1;

    std::cout << "[ROOT_SCHEMA_OK] file=" << path << " mode=" << mode
              << " tree=" << tree_name << " layout=" << layout << " rows=" << rows
              << " branches=" << expected->size() << "\n";
    return 0;
}
