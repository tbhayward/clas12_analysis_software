// split_by_runnum.cpp
// Build:  g++ -O2 split_by_runnum.cpp `root-config --cflags --libs` -o split_by_runnum
// Usage:  ./split_by_runnum /path/to/input.root [tree_name] [threshold]
//
// Defaults: tree_name = "PhysicsEvents", threshold = 17183
//
// Produces: /path/to/input_A.root  (runnum <= threshold)
//           /path/to/input_B.root  (runnum >  threshold)

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <iostream>
#include <string>

static std::string makeOutName(const std::string& in, const std::string& suffix) {
    // insert suffix before .root if present; else append suffix + .root
    const std::string ext = ".root";
    std::size_t pos = in.rfind(ext);
    if (pos == std::string::npos) return in + suffix + ext;
    return in.substr(0, pos) + suffix + ext;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.root> [tree_name] [threshold]\n";
        return 1;
    }

    const std::string inPath   = argv[1];
    const std::string treeName = (argc >= 3) ? argv[2] : "PhysicsEvents";
    const long long   thr      = (argc >= 4) ? std::stoll(argv[3]) : 17183LL;

    // Open input
    TFile* fin = TFile::Open(inPath.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file: " << inPath << "\n";
        return 2;
    }

    // Get tree
    TTree* tin = dynamic_cast<TTree*>(fin->Get(treeName.c_str()));
    if (!tin) {
        std::cerr << "Error: tree not found: " << treeName << "\n";
        fin->Close();
        return 3;
    }

    // Verify the 'runnum' branch exists
    if (!tin->GetBranch("runnum")) {
        std::cerr << "Error: branch 'runnum' not found in tree '" << treeName << "'.\n";
        fin->Close();
        return 4;
    }

    std::cout << "[INFO] Input file: " << inPath << "\n";
    std::cout << "[INFO] Tree name:  " << treeName << "\n";
    std::cout << "[INFO] Threshold:  " << thr << "\n";
    std::cout << "[INFO] Total entries in input: " << tin->GetEntries() << "\n";

    // Build selections
    std::string selA = "runnum<=" + std::to_string(thr);
    std::string selB = "runnum>"  + std::to_string(thr);

    // Prepare outputs
    const std::string outA = makeOutName(inPath, "_A");
    const std::string outB = makeOutName(inPath, "_B");

    // Copy A
    TFile* foutA = TFile::Open(outA.c_str(), "RECREATE");
    if (!foutA || foutA->IsZombie()) {
        std::cerr << "Error: cannot create output file: " << outA << "\n";
        fin->Close();
        return 5;
    }
    std::cout << "[INFO] Writing A (<= threshold) to: " << outA << "\n";
    TTree* tA = tin->CopyTree(selA.c_str());
    if (!tA) {
        std::cerr << "Error: CopyTree for selection '" << selA << "' returned null.\n";
        foutA->Close();
        fin->Close();
        return 6;
    }
    tA->Write();
    std::cout << "[INFO] Entries in A: " << tA->GetEntries() << "\n";
    foutA->Close();

    // Copy B
    TFile* foutB = TFile::Open(outB.c_str(), "RECREATE");
    if (!foutB || foutB->IsZombie()) {
        std::cerr << "Error: cannot create output file: " << outB << "\n";
        fin->Close();
        return 7;
    }
    std::cout << "[INFO] Writing B (> threshold) to: " << outB << "\n";
    TTree* tB = tin->CopyTree(selB.c_str());
    if (!tB) {
        std::cerr << "Error: CopyTree for selection '" << selB << "' returned null.\n";
        foutB->Close();
        fin->Close();
        return 8;
    }
    tB->Write();
    std::cout << "[INFO] Entries in B: " << tB->GetEntries() << "\n";
    foutB->Close();

    fin->Close();
    std::cout << "[DONE] Split complete.\n";
    return 0;
}