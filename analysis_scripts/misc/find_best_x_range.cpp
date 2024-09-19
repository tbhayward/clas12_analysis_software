#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

// Function to find the x range with the highest statistics
std::pair<double, double> findBestXRange(TTree* tree, const char* condition, const char* xBranchName) {
    const double xMin = 0.06;
    const double xMax = 0.40;
    const double binWidth = 0.04;
    
    int numBins = static_cast<int>((xMax - xMin) / binWidth);
    std::vector<int> binCounts(numBins, 0);
    
    double x;
    tree->SetBranchAddress(xBranchName, &x);

    Long64_t nEntries = tree->GetEntries(condition);
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (x >= xMin && x <= xMax) {
            int bin = static_cast<int>((x - xMin) / binWidth);
            if (bin >= 0 && bin < numBins) {
                binCounts[bin]++;
            }
        }
    }
    
    // Find the bin with the highest count
    int maxBin = 0;
    for (int i = 1; i < numBins; ++i) {
        if (binCounts[i] > binCounts[maxBin]) {
            maxBin = i;
        }
    }

    double bestXMin = xMin + maxBin * binWidth;
    double bestXMax = bestXMin + binWidth;

    return std::make_pair(bestXMin, bestXMax);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_root_file>" << std::endl;
        return 1;
    }

    // Open the input ROOT file
    TFile* file = TFile::Open(argv[1]);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << argv[1] << std::endl;
        return 1;
    }

    // Load the PhysicsEvents tree
    TTree* tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: Cannot find PhysicsEvents tree in file " << argv[1] << std::endl;
        file->Close();
        return 1;
    }

    // Kinematic constraints
    const char* baseCondition = "Q2 > 1 && y < 0.80 && W > 2 && Mx2 > 0.16";
    
    // Condition for -(t-tmin) < 1
    const char* condition1 = " && -(t-tmin) < 1";
    std::string fullCondition1 = std::string(baseCondition) + condition1;

    // Condition for -(t-tmin) > 2
    const char* condition2 = " && -(t-tmin) > 2";
    std::string fullCondition2 = std::string(baseCondition) + condition2;

    // Find the best x range for each condition
    std::pair<double, double> bestRange1 = findBestXRange(tree, fullCondition1.c_str(), "x");
    std::pair<double, double> bestRange2 = findBestXRange(tree, fullCondition2.c_str(), "x");

    // Output the results
    std::cout << "Best x range for -(t-tmin) < 1: [" << bestRange1.first << ", " << bestRange1.second << "]" << std::endl;
    std::cout << "Best x range for -(t-tmin) > 2: [" << bestRange2.first << ", " << bestRange2.second << "]" << std::endl;

    // Clean up
    file->Close();
    return 0;
}