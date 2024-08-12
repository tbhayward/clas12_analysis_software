#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <regex>

// Define a struct to hold the asymmetry fit data
struct AsymmetryFit {
    std::string fitType;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> yErr;
};

// Define a struct to hold the kinematic data
struct Kinematics {
    std::vector<std::vector<double>> kinematics;
};

// Function to parse asymmetry fits from the file
std::map<std::string, AsymmetryFit> parseAsymmetries(const std::string &filename) {
    std::map<std::string, AsymmetryFit> asymmetryData;
    std::ifstream file(filename);
    std::string line;

    std::regex pattern(R"((\w+)chi2Fits(\w+) = \{\{(.+?)\}\};)");
    std::smatch match;

    while (std::getline(file, line)) {
        if (std::regex_search(line, match, pattern)) {
            std::string key = match[1].str();
            std::string fitType = match[2].str();
            std::string dataStr = match[3].str();
            
            AsymmetryFit fit;
            fit.fitType = fitType;

            // Split the dataStr into individual data points
            std::regex dataPattern(R"(\{(\d+\.?\d*), (\d+\.?\d*), (\d+\.?\d*)\})");
            std::smatch dataMatch;
            std::string::const_iterator searchStart(dataStr.cbegin());

            while (std::regex_search(searchStart, dataStr.cend(), dataMatch, dataPattern)) {
                fit.x.push_back(std::stod(dataMatch[1].str()));
                fit.y.push_back(std::stod(dataMatch[2].str()));
                fit.yErr.push_back(std::stod(dataMatch[3].str()));
                searchStart = dataMatch.suffix().first;
            }

            asymmetryData[key + fitType] = fit;
        }
    }

    return asymmetryData;
}

// Function to parse kinematics from the file
std::map<std::string, Kinematics> parseKinematics(const std::string &filename) {
    std::map<std::string, Kinematics> kinematicData;
    std::ifstream file(filename);
    std::string line;

    std::regex pattern(R"((\w+Kinematics) = \{\{(.+?)\}\};)");
    std::smatch match;

    while (std::getline(file, line)) {
        if (std::regex_search(line, match, pattern)) {
            std::string key = match[1].str();
            std::string dataStr = match[2].str();

            Kinematics kin;

            // Split the dataStr into individual kinematic points
            std::regex dataPattern(R"(\{(.+?)\})");
            std::smatch dataMatch;
            std::string::const_iterator searchStart(dataStr.cbegin());

            while (std::regex_search(searchStart, dataStr.cend(), dataMatch, dataPattern)) {
                std::vector<double> kinValues;
                std::stringstream ss(dataMatch[1].str());
                std::string value;
                while (std::getline(ss, value, ',')) {
                    kinValues.push_back(std::stod(value));
                }
                kin.kinematics.push_back(kinValues);
                searchStart = dataMatch.suffix().first;
            }

            kinematicData[key] = kin;
        }
    }

    return kinematicData;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <asymmetries.txt> <kinematicPlots.txt>\n";
        return 1;
    }

    std::string asymmetryFile = argv[1];
    std::string kinematicFile = argv[2];

    // Parse the files
    std::map<std::string, AsymmetryFit> asymmetryData = parseAsymmetries(asymmetryFile);
    std::map<std::string, Kinematics> kinematicData = parseKinematics(kinematicFile);

    // Print out the parsed data
    std::cout << "Asymmetry Data:\n";
    for (const auto& [key, fit] : asymmetryData) {
        std::cout << key << ": " << fit.fitType << "\n";
        for (size_t i = 0; i < fit.x.size(); ++i) {
            std::cout << "  (" << fit.x[i] << ", " << fit.y[i] << " ± " << fit.yErr[i] << ")\n";
        }
    }

    std::cout << "\nKinematic Data:\n";
    for (const auto& [key, kin] : kinematicData) {
        std::cout << key << ":\n";
        for (const auto& kinSet : kin.kinematics) {
            std::cout << "  {";
            for (const auto& val : kinSet) {
                std::cout << val << ", ";
            }
            std::cout << "}\n";
        }
    }

    return 0;
}
