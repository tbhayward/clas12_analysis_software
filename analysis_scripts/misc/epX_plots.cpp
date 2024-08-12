#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm> // For std::remove

// Function to read arrays directly from the kinematic file
std::map<std::string, std::vector<std::vector<double>>> readKinematics(const std::string &filename) {
    std::map<std::string, std::vector<std::vector<double>>> kinematicData;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::string::size_type pos = line.find("=");
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos - 1);
            std::string dataStr = line.substr(pos + 1);

            // Remove curly braces and semicolon from the data string
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '{'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '}'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), ';'), dataStr.end());

            // Split dataStr into individual rows of data
            std::vector<std::vector<double>> values;
            std::stringstream ss(dataStr);
            std::string rowStr;
            while (std::getline(ss, rowStr, '}')) {
                rowStr.erase(std::remove(rowStr.begin(), rowStr.end(), '{'), rowStr.end());
                rowStr.erase(std::remove(rowStr.begin(), rowStr.end(), ','), rowStr.end());

                std::stringstream rowStream(rowStr);
                std::string value;
                std::vector<double> row;
                while (std::getline(rowStream, value, ',')) {
                    row.push_back(std::stod(value));
                }
                if (!row.empty()) values.push_back(row);
            }

            kinematicData[key] = values;
        }
    }

    return kinematicData;
}

// Function to read asymmetry fits directly from the file
std::map<std::string, std::vector<std::vector<double>>> readAsymmetries(const std::string &filename) {
    std::map<std::string, std::vector<std::vector<double>>> asymmetryData;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::string::size_type pos = line.find("=");
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos - 1);
            std::string dataStr = line.substr(pos + 1);

            // Remove curly braces and semicolon from the data string
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '{'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '}'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), ';'), dataStr.end());

            // Split dataStr into individual rows of data
            std::vector<std::vector<double>> values;
            std::stringstream ss(dataStr);
            std::string rowStr;
            while (std::getline(ss, rowStr, '}')) {
                rowStr.erase(std::remove(rowStr.begin(), rowStr.end(), '{'), rowStr.end());
                rowStr.erase(std::remove(rowStr.begin(), rowStr.end(), ','), rowStr.end());

                std::stringstream rowStream(rowStr);
                std::string value;
                std::vector<double> row;
                while (std::getline(rowStream, value, ',')) {
                    row.push_back(std::stod(value));
                }
                if (!row.empty()) values.push_back(row);
            }

            asymmetryData[key] = values;
        }
    }

    return asymmetryData;
}

// Function to print the data for verification
void printData(const std::map<std::string, std::vector<std::vector<double>>> &data) {
    for (const auto &entry : data) {
        std::cout << entry.first << " = {\n";
        for (const auto &vec : entry.second) {
            std::cout << "  {";
            for (size_t i = 0; i < vec.size(); ++i) {
                std::cout << vec[i];
                if (i < vec.size() - 1) std::cout << ", ";
            }
            std::cout << "},\n";
        }
        std::cout << "};\n";
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <asymmetries.txt> <kinematicPlots.txt>\n";
        return 1;
    }

    std::string asymmetryFile = argv[1];
    std::string kinematicFile = argv[2];

    // Read the asymmetry data from the file
    std::map<std::string, std::vector<std::vector<double>>> asymmetryData = readAsymmetries(asymmetryFile);

    // Read the kinematic data from the file
    std::map<std::string, std::vector<std::vector<double>>> kinematicData = readKinematics(kinematicFile);

    // Print out the parsed data
    std::cout << "Asymmetry Data:\n";
    printData(asymmetryData);

    std::cout << "\nKinematic Data:\n";
    printData(kinematicData);

    return 0;
}
