#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm> // For std::remove

// Function to read arrays directly from the file
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
            while (std::getline(ss, rowStr, ',')) {
                std::vector<double> row;
                std::stringstream rowStream(rowStr);
                std::string value;
                while (std::getline(rowStream, value, ',')) {
                    row.push_back(std::stod(value));
                }
                values.push_back(row);
            }

            kinematicData[key] = values;
        }
    }

    return kinematicData;
}

// Function to print the data for verification
void printKinematicData(const std::map<std::string, std::vector<std::vector<double>>> &kinematicData) {
    for (const auto &entry : kinematicData) {
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
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <kinematicPlots.txt>\n";
        return 1;
    }

    std::string kinematicFile = argv[1];

    // Read the kinematic data from the file
    std::map<std::string, std::vector<std::vector<double>>> kinematicData = readKinematics(kinematicFile);

    // Print out the parsed data
    printKinematicData(kinematicData);

    return 0;
}
