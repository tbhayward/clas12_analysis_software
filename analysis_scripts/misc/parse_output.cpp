#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>

int main() {
    std::ifstream inFile("/volatile/clas12/thayward/UU_validation/packages/1/clas12_analysis_software/analysis_scripts/asymmetry_extraction/output/results/kinematics_rga_fa18_inb_epX_skimmed_timeStamp_03_06_104135.txt"); // Adjust the path as necessary
    if (!inFile) {
        std::cerr << "Unable to open file" << std::endl;
        return 1;
    }

    std::ofstream outFile("output2/analysis_statistics.txt");
    if (!outFile) {
        std::cerr << "Unable to create output file" << std::endl;
        return 1;
    }

    std::string line;
    bool inTable = false;
    std::regex tableStartRegex(R"(\\begin{tabular})"), tableEndRegex(R"(\\end{tabular})"), 
        captionRegex(R"(\\caption\{.*Q2y(\d)z(\d).*\})"), 
        dataRegex(R"(^1~&~([0-9\.]+)~&~([0-9\.]+)~&~([0-9\.]+)~&~([0-9\.]+)~&~([0-9\.]+)~&~([0-9\.]+)~&~([0-9\.]+)~&~([0-9\.]+)~&~([0-9\.]+)~&~([0-9\.]+)\\ \\hline)");

    // Assuming we are dealing with only one set of kinematic variables for simplicity
    std::vector<float> Q2(25, -100), W(25, -100), xB(25, -100), y(25, -100), z(25, -100), zeta(25, -100), PT(25, -100), xF(25, -100), t(25, -100);

    while (getline(inFile, line)) {
        if (std::regex_search(line, tableStartRegex)) {
            inTable = true;
        } else if (inTable && std::regex_search(line, tableEndRegex)) {
            inTable = false;
        } else if (inTable) {
            std::smatch matches;
            if (std::regex_search(line, matches, dataRegex)) {
                for (size_t i = 1; i < matches.size(); ++i) {
                    // Convert match to float and store in appropriate vector
                    // This example just prints the matches to demonstrate parsing
                    std::cout << "Match " << i << ": " << matches[i] << std::endl;
                }
            }
        }
    }

    // Example of outputting the Q2 vector values to the file
    for (size_t i = 0; i < Q2.size(); ++i) {
        outFile << "Q2[" << i << "] = " << Q2[i] << std::endl;
        // Repeat for other variables
    }

    inFile.close();
    outFile.close();

    return 0;
}
