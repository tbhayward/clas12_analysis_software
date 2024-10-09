#include "bin_boundaries.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<BinBoundary> read_bin_boundaries(const std::string& filename) {
    std::vector<BinBoundary> bin_boundaries;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return bin_boundaries;
    }

    std::string line;
    int line_num = 0;
    int sequential_bin_number = 1;  // Initialize sequential bin number starting from 1

    // Skip the first two lines (header and descriptive line)
    while (line_num < 2 && std::getline(file, line)) {
        ++line_num;
    }

    // Read data from the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        BinBoundary bin;

        // Assign bin number sequentially (this is an integer)
        bin.bin_number = sequential_bin_number++;

        // Second column: bin name (read it as a string)
        std::getline(ss, token, '\t');
        bin.bin_label = token;

        try {
            // xB_low (min)
            std::getline(ss, token, '\t');
            bin.xB_low = std::stod(token);

            // xB_high (max)
            std::getline(ss, token, '\t');
            bin.xB_high = std::stod(token);

            // xB_avg (average)
            std::getline(ss, token, '\t');
            bin.xB_avg = std::stod(token);

            // Q2_low (min)
            std::getline(ss, token, '\t');
            bin.Q2_low = std::stod(token);

            // Q2_high (max)
            std::getline(ss, token, '\t');
            bin.Q2_high = std::stod(token);

            // Q2_avg (average)
            std::getline(ss, token, '\t');
            bin.Q2_avg = std::stod(token);

            // t_low (min)
            std::getline(ss, token, '\t');
            bin.t_low = std::stod(token);

            // t_high (max)
            std::getline(ss, token, '\t');
            bin.t_high = std::stod(token);

            // t_avg (average)
            std::getline(ss, token, '\t');
            bin.t_avg = std::stod(token);

        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument while parsing token: '" << token << "' in bin " << bin.bin_label << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range error while parsing token: '" << token << "' in bin " << bin.bin_label << std::endl;
        }

        bin_boundaries.push_back(bin);
    }

    file.close();
    return bin_boundaries;
}