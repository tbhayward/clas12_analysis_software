#include "load_binning_scheme.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
    // Split on any whitespace (like Python's .split()).
    static inline std::vector<std::string> split_ws(const std::string& line) {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string tok;
        while (iss >> tok) tokens.push_back(tok);
        return tokens;
    }
}

std::vector<Binning> load_binning_scheme(const std::string& csv_file_path) {
    std::vector<Binning> bins;

    std::ifstream in(csv_file_path);
    if (!in) {
        std::cerr << "[load_binning_scheme] ERROR: cannot open file: "
                  << csv_file_path << std::endl;
        return bins;
    }

    // Read all lines first so we can easily skip two header lines
    std::vector<std::string> lines;
    std::string line;
    lines.reserve(1024);
    while (std::getline(in, line)) {
        // Handle possible CRLF
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    in.close();

    if (lines.size() <= 2) {
        std::cerr << "[load_binning_scheme] WARNING: file has <= 2 lines; "
                     "no data rows found after headers."
                  << std::endl;
        return bins;
    }

    // Skip first two header lines; parse the rest
    for (std::size_t i = 2; i < lines.size(); ++i) {
        const std::string& ln = lines[i];
        if (ln.empty()) continue;                  // skip blank lines
        if (!ln.empty() && ln[0] == '#') continue; // skip comment lines

        auto tokens = split_ws(ln);
        // Need at least up to index 11 -> size >= 12
        if (tokens.size() < 12) {
            // Not enough tokens to parse a valid row
            continue;
        }

        try {
            // Columns per your Python docstring:
            // 4: xBmin, 5: xBmax, 7: Q2min, 8: Q2max, 10: |tmin|, 11: |tmax|
            double xBmin = std::stod(tokens[4]);
            double xBmax = std::stod(tokens[5]);
            double Q2min = std::stod(tokens[7]);
            double Q2max = std::stod(tokens[8]);
            double tmin  = std::fabs(std::stod(tokens[10]));
            double tmax  = std::fabs(std::stod(tokens[11]));

            bins.push_back(Binning{xBmin, xBmax, Q2min, Q2max, tmin, tmax});
        } catch (const std::exception& e) {
            // i is 0-based in `lines`; +1 to get file line, +1 again for human-friendly
            std::size_t file_line = i + 1;
            std::cerr << "[load_binning_scheme] Parse error on line "
                      << file_line << ": \"" << ln << "\" -> " << e.what()
                      << std::endl;
            continue;
        }
    }

    return bins;
}