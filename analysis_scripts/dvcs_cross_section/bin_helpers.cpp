// bin_helpers.cpp

#include "bin_helpers.h"
#include <algorithm>
#include <cctype>
#include <cmath> 
#include <iostream>  // Added for debugging

// Function implementation
std::string clean_bin_label(const std::string& label) {
    std::string clean_label = label;
    clean_label.erase(std::remove_if(clean_label.begin(), clean_label.end(), [](unsigned char c) {
        return std::isspace(c) || c == '(' || c == ')';
    }), clean_label.end());
    return clean_label;
}

std::vector<int> precompute_relevant_bins(int xB_bin, const std::vector<BinBoundary>& bin_boundaries) {
    std::vector<int> relevant_bins;
    std::cout << "Entering precompute_relevant_bins with xB_bin = " << xB_bin << std::endl;
    
    for (size_t bin_idx = 0; bin_idx < bin_boundaries.size(); ++bin_idx) {
        std::string bin_label = clean_bin_label(bin_boundaries[bin_idx].bin_label);
        size_t first_comma = bin_label.find(',');

        std::cout << "Processing bin_idx: " << bin_idx << ", bin_label: '" << bin_label << "'" << std::endl;

        if (first_comma != std::string::npos) {
            try {
                std::string xB_substr = bin_label.substr(0, first_comma);
                int xB_label = std::stoi(xB_substr);
                std::cout << "Extracted xB_label: " << xB_label << std::endl;
                
                if (xB_label == xB_bin) {
                    relevant_bins.push_back(bin_idx);
                    std::cout << "Bin_idx " << bin_idx << " added to relevant_bins." << std::endl;
                }
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Cannot convert '" << bin_label.substr(0, first_comma) << "' to integer for bin_idx: " << bin_idx << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Error: Integer out of range when processing bin_idx: " << bin_idx << std::endl;
            }
        } else {
            std::cerr << "Warning: No comma found in bin_label for bin_idx: " << bin_idx << std::endl;
        }
    }

    std::cout << "Total relevant bins found: " << relevant_bins.size() << std::endl;
    return relevant_bins;
}

// Find the next perfect square greater than or equal to the given number
int next_perfect_square(int n) {
    int square_root = std::ceil(std::sqrt(n));
    return square_root * square_root;
}