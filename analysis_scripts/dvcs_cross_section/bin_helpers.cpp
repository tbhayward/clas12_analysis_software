#include "bin_helpers.h"
#include <algorithm>
#include <cctype>
#include <cmath> 

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
    for (size_t bin_idx = 0; bin_idx < bin_boundaries.size(); ++bin_idx) {
        std::string bin_label = clean_bin_label(bin_boundaries[bin_idx].bin_label);
        size_t first_comma = bin_label.find(',');

        if (first_comma != std::string::npos) {
            int xB_label = std::stoi(bin_label.substr(0, first_comma));
            if (xB_label == xB_bin) {
                relevant_bins.push_back(bin_idx);
            }
        }
    }
    return relevant_bins;
}

// Find the next perfect square greater than or equal to the given number
int next_perfect_square(int n) {
    int square_root = std::ceil(std::sqrt(n));
    return square_root * square_root;
}