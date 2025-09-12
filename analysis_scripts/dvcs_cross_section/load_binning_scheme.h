#ifndef LOAD_BINNING_SCHEME_H
#define LOAD_BINNING_SCHEME_H

#include <string>
#include <vector>

// Matches your Python namedtuple fields and order.
struct Binning {
    double xBmin;
    double xBmax;
    double Q2min;
    double Q2max;
    double tmin;   // |tmin|
    double tmax;   // |tmax|
};

// Parse the CSV at `csv_file_path` (expected: imports/integrated_bin_v2.csv).
// Skips the first two header lines. Splits on whitespace, and extracts
// columns 4,5,7,8,10,11 as doubles (with t made positive).
std::vector<Binning> load_binning_scheme(const std::string& csv_file_path);

#endif // LOAD_BINNING_SCHEME_H