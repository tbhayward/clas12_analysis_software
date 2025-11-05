// load_csv.h
#ifndef LOAD_CSV_H
#define LOAD_CSV_H

#include <map>
#include <string>
#include <vector>

struct LeeRow {
    // bin edges (use absolute t)
    double xBmin = 0.0, xBmax = 0.0;
    double Q2min = 0.0, Q2max = 0.0;
    double t_abs_min = 0.0, t_abs_max = 0.0;
    double phimin = 0.0, phimax = 0.0, phiavg = 0.0;

    // raw yields (Fa18 inb/out, three topologies, sum later when comparing)
    double raw_inb_fd_fd = 0.0;
    double raw_inb_cd_fd = 0.0;
    double raw_inb_cd_ft = 0.0;
    double raw_out_fd_fd = 0.0;
    double raw_out_cd_fd = 0.0;
    double raw_out_cd_ft = 0.0;

    // pi0 contamination / signal yields
    double contam_ratio_inb = 0.0;
    double contam_ratio_out = 0.0;
    double signal_yield_inb = 0.0;
    double signal_yield_out = 0.0;

    // acceptance (from full_acc.csv)
    double acceptance_inb = 0.0;
    double acceptance_out = 0.0;

    // acceptance-corrected yield (combined, from all_bin_v3.csv)
    double acc_corr_yield_combined = 0.0;

    // other systematics / factors
    double Frad = 0.0;
    double Fbin = 0.0;
    double bin_volume = 0.0;

    // optional cross section (not used yet)
    double xsec = 0.0;

    // derived helpers
    double raw_inb_sum() const {
        return raw_inb_fd_fd + raw_inb_cd_fd + raw_inb_cd_ft;
    }
    double raw_out_sum() const {
        return raw_out_fd_fd + raw_out_cd_fd + raw_out_cd_ft;
    }
    double raw_combined() const {
        return raw_inb_sum() + raw_out_sum();
    }
};

struct LeeData {
    std::vector<LeeRow> rows;
    std::map<std::string, size_t> index_by_key; // bin key -> row index
};

// Build a unique key from bin edges (rounded to keep string stable)
std::string make_bin_key(double xBmin,double xBmax,
                         double Q2min,double Q2max,
                         double tmin,double tmax,
                         double phimin,double phimax);

// Load and merge the two CSVs; if verbose=true, print row-by-row.
LeeData load_lee_csvs(const std::string& all_csv,
                      const std::string& full_acc_csv,
                      bool verbose);

#endif // LOAD_CSV_H