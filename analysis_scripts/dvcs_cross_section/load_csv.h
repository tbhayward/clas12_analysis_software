// load_csv.h
#ifndef LOAD_CSV_H
#define LOAD_CSV_H

#include <string>
#include <vector>

struct LeeRow {
    // Bin definitions
    double xBmin=0, xBmax=0;
    double Q2min=0, Q2max=0;
    double tmin=0, tmax=0;          // absolute t limits (positive)
    double phimin=0, phimax=0;
    double phiavg=0;

    // Valid flag from all_bin_v3.csv (rows with valid==1 only are kept)
    bool   valid=true;

    // Raw yields (from all_bin_v3.csv)
    double raw_inb_fd_fd=0, raw_inb_cd_fd=0, raw_inb_cd_ft=0;
    double raw_out_fd_fd=0, raw_out_cd_fd=0, raw_out_cd_ft=0;

    // Convenience sums
    double raw_inb_sum=0, raw_out_sum=0, raw_combined=0;

    // Contamination ratios and signal yields (from all_bin_v3.csv)
    double contam_inb=0, contam_out=0;
    double signal_inb=0, signal_out=0;

    // Acceptance (from full_acc.csv; matched by bin)
    double acc_inb=0, acc_out=0;

    // Acceptance corrected yield (combined) and systematic factors
    double acc_corr_yld=0;     // "acceptance corrected yield, ep->epg, exp"
    double Frad=1.0, Fbin=1.0, Vbin=1.0; // Frad, Fbin, bin_volume
};

struct LeeData {
    std::vector<LeeRow> rows;
};

LeeData load_lee_csvs(const std::string& all_bin_csv,
                      const std::string& full_acc_csv,
                      bool verbose);

#endif // LOAD_CSV_H