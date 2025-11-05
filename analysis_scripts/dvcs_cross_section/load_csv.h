#ifndef LOAD_CSV_H
#define LOAD_CSV_H

#include <string>
#include <vector>

struct LeeRow {
    // Bin coordinates
    double xBmin, xBmax;
    double Q2min, Q2max;
    double tmin,  tmax;           // |t| bin (abs(t))
    double phimin, phimax;        // from all_bin_v3
    double phiavg;                // center used for plotting

    // Raw yields (sum across three topologies, per your request)
    double raw_inb_sum;           // FD,FD + CD,FD + CD,FT (inbending)
    double raw_out_sum;           // FD,FD + CD,FD + CD,FT (outbending)
    double raw_combined;          // sum of the above

    // Pi0 contamination & signal (from all_bin_v3)
    double contam_inb;
    double contam_out;
    double signal_inb;
    double signal_out;

    // Acceptances (from full_acc.csv, matched by bin)
    double acc_inb;
    double acc_out;

    // Acceptance-corrected yield (combined)
    double acc_corr_yield;

    // Systematics / factors
    double Frad;
    double Fbin;
    double Vbin; // bin_volume

    // Debug
    bool   has_acceptance;
};

std::vector<LeeRow> load_lee_csvs(const std::string& all_bin_v3_path,
                                  const std::string& full_acc_path,
                                  int& matched_acc,
                                  int& unmatched_acc);

#endif // LOAD_CSV_H