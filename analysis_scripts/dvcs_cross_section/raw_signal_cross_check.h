#ifndef RAW_SIGNAL_CROSS_CHECK_H
#define RAW_SIGNAL_CROSS_CHECK_H

#include <string>

// Make raw-yield comparison plots (Hayward vs. Lee) and the ratio (Hayward/Lee).
//
// This version reads *directly* from the two CSVs:
//
//   - lee_csv_path:
//       Lee pass-1 CSV, e.g. "imports/all_bin_v3.csv"
//       Uses:
//         "bin index", "valid bin",
//         "xBmin", "xBmax", "Q2min", "Q2max", "tmin", "tmax", "phiavg",
//         "raw yield, ep->epg, (FD, FD), exp, inbending",
//         "raw yield, ep->epg, (CD, FD), exp, inbending",
//         "raw yield, ep->epg, (CD, FT), exp, inbending",
//         "raw yield, ep->epg, (FD, FD), exp, outbending",
//         "raw yield, ep->epg, (CD, FD), exp, outbending",
//         "raw yield, ep->epg, (CD, FT), exp, outbending"
//
//   - hayward_csv_path:
//       Hayward pass-2 CSV, e.g. "output/csvs/dvcs_pass2_analysis.csv"
//       Uses (unpolarized only):
//         "bin index", "valid bin",
//         "raw yield, ep->epg, (FD, FD), exp, Fa18 Inb, unpol",
//         "raw yield, ep->epg, (CD, FD), exp, Fa18 Inb, unpol",
//         "raw yield, ep->epg, (CD, FT), exp, Fa18 Inb, unpol",
//         "raw yield, ep->epg, (FD, FD), exp, Fa18 Out, unpol",
//         "raw yield, ep->epg, (CD, FD), exp, Fa18 Out, unpol",
//         "raw yield, ep->epg, (CD, FT), exp, Fa18 Out, unpol"
//
// For each xB bin it produces in output_base_dir:
//
//   raw_counts_fa18_inb_xB_<ix>.png
//   raw_ratio_fa18_inb_xB_<ix>.png
//   raw_counts_fa18_out_xB_<ix>.png
//   raw_ratio_fa18_out_xB_<ix>.png
//
void plot_raw_yield_cross_checks(const std::string& lee_csv_path,
                                 const std::string& hayward_csv_path,
                                 const std::string& output_base_dir);

#endif // RAW_SIGNAL_CROSS_CHECK_H