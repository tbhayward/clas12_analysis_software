#ifndef ALL_BIN_DATA_H
#define ALL_BIN_DATA_H

#include <vector>
#include <string>

// Struct to hold bin data
struct AllBinData {
    int bin_number;
    int phi_min;
    int phi_max;
    double phi_avg;
    int yield_epg_FD_FD_inb;
    int yield_epg_CD_FD_inb;
    int yield_epg_CD_FT_inb;
    int yield_epg_FD_FD_outb;
    int yield_epg_CD_FD_outb;
    int yield_epg_CD_FT_outb;
    int yield_eppi0_FD_FD_inb;
    int yield_eppi0_CD_FD_inb;
    int yield_eppi0_CD_FT_inb;
    int yield_eppi0_FD_FD_outb;
    int yield_eppi0_CD_FD_outb;
    int yield_eppi0_CD_FT_outb;
};

// Function to read the all_bin_v3.csv file
std::vector<AllBinData> read_all_bin_v3(const std::string& filename);

// Function to print the bin data for debugging
void print_all_bin_data(const std::vector<AllBinData>& all_bin_data);

#endif  // ALL_BIN_DATA_H