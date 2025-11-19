#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include <string>
#include <map>

// Simple triple for (value, stat, sys)
struct Triple {
    double value;
    double stat;
    double sys;
};

// Map from period / group label to luminosity triple.
// Keys should be things like "Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Fa18", "Sp18", "10.6 GeV".
using LumiMap = std::map<std::string, Triple>;

/**
 * Build the luminosity map used to fill the integrated luminosity columns.
 *
 * IMPORTANT:
 *   - Fill the 'value' field with your actual integrated luminosities in (nC).
 *   - 'stat' can be 0 (treat lumi as systematic later if desired).
 *   - 'sys' can be 0 for now.
 */
LumiMap build_lumi_map();

/**
 * Update dvcs_pass2_analysis.csv in-place:
 *   - Fills the integrated luminosity columns for all rows:
 *       integrated luminosity, Fa18 Inb (nC)
 *       integrated luminosity, Fa18 Out (nC)
 *       integrated luminosity, Sp19 Inb (nC)
 *       integrated luminosity, Sp18 Inb (nC)
 *       integrated luminosity, Sp18 Out (nC)
 *       integrated luminosity, Fa18 (nC)
 *       integrated luminosity, Sp18 (nC)
 *       integrated luminosity, 10.6 GeV (nC)
 *     (and optionally "integrated luminosity, Fa18 Inb Supp (nC)" if that column exists).
 *
 *   - Computes cross sections for all labels/helicities where columns exist:
 *        cross section, ep->epg, exp, <label>, <helicity>
 *     using:
 *        unfolded yield, ep->epg, exp, <label>, <helicity>
 *        Frad, 10.6 GeV
 *        bin volume, 10.6 GeV
 *        integrated luminosity, <label> (nC)
 *
 *   - Writes cross section columns as "(value, stat, sys)" with sys = 0.
 *   - Propagates statistical uncertainties from the unfolded yield and Frad.
 *
 *   - Computes theory curves at the bin-mean kinematics for each row and label,
 *     sampling 12 equally spaced points in phi from 0 to 2*pi, and writes these
 *     into JSON files under:
 *       output/jsons/cross_sections/<PeriodDir>/...
 *
 * This function does no plotting.
 *
 * Returns true on success, false on fatal error.
 */
bool compute_cross_sections(const std::string &csv_main,
                            const LumiMap &lumi_map);

/**
 * Plot cross sections vs phi for a given label (period or combined group).
 *
 * - Reads cross sections from the CSV (three-tuple columns).
 * - Reads the precomputed theory curves from the JSON cache written by
 *   compute_cross_sections().
 * - Produces “xB canvases” with three pads:
 *     (unpolarized, positive helicity, negative helicity)
 *   showing data (points with error bars) and theory curves
 *   (BH, KM, VGG) as lines through the 12 equidistant phi values.
 *
 * Plots are saved in:
 *   output/cross_sections/<PeriodDir>/cross_sections_<PeriodDir>_xB_<idx>.png
 *
 * 'label' should be one of:
 *   "Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp", "Sp18 Inb",
 *   "Sp18 Out", "Sp19 Inb", "Fa18", "Sp18", "10.6 GeV"
 *
 * theory_json_root:
 *   Root directory where compute_cross_sections() stored JSON:
 *     e.g. "output/jsons/cross_sections"
 *
 * out_root_dir:
 *   Root directory for plots, e.g. "output/cross_sections".
 *
 * Returns true on success, false on fatal error.
 */
bool plot_cross_sections_for_label(const std::string &csv_main,
                                   const std::string &label,
                                   const std::string &theory_json_root,
                                   const std::string &out_root_dir);

#endif // CROSS_SECTIONS_H