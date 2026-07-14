#ifndef CUT_VARIATION_SYSTEMATICS_H
#define CUT_VARIATION_SYSTEMATICS_H

#include <string>

struct CutVariationSystematicsOptions {
    bool enabled = false;
    bool apply_barlow = true;
    double barlow_threshold = 1.0;

    // Pass-1 cut-variation prescription: when the tight cross section differs
    // from nominal by more than this fraction, regard the tight result as
    // statistically unstable and estimate the systematic from the loose
    // variation only.  Otherwise use the arithmetic mean of the absolute
    // loose and tight deviations.
    bool use_pass1_tight_instability_rule = true;
    double tight_relative_difference_threshold = 0.50;
    bool make_plots = true;
    bool write_diagnostic_csv = true;

    std::string nominal_csv = "output/csvs/dvcs_pass2_analysis.csv";
    std::string exclusivity_loose_csv = "output/cut_variation_systematics/csv/exclusivity_95.csv";
    std::string exclusivity_tight_csv = "output/cut_variation_systematics/csv/exclusivity_99p99.csv";
    std::string fiducial_loose_csv = "output/cut_variation_systematics/csv/fiducial_loose.csv";
    std::string fiducial_tight_csv = "output/cut_variation_systematics/csv/fiducial_tight.csv";
    std::string output_dir = "output/cut_variation_systematics";

    std::string cross_section_column =
        "normed cross sections, ep->epg, exp, 10.6 GeV, unpol";
};

// Computes the two point-to-point cut systematics from five completed CSVs.
// The stored values are absolute uncertainties in nb/GeV^4.  By default the
// pass-1 prescription is used: 0.5*(|loose-nominal|+|tight-nominal|), except
// when |tight-nominal|/|nominal| > 0.5, where only |loose-nominal| is used.
// Barlow filtering is applied to the individual deviations before evaluating
// the final value; the unfiltered prescription is retained in the raw column.
// Cross-section cells are expected to be tuples: (value, statistical, systematic).
// The two varied exclusivity CSVs correspond to 95% and 99.99% upper-tail
// quantiles for one-sided cuts, together with the collaboration-approved loose
// and tight settings for all symmetric cuts.  The fiducial CSVs correspond to
// coherent loose and tight +/-2 degree angular-boundary variations.
bool update_cut_variation_systematics(const CutVariationSystematicsOptions& options);

#endif
