#ifndef CORRELATED_SCALE_SYSTEMATICS_H
#define CORRELATED_SCALE_SYSTEMATICS_H

#include <string>

struct CorrelatedScaleSystematicsOptions {
    double theta_bin_width_deg = 4.0;
    int min_ratio_points_per_period = 25;
    // Pre-existing target-thickness ⊕ accumulated-charge normalization.
    double uncorrelated_normalization_fraction = 0.021633307652784;

    // Proton-efficiency overall-normalization prescription.
    double proton_efficiency_fa18_fraction = 0.0283;
    double proton_efficiency_sp18_fraction = 0.1132;
    double proton_efficiency_sp19_fraction = 0.0283;

    // Exact accumulated charges for the final Pass-2 10.6-GeV run selection.
    double charge_sp18_inb_mC = 51.248191;
    double charge_sp18_out_mC = 11.435592;
    double charge_fa18_inb_mC = 29.407050;
    double charge_fa18_out_mC = 31.900540;
    std::string output_dir =
        "output/systematics/correlated_scale";
};

bool correlated_scale_systematics(
    const std::string& csv_path,
    const CorrelatedScaleSystematicsOptions& options =
        CorrelatedScaleSystematicsOptions());

#endif
