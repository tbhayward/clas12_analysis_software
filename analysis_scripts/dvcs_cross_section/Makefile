.RECIPEPREFIX := >

# Compiler and flags
CXX = g++
CXXFLAGS = $(shell root-config --cflags) -O2 -std=c++17 -pthread -fopenmp
LDFLAGS  = $(shell root-config --libs) -lMinuit2 -pthread -fopenmp

# ---------- dvcs_analysis ----------
DVCS_SRC = main.cpp initialize_pass2_csv.cpp load_trees.cpp make_dirs.cpp exclusivity_cuts.cpp \
           global_cuts.cpp load_binning_scheme.cpp \
           bin_means.cpp total_counts.cpp current_dependence.cpp eppi0_normalization.cpp yield_totals.cpp pi0_contamination.cpp pi0_corrected_counts.cpp \
           bsa.cpp radiative_corrections.cpp acceptance.cpp unfolding.cpp bin_volume.cpp \
           cross_sections.cpp model_predictions.cpp overall_normalization_study.cpp \
           bin_centering_corrections.cpp models_vs_data_plots.cpp propagator_study.cpp \
           norm_cross_sections.cpp pass1_paper_plots.cpp branch_data_mc_comparison.cpp \
           branch_plot_config.cpp pi0_subtracted_kinematics.cpp external_scripts_runner.cpp cut_variation_systematics.cpp cut_variation_runner.cpp

DVCS_OBJ = $(DVCS_SRC:.cpp=.o)
DVCS_BIN = dvcs_analysis

# ---------- cross_check_lee ----------
LEE_SRC  = cross_check_lee_main.cpp load_csv.cpp raw_signal_cross_check.cpp \
           pi0_contamination_cross_check.cpp rad_correction_cross_check.cpp \
           load_binning_scheme.cpp unfolded_yields_cross_check.cpp bin_centering_cross_check.cpp \
           bin_volume_cross_check.cpp acceptance_cross_check.cpp cross_section_cross_check.cpp \
           cross_section_hayward_cross_check.cpp cross_section_run_period_consistency.cpp

LEE_OBJ  = $(LEE_SRC:.cpp=.o)
LEE_BIN  = cross_check_lee

# ---------- norm_scenario_comparison ----------
NORM_SCENARIO_SRC = norm_scenario_comparison.cpp
NORM_SCENARIO_OBJ = $(NORM_SCENARIO_SRC:.cpp=.o)
NORM_SCENARIO_BIN = norm_scenario_comparison

# ---------- main_systematics ----------
SYSTEMATICS_SRC = main_systematics.cpp pass1_systematics_import.cpp run_period_consistency_systematics.cpp combination_systematics.cpp sp19_inb_energy_scaling_systematics.cpp combination_point_to_point_systematics.cpp model_predictions.cpp
SYSTEMATICS_OBJ = $(SYSTEMATICS_SRC:.cpp=.o)
SYSTEMATICS_BIN = main_systematics

# Build all by default
all: $(DVCS_BIN) $(LEE_BIN) $(NORM_SCENARIO_BIN) $(SYSTEMATICS_BIN)

$(DVCS_BIN): $(DVCS_OBJ)
>$(CXX) $(DVCS_OBJ) -o $@ $(LDFLAGS)

$(LEE_BIN): $(LEE_OBJ)
>$(CXX) $(LEE_OBJ) -o $@ $(LDFLAGS)

$(NORM_SCENARIO_BIN): $(NORM_SCENARIO_OBJ)
>$(CXX) $(NORM_SCENARIO_OBJ) -o $@ $(LDFLAGS)

$(SYSTEMATICS_BIN): $(SYSTEMATICS_OBJ)
>$(CXX) $(SYSTEMATICS_OBJ) -o $@ $(LDFLAGS)

%.o: %.cpp
>$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
>rm -f $(DVCS_OBJ) $(DVCS_BIN) $(LEE_OBJ) $(LEE_BIN) $(NORM_SCENARIO_OBJ) $(NORM_SCENARIO_BIN) $(SYSTEMATICS_OBJ) $(SYSTEMATICS_BIN)
