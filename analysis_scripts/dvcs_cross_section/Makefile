# Compiler and flags
CXX = g++
CXXFLAGS = $(shell root-config --cflags) -O2 -std=c++17 -pthread -fopenmp
LDFLAGS  = $(shell root-config --libs) -lMinuit2 -pthread -fopenmp
# ---------- dvcs_analysis ----------
DVCS_SRC = main.cpp initialize_pass2_csv.cpp load_trees.cpp make_dirs.cpp exclusivity_cuts.cpp \
		   global_cuts.cpp load_binning_scheme.cpp \
           bin_means.cpp total_counts.cpp yield_totals.cpp pi0_contamination.cpp pi0_corrected_counts.cpp \
           bsa.cpp radiative_corrections.cpp acceptance.cpp unfolding.cpp bin_volume.cpp \
           cross_sections.cpp model_predictions.cpp \
           bin_centering_corrections.cpp models_vs_data_plots.cpp
DVCS_OBJ = $(DVCS_SRC:.cpp=.o)
DVCS_BIN = dvcs_analysis
# ---------- cross_check_lee ----------
LEE_SRC  = cross_check_lee_main.cpp load_csv.cpp raw_signal_cross_check.cpp \
           pi0_contamination_cross_check.cpp rad_correction_cross_check.cpp \
           load_binning_scheme.cpp
LEE_OBJ  = $(LEE_SRC:.cpp=.o)
LEE_BIN  = cross_check_lee
# Build both by default
all: $(DVCS_BIN) $(LEE_BIN)
$(DVCS_BIN): $(DVCS_OBJ)
	$(CXX) $(DVCS_OBJ) -o $@ $(LDFLAGS)
$(LEE_BIN): $(LEE_OBJ)
	$(CXX) $(LEE_OBJ) -o $@ $(LDFLAGS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f $(DVCS_OBJ) $(DVCS_BIN) $(LEE_OBJ) $(LEE_BIN)