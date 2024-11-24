# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -g -Wno-unused-but-set-variable -std=c++17
ROOTLIBS = $(shell root-config --cflags --libs)
LFLAGS = $(ROOTLIBS) -lMinuit -lstdc++fs  # Add -lstdc++fs for filesystem

# Target executable
TARGET = dvcs_cross_section_extraction

# Source files (added calculate_contamination.cpp)
SRC = dvcs_cross_section_extraction.cpp \
      determine_exclusivity.cpp \
      histConfigs.cpp \
      formatLabelName.cpp \
      create_directories.cpp \
      kinematic_cuts.cpp \
      plot_pi0_mass.cpp \
      bin_boundaries.cpp \
      plot_dvcs_data_mc_comparison.cpp \
      all_bin_data.cpp \
      plot_yield_comparison.cpp \
      plot_unfolding.cpp \
      bin_helpers.cpp \
      write_csv.cpp \
      calculate_contamination.cpp  \
      plot_comparison.cpp 

# Object files
OBJ = $(SRC:.cpp=.o)

# Main target
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LFLAGS)

# Generic rule for compiling source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< $(ROOTLIBS)

# Clean target
clean:
	rm -f $(OBJ) $(TARGET)