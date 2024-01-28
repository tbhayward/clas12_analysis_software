# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -g -Wno-unused-but-set-variable
ROOTLIBS = $(shell root-config --cflags --libs)
LFLAGS = -L/site/12gev_phys/2.4/Linux_CentOS7.9.2009-gcc9.2.0/root/6.20.04/lib $(ROOTLIBS) -lMinuit

# Target executable
TARGET = asymmetry_extraction

# Source files
SRC = asymmetry_extraction.cpp load_bins_from_csv.cpp load_run_info_from_csv.cpp dilution_factor.cpp asymmetry_fits.cpp KinematicCuts.cpp inclusiveKinematicCuts.cpp singleHadronKinematicCuts.cpp b2bDihadronKinematicCuts.cpp dihadronKinematicCuts.cpp formatLabelName.cpp readChi2Fits.cpp histConfigs.cpp charge_accumulation.cpp plot_data.cpp modifyTree.cpp

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
