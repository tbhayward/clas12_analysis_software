# Variables
CXX = g++
ROOTLIBS = $(shell root-config --cflags --libs)
LFLAGS = -L/site/12gev_phys/2.4/Linux_CentOS7.9.2009-gcc9.2.0/root/6.20.04/lib $(ROOTLIBS) -lMinuit

# Target executable
TARGET = analysis_scripts/asymmetry_extraction/asymmetry_extraction


# Source files
SRC = BSA_rgc_fits_root_tree.cpp load_bins_from_csv.cpp

# Object files
OBJ = $(SRC:.cpp=.o)

# Main target
$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LFLAGS)

# Generic rule for compiling source files
%.o: %.cpp
	$(CXX) -c $< $(ROOTLIBS)

# Clean target
clean:
	rm -f $(OBJ) $(TARGET)
