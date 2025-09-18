# Compiler and flags
CXX = g++
CXXFLAGS = `root-config --cflags` -O2 -std=c++17 -pthread
LDFLAGS = `root-config --libs` -pthread

# Source files
SRC = main.cpp load_trees.cpp make_dirs.cpp exclusivity_cuts.cpp load_binning_scheme.cpp bin_means.cpp total_counts.cpp pi0_contamination.cpp bsa.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = dvcs_analysis

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)