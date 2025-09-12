# Compiler and flags
CXX = g++
CXXFLAGS = `root-config --cflags` -O2 -std=c++17 -pthread
LDFLAGS = `root-config --libs` -pthread

# Source files
SRC = main.cpp load_trees.cpp make_dirs.cpp exclusivity_cuts.cpp load_binning_scheme.cpp bin_means.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = dvcs_analysis

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(LDFLAGS)

# Build object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJ) $(TARGET)