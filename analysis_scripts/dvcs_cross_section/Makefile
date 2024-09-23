# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -g -Wno-unused-but-set-variable -std=c++17
ROOTLIBS = $(shell root-config --cflags --libs)
LFLAGS = -L/path/to/root/lib $(ROOTLIBS) -lMinuit

# Target executable
TARGET = dvcs_cross_section_extraction

# Source files
SRC = dvcs_cross_section_extraction.cpp determine_exclusivity.cpp histConfigs.cpp formatLabelName.cpp create_directories.cpp

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