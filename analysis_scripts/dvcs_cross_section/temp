# Variables
CXX = g++
CXXFLAGS = -Wall -Wextra -g -Wno-unused-but-set-variable -std=c++17
ROOTLIBS = $(shell root-config --cflags --libs)
LFLAGS = -L/site/12gev_phys/2.4/Linux_CentOS7.9.2009-gcc9.2.0/root/6.20.04/lib $(ROOTLIBS) -lMinuit

# Target executable
TARGET = dvcs_cross_section_extraction

# Source files
SRC = dvcs_cross_section_extraction.cpp  # Add new source files as we develop them

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