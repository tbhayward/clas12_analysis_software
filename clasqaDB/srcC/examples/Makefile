CXX = g++ -std=c++11
FLAGS = -g -Wno-deprecated -fPIC -m64 -fno-inline -Wno-write-strings

# QADB and rapidjson
DEPS += -I$(QADB)/srcC/rapidjson/include -I$(QADB)/srcC/include

# clas12root
DEPS += -I$(CLAS12ROOT)/hipo4 -I$(CLAS12ROOT)/Clas12Banks -I$(CLAS12ROOT)/Clas12Root
LIBS += -L$(CLAS12ROOT)/lib -lClas12Banks -lClas12Root -lHipo4 -llz4

# ROOT
DEPS += $(shell root-config --cflags)
LIBS += $(shell root-config --glibs)
#LIBS += -lMinuit -lRooFitCore -lRooFit


#--------------------------------------------

SOURCES := $(basename $(wildcard *.cpp))
EXES := $(addsuffix .exe, $(SOURCES))

#--------------------------------------------

all: 
	make exe

exe: $(EXES)

%.exe: %.o
	@echo "--- make executable $@"
	$(CXX) -o $@ $< $(LIBS)

%.o: %.cpp
	@echo "----- build $@ -----"
	$(CXX) -c $^ -o $@ $(FLAGS) $(DEPS)

clean:
	@echo "----- clean -----"
	$(RM) $(EXES)
