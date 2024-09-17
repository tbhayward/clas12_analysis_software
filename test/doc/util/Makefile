CXX = g++ -std=c++11
FLAGS = -g -Wno-deprecated -fPIC -m64 -fno-inline -Wno-write-strings

# QADB and rapidjson
DEPS += -I$(QADB)/srcC/rapidjson/include -I$(QADB)/srcC/include

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
