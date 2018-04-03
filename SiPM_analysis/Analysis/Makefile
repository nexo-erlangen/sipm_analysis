
CXX = g++ 
CXXFLAGS = -Wall -g -O2 -std=c++11 $(shell root-config --glibs) -lSpectrum
CXXFLAGS2 = -Wall -g -O2 -std=c++11 $(shell root-config --cflags) 
PROG = ./Analysis
SOURCES =  ./main.cpp ./DATA.cpp ./MAININPUT.cpp ./functions.cpp ./CFILE.cpp ./FOLDER.cpp ./WAVE.cpp
OBJECTS=$(SOURCES:.cpp=.o)


all: $(SOURCES) $(PROG)

$(PROG) : $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(CXXFLAGS)
.cpp.o:
	$(CXX) $(CXXFLAGS2) -c $< -o $@

clean:
	rm -rf *.o
