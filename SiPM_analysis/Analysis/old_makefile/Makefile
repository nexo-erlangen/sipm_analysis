
CXX = g++ 
CXXFLAGS = -Wall -g -O3 -std=c++11 -lSpectrum $(shell root-config --glibs) 
CXXFLAGS2 = -Wall -g -std=c++11 -I $ROOTSYS $(shell root-config --cflags) 
PROG = ./Debug/Analysis
OBJS = ./Debug/main.o ./Debug/DATA.o ./Debug/MAININPUT.o ./Debug/functions.o ./Debug/CFILE.o ./Debug/FOLDER.o ./Debug/WAVE.o

$(PROG) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(PROG) $(OBJS)
./Debug/main.o : ./main.cpp
	$(CXX) $(CXXFLAGS2) -c ./main.cpp -o ./Debug/main.o
./Debug/functions.o : ./functions.h ./functions.cpp
	$(CXX) $(CXXFLAGS2) -c ./functions.cpp -o ./Debug/functions.o
./Debug/WAVE.o : ./WAVE.h ./WAVE.cpp
	$(CXX) $(CXXFLAGS2) -c ./WAVE.cpp -o ./Debug/WAVE.o
./Debug/FOLDER.o : ./FOLDER.h ./FOLDER.cpp
	$(CXX) $(CXXFLAGS2) -c ./FOLDER.cpp -o ./Debug/FOLDER.o
./Debug/CFILE.o : ./CFILE.h ./CFILE.cpp
	$(CXX) $(CXXFLAGS2) -c ./CFILE.cpp -o ./Debug/CFILE.o
./Debug/MAININPUT.o : ./MAININPUT.h ./MAININPUT.cpp
	$(CXX) $(CXXFLAGS2) -c ./MAININPUT.cpp -o ./Debug/MAININPUT.o
./Debug/DATA.o : ./DATA.h ./DATA.cpp
	$(CXX) $(CXXFLAGS2) -c ./DATA.cpp -o ./Debug/DATA.o