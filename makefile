MATCHINGSRC=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
MATCHINGOBJ=$(MATCHINGSRC:.cc=.o)

CXXFLAGS=-O3 -std=c++17 -Wall -Wpedantic -Wextra

#INCLUDES=-I./ -I/work/submit/srothman/miniforge3/envs/uproot/include
INCLUDES=-I./ -I/home/simon/miniforge3/include

default: main

cleanMatching: 
	rm -f $(MATCHINGOBJ) libmatching.so

cleanMain: 
	rm -f main

clean: cleanMatching cleanMain

SRothman/Matching/src/%.o: SRothman/Matching/src/%.cc $(MATCHINGHEADERS)
	g++ -c -fPIC -o $@ $< $(INCLUDES) $(CXXFLAGS)

libmatching.so: $(MATCHINGOBJ)
	g++ -shared -o $@ $(MATCHINGOBJ) $(CXXFLAGS)

main: main.cc SRothman/EECs/src/*.h
	g++ -o $@ $< -lMinuit2 -lblas -llapack $(INCLUDES) $(CXXFLAGS)
