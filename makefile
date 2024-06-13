MATCHINGSRC=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
MATCHINGOBJ=$(MATCHINGSRC:.cc=.o)

CXXFLAGS=-O3 -std=c++17 -pipe -pedantic-errors -Wall -W -Wundef -Wpointer-arith -Wcast-align -Wsign-compare -Wsign-compare -Wmissing-noreturn  -Wmissing-format-attribute -Wunreachable-code -Wdisabled-optimization -Werror

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
	g++ -o $@ $< -lblas -llapack $(INCLUDES) $(CXXFLAGS)
