MATCHINGSRC=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
MATCHINGOBJ=$(MATCHINGSRC:.cc=.o)

CXXFLAGS=-O0 -std=c++17 -pipe -pedantic-errors -Wall -W -Wundef -Wpointer-arith -Wcast-align -Wsign-compare -Wsign-compare -Wmissing-noreturn  -Wmissing-format-attribute -Wunreachable-code -Wdisabled-optimization -Werror -ggdb3 -DNDEBUG

#INCLUDES=-I./ -I/work/submit/srothman/miniforge3/envs/uproot/include
INCLUDES=-I./ -I/home/simon/miniforge3/include -I/home/simon/standaloneEEC

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

EEC.o: SRothman/EECs/src/theOnlyHeader.cc SRothman/EECs/src/*.h
	g++ -o $@ $< -lblas -llapack $(INCLUDES) $(CXXFLAGS) -c

main.o: main.cc cov.h SRothman/SimonTools/src/*.h SRothman/EECs/src/theOnlyHeader.h SRothman/EECs/src/fastStructs.h
	g++ -o $@ $< -lblas -llapack $(INCLUDES) $(CXXFLAGS) -c

main: main.o EEC.o
	g++ main.o EEC.o -o main
