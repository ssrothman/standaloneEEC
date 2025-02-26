MATCHINGSRC=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
MATCHINGOBJ=$(MATCHINGSRC:.cc=.o)

OPTIMIZE=-O3 -mtune=native -flto=-1
DEBUG=-ggdb3 -g
CXXFLAGS=--std=c++17 -pipe -DNDEBUG $(OPTIMIZE) $(DEBUG)
WARNINGFLAGS=-pedantic-errors -Wall -W -Wundef -Wpointer-arith -Wcast-align -Wsign-compare -Wsign-compare -Wmissing-noreturn  -Wmissing-format-attribute -Wunreachable-code -Wdisabled-optimization -Werror -Wno-suggest-attribute=format
LIBS=$(root-config --glibs --cflags --libs) -lMathMore -lGenVector -lMinuit2

#INCLUDES=-I./ -I/work/submit/srothman/miniforge3/envs/uproot/include
#INCLUDES=-I./ -I/home/simon/miniforge3/include -I/home/simon/standaloneEEC
INCLUDES=-I./ -I/home/simon/miniforge3/envs/ROOT/include/eigen3 -I/home/simon/miniforge3/envs/ROOT/include/ -I/home/simon/standaloneEEC -I/home/simon/miniforge3/envs/ROOT/include/eigen3

default: bin/res4_benchmark bin/res3_benchmark bin/CAres3_benchmark 

clean: 
	rm -f lib/*
	rm -f bin/*

clean_output:
	rm -f callgrind.out.*
	rm -f massif.out.*

EECOBJ=$(wildcard SRothman/EECs/src/*.cc)
EECHEADERS=$(wildcard SRothman/EECs/src/*.h)
lib/EEC.so: $(EECOBJ) $(EECHEADERS)
	g++ -shared -fPIC -o $@ $(EECOBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

TOOLSOBJ=$(wildcard SRothman/SimonTools/src/*.cc)
TOOLSHEADERS=$(wildcard SRothman/SimonTools/src/*.h)
lib/SimonTools.so: $(TOOLSOBJ) $(TOOLSHEADERS)
	g++ -shared -fPIC -o $@ $(TOOLSOBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

MATCHINGOBJ=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
lib/Matching.so: $(MATCHINGOBJ) $(MATCHINGHEADERS)
	g++ -shared -fPIC -o $@ $(MATCHINGOBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

ALLLIBS=lib/EEC.so lib/SimonTools.so lib/Matching.so
libs: $(ALLLIBS)

bin/res4_benchmark: res4_benchmark.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res3_benchmark: res3_benchmark.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/CAres3_benchmark: CAres3_benchmark.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/test_matching_v2: test_matching_v2.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)
