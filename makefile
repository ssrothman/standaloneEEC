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

default: all_res3 all_CAres3 all_res4 all_CAres4 all_proj

all_res3: bin/res3_benchmark bin/res3_unittest bin/res3_checkbyhand

all_CAres3: bin/CAres3_benchmark bin/CAres3_unittest bin/CAres3_checkbyhand

all_res4: bin/res4_benchmark bin/res4_unittest bin/res4_checkbyhand

all_CAres4: bin/CAres4_benchmark bin/CAres4_unittest bin/CAres4_checkbyhand

all_proj: bin/proj_unittest

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

lib/EEC_byhand.so: $(EECOBJ) $(EECHEADERS)
	g++ -shared -fPIC -o $@ $(EECOBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) -DCHECK_BY_HAND

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

ALLLIBS_BYHAND=lib/EEC_byhand.so lib/SimonTools.so lib/Matching.so
libs_byhand: $(ALLLIBS_BYHAND)

bin/res4_benchmark: res4_benchmark.cc $(ALLLIBS) testResCalculator.h
	g++ res4_benchmark.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_unittest: res4_unittest.cc $(ALLLIBS) testResCalculator.h
	g++ res4_unittest.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_checkbyhand: res4_checkbyhand.cc $(ALLLIBS_BYHAND) testResCalculator.h
	g++ res4_checkbyhand.cc $(ALLLIBS_BYHAND) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS) -DCHECK_BY_HAND

bin/CAres4_benchmark: CAres4_benchmark.cc $(ALLLIBS) testResCalculator.h
	g++ CAres4_benchmark.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/CAres4_unittest: CAres4_unittest.cc $(ALLLIBS) testResCalculator.h
	g++ CAres4_unittest.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS) 

bin/CAres4_checkbyhand: CAres4_checkbyhand.cc $(ALLLIBS_BYHAND) testResCalculator.h
	g++ CAres4_checkbyhand.cc $(ALLLIBS_BYHAND) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS) -DCHECK_BY_HAND

bin/res3_benchmark: res3_benchmark.cc $(ALLLIBS) testResCalculator.h
	g++ res3_benchmark.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res3_unittest: res3_unittest.cc $(ALLLIBS) testResCalculator.h
	g++ res3_unittest.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res3_checkbyhand: res3_checkbyhand.cc $(ALLLIBS_BYHAND) testResCalculator.h
	g++ res3_checkbyhand.cc $(ALLLIBS_BYHAND) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS) -DCHECK_BY_HAND

bin/CAres3_benchmark: CAres3_benchmark.cc $(ALLLIBS) testResCalculator.h
	g++ CAres3_benchmark.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/CAres3_unittest: CAres3_unittest.cc $(ALLLIBS) testResCalculator.h
	g++ CAres3_unittest.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/CAres3_checkbyhand: CAres3_checkbyhand.cc $(ALLLIBS_BYHAND) testResCalculator.h
	g++ CAres3_checkbyhand.cc $(ALLLIBS_BYHAND) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS) -DCHECK_BY_HAND

bin/proj_unittest: proj_unittest.cc $(ALLLIBS) testProjCalculator.h
	g++ proj_unittest.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/test_matching_v2: test_matching_v2.cc $(ALLLIBS)
	g++ test_matching_v2.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/testshower: testshower.cc $(ALLLIBS)
	g++ testshower.cc $(ALLLIBS) -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)
