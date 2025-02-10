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
INCLUDES=-I./ -I/home/simon/miniforge3/envs/ROOT/include/eigen3 -I/home/simon/miniforge3/envs/ROOT/include/ -I/home/simon/standaloneEEC

default: bin/res4_standalone_benchmark bin/test_res4_standalone bin/res4_standalone_fastEEC bin/res4_standalone_multi_array bin/res4_standalone_multi_array_precomputed bin/res4_standalone_vector bin/res4_standalone_vector_precomputed bin/res4_standalone_transfer_benchmark bin/res4_standalone_transfer_multi_array bin/res4_standalone_transfer_vector bin/res4_standalone_transfer_fastEEC bin/res4_standalone_validate

clean: 
	rm -f lib/*
	rm -f bin/*

clean_output:
	rm -f callgrind.out.*
	rm -f massif.out.*

STANDALONEOBJ=$(wildcard SRothman/EECs/src/standalones/*.cc)
STANDALONEHEADERS=$(wildcard SRothman/EECs/src/standalones/*.h)
lib/EEC_standalone.so: $(STANDALONEOBJ) $(STANDALONEHEADERS)
	g++ -shared -fPIC -o $@ $(STANDALONEOBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

EECOLDOBJ=$(wildcard SRothman/EECs/src/*.cc)
EECOLDHEADERS=$(wildcard SRothman/EECs/src/*.h)
lib/EEC_old.so: $(EECOLDOBJ) $(EECOLDHEADERS)
	g++ -shared -fPIC -o $@ $(EECOLDOBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

TOOLSOBJ=$(wildcard SRothman/SimonTools/src/*.cc)
TOOLSHEADERS=$(wildcard SRothman/SimonTools/src/*.h)
lib/SimonTools.so: $(TOOLSOBJ) $(TOOLSHEADERS)
	g++ -shared -fPIC -o $@ $(TOOLSOBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

MATCHINGOBJ=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
lib/Matching.so: $(MATCHINGOBJ) $(MATCHINGHEADERS)
	g++ -shared -fPIC -o $@ $(MATCHINGOBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

MATCHINGV2OBJ=$(wildcard SRothman/Matching/src/v2/*.cc)
MATCHINGV2HEADERS=$(wildcard SRothman/Matching/src/v2/*.h)
lib/MatchingV2.so: $(MATCHINGV2OBJ) $(MATCHINGV2HEADERS)
	g++ -shared -fPIC -o $@ $(MATCHINGV2OBJ) $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

ALLLIBS=lib/EEC_standalone.so lib/EEC_old.so lib/SimonTools.so lib/Matching.so lib/MatchingV2.so
libs: $(ALLLIBS)

bin/test_res4_standalone: test_res4_standalone.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_benchmark: res4_standalone_benchmark.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_fastEEC: res4_standalone_fastEEC.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_multi_array: res4_standalone_multi_array.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_multi_array_precomputed: res4_standalone_multi_array_precomputed.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_vector: res4_standalone_vector.cc  $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_vector_precomputed: res4_standalone_vector_precomputed.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_transfer_benchmark: res4_standalone_transfer_benchmark.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_transfer_multi_array: res4_standalone_transfer_multi_array.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_transfer_vector: res4_standalone_transfer_vector.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_transfer_fastEEC: res4_standalone_transfer_fastEEC.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)

bin/res4_standalone_validate: res4_standalone_validate.cc $(ALLLIBS)
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS) $(LIBS)
