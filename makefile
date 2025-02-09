MATCHINGSRC=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
MATCHINGOBJ=$(MATCHINGSRC:.cc=.o)

OPTIMIZE=-O3 -mtune=native -flto
DEBUG=-ggdb3 -g
CXXFLAGS=--std=c++17 -pipe -DNDEBUG $(OPTIMIZE) $(DEBUG)
WARNINGFLAGS=-pedantic-errors -Wall -W -Wundef -Wpointer-arith -Wcast-align -Wsign-compare -Wsign-compare -Wmissing-noreturn  -Wmissing-format-attribute -Wunreachable-code -Wdisabled-optimization -Werror 

#INCLUDES=-I./ -I/work/submit/srothman/miniforge3/envs/uproot/include
#INCLUDES=-I./ -I/home/simon/miniforge3/include -I/home/simon/standaloneEEC
INCLUDES=-I./ -I/home/simon/miniforge3/envs/ROOT/include/eigen3 -I/home/simon/miniforge3/envs/ROOT/include/ -I/home/simon/standaloneEEC

default: bin/res4_standalone_benchmark bin/test_res4_standalone bin/res4_standalone_fastEEC bin/res4_standalone_multi_array bin/res4_standalone_multi_array_precomputed bin/res4_standalone_vector bin/res4_standalone_vector_precomputed bin/res4_standalone_transfer_benchmark bin/res4_standalone_transfer_multi_array bin/res4_standalone_transfer_vector bin/res4_standalone_transfer_fastEEC

clean: 
	rm -f lib/*
	rm -f bin/*

clean_output:
	rm -f callgrind.out.*
	rm -f massif.out.*

STANDALONEOBJ=$(wildcard SRothman/EECs/src/standalones/*.cc)
lib/EEC_standalone.so: $(STANDALONEOBJ)
	g++ -shared -fPIC -o $@ $^ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

EECOLDOBJ=$(wildcard SRothman/EECs/src/*.cc)
lib/EEC_old.so: $(EECOLDOBJ)
	g++ -shared -fPIC -o $@ $^ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/test_res4_standalone: test_res4_standalone.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_benchmark: res4_standalone_benchmark.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_fastEEC: res4_standalone_fastEEC.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_multi_array: res4_standalone_multi_array.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_multi_array_precomputed: res4_standalone_multi_array_precomputed.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_vector: res4_standalone_vector.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_vector_precomputed: res4_standalone_vector_precomputed.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_transfer_benchmark: res4_standalone_transfer_benchmark.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_transfer_multi_array: res4_standalone_transfer_multi_array.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_transfer_vector: res4_standalone_transfer_vector.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

bin/res4_standalone_transfer_fastEEC: res4_standalone_transfer_fastEEC.cc lib/EEC_standalone.so lib/EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)
