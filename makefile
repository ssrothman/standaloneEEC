MATCHINGSRC=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
MATCHINGOBJ=$(MATCHINGSRC:.cc=.o)

OPTIMIZE=-O3 -mtune=native -flto
DEBUG=-ggdb3 -g
CXXFLAGS=--std=c++17 -pipe -DNDEBUG $(OPTIMIZE) $(DEBUG)
WARNINGFLAGS=-pedantic-errors -Wall -W -Wundef -Wpointer-arith -Wcast-align -Wsign-compare -Wsign-compare -Wmissing-noreturn  -Wmissing-format-attribute -Wunreachable-code -Wdisabled-optimization -Werror 

#INCLUDES=-I./ -I/work/submit/srothman/miniforge3/envs/uproot/include
#INCLUDES=-I./ -I/home/simon/miniforge3/include -I/home/simon/standaloneEEC
INCLUDES=-I./ -I/home/simon/miniforge3/envs/ROOT/include/eigen3 -I/home/simon/miniforge3/envs/ROOT/include/

default: test_res4_standalone res4_standalone_benchmark res4_standalone_fastEEC res4_standalone_multi_array res4_standalone_multi_array_precomputed res4_standalone_vector res4_standalone_vector_precomputed res4_standalone_transfer_benchmark res4_standalone_transfer_multi_array res4_standalone_transfer_vector res4_standalone_transfer_fastEEC

clean: 
	rm -f *.o
	rm -f *.so
	rm -f main
	rm -f run_tests
	rm -f test_res4_standalone
	rm -f res4_standalone_benchmark

STANDALONEOBJ=$(wildcard SRothman/EECs/src/standalones/*.cc)
EEC_standalone.so: $(STANDALONEOBJ)
	g++ -shared -fPIC -o $@ $^ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

EECOLDOBJ=$(wildcard SRothman/EECs/src/*.cc)
EEC_old.so: $(EECOLDOBJ)
	g++ -shared -fPIC -o $@ $^ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

test_res4_standalone: test_res4_standalone.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_benchmark: res4_standalone_benchmark.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_fastEEC: res4_standalone_fastEEC.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_multi_array: res4_standalone_multi_array.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_multi_array_precomputed: res4_standalone_multi_array_precomputed.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_vector: res4_standalone_vector.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_vector_precomputed: res4_standalone_vector_precomputed.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_transfer_benchmark: res4_standalone_transfer_benchmark.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_transfer_multi_array: res4_standalone_transfer_multi_array.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_transfer_vector: res4_standalone_transfer_vector.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)

res4_standalone_transfer_fastEEC: res4_standalone_transfer_fastEEC.cc EEC_standalone.so EEC_old.so
	g++ $^ -o $@ $(INCLUDES) $(CXXFLAGS) $(WARNINGFLAGS)







SRothman/Matching/src/%.o: SRothman/Matching/src/%.cc $(MATCHINGHEADERS)
	g++ -c -fPIC -o $@ $< $(INCLUDES) $(CXXFLAGS)

libmatching.so: $(MATCHINGOBJ)
	g++ -shared -o $@ $(MATCHINGOBJ) $(CXXFLAGS)

testRes4CA: res4CA.o testRes4CA.o
	g++ res4CA.o testRes4CA.o -o testRes4CA

testRes4CA.o: SRothman/EECs/src/testCA.cc SRothman/EECs/src/res4CA.h
	g++ -o $@ $< $(INCLUDES) $(CXXFLAGS) -c 

res4CA.o : SRothman/EECs/src/res4CA.cc SRothman/EECs/src/res4CA.h
	g++ -o $@ $< $(INCLUDES) $(CXXFLAGS) -c 

EEC.o: SRothman/EECs/src/theOnlyHeader.cc SRothman/EECs/src/*.h
	g++ -o $@ $< -larmadillo -llpack -lblas $(INCLUDES) $(CXXFLAGS) -c 

ToyShower.o : SRothman/SimonTools/src/ToyShowerer.cc SRothman/SimonTools/src/ToyShowerer.h
	g++ -o $@ $< -larmadillo -llpack -lblas $(INCLUDES) $(CXXFLAGS) -c 

main.o: main.cc cov.h SRothman/SimonTools/src/*.h SRothman/EECs/src/theOnlyHeader.h SRothman/EECs/src/fastStructs.h
	g++ -o $@ $< -larmadillo -llapack -lblas $(INCLUDES) $(CXXFLAGS) -c 

run_tests.o: run_tests.cc test_proj.h test_res3.h test_res4.h example_jet.h 
	g++ -o $@ $< -larmadillo -llapack -lblas $(INCLUDES) $(CXXFLAGS) -c 

run_tests: run_tests.o EEC.o ToyShower.o
	g++ $(root-config --glibs --cflags --libs) run_tests.o EEC.o ToyShower.o -o run_tests -lMathMore -lGenVector

main: main.o EEC.o
	g++ main.o EEC.o -o main
