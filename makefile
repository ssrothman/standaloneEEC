MATCHINGSRC=$(wildcard SRothman/Matching/src/*.cc)
MATCHINGHEADERS=$(wildcard SRothman/Matching/src/*.h)
MATCHINGOBJ=$(MATCHINGSRC:.cc=.o)

CXXFLAGS=-O3 -std=c++17 -pipe -pedantic-errors -Wall -W -Wundef -Wpointer-arith -Wcast-align -Wsign-compare -Wsign-compare -Wmissing-noreturn  -Wmissing-format-attribute -Wunreachable-code -Wdisabled-optimization -Werror -ggdb3 -DNDEBUG

#INCLUDES=-I./ -I/work/submit/srothman/miniforge3/envs/uproot/include
#INCLUDES=-I./ -I/home/simon/miniforge3/include -I/home/simon/standaloneEEC
INCLUDES=-I./ -I/home/simon/miniforge3/envs/ROOT/include/eigen3 -I/home/simon/miniforge3/envs/ROOT/include/

default: run_tests

clean: 
	rm -f *.o
	rm -f main
	rm -f run_tests

SRothman/Matching/src/%.o: SRothman/Matching/src/%.cc $(MATCHINGHEADERS)
	g++ -c -fPIC -o $@ $< $(INCLUDES) $(CXXFLAGS)

libmatching.so: $(MATCHINGOBJ)
	g++ -shared -o $@ $(MATCHINGOBJ) $(CXXFLAGS)

EEC.o: SRothman/EECs/src/theOnlyHeader.cc SRothman/EECs/src/*.h
	g++ -o $@ $< -larmadillo -llpack -lblas $(INCLUDES) $(CXXFLAGS) -c -DARMA_DONT_USE_WRAPPER

main.o: main.cc cov.h SRothman/SimonTools/src/*.h SRothman/EECs/src/theOnlyHeader.h SRothman/EECs/src/fastStructs.h
	g++ -o $@ $< -larmadillo -llapack -lblas $(INCLUDES) $(CXXFLAGS) -c -DARMA_DONT_USE_WRAPPER

run_tests.o: run_tests.cc test_proj.h test_res3.h test_res4.h example_jet.h
	g++ -o $@ $< -larmadillo -llapack -lblas $(INCLUDES) $(CXXFLAGS) -c -DARMA_DONT_USE_WRAPPER

run_tests: run_tests.o EEC.o
	g++ run_tests.o EEC.o -o run_tests -DARMA_DONT_USE_WRAPPER

main: main.o EEC.o
	g++ main.o EEC.o -o main
