MKDIR=mkdir -p
OS=$(shell uname)
CLR=clear && clear && clear

CC=g++
CCFLAGS=-std=c++11 -O3 -fopenmp
CU=nvcc
CUFLAGS=-std=c++11 -O3 -arch=sm_62 -rdc=true -DCUDA_ON --compiler-options -std=c++11

AR=ar rcs

CCSRC=./AsyPlot.cc
CUSRC=./AsyPlot.cc
CCOBJ=$(subst ./,./obj/cc/,$(patsubst %.cc,%.o, $(CCSRC)))
CUOBJ=$(subst ./,./obj/cu/,$(patsubst %.cc,%.o, $(CUSRC)))

TESTCCSRC=main.cc
TESTCUSRC=main.cu
TESTCCEXEC=mainCC.out
TESTCUEXEC=mainCU.out

INCLUDECC=srcCC/include
INCLUDECU=srcCU/include
INCCC=-I./inc
INCCU=-I./inc
LIBSCC=-L./lib -lMPMDCC
LIBSCU=-L./lib -lMPMDCU
MORE_FLAGS_OBJ_CC=
MORE_FLAGS_OBJ_CU=
MORE_FLAGS_BIN_CC=
MORE_FLAGS_BIN_CU=

obj/cc/%.o: ./%.cc
	$(CC) $(CCFLAGS) $(MORE_FLAGS_OBJ_CC) $(INCCC) -c -o $@ $< $(LIBSCC)

obj/cu/%.o: AsyPlot.cc
	$(CU) $(MORE_FLAGS_OBJ_CU) $(CUFLAGS) $(INCCU) -c -o $@ $< $(LIBSCU)

mainCC.out: main.cc
	$(CC) $(CCFLAGS) $(MORE_FLAGS_BIN_CC) $(INCCC) -o $@ $(CCOBJ) $< $(LIBSCC)

# exec/%.cu.o: exec/%.cu
# 	$(CU) $(MORE_FLAGS_OBJ_CU) $(CUFLAGS) -c -o $@ $< $(LIBSCU)

mainCU.out: main.cu
	$(CU) $(CUFLAGS) $(MORE_FLAGS_BIN_CU) $(INCCU) -o $@ $(CUOBJ) $< $(LIBSCU)


all: echo CPU GPU 

CPU: obj/ $(CCOBJ) $(TESTCCEXEC)

GPU: obj/ $(CUOBJ) $(TESTCUEXEC)

echo:
	@echo "CCSRC: " $(CCSRC)
	@echo "CUSRC: " $(CUSRC)
	@echo "CCOBJ: " $(CCOBJ)
	@echo "CUOBJ: " $(CUOBJ)
	@echo "TESTCCSRC: " $(TESTCCSRC)
	@echo "TESTCUSRC: " $(TESTCUSRC)
	@echo "TESTCCEXEC: " $(TESTCCEXEC)
	@echo "TESTCUEXEC: " $(TESTCUEXEC)

clean_obj:
	rm -rf obj/

clean_bin:
	rm -f mainCC.out mainCU.out

clean: clean_bin clean_obj

obj/:
	$(MKDIR) ./obj/cc
	$(MKDIR) ./obj/cu

