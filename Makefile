MKDIR=mkdir -p
OS=$(shell uname)
CLR=clear && clear && clear

CC=g++
CCFLAGS=-std=c++11 -O3
CU=nvcc
CUFLAGS=-std=c++11 -O3 -arch=sm_62 -rdc=true -DCUDA_ON --compiler-options -std=c++11

AR=ar rcs

CCSRC=$(wildcard srcCC/*.cc)
CUSRC=$(wildcard srcCU/*.cu)
CCOBJ=$(subst srcCC/,srcCC/obj/cc/,$(patsubst %.cc,%.o, $(CCSRC)))
CUOBJ=$(subst srcCU/,srcCU/obj/cu/,$(patsubst %.cu,%.o, $(CUSRC)))

TESTCCSRC=$(filter-out exec/generate.cc, $(wildcard exec/*.cc))
TESTCUSRC=$(wildcard exec/*.cu)
TESTCCEXEC=$(subst exec/,bin/cc/,$(patsubst %.cc,%.out, $(TESTCCSRC)))
TESTCUEXEC=$(subst exec/,bin/cu/,$(patsubst %.cu,%.out, $(TESTCUSRC)))

INCLUDECC=srcCC/include
INCLUDECU=srcCU/include
INCCC=-I./lib/includeCC
INCCU=-I./lib/includeCU
LIBSCC=-L./lib -lMPMDCC
LIBSCU=-L./lib -lMPMDCU
LIBCC=libMPMDCC.a 
LIBCU=libMPMDCU.a 
MORE_FLAGS_OBJ_CC=
MORE_FLAGS_OBJ_CU=
MORE_FLAGS_BIN_CC=
MORE_FLAGS_BIN_CU=

srcCC/obj/cc/%.o: srcCC/%.cc
	$(CC) $(CCFLAGS) $(MORE_FLAGS_OBJ_CC) $(INCCC) -c -o $@ $<

srcCU/obj/cu/%.o: srcCU/%.cu
	$(CU) $(MORE_FLAGS_OBJ_CU) $(CUFLAGS) $(INCCU) -c -o $@ $<

bin/cc/%.out: exec/%.cc
	$(CC) $(CCFLAGS) $(MORE_FLAGS_BIN_CC) $(INCCC) -o $@ $< $(LIBSCC)

# exec/%.cu.o: exec/%.cu
# 	$(CU) $(MORE_FLAGS_OBJ_CU) $(CUFLAGS) -c -o $@ $< $(LIBSCU)

bin/cu/%.out: exec/%.cu
	$(CU) $(CUFLAGS) $(MORE_FLAGS_BIN_CU) $(INCCU) -o $@ $< $(LIBSCU)


all: echo CPU GPU 

norProj: lib/$(LIBCC) lib/$(LIBCU)
	rm -rf norProject/lib
	rm -rf norProject/inc
	$(MKDIR) norProject/lib
	$(MKDIR) norProject/inc
	cp -f include/*.hh norProject/inc
	cp -f $(INCLUDECC)/*.hh norProject/inc
	cp -f $(INCLUDECU)/*.cuh norProject/inc
	cp -f lib/$(LIBCC) norProject/lib/$(LIBCC)
	cp -f lib/$(LIBCU) norProject/lib/$(LIBCU)

CPU: lib/$(LIBCC) $(TESTCCEXEC)

GPU: lib/$(LIBCU) $(TESTCUEXEC)

echo:
	@echo "CCSRC: " $(CCSRC)
	@echo "CUSRC: " $(CUSRC)
	@echo "CCOBJ: " $(CCOBJ)
	@echo "CUOBJ: " $(CUOBJ)
	@echo "TESTCCSRC: " $(TESTCCSRC)
	@echo "TESTCUSRC: " $(TESTCUSRC)
	@echo "TESTCCEXEC: " $(TESTCCEXEC)
	@echo "TESTCUEXEC: " $(TESTCUEXEC)

lib: lib/$(LIBCC) lib/$(LIBCU)

mvlibCC:
	@rm -rf lib/include
	$(MKDIR) lib/includeCC
	cp -f include/*.hh lib/includeCC
	cp -f $(INCLUDECC)/*.hh lib/includeCC

mvlibCU:
	@rm -rf lib/includeCU
	$(MKDIR) lib/includeCU
	cp -f include/*.hh lib/includeCU
	cp -f $(INCLUDECU)/*.cuh lib/includeCU

lib/$(LIBCC): mvlibCC obj/ bin/ $(CCOBJ) 
	$(AR) lib/$(LIBCC) $(CCOBJ)

lib/$(LIBCU): mvlibCU obj/ bin/ $(CUOBJ) 
	$(AR) lib/$(LIBCU) $(CUOBJ)

clean_lib:
	rm -rf lib/

clean_obj:
	rm -rf srcCC/obj/
	rm -rf srcCU/obj/

clean_bin:
	rm -rf bin

clean: clean_lib clean_bin clean_obj

run:
	bin/cc/main.out

bin/:
	$(MKDIR) bin/cc
	$(MKDIR) bin/cu
obj/:
	$(MKDIR) srcCC/obj/cc
	$(MKDIR) srcCU/obj/cu

