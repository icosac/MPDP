MKDIR=mkdir -p
OS=$(shell uname)
CLR=clear && clear && clear

CC=g++
CCFLAGS=-std=c++11 -Wall -O3
NV=nvcc
NVFLAGS=-std=c++11 -arch=sm_75 -rdc=true --default-stream per-thread

AR=ar rcs

CCSRC=$(wildcard src/*.cc)
NVSRC=$(wildcard src/*.cu)
CCOBJ=$(subst src/,src/obj/cc/,$(patsubst %.cc,%.o, $(CCSRC)))
NVOBJ=$(subst src/,src/obj/cu/,$(patsubst %.cu,%.o, $(NVSRC)))

TESTCCSRC=$(wildcard exec/*.cc)
TESTNVSRC=$(wildcard exec/*.cu)
TESTCCEXEC=$(subst exec/,bin/cc/,$(patsubst %.cc,%.out, $(TESTCCSRC)))
TESTNVEXEC=$(subst exec/,bin/cu/,$(patsubst %.cu,%.out, $(TESTNVSRC)))

INCLUDE=src/include
INC=-I./lib/include
LIBS=-L./lib $(INC) -lClothoids
LIB=libClothoids.a #LIB_DUBINS
MORE_FLAGS=

src/obj/cc/%.o: src/%.cc
	$(CC) $(CCFLAGS) $(MORE_FLAGS) -c -o $@ $< $(LIBS)

src/obj/cu/%.o: src/%.cu
	$(NV) $(NVFLAGS) $(MORE_FLAGS) -c -o $@ $< $(LIBS)

bin/cc/%.out: exec/%.cc
	$(CC) $(CCFLAGS) $(MORE_FLAGS) -o $@ $< $(LIBS)

bin/cu/%.out: exec/%.cu
	$(NV) $(NVFLAGS) $(MORE_FLAGS) -o $@ $< $(LIBS)

all: echo lib $(TESTCCEXEC) $(TESTNVEXEC)

echo:
	@echo "CCSRC: " $(CCSRC)
	@echo "NVSRC: " $(NVSRC)
	@echo "CCOBJ: " $(CCOBJ)
	@echo "NVOBJ: " $(NVOBJ)
	@echo "TESTCCSRC: " $(TESTCCSRC)
	@echo "TESTNVSRC: " $(TESTNVSRC)
	@echo "TESTCCEXEC: " $(TESTCCEXEC)
	@echo "TESTNVEXEC: " $(TESTNVEXEC)

lib: lib/$(LIB)

mvlib:
	@rm -rf lib/include
	$(MKDIR) lib
	$(MKDIR) lib/include
	cp -f $(INCLUDE)/*.hh lib/include
	cp -f $(INCLUDE)/*.tt lib/include

lib/$(LIB): mvlib obj/ bin/ $(CCOBJ) #TODO add CUDA support
	$(AR) lib/$(LIB) $(CCOBJ)

clean_lib:
	rm -rf lib/

clean_obj:
	rm -rf src/obj/

clean_bin:
	rm -rf bin

clean: clean_lib clean_bin clean_obj

run:
	bin/cc/main.out

bin/:
	$(MKDIR) bin/cc
	$(MKDIR) bin/cu
obj/:
	$(MKDIR) src/obj/cc
	$(MKDIR) src/obj/cu

