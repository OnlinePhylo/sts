VERSION := $(shell git describe --always)

LIBS = -lbpp-core -lbpp-seq -lbpp-phyl -lsmctc -lgsl -lgslcblas

CC ?= gcc
CXX ?= g++

BEAGLE_FLAGS := $(shell pkg-config --libs --cflags hmsbeagle-1)
SMCTC ?= ../lib/smctc
BEAGLE_CFLAGS := $(shell pkg-config --cflags hmsbeagle-1)
BEAGLE_LFLAGS := $(shell pkg-config --libs hmsbeagle-1)
CFLAGS := $(CCFLAGS) -I$(SMCTC)/include -g -std=c++0x $(BEAGLE_CFLAGS) $ -I. $(CFLAGS)
LFLAGS := -L$(SMCTC)/lib --std=c++0x $(BEAGLE_LFLAGS) $(LFLAGS)

OBJ = bpp_shim.o hmsbeagle.o phylofunc.o phylomoves.o

all: sts

%.o: %.cc %.hh $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

sts: $(OBJ) phylomain.o
	$(CXX) -o $@ $^ -g $(LFLAGS) $(LIBS)

test: $(OBJ) test.o
	$(CXX) -o $@ $^ -g $(LFLAGS) $(LIBS)
	./$@

clean:
	rm -f *.o sts
