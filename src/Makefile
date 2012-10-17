VERSION := $(shell git describe --always)

LIBS = -lbpp-core -lbpp-seq -lbpp-phyl -lsmctc -lgsl -lgslcblas

CC ?= gcc
CXX ?= g++

BEAGLE_FLAGS := $(shell pkg-config --libs --cflags hmsbeagle-1)
SMCTC ?= ../lib/smctc
CFLAGS := $(CCFLAGS) -I$(SMCTC)/include -L$(SMCTC)/lib -I. -DSTS_VERSION=$(VERSION) -g -std=c++0x
LFLAGS := -I$(SMCTC)/include -L$(SMCTC)/lib $(LFLAGS)

OBJ = bpp_shim.o hmsbeagle.o phylofunc.o phylomain.o

all: sts

%.o: %.cc %.hh $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LFLAGS) $(LIBS) $(BEAGLE_FLAGS)

sts: $(OBJ)
	$(CXX) -o $@ $^ -g $(CFLAGS) $(LFLAGS) $(LIBS) $(BEAGLE_FLAGS)

clean:
	rm -f *.o sts
