LIBS = -lbpp-core -lbpp-seq -lbpp-phyl -lsmctc -lgsl -lgslcblas

CC ?= gcc
CXX ?= g++

SMCTC ?= ../lib/smctc
BEAGLE_CFLAGS := $(shell pkg-config --cflags hmsbeagle-1)
BEAGLE_LFLAGS := $(shell pkg-config --libs hmsbeagle-1)
CFLAGS := $(CCFLAGS) -I$(SMCTC)/include -g -std=c++0x $(BEAGLE_CFLAGS) $(CFLAGS)
LFLAGS := -L$(SMCTC)/lib -lgsl -lgslcblas -lsmctc $(BEAGLE_LFLAGS) $(LFLAGS)

OBJ = bpp_shim.o hmsbeagle.o phylofunc.o phylomain.o phylomoves.o

all: phylo

%.o: %.cc %.hh $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

phylo: $(OBJ)
	$(CXX) -o $@ $^ -g $(LFLAGS) -o phylo $(LIBS)

clean:
	rm -f *.o phylo
