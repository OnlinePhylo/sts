PHYLOC = bpp_shim.cc phylomain.cc phylofunc.cc hmsbeagle.cc
PHYLOH = bpp_shim.hh phylofunc.hh hmsbeagle.hh

CC ?= gcc
CXX ?= g++

SMCTC ?= ../lib/smctc
CCFLAGS := $(CCFLAGS) -I$(SMCTC)/include -L$(SMCTC)/lib -g
LFLAGS := -I$(SMCTC)/include -L$(SMCTC)/lib $(LFLAGS)

all: phylo

phylo: phylomain.cc phylofunc.cc phylofunc.hh hmsbeagle.cc hmsbeagle.hh bpp_shim.cc bpp_shim.hh
	$(CXX) -std=c++0x -O0 -g $(CCFLAGS) $(PHYLOC) -lbpp-core -lbpp-seq -lbpp-phyl -lsmctc -lgsl -lgslcblas -o phylo `pkg-config --libs --cflags hmsbeagle-1`

clean:
	rm -f *.o phylo
