PHYLOC = phylomain.cc phylofunc.cc
PHYLOH = phylofunc.hh

CC ?= gcc
CXX ?= g++

SMCTC ?= ../lib/smctc
CCFLAGS := $(CCFLAGS) -I$(SMCTC)/include -L$(SMCTC)/lib -g
LFLAGS := -I$(SMCTC)/include -L$(SMCTC)/lib $(LFLAGS)

all: phylo

phylo: phylomain.cc phylofunc.cc phylofunc.hh hmsbeagle.hh
	$(CXX) -std=c++0x -O0 -g $(CCFLAGS) $(PHYLOC) -lbpp-core -lbpp-seq -lbpp-phyl -lsmctc -lgsl -lgslcblas -o phylo `pkg-config --libs --cflags hmsbeagle-1`

clean:
	rm -f *.o phylo
