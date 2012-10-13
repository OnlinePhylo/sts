LIBS = -lbpp-core -lbpp-seq -lbpp-phyl -lsmctc -lgsl -lgslcblas

CC ?= gcc
CXX ?= g++

SMCTC ?= ../lib/smctc
CFLAGS := $(CCFLAGS) -I$(SMCTC)/include -L$(SMCTC)/lib -g -std=c++0x
LFLAGS := -I$(SMCTC)/include -L$(SMCTC)/lib $(LFLAGS)

OBJ = bpp_shim.o hmsbeagle.o phylofunc.o phylomain.o

all: phylo

%.o: %.cc %.hh $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) $(LFLAGS) $(LIBS) `pkg-config --libs --cflags hmsbeagle-1`

phylo: $(OBJ)
	$(CXX) -o $@ $^ -g $(CFLAGS) $(LFLAGS) -o phylo $(LIBS) `pkg-config --libs --cflags hmsbeagle-1`

clean:
	rm -f *.o phylo
