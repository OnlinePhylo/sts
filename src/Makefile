VERSION := $(shell git describe --always)

LIBS = -lbpp-core -lbpp-seq -lbpp-phyl -lsmctc -lgsl -lgslcblas

CC ?= gcc
CXX ?= g++

BEAGLE_FLAGS := $(shell pkg-config --libs --cflags hmsbeagle-1)
SMCTC ?= ../lib/smctc
BEAGLE_CFLAGS := $(shell pkg-config --cflags hmsbeagle-1)
BEAGLE_LFLAGS := $(shell pkg-config --libs hmsbeagle-1)
JSONCPP_CFLAGS := $(shell pkg-config --cflags jsoncpp)
JSONCPP_LFLAGS := $(shell pkg-config --libs jsoncpp)
CFLAGS := $(CCFLAGS) -I$(SMCTC)/include -g -std=c++0x $(BEAGLE_CFLAGS) $(JSONCPP_CFLAGS) -DSTS_VERSION=$(VERSION) $ -I../include $(CFLAGS)
LFLAGS := -L$(SMCTC)/lib --std=c++0x $(BEAGLE_LFLAGS) $(JSONCPP_LFLAGS) $(LFLAGS)

DEPS := $(shell find ../include -name "*.hpp")

all: sts

%.o: %.cc $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

sts: sts.o
	$(CXX) -o $@ $< -g $(LFLAGS) $(LIBS)

test: $(OBJ) test.o
	$(CXX) -o $@ $^ -g $(LFLAGS) $(LIBS)

run-test: test
	./test

clean:
	rm -f *.o sts test
