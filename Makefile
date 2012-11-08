.PHONY: all smctc_lib style clean continuous doc jsoncpp_lib

all: smctc_lib jsoncpp_lib
	$(MAKE) -Csrc all

doc:
	doxygen Doxyfile

smctc_lib:
	+make -Clib/smctc libraries

jsoncpp_lib:
	+make -Clib/jsoncpp

style:
	astyle  -A3 \
		--pad-oper \
		--unpad-paren \
		--keep-one-line-blocks \
		--keep-one-line-statements \
		--suffix=none \
		--formatted \
		--lineend=linux \
		`find src include/sts -regextype posix-extended -regex ".*\.(cc|hh|cpp|hpp)$$"`

clean:
	$(MAKE) -Csrc clean
	$(MAKE) -Clib/smctc clean
	$(MAKE) -Clib/jsoncpp clean

test: smctc_lib jsoncpp_lib
	$(MAKE) -Csrc run-test

continuous:
	while :; do inotifywait -q -e modify -r src include @src/sts; $(MAKE) all; done
