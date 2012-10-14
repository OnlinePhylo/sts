.PHONY: all smctc_lib style clean continuous doc

all: smctc_lib
	$(MAKE) -Csrc all

doc:
	doxygen Doxyfile

smctc_lib:
	$(MAKE) -Clib/smctc libraries

style:
	astyle  -A3 \
		--pad-oper \
		--unpad-paren \
		--keep-one-line-blocks \
		--keep-one-line-statements \
		--suffix=none \
		--formatted \
		--lineend=linux \
		`find src -regextype posix-extended -regex ".*\.(cc|hh|c|h)"`

clean:
	$(MAKE) -Csrc clean

continuous:
	while :; do inotifywait -q -e modify -r src @src/phylo; $(MAKE) all; done

