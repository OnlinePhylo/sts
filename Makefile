.PHONY: all smctc_lib style clean

all: smctc_lib
	$(MAKE) -Csrc all

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
		`find src -regextype posix-extended -regex ".*\.(cc|hh)"`

clean:
	$(MAKE) -Csrc clean
