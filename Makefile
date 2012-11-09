.PHONY: all sts setup-cmake clean

BUILD := _build

all: sts

sts: setup-cmake
	+make -C$(BUILD) $@

test: setup-cmake
	+make -C$(BUILD) run-tests
	$(BUILD)/run-tests

setup-cmake:
	mkdir -p $(BUILD)
	cd $(BUILD) && cmake ..

clean:
	rm -rf $(BUILD)

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
