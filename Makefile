.PHONY: all sts-online setup-cmake clean doc test debug release

BUILD = _build/debug

EXE = sts-online

all: debug

release: BUILD=_build/release
release: CMAKE_BUILD_TYPE=Release
release: setup-cmake
release: $(EXE)

debug: BUILD=_build/debug
debug: CMAKE_BUILD_TYPE = Debug
debug: setup-cmake
debug: $(EXE)

$(EXE): setup-cmake
	+make -C$(BUILD) $@

test: BUILD=_build/debug
test: CMAKE_BUILD_TYPE = Debug
test: setup-cmake
	+make -C$(BUILD) run-tests
	$(BUILD)/test/run-tests

setup-cmake: CMakeLists.txt
	mkdir -p $(BUILD)
	cd $(BUILD) && cmake ../..

doc:
	doxygen Doxyfile

clean:
	rm -rf $(BUILD) doc/*

style:
	astyle  -A3 \
	        --pad-oper \
	        --unpad-paren \
	        --keep-one-line-blocks \
	        --keep-one-line-statements \
	        --suffix=none \
	        --formatted \
	        --lineend=linux \
	        `find src -regextype posix-extended -regex ".*\.(cc|h|hpp)$$"`
