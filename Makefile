.PHONY: all sts-online setup-cmake clean doc test debug release

BUILD = _build/release

EXE = sts-online

all: release

release: BUILD=_build/release
release: BUILD_TYPE=Release
release: setup-cmake
release: $(EXE)

debug: BUILD=_build/debug
debug: BUILD_TYPE = Debug
debug: setup-cmake
debug: $(EXE)

$(EXE): setup-cmake
	+make -C$(BUILD) $@

test: BUILD=_build/debug
test: BUILD_TYPE = Debug
test: BUILD_TESTING = -DBUILD_TESTING=ON
test: setup-cmake
	+make -C$(BUILD) run-tests
	$(BUILD)/test/run-tests

setup-cmake: CMakeLists.txt
	mkdir -p $(BUILD)
	cd $(BUILD) && cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ${BUILD_TESTING} ../..

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
