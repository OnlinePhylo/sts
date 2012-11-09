.PHONY: all sts setup-cmake clean

BUILD := _build

all: sts

sts: setup-cmake
	+make -C$(BUILD) $@

test:
	+make -C$(BUILD) run-tests
	$(BUILD)/run-tests

setup-cmake:
	mkdir -p $(BUILD)
	cd $(BUILD) && cmake ..

clean:
	rm -rf $(BUILD)
