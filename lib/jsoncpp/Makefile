CPP_FILES = $(addprefix src/lib_json/, json_reader.cpp json_value.cpp json_writer.cpp)
OBJ = $(CPP_FILES:.cpp=.o)

CPPFLAGS = -Iinclude $(CFLAGS)
LIB := libjsoncpp.a

.PHONY: all clean

all: $(LIB)

clean:
	rm -f $(LIB) $(OBJ)

$(LIB): $(OBJ)
	ar rvs $@ $^


