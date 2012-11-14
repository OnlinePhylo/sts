cmake_minimum_required (VERSION 2.8.1)
cmake_policy(SET CMP0015 NEW)

project (jsoncpp CXX)

if(NOT DEFINED CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "-g -std=c++0x")
endif()

include_directories(include)

set(JSONCPP_SOURCE_FILES
  src/lib_json/json_reader.cpp
  src/lib_json/json_value.cpp
  src/lib_json/json_writer.cpp)

ADD_LIBRARY(jsoncpp STATIC
            ${JSONCPP_SOURCE_FILES})
