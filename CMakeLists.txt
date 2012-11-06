cmake_minimum_required (VERSION 2.8.1)
cmake_policy(SET CMP0015 NEW)
enable_testing()

project (STS CXX)

set(CMAKE_CXX_FLAGS "-g -std=c++0x -Wall")

# Subdirectories
add_subdirectory(lib/smctc)
add_subdirectory(lib/jsoncpp)

include_directories(include)
include_directories(src)

# SMCTC
include_directories ("${PROJECT_SOURCE_DIR}/Project")
include_directories(lib/smctc/include)

# JSONCPP
include_directories(lib/jsoncpp/include)

# BEAGLE
find_package(PkgConfig REQUIRED)
pkg_check_modules(HMS_BEAGLE hmsbeagle-1 REQUIRED)
include_directories(${HMS_BEAGLE_INCLUDE_DIRS})
link_directories(${HMS_BEAGLE_LIBRARY_DIRS})

set(PHYLO_LIBS
  bpp-core
  bpp-seq
  bpp-phyl
  smctc
  gsl
  gslcblas
  smctc
  jsoncpp
  ${HMS_BEAGLE_LIBRARIES})

add_executable(sts ${CMAKE_SOURCE_DIR}/src/sts.cc)
target_link_libraries(sts ${PHYLO_LIBS})

add_executable(run-tests ${CMAKE_SOURCE_DIR}/src/test.cc)
target_link_libraries(run-tests ${PHYLO_LIBS})
add_test(all_tests run-tests)

# Ensure dependency on run-tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS run-tests)
