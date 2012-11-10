cmake_minimum_required (VERSION 2.8.1)
cmake_policy(SET CMP0015 NEW)
enable_testing()

project(STS CXX)

set(CMAKE_CXX_FLAGS "-g -std=c++0x -Wall -pedantic")

# Subdirectories
add_subdirectory(lib/smctc)
add_subdirectory(lib/jsoncpp)
include_directories(lib/include)
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

set(CPP_FILES
  src/base_branch_length_proposer.cc
  src/bpp_shim.cc
  src/delta_branch_length_proposer.cc
  src/edge.cc
  src/exponential_branch_length_proposer.cc
  src/forest_likelihood.cc
  src/gamma_branch_length_proposer.cc
  src/gsl.cc
  src/log/json_logger.cc
  src/log/sampler.cc
  src/mcmc_move.cc
  src/node.cc
  src/online_calculator.cc
  src/rooted_merge.cc
  src/smc_init.cc
  src/smc_move.cc
  src/state.cc
  src/uniform_bl_mcmc_move.cc
  src/uniform_branch_length_proposer.cc
  src/util.cc)
add_library(sts-static STATIC ${CPP_FILES})

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
target_link_libraries(sts sts-static ${PHYLO_LIBS})
install(TARGETS sts
        RUNTIME DESTINATION bin)

add_executable(run-tests EXCLUDE_FROM_ALL ${CMAKE_SOURCE_DIR}/src/test.cc)
target_link_libraries(run-tests sts-static ${PHYLO_LIBS})
add_test(all_tests run-tests)

# Ensure dependency on run-tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS run-tests)
