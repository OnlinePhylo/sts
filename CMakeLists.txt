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
include_directories(lib/smctc/include)

# JSONCPP
include_directories(lib/jsoncpp/include)

# BEAGLE
find_package(PkgConfig REQUIRED)
pkg_check_modules(HMS_BEAGLE hmsbeagle-1 REQUIRED)
include_directories(${HMS_BEAGLE_INCLUDE_DIRS})
link_directories(${HMS_BEAGLE_LIBRARY_DIRS})

set(STS_CPP_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/base_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/bpp_shim.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/child_swap_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/delta_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/edge.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/exponential_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/forest_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gamma_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gsl.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/log/json_logger.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/log/sampler.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/metropolis_hastings_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/node.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/node_deleter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online_calculator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rooted_merge.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/smc_init.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/smc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/state.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/uniform_bl_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/uniform_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/util.cc)
add_library(sts-static STATIC ${STS_CPP_FILES})

set(STS_PHYLO_LIBS
  bpp-core
  bpp-seq
  bpp-phyl
  smctc
  gsl
  gslcblas
  smctc
  jsoncpp
  ${HMS_BEAGLE_LIBRARIES})
set(STS_PHYLO_LIBS ${STS_PHYLO_LIBS} PARENT_SCOPE)

add_executable(sts ${CMAKE_CURRENT_SOURCE_DIR}/src/sts.cc)
target_link_libraries(sts sts-static ${STS_PHYLO_LIBS})
install(TARGETS sts
        RUNTIME DESTINATION bin)

add_executable(run-tests EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/src/test.cc)
target_link_libraries(run-tests sts-static ${STS_PHYLO_LIBS})
add_test(all_tests run-tests)

# Ensure dependency on run-tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS run-tests)
