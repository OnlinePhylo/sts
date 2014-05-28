cmake_minimum_required (VERSION 2.8.1)
cmake_policy(SET CMP0015 NEW)
project(STS CXX)

set(CMAKE_CXX_FLAGS "-g -std=c++0x -Wall -pedantic")

# Versioning
set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake
    ${CMAKE_MODULE_PATH})
include(GetGitRevisionDescription)
git_describe(GIT_DESCRIPTION --always --long)
configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

# Subdirectories
# boost
include_directories(SYSTEM lib/boost_1_54_0)
add_definitions(-DSMCTC_HAVE_BGL)
add_subdirectory(lib/smctc)
add_subdirectory(lib/jsoncpp)
include_directories(lib/include)
include_directories(src)
include_directories(src/online)

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
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online_calculator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rooted_merge.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/smc_init.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/smc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/state.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/uniform_bl_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/uniform_branch_length_proposer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/util.cc

  # ONLINE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/beagle_tree_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/branch_length_prior.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/composite_tree_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/guided_online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/likelihood_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/multiplier_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/multiplier_proposal.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/multiplier_smc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/node_slider_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/node_slider_smc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/online_mcmc_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/online_smc_init.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/online_util.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/tree_particle.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/uniform_online_add_sequence_move.cc
  )
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

add_executable(sts ${CMAKE_CURRENT_SOURCE_DIR}/src/sts.cc)
target_link_libraries(sts sts-static ${STS_PHYLO_LIBS})
install(TARGETS sts
        RUNTIME DESTINATION bin)

add_executable(sts-online ${CMAKE_CURRENT_SOURCE_DIR}/src/online/sts_online.cc)
target_link_libraries(sts-online sts-static ${STS_PHYLO_LIBS})


# TESTS
enable_testing()
find_package(Threads)
include(ExternalProject)
ExternalProject_Add(
    googletest
    URL https://googletest.googlecode.com/files/gtest-1.7.0.zip
    URL_MD5 2d6ec8ccdf5c46b05ba54a9fd1d130d7
    TIMEOUT 10
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON)
ExternalProject_Get_Property(googletest source_dir)
include_directories(${source_dir}/fused-src)
ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIBRARIES ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a
                    ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a
                    ${CMAKE_THREAD_LIBS_INIT})

set(STS_TEST_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/test_sts_beagle_tree_likelihood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/test_sts_gsl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/test_sts_likelihood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/test_sts_parsing.cpp
  )
add_executable(run-tests EXCLUDE_FROM_ALL ${CMAKE_CURRENT_SOURCE_DIR}/src/test.cc ${STS_TEST_FILES})
add_dependencies(run-tests googletest)
target_link_libraries(run-tests sts-static ${STS_PHYLO_LIBS} ${GTEST_LIBRARIES})
add_test(all_tests run-tests)

# Ensure dependency on run-tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS run-tests)
