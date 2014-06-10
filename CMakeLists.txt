cmake_minimum_required (VERSION 2.8.1)
cmake_policy(SET CMP0015 NEW)
project(sts)


# Versioning
set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake/cxx11/Modules
    ${CMAKE_MODULE_PATH})

include(GetGitRevisionDescription)
git_describe(GIT_DESCRIPTION --always --long)

# Check for C++11
include(FindCXXFeatures)

if(NOT _HAS_CXX11_FLAG)
   message( FATAL_ERROR "-std=c++11 flag not supported." )
endif()

set(CMAKE_CXX_FLAGS "-g -Wall -Wextra -pedantic ${CXX11_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")

configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

# Subdirectories
# boost
include_directories(SYSTEM lib/boost_1_55_0)
include_directories(SYSTEM lib/boostable)

add_definitions(-DSMCTC_HAVE_BGL)
add_subdirectory(lib/smctc)
add_subdirectory(lib/jsoncpp)
add_subdirectory(lib/lcfit)
include_directories(lib/include)
include_directories(src)
include_directories(src/online)

# SMCTC
include_directories(lib/smctc/include)

# JSONCPP
include_directories(lib/jsoncpp/include)

# lcfit
include_directories(lib/lcfit/lcfit_src)
include_directories(lib/lcfit/lcfit_cpp_src)

# BEAGLE
find_package(PkgConfig REQUIRED)
pkg_check_modules(HMS_BEAGLE hmsbeagle-1 REQUIRED)
include_directories(${HMS_BEAGLE_INCLUDE_DIRS})
link_directories(${HMS_BEAGLE_LIBRARY_DIRS})

# Bio++
include(FindBpp)
find_package(Bpp REQUIRED)
include_directories(${BPP_INCLUDES})

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
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/attachment_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/beagle_tree_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/branch_length_prior.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/composite_tree_likelihood.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/guided_online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/lcfit_online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/lcfit_rejection_sampler.cc
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
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/tripod_optimizer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/online/uniform_online_add_sequence_move.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lib/lcfit/lcfit_cpp_src/lcfit_cpp.cc
  )
add_library(sts-static STATIC ${STS_CPP_FILES})

set(STS_PHYLO_LIBS
  ${BPP_LIBRARIES}
  smctc
  gsl
  gslcblas
  smctc
  jsoncpp
  lcfit-static
  ${HMS_BEAGLE_LIBRARIES})

add_executable(sts ${CMAKE_CURRENT_SOURCE_DIR}/src/sts.cc)
target_link_libraries(sts sts-static ${STS_PHYLO_LIBS})
install(TARGETS sts
        RUNTIME DESTINATION bin)

add_executable(sts-online ${CMAKE_CURRENT_SOURCE_DIR}/src/online/sts_online.cc)
target_link_libraries(sts-online sts-static ${STS_PHYLO_LIBS})

# Tests
enable_testing()
add_subdirectory(test)
