cmake_minimum_required(VERSION 2.8.1)
cmake_policy(SET CMP0015 NEW)
project(sts CXX)

OPTION(USE_OPENMP "Use OpenMP if available." OFF)

#
# cmake modules
#

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake/cxx11/Modules)

#
# Versioning
#

include(GetGitRevisionDescription)
git_describe(GIT_DESCRIPTION --always --long)

#
# Compiler settings
#

include(FindCXXFeatures)

if(NOT _HAS_CXX11_FLAG)
   message( FATAL_ERROR "-std=c++11 flag not supported." )
endif()

set(CMAKE_CXX_FLAGS "-g -Wall -Wextra -pedantic ${CXX11_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")

IF(USE_OPENMP)
	find_package(OpenMP)
	IF (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	ENDIF(OPENMP_FOUND)
ENDIF(USE_OPENMP)

#
# In-tree dependencies
#

OPTION(USE_BEAGLE "Use BEAGLE library." OFF)

# boostable
include_directories(SYSTEM lib/boostable)

# lcfit
include_directories(lib/lcfit/lcfit_src)
include_directories(lib/lcfit/lcfit_cpp_src)
add_subdirectory(lib/lcfit)

# SMCTC
add_definitions(-DSMCTC_HAVE_BGL)
include_directories(SYSTEM lib/smctc/include)
add_subdirectory(lib/smctc)

# JsonCpp
include_directories(lib/jsoncpp/include)
add_subdirectory(lib/jsoncpp)

# tclap
include_directories(lib/include)


#
# External dependencies
#

# pkg-config
find_package(PkgConfig REQUIRED)

# BEAGLE
IF(USE_BEAGLE)
	pkg_check_modules(HMS_BEAGLE hmsbeagle-1 REQUIRED)
	include_directories(${HMS_BEAGLE_INCLUDE_DIRS})
	link_directories(${HMS_BEAGLE_LIBRARY_DIRS})
ELSE(USE_BEAGLE)
	SET(NO_BEAGLE "1")
	add_definitions("-DNO_BEAGLE")
ENDIF(USE_BEAGLE)

# Bio++
find_package(Bpp REQUIRED)
include_directories(SYSTEM ${BPP_INCLUDES})

# GSL
include(FindGsl)
find_package(Gsl REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

#
# sts proper
#

add_subdirectory(src)


#
# Tests
#

IF (BUILD_TESTING)
	enable_testing()
	add_subdirectory(test)
ENDIF (BUILD_TESTING)
