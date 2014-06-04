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

# Subdirectories
# boost
include_directories(SYSTEM lib/boost_1_54_0)

# SMCTC
add_definitions(-DSMCTC_HAVE_BGL)
include_directories(lib/smctc/include)
add_subdirectory(lib/smctc)


# JSONCPP
include_directories(lib/jsoncpp/include)
add_subdirectory(lib/jsoncpp)

include_directories(lib/include)

# BEAGLE
find_package(PkgConfig REQUIRED)
pkg_check_modules(HMS_BEAGLE hmsbeagle-1 REQUIRED)
include_directories(${HMS_BEAGLE_INCLUDE_DIRS})
link_directories(${HMS_BEAGLE_LIBRARY_DIRS})

# Bio++
include(FindBpp)
find_package(Bpp REQUIRED)
include_directories(SYSTEM ${BPP_INCLUDES})

add_subdirectory(src)

# Tests
enable_testing()
add_subdirectory(test)
