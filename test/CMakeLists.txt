# TESTS
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
include_directories(${source_dir}/include)
ExternalProject_Get_Property(googletest binary_dir)
set(GTEST_LIBRARIES ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a
                    ${binary_dir}/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a
                    ${CMAKE_THREAD_LIBS_INIT})

set(STS_TEST_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/test_sts_beagle_tree_likelihood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_sts_gsl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_sts_likelihood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_sts_parsing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_sts_log_tricks.cpp
  )
add_executable(run-tests EXCLUDE_FROM_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/test.cc ${STS_TEST_FILES})
add_dependencies(run-tests googletest)
target_link_libraries(run-tests sts-static ${STS_PHYLO_LIBS} ${GTEST_LIBRARIES})
add_test(all-tests ${CMAKE_CURRENT_BINARY_DIR}/run-tests)
