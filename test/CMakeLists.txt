# TESTS

enable_testing()
find_package(GTest REQUIRED)

find_package(Threads)
include(ExternalProject)

IF(GTEST_FOUND)
	include_directories(${GTEST_INCLUDE_DIRS})
ELSEIF()
	ExternalProject_Add(
		googletest
		#GIT_REPOSITORY https://github.com/google/googletest.git
		URL https://github.com/google/googletest/archive/release-1.7.0.zip
		URL_HASH SHA256=b58cb7547a28b2c718d1e38aee18a3659c9e3ff52440297e965f5edffe34b6d0
		TIMEOUT 10
		CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		# Disable install step
		INSTALL_COMMAND ""
		# Wrap download, configure and build steps in a script to log output
		LOG_DOWNLOAD ON
		LOG_CONFIGURE ON
		LOG_BUILD ON)

	ExternalProject_Get_Property(googletest source_dir)
	include_directories(${source_dir}/googletest/include)
	ExternalProject_Get_Property(googletest binary_dir)

	set(GTEST gtest)
	set(GTEST_MAIN gtest_main)

	set(GTEST_BOTH_LIBRARIES ${binary_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}${GTEST}${CMAKE_STATIC_LIBRARY_SUFFIX}
						${binary_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}${GTEST_MAIN}${CMAKE_STATIC_LIBRARY_SUFFIX})            
ENDIF()


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src/online)

set(STS_TEST_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/test_sts_log_tricks.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_parsimony.cpp
  )
  
IF(USE_BEAGLE)
  set(STS_TEST_FILES
    ${STS_TEST_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/test_sts_beagle_tree_likelihood.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_sts_flexible_tree_likelihood.cpp
  )
ENDIF(USE_BEAGLE)

add_executable(run-tests EXCLUDE_FROM_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/test.cc ${STS_TEST_FILES})

#add_dependencies(run-tests googletest)

target_link_libraries(run-tests sts-static ${STS_PHYLO_LIBS} ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_test(all-tests ${CMAKE_CURRENT_BINARY_DIR}/run-tests)
