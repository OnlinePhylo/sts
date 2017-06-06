configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/online/sts_config.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/online/sts_config.h"
  )

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/online")
file(GLOB STS_CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/online/*.cc)

IF(NO_BEAGLE)
	list(REMOVE_ITEM STS_CPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/online/beagle_flexible_tree_likelihood.cc")
	list(REMOVE_ITEM STS_CPP_FILES "${CMAKE_CURRENT_SOURCE_DIR}/online/beagle_tree_likelihood.cc")
ENDIF(NO_BEAGLE)

add_library(sts-static STATIC ${STS_CPP_FILES})

set(STS_PHYLO_LIBS
  ${BPP_LIBRARIES}
  smctc
  ${GSL_LIBRARIES}
  smctc
  jsoncpp
  lcfit_cpp-static
  ${HMS_BEAGLE_LIBRARIES})

# Export
set(STS_PHYLO_LIBS ${STS_PHYLO_LIBS} PARENT_SCOPE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_executable(sts-online ${CMAKE_CURRENT_SOURCE_DIR}/online/sts_online.cc)
target_link_libraries(sts-online sts-static ${STS_PHYLO_LIBS})
install(TARGETS sts-online
        RUNTIME DESTINATION bin)
