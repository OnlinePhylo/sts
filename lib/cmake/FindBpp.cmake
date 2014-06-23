if (BPP_INCLUDES)
  # Already in cache, be silent
  set (BPP_FIND_QUIETLY TRUE)
endif (BPP_INCLUDES)

# Includes
find_path(BPP_CORE_INCLUDES "Bpp/Numeric.all"
  PATHS ENV CPATH)
find_path(BPP_SEQ_INCLUDES "Bpp/Seq.all"
  PATHS ENV CPATH)
find_path(BPP_PHYL_INCLUDES "Bpp/Phyl.all"
  PATHS ENV CPATH)

# Handle static libaries
if(STATIC_BPP)
  find_library(BPP_CORE_LIBRARIES
    NAMES libbpp-core.a
    PATHS ENV LIBS ENV LIBRARY_PATH)
  find_library(BPP_SEQ_LIBRARIES
    NAMES libbpp-seq.a
    PATHS ENV LIBS ENV LIBRARY_PATH)
  find_library(BPP_PHYL_LIBRARIES
    NAMES libbpp-phyl.a
    PATHS ENV LIBS ENV LIBRARY_PATH)
else()
  find_library(BPP_CORE_LIBRARIES bpp-core
    PATHS ENV LIBS ENV LIBRARY_PATH)
  find_library(BPP_SEQ_LIBRARIES bpp-seq
    PATHS ENV LIBS ENV LIBRARY_PATH)
  find_library(BPP_PHYL_LIBRARIES bpp-phyl
    PATHS ENV LIBS ENV LIBRARY_PATH)
endif()

set(BPP_LIBRARIES
  ${BPP_PHYL_LIBRARIES}
  ${BPP_SEQ_LIBRARIES}
  ${BPP_CORE_LIBRARIES})
set(BPP_INCLUDES
  ${BPP_CORE_INCLUDES}
  ${BPP_SEQ_INCLUDES}
  ${BPP_PHYL_INCLUDES})

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Bpp
                                   DEFAULT_MSG
                                   BPP_LIBRARIES
                                   BPP_INCLUDES)
mark_as_advanced(BPP_LIBRARIES BPP_INCLUDES)
