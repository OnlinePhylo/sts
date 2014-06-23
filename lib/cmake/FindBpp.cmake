if (BPP_INCLUDES)
  # Already in cache, be silent
  set (BPP_FIND_QUIETLY TRUE)
endif (BPP_INCLUDES)

find_path(BPP_CORE_INCLUDES "Bpp/Numeric.all"
  PATHS ENV CPATH)
find_library(BPP_CORE_LIBRARIES bpp-core
  PATHS ENV LIBS ENV LIBRARY_PATH)
find_path(BPP_SEQ_INCLUDES "Bpp/Seq.all"
  PATHS ENV CPATH)
find_library(BPP_SEQ_LIBRARIES bpp-seq
  PATHS ENV LIBS ENV LIBRARY_PATH)
find_path(BPP_PHYL_INCLUDES "Bpp/Phyl.all"
  PATHS ENV CPATH)
find_library(BPP_PHYL_LIBRARIES bpp-phyl
  PATHS ENV LIBS ENV LIBRARY_PATH)
set(BPP_LIBRARIES
  ${BPP_CORE_LIBRARIES}
  ${BPP_SEQ_LIBRARIES}
  ${BPP_PHYL_LIBRARIES})
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
