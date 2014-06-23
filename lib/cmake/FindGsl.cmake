if (GSL_INCLUDES)
  # Already in cache, be silent
  set (GSL_FIND_QUIETLY TRUE)
endif (GSL_INCLUDES)

find_path(GSL_INCLUDES NAMES "gsl/gsl_rng.h" PATHS ENV CPATH)
find_library(GSL_LIBRARY
  NAMES gsl
  PATHS ENV LIBS ENV LIBRARY_PATH)
find_library(GSL_CBLAS_LIBRARY
  NAMES gslcblas
  PATHS ENV LIBS ENV LIBRARY_PATH)
set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Gsl
                                   DEFAULT_MSG
                                   GSL_LIBRARIES
                                   GSL_INCLUDES)
mark_as_advanced(GSL_LIBRARIES GSL_INCLUDES)
