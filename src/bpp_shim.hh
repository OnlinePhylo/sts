/// \file bpp_shim.hh
/// \author metatangle, inc.
/// \brief Helpers to get things in and out of bpp.

#include <Bpp/Numeric/Matrix/Matrix.h>

static void blit_vector_to_array(double *, const std::vector<double> &);
static void blit_matrix_to_array(double *, const bpp::Matrix<double> &);
