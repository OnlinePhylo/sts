/// \file bpp_shim.hh
/// \author metatangle, inc.
/// \brief Helpers to get things in and out of bpp.

#ifndef __bpp_shim__
#define __bpp_shim__

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Sequence.h>

void blit_vector_to_array(double *, const std::vector<double> &);
void blit_matrix_to_array(double *, const bpp::Matrix<double> &);
std::vector<double> get_partials(const std::string&, const bpp::SubstitutionModel &, const bpp::Alphabet *);

#endif //  __bpp_shim__
