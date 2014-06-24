/// \file bpp_shim.h
/// \author metatangle, inc.
/// \brief Helpers to get things in and out of bpp.
#ifndef STS_LIKELIHOOD_BPP_SHIM_H
#define STS_LIKELIHOOD_BPP_SHIM_H

#include <iostream>
#include <vector>

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Sequence.h>

namespace sts
{
namespace likelihood
{

void blit_vector_to_array(double *, const std::vector<double> &);
void blit_matrix_to_array(double *, const bpp::Matrix<double> &);
void blit_transpose_matrix_to_array(double *, const bpp::Matrix<double> &);
std::vector<double> get_partials(const bpp::Sequence& sequence, const bpp::SubstitutionModel& model, const size_t n_categories=1);

// Debug functions
void print_matrix(const double*, int, int);
void print_matrix(const bpp::Matrix<double> *);

} // namespace likelihood
} // namespace sts

#endif // STS_LIKELIHOOD_BPP_SHIM_H
