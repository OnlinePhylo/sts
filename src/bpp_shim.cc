/// \file bpp_shim.cc
/// \author metatangle, inc.
/// \brief Helpers to get things in and out of bpp.

#include "bpp_shim.h"
#include <algorithm>

namespace sts
{
namespace likelihood
{
/// Copy the contents of vec into arr
/// \param arr destination array, with length at least vec.size()
/// \param vec Vector to copy from
void blit_vector_to_array(double *arr, const std::vector<double> &vec)
{
    for(std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it)
        *arr++ = *it;
}

/// Copy the contents of matrix into arr, in row-major order
/// \param arr destination array, with length at least nrows x ncols in length
/// \param vec Vector to copy from
void blit_matrix_to_array(double *arr, const bpp::Matrix<double> &matrix)
{
    int cols = matrix.getNumberOfColumns(), rows = matrix.getNumberOfRows();
    for(int i = 0; i < rows; ++i) {
        blit_vector_to_array(arr, matrix.row(i));
        arr += cols;
    }
}

/// Copy and transpose matrix into arr
/// \param arr destination array, with length at least nrows x ncols in length
/// \param vec Vector to copy from
void blit_transpose_matrix_to_array(double *arr, const bpp::Matrix<double> &matrix)
{
    int cols = matrix.getNumberOfColumns(), rows = matrix.getNumberOfRows();
    for(int j = 0; j < cols; ++j) {
        for(int i = 0; i < rows; ++i) {
            *arr++ = matrix(i, j);
        }
    }
}

/// Get a vector of partial states from a sequence, substitution model, and alphabet.
/// \param sequence input sequence
/// \param model Substitution model
/// \param alphabet Sequence alphabet
/// \param n_categories Number of rate categories
/// \returns vector with length \c{model.getNumberOfStates()*sequence.size()}
std::vector<double> get_partials(const bpp::Sequence& sequence, const bpp::SubstitutionModel &model, const size_t n_categories)
{
    unsigned int n_states = model.getNumberOfStates(), n_sites = sequence.size();

    std::vector<double> partials(n_sites * n_states * n_categories);
                //size_t idx = n_states * n_sites * cat + n_states * site + i;

    for(unsigned int site = 0; site < n_sites; site++) {
        for(unsigned int i = 0; i < n_states; i++) {
            size_t idx = n_states * site + i;
            partials[idx] = model.getInitValue(i, sequence.getValue(site));
        }
    }

    // Make `n_categories` copies of the partials
    size_t per_category = n_sites * n_states;
    for(size_t c = 1; c < n_categories; c++) {
        std::copy(partials.begin(),
                  partials.begin() + per_category,
                  partials.begin() + (per_category * c));
    }

    return partials;
}

/// Some debugging functions
void print_matrix(const bpp::Matrix<double> *m)
{
    for(unsigned int i = 0; i < m->getNumberOfRows(); i++) {
        for(unsigned int j = 0; j < m->getNumberOfColumns(); j++) {
            std::cerr << (*m)(i, j) << '\t';

        }
        std::cerr << std::endl;
    }
}

void print_matrix(const double* d, int rows, int cols)
{
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            std::cerr << d[i * cols + j] << '\t';
        }
        std::cerr << std::endl;
    }
}

} // namespace likelihood
} // namespace sts
