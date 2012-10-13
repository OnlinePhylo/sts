#include "bpp_shim.hh"

static void
blit_vector_to_array(double *arr, const std::vector<double> &vec)
{
    for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it)
        *arr++ = *it;
}

static void
blit_matrix_to_array(double *arr, const bpp::Matrix<double> &matrix)
{
    int cols = matrix.getNumberOfColumns(), rows = matrix.getNumberOfRows();
    for (int i = 0; i < rows; ++i) {
        blit_vector_to_array(arr, matrix.row(i));
        arr += cols;
    }
}
