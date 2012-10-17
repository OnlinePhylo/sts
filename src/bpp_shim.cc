#include "bpp_shim.hh"
#include <iostream>


void
blit_vector_to_array(double *arr, const std::vector<double> &vec)
{
    for(std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it)
        *arr++ = *it;
}

void
blit_matrix_to_array(double *arr, const bpp::Matrix<double> &matrix)
{
    int cols = matrix.getNumberOfColumns(), rows = matrix.getNumberOfRows();
    for(int i = 0; i < rows; ++i) {
        blit_vector_to_array(arr, matrix.row(i));
        arr += cols;
    }
}

void blit_transpose_matrix_to_array(double *arr, const bpp::Matrix<double> &matrix)
{
    int cols = matrix.getNumberOfColumns(), rows = matrix.getNumberOfRows();
    for(int j = 0; j < rows; ++j) {
        for(int i = 0; i < rows; ++i) {
            *arr++ = matrix(i, j);
        }
    }
}

/// Get a vector of partial states from a sequence, substitution model, and alphabet.
/// This does not have site compression, which would require buying into Bio++'s sitepatterns object.
std::vector<double> get_partials(const bpp::Sequence& sequence, const bpp::SubstitutionModel &model, const bpp::Alphabet *alphabet)
{
    unsigned int n_states = model.getNumberOfStates(), n_sites = sequence.size();

    std::vector<double> partials(n_sites * n_states);

    for(unsigned int site = 0; site < n_sites; site++) {
        for(unsigned int i = 0; i < n_states; i++) {
            partials[n_states * site + i] = model.getInitValue(i, sequence.getValue(site));
        }
    }

    return partials;
}

/// Some debugging functions
void print_matrix(const bpp::Matrix<double> *m)
{
    for(int i = 0; i < m->getNumberOfRows(); i++) {
        for(int j = 0; j < m->getNumberOfColumns(); j++) {
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

/*

A hypothesized replacement for AD's get_partials but for arbitrary alphabets and models, something like
DRASDRTreeLikelihoodData::initLikelihoods

Note that they do the nice thing of only allocating memory for the extant patterns.

// Turn a sequence into partials.
std::vector<double> get_partials(const std::string& sequence, const SubstitutionModel & model) {

    std::vector< double > partials(n * N_STATES);

    printf("%g\n",model.getInitValue(2,3));

    // Init leaves likelihoods:
    const Sequence* seq;
    try
    {
      seq = &sites.getSequence(node->getName());
    }
    leavesLikelihoods_leaf->resize(nbDistinctSites_);
    for (unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      Vdouble* leavesLikelihoods_leaf_i = &(*leavesLikelihoods_leaf)[i];
      leavesLikelihoods_leaf_i->resize(nbStates_);
      int state = seq->getValue(i);
      double test = 0.;
      for (unsigned int s = 0; s < nbStates_; s++)
      {
        //Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
        //otherwise value set to 0:
        ( * leavesLikelihoods_leaf_i)[s] = model.getInitValue(s, state);
        test += ( * leavesLikelihoods_leaf_i)[s];
      }
      if (test < 0.000001) std::cerr << "WARNING!!! Likelihood will be 0 for this site." << std::endl;
    }
  }

}

*/
