#ifndef STS_MOVES_GAMMA_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_GAMMA_BRANCH_LENGTH_PROPOSER_HPP

#include <cassert>
#include <cmath>
#include <utility>
#include "gsl/gsl_randist.h"
#include "smctc.hh"

#include "sts/particle/state.hpp"
#include "sts/moves/base_branch_length_proposer.hpp"

namespace sts
{
namespace moves
{

/// \class Gamma_branch_length_proposer
/// \brief Propose branch lengths from an gamma distribution.
class Gamma_branch_length_proposer : public Base_branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an gamma distribution with mean
    /// \c mean.
    /// \param mean Mean of gamma distribution
    explicit Gamma_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double);

    /// Mean of gamma distribution
    double mean;
protected:
    double propose_bl(smc::rng *rng);
};

/// Propose a branch length

/// \returns A pair consisting of: (branch_length, likelihood)
double Gamma_branch_length_proposer::propose_bl(smc::rng *rng)
{
    return rng->Gamma(2.0, this->mean / 2.0);
}

double Gamma_branch_length_proposer::log_proposal_density(double d)
{
    double gamma_p = gsl_ran_gamma_pdf(d, 2.0, this->mean / 2.0);
    return std::log(gamma_p);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_GAMMA_BRANCH_LENGTH_PROPOSER_HPP
