#include "gamma_branch_length_proposer.h"

#include <cmath>
#include "gsl/gsl_randist.h"
#include "smctc.hh"

namespace sts
{
namespace moves
{
/// Propose a branch length

/// \returns Branch length
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
