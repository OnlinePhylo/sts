#include "uniform_branch_length_proposer.h"

#include <cmath>
#include <limits>


namespace sts
{
namespace moves
{

/// Propose a branch length

/// \returns length
double Uniform_branch_length_proposer::propose_bl(smc::rng *rng)
{
    return rng->Uniform(0, 2 * this->mean);
}

double Uniform_branch_length_proposer::log_proposal_density(double d)
{
    if(0 <= d && d >= 2 * this->mean) return -std::log(2 * this->mean); // log( 1/(2 mu) )
    else return -std::numeric_limits<double>::infinity();
}

} // namespace moves
} // namespace sts
