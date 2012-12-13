#include "exponential_branch_length_proposer.h"
#include <cassert>
#include <cmath>
#include <utility>

namespace sts
{
namespace moves
{
/// Propose a branch length

/// \returns A pair consisting of: (branch_length, likelihood)
double Exponential_branch_length_proposer::propose_single_branch_length(smc::rng *rng)
{
    return rng->Exponential(this->mean);
}

double Exponential_branch_length_proposer::log_proposal_density(double d)
{
    return -(std::log(mean) + (d / mean));
}

} // namespace moves
} // namespace sts
