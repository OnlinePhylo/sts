#include "delta_branch_length_proposer.h"

#include <limits>

namespace sts
{
namespace moves
{
/// Propose a branch length

/// \returns branch length
double Delta_branch_length_proposer::propose_single_branch_length(smc::rng *rng)
{
    return this->mean;
}

double Delta_branch_length_proposer::log_proposal_density(double d)
{
    if(d == this->mean) return 0;
    else return -std::numeric_limits<double>::infinity();
}

} // namespace moves
} // namespace sts
