#ifndef STS_MOVES_UNIFORM_BRANCH_LENGTH_PROPOSER_HPP
#define STS_MOVES_UNIFORM_BRANCH_LENGTH_PROPOSER_HPP

#include "base_branch_length_proposer.h"

namespace sts
{
namespace moves
{

/// \class Uniform_branch_length_proposer
/// \brief Propose branch lengths from a uniform distribution.
class Uniform_branch_length_proposer : public Base_branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an uniform distribution on [0, 2 \c mean].
    /// \param mean Mean of uniform distribution
    explicit Uniform_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double);

    /// Mean of uniform distribution
    double mean;
protected:
    double propose_bl(smc::rng *rng);
};

} // namespace moves
} // namespace sts

#endif // STS_MOVES_UNIFORM_BRANCH_LENGTH_PROPOSER_HPP
