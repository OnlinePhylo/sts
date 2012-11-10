#ifndef STS_MOVES_DELTA_BRANCH_LENGTH_PROPOSER_H
#define STS_MOVES_DELTA_BRANCH_LENGTH_PROPOSER_H

#include "smctc.hh"

#include "state.h"
#include "base_branch_length_proposer.h"

namespace sts
{
namespace moves
{

/// \class Delta_branch_length_proposer
/// \brief "Propose" branch lengths from a delta distribution.
class Delta_branch_length_proposer : public Base_branch_length_proposer
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an delta distribution with mean
    /// \c mean.
    /// \param mean Mean of delta distribution
    explicit Delta_branch_length_proposer(double mean) : mean(mean) {};
    double log_proposal_density(double);

    /// Mean of delta distribution
    double mean;
protected:
    double propose_bl(smc::rng *);
};

} // namespace moves
} // namespace sts

#endif // STS_MOVES_DELTA_BRANCH_LENGTH_PROPOSER_H
