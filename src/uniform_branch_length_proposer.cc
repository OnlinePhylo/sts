#include "uniform_branch_length_proposer.h"

#include <cmath>
#include <limits>


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
/// Propose a branch length

/// \returns length
double Uniform_branch_length_proposer::propose_bl(smc::rng *rng)
{
    return rng->Uniform(0, 2 * this->mean);
}

double Uniform_branch_length_proposer::log_proposal_density(double d)
{
    if(0 <= d && d >= 2 * this->mean) return -std::log(2*this->mean); // log( 1/(2 mu) )
    else return -std::numeric_limits<double>::infinity();
}

} // namespace moves
} // namespace sts
