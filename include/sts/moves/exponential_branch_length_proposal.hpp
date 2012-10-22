#ifndef STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSAL_HPP
#define STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSAL_HPP

#include <cassert>
#include <cmath>
#include <utility>
#include "gsl/gsl_randist.h"
#include "smctc.hh"

#include "sts/particle/phylo_node.hpp"

namespace sts
{
namespace moves
{

/// \class exponential_branch_length_proposal
/// \brief Propose branch lengths from an exponential distribution.
class exponential_branch_length_proposal
{
public:
    /// Instantiate a new BL proposer where branch lengths are drawn from an exponential distribution with mean
    /// \c mean.
    /// \param mean Mean of exponential distribution
    explicit exponential_branch_length_proposal(double mean) : mean(mean) {};

    double operator()(std::shared_ptr<particle::phylo_node>, smc::rng*);
protected:
    typedef std::pair<double, double> doubles;
    /// Mean of exponential distribution
    double mean;
    doubles propose_bl(smc::rng *rng);
};

/// Propose branch lengths on \c node.

/// \param node Phylo node to operate on. Child edges of \c node must be initialized. <b>This function changes child edge branch
///  lengths.</b>
/// \param rng Random number generator
/// \returns The log-likelihood of the proposal
double exponential_branch_length_proposal::operator()(std::shared_ptr<particle::phylo_node> node, smc::rng *rng)
{
    // TODO: different BLs for child1 and child2
    doubles d1 = propose_bl(rng); //, d2 = propose_bl(rng);

    // Children should be initialized
    assert(node->child1);
    assert(node->child2);
    node->child1->length = d1.first;
    node->child2->length = d1.first;
    //node->child2->length = d2.first;
    //return d1.second + d2.second;
    return d1.second;
}

/// Propose a branch length

/// \returns A pair consisting of: (branch_length, likelihood)
exponential_branch_length_proposal::doubles exponential_branch_length_proposal::propose_bl(smc::rng *rng)
{
    double length = rng->Exponential(this->mean);
    double log_like = std::log(gsl_ran_exponential_pdf(length, this->mean));
    return doubles(length, log_like);
}

} // namespace moves
} // namespace sts

#endif // STS_MOVES_EXPONENTIAL_BRANCH_LENGTH_PROPOSAL_HPP
