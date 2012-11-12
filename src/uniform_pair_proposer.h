#ifndef STS_UNIFORM_PAIR_PROPOSER_H
#define STS_UNIFORM_PAIR_PROPOSER_H

#include <string>
#include <vector>
#include "smctc.hh"
#include <forest_likelihood.h>
#include "node.h"
#include "state.h"
#include "edge.h"

namespace sts
{
namespace moves
{

/// \class Uniform_pair_proposer
/// \brief "Propose" pairs of nodes to join uniformly at random from the forest roots.
class Uniform_pair_proposer
{
public:
    /// Instantiate 
    explicit Uniform_pair_proposer(sts::likelihood::Forest_likelihood& log_likelihood) : log_likelihood(log_likelihood) {};
    void operator()(particle::Particle, smc::rng*, particle::Node_ptr& a, particle::Node_ptr& b, double& fwd_density, double& back_density);    
private:    
    sts::likelihood::Forest_likelihood& log_likelihood;
};


} // namespace moves
} // namespace sts

#endif // STS_UNIFORM_PAIR_PROPOSER_H
