/// \file smc_init.hpp
/// \brief smc_init class
/// \author metatangle, inc.

#ifndef STS_MOVES_SMC_INIT_HPP
#define STS_MOVES_SMC_INIT_HPP

#include <memory>
#include "smctc.hh"
#include "sts/likelihood/forest_likelihood.hpp"
#include "sts/particle/phylo_particle.hpp"

namespace sts
{
namespace moves
{
/// Particle initializer.

/// Uses as log_likelihood the \\perp ll
class smc_init
{
public:
    explicit smc_init(sts::likelihood::forest_likelihood& log_likelihood) : log_likelihood(log_likelihood) {};
    virtual smc::particle<particle::particle> operator()(smc::rng*);
protected:
    sts::likelihood::forest_likelihood log_likelihood;
};

///A function to initialise particles

/// \param rng A pointer to the random number generator which is to be used
smc::particle<particle::particle> smc_init::operator()(smc::rng* rng)
{
    particle::particle value;
    // initial particles have all sequences uncoalesced
    value = std::make_shared<particle::phylo_particle>();
    // Note that the likelihood of the equivalent \perp particles doesn't matter. We set it to zero.
    return smc::particle<particle::particle>(value, 0.);
}

} // namespace moves
} // namespace sts
#endif // STS_MOVES_SMC_INIT_HPP
