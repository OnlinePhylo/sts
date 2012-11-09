/// \file smc_init.hpp
/// \brief Smc_init class
/// \author metatangle, inc.

#ifndef STS_MOVES_SMC_INIT_HPP
#define STS_MOVES_SMC_INIT_HPP

#include <memory>
#include "smctc.hh"
#include "sts/likelihood/forest_likelihood.hpp"
#include "sts/particle/state.hpp"

namespace sts
{
namespace moves
{
/// Particle initializer.

/// Uses as log_likelihood the \\perp ll
class Smc_init
{
public:
    explicit Smc_init(sts::likelihood::Forest_likelihood& log_likelihood) : log_likelihood(log_likelihood) {};
    virtual smc::particle<particle::Particle> operator()(smc::rng*);
    virtual ~Smc_init() {};
protected:
    sts::likelihood::Forest_likelihood log_likelihood;
};

///A function to initialise particles

/// \param rng A pointer to the random number generator which is to be used
smc::particle<particle::Particle> Smc_init::operator()(smc::rng* rng)
{
    particle::Particle value;
    // initial particles have all sequences uncoalesced
    value = std::make_shared<particle::State>();
    // Note that the likelihood of the equivalent \perp particles doesn't matter. We set it to zero.
    return smc::particle<particle::Particle>(value, 0.);
}

} // namespace moves
} // namespace sts
#endif // STS_MOVES_SMC_INIT_HPP
