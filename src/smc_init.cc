/// \file smc_init.hpp
/// \brief Smc_init class
/// \author metatangle, inc.

#include "smc_init.h"

#include "state.h"

#include <memory>

namespace sts
{
namespace moves
{
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
