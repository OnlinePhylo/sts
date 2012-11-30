/// \file smc_init.hpp
/// \brief Smc_init class
/// \author metatangle, inc.

#ifndef STS_MOVES_SMC_INIT_H
#define STS_MOVES_SMC_INIT_H

#include "forest_likelihood.h"
#include "particle.h"

#include "smctc.hh"

namespace sts
{
namespace moves
{
/// Particle initializer.

/// Uses as log_likelihood the \\perp ll.
class Smc_init
{
public:
    virtual smc::particle<particle::Particle> operator()(smc::rng*);
    virtual ~Smc_init() {};
protected:
};

} // namespace moves
} // namespace sts
#endif // STS_MOVES_SMC_INIT_H
