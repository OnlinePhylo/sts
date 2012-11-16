#ifndef STS_MOVES_BL_PROPOSAL_FN_H
#define STS_MOVES_BL_PROPOSAL_FN_H

#include <functional>
#include "particle.h"
#include "smctc.hh"

namespace sts
{
namespace moves
{
/// Branch length proposal function.
/// Accepts two parameters: a Node with initialized edges and a random source; returns the log-likelihood.
typedef std::function<double(particle::Particle, smc::rng*)> Bl_proposal_fn;
}
}

#endif