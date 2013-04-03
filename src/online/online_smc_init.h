/// \file smc_init.h
/// \brief Smc_init class
/// \author metatangle, inc.

#ifndef STS_MOVES_ONLINE_SMC_INIT_H
#define STS_MOVES_ONLINE_SMC_INIT_H

#include "smctc.hh"
#include <vector>

namespace sts { namespace online {

// Forwards
class TreeParticle;

/// Particle initializer.

/// Uses as log_likelihood the \\perp ll.
class OnlineSMCInit
{
public:
    OnlineSMCInit(const std::vector<TreeParticle>& p) :
        particles(p),
        i(0) { };

    smc::particle<TreeParticle> operator()(smc::rng*);
private:
    const std::vector<TreeParticle> particles;
    size_t i;
};

}} // namespaces
#endif // STS_MOVES_ONLINE_SMC_INIT
