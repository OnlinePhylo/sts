/// \file smc_init.h
/// \brief Smc_init class
/// \author metatangle, inc.

#ifndef STS_MOVES_ONLINE_SMC_INIT_H
#define STS_MOVES_ONLINE_SMC_INIT_H

#include "smctc.hh"
#include <vector>

namespace sts { namespace online {

// Forwards
class Tree_particle;

/// Particle initializer.

/// Uses as log_likelihood the \\perp ll.
class Online_smc_init
{
public:
    Online_smc_init(const std::vector<Tree_particle>& p) :
        particles(p),
        i(0) { };

    smc::particle<Tree_particle> operator()(smc::rng*);
private:
    const std::vector<Tree_particle> particles;
    size_t i;
};

}} // namespaces
#endif // STS_MOVES_ONLINE_SMC_INIT
