/// \file online_smc_init.cc

#include "online_smc_init.h"
#include "tree_particle.h"

using namespace std;

namespace sts { namespace online {

///A function to initialise particles
smc::particle<TreeParticle> OnlineSMCInit::operator()(smc::rng*)
{

    TreeParticle value = particles[i++ % particles.size()];
    return smc::particle<TreeParticle>(value, 0.);
}

}} // namespaces
