/// \file online_smc_init.cc

#include "online_smc_init.h"
#include "tree_particle.h"

using namespace std;

namespace sts { namespace moves {

///A function to initialise particles
smc::particle<sts::Tree_particle> Online_smc_init::operator()(smc::rng*)
{

    sts::Tree_particle value = particles[i++ % particles.size()];
    return smc::particle<sts::Tree_particle>(value, 0.);
}

}} // namespaces
