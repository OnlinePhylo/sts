/// \file online_smc_init.cc

#include "online_smc_init.h"
#include "tree_particle.h"

using namespace std;

namespace sts { namespace online {

///A function to initialise particles
smc::particle<Tree_particle> Online_smc_init::operator()(smc::rng*)
{

    Tree_particle value = particles[i++ % particles.size()];
    return smc::particle<Tree_particle>(value, 0.);
}

}} // namespaces
