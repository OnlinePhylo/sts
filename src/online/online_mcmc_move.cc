#include "online_mcmc_move.h"
#include "tree_particle.h"

namespace sts { namespace online {


Online_mcmc_move::Online_mcmc_move() :
    n_attempted(0),
    n_accepted(0)
{}

int Online_mcmc_move::operator()(long time, smc::particle<Tree_particle>& particle, smc::rng* rng)
{
    ++n_attempted;
    const int result = propose_move(time, particle, rng);
    if(result)
        ++n_accepted;
    return result;
}

double Online_mcmc_move::acceptance_probability() const
{
    if(!n_attempted)
        return 0.0;

    return static_cast<double>(n_accepted) / n_attempted;
}

}} // Namespaces
