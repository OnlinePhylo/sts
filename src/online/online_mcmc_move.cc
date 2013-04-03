#include "online_mcmc_move.h"
#include "tree_particle.h"

namespace sts { namespace online {


OnlineMCMCMove::OnlineMCMCMove() :
    n_attempted(0),
    n_accepted(0)
{}

int OnlineMCMCMove::operator()(long time, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    ++n_attempted;
    const int result = propose_move(time, particle, rng);
    if(result)
        ++n_accepted;
    return result;
}

double OnlineMCMCMove::acceptance_probability() const
{
    if(!n_attempted)
        return 0.0;

    return static_cast<double>(n_accepted) / n_attempted;
}

}} // Namespaces
