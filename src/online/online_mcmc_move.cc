#include "online_mcmc_move.h"
#include "tree_particle.h"

namespace sts { namespace online {


OnlineMCMCMove::OnlineMCMCMove(double lambda) :
    n_attempted(0),
    n_accepted(0),
    _lambda(lambda),
    _target(0.234),
    _min(0.001),
    _max(100)
{}

int OnlineMCMCMove::operator()(long time, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    ++n_attempted;
    const int result = proposeMove(time, particle, rng);
    if(result){
        ++n_accepted;
    }
    _lambda = tune();
    TreeParticle* value = particle.GetValuePointer();
    return result;
}

double OnlineMCMCMove::acceptanceProbability() const
{
    if(!n_attempted)
        return 0.0;

    return static_cast<double>(n_accepted) / n_attempted;
}
    
double OnlineMCMCMove::tune()
{
    double delta = 1./sqrt(n_attempted);
    delta = fmin(0.01, delta);
    double logTuner = log(_lambda);
    double accRatio = acceptanceProbability();
    if(accRatio > _target){
        logTuner += delta;
    }
    else {
        logTuner -= delta;
    }
    double t = exp(logTuner);
    
    if (t <= _min || t >= _max)
        t = _lambda;
    return t;
}

}} // Namespaces
