#ifndef STS_ONLINE_MULTIPLIER_MCMC_MOVE_H
#define STS_ONLINE_MULTIPLIER_MCMC_MOVE_H

#include "online_mcmc_move.h"
#include "tree_particle.h"

#include "smctc.hh"

namespace sts { namespace online {

// Forwards
class CompositeTreeLikelihood;

class MultiplierMCMCMove : public OnlineMCMCMove
{
public:
	MultiplierMCMCMove(CompositeTreeLikelihood& calculator, const std::vector<std::string>& parameters={},
                         const double lambda=3.0);
    ~MultiplierMCMCMove();
    int proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
    int proposeMove(TreeParticle& particle, smc::rng* rng);
private:
    CompositeTreeLikelihood& calculator;
};

}}

#endif
