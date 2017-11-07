#ifndef STS_ONLINE_NODE_SLIDER_MCMC_MOVE_H
#define STS_ONLINE_NODE_SLIDER_MCMC_MOVE_H

#include "online_mcmc_move.h"
#include "tree_particle.h"
#include "smctc.hh"

namespace sts { namespace online {

// Forwards
class CompositeTreeLikelihood;

class NodeSliderMCMCMove : public OnlineMCMCMove
{
public:
    NodeSliderMCMCMove(CompositeTreeLikelihood& calculator,
                          const double lambda=3.0);
    ~NodeSliderMCMCMove();
    std::pair<int, double> proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
private:
    CompositeTreeLikelihood& calculator;
};

}}

#endif
