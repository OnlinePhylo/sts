#ifndef STS_ONLINE_NODE_SLIDER_SMC_MOVE_H
#define STS_ONLINE_NODE_SLIDER_SMC_MOVE_H

#include "smctc.hh"

namespace sts { namespace online {

// Forwards
class CompositeTreeLikelihood;
class TreeParticle;

class NodeSliderSMCMove
{
public:
    NodeSliderSMCMove(CompositeTreeLikelihood& calculator,
                      const double lambda=0.191);
    void operator()(long, smc::particle<TreeParticle>&, smc::rng*);
private:
    CompositeTreeLikelihood& calculator;
    double lambda;
};

}}

#endif
