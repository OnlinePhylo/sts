#ifndef STS_ONLINE_NODE_SLIDER_MCMC_MOVE_H
#define STS_ONLINE_NODE_SLIDER_MCMC_MOVE_H

#include "tree_particle.h"
#include "smctc.hh"

namespace sts { namespace online {

// Forwards
class Beagle_tree_likelihood;

class Node_slider_mcmc_move
{
public:
    Node_slider_mcmc_move(Beagle_tree_likelihood& calculator,
                         const double lambda=3.0);
    // debug - log
    ~Node_slider_mcmc_move();
    int operator()(long, smc::particle<Tree_particle>&, smc::rng*);
private:
    Beagle_tree_likelihood& calculator;
    double lambda;
    size_t proposed;
    size_t accepted;
};

}}

#endif
