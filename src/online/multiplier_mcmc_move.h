#ifndef STS_ONLINE_MULTIPLIER_MCMC_MOVE_H
#define STS_ONLINE_MULTIPLIER_MCMC_MOVE_H

#include "tree_particle.h"
#include "smctc.hh"


namespace sts { namespace online {

// Forwards
class Beagle_tree_likelihood;

class Multiplier_mcmc_move
{
public:
    Multiplier_mcmc_move(Beagle_tree_likelihood& calculator,
                         const double lambda=3.0);
    // debug - log
    ~Multiplier_mcmc_move();
    int operator()(long, smc::particle<Tree_particle>&, smc::rng*);
    double acceptance_rate() const;
private:
    Beagle_tree_likelihood& calculator;
    double lambda;
    size_t proposed;
    size_t accepted;
};

}}

#endif
