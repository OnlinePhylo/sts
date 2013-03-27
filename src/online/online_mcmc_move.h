#ifndef STS_ONLINE_ONLINE_MCMC_MOVE_H
#define STS_ONLINE_ONLINE_MCMC_MOVE_H

#include <smctc.hh>

namespace sts { namespace online {

// Forwards
class Tree_particle;

class Online_mcmc_move
{
public:
    Online_mcmc_move();
    virtual ~Online_mcmc_move() {};

    double acceptance_probability() const;

    int operator()(long, smc::particle<Tree_particle>&, smc::rng*);
protected:
    virtual int propose_move(long time, smc::particle<Tree_particle>& particle, smc::rng* rng) = 0;

    /// Number of times the move was attempted
    unsigned int n_attempted;
    /// Number of times the move was accepted
    unsigned int n_accepted;
};

}}

#endif
