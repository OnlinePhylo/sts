#ifndef STS_ONLINE_ONLINE_MCMC_MOVE_H
#define STS_ONLINE_ONLINE_MCMC_MOVE_H

#include <smctc.hh>

#include <Bpp/Numeric/Parameter.h>

namespace sts { namespace online {

// Forwards
class TreeParticle;

class OnlineMCMCMove
{
public:
    OnlineMCMCMove(double lambda=3);
    virtual ~OnlineMCMCMove() {};

    double acceptanceProbability() const;

    int operator()(long, smc::particle<TreeParticle>&, smc::rng*);
    
    int operator()(TreeParticle&, smc::rng*);
    
protected:
    virtual int proposeMove(long time, smc::particle<TreeParticle>& particle, smc::rng* rng) = 0;
    
    virtual int proposeMove(TreeParticle& particle, smc::rng* rng) = 0;

    virtual double tune();
    
    /// Number of times the move was attempted
    unsigned int n_attempted;
    /// Number of times the move was accepted
    unsigned int n_accepted;
    
    double _lambda;
    
    double _target;
    double _min;
    double _max;
};

}}

#endif
