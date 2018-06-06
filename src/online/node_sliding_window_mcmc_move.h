#ifndef node_sliding_window_mcmc_move_hpp
#define node_sliding_window_mcmc_move_hpp

#include "online_mcmc_move.h"
#include "tree_particle.h"

#include "smctc.hh"

namespace sts { namespace online {
    
    // Forwards
    class CompositeTreeLikelihood;
    
    class SlidingWindowMCMCMove : public OnlineMCMCMove
    {
    public:
		SlidingWindowMCMCMove(CompositeTreeLikelihood& calculator, const std::vector<std::string>& parameters={},
                           const double lambda=3.0);
        ~SlidingWindowMCMCMove();
        int proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
        int proposeMove(TreeParticle&, smc::rng*);
        
    private:
        CompositeTreeLikelihood& calculator;
    };
    
}}

#endif /* node_sliding_window_mcmc_move_hpp */
