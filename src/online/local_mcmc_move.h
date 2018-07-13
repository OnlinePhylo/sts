//
//  local_mcmc_move.hpp
//  sts
//
//  Created by Mathieu Fourment on 19/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef local_mcmc_move_hpp
#define local_mcmc_move_hpp

#include <stdio.h>

#include "online_mcmc_move.h"
#include "tree_particle.h"
#include "composite_tree_likelihood.h"

namespace sts { namespace online {
  
    class LocalMCMCMove : public OnlineMCMCMove{
    public:
        LocalMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator, const std::vector<std::string>& parameters={}, const double lambda=3.0);
        
        virtual ~LocalMCMCMove();
        
        int proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
        int proposeMove(TreeParticle&, smc::rng*);
    };
    
}}
#endif /* local_mcmc_move_hpp */
