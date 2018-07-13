//
//  snni_mcmc_move.h
//  sts
//
//  Created by Mathieu Fourment on 6/07/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef snni_mcmc_move_h
#define snni_mcmc_move_h

#include <stdio.h>

#include "online_mcmc_move.h"
#include "tree_particle.h"
#include "composite_tree_likelihood.h"

namespace sts { namespace online {
	
	class NNIMCMCMove : public OnlineMCMCMove{
	public:
		NNIMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator, const std::vector<std::string>& parameters={}, const double lambda=3.0);
		
		virtual ~NNIMCMCMove();
		
		int proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
		int proposeMove(TreeParticle&, smc::rng*);
		
		virtual double tune(){return _lambda;}
	};
	
}}

#endif /* snni_mcmc_move_h */
