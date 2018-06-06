//
//  delta_exchange_mcmc_move.h
//  sts
//
//  Created by Mathieu Fourment on 4/06/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef delta_exchange_mcmc_move_h
#define delta_exchange_mcmc_move_h

#include "online_mcmc_move.h"
#include "tree_particle.h"

#include "smctc.hh"

namespace sts { namespace online {

	// Forwards
	class CompositeTreeLikelihood;
	
	class DeltaExchangeMCMCMove : public OnlineMCMCMove
	{
	public:
		DeltaExchangeMCMCMove(CompositeTreeLikelihood& calculator, const std::vector<std::string>& parameters={},
						   const double lambda=3.0, bool transformed=true);
		~DeltaExchangeMCMCMove();
		int proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
		int proposeMove(TreeParticle& particle, smc::rng* rng);
	private:
		CompositeTreeLikelihood& calculator;
		bool _transformed;
	};
	
}}

#endif /* delta_exchange_mcmc_move_h */
