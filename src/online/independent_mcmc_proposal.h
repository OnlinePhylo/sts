//
//  independent_mcmc_proposal.h
//  sts
//
//  Created by Mathieu Fourment on 20/08/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef independent_mcmc_proposal_h
#define independent_mcmc_proposal_h

#include "online_mcmc_move.h"
#include "tree_particle.h"
#include "transform.h"
#include "smctc.hh"

namespace sts { namespace online {
	
	// Forwards
	class CompositeTreeLikelihood;
	
	class IndependentMCMCMove : public OnlineMCMCMove
	{
	public:
		IndependentMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
						 std::vector<double> params,
						 std::vector<std::unique_ptr<Transform>> transforms,
						 std::function<double(smc::rng*, std::vector<double> parameters)> random,
						 std::function<double(double x, std::vector<double> parameters)> logP,
						 const double lambda=3.0);
		
		IndependentMCMCMove(const IndependentMCMCMove& proposal);
		
		~IndependentMCMCMove();
		int proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
		int proposeMove(TreeParticle& particle, smc::rng* rng);
	private:
		std::vector<double> _params;
		std::vector<std::unique_ptr<Transform>> _transforms;
		std::function<double(smc::rng*, std::vector<double> parameters)> _random;
		std::function<double(double x, std::vector<double> parameters)> _logP;
	};
	
}}


#endif /* independent_mcmc_proposal_h */
