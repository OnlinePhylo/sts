//
//  adaptive_mcmc_move.h
//  sts
//
//  Created by Mathieu Fourment on 5/06/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef adaptive_mcmc_move_h
#define adaptive_mcmc_move_h

#include "online_mcmc_move.h"
#include "tree_particle.h"
#include "Transform.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>

#include "smctc.hh"

namespace sts { namespace online {
	
	// Forwards
	class CompositeTreeLikelihood;
	
	class AdaptiveMCMCMove : public OnlineMCMCMove
	{
	public:
		AdaptiveMCMCMove(CompositeTreeLikelihood& calculator,
						 gsl_matrix& L,
						 gsl_vector& mu,
						 std::vector<Transform*> transforms,
						 const double lambda=3.0);
		
		AdaptiveMCMCMove(const AdaptiveMCMCMove& adapt);
		
		~AdaptiveMCMCMove();
		int proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
		int proposeMove(TreeParticle& particle, smc::rng* rng);
	private:
		CompositeTreeLikelihood& calculator;
		gsl_matrix& _L;
		gsl_vector& _mu;
		gsl_vector* _result;
		gsl_vector* _work;
		std::vector<Transform*> _transforms;
	};
	
}}

#endif /* adatptive_mcmc_mode_h */