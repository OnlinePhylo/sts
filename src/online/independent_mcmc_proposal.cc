//
//  independent_mcmc_proposal.cc
//  sts
//
//  Created by Mathieu Fourment on 20/08/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "independent_mcmc_proposal.h"
#include "composite_tree_likelihood.h"
#include "online_util.h"
#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace bpp;

namespace sts { namespace online {
	
	IndependentMCMCMove::IndependentMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator,
									   std::vector<double> params,
									   std::vector<std::unique_ptr<Transform>> transforms,
									   std::function<double(smc::rng*, std::vector<double> parameters)> random,
									   std::function<double(double x, std::vector<double> parameters)> logP,
									   const double lambda) :
	OnlineMCMCMove(calculator, {}, lambda),
	_params(params),
	_transforms(std::move(transforms)),
	_random(random),
	_logP(logP){
	}
	
	IndependentMCMCMove::IndependentMCMCMove(const IndependentMCMCMove& proposal):
	OnlineMCMCMove(proposal._calculator, proposal._parameters, proposal._lambda),
	_params(proposal._params), _random(proposal._random), _logP(proposal._logP){
		for(int i = 0; i < proposal._transforms.size(); i++){
			_transforms.emplace_back(proposal._transforms[i]->clone());
		}
	}
	
	IndependentMCMCMove::~IndependentMCMCMove()
	{
		// Debug bits
		if(n_attempted > 0) {
			std::clog << "indepedent_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << " " <<_lambda << std::endl;
		}
	}
	
	int IndependentMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
	{
		// Choose an edge at random
		TreeParticle* value = particle.GetValuePointer();
		return proposeMove(*value, rng);
	}
	
	int IndependentMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng){
		
		size_t indexCalculator = 0;
#if defined(_OPENMP)
		indexCalculator = omp_get_thread_num();
#endif
		double orig_ll = particle.logP;

		double new_transformed_value = _random(rng, _params);
		double new_logQ = _logP(new_transformed_value, _params);
		
		std::unique_ptr<Transform>& transform = _transforms[0];
		
		std::string name = transform->getNames()[0];
		double new_value = transform->inverse_transform(new_transformed_value);
		
		double orig_value;
		if (particle.model->hasParameter(name)) {
			orig_value = particle.model->getParameterValue(name);
			particle.model->setParameterValue(name, new_value);
		}
		else if (particle.rateDist->hasParameter(name)) {
			orig_value = particle.rateDist->getParameterValue(name);
			particle.rateDist->setParameterValue(name, new_value);
		}
		const double orig_transformed_value = transform->transform(orig_value);
		
		double logJacobian = transform->logJacobian(orig_value) - transform->logJacobian(new_value);
		double logQ = _logP(orig_transformed_value, _params);
		
		_calculator[indexCalculator]->initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double new_ll = _calculator[indexCalculator]->operator()();
		
		particle.logP = new_ll;
		
		double mh_ratio = std::exp(new_ll + logQ - orig_ll - new_logQ + logJacobian);
		if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
			return 1;
		} else {
			// Rejected
			particle.logP = orig_ll;
			if (particle.model->hasParameter(name)) {
				particle.model->setParameterValue(name, orig_value);
			}
			else if (particle.rateDist->hasParameter(name)) {
				particle.rateDist->setParameterValue(name, orig_value);
			}

			return 0;
		}
		
	}
	
}}
