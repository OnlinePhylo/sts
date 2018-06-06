//
//  delta_exchange_mcmc_move.cc
//  sts
//
//  Created by Mathieu Fourment on 4/06/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "delta_exchange_mcmc_move.h"
#include "composite_tree_likelihood.h"
#include "online_util.h"
#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>

using namespace bpp;

namespace sts { namespace online {
	
	DeltaExchangeMCMCMove::DeltaExchangeMCMCMove(CompositeTreeLikelihood& calculator, const std::vector<std::string>& parameters,
										   const double lambda, bool transformed) :
	OnlineMCMCMove(parameters, lambda),
	calculator(calculator),
	_transformed(transformed)
	{}
	
	DeltaExchangeMCMCMove::~DeltaExchangeMCMCMove()
	{
		// Debug bits
		if(n_attempted > 0) {
			std::clog << "Delta_exchange_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << " " <<_lambda << std::endl;
		}
	}
	
	int DeltaExchangeMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
	{
		// Choose an edge at random
		TreeParticle* value = particle.GetValuePointer();
		return proposeMove(*value, rng);
	}
	
	int DeltaExchangeMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng){
		size_t parameterCount = _parameters.size();
		if(_transformed) parameterCount++;
		const long idx1 = rng->UniformDiscrete(0, parameterCount-1);
		long idx2 = idx1;
		while(idx1 == idx2){
			idx2 = rng->UniformDiscrete(0, parameterCount-1);
		}
		
		std::vector<double> values = particle.model->getFrequencies();
		std::vector<double> backup = values;

		double w = rng->UniformS()*_lambda;
		// watch for infinite loop
		while(w >= values[idx1] || w >= values[idx2]){
			w = rng->UniformS()*_lambda;
		}
		double orig1 = values[idx1];
		double orig2 = values[idx2];
		values[idx1] += w;
		values[idx2] -= w;
		
		calculator.initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double orig_ll = calculator();
		//assert(orig_ll==particle.logP);
		
		std::map<int, double> a;
		for(int i = 0; i < 4; i++){
			a[i] = values[i];
		}
		particle.model->setFreq(a);
		
		// need to reinitialize
		calculator.initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double new_ll = calculator();
		
		particle.logP = new_ll;
		double mh_ratio = std::exp(new_ll - orig_ll);
		if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
			return 1;
		} else {
			// Rejected
			particle.logP = orig_ll;
			if(_transformed){
				a.clear();
				for(int i = 0; i < 4; i++){
					a[i] = backup[i];
				}
				particle.model->setFreq(a);
			}
			return 0;
		}
	}
	
}}
