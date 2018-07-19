//
//  snni_mcmc_move.cc
//  sts
//
//  Created by Mathieu Fourment on 6/07/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "snni_mcmc_move.h"

#include "online_util.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace sts { namespace online {
	
	NNIMCMCMove::NNIMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator, const std::vector<std::string>& parameters, const double lambda) : OnlineMCMCMove(calculator, parameters, lambda){}
	
	NNIMCMCMove::~NNIMCMCMove()
	{
		// Debug bits
		if(n_attempted > 0) {
			std::clog << "snni_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << std::endl;
		}
	}
	
	int NNIMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
	{
		// Choose an edge at random
		TreeParticle* value = particle.GetValuePointer();
		return proposeMove(*value, rng);
	}
	
	int NNIMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng){
		
		size_t index = 0;
#if defined(_OPENMP)
		index = omp_get_thread_num();
#endif
		double orig_ll = particle.logP;
		
		std::vector<bpp::Node*> nodes = onlineAvailableInternalEdges(*particle.tree);
		size_t idx_j = rng->UniformDiscrete(0, nodes.size() - 1); // index of branch ij
		size_t idx_d = rng->UniformDiscrete(0, 1); // choose node d to move
		bpp::Node* j = nodes[idx_j];
		// we don't want a child of the root
		while (!j->getFather()->hasFather()) {
			idx_j = rng->UniformDiscrete(0, nodes.size() - 1);
			j = nodes[idx_j];
		}
		bpp::Node* i = j->getFather();
		bpp::Node* d = j->getSon(idx_d);
		bpp::Node* c = j->getSon(1-idx_d);
		size_t posj = i->getSonPosition(j);
		bpp::Node* b = i->getSon(1-posj);
		
		// swap b and d
		i->setSon(1-posj, d);
		j->setSon(idx_d, b);
		
		_calculator[index]->initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double new_ll = _calculator[index]->operator()();
		particle.logP = new_ll;
		double mh_ratio = std::exp(new_ll - orig_ll);
		if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
			return 1;
		} else {
			// Rejected
			particle.logP = orig_ll;
			
			i->setSon(1-posj, b);
			j->setSon(idx_d, d);
			
			return 0;
		}
	}
	
	
}}