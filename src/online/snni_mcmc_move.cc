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
	
	NNIMCMCMove::NNIMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator, const std::vector<std::string>& parameters, const double lambda, bool delayed) : OnlineMCMCMove(calculator, parameters, lambda), _delayed(delayed){}
	
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
			
			if(_delayed){
				// swap b and c
				i->setSon(1-posj, c);
				j->setSon(1-idx_d, b);

				_calculator[index]->initialize(*particle.model, *particle.rateDist, *particle.tree);

				double new_ll2 = _calculator[index]->operator()();
				particle.logP = new_ll2;
				double mh_ratio1 = std::min(1.0, std::exp(new_ll - new_ll2));
				double mh_ratio2 = std::exp(new_ll2 + std::log(1.0 - mh_ratio1) - orig_ll - std::log1p(-exp(new_ll-orig_ll)));
				//std::cout << mh_ratio << " " << mh_ratio2 <<std::endl;
				if(mh_ratio2 >= 1.0 || rng->UniformS() < mh_ratio2) {
					return 1;
				} else {

					// Rejected
					particle.logP = orig_ll;

					i->setSon(1-posj, b);
					j->setSon(1-idx_d, c);
				}
			}
			return 0;
		}
	}
	
	
}}