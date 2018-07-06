//
//  snni_mcmc_move.cc
//  sts
//
//  Created by Mathieu Fourment on 6/07/2018.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "snni_mcmc_move.h"

#include "online_util.h"

namespace sts { namespace online {
	
	NNIMCMCMove::NNIMCMCMove(CompositeTreeLikelihood& calculator, const std::vector<std::string>& parameters, const double lambda) : OnlineMCMCMove(parameters, lambda), _calculator(calculator){}
	
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
		
		
		_calculator.initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double orig_ll = _calculator();
		assert(orig_ll==particle.logP);
		
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
		bpp::Node* c = j->getSon(1-idx_d);
		bpp::Node* a = i->getFather();
		const double wai = i->getDistanceToFather();
		const double wij = j->getDistanceToFather();
		
		// detach node j
		size_t posj = i->getSonPosition(j);
		i->setSon(posj, c);
		
		// insert node j between a and i
		size_t posi = a->getSonPosition(i);
		a->setSon(posi, j);
		j->setSon(j->getSonPosition(c), i);
		
		j->setDistanceToFather(wai);
		i->setDistanceToFather(wij);
		
		_calculator.initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double new_ll = _calculator();
		particle.logP = new_ll;
		double mh_ratio = std::exp(new_ll - orig_ll);
		if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
			return 1;
		} else {
			// Rejected
			particle.logP = orig_ll;
			
			// detach node i
			size_t posi = j->getSonPosition(i);
			j->setSon(posi, c);
			
			// insert node i between a and j
			size_t posj = a->getSonPosition(j);
			a->setSon(posj, i);
			i->setSon(i->getSonPosition(c), j);

			j->setDistanceToFather(wij);
			i->setDistanceToFather(wai);
			
			return 0;
		}
	}
	
	
}}