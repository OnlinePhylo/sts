//
//  local_mcmc_move.cpp
//  sts
//
//  Created by Mathieu Fourment on 19/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "local_mcmc_move.h"

#include "online_util.h"

namespace sts { namespace online {
    
LocalMCMCMove::LocalMCMCMove(std::vector<std::unique_ptr<CompositeTreeLikelihood>>& calculator, const std::vector<std::string>& parameters, const double lambda) : OnlineMCMCMove(calculator, parameters, lambda){}
    
    LocalMCMCMove::~LocalMCMCMove()
    {
        // Debug bits
        if(n_attempted > 0) {
            std::clog << "Local_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << " "<<  _lambda << std::endl;
        }
    }
    
    int LocalMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
    {
        // Choose an edge at random
        TreeParticle* value = particle.GetValuePointer();
        return proposeMove(*value, rng);
    }
    
    int LocalMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng){
		
		size_t index = 0;
#if defined(_OPENMP)
		index = omp_get_thread_num();
#endif
		_calculator[index]->initialize(*particle.model, *particle.rateDist, *particle.tree);
		
		double orig_ll = _calculator[index]->operator()();
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
		bpp::Node* d = j->getSon(idx_d);
        //bpp::Node* b = i->getSon(1-i->getSonPosition(j));
		bpp::Node* a = i->getFather();
		
        const double wai = i->getDistanceToFather();
        const double wij = j->getDistanceToFather();
        const double wjc = c->getDistanceToFather();
        
		double wac = wjc + wij + wai;
        double u1 = rng->UniformS();
        double u2 = rng->UniformS();
        double scaler = exp((u1 - 0.5)*_lambda);
        
        const double wpac = wac*scaler;
        const double wpai = wai*scaler;
        const double wpaj = u2*wac*scaler;
        
        // topology changed
        if(wpaj < wpai){
            const double wpij = wpai - wpaj;
            const double wpic = wpac - wpai;
			// detach node j
			size_t posj = i->getSonPosition(j);
			i->setSon(posj, c);
			
			// insert node j between a and i
			size_t posi = a->getSonPosition(i);
			a->setSon(posi, j);
			j->setSon(j->getSonPosition(c), i);
			
			j->setDistanceToFather(wpaj);
			i->setDistanceToFather(wpij);
			c->setDistanceToFather(wpic);
        }
        else{
            const double wpij = wpaj - wpai;
            const double wpjc = wpac - wpaj;
            j->setDistanceToFather(wpij);
            c->setDistanceToFather(wpjc);
            i->setDistanceToFather(wpai);
        }
        
        _calculator[index]->initialize(*particle.model, *particle.rateDist, *particle.tree);
		
        double new_ll = _calculator[index]->operator()();
        particle.logP = new_ll;
        double mh_ratio = std::exp(new_ll + 3.0*std::log(scaler) - orig_ll);
		if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
            return 1;
		} else {
            // Rejected
            particle.logP = orig_ll;
			
			// topology changed
			if(wpaj < wpai){
				// detach node i
				size_t posi = j->getSonPosition(i);
				j->setSon(posi, c);
				
				// insert node i between a and j
				size_t posj = a->getSonPosition(j);
				a->setSon(posj, i);
				i->setSon(i->getSonPosition(c), j);
			}
			i->setDistanceToFather(wai);
			j->setDistanceToFather(wij);
			c->setDistanceToFather(wjc);
			
            return 0;
        }
    }


}}
