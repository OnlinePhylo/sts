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
    
LocalMCMCMove::LocalMCMCMove(CompositeTreeLikelihood& calculator, const double lambda) : OnlineMCMCMove(lambda), _calculator(calculator){}
    
    LocalMCMCMove::~LocalMCMCMove()
    {
        // Debug bits
        if(n_attempted > 0) {
            std::clog << "Multiplier_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << std::endl;
        }
    }
    
    int LocalMCMCMove::proposeMove(long, smc::particle<TreeParticle>& particle, smc::rng* rng)
    {
        // Choose an edge at random
        TreeParticle* value = particle.GetValuePointer();
        return proposeMove(*value, rng);
    }
    
    int LocalMCMCMove::proposeMove(TreeParticle& particle, smc::rng* rng){
        std::vector<bpp::Node*> nodes = onlineAvailableEdges(*particle.tree);
        size_t idx_central = rng->UniformDiscrete(0, nodes.size() - 1); // index of central branch
        size_t idx_up = rng->UniformDiscrete(0, 1); // choose branch to switch
        size_t idx_down = rng->UniformDiscrete(0, 1); // choose other branch to switch
        bpp::Node* central = nodes[idx_central];
        bpp::Node* down = central->getSon(idx_down);
        bpp::Node* up = (idx_up == 1 ? central->getFather()->getSon(1-central->getFather()->getSonPosition(central)) : central->getFather());
        
        const double orig_central_dist = central->getDistanceToFather();
        const double orig_down_dist = down->getDistanceToFather();
        const double orig_up_dist = up->getDistanceToFather(); // up does not change relative to total length
        
        double orig_total_dist = orig_central_dist + orig_down_dist + orig_up_dist;
        double u1 = rng->UniformS();
        double u2 = rng->UniformS();
        double scaler = exp((u1 - 0.5)*_lambda);
        
        const double new_total_dist = orig_total_dist*scaler;
        const double new_up_dist = orig_up_dist*scaler;
        const double new_up_central_dist = u2*orig_total_dist*scaler;
        
        double new_central_dist;
        double new_down_dist;
        
        // topology changed
        if(new_up_central_dist < new_up_dist){
            new_central_dist = new_up_dist - new_up_central_dist;
            new_down_dist = new_total_dist - new_up_dist;
        }
        else{
            new_central_dist = new_up_central_dist - new_up_dist;
            new_down_dist = new_total_dist - new_up_central_dist;
            central->setDistanceToFather(new_central_dist);
            down->setDistanceToFather(new_down_dist);
            up->setDistanceToFather(new_up_dist);
        }
        
        _calculator.initialize(*particle.model, *particle.rateDist, *particle.tree);
        
        double orig_ll = _calculator();
        
        const Proposal p = positive_real_multiplier(orig_dist, 1e-6, 100.0, _lambda, rng);
        n->setDistanceToFather(p.value);
        double new_ll = _calculator();
        
        particle.logP = new_ll;
        double mh_ratio = std::exp(new_ll + std::log(p.hastingsRatio) - orig_ll);
        if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
            return 1;
        } else {
            // Rejected
            particle.logP = orig_ll;
            n->setDistanceToFather(orig_dist);
            return 0;
        }
    }


}}
