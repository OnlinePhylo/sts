//
//  scale_mcmc_move.cpp
//  sts
//
//  Created by Mathieu Fourment on 19/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "scale_mcmc_move.h"
#include "online_util.h"

namespace sts { namespace online {
    
    
    
    double ScaleMCMCMove::propose(TreeParticle& particle, const std::vector<bpp::Parameter*>* parameters, smc::rng* rng){
        
        _compositeTreeLilkelihood.initialize(*particle.model, *particle.rateDist, *particle.tree);
        
        double orig_logP = particle.logP;// _compositeTreeLilkelihood();
        
        const double scaler = _delta + (rng->UniformS() * ((1.0 / _delta) - _delta));
        
        if(parameters == nullptr){
            double logRatio = -std::log(scaler);
            std::vector<bpp::Node*> nodes = onlineAvailableEdges(*particle.tree);
            size_t idx = rng->UniformDiscrete(0, nodes.size() - 1);
            bpp::Node* n = nodes[idx];
            const double orig_dist = n->getDistanceToFather();
             n->setDistanceToFather(scaler*orig_dist);
            
            double new_logP = _compositeTreeLilkelihood();
            particle.logP = new_logP;
            
            double mh_ratio = std::exp(new_logP + logRatio - orig_logP);
            
            if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
                return 1;
            } else {
                // Rejected
                n->setDistanceToFather(scaler*orig_dist);
            }
        }
        else{
            double logRatio = -std::log(scaler);
            size_t idx = rng->UniformDiscrete(0, parameters->size() - 1);
            
            const double old_v = parameters->at(idx)->getValue();
            parameters->at(idx)->setValue(old_v*scaler);
            
            double new_logP = _compositeTreeLilkelihood();
            particle.logP = new_logP;
            
            double mh_ratio = std::exp(new_logP + logRatio - orig_logP);
            
            if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
                return 1;
            } else {
                // Rejected
                parameters->at(idx)->setValue(old_v);
            }
        }
        particle.logP = orig_logP;
        return 0;
    }
    
}}
