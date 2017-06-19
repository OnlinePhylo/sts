//
//  scale_mcmc_move.cpp
//  sts
//
//  Created by Mathieu Fourment on 19/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "scale_mcmc_move.h"

namespace sts { namespace online {
    
    
    
    double ScaleMCMCMove::propose(TreeParticle& particle, const std::vector<bpp::Parameter*>& parameters, smc::rng* rng){
        
        _compositeTreeLilkelihood.initialize(*particle.model, *particle.rateDist, *particle.tree);
        
        double orig_logP = particle.logP;// _compositeTreeLilkelihood();
        
        const double scaler = (_delta + (rng->UniformS() * ((1.0 / _delta) - _delta)));
        double logRatio = -std::log(scaler);
        size_t idx = rng->UniformDiscrete(0, parameters.size() - 1);
        
        const double old_v = parameters[idx]->getValue();
        parameters[idx]->setValue(old_v*scaler);
        
        double new_logP = _compositeTreeLilkelihood();
        particle.logP = new_logP;
        
        double mh_ratio = std::exp(new_logP + logRatio - orig_logP);
        
        if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
            return 1;
        } else {
            // Rejected
            particle.logP = orig_logP;
            parameters[idx]->setValue(old_v);
            return 0;
        }
    }
    
}}
