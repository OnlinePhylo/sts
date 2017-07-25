//
//  exchange_mcmc_move.cpp
//  sts
//
//  Created by Mathieu Fourment on 19/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "exchange_mcmc_move.h"

namespace sts { namespace online {
    


double ExchangeMCMCMove::propose(TreeParticle& particle, const std::vector<bpp::Parameter*>* parameters, smc::rng* rng){
    
    _compositeTreeLilkelihood.initialize(*particle.model, *particle.rateDist, *particle.tree);
    double orig_logP = _compositeTreeLilkelihood();
    
    double logRatio = 0;
    size_t idx1 = rng->UniformDiscrete(0, parameters->size() - 1);
    size_t idx2 = idx1;
    while(idx1 == idx2){
        idx2 = rng->UniformDiscrete(0, parameters->size() - 1);
    }
    
    const double d = rng->UniformS() * _delta;
    const double old_v1 = parameters->at(idx1)->getValue();
    const double old_v2 = parameters->at(idx2)->getValue();
    parameters->at(idx1)->setValue(old_v1-d);
    parameters->at(idx2)->setValue(old_v2+d);
    
    double new_logP = _compositeTreeLilkelihood();
    
    double mh_ratio = std::exp(new_logP + logRatio - orig_logP);
    particle.logP = new_logP;
    
    if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
        return 1;
    } else {
        // Rejected
        particle.logP = orig_logP;
        parameters->at(idx1)->setValue(old_v1);
        parameters->at(idx2)->setValue(old_v2);
        return 0;
    }
}

}}
