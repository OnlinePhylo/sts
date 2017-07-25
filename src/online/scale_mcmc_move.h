//
//  scale_mcmc_move.hpp
//  sts
//
//  Created by Mathieu Fourment on 19/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef scale_mcmc_move_hpp
#define scale_mcmc_move_hpp

#include <stdio.h>

#include "mcmc_move.h"
#include "composite_tree_likelihood.h"

namespace sts { namespace online {
    
    class ScaleMCMCMove : public MCMCMove{
        
    public:
        ScaleMCMCMove(CompositeTreeLikelihood& compositeTreeLilkelihood, double delta=0.5) : _compositeTreeLilkelihood(compositeTreeLilkelihood), _delta(delta){}
        
        virtual ~ScaleMCMCMove(){}
        
        virtual double propose(TreeParticle& particle, const std::vector<bpp::Parameter*>* parameters, smc::rng* rng);
        
    protected:
        CompositeTreeLikelihood& _compositeTreeLilkelihood;
        double _delta;
    };
    
}}

#endif /* scale_mcmc_move_hpp */
