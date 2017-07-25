//
//  mcmc_move.hpp
//  sts
//
//  Created by Mathieu Fourment on 19/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef mcmc_move_hpp
#define mcmc_move_hpp

#include <stdio.h>

#include "tree_particle.h"

#include "smctc.hh"

#include <Bpp/Numeric/Parameter.h>

namespace sts { namespace online {
    
    class MCMCMove{
        
    public:
        
        virtual ~MCMCMove(){}
        
        virtual double propose(TreeParticle& particle, const std::vector<bpp::Parameter*>* parameters, smc::rng* rng) = 0;
        
        //virtual double tune();
    };
}}

#endif /* mcmc_move_hpp */
