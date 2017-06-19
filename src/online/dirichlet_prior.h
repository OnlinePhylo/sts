//
//  dirichletPrior.hpp
//  sts
//
//  Created by Mathieu Fourment on 14/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef dirichletPrior_hpp
#define dirichletPrior_hpp

#include <stdio.h>

#include "prior.h"

namespace sts {
    namespace online {
        
        class DirichletPrior : public Prior{
            
        public:
            DirichletPrior(std::vector<std::string> paramNames):Prior(paramNames, 0){}
            
            virtual ~DirichletPrior(){}
            
            virtual double calculateLogLikelihood();
        };
    }
}
#endif /* dirichletPrior_hpp */
