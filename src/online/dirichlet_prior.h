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
			DirichletPrior(std::vector<std::string> paramNames, const std::vector<double>& alphas=std::vector<double>()):Prior(paramNames, 0), _alphas(alphas){}
            
            virtual ~DirichletPrior(){}
            
            virtual double calculateLogLikelihood();
		private:
			std::vector<double> _alphas;
        };
    }
}
#endif /* dirichletPrior_hpp */
