//
//  likelihood.hpp
//  sts
//
//  Created by Mathieu Fourment on 14/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef likelihood_h
#define likelihood_h

#include <stdio.h>

namespace sts {
    namespace online {
        
        class Likelihood {
            
        public:
            virtual double calculateLogLikelihood() = 0;
            
            virtual ~Likelihood(){}
        };
    }
}
#endif /* likelihood_h */
