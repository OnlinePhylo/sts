//
//  dirichletPrior.cpp
//  sts
//
//  Created by Mathieu Fourment on 14/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "dirichlet_prior.h"

#include <gsl/gsl_sf.h>

namespace sts { namespace online {
        
    double DirichletPrior::calculateLogLikelihood(){
        return gsl_sf_lngamma(4);
    }

}}
