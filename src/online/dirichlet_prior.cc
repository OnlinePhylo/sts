//
//  dirichletPrior.cpp
//  sts
//
//  Created by Mathieu Fourment on 14/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "dirichlet_prior.h"

#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>

namespace sts { namespace online {
        
    double DirichletPrior::calculateLogLikelihood(){
        size_t dimAlpha = _alphas.size();
        if(dimAlpha == 0){
            return gsl_sf_lngamma(_paramNames.size());
        }
        else{
            std::vector<double> thetas(dimAlpha, 0);
            std::transform(_parameters.begin(), _parameters.end(), thetas.begin(),
                           [](const bpp::Parameter* p) -> double { return p->getValue(); });
            return gsl_ran_dirichlet_lnpdf(dimAlpha, _alphas.data(), thetas.data());
        }
    }

}}
