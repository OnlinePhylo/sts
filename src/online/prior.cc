//
//  Prior.cpp
//  sts
//
//  Created by Mathieu Fourment on 14/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "prior.h"

namespace sts { namespace online {

    double Prior::calculateLogLikelihood(){
        double logP = 0;
        for(const bpp::Parameter* p : _parameters){
            logP += _logDensity(p->getValue());
        }
        return logP;
    }

    const std::vector<std::string>& Prior::getParameterNames() const{
        return _paramNames;
    }
    
    void Prior::setParameters(const std::vector<const bpp::Parameter*> parameters){
        _parameters = parameters;
    }
}}
