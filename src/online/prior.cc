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
        if (_paramNames.size() == 0) {
            return _logDensity(0);
        }
        for( std::string paramString: _parameters.getParameterNames()){
            const bpp::Parameter& p = _parameters.getParameter(paramString);
            logP += _logDensity(p.getValue());
        }
        return logP;
    }

    const std::vector<std::string>& Prior::getParameterNames() const{
        return _paramNames;
    }
    
    void Prior::setParameters(bpp::ParameterList parameters){
        _parameters = parameters;
    }
}}
