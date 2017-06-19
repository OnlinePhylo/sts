//
//  Prior.hpp
//  sts
//
//  Created by Mathieu Fourment on 14/06/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef prior_hpp
#define prior_hpp

#include <stdio.h>

#include "likelihood.h"

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/Function/Functions.h>


namespace sts {
    namespace online {
        
        class Prior : public Likelihood{
            
        public:
            Prior(std::vector<std::string> paramNames, std::function<double(double)> logDensity):
                _logDensity(logDensity), _paramNames(paramNames), _parameters(){}
            
            virtual ~Prior(){}
            
            const std::vector<std::string>& getParameterNames() const;
            
            void setParameters(const std::vector<const bpp::Parameter*> parameters);
            
            virtual double calculateLogLikelihood();
            
        private:
            std::function<double(double)> _logDensity;
            
        protected:
            std::vector<std::string> _paramNames;
            std::vector<const bpp::Parameter*> _parameters;
        };
    }
}

#endif /* prior_h */
