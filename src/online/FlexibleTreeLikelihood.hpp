//
//  FlexibleTreeLikelihood.hpp
//  sts
//
//  Created by Mathieu Fourment on 22/02/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#ifndef FlexibleTreeLikelihood_hpp
#define FlexibleTreeLikelihood_hpp

#include <stdio.h>
#include <vector>

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace sts {
    namespace online {
        
        class FlexibleTreeLikelihood{
            
        public:

            virtual ~FlexibleTreeLikelihood(){}
            
            virtual void initialize(bpp::TreeTemplate<bpp::Node>& tree, bpp::SubstitutionModel &model, bpp::DiscreteDistribution& rateDist) = 0;
            
            virtual double calculateLogLikelihood() = 0;
            
            virtual double calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength) = 0;
            
            // Compute derivatives at node node with current tree
//            virtual void calculateDerivatives(const bpp::Node& node, double* d1, double* d2) = 0;
            
            // Compute derivatives of pendant branch with taxon taxonName
            virtual void calculateDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2) = 0;
            
            virtual void calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2) = 0;
        };
    }
}
    
#endif /* FlexibleTreeLikelihood_hpp */
