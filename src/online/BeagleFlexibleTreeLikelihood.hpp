//
//  BeagleFlexibleTreeLikleihood.hpp
//  sts
//
//  Created by Mathieu Fourment on 26/06/2016.
//  Copyright © 2016 Mathieu Fourment. All rights reserved.
//

#ifndef BeagleFlexibleTreeLikleihood_hpp
#define BeagleFlexibleTreeLikleihood_hpp

#include <stdio.h>
#include <vector>


#include "AbstractFlexibleTreeLikelihood.hpp"

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "libhmsbeagle/beagle.h"

namespace sts {
    namespace online {
        
        class BeagleFlexibleTreeLikelihood : public AbstractFlexibleTreeLikelihood {
            
        public:
            BeagleFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, bpp::SubstitutionModel &model, bpp::DiscreteDistribution& rateDist);
            
            virtual ~BeagleFlexibleTreeLikelihood();
            
            virtual void initialize(bpp::TreeTemplate<bpp::Node>& tree, bpp::SubstitutionModel &model, bpp::DiscreteDistribution& rateDist);
            
            virtual double calculateLogLikelihood();
            
            virtual double calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength);
            
//            virtual void calculateDerivatives(const bpp::Node& node, double* d1, double* d2);
            
            virtual void calculateDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);
            
        protected:
            
            void updateSubstitutionModel();
            
            void updateSiteModel();
            
            bool traverse(const bpp::Node* node);
            
            void traverseUpper(const bpp::Node* node, int& index);
            
        private:
            double _logLnl;
            
            void registerLeaves(const bpp::SiteContainer& sites);
            
            int _beagleInstance;
            BeagleInstanceDetails _instanceDetails;
            
            //size_t _nBeagleUpdateTransitionsCalls;
            
            bool _useAutoScaling;
            
            std::vector<int> _matrixUpdateIndices;
            std::vector<double> _branchLengths;
            std::vector<int> _scaleBufferIndices;
            std::vector<BeagleOperation> _operations;
            std::vector<int> _upperIndices;
            
        };
    }
}

#endif /* BeagleFlexibleTreeLikleihood_hpp */
