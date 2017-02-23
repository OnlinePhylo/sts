//
//  SimpleFlexibleTreeLikelihood.hpp
//  sts
//
//  Created by Mathieu Fourment on 23/06/2016.
//  Copyright © 2016 Mathieu Fourment. All rights reserved.
//

#ifndef SimpleFlexibleTreeLikelihood_hpp
#define SimpleFlexibleTreeLikelihood_hpp

#include <stdio.h>
#include <vector>

#include "AbstractFlexibleTreeLikelihood.hpp"

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace sts {
    namespace online {
    
        class SimpleFlexibleTreeLikelihood : public AbstractFlexibleTreeLikelihood{
            
        public:
            SimpleFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, bpp::SubstitutionModel &model, bpp::DiscreteDistribution& rateDist);
            
            virtual ~SimpleFlexibleTreeLikelihood(){}
            
            virtual double calculateLogLikelihood();
            
            virtual double calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength);
            
            virtual void calculateDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);

            virtual void calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);
            
            void updateNode(const bpp::Node& node);
            
            void updateAllNodes();
    
        protected:
            
            void setStates(const bpp::SiteContainer& sites);
            
            bool traverse(const bpp::Node* node);
            
            void calculatePatternLikelihood( const double *partials, const double *frequencies, double *outLogLikelihoods)const;
            
            void integratePartials( const double *inPartials, const double *proportions, double *outPartials )const;
            
            void updatePartialsKnownKnown( const int *states1, const double *matrices1, const int *states2, const double *matrices2, double *partials )const;
            
            void updatePartialsKnownUndefined( const int *states1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 )const;
            
            void updatePartialsUndefinedUndefined( const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 )const;
            
            void updatePartials( int nodeIndex1, int nodeIndex2, const double* matrix1, const double *matrix2, double *partials );
            
            void updatePartials( int nodeIndex1, int nodeIndex2, int nodeIndex3 );
            
            void updateUpperPartialsKnown( const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const int *states, double *partials ) const;
            
            void updateUpperPartialsUndefined( const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ) const;
            
            void updateUpperPartialsRootUndefined( const double *partials1, const double *matrices1, const double *frequencies, double *partials_upper )const;
            
            void updateUpperPartialsRootKnown(const double *matrices1, const int *states1, const double *frequencies, double *partials_upper ) const;
            
            void updateUpperPartials(const bpp::Node* node);
            
            void calculatePatternLikelihood( const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ) const;
            
            void calculatePatternLikelihoodLeaf( const double *partials_upper, const int* states, const double *matrix_lower, const double *proportions, double *pattern_lk ) const;
            
            
        private:
            void calculateDerivatives(const bpp::Node& distal, std::string taxonName, int index, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2);
            
            std::vector<std::vector<int> > _states;

            double _logLnl;
            
            int _matrixSize;
            
            std::vector<double> _patternWeights;
            
            std::vector<std::vector<double> > _matrices;
            std::vector<std::vector<double> > _partials;
            std::vector<double> _rootPartials;
            std::vector<double> _patternLikelihoods;
            std::vector<std::vector<double> > _upperPartials;
            
            std::vector<double> _temporaryPartials;
            std::vector<std::vector<double> > _temporaryMatrices;
            std::vector<std::vector<double> > _temporaryPatternLikelihood;
        };
    }
}
#endif /* SimpleFlexibleTreeLikelihood_hpp */
