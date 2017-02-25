//
//  AbstractFlexibleTreeLikelihood.cpp
//  sts
//
//  Created by Mathieu Fourment on 22/02/2017.
//  Copyright © 2017 Mathieu Fourment. All rights reserved.
//

#include "AbstractFlexibleTreeLikelihood.hpp"

using namespace std;
using namespace bpp;

namespace sts {
    namespace online {
        
        AbstractFlexibleTreeLikelihood::AbstractFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist):
        _patterns(patterns), _model(&model), _rateDist(&rateDist), _tree(nullptr){
            _stateCount = model.getNumberOfStates();
            std::unique_ptr<bpp::SiteContainer> sc(_patterns.getSites());
            _taxa = sc->getSequencesNames();
            _sequenceCount = _taxa.size();
            _totalNodeCount = (_sequenceCount * 2) - 1;
            _rateCount = rateDist.getNumberOfCategories();
            
            _patternCount = _patterns.getWeights().size();
            
            _useScaleFactors = false;
            _recomputeScaleFactors = false;
            
            _updateSiteModel = true;
            _updateSubstitutionModel = true;
            _needNodeUpdate.assign(_totalNodeCount, true);
            
            _partialCount = _totalNodeCount*2; // includes upperPartials
            _matrixCount = _totalNodeCount + 2; // temporary matrices
            
            _updatePartials = true;
            _updateUpperPartials = true;
        }
        
        void AbstractFlexibleTreeLikelihood::initialize(bpp::TreeTemplate<bpp::Node>& tree, const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist){
            _tree = &tree;
            _model = &model;
            _rateDist = &rateDist;
            
            
            _nodeCount = _tree->getNumberOfNodes();
            _leafCount = _tree->getNumberOfLeaves();
            _internalNodeCount = _nodeCount - _leafCount;
            
            _needNodeUpdate.assign(_totalNodeCount, true);
            _updatePartials = true;
            _updateUpperPartials = true;
        }
    }
}
