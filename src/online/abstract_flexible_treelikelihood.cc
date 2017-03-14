#include "abstract_flexible_treelikelihood.h"

using namespace std;
using namespace bpp;

namespace sts {
    namespace online {
        
        size_t AbstractFlexibleTreeLikelihood::operationCallCount = 0;
        
        AbstractFlexibleTreeLikelihood::AbstractFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bool useAmbiguities):
			_model(&model), _rateDist(&rateDist), _tree(nullptr), _useAmbiguities(useAmbiguities){
            _stateCount = model.getNumberOfStates();
            std::unique_ptr<bpp::SiteContainer> sc(patterns.getSites());
            _taxa = sc->getSequencesNames();
            _sequenceCount = _taxa.size();
            _totalNodeCount = (_sequenceCount * 2) - 1;
            _rateCount = rateDist.getNumberOfCategories();
            _patternCount = patterns.getWeights().size();
            
            _useScaleFactors = false;
            _recomputeScaleFactors = false;
            
            _updateSiteModel = true;
            _updateSubstitutionModel = true;
            _needNodeUpdate.assign(_totalNodeCount, true);
            
            _partialCount = _totalNodeCount*2+1; // includes upperPartials
            _matrixCount = _totalNodeCount + 3; // temporary matrices
            
            _updatePartials = true;
            _updateUpperPartials = true;
        }
        
        
        void AbstractFlexibleTreeLikelihood::initialize(const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bpp::TreeTemplate<bpp::Node>& tree){
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
