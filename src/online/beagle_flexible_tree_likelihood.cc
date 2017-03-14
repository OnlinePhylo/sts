#include "beagle_flexible_tree_likelihood.h"


#include "bpp_shim.h"
#include "util.h"

using sts::util::beagle_check;
using namespace std;

namespace sts {
    namespace online {
        
        BeagleFlexibleTreeLikelihood::BeagleFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, const bpp::SubstitutionModel& model, const bpp::DiscreteDistribution& rateDist, bool useAmbiguities ):
            AbstractFlexibleTreeLikelihood(patterns, model, rateDist){
                
            _matrixCount = _totalNodeCount + 5; // temporary matrices: pendant, proxiximal, distal + 2 derivatives
            _partialCount = _totalNodeCount*2+3; // includes upperPartials, 1 for attachment, 2 for derivatives
            const int scaleFactorCount = _sequenceCount + (_sequenceCount*2-2);
            
            _beagleInstance = beagleCreateInstance(
                                                   0,              // Number of tip data elements (input)
                                                   _partialCount,  // Number of partials buffers to create (input)
                                                   0,              // Number of compact state representation buffers to create (input)
                                                   _stateCount,    // Number of states in the continuous-time Markov chain (input)
                                                   _patternCount,  // Number of site patterns to be handled by the instance (input)
                                                   1,              // Number of rate matrix eigen-decomposition buffers to allocate (input)
                                                   _matrixCount,   // Number of rate matrix buffers (input)
                                                   _rateCount,     // Number of rate categories (input)
                                                   scaleFactorCount,  // Number of scaling buffers (one for each internal node + 1 for accumulation) - 1 extra buffer for prox, distal
                                                   NULL,           // List of potential resource on which this instance is allowed (input, NULL implies no
                                                   // restriction
                                                   0,              // Length of resourceList list (input)
                                                   BEAGLE_FLAG_VECTOR_SSE | BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_SCALING_MANUAL, // Bit-flags indicating
                                                   //preferred implementation charactertistics, see BeagleFlags (input)
                                                   0,              // Bit-flags indicating required implementation characteristics, see BeagleFlags (input)
                                                   &_instanceDetails);
            
            if(_beagleInstance < 0)
                beagle_check(_beagleInstance);
            
			const std::vector<unsigned int> weights = patterns.getWeights();
            std::unique_ptr<const bpp::SiteContainer> sites(patterns.getSites());
			setPartials(*sites);
                
			std::vector<double> w;
			w.reserve(_patternCount);
			auto castit = [](unsigned int w) { return static_cast<double>(w); };
			std::transform(weights.begin(), weights.end(), w.begin(), castit);
			beagle_check(beagleSetPatternWeights(_beagleInstance, w.data()));
            
            _upperPartialsIndexes.resize(_totalNodeCount);
            _logLnl = 0;
        }
        
        
        BeagleFlexibleTreeLikelihood::~BeagleFlexibleTreeLikelihood(){
            if(_beagleInstance >= 0)
                beagleFinalizeInstance(_beagleInstance);
        }
        
        void BeagleFlexibleTreeLikelihood::initialize(const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bpp::TreeTemplate<bpp::Node>& tree){
            AbstractFlexibleTreeLikelihood::initialize(model, rateDist, tree);
            
            _scaleBufferIndices.resize(_internalNodeCount+1); // internal node + 1 for accumulation
            _scaleBufferUpperIndices.resize(_nodeCount-1); // every node except root and its kids + 1 for new node + 1 for accummulation
            
            updateSiteModel();
            updateSubstitutionModel();
        }
        
        void BeagleFlexibleTreeLikelihood::setStates(const bpp::SiteContainer& sites){
            
        }

        void BeagleFlexibleTreeLikelihood::setPartials(const bpp::SiteContainer& sites){
            std::vector<double> partials(_patternCount * _stateCount * _rateCount);
            
            for(int i = 0; i < _sequenceCount; i++){
                const bpp::Sequence& sequence = sites.getSequence(i);
                for(size_t site = 0; site < _patternCount; site++) {
                    for(size_t j = 0; j < _stateCount; j++) {
                        size_t idx = _stateCount * site + j;
                        partials[idx] = _model->getInitValue(j, sequence.getValue(site));
                    }
                }
                
                size_t per_category = _patternCount * _stateCount;
                for(size_t c = 1; c < _rateCount; c++) {
                    std::copy(partials.begin(),
                              partials.begin() + per_category,
                              partials.begin() + (per_category * c));
                }
                
                beagle_check(beagleSetPartials(_beagleInstance, i, partials.data()));
            }
        }
        
        void BeagleFlexibleTreeLikelihood::traverseUpper(const bpp::Node* node){
            
            if(node->hasFather()){
                const bpp::Node* parent = node->getFather();
                const bpp::Node* sibling = parent->getSon(0);
                if (sibling == node) {
                    sibling = parent->getSon(1);
                }
                
                if(parent->hasFather()){
                    const bpp::Node* grandParent = parent->getFather();
                    const int idSibling = sibling->getId();
                    
                    int idMatrix = parent->getId();
                    
                    // The sons of the right node of the root are going to use the lower partials of the left node of the root
                    if(!grandParent->hasFather() && grandParent->getSon(1) == parent){
                        idMatrix = grandParent->getSon(0)->getId();
                    }
                    
                    _upperPartialsIndexes[node->getId()] = node->getId()+_totalNodeCount;

                    _operations.push_back(BeagleOperation(
                                                          {_upperPartialsIndexes[node->getId()],      // Destination buffer
                                                              BEAGLE_OP_NONE,    // (output) scaling buffer index
                                                              BEAGLE_OP_NONE,    // (input) scaling buffer index
                                                              _upperPartialsIndexes[parent->getId()],     // Index of first child partials buffer
                                                              idMatrix,          // Index of first child transition matrix
                                                              idSibling,         // Index of second child partials buffer
                                                              idSibling}));
                    if (_useScaleFactors) {
                        // get the index of this scaling buffer
                        const int indexInternal = node->getId() - _leafCount + _internalNodeCount+1;
                        _scaleBufferUpperIndices[indexInternal] = indexInternal;
                        
                        if (_recomputeScaleFactors) {
                            
                            // store the index
                            _operations.back().destinationScaleWrite = _scaleBufferUpperIndices[indexInternal];
                            _operations.back().destinationScaleRead  = BEAGLE_OP_NONE;
                            
                        }
                        else {
                            _operations.back().destinationScaleWrite  = BEAGLE_OP_NONE;
                            _operations.back().destinationScaleRead = _scaleBufferUpperIndices[indexInternal];// Read existing scaleFactor
                        }
                    }
                }
                // We dont need to calculate upper partials for the left son of the root as it is using the lower partials of its sibling
                // Left node of the root (function should never be called on theright node of the root)
                else if(parent->getSon(0) == node){
                    _upperPartialsIndexes[node->getId()] = parent->getSon(1)->getId();
                    if (_useScaleFactors) {
                        const int indexInternal = node->getId() - _leafCount + _internalNodeCount+1;
                        _scaleBufferUpperIndices[indexInternal] = parent->getSon(1)->getId() - _sequenceCount;
                    }
                }
                else{
                    _upperPartialsIndexes[node->getId()] = parent->getSon(0)->getId();
                    if (_useScaleFactors) {
                        const int indexInternal = node->getId() - _leafCount + _internalNodeCount+1;
                        _scaleBufferUpperIndices[indexInternal] = parent->getSon(0)->getId() - _sequenceCount;
                    }
                }
            }
            
            if(node->getNumberOfSons() > 0){
                traverseUpper(node->getSon(0));
                traverseUpper(node->getSon(1));
            }
        }

        
        bool BeagleFlexibleTreeLikelihood::traverse(const bpp::Node* node){
            bool update = _needNodeUpdate[node->getId()];
            
            if(node->hasFather() && update){
                _matrixUpdateIndices.push_back(node->getId());
                _branchLengths.push_back(node->getDistanceToFather());
            }
            
            if(node->getNumberOfSons() > 0){
                
                bool update1 = traverse(node->getSon(0));
                bool update2 = traverse(node->getSon(1));
                
                if( update1 || update2 ){
                    
                    int id1 = node->getSon(0)->getId();
                    int id2 = node->getSon(1)->getId();
        
                    _operations.push_back(BeagleOperation(
                                                          {node->getId(),      // Destination buffer
                                                              BEAGLE_OP_NONE,    // (output) scaling buffer index
                                                              BEAGLE_OP_NONE,    // (input) scaling buffer index
                                                              id1,               // Index of first child partials buffer
                                                              id1,               // Index of first child transition matrix
                                                              id2,               // Index of second child partials buffer
                                                              id2}));            // Index of second child transition matrix
                    
                    if (_useScaleFactors) {
                        // get the index of this scaling buffer
                        int indexInternal = node->getId() - _leafCount;
                        
                        if (_recomputeScaleFactors) {
                            // store the index
                            _scaleBufferIndices[indexInternal] = indexInternal;
                            _operations.back().destinationScaleWrite = _scaleBufferIndices[indexInternal];
                            _operations.back().destinationScaleRead  = BEAGLE_OP_NONE;
                        }
                        else {
                            _operations.back().destinationScaleWrite  = BEAGLE_OP_NONE;
                            _operations.back().destinationScaleRead = _scaleBufferIndices[indexInternal];// Read existing scaleFactor
                        }
                    }
                    
                    update |= (update1 | update2);
                }
            }
            return update;
        }
        
        double BeagleFlexibleTreeLikelihood::calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength){
            
            _branchLengths.clear();
            _matrixUpdateIndices.clear();
            _operations.clear();
            
            if(_updatePartials){
                traverse(_tree->getRootNode());
                traverseUpper(_tree->getRootNode());
                _updatePartials = _updateUpperPartials = false;
                _needNodeUpdate.assign(_totalNodeCount, false);
            }
            else if(_updateUpperPartials){
            	traverseUpper(_tree->getRootNode());
                _updateUpperPartials = false;
            }
            
            const int category_weight_index = 0;
            const int state_frequency_index = 0;
            int cumulateScaleBufferIndex = BEAGLE_OP_NONE;

            const int distalIndex = distal.getId();
            const int indexTaxon = std::find(_taxa.begin(), _taxa.end(), taxonName) - _taxa.begin();
            const int tmpPartialsIndex = _totalNodeCount*2;
            
            // Temporary transition matrices
            const int tempMatrixPendant = _totalNodeCount;
            const int tempMatrixDistal = _totalNodeCount + 1;
            const int tempMatrixProximal =  _totalNodeCount + 2;
            
            // Distal upper partial index
            const int distalUpperIndex = _upperPartialsIndexes[distalIndex];
            
            _branchLengths.push_back(pendantLength);
            _matrixUpdateIndices.push_back(tempMatrixPendant);
            
            _branchLengths.push_back(distalLength);
            _matrixUpdateIndices.push_back(tempMatrixDistal);
            
            _branchLengths.push_back(proximalLength);
            _matrixUpdateIndices.push_back(tempMatrixProximal);
            
            // Update distal, pendant and proximal transition matrices
            beagle_check(beagleUpdateTransitionMatrices(_beagleInstance,               // instance
                                                        0,                             // eigenIndex
                                                        _matrixUpdateIndices.data(),   // probabilityIndices
                                                        NULL,                          // firstDerivativeIndices
                                                        NULL,                          // secondDerivativeIndices
                                                        _branchLengths.data(),         // edgeLengths
                                                        _matrixUpdateIndices.size())); // count

            // Update partials at ancestral node of distal and proximal                                                            
			_operations.push_back(BeagleOperation({tmpPartialsIndex,       // Destination buffer
													BEAGLE_OP_NONE,        // (output) scaling buffer index
                                                    BEAGLE_OP_NONE,        // (input) scaling buffer index
                                                    distalIndex,           // Index of first child partials buffer
                                                    tempMatrixDistal,      // Index of first child transition matrix
                                                    distalUpperIndex,      // Index of second child partials buffer
                                                    tempMatrixProximal})); // Index of second child transition matrix
            
            if (_useScaleFactors) {
                // get the index of this scaling buffer
                int indexInternal = _scaleBufferUpperIndices.size()-2;
                
                if (_recomputeScaleFactors) {
                    
                    // store the index
                    _scaleBufferUpperIndices[indexInternal] = indexInternal;
                    _operations.back().destinationScaleWrite = _scaleBufferUpperIndices[indexInternal];
                    _operations.back().destinationScaleRead  = BEAGLE_OP_NONE;
                    
                }
                else {
                    _operations.back().destinationScaleWrite  = BEAGLE_OP_NONE;
                    _operations.back().destinationScaleRead = _scaleBufferUpperIndices[indexInternal];// Read existing scaleFactor
                }
            }
            
            beagle_check(beagleUpdatePartials(_beagleInstance,
                                              _operations.data(),
                                              _operations.size(),
                                              BEAGLE_OP_NONE));
            
            operationCallCount += _operations.size();

            double logLnl = 0;
            
            for(int pass = 0; pass < 2; pass++){
                if (_useScaleFactors) {
                    if (_recomputeScaleFactors) {
                        cumulateScaleBufferIndex = _scaleBufferUpperIndices.size()-1;
                        beagle_check(beagleResetScaleFactors(_beagleInstance, cumulateScaleBufferIndex));

                        beagle_check(beagleAccumulateScaleFactors(_beagleInstance,
                                                                  _scaleBufferIndices.data()+_internalNodeCount+1,
                                                                  _nodeCount-2,
                                                                  cumulateScaleBufferIndex));
                    }
                    else {
                        cumulateScaleBufferIndex = _scaleBufferUpperIndices.size()-1;
                    }
                }
                
                logLnl = 0;
                
                // Calculate log likelihood                         
                beagle_check(beagleCalculateEdgeLogLikelihoods(_beagleInstance,
                                          &tmpPartialsIndex,
                                          &indexTaxon,
                                          &tempMatrixPendant,
                                          NULL,
                                          NULL,
                                          &category_weight_index,
                                          &state_frequency_index,
                                          &cumulateScaleBufferIndex,
                                          1,
                                          &logLnl,
                                          NULL,
                                          NULL));
                
                if (std::isnan(_logLnl) || std::isinf(_logLnl)) {
                    _useScaleFactors = true;
                    _recomputeScaleFactors = true;
                    
                    traverse(_tree->getRootNode());
                }
                else{
                    break;
                }
            }
            
            return logLnl;
        }
        
        double BeagleFlexibleTreeLikelihood::calculateLogLikelihood(){
            
            if(!_updatePartials)return _logLnl;
            
            _branchLengths.clear();
            _matrixUpdateIndices.clear();
            _operations.clear();
            
            const int category_weight_index = 0;
            const int state_frequency_index = 0;
            
            int rootIndex = _tree->getRootNode()->getId();
            traverse(_tree->getRootNode());
            
            _logLnl = 0;
            
            if (_updateSubstitutionModel) {
                updateSubstitutionModel();
            }
            
            if (_updateSiteModel) {
                updateSiteModel();
            }
            
            // Register topology, branch lengths; update transition matrices.
            beagle_check(beagleUpdateTransitionMatrices(_beagleInstance,               // instance
                                                        0,                             // eigenIndex
                                                        _matrixUpdateIndices.data(),   // probabilityIndices
                                                        NULL,                          // firstDerivativeIndices
                                                        NULL,                          // secondDerivativeIndices
                                                        _branchLengths.data(),         // edgeLengths
                                                        _matrixUpdateIndices.size())); // count
            
            for(int pass = 0; pass < 2; pass++){
                // Update partials for all traversed nodes
                beagle_check(beagleUpdatePartials(_beagleInstance,
                                                  _operations.data(),
                                                  _operations.size(),
                                                  BEAGLE_OP_NONE));
                
                operationCallCount += _operations.size();
                
                int cumulateScaleBufferIndex = BEAGLE_OP_NONE;
                if (_useScaleFactors) {
                    cumulateScaleBufferIndex = _internalNodeCount;
                    if (_recomputeScaleFactors) {
                        beagle_check(beagleResetScaleFactors(_beagleInstance, cumulateScaleBufferIndex));
                        
                        beagle_check(beagleAccumulateScaleFactors(_beagleInstance,
                                                                  _scaleBufferIndices.data(),
                                                                  _internalNodeCount,
                                                                  cumulateScaleBufferIndex));
                    }
                }
                
                beagle_check(beagleCalculateRootLogLikelihoods(_beagleInstance,
                                                               &rootIndex,
                                                               &category_weight_index,
                                                               &state_frequency_index,
                                                               &cumulateScaleBufferIndex,
                                                               1,  // op count
                                                               &_logLnl));
                
                if (std::isnan(_logLnl) || std::isinf(_logLnl)) {
                    _useScaleFactors = true;
                    _recomputeScaleFactors = true;
                    _operations.clear();
                    
                    traverse(_tree->getRootNode());
                }
                else{
                    break;
                }
                
            }
            
            _needNodeUpdate.assign(_totalNodeCount, false);
            _updatePartials = false;
            _updateUpperPartials = true;
            
            return _logLnl;
        }
        
        // Compute derivatives of pendant branch with taxon taxonName
        void BeagleFlexibleTreeLikelihood::calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
            
            _branchLengths.clear();
            _matrixUpdateIndices.clear();
            _operations.clear();
            
            if(_updatePartials){
                traverse(_tree->getRootNode());
                traverseUpper(_tree->getRootNode());
                _updatePartials = _updateUpperPartials = false;
                _needNodeUpdate.assign(_totalNodeCount, false);
            }
            else if(_updateUpperPartials){
                traverseUpper(_tree->getRootNode());
                _updateUpperPartials = false;
            }
            
            const int distalIndex = distal.getId();
            const int indexTaxon = std::find(_taxa.begin(), _taxa.end(), taxonName) - _taxa.begin();
            const int tmpPartialsIndex = _totalNodeCount*2;
            
            // Temporary transition matrices
            const int tempMatrixPendant = _totalNodeCount;
            const int tempMatrixDistal = _totalNodeCount + 1;
            const int tempMatrixProximal =  _totalNodeCount + 2;
            const int tempMatrixDistald1 =  _totalNodeCount + 3;
            const int tempMatrixDistald2 =  _totalNodeCount + 4;
            
            // Distal upper partial index
            const int distalUpperIndex = _upperPartialsIndexes[distalIndex];
            
            // calculate partial for proximal and pendant
            _operations.push_back(BeagleOperation({tmpPartialsIndex,
                BEAGLE_OP_NONE,
                BEAGLE_OP_NONE,
                indexTaxon,
                tempMatrixPendant,
                distalUpperIndex,
                tempMatrixProximal}));
            
            _branchLengths.push_back(pendantLength);
            _matrixUpdateIndices.push_back(tempMatrixPendant);
            
            _branchLengths.push_back(proximalLength);
            _matrixUpdateIndices.push_back(tempMatrixProximal);

            // update matrices of proximal and pendant
            beagle_check(beagleUpdateTransitionMatrices(_beagleInstance,
                                                        0,
                                                        _matrixUpdateIndices.data(),
                                                        NULL,
                                                        NULL,
                                                        _branchLengths.data(),
                                                        _matrixUpdateIndices.size()));

            // update matrix and derivatives of distal
            beagle_check(beagleUpdateTransitionMatrices(_beagleInstance,
                                                        0,
                                                        &tempMatrixDistal,
                                                        &tempMatrixDistald1,
                                                        &tempMatrixDistald2,
                                                        &distalLength,
                                                        1));

            beagle_check(beagleUpdatePartials(_beagleInstance, _operations.data(), _operations.size(), NULL));

            operationCallCount += _operations.size();

            const int category_weight_index = 0;
            const int state_frequency_index = 0;
            int cumulateScaleBufferIndex = BEAGLE_OP_NONE;
            double logLike=0;
            double dd1=0,dd2=0;
            
            beagle_check(beagleCalculateEdgeLogLikelihoods(_beagleInstance,
                                                           &tmpPartialsIndex,
                                                           &distalIndex,
                                                           &tempMatrixDistal,
                                                           &tempMatrixDistald1,
                                                           &tempMatrixDistald2,
                                                           &category_weight_index,
                                                           &state_frequency_index,
                                                           &cumulateScaleBufferIndex,
                                                           1,
                                                           &logLike,
                                                           &dd1,
                                                           &dd2));

            *d1 = dd1;
            *d2 = dd2;
        }
        
        void BeagleFlexibleTreeLikelihood::calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
            
            _branchLengths.clear();
            _matrixUpdateIndices.clear();
            _operations.clear();
            
            if(_updatePartials){
                traverse(_tree->getRootNode());
                traverseUpper(_tree->getRootNode());
                _updatePartials = _updateUpperPartials = false;
                _needNodeUpdate.assign(_totalNodeCount, false);
            }
            else if(_updateUpperPartials){
                traverseUpper(_tree->getRootNode());
                _updateUpperPartials = false;
            }
            
            const int distalIndex = distal.getId();
            const int indexTaxon = std::find(_taxa.begin(), _taxa.end(), taxonName) - _taxa.begin();
            const int tmpPartialsIndex = _totalNodeCount*2;
            
            // Temporary transition matrices
            const int tempMatrixPendant = _totalNodeCount;
            const int tempMatrixDistal = _totalNodeCount + 1;
            const int tempMatrixProximal =  _totalNodeCount + 2;
            const int tempMatrixPendantd1 =  _totalNodeCount + 3;
            const int tempMatrixPendantd2 =  _totalNodeCount + 4;
            
            // Distal upper partial index
            const int distalUpperIndex = _upperPartialsIndexes[distalIndex];
            
            // calculate partial for proximal and distal
            _operations.push_back(BeagleOperation({tmpPartialsIndex,
                BEAGLE_OP_NONE,
                BEAGLE_OP_NONE,
                distalIndex,
                tempMatrixDistal,
                distalUpperIndex,
                tempMatrixProximal}));
            
            _branchLengths.push_back(distalLength);
            _matrixUpdateIndices.push_back(tempMatrixDistal);
            
            _branchLengths.push_back(proximalLength);
            _matrixUpdateIndices.push_back(tempMatrixProximal);
            
            // update matrices of proximal and pendant
            beagle_check(beagleUpdateTransitionMatrices(_beagleInstance,
                                                        0,
                                                        _matrixUpdateIndices.data(),
                                                        NULL,
                                                        NULL,
                                                        _branchLengths.data(),
                                                        _matrixUpdateIndices.size()));
            
            // update matrix and derivatives of pendant
            beagle_check(beagleUpdateTransitionMatrices(_beagleInstance,
                                                        0,
                                                        &tempMatrixPendant,
                                                        &tempMatrixPendantd1,
                                                        &tempMatrixPendantd2,
                                                        &pendantLength,
                                                        1));
            
            beagle_check(beagleUpdatePartials(_beagleInstance, _operations.data(), _operations.size(), NULL));
            
            operationCallCount += _operations.size();
            
            const int category_weight_index = 0;
            const int state_frequency_index = 0;
            int cumulateScaleBufferIndex = BEAGLE_OP_NONE;
            double logLike=0;
            double dd1=0,dd2=0;
            
            beagle_check(beagleCalculateEdgeLogLikelihoods(_beagleInstance,
                                                           &tmpPartialsIndex,
                                                           &distalIndex,
                                                           &tempMatrixPendant,
                                                           &tempMatrixPendantd1,
                                                           &tempMatrixPendantd2,
                                                           &category_weight_index,
                                                           &state_frequency_index,
                                                           &cumulateScaleBufferIndex,
                                                           1,
                                                           &logLike,
                                                           &dd1,
                                                           &dd2));
            
            *d1 = dd1;
            *d2 = dd2;
        }
        
        void BeagleFlexibleTreeLikelihood::updateSiteModel(){
            const std::vector<double>& categories = _rateDist->getCategories();
            const std::vector<double>& weights = _rateDist->getProbabilities();
        
            beagle_check(beagleSetCategoryRates(_beagleInstance, categories.data()));
            beagle_check(beagleSetCategoryWeights(_beagleInstance, 0, weights.data()));
            
            _updateSiteModel = false;
            _needNodeUpdate.assign(_totalNodeCount, true);
        }
        
        void BeagleFlexibleTreeLikelihood::updateSubstitutionModel(){
            std::vector<double> evec(_stateCount * _stateCount),
            ivec(_stateCount * _stateCount);
            
            const vector<double> eval = _model->getEigenValues();
            const bpp::Matrix<double>& lvec = _model->getRowLeftEigenVectors();
            const bpp::Matrix<double>& rvec = _model->getColumnRightEigenVectors();
            
            for(int i = 0; i < _stateCount; i++){
                const vector<double>& lvecr = lvec.row(i);
                const vector<double>& rvecr = rvec.row(i);
                std::copy(lvecr.begin(), lvecr.end(), ivec.begin()+i*_stateCount);
                std::copy(rvecr.begin(), rvecr.end(), evec.begin()+i*_stateCount);
            }

            // Eigen decomposition
            beagle_check(beagleSetEigenDecomposition(_beagleInstance, 0, evec.data(), ivec.data(), eval.data()));
            // State frequencies
            beagle_check(beagleSetStateFrequencies(_beagleInstance, 0, _model->getFrequencies().data()));
            
            _updateSubstitutionModel = false;
            _needNodeUpdate.assign(_totalNodeCount, true);
        }
    }
}
