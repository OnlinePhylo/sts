#include "simple_flexible_tree_likelihood.h"


using namespace std;

namespace sts {
    namespace online {
        
        SimpleFlexibleTreeLikelihood::SimpleFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist, bool useAmbiguities):
        AbstractFlexibleTreeLikelihood(patterns, model, rateDist, useAmbiguities){

            _matrixSize = _stateCount*_stateCount;
            
            _patternWeights.resize(_patternCount);
            const std::vector<unsigned int>& w = patterns.getWeights();
            auto castit = [](unsigned int w) { return static_cast<double>(w); };
            std::transform(w.begin(), w.end(), _patternWeights.begin(), castit);
            _states.resize(_totalNodeCount);
            for(auto it = _states.begin(); it != _states.end(); ++it){
                it->resize(_patternCount);
            }
            std::unique_ptr<bpp::SiteContainer> sites(patterns.getSites());
            
            // Lower partials 
            _partials.resize(_totalNodeCount*2+1);
            for(auto it = _partials.begin(); it != _partials.end(); ++it){
                it->resize(_rateCount*_stateCount*_patternCount);
            }
            
            // Probability matrices
            _matrices.resize(_totalNodeCount+3);
            for(auto it = _matrices.begin(); it != _matrices.end(); ++it){
                it->resize(_rateCount*_matrixSize);
            }
            
            // Rates are integrated out
            _rootPartials.resize(4);
            for(auto it = _rootPartials.begin(); it != _rootPartials.end(); ++it){
                it->assign(_stateCount*_patternCount, 0.);
            }
            
            // 0 current tree
            // 1 current tree + taxon
            // 2 1st derivative current tree + taxon
            // 3 2nd derivative current tree + taxon
            _patternLikelihoods.resize(4);
            for(auto it = _patternLikelihoods.begin(); it != _patternLikelihoods.end(); ++it){
                it->assign(_patternCount, 0.);
            }
            
            if(!_useAmbiguities){
	            setStates(*sites);
	        }
	        else{
	        	setPartials(*sites);
	        }
            
            _upperPartialsIndexes.resize(_totalNodeCount);
	        
            _updatePartials = true;
            _updateUpperPartials = true;
        }
        
        
        void SimpleFlexibleTreeLikelihood::setStates(const bpp::SiteContainer& sites){
            
            for(int i = 0; i < _sequenceCount; i++){
                const bpp::Sequence& sequence = sites.getSequence(i);
                for(size_t site = 0; site < _patternCount; site++) {
                    _states[i][site] = sequence.getValue(site);
                }
            }
        }
        
        void SimpleFlexibleTreeLikelihood::setPartials(const bpp::SiteContainer& sites){
            for(int i = 0; i < _sequenceCount; i++){
                const bpp::Sequence& sequence = sites.getSequence(i);
                for(size_t site = 0; site < _patternCount; site++) {
                    for(size_t j = 0; j < _stateCount; j++) {
                        size_t idx = _stateCount * site + j;
                        _partials[i][idx] = _model->getInitValue(j, sequence.getValue(site));
                    }
                    for(size_t c = 1; c < _rateCount; c++)
                        std::copy(_partials[i].begin(),
                                  _partials[i].begin() + _stateCount*_patternCount,
                                  _partials[i].begin() + (_stateCount*_patternCount * c));
                }
            }
        }
        
        void SimpleFlexibleTreeLikelihood::traverseUpper(const bpp::Node* node){
        
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
                    updatePartials(_upperPartialsIndexes[node->getId()], _upperPartialsIndexes[parent->getId()], idMatrix, idSibling, idSibling);
				}
                // We dont need to calculate upper partials for the children of the root as it is using the lower partials of its sibling
                // Left node of the root
                else if(parent->getSon(0) == node){
                    _upperPartialsIndexes[node->getId()] = parent->getSon(1)->getId();
                }
                else{
                    _upperPartialsIndexes[node->getId()] = parent->getSon(0)->getId();
                }
        	}

        	if(node->getNumberOfSons() > 0){
				traverseUpper(node->getSon(0));
                traverseUpper(node->getSon(1));
        	}
        }
        
        bool SimpleFlexibleTreeLikelihood::traverse(const bpp::Node* node){
            bool update = _needNodeUpdate[node->getId()];
            
            if(node->hasFather() && update){
                int id = node->getId();
                
                int offset = 0;
                for(int c = 0; c < _rateCount; c++){
                    const bpp::Matrix<double>& m = _model->getPij_t(node->getDistanceToFather() * _rateDist->getCategory(c));
                    
                    for(int i = 0; i < _stateCount; i++){
                        const vector<double>& mi = m.row(i);
                        std::copy(mi.begin(), mi.end(), _matrices[id].data()+offset);
                        offset += _stateCount;
                    }
                }
            }
            
            if(node->getNumberOfSons() > 0){
                
                bool update1 = traverse(node->getSon(0));
                bool update2 = traverse(node->getSon(1));
                
                if( update1 || update2 ){
                    
                    updatePartials(node->getId(), node->getSon(0)->getId(), node->getSon(0)->getId(), node->getSon(1)->getId(), node->getSon(1)->getId());

//                    if (_useScaleFactors) {
//                        // get the index of this scaling buffer
//                        int indexInternal = node->getId() - _leafCount;
//                        
//                        if (_recomputeScaleFactors) {
//                            
//                            // store the index
//                            _scaleBufferIndices[indexInternal] = indexInternal;
//                            _operations.back().destinationScaleWrite = _scaleBufferIndices[indexInternal];
//                            _operations.back().destinationScaleRead  = BEAGLE_OP_NONE;
//                            
//                        }
//                        else {
//                            _operations.back().destinationScaleWrite  = BEAGLE_OP_NONE;
//                            _operations.back().destinationScaleRead = _scaleBufferIndices[indexInternal];// Read existing scaleFactor
//                        }
//                        
//                    }
//                    else if (_useAutoScaling) {
//                        _scaleBufferIndices[node->getId() - _leafCount] = node->getId();
//                    }
                    
                    update |= (update1 | update2);
                }
            }
            return update;
        }
        
        void SimpleFlexibleTreeLikelihood::calculatePatternLikelihood( const double *partials, const double *frequencies, double *outLogLikelihoods)const{
            int v = 0;
            int i = 0;
            
            for ( int k = 0; k < _patternCount; k++ ) {
                
                outLogLikelihoods[k] = 0;
                for ( i = 0; i < _stateCount; i++ ) {
                    outLogLikelihoods[k] += frequencies[i] * partials[v];
                    v++;
                }
                //outLogLikelihoods[k] = log(outLogLikelihoods[k]);
                
                if ( _useScaleFactors ) {
                    //outLogLikelihoods[k] += getLogScalingFactor( tlk, k);
                }
            }
            
        }
        

        void SimpleFlexibleTreeLikelihood::integratePartials( const double *inPartials, const double *proportions, double *outPartials )const{
            int i,k;
            double *pPartials = outPartials;
            const double *pInPartials = inPartials;
            
            if( _rateCount == 1 ){
                memcpy(outPartials, inPartials, sizeof(double)*_patternCount*_stateCount);
                return;
            }
            
            for ( k = 0; k < _patternCount; k++ ) {
                
                for ( i = 0; i < _stateCount; i++ ) {
                    
                    *pPartials++ = *pInPartials++ * proportions[0];
                }
            }
            
            
            for ( int l = 1; l < _rateCount; l++ ) {
                pPartials = outPartials;
                
                for ( k = 0; k < _patternCount; k++ ) {
                    
                    for ( i = 0; i < _stateCount; i++ ) {
                        
                        *pPartials += *pInPartials++ * proportions[l];
                        pPartials++;
                    }
                }
            }
            
        }
        
        void SimpleFlexibleTreeLikelihood::updatePartialsKnownKnown( const int *states1, const double *matrices1, const int *states2, const double *matrices2, double *partials )const{
            int k,i,w;
            int u = 0;
            int state1, state2;
            double *pPartials = partials;
    
            u = 0;
            pPartials = partials;
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( k = 0; k < _patternCount; k++ ) {
                    
                    state1 = states1[k];
                    state2 = states2[k];
                    
                    w = u;
                    
                    if (state1 < _stateCount && state2 < _stateCount) {
                        
                        for ( i = 0; i < _stateCount; i++ ) {
                            
                            *pPartials++ = matrices1[w + state1] * matrices2[w + state2];
                            w += _stateCount;
                        }
                        
                    }
                    else if (state1 < _stateCount) {
                        // child 1 has a gap or unknown state so treat it as unknown
                        
                        for ( i = 0; i < _stateCount ; i++ ) {
                            
                            *pPartials++ = matrices1[w + state1];
                            
                            w += _stateCount;
                        }
                    } else if (state2 < _stateCount ) {
                        // child 2 has a gap or unknown state so treat it as unknown
                        
                        for ( i = 0; i < _stateCount ; i++ ) {
                            
                            *pPartials++ = matrices2[w + state2];
                            
                            w += _stateCount;
                        }
                    } else {
                        // both children have a gap or unknown state so set partials to 1
                        
                        for ( i = 0; i < _stateCount; i++ ) {
                            *pPartials++ = 1.0;
                        }
                    }
                }
                u += _matrixSize;
            }
        }
        
        
        void SimpleFlexibleTreeLikelihood::updatePartialsKnownUndefined( const int *states1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 )const{
            double sum;
            int v = 0;
            int i,j,k;
            int w;
            int state1;
            
            double *pPartials = partials3;
            
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( k = 0; k < _patternCount; k++ ) {
                    
                    state1 = states1[k];
                    
                    w = l * _matrixSize;
                    
                    if ( state1 < _stateCount) {
                        
                        
                        for ( i = 0; i < _stateCount; i++) {
                            
                            *pPartials = matrices1[w + state1];
                            
                            sum = 0.0;
                            for ( j = 0; j < _stateCount; j++) {
                                sum += matrices2[w] * partials2[v + j];
                                w++;
                            }
                            
                            *pPartials *= sum;
                            pPartials++;
                        }
                        
                    }
                    else {
                        // Child 1 has a gap or unknown state so don't use it
                        
                        for ( i = 0; i < _stateCount; i++) {
                            
                            *pPartials = 0.0;
                            for ( j = 0; j < _stateCount; j++) {
                                *pPartials += matrices2[w] * partials2[v + j];
                                w++;
                            }
                            
                            pPartials++;
                        }
                        
                    }
                    
                    v += _stateCount;
                }
            }
        }        
        
        void SimpleFlexibleTreeLikelihood::updatePartialsUndefinedUndefined( const double *partials1, const double *matrices1, const double *partials2, const double *matrices2, double *partials3 )const{
            double sum1, sum2;
            int v = 0;
            int i,j,k;
            int w;
            double *pPartials = partials3;
            
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( k = 0; k < _patternCount; k++ ) {
                    
                    w = l * _matrixSize;
                    
                    for ( i = 0; i < _stateCount; i++ ) {
                        
                        sum1 = sum2 = 0.;
                        
                        for ( j = 0; j < _stateCount; j++) {
                            sum1 += matrices1[w] * partials1[v + j];
                            sum2 += matrices2[w] * partials2[v + j];
                            
                            w++;
                        }
                        
                        *pPartials++ = sum1 * sum2;
                    }
                    v += _stateCount;
                }
            }

        }
        
        void SimpleFlexibleTreeLikelihood::updatePartials(int partialsIndex, int partialsIndex1, int matrixIndex1, int partialsIndex2, int matrixIndex2 ) {
            if( _partials[partialsIndex1].size() > 0 ){
                if(  _partials[partialsIndex2].size() > 0 ){
                    updatePartialsUndefinedUndefined(_partials[partialsIndex1].data(),
                                                     _matrices[matrixIndex1].data(),
                                                     _partials[partialsIndex2].data(),
                                                     _matrices[matrixIndex2].data(),
                                                     _partials[partialsIndex].data());
                }
                else {
                    updatePartialsKnownUndefined(_states[partialsIndex2].data(),
                                                 _matrices[matrixIndex2].data(),
                                                 _partials[partialsIndex1].data(),
                                                 _matrices[matrixIndex1].data(),
                                                 _partials[partialsIndex].data());
                }
                
            }
            else{
                if(  _partials[partialsIndex2].size() > 0 ){
                    updatePartialsKnownUndefined(_states[partialsIndex1].data(),
                                                 _matrices[matrixIndex1].data(),
                                                 _partials[partialsIndex2].data(),
                                                 _matrices[matrixIndex2].data(),
                                                 _partials[partialsIndex].data());
                    
                }
                else{
                    updatePartialsKnownKnown(_states[partialsIndex1].data(),
                                             _matrices[matrixIndex1].data(),
                                             _states[partialsIndex2].data(),
                                             _matrices[matrixIndex2].data(),
                                             _partials[partialsIndex].data());
                }
            }
            
            if ( _useScaleFactors ) {
                //SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
            }
            operationCallCount++;
        }
        
        void SimpleFlexibleTreeLikelihood::calculateBranchLikelihood(double* rootPartials, const double* attachmentPartials, const double* pendantPartials, const double* pendantMatrices, const double* weights){
            memset(rootPartials, 0, sizeof(double)*_stateCount*_patternCount);
            int v = 0;
            for(int l = 0; l < _rateCount; l++) {
                int u = 0;
                const double weight = weights[l];
                for(int k = 0; k < _patternCount; k++) {
                    int w = l * _matrixSize;
                    const double* partialsChildPtr = pendantPartials+v;
                    for(int i = 0; i < _stateCount; i++) {
                        const double* transMatrixPtr = &pendantMatrices[w];
                        double sum = 0.0;
                        for (int j = 0; j < _stateCount; j++) {
                            sum += transMatrixPtr[j] * partialsChildPtr[j];
                        }
                        rootPartials[u] += sum * attachmentPartials[v] * weight;
                        u++;
                        v++;
                        w += _stateCount;
                    }
                }
            }
                
            
        }
        
        void SimpleFlexibleTreeLikelihood::calculateBranchLikelihood(double* rootPartials, const double* attachmentPartials, const int* pendantStates, const double* pendantMatrices, const double* weights){
            memset(rootPartials, 0, sizeof(double)*_stateCount*_patternCount);
            int v = 0;
            for(int l = 0; l < _rateCount; l++) {
                int u = 0; // Index in resulting product-partials (summed over categories)
                const double weight = weights[l];
                for(int k = 0; k < _patternCount; k++) {
                    int w =  l * _matrixSize + pendantStates[k];
                    for(int i = 0; i < _stateCount; i++) {
                        rootPartials[u] += pendantMatrices[w] * attachmentPartials[v] * weight;
                        u++;
                        v++;
                        w += _stateCount;
                    }
                }
            }
            
        }
        
        double SimpleFlexibleTreeLikelihood::calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength){
            
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
            
            // update matrices of pendant, proximal and distal
            int offset = 0;
            for(int c = 0; c < _rateCount; c++){
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixPendant  = _model->getPij_t(pendantLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixPendant.row(i);
                    std::copy(row.begin(), row.end(), _matrices[_totalNodeCount].data()+offset);
                    offset += _stateCount;
                }
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixDistal  = _model->getPij_t(distalLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixDistal.row(i);
                    std::copy(row.begin(), row.end(), _matrices[_totalNodeCount+1].data()+offset);
                    offset += _stateCount;
                }
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixProximal  = _model->getPij_t(proximalLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixProximal.row(i);
                    std::copy(row.begin(), row.end(), _matrices[_totalNodeCount+2].data()+offset);
                    offset += _stateCount;
                }
                
            }
            
            // Distal and Proximal  are attached
            updatePartials(tmpPartialsIndex, distalIndex, _totalNodeCount+1, _upperPartialsIndexes[distalIndex], _totalNodeCount+2);
            
            const vector<double>& weights = _rateDist->getProbabilities();


            if(_partials[_totalNodeCount].size()>0){
                calculateBranchLikelihood(_rootPartials[1].data(), _partials[tmpPartialsIndex].data(), _partials[indexTaxon].data(), _matrices[_totalNodeCount].data(), weights.data());
            }
            else{
                calculateBranchLikelihood(_rootPartials[1].data(), _partials[tmpPartialsIndex].data(), _states[indexTaxon].data(), _matrices[_totalNodeCount].data(), weights.data());
            }
            
            double* patternLikelihood = _patternLikelihoods[1].data();
			calculatePatternLikelihood(_rootPartials[1].data(), _model->getFrequencies().data(), patternLikelihood);
                
            double logLnl = 0.;
            for ( int i = 0; i < _patternCount; i++) {
                logLnl += log(patternLikelihood[i]) * _patternWeights[i];
            }
            
            return logLnl;
        }
        
        void SimpleFlexibleTreeLikelihood::calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
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
            
            const vector<double>& weights = _rateDist->getProbabilities();
            
            const int distalIndex = distal.getId();
            const int indexTaxon = std::find(_taxa.begin(), _taxa.end(), taxonName) - _taxa.begin();
            const int tmpPartialsIndex = _totalNodeCount*2;
            
            // update matrices of pendant, proximal and distal
            int offset = 0;
            for(int c = 0; c < _rateCount; c++){
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixDistal  = _model->getPij_t(distalLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixDistal.row(i);
                    std::copy(row.begin(), row.end(), _matrices[_totalNodeCount+1].data()+offset);
                    offset += _stateCount;
                }
                
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixProximal  = _model->getPij_t(proximalLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixProximal.row(i);
                    std::copy(row.begin(), row.end(), _matrices[_totalNodeCount+2].data()+offset);
                    offset += _stateCount;
                }
                
            }
            
            offset = 0;
            for(int c = 0; c < _rateCount; c++){
                offset = c * _matrixSize;
                const bpp::Matrix<double>& dMatrix  = _model->getdPij_dt(pendantLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = dMatrix.row(i);
                    std::copy(row.begin(), row.end(), _matrices[_totalNodeCount].data()+offset);
                    offset += _stateCount;
                }
                for (int k = 0; k < _matrixSize; k++) {
                    _matrices[_totalNodeCount][c*_matrixSize+k] *= _rateDist->getCategory(c);
                }
            }
            
            // Proximal and distal  are attached
            updatePartials(tmpPartialsIndex, distalIndex, _totalNodeCount+1, _upperPartialsIndexes[distalIndex], _totalNodeCount+2);
            
            if(_partials[_totalNodeCount].size()>0){
                calculateBranchLikelihood(_rootPartials[2].data(), _partials[tmpPartialsIndex].data(), _partials[indexTaxon].data(), _matrices[_totalNodeCount].data(), weights.data());
            }
            else{
                calculateBranchLikelihood(_rootPartials[2].data(), _partials[tmpPartialsIndex].data(), _states[indexTaxon].data(), _matrices[_totalNodeCount].data(), weights.data());
            }
            
            
            double* patternLikelihood = _patternLikelihoods[1].data();
            double* d1PatternLikelihood = _patternLikelihoods[2].data();
            
            calculatePatternLikelihood(_rootPartials[2].data(), _model->getFrequencies().data(), d1PatternLikelihood);
            
            if(d1 != NULL){
                double dd1 = 0.;
                for ( int i = 0; i < _patternCount; i++) {
                    dd1 += (d1PatternLikelihood[i]/patternLikelihood[i]) * _patternWeights[i];
                }
                *d1 = dd1;
            }
            
            if(d2 != NULL){
                offset = 0;
                for(int c = 0; c < _rateCount; c++){
                    const double rate = _rateDist->getCategory(c);
                    offset = c * _matrixSize;
                    const bpp::Matrix<double>& d2Matrix = _model->getd2Pij_dt2(pendantLength*rate);
                    for(int i = 0; i < _stateCount; i++){
                        const vector<double>& row = d2Matrix.row(i);
                        std::copy(row.begin(), row.end(), _matrices[_totalNodeCount].data()+offset);
                        offset += _stateCount;
                    }
                    for (int k = 0; k < _matrixSize; k++) {
                        _matrices[_totalNodeCount][c*_matrixSize+k] *= rate*rate;
                    }
                }
                
                // Proximal and distal  are attached
                updatePartials(tmpPartialsIndex, distalIndex, _totalNodeCount+1, _upperPartialsIndexes[distalIndex], _totalNodeCount+2);
                
                
                if(_partials[_totalNodeCount].size()>0){
                    calculateBranchLikelihood(_rootPartials[3].data(), _partials[tmpPartialsIndex].data(), _partials[indexTaxon].data(), _matrices[_totalNodeCount].data(), weights.data());
                }
                else{
                    calculateBranchLikelihood(_rootPartials[3].data(), _partials[tmpPartialsIndex].data(), _states[indexTaxon].data(), _matrices[_totalNodeCount].data(), weights.data());
                }
                
                
                double* d2PatternLikelihood = _patternLikelihoods[3].data();
                calculatePatternLikelihood(_rootPartials[3].data(), _model->getFrequencies().data(), d2PatternLikelihood);
                
                double dd2 = 0.;
                for ( int i = 0; i < _patternCount; i++) {
                    dd2 += ((d2PatternLikelihood[i]*patternLikelihood[i]  - d1PatternLikelihood[i]*d1PatternLikelihood[i])/(patternLikelihood[i]*patternLikelihood[i])) * _patternWeights[i];
                }
                *d2 = dd2;
            }
        }
        
        void SimpleFlexibleTreeLikelihood::calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
            
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
            
            const vector<double>& weights = _rateDist->getProbabilities();
            
            const int distalIndex = distal.getId();
            const int indexTaxon = std::find(_taxa.begin(), _taxa.end(), taxonName) - _taxa.begin();
            const int tmpPartialsIndex = _totalNodeCount*2;
            
            // Temporary transition matrices
            const int tempMatrixPendant = _totalNodeCount;
            const int tempMatrixDistal = _totalNodeCount + 1;
            const int tempMatrixProximal =  _totalNodeCount + 2;
            
            // update matrices of pendant, proximal and distal
            int offset = 0;
            for(int c = 0; c < _rateCount; c++){
                const double rate = _rateDist->getCategory(c);
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixPendant  = _model->getPij_t(pendantLength*rate);
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixPendant.row(i);
                    std::copy(row.begin(), row.end(), _matrices[tempMatrixPendant].data()+offset);
                    offset += _stateCount;
                }

                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixProximal  = _model->getPij_t(proximalLength*rate);
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixProximal.row(i);
                    std::copy(row.begin(), row.end(), _matrices[tempMatrixProximal].data()+offset);
                    offset += _stateCount;
                }
            }

            offset = 0;
            for(int c = 0; c < _rateCount; c++){
                const double rate = _rateDist->getCategory(c);
                offset = c * _matrixSize;
                const bpp::Matrix<double>& dMatrix  = _model->getdPij_dt(distalLength*rate);
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = dMatrix.row(i);
                    std::copy(row.begin(), row.end(), _matrices[tempMatrixDistal].data()+offset);
                    offset += _stateCount;
                }
                for (int k = 0; k < _matrixSize; k++) {
                    _matrices[tempMatrixDistal][c*_matrixSize+k] *= rate;
                }
            }
            
            // Pendant and proximal  are attached
            updatePartials(tmpPartialsIndex, indexTaxon, tempMatrixPendant, _upperPartialsIndexes[distalIndex], tempMatrixProximal);
            
            if(_partials[distalIndex].size()>0){
                calculateBranchLikelihood(_rootPartials[2].data(), _partials[tmpPartialsIndex].data(), _partials[distalIndex].data(), _matrices[tempMatrixDistal].data(), weights.data());
            }
            else{
                calculateBranchLikelihood(_rootPartials[2].data(), _partials[tmpPartialsIndex].data(), _states[distalIndex].data(), _matrices[tempMatrixDistal].data(), weights.data());
            }
            
            
            double* patternLikelihood = _patternLikelihoods[1].data();
            double* d1PatternLikelihood = _patternLikelihoods[2].data();
            
            calculatePatternLikelihood(_rootPartials[2].data(), _model->getFrequencies().data(), d1PatternLikelihood);
            
            if(d1 != NULL){
                double dd1 = 0.;
                for ( int i = 0; i < _patternCount; i++) {
                    dd1 += (d1PatternLikelihood[i]/patternLikelihood[i]) * _patternWeights[i];
                }
                *d1 = dd1;
            }
            
            if(d2 != NULL){
                offset = 0;
                for(int c = 0; c < _rateCount; c++){
                    const double rate = _rateDist->getCategory(c);
                    offset = c * _matrixSize;
                    const bpp::Matrix<double>& d2MatrixPendant  = _model->getd2Pij_dt2(distalLength*rate);
                    for(int i = 0; i < _stateCount; i++){
                        const vector<double>& row = d2MatrixPendant.row(i);
                        std::copy(row.begin(), row.end(), _matrices[tempMatrixDistal].data()+offset);
                        offset += _stateCount;
                    }
                    for (int k = 0; k < _matrixSize; k++) {
                        _matrices[tempMatrixDistal][c*_matrixSize+k] *= rate*rate;
                    }
                }
                
                // Pendant and proximal  are attached
                updatePartials(tmpPartialsIndex, indexTaxon, tempMatrixPendant, _upperPartialsIndexes[distalIndex], tempMatrixProximal);
                
                
                if(_partials[distalIndex].size()>0){
	                calculateBranchLikelihood(_rootPartials[3].data(), _partials[tmpPartialsIndex].data(), _partials[distalIndex].data(), _matrices[tempMatrixDistal].data(), weights.data());
				}
				else{
					calculateBranchLikelihood(_rootPartials[3].data(), _partials[tmpPartialsIndex].data(), _states[distalIndex].data(), _matrices[tempMatrixDistal].data(), weights.data());
				}
                

                double* d2PatternLikelihood = _patternLikelihoods[3].data();
                calculatePatternLikelihood(_rootPartials[3].data(), _model->getFrequencies().data(), d2PatternLikelihood);
            
                double dd2 = 0.;
                for ( int i = 0; i < _patternCount; i++) {
                    dd2 += ((d2PatternLikelihood[i]*patternLikelihood[i]  - d1PatternLikelihood[i]*d1PatternLikelihood[i])/(patternLikelihood[i]*patternLikelihood[i])) * _patternWeights[i];
                }
                *d2 = dd2;
            }
            
        }
        
        double SimpleFlexibleTreeLikelihood::calculateLogLikelihood(){
            
            if(!_updatePartials)return _logLnl;
            
            traverse(_tree->getRootNode());
            
            _logLnl = 0;
            
            for(int pass = 0; pass < 2; pass++){
                
                if (_useScaleFactors) {
//                    if (_recomputeScaleFactors) {
//                        cumulateScaleBufferIndex = _internalNodeCount;
//                        beagle_check(beagleResetScaleFactors(_beagleInstance, cumulateScaleBufferIndex));
//                        
//                        beagle_check(beagleAccumulateScaleFactors(_beagleInstance,
//                                                                  _scaleBufferIndices.data(),
//                                                                  _internalNodeCount,
//                                                                  cumulateScaleBufferIndex));
//                    }
//                    else {
//                        cumulateScaleBufferIndex = _internalNodeCount;
//                    }
                }
                integratePartials(_partials[_tree->getRootNode()->getId()].data(), _rateDist->getProbabilities().data(), _rootPartials[0].data());
                
                calculatePatternLikelihood(_rootPartials[0].data(), _model->getFrequencies().data(), _patternLikelihoods[0].data());
                
                const double* patternLikelihood = _patternLikelihoods[0].data();
                _logLnl = 0;
                for ( int i = 0; i < _patternCount; i++) {
                    _logLnl += log(patternLikelihood[i]) * _patternWeights[i];
                }

                if (std::isnan(_logLnl) || std::isinf(_logLnl)) {
                    _useScaleFactors = true;
                    _recomputeScaleFactors = true;
                    
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
        
        void SimpleFlexibleTreeLikelihood::updateNode(const bpp::Node& node){
            _needNodeUpdate[node.getId()] = true;
            _updatePartials = true;
            _updateUpperPartials = true;
        }
        
        void SimpleFlexibleTreeLikelihood::updateAllNodes(){
            _needNodeUpdate.assign(_needNodeUpdate.size(), true);
            _updatePartials = true;
            _updateUpperPartials = true;
        }
     
    }
}
