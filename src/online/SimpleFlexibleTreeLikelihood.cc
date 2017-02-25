//
//  SimpleFlexibleTreeLikelihood.cpp
//  sts
//
//  Created by Mathieu Fourment on 23/06/2016.
//  Copyright © 2016 Mathieu Fourment. All rights reserved.
//

#include "SimpleFlexibleTreeLikelihood.hpp"


using namespace std;

// g++ -lhmsbeagle -I/usr/local/include/libhmsbeagle-1/ -std=c++11  -lbpp-seq -lbpp-core -lbpp-phyl -g SimpleFlexibleTreeLikelihood.cpp AbstractFlexibleTreeLikelihood.cpp

namespace sts {
    namespace online {
        
        SimpleFlexibleTreeLikelihood::SimpleFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, const bpp::SubstitutionModel &model, const bpp::DiscreteDistribution& rateDist):
        AbstractFlexibleTreeLikelihood(patterns, model, rateDist){

            _matrixSize = _stateCount*_stateCount;
            
            _patternWeights.resize(_patternCount);
            const std::vector<unsigned int>& w = _patterns.getWeights();
            auto castit = [](unsigned int w) { return static_cast<double>(w); };
            std::transform(w.begin(), w.end(), _patternWeights.begin(), castit);
            _states.resize(_totalNodeCount);
            for(auto it = _states.begin(); it != _states.end(); ++it){
                it->resize(_patternCount);
            }
            std::unique_ptr<bpp::SiteContainer> sites(patterns.getSites());
            setStates(*sites.get());
            
            // Lower partials
            // Leaves don't need partials
            _partials.resize(_totalNodeCount);
            for(auto it = _partials.begin(); it != _partials.end(); ++it){
                it->resize(_rateCount*_stateCount*_patternCount);
            }
            
            // Probability matrices
            _matrices.resize(_totalNodeCount);
            for(auto it = _matrices.begin(); it != _matrices.end(); ++it){
                it->resize(_rateCount*_matrixSize);
            }
            
            // Rates are integrated out
            _rootPartials.assign(_stateCount*_patternCount, 0.);
            _patternLikelihoods.assign(_patternCount, 0.);
            
            // Upper partials
            _upperPartials.resize(_totalNodeCount);
            for(auto it = _upperPartials.begin(); it != _upperPartials.end(); ++it){
                it->resize(_rateCount*_stateCount*_patternCount);
            }
            
            // Temporary storage
            _temporaryMatrices.resize(3);
            for (int i = 0; i < _temporaryMatrices.size(); i++) {
                _temporaryMatrices[i].assign(_rateCount*_matrixSize, 0.);
            }
            
            _temporaryPartials.assign(_rateCount*_stateCount*_patternCount, 0.);
            
            // 1 for LnL with upper likelihood
            // 2 for 1st and 2nd derivatives
            _temporaryPatternLikelihood.resize(3);
            for (int i = 0; i < _temporaryPatternLikelihood.size(); i++) {
                _temporaryPatternLikelihood[i].assign(_patternCount, 0.);
            }
            
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
                    
                    updatePartials(node->getSon(0)->getId(), node->getSon(1)->getId(), node->getId());

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
                outLogLikelihoods[k] = log(outLogLikelihoods[k]);
                
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
        
        void SimpleFlexibleTreeLikelihood::updatePartials( int nodeIndex1, int nodeIndex2, const double* matrix1, const double *matrix2, double *partials ) {
            if( nodeIndex1 >= _sequenceCount ){
                if(  nodeIndex2 >= _sequenceCount ){
                    updatePartialsUndefinedUndefined(_partials[nodeIndex1].data(),
                                                     matrix1,
                                                     _partials[nodeIndex2].data(),
                                                     matrix2,
                                                     partials);
                }
                else {
                    updatePartialsKnownUndefined(_states[nodeIndex2].data(),
                                                 matrix2,
                                                 _partials[nodeIndex1].data(),
                                                 matrix1,
                                                 partials);
                }
                
            }
            else{
                if(  nodeIndex2 >= _sequenceCount ){
                    updatePartialsKnownUndefined(_states[nodeIndex1].data(),
                                                 matrix1,
                                                 _partials[nodeIndex2].data(),
                                                 matrix2,
                                                 partials);
                    
                }
                else{
                    updatePartialsKnownKnown(_states[nodeIndex1].data(),
                                             matrix1,
                                             _states[nodeIndex2].data(),
                                             matrix2,
                                             partials);
                }
            }
            
            if ( _useScaleFactors ) {
                //SingleTreeLikelihood_scalePartials( tlk, nodeIndex3);
            }
        }
        
        void SimpleFlexibleTreeLikelihood::updatePartials( int nodeIndex1, int nodeIndex2, int nodeIndex3 ) {
            updatePartials(nodeIndex1, nodeIndex2, _matrices[nodeIndex1].data(), _matrices[nodeIndex2].data(), _partials[nodeIndex3].data());
        }
        
        // Called by a node whose parent is NOT the root and the node's sibling is a leaf
        void SimpleFlexibleTreeLikelihood::updateUpperPartialsKnown( const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const int *states, double *partials ) const{
            int w,w2,i,j,k;
            int v = 0;
            double sum1,sum2;
            int state;
            double *pPartials = partials;
            
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( k = 0; k < _patternCount; k++ ) {
                    
                    state = states[k];
                    
                    w  = l * _matrixSize;
                    
                    for ( i = 0; i < _stateCount; i++ ) {
                        w2 = l * _matrixSize + i;
                        sum1 = sum2 = 0.;
                        
                        for ( j = 0; j < _stateCount; j++) {
                            sum1 += matrix_upper[w2] * partials_upper[v + j];
                            w2 += _stateCount;
                        }
                        
                        if( state < _stateCount){
                            sum2 = matrix_lower[w+state];
                            w += _stateCount;
                        }
                        else {
                            sum2 = 1.0;
                            w += _stateCount;
                        }
                        
                        *pPartials++ = sum1 * sum2 ;
                    }
                    v += _stateCount;
                }
            }
        }
        
        // Called by a node whose parent is NOT the root and the node's sibling is NOT a leaf
        void SimpleFlexibleTreeLikelihood::updateUpperPartialsUndefined( const double *matrix_upper, const double *partials_upper, const double *matrix_lower, const double *partials_lower, double *partials ) const{
            int w,w2,i,j,k;
            int v = 0;
            double sum1,sum2;
            double *pPartials = partials;
            
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( k = 0; k < _patternCount; k++ ) {
                    
                    w  = l * _matrixSize;
                    
                    for ( i = 0; i < _stateCount; i++ ) {
                        w2 = l * _matrixSize + i;
                        sum1 = sum2 = 0.;
                        
                        for ( j = 0; j < _stateCount; j++) {
                            sum1 += matrix_upper[w2] * partials_upper[v + j];
                            sum2 += matrix_lower[w]  * partials_lower[v + j];
                            w++;
                            w2 += _stateCount;
                        }
                        
                        *pPartials++ = sum1 * sum2 ;
                    }
                    v += _stateCount;
                }
            }
        }
        
        // Called by a child of the root but not if the child's sibling is NOT leaf
        void SimpleFlexibleTreeLikelihood::updateUpperPartialsRootUndefined( const double *partials1, const double *matrices1, const double *frequencies, double *partials_upper )const{
            int w,i,j,k;
            int v = 0;
            double *pPartials = partials_upper;
            
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( k = 0; k < _patternCount; k++ ) {
                    
                    w = l * _matrixSize;
                    
                    for ( i = 0; i < _stateCount; i++ ) {
                        
                        *pPartials = 0.;
                        
                        for ( j = 0; j < _stateCount; j++ ) {
                            *pPartials += matrices1[w] * partials1[v + j];
                            w++;
                        }
                        *pPartials *= frequencies[i];
                        pPartials++;
                    }
                    v += _stateCount;
                }
            }
        }
        
        // Called by a child of the root and the child's sibling is a leaf
        void SimpleFlexibleTreeLikelihood::updateUpperPartialsRootKnown(const double *matrices1, const int *states1, const double *frequencies, double *partials_upper ) const{
            int w;
            double *pPartials = partials_upper;
            int state1;
            
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( int k = 0; k < _patternCount; k++ ) {
                    
                    state1 = states1[k];
                    
                    w = l * _matrixSize;
                    
                    if( state1 < _stateCount ){
                        for ( int i = 0; i < _stateCount; i++ ) {
                            *pPartials++ = matrices1[w+state1] * frequencies[i];
                            w += _stateCount;
                        }
                    }
                    else {
                        memcpy(pPartials, frequencies, sizeof(double)*_stateCount);
                        pPartials += _stateCount;
                    }
                    
                }
            }
        }
        
        void SimpleFlexibleTreeLikelihood::updateUpperPartials(const bpp::Node* node){
            if( node->hasFather() ){
                const bpp::Node *parent = node->getFather();
                const bpp::Node *sibling = (parent->getSon(0) == node ? parent->getSon(1) : parent->getSon(0));
                
                if( !parent->hasFather() ){
                    if( sibling->isLeaf() ){
                        updateUpperPartialsRootKnown(_matrices[sibling->getId()].data(), _states[sibling->getId()].data(), _model->getFrequencies().data(), _upperPartials[node->getId()].data() );
                    }
                    else {
                        updateUpperPartialsRootUndefined(_partials[sibling->getId()].data(),  _matrices[sibling->getId()].data(), _model->getFrequencies().data(), _upperPartials[node->getId()].data() );
                    }
                }
                else if( sibling->isLeaf() ){
                    updateUpperPartialsKnown(_matrices[parent->getId()].data(), _upperPartials[parent->getId()].data(), _matrices[sibling->getId()].data(), _states[sibling->getId()].data(), _upperPartials[node->getId()].data());
                }
                else {
                    updateUpperPartialsUndefined(_matrices[parent->getId()].data(), _upperPartials[parent->getId()].data(), _matrices[sibling->getId()].data(), _partials[sibling->getId()].data(), _upperPartials[node->getId()].data());
                }
            }
            
            if( !node->isLeaf() ){
                updateUpperPartials(node->getSon(0));
                updateUpperPartials(node->getSon(1));
            }
        }
        
        void SimpleFlexibleTreeLikelihood::calculatePatternLikelihood( const double *partials_upper, const double *partials_lower, const double *matrix_lower, const double *proportions, double *pattern_lk ) const{
            int w,i,j,k;
            int v = 0;
            double p,sum;
            
            memset(pattern_lk, 0, sizeof(double)*_patternCount);
            
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( k = 0; k < _patternCount; k++ ) {
                    
                    w = l * _matrixSize;
                    p = 0;
                    
                    for ( i = 0; i < _stateCount; i++ ) {
                        sum = 0;
                        
                        for ( j = 0; j < _stateCount; j++) {
                            sum += matrix_lower[w] * partials_lower[v + j];
                            w++;
                        }
                        
                        p += sum * partials_upper[v + i];
                    }
                    pattern_lk[k] += p * proportions[l];
                    v += _stateCount;
                }
            }
        }
        
        void SimpleFlexibleTreeLikelihood::calculatePatternLikelihoodLeaf( const double *partials_upper, const int* states, const double *matrix_lower, const double *proportions, double *pattern_lk ) const{
            int w,j,k;
            int v = 0;
            double p = 0;
            int state;
            
            memset(pattern_lk, 0, sizeof(double)*_patternCount);
            
            for ( int l = 0; l < _rateCount; l++ ) {
                
                for ( k = 0; k < _patternCount; k++ ) {
                    state = states[k];
                    
                    w = l * _matrixSize;
                    p = 0;
                    if( state < _stateCount ){
                        
                        for ( j = 0; j < _stateCount; j++) {
                            p += matrix_lower[w+state] * partials_upper[v + j];
                            w += _stateCount;
                        }
                    }
                    else {
                        for ( j = 0; j < _stateCount; j++) {
                            p += partials_upper[v + j];
                        }
                    }
                    pattern_lk[k] += p * proportions[l];
                    v += _stateCount;
                }
            }
        }
        
        double SimpleFlexibleTreeLikelihood::calculateLogLikelihood(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength){
            
            if(_updatePartials){
                traverse(_tree->getRootNode());
                updateUpperPartials(_tree->getRootNode());
                _updatePartials = _updateUpperPartials = false;
                _needNodeUpdate.assign(_totalNodeCount, false);
            }
            else if(_updateUpperPartials){
                updateUpperPartials(_tree->getRootNode());
                _updateUpperPartials = false;
            }
            
            int distalIndex = distal.getId();
            int indexTaxon = std::find(_taxa.begin(), _taxa.end(), taxonName) - _taxa.begin();
            int newNodeIndex = _tree->getNumberOfNodes()+2;
            
            // update matrices of pendant, proximal and distal
            int offset = 0;
            for(int c = 0; c < _rateCount; c++){
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixPendant  = _model->getPij_t(pendantLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixPendant.row(i);
                    std::copy(row.begin(), row.end(), _temporaryMatrices[0].data()+offset);
                    offset += _stateCount;
                }
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixDistal  = _model->getPij_t(distalLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixDistal.row(i);
                    std::copy(row.begin(), row.end(), _temporaryMatrices[1].data()+offset);
                    offset += _stateCount;
                }
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixProximal  = _model->getPij_t(proximalLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixProximal.row(i);
                    std::copy(row.begin(), row.end(), _temporaryMatrices[2].data()+offset);
                    offset += _stateCount;
                }
                
            }
            
            // update lower partials at new node
            updatePartials(indexTaxon, distalIndex, _temporaryMatrices[0].data(), _temporaryMatrices[1].data(), _temporaryPartials.data());
            const vector<double>& rates = _rateDist->getProbabilities();
            
            calculatePatternLikelihood(_upperPartials[distal.getId()].data(), _temporaryPartials.data(), _temporaryMatrices[2].data(), rates.data(), _temporaryPatternLikelihood[0].data());
            
            double logLnl = 0.;
            for ( int i = 0; i < _patternCount; i++) {
                logLnl += log(_temporaryPatternLikelihood[0][i]) * _patternWeights[i];
            }
            
            return logLnl;
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
                integratePartials(_partials[_tree->getRootNode()->getId()].data(), _rateDist->getProbabilities().data(), _rootPartials.data());
                
                calculatePatternLikelihood(_rootPartials.data(), _model->getFrequencies().data(), _patternLikelihoods.data());
                
                _logLnl = 0;
                for ( int i = 0; i < _patternCount; i++) {
                    _logLnl += _patternLikelihoods[i] * _patternWeights[i];
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
        
        void SimpleFlexibleTreeLikelihood::calculateDistalDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
            calculateDerivatives(distal, taxonName, distal.getId(), pendantLength, distalLength, proximalLength, d1, d2);
        }
        
        void SimpleFlexibleTreeLikelihood::calculatePendantDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
            int indexTaxon = std::find(_taxa.begin(), _taxa.end(), taxonName) - _taxa.begin();
            calculateDerivatives(distal, taxonName, indexTaxon, pendantLength, distalLength, proximalLength, d1, d2);
//            assert(d1 != NULL || d2 != NULL);
//            
//            int distalIndex = distal.getId();
//            int indexTaxon = _patterns.getSites()->getSequencePosition(taxonName);
//            int newNodeIndex = _tree->getNumberOfNodes();
//            
//            // update matrices of proximal and distal
//            int offset = 0;
//            for(int c = 0; c < _rateCount; c++){
//                offset = c * _matrixSize;
//                const bpp::Matrix<double>& matrixDistal  = _model->getPij_t(distalLength*_rateDist->getCategory(c));
//                for(int i = 0; i < _stateCount; i++){
//                    const vector<double>& row = matrixDistal.row(i);
//                    std::copy(row.begin(), row.end(), _temporaryMatrices[1].data()+offset);
//                    offset += _stateCount;
//                }
//                offset = c * _matrixSize;
//                const bpp::Matrix<double>& matrixProximal  = _model->getPij_t(proximalLength*_rateDist->getCategory(c));
//                for(int i = 0; i < _stateCount; i++){
//                    const vector<double>& row = matrixProximal.row(i);
//                    std::copy(row.begin(), row.end(), _temporaryMatrices[2].data()+offset);
//                    offset += _stateCount;
//                }
//            }
//            
//            const vector<double>& rates = _rateDist->getProbabilities();
//            
//            // Only d1
//            if(d2 == NULL){
//                int offset = 0;
//                for(int c = 0; c < _rateCount; c++){
//                    offset = c * _matrixSize;
//                    const bpp::Matrix<double>& dMatrixPendant  = _model->getdPij_dt(pendantLength*_rateDist->getCategory(c));
//                    for(int i = 0; i < _stateCount; i++){
//                        const vector<double>& row = dMatrixPendant.row(i);
//                        std::copy(row.begin(), row.end(), _temporaryMatrices[0].data()+offset);
//                        offset += _stateCount;
//                    }
//                    for (int k = 0; k < _matrixSize; k++) {
//                        _temporaryMatrices[0][c*_matrixSize+k] *= _rateDist->getCategory(c);
//                    }
//                }
//                
//                updatePartials(indexTaxon, distalIndex, _temporaryMatrices[0].data(), _temporaryMatrices[1].data(), _temporaryPartials.data());
//                
//                calculatePatternLikelihood(_upperPartials[distal.getId()].data(), _temporaryPartials.data(), _temporaryMatrices[2].data(), rates.data(), _temporaryPatternLikelihood[1].data());
//                double dd1 = 0.;
//                for ( int i = 0; i < _patternCount; i++) {
//                    dd1 += (_temporaryPatternLikelihood[1][i]/_temporaryPatternLikelihood[0][i]) * _patternWeights[i];
//                }
//                *d1 = dd1;
//            }
//            // d1 and d2 or only d2
//            else {
//                int offset = 0;
//                for(int c = 0; c < _rateCount; c++){
//                    offset = c * _matrixSize;
//                    const bpp::Matrix<double>& d1MatrixPendant  = _model->getdPij_dt(pendantLength*_rateDist->getCategory(c));
//                    for(int i = 0; i < _stateCount; i++){
//                        const vector<double>& row = d1MatrixPendant.row(i);
//                        std::copy(row.begin(), row.end(), _temporaryMatrices[0].data()+offset);
//                        offset += _stateCount;
//                    }
//                    for (int k = 0; k < _matrixSize; k++) {
//                        _temporaryMatrices[0][c*_matrixSize+k] *= _rateDist->getCategory(c);
//                    }
//                }
//                updatePartials(indexTaxon, distalIndex, _temporaryMatrices[0].data(), _temporaryMatrices[1].data(), _temporaryPartials.data());
//                calculatePatternLikelihood(_upperPartials[distal.getId()].data(), _temporaryPartials.data(), _temporaryMatrices[2].data(), rates.data(), _temporaryPatternLikelihood[1].data());
//                
//                offset = 0;
//                for(int c = 0; c < _rateCount; c++){
//                    offset = c * _matrixSize;
//                    const bpp::Matrix<double>& d2MatrixPendant  = _model->getd2Pij_dt2(pendantLength*_rateDist->getCategory(c));
//                    for(int i = 0; i < _stateCount; i++){
//                        const vector<double>& row = d2MatrixPendant.row(i);
//                        std::copy(row.begin(), row.end(), _temporaryMatrices[0].data()+offset);
//                        offset += _stateCount;
//                    }
//                    for (int k = 0; k < _matrixSize; k++) {
//                        _temporaryMatrices[0][c*_matrixSize+k] *= _rateDist->getCategory(c)*_rateDist->getCategory(c);
//                    }
//                }
//                
//                updatePartials(indexTaxon, distalIndex, _temporaryMatrices[0].data(), _temporaryMatrices[1].data(), _temporaryPartials.data());
//                calculatePatternLikelihood(_upperPartials[distal.getId()].data(), _temporaryPartials.data(), _temporaryMatrices[2].data(), rates.data(), _temporaryPatternLikelihood[2].data());
//                
//                double dd2 = 0.;
//                for ( int i = 0; i < _patternCount; i++) {
//                    dd2 += ((_temporaryPatternLikelihood[2][i]*_temporaryPatternLikelihood[0][i]  - _temporaryPatternLikelihood[1][i]*_temporaryPatternLikelihood[1][i])/(_temporaryPatternLikelihood[0][i]*_temporaryPatternLikelihood[0][i])) * _patternWeights[i];
//                }
//                *d2 = dd2;
//                
//                if(d1 != NULL){
//                    double dd1 = 0.;
//                    for ( int i = 0; i < _patternCount; i++) {
//                        dd1 += (_temporaryPatternLikelihood[1][i]/_temporaryPatternLikelihood[0][i]) * _patternWeights[i];
//                    }
//                    *d1 = dd1;
//                }
//            }
        }
        
        void SimpleFlexibleTreeLikelihood::calculateDerivatives(const bpp::Node& distal, std::string taxonName, int index, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
            assert(d1 != NULL || d2 != NULL);
            
            int distalIndex = distal.getId();
            int indexTaxon = std::find(_taxa.begin(), _taxa.end(), taxonName) - _taxa.begin();
            int newNodeIndex = _tree->getNumberOfNodes()+2;
            
            double siblingLength = pendantLength;
            double theLength = distalLength;
            int siblingNode = indexTaxon;
            int theNode = distal.getId();
            
            // derivatives of pendant
            if (index == indexTaxon) {
                siblingLength = distalLength;
                theLength = pendantLength;
                siblingNode = distal.getId();
                theNode = indexTaxon;
            }
            
            // update matrices of proximal and distal
            int offset = 0;
            for(int c = 0; c < _rateCount; c++){
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrix  = _model->getPij_t(siblingLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrix.row(i);
                    std::copy(row.begin(), row.end(), _temporaryMatrices[1].data()+offset);
                    offset += _stateCount;
                }
                offset = c * _matrixSize;
                const bpp::Matrix<double>& matrixProximal  = _model->getPij_t(proximalLength*_rateDist->getCategory(c));
                for(int i = 0; i < _stateCount; i++){
                    const vector<double>& row = matrixProximal.row(i);
                    std::copy(row.begin(), row.end(), _temporaryMatrices[2].data()+offset);
                    offset += _stateCount;
                }
            }
            
            const vector<double>& rates = _rateDist->getProbabilities();
            
            // Only d1
            if(d2 == NULL){
                int offset = 0;
                for(int c = 0; c < _rateCount; c++){
                    offset = c * _matrixSize;
                    const bpp::Matrix<double>& dMatrix  = _model->getdPij_dt(theLength*_rateDist->getCategory(c));
                    for(int i = 0; i < _stateCount; i++){
                        const vector<double>& row = dMatrix.row(i);
                        std::copy(row.begin(), row.end(), _temporaryMatrices[0].data()+offset);
                        offset += _stateCount;
                    }
                    for (int k = 0; k < _matrixSize; k++) {
                        _temporaryMatrices[0][c*_matrixSize+k] *= _rateDist->getCategory(c);
                    }
                }
                
                updatePartials(theNode, siblingNode, _temporaryMatrices[0].data(), _temporaryMatrices[1].data(), _temporaryPartials.data());
                
                calculatePatternLikelihood(_upperPartials[distal.getId()].data(), _temporaryPartials.data(), _temporaryMatrices[2].data(), rates.data(), _temporaryPatternLikelihood[1].data());
                double dd1 = 0.;
                for ( int i = 0; i < _patternCount; i++) {
                    dd1 += (_temporaryPatternLikelihood[1][i]/_temporaryPatternLikelihood[0][i]) * _patternWeights[i];
                }
                *d1 = dd1;
            }
            // d1 and d2 or only d2
            else {
                int offset = 0;
                for(int c = 0; c < _rateCount; c++){
                    offset = c * _matrixSize;
                    const bpp::Matrix<double>& d1Matrix  = _model->getdPij_dt(theLength*_rateDist->getCategory(c));
                    for(int i = 0; i < _stateCount; i++){
                        const vector<double>& row = d1Matrix.row(i);
                        std::copy(row.begin(), row.end(), _temporaryMatrices[0].data()+offset);
                        offset += _stateCount;
                    }
                    for (int k = 0; k < _matrixSize; k++) {
                        _temporaryMatrices[0][c*_matrixSize+k] *= _rateDist->getCategory(c);
                    }
                }
                updatePartials(theNode, siblingNode, _temporaryMatrices[0].data(), _temporaryMatrices[1].data(), _temporaryPartials.data());
                calculatePatternLikelihood(_upperPartials[distal.getId()].data(), _temporaryPartials.data(), _temporaryMatrices[2].data(), rates.data(), _temporaryPatternLikelihood[1].data());
                
                offset = 0;
                for(int c = 0; c < _rateCount; c++){
                    offset = c * _matrixSize;
                    const bpp::Matrix<double>& d2MatrixPendant  = _model->getd2Pij_dt2(theLength*_rateDist->getCategory(c));
                    for(int i = 0; i < _stateCount; i++){
                        const vector<double>& row = d2MatrixPendant.row(i);
                        std::copy(row.begin(), row.end(), _temporaryMatrices[0].data()+offset);
                        offset += _stateCount;
                    }
                    for (int k = 0; k < _matrixSize; k++) {
                        _temporaryMatrices[0][c*_matrixSize+k] *= _rateDist->getCategory(c)*_rateDist->getCategory(c);
                    }
                }
                
                updatePartials(theNode, siblingNode, _temporaryMatrices[0].data(), _temporaryMatrices[1].data(), _temporaryPartials.data());
                calculatePatternLikelihood(_upperPartials[distal.getId()].data(), _temporaryPartials.data(), _temporaryMatrices[2].data(), rates.data(), _temporaryPatternLikelihood[2].data());
                
                double dd2 = 0.;
                for ( int i = 0; i < _patternCount; i++) {
                    dd2 += ((_temporaryPatternLikelihood[2][i]*_temporaryPatternLikelihood[0][i]  - _temporaryPatternLikelihood[1][i]*_temporaryPatternLikelihood[1][i])/(_temporaryPatternLikelihood[0][i]*_temporaryPatternLikelihood[0][i])) * _patternWeights[i];
                }
                *d2 = dd2;
                
                if(d1 != NULL){
                    double dd1 = 0.;
                    for ( int i = 0; i < _patternCount; i++) {
                        dd1 += (_temporaryPatternLikelihood[1][i]/_temporaryPatternLikelihood[0][i]) * _patternWeights[i];
                    }
                    *d1 = dd1;
                }
            }
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


//#include <Bpp/Phyl/PatternTools.h>
//#include <Bpp/Seq/Container/SequenceContainer.h>
//#include <Bpp/Seq/Container/SiteContainerTools.h>
//#include <Bpp/Seq/Container/VectorSiteContainer.h>
//#include <Bpp/Seq/Io/IoSequenceFactory.h>
//#include <Bpp/Seq/Io/ISequence.h>
//#include <Bpp/Seq/SiteTools.h>
////#include <iostream>
//#include <Bpp/Seq/Alphabet/DNA.h>
//#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
//#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
//
//const bpp::DNA DNA;
//
//int main(int argc, char **argv){
//    std::ifstream in("/Users/mathieu/Desktop/rep0/sequences.fa");
//    bpp::IoSequenceFactory fac;
//    std::unique_ptr<bpp::ISequence> reader = std::unique_ptr<bpp::ISequence>(fac.createReader(bpp::IoSequenceFactory::FASTA_FORMAT));
//    std::unique_ptr<bpp::SequenceContainer> raw_seqs(reader->readSequences(in, &DNA));
//    //bpp::SiteContainer* sequences = new bpp::VectorSiteContainer(*raw_seqs);
//    std::unique_ptr<bpp::SiteContainer> sequences(new bpp::VectorSiteContainer(*raw_seqs));
//    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);
//    in.close();
//    
//    std::unique_ptr<bpp::SitePatterns> _patterns(new bpp::SitePatterns(sequences.get()));
//    
//    bpp::JCnuc model(&DNA);
//    bpp::ConstantRateDistribution rateDist;
//    
//    bpp::Node *root = new bpp::Node("root");
//    
//    bpp::Node* new_leaf1 = new bpp::Node(0, sequences->getSequence(0).getName());
//    bpp::Node* new_leaf2 = new bpp::Node(1, sequences->getSequence(1).getName());
//    
//    root->addSon(new_leaf1);
//    root->addSon(new_leaf2);
//    //((taxon43_1997.0:0.53115658,taxon20_1995.0:0.01274880)node1:0.01368644,taxon1_1997.0:0.00000000);
//    double pendantLength = 0.53115658;
//    double distalLength = 0.01274880;
//    double proximalLength = 0.01368644;
//    root->getSons()[0]->setDistanceToFather(distalLength+proximalLength);
//    root->getSons()[1]->setDistanceToFather(0);
//    
//    std::unique_ptr<bpp::TreeTemplate<bpp::Node> > tree(new bpp::TreeTemplate<bpp::Node>(root));
//    
//    vector<bpp::Node*> r = tree->getNodes();
//    vector<string> names = sequences->getSequencesNames();
//    int counter = sequences->getNumberOfSequences();
//    for(bpp::Node* node : r){
//        if(node->getNumberOfSons() == 0){
//            int pos = find(names.begin(), names.end(), node->getName()) - names.begin();
//            node->setId(pos);
//        }
//        else {
//            node->setId(counter++);
//        }
//        cout<<node->getName()<< " "<<node->getId()<<endl;
//    }
//    
//    sts::online::SimpleFlexibleTreeLikelihood likelihood(*_patterns.get(), model, rateDist);
//    likelihood.initialize(*tree, model, rateDist);
//    likelihood.updateAllNodes();
//    double score = likelihood.calculateLogLikelihood();
//    cout << "LnL2: " << score << endl;
//    
//    bpp::Node* distal = root->getSon(0);
//    double lnl3a = likelihood.calculateLogLikelihood(*distal, sequences->getSequence(2).getName(), pendantLength, distalLength, proximalLength);
//    cout << "LnL3: " << lnl3a << endl;
//    
//    double d1=0;
//    double d2=0;
//    likelihood.calculateDerivatives(*distal, sequences->getSequence(2).getName(), pendantLength, distalLength, proximalLength, &d1, &d2);
//    cout <<"Pendant d1: "<<d1<<" "<<" d2: "<<d2<<endl;
//    
//    double lnl3p1 = likelihood.calculateLogLikelihood(*distal, sequences->getSequence(2).getName(), pendantLength-1e-7, distalLength, proximalLength);
//    double lnl3p2 = likelihood.calculateLogLikelihood(*distal, sequences->getSequence(2).getName(), pendantLength+1e-7, distalLength, proximalLength);
//    cout <<"Distal approximate d1: "<<(lnl3p2-lnl3p1)/2e-7<<endl;
//    
//    d2 = -1;
//    likelihood.calculateDistalDerivatives(*distal, sequences->getSequence(2).getName(), pendantLength, distalLength, proximalLength, &d1, &d2);
//    cout <<"Distal d1: "<<d1<<" "<<" d2: "<<d2<<endl;
//    
//    double lnl3d1 = likelihood.calculateLogLikelihood(*distal, sequences->getSequence(2).getName(), pendantLength, distalLength-1e-7, proximalLength);
//    double lnl3d2 = likelihood.calculateLogLikelihood(*distal, sequences->getSequence(2).getName(), pendantLength, distalLength+1e-7, proximalLength);
//    cout <<"Distal approximate d1: "<<(lnl3d2-lnl3d1)/2e-7<<endl;
//    
//    bpp::Node* new_leaf3 = new bpp::Node(2, sequences->getSequence(2).getName());
//    bpp::Node* temp = new bpp::Node(counter++, "leaf3dad");
//
//    size_t pos = root->getSonPosition(distal);
//    root->setSon(pos, temp);
//    
//    temp->addSon(distal);
//    temp->addSon(new_leaf3);
//    temp->setDistanceToFather(proximalLength);
//    new_leaf3->setDistanceToFather(pendantLength);
//    distal->setDistanceToFather(distalLength);
//    
////    likelihood.updateNode(*new_leaf3);
////    likelihood.updateNode(*temp);
////    likelihood.updateNode(*distal);
//    
//    likelihood.initialize(*tree, model, rateDist);
//    double lnl3b = likelihood.calculateLogLikelihood();
//    cout << "LnL3: " << lnl3b << " = " << lnl3a << endl;
////
////    score = p.getScore(*tree, *temp, sequences->getSequence(3).getName());
////    cout << "Score: " << score << endl;
//    
//    return 0;
//}
