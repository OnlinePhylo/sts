//
//  BeagleFlexibleTreeLikleihood.cpp
//  sts
//
//  Created by Mathieu Fourment on 26/06/2016.
//  Copyright © 2016 Mathieu Fourment. All rights reserved.
//

#include "BeagleFlexibleTreeLikleihood.hpp"


#include "bpp_shim.h"
#include "util.h"

using sts::util::beagle_check;
using namespace std;

namespace sts {
    namespace online {
        
        BeagleFlexibleTreeLikelihood::BeagleFlexibleTreeLikelihood(const bpp::SitePatterns& patterns, bpp::SubstitutionModel& model, bpp::DiscreteDistribution& rateDist ):
            AbstractFlexibleTreeLikelihood(patterns, model, rateDist){
            
            _beagleInstance = beagleCreateInstance(
                                                   0,              // Number of tip data elements (input)
                                                   _partialCount,  // Number of partials buffers to create (input)
                                                   0,              // Number of compact state representation buffers to create (input)
                                                   _stateCount,    // Number of states in the continuous-time Markov chain (input)
                                                   _patternCount,  // Number of site patterns to be handled by the instance (input)
                                                   1,              // Number of rate matrix eigen-decomposition buffers to allocate (input)
                                                   _matrixCount,   // Number of rate matrix buffers (input)
                                                   _rateCount,     // Number of rate categories (input)
                                                   _patternCount,  // Number of scaling buffers (one for each internal node + 1 for accumulation) - 1 extra buffer for prox, distal
                                                   NULL,           // List of potential resource on which this instance is allowed (input, NULL implies no
                                                   // restriction
                                                   0,              // Length of resourceList list (input)
                                                   BEAGLE_FLAG_VECTOR_SSE | BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_SCALING_MANUAL, // Bit-flags indicating
                                                   //preferred implementation charactertistics, see BeagleFlags (input)
                                                   0,              // Bit-flags indicating required implementation characteristics, see BeagleFlags (input)
                                                   &_instanceDetails);
            
            if(_beagleInstance < 0)
                beagle_check(_beagleInstance);
            
			const std::vector<unsigned int> weights = _patterns.getWeights();
			const bpp::SiteContainer* sites = _patterns.getSites();
			registerLeaves(*sites);
			std::vector<double> w;
			w.reserve(_patternCount);
			auto castit = [](unsigned int w) { return static_cast<double>(w); };
			std::transform(weights.begin(), weights.end(), w.begin(), castit);
			beagle_check(beagleSetPatternWeights(_beagleInstance, w.data()));
            
                
            // all point to the last index
            _upperIndices.assign(_totalNodeCount, 2*_totalNodeCount-1);
                
            _useAutoScaling = false;//_instanceDetails.flags & BEAGLE_FLAG_SCALING_AUTO;
            _logLnl = 0;
        }
        
        
        BeagleFlexibleTreeLikelihood::~BeagleFlexibleTreeLikelihood(){
            if(_beagleInstance >= 0)
                beagleFinalizeInstance(_beagleInstance);
        }
        
        void BeagleFlexibleTreeLikelihood::initialize(bpp::TreeTemplate<bpp::Node>& tree, bpp::SubstitutionModel &model, bpp::DiscreteDistribution& rateDist){
            AbstractFlexibleTreeLikelihood::initialize(tree, model, rateDist);
            
            updateSiteModel();
            updateSubstitutionModel();
        }
        
        void BeagleFlexibleTreeLikelihood::registerLeaves(const bpp::SiteContainer& sites){
            std::vector<double> partials(_patternCount * _stateCount * _rateCount);
            
            for(int i = 0; i < _sequenceCount; i++){
                const bpp::Sequence& sequence = sites.getSequence(i);
                for(size_t site = 0; site < _patternCount; site++) {
                    for(size_t i = 0; i < _stateCount; i++) {
                        size_t idx = _stateCount * site + i;
                        partials[idx] = _model->getInitValue(i, sequence.getValue(site));
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
        
        void BeagleFlexibleTreeLikelihood::traverseUpper(const bpp::Node* node, int& index){
        
        	if(node->hasFather()){
                const bpp::Node* parent = node->getFather();
                const bpp::Node* sibling = parent->getSon(0);
                if (sibling == node) {
                    sibling = parent->getSon(1);
                }
        		int idSibling = sibling->getId();
        		
        		if(parent->hasFather()){
                    const bpp::Node* grandParent = parent->getFather();
					int idParent = parent->getId();
					int idGrandParent = grandParent->getId();					
                    _upperIndices[node->getId()] = index++;
                    
					_operations.push_back(BeagleOperation(
													     {_upperIndices[node->getId()],      // Destination buffer
														  BEAGLE_OP_NONE,    // (output) scaling buffer index
														  BEAGLE_OP_NONE,    // (input) scaling buffer index
														  _upperIndices[idGrandParent],     // Index of first child partials buffer
														  idParent,          // Index of first child transition matrix
														  idSibling,         // Index of second child partials buffer
														  idSibling})); 
				}
				else{
                    _upperIndices[node->getId()] = idSibling;
                    
				}
        	}
        	if(node->getNumberOfSons() > 0){
        	
				traverseUpper(node->getSon(0), index);
                traverseUpper(node->getSon(1), index);
        	}
        }
        
        bool BeagleFlexibleTreeLikelihood::traverse(const bpp::Node* node){
            bool update = _needNodeUpdate[node->getId()];
            
            if(node->hasFather() && update){
                _matrixUpdateIndices.push_back(node->getId());
                _branchLengths.push_back(node->getDistanceToFather());//rate
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
                    else if (_useAutoScaling) {
                        _scaleBufferIndices[node->getId() - _leafCount] = node->getId();
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
            
            if(_updateUpperPartials){
                _upperIndices.assign(_totalNodeCount, 2*_totalNodeCount-1 );
                
                int index = _totalNodeCount;
            	traverseUpper(_tree->getRootNode(), index);
                _updateUpperPartials = false;
            }
            
            const int category_weight_index = 0;
            const int state_frequency_index = 0;
            int cumulateScaleBufferIndex = BEAGLE_OP_NONE;
            
            int indexTaxon = _patterns.getSites()->getSequencePosition(taxonName);
            
            // Temporary transition matrices
            int tempMatrixDistal = _totalNodeCount;
            int tempMatrixProximal =  _totalNodeCount + 1;
            
            
            // Temporary partial
            int  tempPartialIndex = _upperIndices[indexTaxon];
            
            // Distal upper partial index
            int distalUpperIndex = _upperIndices[distal.getId()];
            
            
            //cout <<"md: "<< tempMatrixDistal << " mp: "<<tempMatrixProximal<< " taxon: " << indexTaxon<<" p: "<<tempPartialIndex<<" d: "<<distal.getId()<<" du: "<<distalUpperIndex<<endl;
            
            _branchLengths.push_back(pendantLength);
            _matrixUpdateIndices.push_back(indexTaxon);
            
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
			_operations.push_back(BeagleOperation({tempPartialIndex,       // Destination buffer
													BEAGLE_OP_NONE,        // (output) scaling buffer index
                                                    BEAGLE_OP_NONE,        // (input) scaling buffer index
                                                    distal.getId(),       // Index of first child partials buffer
                                                    tempMatrixDistal,      // Index of first child transition matrix
                                                    distalUpperIndex,      // Index of second child partials buffer
                                                    tempMatrixProximal})); // Index of second child transition matrix
                                                
            beagle_check(beagleUpdatePartials(_beagleInstance,
                                              _operations.data(),
                                              _operations.size(),
                                              BEAGLE_OP_NONE));

            double logLnl;
            /*beagle_check(beagleCalculateRootLogLikelihoods(_beagleInstance,
                                                               &rootIndex,
                                                               &category_weight_index,
                                                               &state_frequency_index,
                                                               &cumulateScaleBufferIndex,
                                                               1,  // op count
                                                               &logLnl));*/
            // Calculate log likelihood                         
            beagle_check(beagleCalculateEdgeLogLikelihoods(_beagleInstance,
                                      &tempPartialIndex,
                                      &indexTaxon,
                                      &indexTaxon,
                                      NULL,
                                      NULL,
                                      &category_weight_index,
                                      &state_frequency_index,
                                      &cumulateScaleBufferIndex,
                                      1,
                                      &logLnl,
                                      NULL,
                                      NULL));
            
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
                
                int cumulateScaleBufferIndex = BEAGLE_OP_NONE;
                if (_useScaleFactors) {
                    if (_recomputeScaleFactors) {
                        cumulateScaleBufferIndex = _internalNodeCount;
                        beagle_check(beagleResetScaleFactors(_beagleInstance, cumulateScaleBufferIndex));
                        
                        beagle_check(beagleAccumulateScaleFactors(_beagleInstance,
                                                                  _scaleBufferIndices.data(),
                                                                  _internalNodeCount,
                                                                  cumulateScaleBufferIndex));
                    }
                    else {
                        cumulateScaleBufferIndex = _internalNodeCount;
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
            
            return _logLnl;
        }
        
        void BeagleFlexibleTreeLikelihood::calculateDerivatives(const bpp::Node& distal, std::string taxonName, double pendantLength, double distalLength, double proximalLength, double* d1, double* d2){
            
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


void newick(bpp::Node* node){
    
    if(node->getNumberOfSons() > 0){
        cout <<"("<< node->getSon(0)->getName()<<"["<<node->getSon(0)->getId()<<"]:"<<node->getSon(0)->getDistanceToFather();
        newick(node->getSon(0));
        cout <<","<< node->getSon(1)->getName()<<"["<<node->getSon(1)->getId()<<"]:"<<node->getSon(1)->getDistanceToFather();
        newick(node->getSon(1));
        cout <<")";
    }
    if(!node->hasFather()){
         cout<<node->getName()<<"["<<node->getId()<<"]" <<endl;
    }
}


#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Io/IoSequenceFactory.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/SiteTools.h>
//#include <iostream>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

const bpp::DNA DNA;

int main(int argc, char **argv){
    std::ifstream in("/Users/mathieu/Desktop/rep0/sequences.fa");
    bpp::IoSequenceFactory fac;
    std::unique_ptr<bpp::ISequence> reader = std::unique_ptr<bpp::ISequence>(fac.createReader(bpp::IoSequenceFactory::FASTA_FORMAT));
    std::unique_ptr<bpp::SequenceContainer> raw_seqs(reader->readSequences(in, &DNA));
    //bpp::SiteContainer* sequences = new bpp::VectorSiteContainer(*raw_seqs);
    std::unique_ptr<bpp::SiteContainer> sequences(new bpp::VectorSiteContainer(*raw_seqs));
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sequences);
    in.close();
    
    cout << "Number of sequences: " << sequences->getNumberOfSequences()<<endl;
    
    std::unique_ptr<bpp::SitePatterns> _patterns(new bpp::SitePatterns(sequences.get()));
    
    bpp::JCnuc model(&DNA);
    bpp::ConstantRateDistribution rateDist;
    
    bpp::Node *root = new bpp::Node("root");
    
    bpp::Node* new_leaf1 = new bpp::Node(0, sequences->getSequence(0).getName());
    bpp::Node* new_leaf2 = new bpp::Node(1, sequences->getSequence(1).getName());
    
    root->addSon(new_leaf1);
    root->addSon(new_leaf2);
    root->getSons()[0]->setDistanceToFather(0.3);
    root->getSons()[1]->setDistanceToFather(0);
    
    std::unique_ptr<bpp::TreeTemplate<bpp::Node> > tree(new bpp::TreeTemplate<bpp::Node>(root));
    
    vector<bpp::Node*> r = tree->getNodes();
    vector<string> names = sequences->getSequencesNames();
    int counter = sequences->getNumberOfSequences();
    for(bpp::Node* node : r){
        if(node->getNumberOfSons() == 0){
            int pos = find(names.begin(), names.end(), node->getName()) - names.begin();
            node->setId(pos);
        }
        else {
            node->setId(counter++);
        }
        cout<<node->getName()<< " "<<node->getId()<<endl;
    }
    newick(tree->getRootNode());
    cout<<endl;
    
    sts::online::BeagleFlexibleTreeLikelihood likelihood(*_patterns.get(), model, rateDist);
    likelihood.initialize(*tree, model, rateDist);
    double lnl = likelihood.calculateLogLikelihood();
    cout << "LnL: " << lnl << endl;
    
    bpp::Node* distal = root->getSon(0);
    double pendantLength = 0.1;
    double distalLength = 0.1;
    double proximalLength = 0.2;
    double lnl2 = likelihood.calculateLogLikelihood(*distal, sequences->getSequence(2).getName(), pendantLength, distalLength, proximalLength);
    cout << "LnL2: " << lnl2 << endl;
    
    //    score = p.getScore(*tree, *root->getSon(0), sequences->getSequence(2).getName());
    //    cout << "Score: " << score << endl;
    //
    bpp::Node* new_leaf3 = new bpp::Node(2, sequences->getSequence(2).getName());
    bpp::Node* temp = new bpp::Node(counter++, "leaf3dad");
    
    size_t pos = root->getSonPosition(distal);
    root->setSon(pos, temp);
    
    temp->addSon(distal);
    temp->addSon(new_leaf3);
    temp->setDistanceToFather(proximalLength);
    distal->setDistanceToFather(distalLength);
    new_leaf3->setDistanceToFather(0.1);
    
    
    newick(tree->getRootNode());
    cout<<endl;
    
    likelihood.initialize(*tree, model, rateDist);
    double lnl3 = likelihood.calculateLogLikelihood();
    cout << "Lnl3: " << lnl3 << endl;
    //
    //    score = p.getScore(*tree, *temp, sequences->getSequence(3).getName());
    //    cout << "Score: " << score << endl;
    
    return 0;
}
