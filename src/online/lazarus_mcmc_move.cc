#include "lazarus_mcmc_move.h"
#include "composite_tree_likelihood.h"
#include "multiplier_proposal.h"
#include "online_util.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_set>
#include <tuple>
#include <unordered_map>
#include <string>
#include <Bpp/Phyl/TreeTemplateTools.h>

using namespace bpp;

namespace sts { namespace online {

LazarusMCMCMove::LazarusMCMCMove(CompositeTreeLikelihood& calculator,
                                       smc::sampler<TreeParticle>& sampler,
                                        std::unordered_map<std::string, size_t>& leaf_ids) :
    calculator(calculator),
    sampler(sampler),
    leaf_ids(leaf_ids)
{}

LazarusMCMCMove::~LazarusMCMCMove()
{
    // Debug bits
    if(n_attempted > 0) {
        std::clog << "Lazarus_mcmc_move: " << n_accepted << '/' << n_attempted << ": " << acceptanceProbability() << std::endl;
    }
}

// compute the splits in the tree
const std::vector< std::pair< std::vector<bool>, bpp::Node* > > LazarusMCMCMove::getSplits( bpp::TreeTemplate<bpp::Node>* tree, bool remove_last )
{
    std::vector< std::pair< std::vector<bool>, bpp::Node* > > splits;

    std::vector<bpp::Node*> nodes = postorder(tree->getRootNode());
    nodes.pop_back();// ignore root split for now, TODO
    //std::unordered_map<bpp::Node*, double> hashkeys;
    std::unordered_map<bpp::Node* , std::vector<bool>> split_map;
    bpp::Node* unwantedNode = tree->getRootNode()->getSon(1);
    const size_t numberOfLeaves = tree->getNumberOfLeaves();
    
    for(auto n : nodes){
        // In sts every tree is rooted using the same leaf as outgroup.
        // The root child with index 1 is always ignored when adding sequences or in other mcmc moves. Its branch length is set to 0 so no need to exchange this one.
        // If we use the other child it is like duplicating the historic tree
        std::vector<bool> split(numberOfLeaves, false);
        if(n->isLeaf()){
            split[leaf_ids[n->getName()]]=true;
//            hashkeys[n] = n->getDistanceToFather();
        }else{
            std::transform(split_map[n->getSon(0)].begin(), split_map[n->getSon(0)].end(), split_map[n->getSon(1)].begin(), split.begin(), std::logical_or<bool>());
//            hashkeys[n] = n->getDistanceToFather() + hashkeys[n->getSon(0)] + hashkeys[n->getSon(1)];
        }
        // always start with false (1100 = 0011)
        if(split[0] == true){
            std::for_each(split.begin(), split.end(), std::logical_not<bool>());
        }
        split_map[n] = split;
        
        // avoid redundant bipartition
        if(n != unwantedNode){
            splits.push_back(std::make_pair(split_map[n], n));
        }
    }

    if(remove_last){ // remove the last leaf from the split system
        for(size_t n = 0; n<splits.size(); n++){
            splits[n].first.pop_back();
        }
    }

    return splits;
}

    void deleteClade(bpp::Node* node){
        for(bpp::Node* n : node->getSons()){
            deleteClade(n);
        }
        delete node;
    }
    
int LazarusMCMCMove::proposeMove(long lTime, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    if(lTime<2) return 0; // no history yet
    bpp::TreeTemplate<bpp::Node>* tree = particle.GetValue().tree.get();
    
    
    // Find the new leaf
    const size_t treeLeafCount = tree->getNumberOfLeaves();
    std::string new_id = "";
    for (auto key = leaf_ids.cbegin(); key != leaf_ids.cend(); ++key) {
        if(key->second == treeLeafCount-1){
            new_id = key->first;
            break;
        }
    }
    bpp::Node* new_leaf = tree->getNode(new_id);
    
    const std::vector< std::pair< std::vector<bool>, bpp::Node* > >& splits = getSplits(tree, true);
    // Remove splits that have the new leaf below the edge defining the split
    std::vector<int> ancs = tree->getAncestorsId(new_leaf->getId());
    
    
    // get the history
    const smc::history< smc::particle<TreeParticle> >* h = sampler.GetHistory();
    
    // get the particles from one generation back
    const smc::historyelement<smc::particle<TreeParticle> >* ele = h->GetElement();
    while(ele->GetNext() != NULL){
        ele = ele->GetNext();
    }
    const smc::particle<TreeParticle>* particles = ele->GetValues();
    const long particleCount = ele->GetNumber();
    
    // make a list of all bipartitions in the previous generation
    std::unordered_map< std::vector<bool>, std::vector< bpp::Node* > > all_splits;
    std::unordered_map< size_t, bool > done;
    for( int p=0; p < particleCount; p++ ){
        bpp::TreeTemplate<bpp::Node>* htree = particles[p].GetValue().tree.get();
        //std::cout << htree->getNumberOfLeaves() << std::endl;
        if(done.find(particles[p].GetValue().particleID) == done.end()){
            const std::vector< std::pair< std::vector<bool>, bpp::Node*> >& splits = getSplits( htree, false );
            for( auto s : splits ){
                all_splits[s.first].push_back( s.second );
            }
            done[particles[p].GetValue().particleID] = true;
        }
    }
    
    done.clear();
    // map all node pointers to their trees
    std::unordered_map< bpp::Node*, bpp::TreeTemplate<bpp::Node>* > node_tree_map;
    for( int p=0; p < particleCount; p++ ){
        if(done.find(particles[p].GetValue().particleID) == done.end()){
            bpp::TreeTemplate<bpp::Node>* old_tree = particles[p].GetValue().tree.get();
            for(auto n : old_tree->getNodes()){
                node_tree_map[n]=old_tree;
            }
            done[particles[p].GetValue().particleID] = true;
        }
    }
    
    // find common splits
    std::vector< std::tuple< std::vector<bool>, bpp::Node*, bpp::Node*> > common_splits;
    for( auto s : splits ){
        // if the leaf not present in the old trees is below the node defining the split then
        // the split is not used since we would lose the leaf in the exchange
        if ( std::find(ancs.begin(), ancs.end(), s.second->getId()) == ancs.end() ){
            for( auto a : all_splits[s.first] ){
                common_splits.push_back( std::make_tuple( s.first, s.second, a ) );
            }
        }
    }
    
    assert(common_splits.size() != 0);

    // calculate likelihood before subtree resurrection
    calculator.initialize(*particle.GetValuePointer()->model,
                          *particle.GetValuePointer()->rateDist,
                          *tree);
    double orig_ll = calculator();


    // choose one of the common subtrees to resurrect
    size_t idx = rng->UniformDiscrete(0, common_splits.size() - 1);

    // create a new tree that is a hybrid of the current tree and the historical tree
    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> new_tree(tree->clone());
    
    bpp::Node* n_cur = get<1>(common_splits[idx]);
    bpp::Node* n_old = get<2>(common_splits[idx]);
    
    bpp::Node* n_cur_new = new_tree->getNode(static_cast<int>(n_cur->getId()));
    
    
    // Uses only one branch
    if (n_cur->getFather() == NULL) {
        n_cur_new->setDistanceToFather(n_old->getDistanceToFather());
    }
    else{
        bpp::TreeTemplate<bpp::Node>* donor_tree = node_tree_map[n_old];
        bpp::TreeTemplate<bpp::Node>* new_subtree = donor_tree->cloneSubtree(n_old->getId());
        
        bpp::Node* f_cur_new = n_cur_new->getFather();
        size_t new_son_p = f_cur_new->getSonPosition(n_cur_new);
        f_cur_new->setSon(new_son_p, new_subtree->getRootNode());
        f_cur_new->getSon(new_son_p)->setDistanceToFather(n_old->getDistanceToFather());
        deleteClade(n_cur_new);
    }
    
    // Leaf IDs reflect their ranking in the alignment (starting from 0)
    // Internal nodes ID are greater or equal than the total number of sequences
    size_t nameCounter = leaf_ids.size();
    for(bpp::Node* node : new_tree->getNodes()){
        if(!node->isLeaf()){
            node->setId(static_cast<int>(nameCounter));
            nameCounter++;
        }
    }

    // calculate new tree likelihood
    calculator.initialize(*particle.GetValuePointer()->model,
                          *particle.GetValuePointer()->rateDist,
                          *new_tree);
    double new_ll = calculator();

    // the hastings ratio is the probability of picking the given subtree
    // divided by the probability of picking the original subtree from the previous generation
    double hastingsRatio=1;
    double mh_ratio = std::exp(new_ll + std::log(hastingsRatio) - orig_ll);
    if(mh_ratio >= 1.0 || rng->UniformS() < mh_ratio) {
        particle.GetValuePointer()->tree = std::move(new_tree);
        return 1;
    } else {
        // Rejected
        return 0;
    }

}

}}
