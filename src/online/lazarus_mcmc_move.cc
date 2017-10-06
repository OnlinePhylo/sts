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
std::vector< std::pair< std::vector<bool>, std::pair<bpp::Node*,bpp::Node*> > > LazarusMCMCMove::getSplits( bpp::TreeTemplate<bpp::Node>* tree, bool remove_last ) 
{
    std::vector< std::pair< std::vector<bool>, std::pair<bpp::Node*,bpp::Node*> > > splits;

    std::vector<bpp::Node*> nodes = postorder(tree->getRootNode());
    std::unordered_map<bpp::Node*, double> hashkeys;
    std::unordered_map<bpp::Node* , std::vector<bool>> split_map;
    for(auto n : nodes){
        if(n->getFather() == NULL)  continue;   // ignore root split for now, TODO
        std::vector<bool> split;
        split.resize(tree->getNumberOfLeaves());
        if(n->isLeaf()){
            split[leaf_ids[n->getName()]]=true;
//            hashkeys[n] = n->getDistanceToFather();
        }else{
            std::transform(split_map[n->getSon(0)].begin(), split_map[n->getSon(0)].end(), split_map[n->getSon(1)].begin(), split.begin(), std::logical_or<bool>());
//            hashkeys[n] = n->getDistanceToFather() + hashkeys[n->getSon(0)] + hashkeys[n->getSon(1)];
        }
        split_map[n] = split;

        splits.push_back(std::make_pair(split, std::make_pair(n,n->getFather())));
        std::vector<bool> notsplit;
        notsplit.resize(tree->getNumberOfLeaves());
        std::transform(split.begin(), split.end(), notsplit.begin(), std::logical_not<bool>());
        splits.push_back(std::make_pair(notsplit, std::make_pair(n->getFather(),n)));
    }

    if(remove_last){ // remove the last leaf from the split system
        for(size_t n = 0; n<splits.size(); n++){
            splits[n].first.resize(tree->getNumberOfLeaves()-1);
        }
    }

    return splits;
}

int LazarusMCMCMove::proposeMove(long lTime, smc::particle<TreeParticle>& particle, smc::rng* rng)
{
    if(lTime<2) return 0; // no history yet

    bpp::TreeTemplate<bpp::Node>* tree = particle.GetValue().tree.get();

    // get the history
    const smc::history< smc::particle<TreeParticle> >* h = sampler.GetHistory();
    
    // get the particles from one generation back
    const smc::particle<TreeParticle>* particles = h->GetLastElement()->GetValues();

    // make a list of all bipartitions in the previous generation
    std::unordered_map< std::vector<bool>, std::vector< std::pair<bpp::Node*,bpp::Node*> > > all_splits;
    for( int p=0; p < h->GetLastElement()->GetNumber(); p++ ){
        
        bpp::TreeTemplate<bpp::Node>* htree = particles[p].GetValue().tree.get();
        if(htree->getNumberOfLeaves() != tree->getNumberOfLeaves()-1) continue;
        std::vector< std::pair< std::vector<bool>, std::pair< bpp::Node*, bpp::Node* > > > splits = getSplits( htree, false );
        for( auto s : splits ){
            all_splits[s.first].push_back( s.second );
        }
    }

    // map all node pointers to their trees
    std::unordered_map< bpp::Node*, bpp::TreeTemplate<bpp::Node>* > node_tree_map;
    for( int p=0; p < h->GetLastElement()->GetNumber(); p++ ){
        bpp::TreeTemplate<bpp::Node>* old_tree = particles[p].GetValue().tree.get();
        for(auto n : old_tree->getNodes()){
            node_tree_map[n]=old_tree;
        }
    }

    // find common splits
    std::vector< std::pair< std::vector<bool>, std::pair<bpp::Node*, bpp::Node*> > > splits = getSplits(tree, true);
    std::vector< std::tuple< std::vector<bool>, std::pair< bpp::Node*, bpp::Node* >, std::pair< bpp::Node*, bpp::Node* > > > common_splits;
    for( auto s : splits ){
        for( auto a : all_splits[s.first] ){
            if(s.second.first->getFather() == s.second.second){
                common_splits.push_back( std::make_tuple( s.first, s.second, a ) );
            }
        }
    }

    // calculate likelihood before subtree resurrection
    calculator.initialize(*particle.GetValuePointer()->model,
                          *particle.GetValuePointer()->rateDist,
                          *tree);
    double orig_ll = calculator();


    // choose one of the common subtrees to resurrect
    size_t idx = rng->UniformDiscrete(0, common_splits.size() - 1);

    // create a new tree that is a hybrid of the current tree and the historical tree
    bpp::TreeTemplate<bpp::Node>* new_tree = tree->clone();
    size_t nid_cur = get<1>(common_splits[idx]).first->getId();
    size_t pid_cur = get<1>(common_splits[idx]).second->getId();
    size_t nid_old = get<2>(common_splits[idx]).first->getId();
    bpp::TreeTemplate<bpp::Node>* donor_tree = node_tree_map[get<2>(common_splits[idx]).first];
    bpp::TreeTemplate<bpp::Node>* new_subtree = donor_tree->cloneSubtree(nid_old);
    // find the node in the current tree to replace
    for(auto n : new_tree->getNodes()){
        if(n->getId() != nid_cur) continue;
        if(n->getFather()->getId() == pid_cur){
            // connect the father to the new subtree
            double fdist = n->getDistanceToFather();
            bpp::Node* f = n->getFather();
            size_t son_p = f->getSonPosition(n);
            f->setSon(son_p, new_subtree->getRootNode());
            f->getSon(son_p)->setDistanceToFather(fdist); // maintain distance to father on this 
            n->removeFather();
            if(new_tree->getNumberOfLeaves() != tree->getNumberOfLeaves()){
                std::cerr << "broken\n";
            }
        }else{
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
        particle.GetValuePointer()->tree.reset(new_tree);
        return 1;
    } else {
        // Rejected
        // TODO: delete the new tree
        return 0;
    }

}

}}
