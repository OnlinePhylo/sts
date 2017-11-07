#ifndef STS_ONLINE_LAZARUS_MCMC_MOVE_H
#define STS_ONLINE_LAZARUS_MCMC_MOVE_H

#include "online_mcmc_move.h"
#include "tree_particle.h"
#include "smctc.hh"

namespace sts { namespace online {

// Forwards
class CompositeTreeLikelihood;

class LazarusMCMCMove : public OnlineMCMCMove
{
public:
    LazarusMCMCMove(CompositeTreeLikelihood& calculator,
                        smc::sampler<TreeParticle>& sampler,
                        std::unordered_map<std::string, size_t>& leaf_ids);
    ~LazarusMCMCMove();
    std::pair<int, double> proposeMove(long, smc::particle<TreeParticle>&, smc::rng*);
protected:
    const std::vector< std::pair< std::vector<bool>, bpp::Node* > > getSplits( bpp::TreeTemplate<bpp::Node>* tree, bool remove_last = false );

private:
    CompositeTreeLikelihood& calculator;
    smc::sampler<TreeParticle>& sampler;
    std::unordered_map<std::string, size_t>& leaf_ids;
    long previous_lTime;
    std::unordered_map< bpp::Node*, bpp::TreeTemplate<bpp::Node>* > _node_tree_map;
    //std::unordered_map< std::vector<bool>, std::vector< bpp::Node* > > _all_splits;
    std::unordered_map<std::vector<bool>, std::unordered_map<size_t, bpp::Node*> > _all_splits2;
    std::unordered_map<size_t, int> _counts;
    int _counter;
    std::vector<int> _debug;
    std::vector<int> _same_tree;
    std::unordered_map<bpp::TreeTemplate<bpp::Node>*, const TreeParticle*> _particle_map;
};

}}

#endif
