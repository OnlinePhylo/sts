#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <vector>
#include <stack>
#include <unordered_set>
#include <algorithm>
#include <memory>
#include <assert.h>
#include "smctc.hh"
#include "phylofunc.hh"
#include "hmsbeagle.hh"

using namespace std;

vector< shared_ptr< phylo_node > > leaf_nodes;
std::vector< std::pair< std::string, std::string > > aln;
std::unordered_map< std::shared_ptr< phylo_node >, int > leaf_sequence_ids;
OnlineCalculator calc;


phylo_node::phylo_node() : id(-1) {}
phylo_node::~phylo_node()
{
    if(id >= 0) calc.free_id(id);
}

bool phylo_node::is_leaf()
{
    return this->child1 == NULL && this->child2 == NULL;
}

/// The function corresponding to the log likelihood at specified time and position (up to normalisation)
/// \param lTime The current time (i.e. the number of coalescence events so far)
/// \param X     The state to consider
/// \return the log likelihood.
double logLikelihood(long lTime, const particle& X)
{

    // Walk backwards through the forest to calculate likelihoods of each tree.
    vector<bool> visited;
    double ll_sum = 0;
    shared_ptr< phylo_particle > cur = X.pp;
    for(shared_ptr< phylo_particle > cur = X.pp; cur != NULL; cur = cur->predecessor) {
        if(cur->node == NULL) continue; // this particle's join was made obsolete
        if(visited.size() < cur->node->id || !visited[ cur->node->id ]) {
            ll_sum += calc.calculate_ll(cur->node, visited);
        }
    }
    // add background freqs for all uncoalesced leaves
    for(int i = 0; i < leaf_nodes.size(); i++) {
        if(visited.size() > i && visited[i]) continue;
        ll_sum += calc.calculate_ll(leaf_nodes[i], visited);
    }

    return ll_sum;
}

/// A function to give an initial particle.
/// Because our starting particle is \perp, we don't use the rng.
/// \param pRng A pointer to the random number generator which is to be used.
/// \return An initial particle.
smc::particle<particle> fInitialise(smc::rng *pRng)
{
    particle value;
    value.pp = make_shared< phylo_particle >();
    // The loglike should just be the background distribution on character state frequencies.
    return smc::particle<particle>(value, logLikelihood(0, value));
}

/// Find the number of trees (that is, trees consisting of more than one node) from a collection of uncoalesced nodes.
/// \param uncoalesced The uncoalesced nodes.
/// \return The count.
int tree_count(const vector< shared_ptr< phylo_node > > &uncoalesced)
{
    int result = 0;
    for(auto i = uncoalesced.begin(), j = uncoalesced.end(); i != j; ++i) {
        if(!i->get()->is_leaf())
            result++;
    }
    return result;
}

/// Collect all of the subtree phylo_nodes into a vector.
/// \param node Desired root node.
/// \return vector of nodes below.
vector< shared_ptr< phylo_node > > subtree_nodes(const shared_ptr<phylo_node> node)
{
    vector< shared_ptr< phylo_node > > subtree_node_v;
    stack< shared_ptr< phylo_node > > s;
    s.push(node);
    while(s.size() > 0) {
        shared_ptr< phylo_node > n = s.top();
        s.pop();
        if(n->child1 == NULL) continue;	// leaf node, nothing more to do.
        subtree_node_v.push_back(n->child1);
        subtree_node_v.push_back(n->child2);
        s.push(n->child1);
        s.push(n->child2);
    }
    return subtree_node_v;
}

/// Find the uncoalesced nodes for a particle.
/// \param pp Input particle
/// \return vector of uncoalesced phylo_nodes.
vector< shared_ptr< phylo_node > > uncoalesced_nodes(const shared_ptr<phylo_particle> pp)
{
    // Our set of phylo nodes that can be used in proposal.
    unordered_set< shared_ptr< phylo_node > > proposal_set;
    // The nodes that have already been coalesced, to be removed later.
    unordered_set< shared_ptr< phylo_node > > coalesced;
    // Insert all of the leaf nodes into the proposal set.
    proposal_set.insert(leaf_nodes.begin(), leaf_nodes.end());
    // Walk back to predecessor particles, adding root nodes to
    // proposal_set and collecting coalesced nodes in `coalesced`.
    for(shared_ptr< phylo_particle > cur = pp->predecessor; cur != NULL; cur = cur->predecessor) {
        // Skip if the particle is \perp.
        if(cur->node == NULL) continue;
        // Skip if we've already processed this subtree, such that it's already found in coalesced.
        if(coalesced.find(cur->node) != coalesced.end()) continue;
        // Insert this active root node to the proposal set.
        proposal_set.insert(cur->node);
        vector< shared_ptr< phylo_node > > subtree = subtree_nodes(cur->node);
        coalesced.insert(subtree.begin(), subtree.end());
    }

    vector< shared_ptr< phylo_node > > pvec(proposal_set.begin(), proposal_set.end());
    vector< shared_ptr< phylo_node > > cvec(coalesced.begin(), coalesced.end());
    sort(pvec.begin(), pvec.end());
    sort(cvec.begin(), cvec.end());

    // The set difference of available (i.e. proposal_set) and coalesced nodes yields the final proposal set; store it
    // in prop_vector.
    vector< shared_ptr< phylo_node > > prop_vector(proposal_set.size() + coalesced.size());
    // UGH: std::set_difference requires an ordered container class
    // AG: that's the only way to do a set difference efficiently, right?
    vector< shared_ptr< phylo_node > >::iterator last_ins = set_difference(pvec.begin(), pvec.end(), cvec.begin(),
            cvec.end(), prop_vector.begin());
    prop_vector.resize(last_ins - prop_vector.begin());

    return prop_vector;
}

/// Get all phylo_nodes for a particle.
/// \param pp Input particle.
/// \return vector of all phylo_nodes.
vector< shared_ptr< phylo_node > > all_nodes(const shared_ptr<phylo_particle> pp)
{
    unordered_set< shared_ptr< phylo_node > > nodes;
    for(shared_ptr< phylo_particle > cur = pp->predecessor; cur != NULL; cur = cur->predecessor) {
        if(cur->node == NULL) continue;
        // Skip if we've already processed this subtree.
        if(nodes.find(cur->node) != nodes.end()) continue;
        // Recursively add all descendants of the root nodes to the node set using a stack.
        vector< shared_ptr< phylo_node > > subtree = subtree_nodes(cur->node);
        nodes.insert(subtree.begin(), subtree.end());
        nodes.insert(cur->node);
    }
    // add in any leaf nodes we may have missed
    nodes.insert(leaf_nodes.begin(), leaf_nodes.end());
    vector< shared_ptr< phylo_node > > nvec( nodes.begin(), nodes.end() );
    return nvec;
}

/// Determine the path from a root to target.
/// \param target Desired destination.
/// \param uncoalesced The uncoalesced phylo_nodes.
/// \return vector of phylo_nodes going from a root to the target.
vector< shared_ptr< phylo_node > > find_root_path(const shared_ptr<phylo_node> target, const vector< shared_ptr< phylo_node > >& uncoalesced)
{
    vector< shared_ptr< phylo_node > > path;
    for(auto n : uncoalesced){
        bool found = false;
        stack< shared_ptr< phylo_node > > s;
        unordered_set< shared_ptr< phylo_node > > visited;
        s.push(n);
        while(s.size()>0){
            shared_ptr< phylo_node > cur = s.top();
            s.pop();
            if(found && (path.back() == cur->child1 || path.back() == cur->child2)){
                // We have found target, and are unwinding to get all entries behind us for path.
                path.push_back(cur);
            }
            else if(cur == target){
                // Success!
                path.push_back(cur);
                found = true;
            }
            else if(visited.find(cur) != visited.end()){
                // If cur has already been visited do nothing.
            }
            else if(cur->child1 != NULL && !found){
                // There is more to explore: we haven't found target, and we haven't visited this node.
                assert(cur->child2 != NULL);
                s.push(n); // Keep current node around so we can visit again post-order.
                s.push(n->child1);
                s.push(n->child2);
            }
            else{
                // We haven't found target or visited this leaf node. Just continue.
            }
            visited.insert(cur);
        }
        if(found) break;
    }
    assert(!path.empty());
    return path;
}

/// Pick an uncoalesced node uniformly at random, then pick an edge to join it to, also uniformly at random.
/// Find the path to the root from the join location and create new nodes along that path.
/// Any new nodes will get likelihood peeled.
/// \param lTime The current time.
/// \param pFrom The current particle. Gets modified in place.
/// \param pRng Our rng.
void move_root_join_anywhere(long lTime, smc::particle<particle>& pFrom, smc::rng *pRng)
{
    // Set pp up as the next phylo_particle and make pFrom point at it.
    particle* part = pFrom.GetValuePointer();
    shared_ptr< phylo_particle > pp = make_shared< phylo_particle >();
    pp->predecessor = part->pp;
    part->pp = pp;
    // Pick a root node uniformly.
    vector< shared_ptr< phylo_node > > uncoalesced = uncoalesced_nodes(pp);
    int n1 = pRng->UniformDiscrete(0, uncoalesced.size() - 1);
    // Uniformly pick a node, excluding that uniformly picked root node.
    vector< shared_ptr< phylo_node > > subtree = subtree_nodes(uncoalesced[n1]);
    subtree.push_back(uncoalesced[n1]);
    sort( subtree.begin(), subtree.end() );
    vector< shared_ptr< phylo_node > > all = all_nodes(pp);
    sort( all.begin(), all.end() );
    vector< shared_ptr< phylo_node > > dest_set( all.size() - subtree.size() );
    set_difference(all.begin(), all.end(), subtree.begin(), subtree.end(), dest_set.begin() );
    int n2 = pRng->UniformDiscrete(0, dest_set.size() - 1);
    // We will propose to attach above n2.
    // Determine the path from a root to n2.
    vector< shared_ptr< phylo_node > > path = find_root_path( dest_set[n2], uncoalesced );
    // Walk through the particle history and eliminate references to the deleted nodes. These no longer exist because
    // the merges they represent no longer happened... they now have one more taxon in them.
    int i=path.size()-1;
    for(shared_ptr< phylo_particle > cur = pp->predecessor; cur != NULL && i > 0; cur = cur->predecessor) {
        if(cur->node == path[i]){
            cur->node = NULL;
            // AD: why wouldn't we skip over this one via a predecessor re-assignment rather than makeing cur-> node = NULL?
            i--;
        }
    }
    // Create a new node.
    shared_ptr< phylo_node > join_node = make_shared< phylo_node >();
    join_node->id = calc.get_id();
    join_node->child1 = uncoalesced[n1];
    join_node->child2 = dest_set[n2];
    // Propose distances to child nodes from Exponential(1)
    join_node->dist1 = pRng->Exponential(1.0);
    join_node->dist2 = pRng->Exponential(1.0);
    double h_prob1 = exp(-join_node->dist1);
    double h_prob2 = exp(-join_node->dist2);
    shared_ptr< phylo_node > cur_node = join_node;
    for(int i=1; i<path.size(); i++){
        shared_ptr< phylo_node > new_node = make_shared< phylo_node >();
        new_node->id = calc.get_id();
        new_node->child1 = cur_node;
        new_node->child2 = path[i]->child2 == path[i-1] ? path[i]->child1 : path[i]->child2;
        new_node->dist1 = path[i]->dist1;
        new_node->dist1 = path[i]->child2 == path[i-1] ? path[i]->dist2 : path[i]->dist1;
        new_node->dist2 = path[i]->child2 == path[i-1] ? path[i]->dist1 : path[i]->dist2;
        cur_node = new_node;
    }
    // cur_node is the new root for this particle.
    pp->node = cur_node;

    //XXX: need to add in reverse proposal probability. This should be the number of edges in the forest.
    pFrom.AddToLogWeight(logLikelihood(lTime, *part) - log(h_prob1) - log(h_prob1));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<particle>& pFrom, smc::rng *pRng)
{
    particle* part = pFrom.GetValuePointer();
    shared_ptr< phylo_particle > pp = make_shared< phylo_particle >();
    pp->predecessor = part->pp;
    part->pp = pp;
    vector< shared_ptr< phylo_node > > prop_vector = uncoalesced_nodes(pp);

    // Pick two nodes from the prop_vector to join.
    int n1 = pRng->UniformDiscrete(0, prop_vector.size() - 1);
    int n2 = pRng->UniformDiscrete(0, prop_vector.size() - 2);
    if(n2 >= n1) n2++;
    pp->node = make_shared< phylo_node >();
    pp->node->id = calc.get_id();
    pp->node->child1 = prop_vector[n1];
    pp->node->child2 = prop_vector[n2];
    // Propose child node distances from a normal distribution centered on the average
    // distance in partial probability vectors
    if(pp->node->child1->is_leaf() && pp->node->child2->is_leaf()) {
        int c1id = leaf_sequence_ids[pp->node->child1];
        int c2id = leaf_sequence_ids[pp->node->child2];
        double dist = 0;
        for(int i = 0; i < aln[c1id].second.size(); i++) {
            if(aln[c1id].second[i] != '-' && aln[c2id].second[i] != '-' && aln[c1id].second[i] != aln[c2id].second[i]) dist++;
        }
        dist /= aln[c1id].second.size();
        double gg = -3.0 / 4.0 * log(((1.0 - dist) - 1 / 4.0) / (3 / 4.0));
        double ss = gg / 10.0;
        pp->node->dist1 = pRng->Normal(gg / 2.0, ss * ss);
        pp->node->dist2 = pRng->Normal(gg / 2.0, ss * ss);
        const double PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982;
        // Calculate normal probability density
        double d_prob1 = exp(-pow(pp->node->dist1 - gg, 2) / (2 * ss * ss)) / (ss * sqrt(2.0 * PI));
        double d_prob2 = exp(-pow(pp->node->dist2 - gg, 2) / (2 * ss * ss)) / (ss * sqrt(2.0 * PI));
        pFrom.AddToLogWeight(logLikelihood(lTime, *part) - log(d_prob1) - log(d_prob1));
    } else {
        // Propose a coalescence time.
        pp->node->dist1 = pRng->Exponential(1.0);
        pp->node->dist2 = pRng->Exponential(1.0);
        double h_prob1 = exp(-pp->node->dist1);
        double h_prob2 = exp(-pp->node->dist2);
        pFrom.AddToLogWeight(logLikelihood(lTime, *part) - log(h_prob1) - log(h_prob1));
    }

    // Add reverse transition probability q(r' -> r)
    // 1/(# of trees in forest), omitting trees consisting of a single leaf
    // This prevents over-counting
    const int tc = tree_count(prop_vector);
    if(tc > 1)
        pFrom.AddToLogWeight(-log(tc));
}

int moves = 0;
int moves_accepted = 0;

int fMoveNodeAgeMCMC(long lTime, smc::particle<particle>& pFrom, smc::rng *pRng)
{
    moves++;
    particle* part = pFrom.GetValuePointer();
    shared_ptr< phylo_node > cur_node = part->pp->node;
    shared_ptr< phylo_node > new_node = make_shared< phylo_node >();
    new_node->child1 = cur_node->child1;
    new_node->child2 = cur_node->child2;
    new_node->id = calc.get_id();

    double cur_ll = logLikelihood(lTime, *part);
    // Choose an amount to shift each child distance uniformly at random.
    double shift1 = pRng->Uniform(-0.1, 0.1);
    double shift2 = pRng->Uniform(-0.1, 0.1);
    // If the shift amount would create a negative child distance we will reflect it back to a positive number.
    // This means the proposal density for the reflection zone is double the rest of the area, but the back-proposal
    // probability is also double in the same area so these terms cancel in the Metropolis-Hastings ratio.
    new_node->dist1 = abs(cur_node->dist1 + shift1);
    new_node->dist2 = abs(cur_node->dist2 + shift2);
    part->pp->node = new_node;

    double alpha = exp(logLikelihood(lTime, *part) - cur_ll);
    if(alpha < 1 && pRng->UniformS() > alpha) {
        // Move rejected, restore the original node.
        part->pp->node = cur_node;
        return false;
    }
    // Accept the new state.
    moves_accepted++;
    return true;
}
