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
#include "hmsbeagle.h"

using namespace std;

vector< shared_ptr< phylo_node > > leaf_nodes;
std::vector< std::pair< std::string, std::string > > aln;
std::vector< std::string > just_the_seqs_maam;
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

///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the number of coalescence events so far)
///  \param X     The state to consider
double logLikelihood(long lTime, const particle& X)
{
    if(just_the_seqs_maam.size() == 0) {
        for(int i = 0; i < aln.size(); i++)
            just_the_seqs_maam.push_back(aln[i].second);
        calc.initialize(just_the_seqs_maam);
    }

    // Walk backwards through the forest to calculate likelihoods of each tree.
    vector<bool> visited;
    double ll_sum = 0;
    shared_ptr< phylo_particle > cur = X.pp;
    while(cur != NULL && cur->node != NULL) {
        if(visited.size() < cur->node->id || !visited[ cur->node->id ]) {
            ll_sum += calc.calculate_ll(cur->node, visited);
        }
        cur = cur->predecessor;
    }
    // add background freqs for all uncoalesced leaves
    for(int i = 0; i < leaf_nodes.size(); i++) {
        if(visited.size() > i && visited[i]) continue;
        ll_sum += calc.calculate_ll(leaf_nodes[i], visited);
    }

    return ll_sum;
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<particle> fInitialise(smc::rng *pRng)
{
    particle value;
    // initial particles have all sequences uncoalesced
    value.pp = make_shared< phylo_particle >();
    // loglike should just be the background distribution on character state frequencies
    return smc::particle<particle>(value, logLikelihood(0, value));
}

/// Find the number of trees (that is, trees consisting of more than one node) from a collection of uncoalesced nodes.
int tree_count(const vector< shared_ptr< phylo_node > > &uncoalesced)
{
    int result = 0;
for(const shared_ptr<phylo_node> &i : uncoalesced) {
        if(!i->is_leaf())
            result++;
    }
    return result;
}

/// Find the uncoalesced nodes for a particle.
/// \param pp Input particle
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
        // Recursively add all descendants of the root nodes to the coalesced set using a stack.
        stack< shared_ptr< phylo_node > > s;
        s.push(cur->node);
        while(s.size() > 0) {
            shared_ptr< phylo_node > n = s.top();
            s.pop();
            if(n->child1 == NULL) continue;	// leaf node, nothing more to do.
            coalesced.insert(n->child1);
            coalesced.insert(n->child2);
            s.push(n->child1);
            s.push(n->child2);
        }
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
    if(pp->node->child1->id < aln.size() && pp->node->child2->id < aln.size()) {
        int c1id = pp->node->child1->id;
        int c2id = pp->node->child2->id;
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
