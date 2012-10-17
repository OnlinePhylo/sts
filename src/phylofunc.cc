#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <vector>
#include <stack>
#include <unordered_set>
#include <algorithm>
#include <memory>
#include <assert.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include "smctc.hh"

#include "phylofunc.hh"
#include "hmsbeagle.hh"

using namespace std;

edge::edge(std::shared_ptr<phylo_node> node, double length) : node(node), length(length) {}

phylo_node::phylo_node(std::shared_ptr<online_calculator> calc) : calc(calc), id(-1) {};
phylo_node::phylo_node(const phylo_node &other) : calc(other.calc), id(other.id) {};

phylo_node::~phylo_node()
{
    if(id >= 0) {
        auto p = calc.lock();
        if(p)
            p->free_id(id);
    }
}

bool phylo_node::is_leaf()
{
    return this->child1 == NULL && this->child2 == NULL;
}

void phylo_node::calc_height()
{
    if(is_leaf())
        this->height = 0.0;
    else
        this->height = max(child1->node->height + 2 * child1->length, child2->node->height + 2 * child2->length);
}

/// Find the number of trees (that is, trees consisting of more than one node) from a collection of uncoalesced nodes.
int tree_count(const vector< shared_ptr< phylo_node > > &uncoalesced)
{
    int result = 0;
    for(auto i = uncoalesced.begin(), j = uncoalesced.end(); i != j; ++i) {
        if(!i->get()->is_leaf())
            result++;
    }
    return result;
}

/// Find the uncoalesced nodes for a particle.
/// \param pp Input particle
vector< shared_ptr< phylo_node > > uncoalesced_nodes(const shared_ptr<phylo_particle> pp, const vector<shared_ptr<phylo_node>> leaf_nodes)
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
            if(n->is_leaf()) continue;	// leaf node, nothing more to do.
            coalesced.insert(n->child1->node);
            coalesced.insert(n->child2->node);
            s.push(n->child1->node);
            s.push(n->child2->node);
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
