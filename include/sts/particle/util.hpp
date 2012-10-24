#ifndef STS_PARTICLE_UTIL_HPP
#define STS_PARTICLE_UTIL_HPP

#include <iostream>
#include <vector>
#include <memory>
#include <stack>
#include <string>
#include <unordered_set>

#include "sts/likelihood/online_calculator.hpp"
#include "sts/particle/phylo_node.hpp"

namespace sts
{
namespace particle
{

int tree_count(const std::vector<node> &);
std::vector<node> uncoalesced_nodes(particle pp, std::vector<node> leaf_nodes);

void write_tree(std::ostream &out, const node root, const std::unordered_map<node, std::string> &names);

/// Find the number of trees (that is, trees consisting of more than one node) from a collection of uncoalesced nodes.
/// \param uncoalesced The uncoalesced nodes.
/// \return The count.
int tree_count(const std::vector<node> &uncoalesced)
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
/// \return vector of uncoalesced phylo_nodes.
std::vector<node> uncoalesced_nodes(const particle pp, const std::vector<node> leaf_nodes)
{
    // Our set of phylo nodes that can be used in proposal.
    std::unordered_set<node> proposal_set;
    // The nodes that have already been coalesced, to be removed later.
    std::unordered_set<node> coalesced;
    // Insert all of the leaf nodes into the proposal set.
    proposal_set.insert(leaf_nodes.begin(), leaf_nodes.end());
    // Walk back to predecessor particles, adding root nodes to
    // proposal_set and collecting coalesced nodes in `coalesced`.
    for(particle cur = pp->predecessor; cur != NULL; cur = cur->predecessor) {
        // Skip if the particle is \perp.
        if(cur->node == NULL) continue;
        // Skip if we've already processed this subtree, such that it's already found in coalesced.
        if(coalesced.find(cur->node) != coalesced.end()) continue;
        // Insert this active root node to the proposal set.
        proposal_set.insert(cur->node);
        // Recursively add all descendants of the root nodes to the coalesced set using a std::stack.
        std::stack<node> s;
        s.push(cur->node);
        while(s.size() > 0) {
            node n = s.top();
            s.pop();
            if(n->is_leaf()) continue;	// leaf node, nothing more to do.
            coalesced.insert(n->child1->node);
            coalesced.insert(n->child2->node);
            s.push(n->child1->node);
            s.push(n->child2->node);
        }
    }

    std::vector<node> pvec(proposal_set.begin(), proposal_set.end());
    std::vector<node> cvec(coalesced.begin(), coalesced.end());
    sort(pvec.begin(), pvec.end());
    sort(cvec.begin(), cvec.end());

    // The set difference of available (i.e. proposal_set) and coalesced nodes yields the final proposal set; store it
    // in prop_vector.
    std::vector<node> prop_vector(proposal_set.size() + coalesced.size());
    // UGH: std::set_difference requires an ordered container class
    // AG: that's the only way to do a set difference efficiently, right?
    auto last_ins = set_difference(pvec.begin(), pvec.end(), cvec.begin(),
                                   cvec.end(), prop_vector.begin());
    prop_vector.resize(last_ins - prop_vector.begin());

    return prop_vector;
}

void write_tree(std::ostream &out, const node root, const std::unordered_map<node, std::string>& names)
{
    std::unordered_set<node> visited;
    std::stack<node> s;
    s.push(root);
    while(!s.empty()) {
        node cur = s.top();
        if(cur->is_leaf()) {
            auto iter = names.find(cur);
            out << iter->second;
            visited.insert(cur);
            s.pop();
            continue;
        }
        if(!visited.count(cur->child1->node)) {
            out << "(";
            s.push(cur->child1->node);
            continue;
        } else if(!visited.count(cur->child2->node)) {
            out << ":" << cur->child1->length << ",";
            s.push(cur->child2->node);
            continue;
        }
        out << ":" << cur->child2->length << ")";
        visited.insert(cur);
        s.pop();
    }
    out << ";\n";
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_UTIL_HPP
