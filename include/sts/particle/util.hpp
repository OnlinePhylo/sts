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

int tree_count(const std::vector<std::shared_ptr<phylo_node>> &);
std::vector<std::shared_ptr<phylo_node>> uncoalesced_nodes(std::shared_ptr<phylo_particle> pp,
                                      std::vector<std::shared_ptr<phylo_node>> leaf_nodes);

void write_tree(std::ostream &out, const std::shared_ptr<phylo_node> root, const std::vector<std::string> &names);

/// Find the number of trees (that is, trees consisting of more than one node) from a collection of uncoalesced nodes.
/// \param uncoalesced The uncoalesced nodes.
/// \return The count.
int tree_count(const std::vector<std::shared_ptr<phylo_node>> &uncoalesced)
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
std::vector<std::shared_ptr<phylo_node>> uncoalesced_nodes(const std::shared_ptr<phylo_particle> pp, const std::vector<std::shared_ptr<phylo_node>> leaf_nodes)
{
    // Our set of phylo nodes that can be used in proposal.
    std::unordered_set<std::shared_ptr<phylo_node>> proposal_set;
    // The nodes that have already been coalesced, to be removed later.
    std::unordered_set<std::shared_ptr<phylo_node>> coalesced;
    // Insert all of the leaf nodes into the proposal set.
    proposal_set.insert(leaf_nodes.begin(), leaf_nodes.end());
    // Walk back to predecessor particles, adding root nodes to
    // proposal_set and collecting coalesced nodes in `coalesced`.
    for(std::shared_ptr<phylo_particle> cur = pp->predecessor; cur != NULL; cur = cur->predecessor) {
        // Skip if the particle is \perp.
        if(cur->node == NULL) continue;
        // Skip if we've already processed this subtree, such that it's already found in coalesced.
        if(coalesced.find(cur->node) != coalesced.end()) continue;
        // Insert this active root node to the proposal set.
        proposal_set.insert(cur->node);
        // Recursively add all descendants of the root nodes to the coalesced set using a std::stack.
        std::stack<std::shared_ptr<phylo_node>> s;
        s.push(cur->node);
        while(s.size() > 0) {
            std::shared_ptr<phylo_node> n = s.top();
            s.pop();
            if(n->is_leaf()) continue;	// leaf node, nothing more to do.
            coalesced.insert(n->child1->node);
            coalesced.insert(n->child2->node);
            s.push(n->child1->node);
            s.push(n->child2->node);
        }
    }

    std::vector<std::shared_ptr<phylo_node>> pvec(proposal_set.begin(), proposal_set.end());
    std::vector<std::shared_ptr<phylo_node>> cvec(coalesced.begin(), coalesced.end());
    sort(pvec.begin(), pvec.end());
    sort(cvec.begin(), cvec.end());

    // The set difference of available (i.e. proposal_set) and coalesced nodes yields the final proposal set; store it
    // in prop_vector.
    std::vector<std::shared_ptr<phylo_node>> prop_vector(proposal_set.size() + coalesced.size());
    // UGH: std::set_difference requires an ordered container class
    // AG: that's the only way to do a set difference efficiently, right?
    auto last_ins = set_difference(pvec.begin(), pvec.end(), cvec.begin(),
                                   cvec.end(), prop_vector.begin());
    prop_vector.resize(last_ins - prop_vector.begin());

    return prop_vector;
}

std::shared_ptr<phylo_particle> phylo_of_tree(std::shared_ptr<likelihood::online_calculator> calc, bpp::TreeTemplate<bpp::Node> &tree)
{
    std::shared_ptr<phylo_particle> particle = std::make_shared<phylo_particle>();
    particle->node = phylo_node::of_tree(calc, tree);
    if(particle->node->is_leaf())
        return particle;

    std::shared_ptr<phylo_particle> prev = particle;
    std::stack<std::shared_ptr<phylo_node>> node_stack;
    node_stack.push(particle->node->child1->node);
    node_stack.push(particle->node->child2->node);
    while(!node_stack.empty()) {
        std::shared_ptr<phylo_particle> cur = std::make_shared<phylo_particle>();
        cur->node = node_stack.top();
        node_stack.pop();
        prev->predecessor = cur;
        if(!cur->node->is_leaf()) {
            node_stack.push(cur->node->child1->node);
            node_stack.push(cur->node->child2->node);
        }
        prev = cur;
    }

    return particle;
}

std::shared_ptr <phylo_particle> phylo_of_newick_string(std::shared_ptr<likelihood::online_calculator> calc, std::string &tree_string)
{
    bpp::TreeTemplate<bpp::Node> *tree = bpp::TreeTemplateTools::parenthesisToTree(tree_string);
    std::shared_ptr<phylo_particle> node = phylo_of_tree(calc, *tree);
    delete tree;
    return node;
}

static void check_visited(std::vector<bool> &visited, int id)
{
    // ensure visited has enough space allocated to store the id
    // if not, resize it large enough and leave some wiggle to prevent frequent resizings
    if(id >= visited.size()) {
        visited.resize(id + 100);
    }
}

static bool visited_id(std::vector<bool> &visited, int id)
{
    check_visited(visited, id);
    return visited[id];
}

static void set_visited_id(std::vector<bool> &visited, int id)
{
    check_visited(visited, id);
    visited[id] = true;
}

void write_tree(std::ostream &out, const std::shared_ptr<phylo_node> root, const std::vector<std::string> &names)
{
    std::vector<bool> visited;
    std::stack<std::shared_ptr<phylo_node>> s;
    s.push(root);
    while(!s.empty()) {
        std::shared_ptr<phylo_node> cur = s.top();
        if(cur->is_leaf()) {
            out << names[cur->id];
            set_visited_id(visited, cur->id);
            s.pop();
            continue;
        }
        if(!visited_id(visited, cur->child1->node->id)) {
            out << "(";
            s.push(cur->child1->node);
            continue;
        } else if(!visited_id(visited, cur->child2->node->id)) {
            out << ":" << cur->child1->length << ",";
            s.push(cur->child2->node);
            continue;
        }
        out << ":" << cur->child2->length << ")";
        set_visited_id(visited, cur->id);
        s.pop();
    }
    out << ";\n";
}

} // namespace particle
} // namespace sts

#endif // STS_PARTICLE_UTIL_HPP
