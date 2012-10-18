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
#include <Bpp/Phyl/TreeTemplateTools.h>
#include "smctc.hh"

#include "phylofunc.hh"
#include "hmsbeagle.hh"

using namespace std;

edge::edge(std::shared_ptr<phylo_node> node, double length) : length(length), node(node) {}

phylo_node::phylo_node(std::shared_ptr<online_calculator> calc) : id(-1), calc(calc) {};
phylo_node::phylo_node(const phylo_node &other) : id(other.id), calc(other.calc) {};

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

shared_ptr< phylo_node >
phylo_node::of_tree(shared_ptr< online_calculator > calc, bpp::TreeTemplate< bpp::Node > &tree, int node_number)
{
    shared_ptr< phylo_node > node = make_shared< phylo_node >(calc);
    node->id = node_number;
    if (tree.isLeaf(node_number))
        return node;
    vector< int > children = tree.getSonsId(node_number);
    assert(children.size() == 2);
    node->child1 = edge::of_tree(calc, tree, children[0]);
    node->child2 = edge::of_tree(calc, tree, children[1]);
    return node;
}

shared_ptr< edge >
edge::of_tree(shared_ptr< online_calculator > calc, bpp::TreeTemplate< bpp::Node > &tree, int node_number)
{
    return make_shared< edge >(
        phylo_node::of_tree(calc, tree, node_number),
        tree.getDistanceToFather(node_number));
}

shared_ptr< phylo_node >
phylo_node::of_tree(shared_ptr< online_calculator > calc, bpp::TreeTemplate< bpp::Node > &tree)
{
    return phylo_node::of_tree(calc, tree, tree.getRootId());
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

shared_ptr< phylo_particle >
phylo_particle::of_tree(shared_ptr< online_calculator > calc, bpp::TreeTemplate< bpp::Node > &tree)
{
    shared_ptr< phylo_particle > particle = make_shared< phylo_particle >();
    particle->node = phylo_node::of_tree(calc, tree);
    if (particle->node->is_leaf())
        return particle;

    shared_ptr< phylo_particle > prev = particle;
    stack< shared_ptr< phylo_node > > node_stack;
    node_stack.push(particle->node->child1->node);
    node_stack.push(particle->node->child2->node);
    while (!node_stack.empty()) {
        shared_ptr< phylo_particle > cur = make_shared< phylo_particle >();
        cur->node = node_stack.top();
        node_stack.pop();
        prev->predecessor = cur;
        if (!cur->node->is_leaf()) {
            node_stack.push(cur->node->child1->node);
            node_stack.push(cur->node->child2->node);
        }
        prev = cur;
    }

    return particle;
}

shared_ptr < phylo_particle >
phylo_particle::of_newick_string(shared_ptr< online_calculator > calc, string &tree_string)
{
    bpp::TreeTemplate< bpp::Node > *tree = bpp::TreeTemplateTools::parenthesisToTree(tree_string);
    shared_ptr< phylo_particle > node = phylo_particle::of_tree(calc, *tree);
    delete tree;
    return node;
}

static void
check_visited(vector< bool > &visited, int id)
{
    // ensure visited has enough space allocated to store the id
    // if not, resize it large enough and leave some wiggle to prevent frequent resizings
    if(id >= visited.size()) {
        visited.resize(id + 100);
    }
}

static bool
visited_id(vector< bool > &visited, int id)
{
    check_visited(visited, id);
    return visited[id];
}

static void
set_visited_id(vector< bool > &visited, int id)
{
    check_visited(visited, id);
    visited[id] = true;
}

void
write_tree(ostream &out, const shared_ptr< phylo_node > root, const vector< string > &names)
{
    vector< bool > visited;
    stack< shared_ptr< phylo_node > > s;
    s.push(root);
    while (!s.empty()) {
        shared_ptr< phylo_node > cur = s.top();
        if (cur->is_leaf()) {
            out << names[cur->id];
            set_visited_id(visited, cur->id);
            s.pop();
            continue;
        }
        if (!visited_id(visited, cur->child1->node->id)) {
            out << "(";
            s.push(cur->child1->node);
            continue;
        } else if (!visited_id(visited, cur->child2->node->id)) {
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
