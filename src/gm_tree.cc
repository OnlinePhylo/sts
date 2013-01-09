/// \file gm_tree.cc
/// \brief GM_tree implementation

#include "gm_tree.h"
#include "node.h"
#include <algorithm>
#include <cassert>
#include <stack>
#include <stdexcept>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

using namespace std;
using bpp::Node;
using bpp::TreeTemplate;
using sts::particle::Node_ptr;

namespace sts
{
namespace guidedmerge
{

template<typename K, typename V>
V get_with_default(const unordered_map<K, V>& m, const K& k)
{
    if(m.find(k) == end(m))
        return V();
    else
        return m.at(k);
}

// GM_tree

/// \brief Add \c node to the set of nodes adjacent to \c other
void GM_tree::add_node_to(GM_node_ptr node, GM_node_ptr other)
{
    assert(node != nullptr);
    assert(other != nullptr);
    if(node->is_leaf())
        leaves.insert(node);
    adjacent_nodes[other].insert(node);
}

/// \brief add an edge from \c n1 to \c n2
void GM_tree::add_edge(GM_node_ptr n1, GM_node_ptr n2)
{
    assert(n1 != nullptr);
    assert(n2 != nullptr);
    add_node_to(n1, n2);
    add_node_to(n2, n1);
}

/// \brief Remove \c node from set of nodes adjacent to \c other
void GM_tree::remove_node_from(GM_node_ptr node, GM_node_ptr other)
{
    assert(node != nullptr);
    assert(other != nullptr);

    adjacent_nodes[other].erase(node);
    if(adjacent_nodes[other].empty()) {
        if(other->is_leaf()) leaves.erase(other);
        adjacent_nodes.erase(other);
    }
}

/// \brief Remove the edge between \c n1 and \c n2
void GM_tree::remove_edge(GM_node_ptr n1, GM_node_ptr n2)
{
    assert(n1 != nullptr);
    assert(n2 != nullptr);
    remove_node_from(n1, n2);
    remove_node_from(n2, n1);
}

/// \brief Remove \c node from the tree
void GM_tree::remove_node(GM_node_ptr node)
{
    assert(node != nullptr);
    if(node->is_leaf())
        leaves.erase(node);
    for(auto &o : adjacent_nodes[node])
        remove_node_from(node, o);
    adjacent_nodes.erase(node);
}

/// \brief Merge two leaves, eliminating all nodes between.
GM_node_ptr GM_tree::merge(GM_node_ptr n1, GM_node_ptr n2, Node_ptr payload) {
    if(payload == nullptr) payload = make_shared<sts::particle::Node>(nullptr);
    if(!(n1->is_leaf() && n2->is_leaf()))
        throw std::runtime_error("Can only merge leaves");
    auto to_remove = find_path(n1, n2);
    GM_node_ptr merged = make_shared<GM_node>(payload);
    unordered_set<GM_node_ptr> new_adjacent; // Nodes that will be adjacent
                                             // to `merged`
    for(auto n : to_remove) {
        auto a = get_with_default(adjacent_nodes, n);
        new_adjacent.insert(begin(a), end(a));
        remove_node(n);
    }

    // Remove all removed nodes from new_adjacent
    for(auto n : to_remove)
        new_adjacent.erase(n);

    // Add edges from the new merge node to all in new_adjacent
    for(auto n : new_adjacent)
        add_edge(merged, n);

    return merged;
}

/// \brief Find a path between nodes \c n1 and \c n2

/// \param n1 First node
/// \param n2 Second node
/// \return Path from \c n1 to \c n2, including both
/// \throw No_path
std::vector<GM_node_ptr> GM_tree::find_path(GM_node_ptr n1, GM_node_ptr n2) const throw (No_path)
{
    typedef vector<GM_node_ptr> Path;
    typedef pair<GM_node_ptr,Path> PPath;
    unordered_set<GM_node_ptr> seen;
    stack<PPath> s;
    s.push(PPath(n1, {n1}));
    while(!s.empty()) {
        PPath i = s.top();
        s.pop();
        for(auto &n : get_with_default(adjacent_nodes, i.first)) {
            if(seen.find(n) != seen.end()) continue;
            seen.insert(n);
            Path p = i.second;
            p.push_back(n);
            if(n == n2)
                return p;
            s.emplace(n, p);
        }
    }
    throw No_path("no path found between nodes");
}

unordered_set<GM_node_ptr> GM_tree::adjacent_via(const GM_node_ptr& node, const GM_node_ptr& via) const
{
    unordered_set<GM_node_ptr> r = get_with_default(adjacent_nodes, node);
    r.erase(via);
    return r;
}

/// \brief Find all merges between leaves within an RF distance of \c k from the current GM_tree
///
/// \param k RF distance
/// \return Set of nodes which may be merged
unordered_set<pair<GM_node_ptr,GM_node_ptr>> GM_tree::find_k_distance_merges(const size_t k) const
{
    assert(!leaves.empty());
    assert(k > 1);
    GM_node_ptr f = *begin(leaves);
    unordered_set<GM_node_ptr> seen;
    seen.insert(f);
    stack<GM_node_ptr> s;
    s.push(f);

    unordered_map<GM_node_ptr, unordered_map<GM_node_ptr, int>> distances;
    unordered_set<pair<GM_node_ptr, GM_node_ptr>> merges;

    while(!s.empty())
    {
        GM_node_ptr cur = s.top();
        s.pop();

        for(auto &n : get_with_default(adjacent_nodes, cur)) {
            if(seen.count(n)) continue;

            for(auto &p : distances[cur]) {
                if(p.second >= k) continue;

                distances[n][p.first] = distances[p.first][n] = p.second + 1;
                if(n->is_leaf() && p.first->is_leaf() && p.second + 1 == k)
                    merges.emplace(n, p.first);
            }

            distances[cur][n] = distances[n][cur] = 1;
            seen.insert(n);
            s.push(n);
        }
    }

    return merges;
}

/// \brief Convert a GM_tree to a newick-format string
///
/// \return a newick-format string representing the GM_tree.
/// Branch lengths are not specified; tree is arbitrary rooted on the leaf lexicographic min name.
std::string GM_tree::to_newick_string(const unordered_map<Node_ptr,string>& name_map) const
{
    auto get_name = [&name_map](const GM_node_ptr p) -> string {
        if(p->node && name_map.count(p->node)) return name_map.at(p->node);
        return std::to_string((size_t)p.get());
    };

    auto ptr_comp = [](const GM_node_ptr& n1, const GM_node_ptr& n2) {
        return n1->node < n2->node;
    };

    // The GM_tree is stored unrooted.
    // We root using the leaf with lexicographic min taxon_name as an outgroup
    GM_node_ptr m = *min_element(begin(leaves), end(leaves), ptr_comp);

    // Start building the tree
    Node* root = new Node();
    std::unique_ptr<TreeTemplate<Node>> tree(new TreeTemplate<Node>(root));
    root->addSon(new Node(get_name(m)));

    stack<GM_node_ptr> gm_nodes;      // Nodes that have yet to be added to the Bio++ tree
    stack<Node*> parents;        // bpp parent corresponding to node in `s`
    unordered_set<GM_node_ptr> seen;  // Nodes added to the tree
    seen.insert(m);

    // Add all nodes adjacent to `m` as children of `root`
    for(auto &child : get_with_default(adjacent_nodes, m)) {
        gm_nodes.push(child);
        parents.push(root);
    }

    // Add remaining nodes
    while(!gm_nodes.empty()) {
        assert(parents.size() == gm_nodes.size());
        GM_node_ptr n = gm_nodes.top();
        gm_nodes.pop();
        Node* parent = parents.top();
        parents.pop();

        if(seen.count(n))
            continue;
        else
            seen.insert(n);

        Node* new_node = n->is_leaf() ? new Node(get_name(n)) : new Node();
        parent->addSon(new_node);

        for(auto &child : get_with_default(adjacent_nodes, n)) {
            gm_nodes.push(child);
            parents.push(new_node);
        }
    }

    return bpp::TreeTemplateTools::treeToParenthesis(*tree);
}

/// Generate a GM_tree from a TreeTemplate
///
/// \param tree Tree
/// \param name_map Map to fill with node names
GM_tree of_treetemplate(const TreeTemplate<Node>* tree, unordered_map<string,Node_ptr>& name_map)
{
    GM_tree gm;
    const vector<const Node*> nodes = tree->getNodes();
    unordered_map<const Node*,GM_node_ptr> clade_map;
    auto get_node = [&clade_map,&name_map](const Node* n) {
        if(!clade_map.count(n)) {
            Node_ptr np = nullptr;
            if(n->isLeaf()) {
                np = make_shared<sts::particle::Node>(nullptr);
                name_map[n->getName()] = np;
            }
            clade_map[n] = make_shared<GM_node>(np);
        }
        return clade_map[n];
    };
    for(const Node* n : nodes) {
        GM_node_ptr p = get_node(n);
        for(int i = 0; i < n->getNumberOfSons(); ++i) {
            const Node* son = n->getSon(i);
            gm.add_edge(p, get_node(son));
        }
    }
    return gm;
}

/// \brief Generate a GM_tree from a newick string.
///
/// \param nwk Tree in newick format
/// \param name_map Map to fill with node names
GM_tree GM_tree::of_newick_string(const string& nwk, unordered_map<string,Node_ptr>& name_map)
{
    const unique_ptr<const TreeTemplate<Node>> tree(bpp::TreeTemplateTools::parenthesisToTree(nwk));
    return of_treetemplate(tree.get(), name_map);
}

/// \brief Generate a GM_tree from a file containing a newick tree
///
/// \param path Path to tree, in newick format
/// \param name_map Map to fill with node names
GM_tree GM_tree::of_newick_path(const string& path, unordered_map<string,Node_ptr>& name_map)
{
    bpp::Newick newick;
    const unique_ptr<const TreeTemplate<Node>> tree(newick.read(path));
    return of_treetemplate(tree.get(), name_map);
}

}
}
