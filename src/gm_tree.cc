/// \file gm_tree.cc
/// \brief GM_tree implementation

#include "gm_tree.h"
#include "node.h"
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include <algorithm>
#include <cassert>
#include <stack>
#include <stdexcept>

#include <iostream>

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
/// \brief Constructor
GM_tree::GM_tree()
{
}

/// \brief Copy constructor
GM_tree::GM_tree(const GM_tree& other)
{
    std::unordered_map<GM_node*,GM_node*> translation_table;

    // Copy nodes
    for(auto &p : other.all_nodes)
        translation_table[p.first] = create_node(p.first->node);
    for(auto p : other.node_gmnode)
        node_gmnode.emplace(p.first, translation_table.at(p.second));
    for(auto p : other.gmnode_node)
        gmnode_node.emplace(translation_table.at(p.first), p.second);
    for(auto p : other.adjacent_nodes) {
        GM_node* k = translation_table.at(p.first);
        for(auto n : p.second)
            adjacent_nodes[k].insert(translation_table.at(n));
    }

    assert(get_leaf_count() == other.get_leaf_count());
    assert(adjacent_nodes.size() == other.adjacent_nodes.size());
}

/// \brief Assignment operator
GM_tree& GM_tree::operator=(const GM_tree& other)
{
    // Leverage the copy constructor, then steal bits from the copy
    GM_tree tmp(other);
    all_nodes = std::move(tmp.all_nodes);
    adjacent_nodes = std::move(tmp.adjacent_nodes);
    node_gmnode = std::move(tmp.node_gmnode);
    gmnode_node = std::move(tmp.gmnode_node);

    return *this;
}

/// \brief Create a node with given payload
GM_node* GM_tree::create_node(Node_ptr np)
{
    GM_node* n = new GM_node(np);
    all_nodes.emplace(n, unique_ptr<GM_node>(n));
    if(np != nullptr) {
        add_leaf(n);
    }
    return n;
}


/// Add a leaf
///
/// Idempotent
/// \pre \c gm_node is a leaf
/// \param gm_node Node to add
void GM_tree::add_leaf(GM_node* gm_node)
{
    assert(gm_node->is_leaf());
    assert(all_nodes.count(gm_node) == 1);

    node_gmnode[gm_node->node] = gm_node;
    gmnode_node[gm_node] = gm_node->node;
}

/// \brief Add \c node to the set of nodes adjacent to \c other
void GM_tree::add_node_to(GM_node* node, GM_node* other)
{
    assert(node != nullptr);
    assert(other != nullptr);
    if(node->is_leaf())
        add_leaf(node);
    adjacent_nodes[other].insert(node);
}

/// \brief add an edge from \c n1 to \c n2
void GM_tree::add_edge(GM_node* n1, GM_node* n2)
{
    assert(n1 != nullptr);
    assert(n2 != nullptr);
    add_node_to(n1, n2);
    add_node_to(n2, n1);
}

/// \brief Remove \c node from set of nodes adjacent to \c other
///
/// \param remove_lonely When \c true: If \c other is no longer adjacent to any nodes after removing \c node, delete it.
void GM_tree::remove_node_from(GM_node* node, GM_node* other, bool remove_lonely)
{
    assert(node != nullptr);
    assert(other != nullptr);

    adjacent_nodes[other].erase(node);
    if(remove_lonely && adjacent_nodes[other].empty()) {
        remove_node(other);
    }
}

/// \brief Remove the edge between \c n1 and \c n2
///
/// \param remove_lonely When \c true: if removing the edge between \c n1 and \c n2 causes either to no longer have any
/// adjacent nodes, delete the lonely node(s).
void GM_tree::remove_edge(GM_node* n1, GM_node* n2, bool remove_lonely)
{
    assert(n1 != nullptr);
    assert(n2 != nullptr);
    remove_node_from(n1, n2, remove_lonely);
    remove_node_from(n2, n1, remove_lonely);
}

/// \brief Remove \c node from the tree
///
/// \param remove_lonely When \c true : If any of the nodes adjacent to \c node have no additional adjacent nodes, remove them.
void GM_tree::remove_node(GM_node* node, bool remove_lonely)
{
    assert(node != nullptr);
    if(node->is_leaf()) {
        node_gmnode.erase(node->node);
        gmnode_node.erase(node);
    }
    for(auto &o : get_with_default(adjacent_nodes, node)) {
        remove_node_from(node, o, remove_lonely);
    }
    adjacent_nodes.erase(node);
    all_nodes.erase(node);
}

/// \brief Merge two leaves, eliminating all nodes between.
void GM_tree::merge(Node_ptr n1, Node_ptr n2, Node_ptr payload) throw (No_path)
{
    size_t lcount = get_leaf_count();
    GM_node *gm1 = node_gmnode.at(n1), *gm2 = node_gmnode.at(n2);
    if(!(gm1->is_leaf() && gm2->is_leaf()))
        throw std::runtime_error("Can only merge leaves");
    auto to_remove = find_path(gm1, gm2);
    GM_node* merged = create_node(payload);
    assert(merged->is_leaf());
    unordered_set<GM_node*> new_adjacent; // Nodes that will be adjacent
                                          // to `merged`

    for(GM_node* n : to_remove) {
        auto a = get_with_default(adjacent_nodes, n);
        new_adjacent.insert(begin(a), end(a));
        remove_node(n, false);
    }

    // Remove all removed nodes from new_adjacent
    for(GM_node* n : to_remove)
        new_adjacent.erase(n);

    // Add edges from the new merge node to all in new_adjacent
    for(GM_node* n : new_adjacent) {
        assert(all_nodes.count(n) == 1);
        assert(all_nodes.count(merged) == 1);
        add_edge(merged, n);
    }

    assert(get_leaf_count() == lcount - 1);
}

/// \brief Find a path between nodes \c n1 and \c n2

/// \param n1 First node
/// \param n2 Second node
/// \return Path from \c n1 to \c n2, including both
/// \throw No_path
std::vector<GM_node*> GM_tree::find_path(GM_node* n1, GM_node* n2) const throw (No_path)
{
    typedef vector<GM_node*> Path;
    typedef pair<GM_node*,Path> PPath;
    unordered_set<GM_node*> seen;
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

/// \pre \c n1 and \c n2 are nodes present in the GM_tree
/// \return Whether a path exists between \c n1 and \c n2
bool GM_tree::path_exists(Node_ptr n1, Node_ptr n2) const
{
    try {
        return find_path(node_gmnode.at(n1), node_gmnode.at(n2)).size() > 2;
    } catch (No_path& e) {
        return false;
    }
}

/// RF distance between \c n1 and \c n2 in the GM_tree
///
/// \pre \c n1 and \c n2 are nodes present in the GM_tree
///
/// \throw No_path No path exists between n1 and n2
size_t GM_tree::rf_distance(Node_ptr n1, Node_ptr n2) const throw (No_path)
{
    return find_path(node_gmnode.at(n1), node_gmnode.at(n2)).size() - 2;
}

/// \brief Find all merges between leaves within an RF distance of \c k from the current GM_tree
///
/// \note The python implementation \c sts.py includes pendant edges in \c k, while this function does not.
///
/// \param k Number of splits that may disagree with this GM_tree
/// \return Node pairs which, when joined, disagree with \c k splits in the GM_tree.
unordered_set<pair<Node_ptr,Node_ptr>> GM_tree::find_k_distance_merges(const size_t k) const
{
    const size_t k_adj = k + 2;
    assert(!node_gmnode.empty());
    GM_node* f = begin(gmnode_node)->first;
    unordered_set<GM_node*> seen;
    seen.insert(f);
    stack<GM_node*> s;
    s.push(f);

    unordered_map<GM_node*, unordered_map<GM_node*, size_t>> distances;
    unordered_set<pair<Node_ptr,Node_ptr>> merges;

    while(!s.empty())
    {
        GM_node* cur = s.top();
        s.pop();

        for(auto &n : get_with_default(adjacent_nodes, cur)) {
            if(seen.count(n)) continue;

            for(auto &p : distances[cur]) {
                if(p.second >= k_adj) continue;

                distances[n][p.first] = distances[p.first][n] = p.second + 1;
                if(n->is_leaf() && p.first->is_leaf() && p.second + 1 == k_adj)
                    merges.emplace(gmnode_node.at(n), gmnode_node.at(p.first));
            }

            distances[cur][n] = distances[n][cur] = 1;
            seen.insert(n);
            s.push(n);
        }
    }

    return merges;
}

/// \brief Convert a GM_tree to a Bio++ TreeTemplate
///
/// Branch lengths are not specified; tree is arbitrary rooted on the leaf lexicographic min name.
TreeTemplate<Node>* GM_tree::to_treetemplate(const unordered_map<Node_ptr,string>& name_map) const
{
    auto get_name = [&name_map](const GM_node* p) -> string {
        if(p->node && name_map.count(p->node)) return name_map.at(p->node);
        return std::to_string((size_t)p);
    };

    auto ptr_comp = [](const pair<GM_node*,Node_ptr>& n1, const pair<GM_node*,Node_ptr>& n2) {
        return n1.first->node < n2.first->node;
    };

    // The GM_tree is stored unrooted.
    // We root using the leaf with lexicographic min taxon_name as an outgroup
    GM_node* m = min_element(begin(gmnode_node), end(gmnode_node), ptr_comp)->first;

    // Start building the tree
    Node* root = new Node();
    std::unique_ptr<TreeTemplate<Node>> tree(new TreeTemplate<Node>(root));
    root->addSon(new Node(get_name(m)));

    stack<GM_node*> gm_nodes;      // Nodes that have yet to be added to the Bio++ tree
    stack<Node*> parents;        // bpp parent corresponding to node in `s`
    unordered_set<GM_node*> seen;  // Nodes added to the tree
    seen.insert(m);

    // Add all nodes adjacent to `m` as children of `root`
    for(auto &child : get_with_default(adjacent_nodes, m)) {
        gm_nodes.push(child);
        parents.push(root);
    }

    // Add remaining nodes
    while(!gm_nodes.empty()) {
        assert(parents.size() == gm_nodes.size());
        GM_node* n = gm_nodes.top();
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

    return tree.release();
}

/// \brief Convert a GM_tree to a newick-format string
///
/// \return a newick-format string representing the GM_tree.
/// Branch lengths are not specified; tree is arbitrary rooted on the leaf lexicographic min name.
std::string GM_tree::to_newick_string(const unordered_map<Node_ptr,string>& name_map) const
{
    unique_ptr<TreeTemplate<Node>> tree(to_treetemplate(name_map));
    return bpp::TreeTemplateTools::treeToParenthesis(*tree);
}

/// Generate a GM_tree from a TreeTemplate
///
/// \param tree Tree
/// \param name_map Map to fill with node names
GM_tree GM_tree::of_treetemplate(const TreeTemplate<Node>* tree, unordered_map<string,Node_ptr>& name_map)
{
    GM_tree gm;
    const vector<const Node*> nodes = tree->getNodes();
    unordered_map<const Node*,GM_node*> clade_map;
    auto get_node = [&clade_map,&name_map,&gm](const Node* n) -> GM_node* {
        if(!clade_map.count(n)) {
            Node_ptr np = nullptr;
            if(n->isLeaf()) {
                np = make_shared<sts::particle::Node>(nullptr);
                name_map[n->getName()] = np;
            }
            clade_map[n] = gm.create_node(np);
        }
        return clade_map[n];
    };
    for(const Node* n : nodes) {
        GM_node* p = get_node(n);
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
