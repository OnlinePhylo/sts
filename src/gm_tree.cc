#include "gm_tree.h"
#include <algorithm>
#include <cassert>
#include <stack>
#include <stdexcept>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

using namespace std;

namespace sts
{
namespace guidedmerge
{

// GM_tree
void GM_tree::add_node_to(GM_node_ptr node, GM_node_ptr other)
{
    assert(node != nullptr);
    assert(other != nullptr);
    if(node->is_leaf)
        leaves.insert(node);
    if(adjacent_nodes.find(other) == adjacent_nodes.end())
        adjacent_nodes[other] = {node};
    else
        adjacent_nodes[other].insert(node);
}

void GM_tree::add_edge(GM_node_ptr n1, GM_node_ptr n2)
{
    assert(n1 != nullptr);
    assert(n2 != nullptr);
    add_node_to(n1, n2);
    add_node_to(n2, n1);
}

void GM_tree::remove_node_from(GM_node_ptr node, GM_node_ptr other)
{
    assert(node != nullptr);
    assert(other != nullptr);

    adjacent_nodes[other].erase(node);
    if(adjacent_nodes[other].empty()) {
        if(other->is_leaf) leaves.erase(other);
        adjacent_nodes.erase(other);
    }
}

void GM_tree::remove_edge(GM_node_ptr n1, GM_node_ptr n2)
{
    assert(n1 != nullptr);
    assert(n2 != nullptr);
    remove_node_from(n1, n2);
    remove_node_from(n2, n1);
}

void GM_tree::remove_node(GM_node_ptr node)
{
    assert(node != nullptr);
    if(node->is_leaf)
        leaves.erase(node);
    for(auto &o : adjacent_nodes[node])
        remove_node_from(node, o);
    adjacent_nodes.erase(node);
}

GM_node_ptr GM_tree::merge(GM_node_ptr n1, GM_node_ptr n2)
{
    if(!(n1->is_leaf && n2->is_leaf))
        throw std::runtime_error("Can only merge leaves");
    auto to_remove = find_path(n1, n2);
    unordered_set<GM_node_ptr> new_adjacent;
    GM_node_ptr merged = make_shared<GM_node>();
    for(auto n : to_remove) {
        auto a = adjacent_nodes[n];
        new_adjacent.insert(begin(a), end(a));
        remove_node(n);
    }

    return merged;
}

std::vector<GM_node_ptr> GM_tree::find_path(GM_node_ptr n1, GM_node_ptr n2)
{
    typedef vector<GM_node_ptr> Path;
    typedef pair<GM_node_ptr,Path> PPath;
    unordered_set<GM_node_ptr> seen;
    stack<PPath> s;
    s.push(PPath(n1, {n1}));
    while(!s.empty()) {
        PPath i = s.top();
        s.pop();
        for(auto &n : adjacent_nodes[i.first]) {
            if(seen.find(n) != seen.end()) continue;
            seen.insert(n);
            Path p = i.second;
            p.push_back(n);
            if(n == n2)
                return p;
            s.emplace(n, p);
        }
    }
    assert(false);
}
unordered_set<GM_node_ptr> GM_tree::adjacent_via(GM_node_ptr node, GM_node_ptr via)
{
    unordered_set<GM_node_ptr> r;
    if(!adjacent_nodes.count(node)) return r;
    r = adjacent_nodes[node];
    r.erase(via);
    return r;
}

unordered_set<pair<GM_node_ptr,GM_node_ptr>> GM_tree::find_k_distance_merges(const size_t k)
{
    assert(!leaves.empty());
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

        for(auto &n : adjacent_nodes[cur]) {
            if(seen.count(n)) continue;

            for(auto &p : distances[cur]) {
                if(p.second >= k) continue;

                distances[n][p.first] = distances[p.first][n] = p.second + 1;
                if(n->is_leaf && p.first->is_leaf && p.second + 1 == k)
                    merges.emplace(n, p.first);
            }

            distances[cur][n] = distances[n][cur] = 1;
            seen.insert(n);
            s.push(n);
        }
    }

    return merges;
}

std::string GM_tree::to_newick_string() const
{
    auto comp = [](const GM_node_ptr& n1, const GM_node_ptr& n2) {
        return n1->taxon_name < n2 ->taxon_name;
    };

    // Root on an "arbitrary" leaf - lexicographic min taxon_name
    GM_node_ptr m = *min_element(begin(leaves), end(leaves), comp);

    // Build the tree
    std::unique_ptr<bpp::Node> root(new bpp::Node());
    root->addSon(new bpp::Node(m->taxon_name));

    stack<GM_node_ptr> gm_nodes;      // Nodes that have yet to be added to the Bio++ tree
    stack<bpp::Node*> parents;        // bpp parent corresponding to node in `s`
    unordered_set<GM_node_ptr> seen;  // Nodes added to the tree
    seen.insert(m);

    // Add all nodes adjacent to `m` as children of `root`
    if(adjacent_nodes.count(m)) {
        for(auto &child : adjacent_nodes.at(m)) {
            gm_nodes.push(child);
            parents.push(root.get());
        }
    }

    while(!gm_nodes.empty()) {
        assert(parents.size() == gm_nodes.size());
        GM_node_ptr n = gm_nodes.top();
        gm_nodes.pop();
        bpp::Node* parent = parents.top();
        parents.pop();

        if(seen.count(n))
            continue;
        else
            seen.insert(n);

        bpp::Node* new_node = n->is_leaf ? new bpp::Node(n->taxon_name) : new bpp::Node("");
        parent->addSon(new_node);

        if(adjacent_nodes.count(n)) {
            for(auto &child : adjacent_nodes.at(n)) {
                gm_nodes.push(child);
                parents.push(new_node);
            }
        }
    }

//    std::unique_ptr<bpp::TreeTemplate<bpp::Node>> tree(new bpp::TreeTemplate<bpp::Node>(root.release());

    return bpp::TreeTemplateTools::nodeToParenthesis(*root);
}

GM_tree of_treetemplate(const bpp::TreeTemplate<bpp::Node>* tree)
{
    GM_tree gm;
    const vector<const bpp::Node*> nodes = tree->getNodes();
    unordered_map<const bpp::Node*,GM_node_ptr> clade_map;
    auto get_node = [&clade_map](const bpp::Node* n) {
        if(!clade_map.count(n))
            clade_map[n] = make_shared<GM_node>(n->hasName() ? n->getName() : "",
                    n->isLeaf());
        return clade_map[n];
    };
    for(const bpp::Node* n : nodes) {
        GM_node_ptr p = get_node(n);
        for(int i = 0; i < n->getNumberOfSons(); ++i) {
            const bpp::Node* son = n->getSon(i);
            gm.add_edge(p, get_node(son));
        }
    }
    return gm;
}

GM_tree GM_tree::of_newick_string(const std::string& nwk)
{
    const unique_ptr<const bpp::TreeTemplate<bpp::Node>> tree(bpp::TreeTemplateTools::parenthesisToTree(nwk));
    return of_treetemplate(tree.get());
}

GM_tree GM_tree::of_newick_path(const std::string& path)
{
    bpp::Newick newick;
    const unique_ptr<const bpp::TreeTemplate<bpp::Node>> tree(newick.read(path));
    return of_treetemplate(tree.get());
}

}
}
