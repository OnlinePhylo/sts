#include "gm_tree.h"
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

    throw std::runtime_error("Not implemented.");
}

GM_tree GM_tree::of_newick_path(const std::string& path)
{
    GM_tree gm;
    bpp::Newick newick;
    const unique_ptr<bpp::TreeTemplate<bpp::Node>> tree(newick.read(path));
    const vector<bpp::Node*> nodes = tree->getNodes();
    unordered_map<const bpp::Node*,GM_node_ptr> clade_map;
    auto get_node = [&clade_map](const bpp::Node* n) {
        if(!clade_map.count(n))
            clade_map[n] = make_shared<GM_node>(n->getName(), n->isLeaf());
        return clade_map[n];
    };
    for(const bpp::Node* n : nodes) {
        GM_node_ptr p = get_node(n);
        for(auto &i : n->getSonsId()) {
            const bpp::Node* son = n->getSon(i);
            gm.add_edge(p, get_node(son));
        }
    }
    return gm;
}

}
}
