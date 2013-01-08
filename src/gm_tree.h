#ifndef STS_GUIDEDMERGE_GM_NODE_H
#define STS_GUIDEDMERGE_GM_NODE_H

#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace sts
{
namespace guidedmerge
{

// Boost
template<typename T>
void hash_combine(size_t& seed, const T& t)
{
    std::hash<T> h;
    seed ^= h(t) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct GM_node
{
    GM_node(const std::string& taxon_name, bool is_leaf) :
        taxon_name(taxon_name), is_leaf(is_leaf) {};
    std::string taxon_name;
    bool is_leaf;
};

typedef std::shared_ptr<GM_node> GM_node_ptr;

class GM_tree
{
public:
    void add_edge(GM_node_ptr n1, GM_node_ptr n2);
    void remove_edge(GM_node_ptr n1, GM_node_ptr n2);
    void remove_node(GM_node_ptr node);
    std::vector<GM_node_ptr> find_path(GM_node_ptr n1, GM_node_ptr n2);
    std::unordered_set<GM_node_ptr> adjacent_via(GM_node_ptr node, GM_node_ptr via);
    std::unordered_set<std::pair<GM_node_ptr,GM_node_ptr>> find_k_distance_merges(const size_t k);

    static GM_tree of_newick_path(const std::string& path);
private:
    void add_node_to(GM_node_ptr node, GM_node_ptr other);
    void remove_node_from(GM_node_ptr node, GM_node_ptr other);

    std::unordered_map<GM_node_ptr,std::unordered_set<GM_node_ptr>> adjacent_nodes;

    std::unordered_set<GM_node_ptr> leaves;
};

}
}

namespace std
{
template<>
struct hash<sts::guidedmerge::GM_node>
{
    size_t operator()(const sts::guidedmerge::GM_node& node)
    {
        return (size_t)&node;
    };
};

template<>
struct less<sts::guidedmerge::GM_node>
{
    bool operator()(const sts::guidedmerge::GM_node& n1,
                    const sts::guidedmerge::GM_node& n2)
    {
        hash<sts::guidedmerge::GM_node> h;
        return h(n1) < h(n2);
    };

};

template<typename A, typename B>
struct hash<std::pair<A, B>> {
private:
   const hash<A> ah;
   const hash<B> bh;
public:
   hash() : ah(), bh() {}
   size_t operator()(const std::pair<A, B>& p) const {
        size_t result = 0;
        sts::guidedmerge::hash_combine(result, ah(p.first));
        sts::guidedmerge::hash_combine(result, bh(p.second));
        return result;
   }
};
}

#endif // STS_GUIDEDMERGE_GM_NODE_H
