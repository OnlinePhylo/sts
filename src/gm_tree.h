/// \file gm_tree.h
/// \brief Guided merge tree declarations
#ifndef STS_GUIDEDMERGE_GM_NODE_H
#define STS_GUIDEDMERGE_GM_NODE_H

#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "node_ptr.h"
#include "util.h"

namespace sts
{
namespace guidedmerge
{

/// Exception raised when there is no path between two nodes
class No_path: virtual public std::exception
{
protected:
    std::string message;
public:
    No_path(const std::string& text): message(text) {};
    const char* what() const throw() { return message.c_str(); };
};

/// \brief Node in a GM_tree
struct GM_node
{
    /// \brief Default constructor. Creates an internal node
    GM_node() {};
    GM_node(const sts::particle::Node_ptr& n) : node(n) {};

    /// The STS node associated with this GM_node
    /// Only present for leaves
    sts::particle::Node_ptr node;
    /// Is this node a leaf?
    inline bool is_leaf() const { return static_cast<bool>(node); };
};

/// \brief Shared pointer to a GM_node
typedef std::shared_ptr<GM_node> GM_node_ptr;

/// \brief Guided merge tree
///
/// Stores an unrooted tree of guided merges, in the form of adjacency lists between nodes.
class GM_tree
{
public:
    void add_edge(GM_node_ptr n1, GM_node_ptr n2);
    void remove_edge(GM_node_ptr n1, GM_node_ptr n2);
    void remove_node(GM_node_ptr node);
    GM_node_ptr merge(GM_node_ptr n1, GM_node_ptr n2);
    std::vector<GM_node_ptr> find_path(GM_node_ptr n1, GM_node_ptr n2) const throw (No_path);
    std::unordered_set<GM_node_ptr> adjacent_via(const GM_node_ptr& node, const GM_node_ptr& via) const;
    std::unordered_set<std::pair<GM_node_ptr,GM_node_ptr>> find_k_distance_merges(const size_t k) const;

    std::string to_newick_string(const std::unordered_map<sts::particle::Node_ptr,std::string>& name_map) const;
    static GM_tree of_newick_path(const std::string& path, std::unordered_map<std::string,sts::particle::Node_ptr>& name_map);
    static GM_tree of_newick_string(const std::string& s, std::unordered_map<std::string,sts::particle::Node_ptr>& name_map);
private:
    void add_node_to(GM_node_ptr node, GM_node_ptr other);
    void remove_node_from(GM_node_ptr node, GM_node_ptr other);

    /// Adjacency list. Maps from a node pointer to adjacent nodes
    std::unordered_map<GM_node_ptr,std::unordered_set<GM_node_ptr>> adjacent_nodes;

    /// All leaves
    std::unordered_set<GM_node_ptr> leaves;
};

}
}

namespace std
{

/// Generic hash function for pairs
template<typename A, typename B>
struct hash<std::pair<A, B>> {
private:
   const hash<A> ah;
   const hash<B> bh;
public:
   hash() : ah(), bh() {}
   size_t operator()(const std::pair<A, B>& p) const {
        size_t result = 0;
        sts::util::hash_combine(result, ah(p.first));
        sts::util::hash_combine(result, bh(p.second));
        return result;
   }
};

}

#endif // STS_GUIDEDMERGE_GM_NODE_H
