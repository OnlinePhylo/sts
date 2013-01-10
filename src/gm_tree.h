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

#include <Bpp/Phyl/TreeTemplate.h>

#include "node_ptr.h"
#include "util.h"

namespace sts
{
namespace guidedmerge
{

/// Exception raised when there is no path between two nodes
class No_path
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

/// \brief Guided merge tree
///
/// Stores an unrooted tree of guided merges, in the form of adjacency lists between nodes.
class GM_tree
{
public:
    GM_tree();
    GM_tree(const GM_tree& other);
    GM_tree& operator=(const GM_tree& other);

    void merge(sts::particle::Node_ptr n1, sts::particle::Node_ptr n2, sts::particle::Node_ptr payload);
    bool path_exists(sts::particle::Node_ptr n1, sts::particle::Node_ptr n2) const;
    size_t rf_distance(sts::particle::Node_ptr n1, sts::particle::Node_ptr n2) const;
    std::unordered_set<std::pair<sts::particle::Node_ptr,sts::particle::Node_ptr>> find_k_distance_merges(const size_t k) const;

    /// \brief Number of leaves is the GM_tree
    size_t get_leaf_count() const { return node_gmnode.size(); };
    std::vector<sts::particle::Node_ptr> leaves() const;

    // Newick-related
    bpp::TreeTemplate<bpp::Node>* to_treetemplate(const std::unordered_map<sts::particle::Node_ptr,std::string>& name_map) const;
    std::string to_newick_string(const std::unordered_map<sts::particle::Node_ptr,std::string>& name_map) const;
    static GM_tree of_newick_path(const std::string& path, const std::unordered_map<std::string,sts::particle::Node_ptr>& name_map);
    static GM_tree of_newick_string(const std::string& nwk, const std::unordered_map<std::string,sts::particle::Node_ptr>& name_map);
    static GM_tree of_treetemplate(const bpp::TreeTemplate<bpp::Node>& tree, const std::unordered_map<std::string,sts::particle::Node_ptr>& name_map);

private:
    void add_leaf(GM_node* node);
    void add_edge(GM_node* n1, GM_node* n2);
    void remove_edge(GM_node* n1, GM_node* n2, bool remove_lonely=true);
    void remove_node(GM_node* node, bool remove_lonely=true);
    GM_node* create_node(sts::particle::Node_ptr np=nullptr);
    std::vector<GM_node*> find_path(GM_node* n1, GM_node* n2) const;
    void add_node_to(GM_node* node, GM_node* other);
    void remove_node_from(GM_node* node, GM_node* other, bool remove_lonely=true);

    /// Adjacency list. Maps from a node pointer to adjacent nodes
    std::unordered_map<GM_node*,std::unordered_set<GM_node*>> adjacent_nodes;

    // Maps between GM nodes and Node_ptrs
    std::unordered_map<sts::particle::Node_ptr,GM_node*> node_gmnode;
    std::unordered_map<GM_node*,sts::particle::Node_ptr> gmnode_node;

    /// All nodes
    std::unordered_map<GM_node*,std::unique_ptr<GM_node>> all_nodes;
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
